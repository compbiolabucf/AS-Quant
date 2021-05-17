import matplotlib
#matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import bisect
from bisect import bisect_left
import csv
import methods
import time
import sys
import os, glob
import configparser
import xlsxwriter
from matplotlib import rcParams
plt.rc('legend',**{'fontsize':8})

rcParams.update({
    'font.family':'arial',
    })

y_limit = 0

def bi_contains(lst, item):
    return bisect_left(lst, item)

def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, '*')))]

def Generate_read_coverate_plot(ax, pathin, sample, chrom, geneID, start, end, startAll, endAll, group_name, number):
	bam_file_reader= open(pathin+'/'+sample+'/'+chrom+".txt", "rt")
	bam_read = csv.reader(bam_file_reader, delimiter="\t")
	bam_list = list(bam_read)
	position_row = [int(bam_list[i][1]) for i in range(len(bam_list))]

	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)

	pos1 = bi_contains(position_row, startAll)
	pos2 = bi_contains(position_row, endAll)
	if(int(bam_list[pos2][1]) != endAll):
		pos2 = pos2 - 1

	p = []
	c = []
	read = 0
	length = endAll - startAll + 1
	
	for t in range(length):
		p.append(t+startAll)
		c.append(0)
		
	for t in range(pos1, pos2+1):
		position = int(bam_list[t][1])
		read = int(bam_list[t][2])
		index = p.index(position)
		c[index] = read

	p = np.array(p)
	c = np.array(c)

	pos3 = bi_contains(p,start)
	pos4 = bi_contains(p,end)

	global y_limit
	m = max(c)
	if m > y_limit:
		y_limit = m

	if number == 1:
		#caption = ax.fill_between(p,c, color="midnightblue", alpha=0.9)
		caption = ax.fill_between(p,c, color="midnightblue", alpha=0.9, label = group_name)
	else:
		#caption = ax.fill_between(p,c, color="crimson", alpha=0.9)
		caption = ax.fill_between(p,c, color="crimson", alpha=0.9, label = group_name)

	ax.legend(handles = [caption])
	ax.fill_between(p[pos3:pos4+1],c[pos3:pos4+1], color="orange", alpha=0.9)
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.set_xlim(startAll, endAll)
	ax.autoscale(enable = True)
	ax.set_xticklabels([])
	ax.tick_params(axis='both', bottom=False, which='major', labelsize=8)

	return y_limit

def Generate_annotation_plot(ax, isoforms, exonCountList, exonStartList, exonEndList, start, end, startAll, endAll):
	#print(startAll, endAll)
	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)
	print("isoforms are: ",isoforms)

	ystart = 0
	height = 3
	for i in range(isoforms):
		if i>=15:
			print("15 isoforms of this Gene is plotted.")
			break;
		else:
			ax.hlines(y=(ystart+ystart+height)/2, xmin=startAll, xmax=endAll, linewidth=1, color='skyblue', linestyle = '--')
			
			ecount = int(exonCountList[i])
			stList = exonStartList[i]
			enList = exonEndList[i]
			for p in range(ecount):
				ex_s = int(stList[p])
				ex_e = int(enList[p])
				width = int(enList[p]) - int(stList[p]) + 1
				
				rect = patches.Rectangle((ex_s,ystart), width, height, color = 'skyblue', fill = True)
				ax.add_patch(rect)
				
				if (start >= ex_s and end <= ex_e):
					width = end - start + 1
					rect1 = patches.Rectangle((start,ystart), width, height, color = 'orange', alpha=0.9, fill = True)
					ax.add_patch(rect1)
				elif(start <= ex_s and end >= ex_e):
					width = ex_e - ex_s + 1
					rect1 = patches.Rectangle((ex_s,ystart), width, height, color = 'orange', alpha=0.9, fill = True)
					ax.add_patch(rect1)
			
			ystart +=5

	ax.set_xlim(startAll, endAll)
	ax.autoscale(enable = True)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ax.set_yticklabels([])
	ax.tick_params(left=False, axis='both', which='major', labelsize=10)

	return

def Process_user_inputs(region, input1_dir, input2_dir, s1_namelist, s2_namelist, g1_name, g2_name, output_dir, ann_list):
	chrom, geneID, rng = region.split(':')
	start, end = rng.split('-')
	ChromDict = methods.MakeFullDictionary(ann_list, chromosomes)
	GeneDict = ChromDict[chrom]
	exList = GeneDict[geneID]

	mergedExList = methods.MergeIntervals(exList)
	startAll = int(exList[0].st)
	endAll = int(exList[-1].en)

	df = pd.DataFrame(ann_list)
	ann_tt = df.loc[(df[1].str.upper()==geneID) & (df[2]==chrom)]

	exonStartList = {}
	exonEndList = {}
	exonCountList = {}

	isoforms = 0
	for a_row in ann_tt.itertuples():
		exonCount = int(a_row[9])
		exonCountList[isoforms] = exonCount
		exonStartList[isoforms] = a_row[10].split(',')
		exonEndList[isoforms] = a_row[11].split(',')
		isoforms+=1

	title = ""+chrom+":"+start+"-"+end+"("+geneID+")"
	print("Number of isoforms found for gene "+geneID+": "+str(isoforms))
	number_of_subplots = len(s1_namelist)+len(s2_namelist)+1
	fig, axes = plt.subplots(nrows=number_of_subplots, ncols=1)

	print("Generating read coverage plots...")
	for i, ax1 in enumerate(axes[0:len(s1_namelist)]):
		y_limit = Generate_read_coverate_plot(ax1, input1_dir, s1_namelist[i], chrom, geneID, int(start), int(end), startAll, endAll, g1_name, 1)
		if i == 0:
			ax1.set_title(title, color = "black", fontsize = 14)
	for i, ax2 in enumerate(axes[len(s1_namelist):number_of_subplots-1]):
		j = i - len(s1_namelist)
		y_limit = Generate_read_coverate_plot(ax2, input2_dir, s2_namelist[j], chrom, geneID, int(start), int(end), startAll, endAll, g2_name, 2)

	print("Generating annotation plots...")
	ax3 = axes[number_of_subplots-1]
	Generate_annotation_plot(ax3, isoforms, exonCountList, exonStartList, exonEndList, int(start), int(end), startAll, endAll)
	ax3.set_ylabel('Annotation', fontsize = 12)		
	ax3.set_xlabel('Position', fontsize = 12)
	ax3.spines['top'].set_color('none')
	ax3.spines['bottom'].set_color('none')
	ax3.spines['left'].set_color('none')
	ax3.spines['right'].set_color('none')
	ax3.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

	plt.savefig(output_dir+title+'.png')
	plt.savefig(output_dir+title+'.eps', format = 'eps', dpi = 100)
	print(title+" Plotted successfully.")
	y_limit = 0


######### Main program ###########
def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, '*')))]

startTime = time.time()
chromosomes_h = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21', 'chr22','chrX','chrY']
chromosomes_m = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
target_AS = ['SE', 'RI', 'MXE', 'A3SS', 'A5SS']

if(len(sys.argv)<6):
	print("Please provide all mandatory arguments. Example: $ python3 as_quant.py -s mouse -i dir1 dir2")
	sys.exit()

for ii in range(len(sys.argv)):
	if sys.argv[ii] == '-s' or sys.argv[ii] == '-S':
		species = sys.argv[ii+1]
	if sys.argv[ii] == '-o' or sys.argv[ii] == '-O':
			output_dir = sys.argv[ii+1]
	if sys.argv[ii] == '-i' or sys.argv[ii] == '-I':
		input1_dir = sys.argv[ii+1]
		input2_dir = sys.argv[ii+2]

if "-o" not in sys.argv and "-O" not in sys.argv:
	output_dir = 'Output/'

if output_dir[-1] != "/":
	output_dir += "/"
os.makedirs(output_dir, exist_ok=True)

if input1_dir[-1] == "/":
	input1_dir = input1_dir[:-1]
if input2_dir[-1] == "/":
	input2_dir = input2_dir[:-1]

g1_name = input1_dir.split("/")[-1]
g2_name = input2_dir.split("/")[-1]

if species =='human':
	chromosomes = chromosomes_h
	species_folder = 'hg19/'
elif species == 'mouse':
	chromosomes = chromosomes_m
	species_folder = 'mm10/'

s1_namelist = list_dirs(input1_dir)
s2_namelist = list_dirs(input2_dir)

ann_file_reader= open(species_folder+'annotation.csv', "rt")
ann_read = csv.reader(ann_file_reader, delimiter="\t")
ann_list = list(ann_read)

#region = input("Enter the range: (chr:gene:start-end): ")
region = "chr11:ACOX1:116183463-116183624"
Process_user_inputs(region, input1_dir, input2_dir, s1_namelist, s2_namelist, g1_name, g2_name, output_dir, ann_list)
totalTime = time.time() - startTime
print("Time elapsed : ",totalTime)


# python3 as_quant.py -s mouse -i /home/naima/input/mouse_M-_M+/RNA-seq_bam/Minus_M /home/naima/input/mouse_M-_M+/RNA-seq_bam/Plus_M

