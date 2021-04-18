import matplotlib
matplotlib.use("agg")
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
import glob, os
import configparser
import xlsxwriter
from matplotlib import rcParams
plt.rc('legend',**{'fontsize':14})

rcParams.update({
    'font.family':'arial',
    })

y_limit = 0

def bi_contains(lst, item):
    return bisect_left(lst, item)

def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, '*')))]

def Generate_read_coverate_plot(ax, pathin, sample, chrom, geneID, start, end, startAll, endAll, number):
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
		caption = ax.fill_between(p,c, color="midnightblue", alpha=0.9)
	else:
		caption = ax.fill_between(p,c, color="crimson", alpha=0.9)

	ax.legend(handles = [caption])
	ax.fill_between(p[pos3:pos4+1],c[pos3:pos4+1], color="orange", alpha=0.9)
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.set_xlim(startAll, endAll)
	ax.autoscale(enable = True)
	ax.set_xticklabels([])
	ax.tick_params(axis='both', bottom=False, which='major', labelsize=10)

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

def Take_user_inputs(chrom, geneID, start, end, input1_dir, input2_dir, s1_namelist, s2_namelist, output_dir, ann_list):
	ChromDict = methods.MakeFullDictionary(ann_list, chromosomes)
	GeneDict = ChromDict[chrom]
	exList = GeneDict[geneID]

	mergedExList = methods.MergeIntervals(exList)
	startAll = int(exList[0].st)
	endAll = int(exList[-1].en)

	df = pd.DataFrame(ann_list)
	ann_tt = df.loc[(df[1]==geneID) & (df[2]==chrom)]

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

	number_of_subplots = len(s1_namelist)+len(s2_namelist)+1
	fig, axes = plt.subplots(nrows=number_of_subplots, ncols=1, figsize=(12,8))

	for i, ax1 in enumerate(axes[0:len(s1_namelist)]):
		y_limit = Generate_read_coverate_plot(ax1, input1_dir, s1_namelist[i], chrom, geneID, int(start), int(end), startAll, endAll, 1)

	for i, ax2 in enumerate(axes[len(s1_namelist):number_of_subplots-1]):
		y_limit = Generate_read_coverate_plot(ax2, input2_dir, s2_namelist[i], chrom, geneID, int(start), int(end), startAll, endAll, 2)

	ax3 = axes[number_of_subplots-1]
	Generate_annotation_plot(ax3, isoforms, exonCountList, exonStartList, exonEndList, int(start), int(end), startAll, endAll)
	ax3.set_ylabel('Counts')		
	ax3.spines['top'].set_color('none')
	ax3.spines['bottom'].set_color('none')
	ax3.spines['left'].set_color('none')
	ax3.spines['right'].set_color('none')
	ax3.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

	"""
    ax.set_title(title, color = "black")
	
	ax.set_title(title, color = "black", fontsize = 20)
	ax.set_ylabel('Read Coverage', fontsize = 16)
	"""

	"""
	ax1 = fig.add_subplot(7,1,1)
	y_limit = Generate_read_coverate_plot(ax1, input1_dir, s1_namelist[0], chrom, geneID, int(start), int(end), startAll, endAll, labels, 1)
	ax2 = fig.add_subplot(7,1,2)
	y_limit = Generate_read_coverate_plot(ax2, input1_dir, s1_namelist[1], chrom, geneID, int(start), int(end), startAll, endAll, labels, 1)
	ax3 = fig.add_subplot(7,1,3)
	y_limit = Generate_read_coverate_plot(ax3, input1_dir, s1_namelist[2], chrom, geneID, int(start), int(end), startAll, endAll, labels, 1)
	ax4 = fig.add_subplot(7,1,4)
	y_limit = Generate_read_coverate_plot(ax4, input2_dir, s2_namelist[0], chrom, geneID, int(start), int(end), startAll, endAll, labels, 2)
	ax5 = fig.add_subplot(7,1,5)
	y_limit = Generate_read_coverate_plot(ax5, input2_dir, s2_namelist[1], chrom, geneID, int(start), int(end), startAll, endAll, labels, 2)
	ax6 = fig.add_subplot(7,1,6)
	y_limit = Generate_read_coverate_plot(ax6, input2_dir, s2_namelist[2], chrom, geneID, int(start), int(end), startAll, endAll, labels, 2)
	ax7 = fig.add_subplot(7,1,7)
	ax7.set_ylabel('Annotation', fontsize="16")
	Generate_annotation_plot(ax7, isoforms, exonCountList, exonStartList, exonEndList, int(start), int(end), startAll, endAll)
	"""
	"""
	ax1 = fig.add_subplot(3,1,1)
	y_limit = Generate_read_coverate_plot(ax1, input1_dir, s1_namelist[0], chrom, geneID, int(start), int(end), startAll, endAll, labels, 1)
	ax2 = fig.add_subplot(3,1,2)
	y_limit = Generate_read_coverate_plot(ax2, input2_dir, s2_namelist[0], chrom, geneID, int(start), int(end), startAll, endAll, labels, 2)
	ax3 = fig.add_subplot(3,1,3)
	ax3.set_ylabel('Annotation', fontsize="16")
	Generate_annotation_plot(ax3, isoforms, exonCountList, exonStartList, exonEndList, int(start), int(end), startAll, endAll)
	"""
	title = title+"_"+start
	os.makedirs(output_dir, exist_ok=True)
	plt.savefig(output_dir+title+'.png')
	plt.savefig(output_dir+title+'.eps', format = 'eps', dpi = 1000)
	print("Plotted successfully.")
	y_limit = 0



######### Main starts here #################
startTime = time.time()

chromosomes_h = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21', 'chr22','chrX','chrY']
chromosomes_m = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']

config = configparser.ConfigParser()
config.read('configuration_plot.ini')
species = config['SPECIES']['species']
input1_dir = config['INPUT_SAMPLES']['input1']
input2_dir = config['INPUT_SAMPLES']['input2']
if input1_dir[-1] == "/":
	input1_dir = input1_dir[:-1]
if input2_dir[-1] == "/":
	input2_dir = input2_dir[:-1]

output_dir = config['OUTPUT_FOLDER']['output_dir']
if output_dir[-1] != "/":
	output_dir += "/"
os.makedirs(output_dir, exist_ok=True)

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
chrom, geneID, rng = region.split(':')
start, end = rng.split('-')
print(chrom, geneID, start, end)
Take_user_inputs(chrom, geneID, start, end, input1_dir, input2_dir, s1_namelist, s2_namelist, output_dir, ann_list)

totalTime = time.time() - startTime
print("Total program time is : ",totalTime)


# "chr11:ACOX1:116183463-116183624"