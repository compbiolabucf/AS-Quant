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
import os

y_limit = 0

def bi_contains(lst, item):
    return bisect_left(lst, item)

def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, '*')))]

def Generate_read_coverate_plot(ax, pathin, sample, chrom, geneID, start, end, startAll, endAll, number):
	bam_file_reader= open(pathin+'/'+chrom+".txt", "rt")
	bam_read = csv.reader(bam_file_reader, delimiter="\t")
	bam_list = list(bam_read)
	position_row = [int(bam_list[i][1]) for i in range(len(bam_list))]

	ax = ax or plt.gca()
	x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
	x_formatter.set_scientific(False)
	ax.xaxis.set_major_formatter(x_formatter)
	ax.tick_params(axis='both', which='major', labelsize=6)

	pos1 = bi_contains(position_row, startAll)
	pos2 = bi_contains(position_row, endAll)
	if(int(bam_list[pos2][1]) != endAll):
		pos2 = pos2 - 1

	p = []
	c = []
	read = 0
	length = endAll - startAll + 1
	#print(length)
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

	# label = 'Sample'+str(number)
	caption = ax.fill_between(p,c, color="skyblue", alpha=0.9, label = sample)
	ax.fill_between(p[pos3:pos4+1],c[pos3:pos4+1], color="crimson", alpha=0.9)
	ax.legend(handles = [caption])
	ax.set_xlim(startAll, endAll)
	ax.autoscale(enable = True)

	return y_limit

def Generate_annotation_plot(ax, isoforms, exonCountList, exonStartList, exonEndList, start, end, startAll, endAll):
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
			ax.hlines(y=(ystart+ystart+height)/2, xmin=startAll, xmax=endAll, linewidth=1, color='slateblue', linestyle = '--')
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
					rect1 = patches.Rectangle((start,ystart), width, height, color = 'crimson', alpha=0.9, fill = True)
					ax.add_patch(rect1)
				elif(start <= ex_s and end >= ex_e):
					width = ex_e - ex_s + 1
					rect1 = patches.Rectangle((ex_s,ystart), width, height, color = 'crimson', alpha=0.9, fill = True)
					ax.add_patch(rect1)

			ystart +=5

	ax.set_xlim(startAll, endAll)
	ax.autoscale(enable = True)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["left"].set_visible(False)
	ax.set_yticklabels([])
	ax.tick_params(left=False, axis='both', which='major', labelsize=6)

	return

def Take_user_inputs(samplenames, input1_dir, input2_dir, s1_namelist, s2_namelist, g1_name, g2_name, output_dir, ann_list):
	region = input("Enter the range: ")
	chrom, geneID, rng = region.split(':')
	start, end = rng.split('-')

	ChromDict = methods.MakeFullDictionary(ann_list, chromosomes)
	GeneDict = ChromDict[chrom]
	exList = GeneDict[geneID]

	mergedExList = methods.MergeIntervals(exList)
	startAll = int(exList[0].st)
	endAll = int(exList[-1].en)

	df = pd.DataFrame(ann_list)
	ann_tt = df.loc[df[2]==chrom]

	exonStartList = {}
	exonEndList = {}
	exonCountList = {}

	isoforms = 0
	for a_row in ann_tt.itertuples():
		geneName = a_row[2].strip()
		if geneName == geneID:
			exonCount = int(a_row[9])
			exonCountList[isoforms] = exonCount
			exonStartList[isoforms] = a_row[10].split(',')
			exonEndList[isoforms] = a_row[11].split(',')
			isoforms+=1

	y_axis_height = 20
	title = ""+chrom+":"+start+"-"+end+"("+geneID+")"

	fig = plt.figure()

	ax1 = fig.add_subplot(3,1,1)

	subplots_adjust(hspace=0.000)
	number_of_subplots=len(s1_namelist)+len(s2_namelist)

	for i,v in enumerate(xrange(number_of_subplots)):
	    v = v+1
	    ax1 = subplot(number_of_subplots,1,v)
	    ax1.plot(x,y)
		ax1.set_title(title, color = "black")
		ax1.set_ylabel('Counts')
		y_limit = Generate_read_coverate_plot(ax1, input1_dir, samplenames[0], chrom, geneID, int(start), int(end), startAll, endAll, 1)

	ax2 = fig.add_subplot(3,1,2)
	ax2.set_xlabel('Position')
	ax2.set_ylabel('Counts')
	y_limit = Generate_read_coverate_plot(ax2, input2_dir, samplenames[1], chrom, geneID, int(start), int(end), startAll, endAll, 2)

	# set ylim for both axes after getting max ylim values
	ax1.set_ylim(0, y_limit+y_limit*0.2)
	ax2.set_ylim(0, y_limit+y_limit*0.2)
	
	ax3 = fig.add_subplot(3,1,3)
	ax3.set_ylabel('Isoforms')
	Generate_annotation_plot(ax3, isoforms, exonCountList, exonStartList, exonEndList, int(start), int(end), startAll, endAll)

	os.makedirs(output_dir, exist_ok=True)
	plt.savefig(output_dir+title+'.png')
	plt.savefig(output_dir+title+'.eps', format = 'eps', dpi = 1000)
	print("Plotted successfully.")
	y_limit = 0



######### Main starts here #################
startTime = time.time()

chromosomes_h = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21', 'chr22','chrX','chrY']
chromosomes_m = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']

if(len(sys.argv)<5):
	print("Please provide all of the mandatory arguments. Example: $ python3 make_plots.py -s human s1 s2")
	sys.exit()

if len(sys.argv)==7:
	for ii in range(len(sys.argv)):
		if sys.argv[ii] == '-s' or sys.argv[ii] == '-S':
			species = sys.argv[ii+1]
		if sys.argv[ii] == '-o' or sys.argv[ii] == '-O':
			output_dir = sys.argv[ii+1]+'/'

elif len(sys.argv)==5:
	species = sys.argv[2]
	output_dir = 'Plots/'
os.makedirs(output_dir, exist_ok=True)

input1_dir = sys.argv[-2]
input2_dir = sys.argv[-1]
g1_name = input1_dir.split("/")[-1]
g2_name = input2_dir.split("/")[-1]

s1_namelist = list_dirs(input1_dir)
s2_namelist = list_dirs(input2_dir)
print(s1_namelist, s2_namelist)

if species =='human':
	chromosomes = chromosomes_h
	inp = 'hg19/'
elif species == 'mouse':
	chromosomes = chromosomes_m
	inp = 'mm10/'

ann_file_reader= open(inp+'annotation.csv', "rt")
ann_read = csv.reader(ann_file_reader, delimiter="\t")
ann_list = list(ann_read)

Take_user_inputs(samplenames, input1_dir, input2_dir, s1_namelist, s2_namelist, g1_name, g2_name, output_dir, ann_list)
totalTime = time.time() - startTime
print("Total program time is : ",totalTime)


