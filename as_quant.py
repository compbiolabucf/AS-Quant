import csv
import glob, os
from operator import attrgetter
import pandas as pd
import time
import bisect
import count_pvalue
import methods
import preprocess
import sys
import configparser

def list_dirs(path):
    return [os.path.basename(x) for x in filter(os.path.isdir, glob.glob(os.path.join(path, '*')))]

startTime = time.time()
chromosomes_h = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21', 'chr22','chrX','chrY']
chromosomes_m = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
target_AS = ['SE', 'RI', 'MXE', 'A3SS', 'A5SS']


if(len(sys.argv)<6):
	print("Please provide all mandatory arguments. Example: $ python3 as_quant.py -s human -i dir1 dir2")
	sys.exit()

for ii in range(len(sys.argv)):
	if sys.argv[ii] == '-s' or sys.argv[ii] == '-S':
		species = sys.argv[ii+1]
	if sys.argv[ii] == '-o' or sys.argv[ii] == '-O':
			output_dir = sys.argv[ii+1]
	if sys.argv[ii] == '-i' or sys.argv[ii] == '-I':
		input1_dir = sys.argv[ii+1]
		input2_dir = sys.argv[ii+2]
	if sys.argv[ii] == '-method':
		method = sys.argv[ii+1]

if "-o" not in sys.argv and "-O" not in sys.argv:
	output_dir = 'Output/'

if output_dir[-1] != "/":
	output_dir += "/"
os.makedirs(output_dir, exist_ok=True)

if input1_dir[-1] == "/":
	input1_dir = input1_dir[:-1]
if input2_dir[-1] == "/":
	input2_dir = input2_dir[:-1]

if '-novel' in sys.argv:
	novel = 'yes'
else:
	novel = 'no'

if method == 'ranksum':
	count1 = len(glob.glob1(input1_dir,"*.bam"))
	count2 = len(glob.glob1(input2_dir,"*.bam"))
	if count1 < 2 or count2 < 2:
		print("Please provide multiple samples/replicates in each group to run ranksum test, otherwise select chisquare.")
		sys.exit()

if '-method' not in sys.argv:
	method = 'chisquare'

g1_name = input1_dir.split("/")[-1]
g2_name = input2_dir.split("/")[-1]

if species =='human':
	chromosomes = chromosomes_h
	species_folder = 'hg19/'
elif species == 'mouse':
	chromosomes = chromosomes_m
	species_folder = 'mm10/'

print("Generating read coverage files for each chromosome...")
current = os.getcwd()
os.chdir(input1_dir)
for file1 in glob.glob("*.bam"):
    preprocess.SamtoText(input1_dir, file1, chromosomes)
os.chdir(input2_dir)
for file2 in glob.glob("*.bam"):
    preprocess.SamtoText(input2_dir, file2, chromosomes)
os.chdir(current)

s1_namelist = list_dirs(input1_dir)
s2_namelist = list_dirs(input2_dir)

ann_file_reader= open(species_folder+'annotation.csv', "rt")
ann_read = csv.reader(ann_file_reader, delimiter="\t")
ann_list = list(ann_read)

ChromDict = methods.MakeFullDictionary(ann_list, chromosomes)
for chrom in ChromDict.keys():
	geneDict = ChromDict[chrom]
	for gene in geneDict.keys():
		exonList = geneDict[gene]
		mergedExonList = methods.MergeIntervals(exonList)
		geneDict[gene] = mergedExonList
	ChromDict[chrom] = geneDict

if novel.upper() == 'YES':
	print("Running AS-Quant for detecting significant novel (unannotated) splicing events...")
	target_AS = ['All']
	ChromDict_exon = methods.MakeFullDictionary(ann_list, chromosomes)
	for AS in target_AS:
		for sample in s1_namelist:
			print("Executing: ",sample, "in group 1")
			methods.Generate_Novel(ChromDict, ChromDict_exon, chromosomes[0:5], AS, input1_dir, species_folder, sample, output_dir)
		for sample in s2_namelist:
			print("Executing: ",sample, "in group 2")
			methods.Generate_Novel(ChromDict, ChromDict_exon, chromosomes[0:5], AS, input2_dir, species_folder, sample, output_dir)
			
else:
	print("Running AS-Quant for detecting significant annotated splicing events...")
	for AS in target_AS:
		print("Spliced Exon type: ",AS)
		for sample in s1_namelist:
			print("Executing: ",sample, "in group 1")
			methods.Generate(ChromDict, chromosomes, AS, input1_dir, species_folder, sample, output_dir)

		for sample in s2_namelist:
			print("Executing: ",sample, "in group 2")
			methods.Generate(ChromDict, chromosomes, AS, input2_dir, species_folder, sample, output_dir)

writer_out = pd.ExcelWriter(output_dir+g1_name+"_Vs_"+g2_name+".xlsx", engine='xlsxwriter')
for AS in target_AS:
	if method == 'ranksum':
		count_pvalue.Count_pvalue_replicates(AS, output_dir, s1_namelist, s2_namelist, g1_name, g2_name)
	elif method == 'chisquare':
		count_pvalue.Count_pvalue(AS, output_dir, s1_namelist, s2_namelist, g1_name, g2_name)
	print(AS, "counting p-val Complete")

	df = pd.read_csv(output_dir+AS+"_Output.csv", delimiter = '\t')
	df.to_excel(writer_out, sheet_name=AS, index = None, header=True)
	os.remove(output_dir+AS+"_Output.csv")
writer_out.save()


totalTime = time.time() - startTime
print("Total AS-Quant time is : ",round((totalTime/60),2), "minutes")