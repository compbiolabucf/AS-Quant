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
#target_AS = ['MXE']

config = configparser.ConfigParser()
config.read('configuration.ini')
#print(config.sections())

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
#print(species, input1_dir, input2_dir, output_dir)

g1_name = input1_dir.split("/")[-1]
g2_name = input2_dir.split("/")[-1]

if species =='human':
	chromosomes = chromosomes_h
	species_folder = 'hg19/'
elif species == 'mouse':
	chromosomes = chromosomes_m
	species_folder = 'mm10/'


#### Listing out the bamfile's names

current = os.getcwd()
os.chdir(input1_dir)
for file1 in glob.glob("*.bam"):
    preprocess.SamtoText(input1_dir, file1, chromosomes)
os.chdir(input2_dir)
for file2 in glob.glob("*.bam"):
    preprocess.SamtoText(input2_dir, file2, chromosomes)
os.chdir(current)
#print("################", s1_namelist, s2_namelist, "#########", g1_name, g2_name)


######## Listing out the folders
s1_namelist = list_dirs(input1_dir)
s2_namelist = list_dirs(input2_dir)
print(s1_namelist, s2_namelist, g1_name, g2_name)

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

"""
###################### For printing ChromDict to Debug #######################################
for chrom in ChromDict.keys():
	#print(chrom)
	if chrom == 'chr5':
		geneDict = ChromDict[chrom]
		#print(len(geneDict))
		for gene in geneDict.keys():
			print(gene)
			if gene == 'Wdr66':
				exlist = geneDict[gene]
				print(len(exlist))
				print("#########################")
				sys.exit()
###################### For printing ChromDict to Debug #######################################
"""
time_s = time.time()
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
	count_pvalue.Count_pvalue(AS, output_dir, s1_namelist, s2_namelist, g1_name, g2_name)
	print(AS, "counting p-val Complete")

	df = pd.read_csv(output_dir+AS+"_Output.csv", delimiter = '\t')
	df.to_excel(writer_out, sheet_name=AS, index = None, header=True)
	os.remove(output_dir+AS+"_Output.csv")
writer_out.save()


totalTime = time.time() - startTime
method_time = time.time() - time_s
print("Total AS-Quant time is : ",round((totalTime/60),2), "minutes")
print("method time: ", method_time)