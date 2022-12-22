import csv
from operator import itemgetter
import pandas as pd
import time
import bisect
from bisect import bisect_left
import sys, os

class Stack:
	def __init__(self):
		self.items = []

	def size(self):
		return len(self.items)

	def isEmpty(self):
		return self.items == []

	def push(self, val):
		self.items.append(val)

	def top(self):
		if self.isEmpty():
			return None
		else:
			return self.items[self.size()-1]

	def pop(self):
		if self.isEmpty():
			return None
		else:
			return self.items.pop()

def bi_contains(lst, item):
    return bisect_left(lst, item)

def MergeIntervals(inputlist):
	n = len(inputlist)
	inputlist.sort(key = itemgetter(1), reverse = False)

	st = Stack()
	st.push(inputlist[0])

	for i in range(1,n):
		stacktop = st.top()
		if inputlist[i][0] <= stacktop[1]:
			st.pop()
			st_st = stacktop[0]
			st_en = max(stacktop[1],inputlist[i][1])
			st.push((st_st, st_en))
		else:
			st.push(inputlist[i])

	mergedExList = []
	while(True):
		if st.size() == 0:
			break;
		stacktop = st.top()
		mergedExList.append(stacktop)
		st.pop()

	return mergedExList


def CountTotalReadCount(chrom, exList, bam_list, position_row):
	totalCount = 0
	for p in range(len(exList)):
		start = int(exList[p][0])
		end = int(exList[p][1])

		pos1 = bi_contains(position_row, start)
		pos2 = bi_contains(position_row, end)

		if(pos1 < len(position_row) and pos2 < len(position_row)):
			if(int(bam_list[pos2][0]) != end):
				pos2 = pos2 - 1

			for t in range(pos1, pos2+1):
				read = int(bam_list[t][1])
				totalCount += read
		
	return totalCount

def writeResult(chrom, gene, start, end, bam_list, position_row, RC, mergedExListLength, writer_list):
	targetRC = CountTotalReadCount(chrom, [(start, end)], bam_list, position_row)
	targetLength = end - start + 1

	#### Avoiding divide by zero error ##### 
	if targetLength == 0:
		averageTargetRC = 0
	else:
		averageTargetRC = targetRC/targetLength

	if mergedExListLength==targetLength:
		averageRCothers = 0
	else:
		averageRCothers = (RC-targetRC)/(mergedExListLength-targetLength)

	writer_list.append((chrom, gene, start, end, targetRC, targetLength, RC, mergedExListLength, RC-targetRC, mergedExListLength-targetLength,  averageTargetRC, averageRCothers))
	
	return writer_list


def Find_Novel_splicing_events(ChromDict_merged, ChromDict, chromosomes, AS, input_dir, species_folder, sample, output_dir):
	tt = time.time()
	AS_flag = []
	as_df = pd.read_csv(os.path.join(species, AS+'.csv'), delimiter='\t')
	writer_list = []
	output_columns = ['chrom', 'geneName', 'splicedExonStart', 'splicedExonEnd', 'splicedExonReadCount(rc)', 'splicedExonLength(sl)', 'othersExonsReadCount(RC)', 'othersExonsLength(L)', 'RC - rc', 'L - sl', 'splicedExonAverageReadCoverage(n)', 'otherExonsAverageReadCoverage(N)']

	for chrom in chromosomes:
		print("Starting:",chrom)
		tts = time.time()
		GeneDict = ChromDict[chrom]
		GeneDict_merged = ChromDict_merged[chrom]
		txtfile = os.path.join(input_dir, sample, chrom+".txt")
		if os.path.getsize(txtfile) > 0:
			bam_df = pd.read_csv(txtfile, delimiter='\t')
			
			position_row = bam_df.iloc[:, 0].tolist()
			bam_list = bam_df.values.tolist()

			for gene in GeneDict.keys():
				exonList = list(set(GeneDict[gene.upper()]))
				mergedExList = GeneDict_merged[gene.upper()]

				mergedExListLength = 0
				for p in range(len(mergedExList)):
					mergedExListLength += mergedExList[p][1] - mergedExList[p][0] + 1
				RC = CountTotalReadCount(chrom, mergedExList, bam_list, position_row)

				for ex in range(len(exonList)):
					start, end = int(exonList[ex][0]), int(exonList[ex][1])

					if (chrom, gene, start, end) not in AS_flag:
						writer_list = writeResult(chrom, gene, start, end, bam_list, position_row, RC, mergedExListLength, writer_list)
						AS_flag.append((chrom, gene, start, end))

	df_out = pd.DataFrame(writer_list, columns = output_columns)
	df_out.to_csv(os.path.join(output_dir, sample+"_"+AS+".csv"), sep='\t', index=False)

	print("Elapsed time: ",round(((time.time()-tt)/60),2), "minutes")


def Find_splicing_events(ChromDict_merged, chromosomes, AS, input_dir, species, sample, output_dir):
	tt = time.time()
	AS_flag = []
	as_df = pd.read_csv(os.path.join(species, AS+'.csv'), delimiter='\t')
	
	writer_list = []
	output_columns = ['chrom', 'geneName', 'splicedExonStart', 'splicedExonEnd', 'splicedExonReadCount(rc)', 'splicedExonLength(sl)', 'othersExonsReadCount(RC)', 'othersExonsLength(L)', 'RC - rc', 'L - sl', 'splicedExonAverageReadCoverage(n)', 'otherExonsAverageReadCoverage(N)']

	for chrom in chromosomes:
		print("Starting:",chrom)
		tts = time.time()
		GeneDict = ChromDict_merged[chrom]
		txtfile = os.path.join(input_dir, sample, chrom+".txt")
		if os.path.getsize(txtfile) > 0:
			bam_df = pd.read_csv(txtfile, delimiter='\t')
			position_row = bam_df.iloc[:, 0].tolist()
			bam_list = bam_df.values.tolist()

			as_chr_rows = as_df[as_df['chrom']==chrom]
			for ind1, t_row in as_chr_rows.iterrows():
				gene = t_row['gene'].strip().upper()
				mergedExList = GeneDict[gene]

				mergedExListLength = 0
				for p in range(len(mergedExList)):
					mergedExListLength += mergedExList[p][1] - mergedExList[p][0] + 1

				RC = CountTotalReadCount(chrom, mergedExList, bam_list, position_row)

				if AS in ['SE', 'RI']:
					exonStart, exonEnd = t_row['exonStart'], t_row['exonEnd']
					if (chrom, gene, exonStart, exonEnd) not in AS_flag:
						writer_list = writeResult(chrom, gene, exonStart, exonEnd, bam_list, position_row, RC, mergedExListLength, writer_list)
						AS_flag.append((chrom, gene, exonStart, exonEnd))

				elif AS == 'MXE':
					exon1Start, exon1End = t_row['exon1Start'], t_row['exon1End']
					exon2Start, exon2End = t_row['exon2Start'], t_row['exon2End']
					if (chrom, gene, exon1Start, exon1End) not in AS_flag:
						writer_list = writeResult(chrom, gene, exon1Start, exon1End, bam_list, position_row, RC, mergedExListLength, writer_list)
						AS_flag.append((chrom, gene, exon1Start, exon1End))

					if (chrom, gene, exon2Start, exon2End) not in AS_flag:
						writer_list = writeResult(chrom, gene, exon2Start, exon2End, bam_list, position_row, RC, mergedExListLength, writer_list)
						AS_flag.append((chrom, gene, exon2Start, exon2End))
				
				else:
					longExonStart, longExonEnd, shortExonStart, shortExonEnd, strand = t_row['longExonStart'], t_row['longExonEnd'], t_row['shortExonStart'], t_row['shortExonEnd'], t_row['strand']
					start, end = 0, 0
					if AS == 'A5SS':
						if strand == '+':
							start, end = longExonEnd+1, shortExonEnd
						else:
							start, end = shortExonStart, longExonStart-1

					elif AS == 'A3SS':
						if strand == '+':
							start, end = shortExonStart, longExonStart-1
						else:
							start, end = longExonEnd+1, shortExonEnd

					if (chrom, gene, start, end) not in AS_flag:
						writer_list = writeResult(chrom, gene, start, end, bam_list, position_row, RC, mergedExListLength, writer_list)
						AS_flag.append((chrom, gene, start, end))

	df_out = pd.DataFrame(writer_list, columns = output_columns)
	df_out.to_csv(os.path.join(output_dir, sample+"_"+AS+".csv"), sep='\t', index=False)

	print("Elapsed time: ",round(((time.time()-tt)/60),2), "minutes")


def MakeFullDictionary(ann_df, chromosomes):
	ChromDict = {}
	for chrom in chromosomes:
		GeneDict = {}
		chr_rows = ann_df[ann_df['chrom']==chrom]
		gene_list = list(set(chr_rows['gene']))
		for gene in gene_list:
			gene_rows = chr_rows[chr_rows['gene']==gene]
			exList = []
			for index, row in gene_rows.iterrows():
				exonCount = row['exonCount']
				exonStarts = list(filter(None, row['exonStarts'].split(',')))
				exonEnds = list(filter(None, row['exonEnds'].split(',')))
				for i in range(exonCount):
					st, en = int(exonStarts[i]), int(exonEnds[i])
					if (st, en) not in exList:
						exList.append((st, en))

			GeneDict[gene.strip().upper()] = exList

		ChromDict[chrom] = GeneDict

	return ChromDict

def merge_ChromDict(ChromDict, chromosomes):
	ChromDict_merged = {}
	for chrom in chromosomes:
		GeneDict_merged = {}
		GeneDict = ChromDict[chrom]
		for gene in GeneDict.keys():
			exonList = GeneDict[gene.upper()]
			mergedExonList = MergeIntervals(exonList)
			GeneDict_merged[gene.upper()] = mergedExonList
		ChromDict_merged[chrom] = GeneDict_merged

	return ChromDict_merged