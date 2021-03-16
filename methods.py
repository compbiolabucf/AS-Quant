import csv
from operator import attrgetter
import pandas as pd
import time
import bisect
from bisect import bisect_left
import sys

class EXON:
	def __init__(self):
		self.st = 0
		self.en = 0

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


def InsertIntoOldChromDict(ChromDict, chrom, exListNew, geneName):
	GeneDict = ChromDict[chrom]
	if geneName not in GeneDict.keys():
		GeneDict[geneName] = exListNew
	else:
		exList = GeneDict[geneName]
		exList.extend(exListNew)
		GeneDict[geneName] = exList
	return GeneDict


def bi_contains(lst, item):
    return bisect_left(lst, item)

def MergeIntervals(inputlist):
	n = len(inputlist)
	inputlist.sort(key = attrgetter('st'), reverse = False)

	st = Stack()
	st.push(inputlist[0])

	for i in range(1,n):
		stacktop = st.top()
		if inputlist[i].st <= stacktop.en:
			st.pop()
			stacktop.en = max(stacktop.en,inputlist[i].en)
			st.push(stacktop)
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
		start = int(exList[p].st)
		end = int(exList[p].en)

		pos1 = bi_contains(position_row, start)
		pos2 = bi_contains(position_row, end)

		if(pos1 < len(bam_list) and pos2 < len(bam_list)):
			if(int(bam_list[pos2][1]) != end):
				pos2 = pos2 - 1

			for t in range(pos1, pos2+1):
				read = int(bam_list[t][2])
				totalCount += read
		
	return totalCount

def writeResult(chrom, geneID, start, end, newList, bam_list, position_row, RC, mergedExListLength, writer):
	#tt = time.time()
	targetRC = CountTotalReadCount(chrom, newList, bam_list, position_row)
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

	writer.writerow([chrom, geneID, start, end, targetRC, targetLength, RC, mergedExListLength, RC-targetRC, mergedExListLength-targetLength,  averageTargetRC, averageRCothers])
	#print("Write result exact time: ", time.time() - tt)

def callA3SS(t_row):
	newList = []
	e = EXON()
	strand = t_row[7]
	if strand == '+':
		e.st = int(t_row[5])
		e.en = int(t_row[3])-1
	else:
		e.st = int(t_row[4])+1
		e.en = int(t_row[6])
	newList.append(e)
	return e.st, e.en, newList

def callA5SS(t_row):
	newList = []
	e = EXON()
	strand = t_row[7]
	if strand == '+':
		e.st = int(t_row[4])+1
		e.en = int(t_row[6])
	else:
		e.st = int(t_row[5])-1
		e.en = int(t_row[4])
	newList.append(e)
	return e.st, e.en, newList

def callMXE(t_row):
	newList1 = []
	e1 = EXON()
	e1.st = int(t_row[3])
	e1.en = int(t_row[4])
	newList1.append(e1)
	#writeResult(chrom, geneID, e1.st, e1.en, newList, bam_list, position_row, RC, mergedExListLength, writer)

	newList2 = []
	e2 = EXON()
	e2.st = int(t_row[5])
	e2.en = int(t_row[6])
	newList2.append(e2)
	#writeResult(chrom, geneID, e2.st, e2.en, newList, bam_list, position_row, RC, mergedExListLength, writer)
	return e1.st, e1.en, newList1, e2.st, e2.en, newList2

def callSE_RI(t_row):
	newList = []
	e = EXON()
	e.st = int(t_row[3])
	e.en = int(t_row[4])
	newList.append(e)
	return e.st, e.en, newList

def Generate(ChromDict, chromosomes, AS, input_dir, species_folder, sample, output_dir):
	tt = time.time()
	AS_flag = []
	target_file_reader= open(species_folder+AS+'.csv', "rt")
	target_read = csv.reader(target_file_reader, delimiter="\t")
	target_list = list(target_read)
	df = pd.DataFrame(target_list)
	
	with open(output_dir+sample+"_"+AS+".csv",'w') as f:
		writer = csv.writer(f, delimiter='\t')
		writer.writerow(['Chrom', 'Gene Name', 'Exon Start', 'Exon End', 'Target Read Count', 'Target Length', 'Others: Read Count', 'Others: Length', 'Others: Read Count- target RC', 'Others: Length - exon length', 'Average Target Read Count(n)', 'Average Read Count All(N)'])

		for chrom in chromosomes:
			print("Starting:",chrom)
			tts = time.time()
			bam_file_reader= open(input_dir+'/'+sample+'/'+chrom+".txt", "rt")
			bam_read = csv.reader(bam_file_reader, delimiter="\t")
			bam_list = list(bam_read)
			position_row = [int(bam_list[i][1]) for i in range(len(bam_list))]
			target_tt = df.loc[df[0]==chrom]
			
			for t_row in target_tt.itertuples():
				geneID = t_row[2].strip().upper()
				GeneDict = ChromDict[chrom]
				mergedExList = GeneDict[geneID]

				mergedExListLength = 0
				for p in range(len(mergedExList)):
					mergedExListLength += mergedExList[p].en - mergedExList[p].st + 1

				RC = CountTotalReadCount(chrom, mergedExList, bam_list, position_row)
				if(AS == 'SE' or AS == 'RI'):
					start, end, newList = callSE_RI(t_row)
					if (chrom, geneID, start, end) not in AS_flag:
						writeResult(chrom, geneID, start, end, newList, bam_list, position_row, RC, mergedExListLength, writer)
						AS_flag.append((chrom, geneID, start, end))

				elif(AS == 'MXE'):
					start1, end1, newList1, start2, end2, newList2 = callMXE(t_row)
					if (chrom, geneID, start1, end1) not in AS_flag:
						writeResult(chrom, geneID, start1, end1, newList1, bam_list, position_row, RC, mergedExListLength, writer)
						AS_flag.append((chrom, geneID, start1, end1))
					if (chrom, geneID, start2, end2) not in AS_flag:
						writeResult(chrom, geneID, start2, end2, newList2, bam_list, position_row, RC, mergedExListLength, writer)
						AS_flag.append((chrom, geneID, start2, end2))
				
				elif(AS == 'A5SS'):
					start, end, newList = callA5SS(t_row)
					if (chrom, geneID, start, end) not in AS_flag:
						writeResult(chrom, geneID, start, end, newList, bam_list, position_row, RC, mergedExListLength, writer)
						AS_flag.append((chrom, geneID, start, end))

				elif(AS == 'A3SS'):
					start, end, newList = callA3SS(t_row)
					if (chrom, geneID, start, end) not in AS_flag:
						writeResult(chrom, geneID, start, end, newList, bam_list, position_row, RC, mergedExListLength, writer)
						AS_flag.append((chrom, geneID, start, end))
			print("Time for chrom: ", chrom, "is: ", time.time() - tts)

		f.close()
	print("Elapsed time: ",round(((time.time()-tt)/60),2), "minutes")


def MakeFullDictionary(ann_list, chromosomes):
	ChromDict = {}
	
	for a_row in ann_list:
		chrom = a_row[2]
		if chrom in chromosomes:
			geneName = a_row[1].strip().upper()
			exonCount = int(a_row[8])
			exonStartList = a_row[9].split(',')
			exonEndList = a_row[10].split(',')

			exList = []
			for i in range(exonCount):
				exonStart = int(exonStartList[i])
				exonEnd = int(exonEndList[i])
				newExon = EXON()
				newExon.st = exonStart
				newExon.en = exonEnd
				exList.append(newExon)
			
			if chrom not in ChromDict.keys():
				GeneDict = {}
				GeneDict[geneName] = exList
			else:
				GeneDict = InsertIntoOldChromDict(ChromDict, chrom, exList, geneName)

			ChromDict[chrom] = GeneDict


	return ChromDict