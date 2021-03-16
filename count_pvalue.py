import csv
from operator import attrgetter
import pandas as pd
import time
from scipy.stats import chisquare
import sys
import os
import xlsxwriter

def Count_pvalue(AS, output_dir, s1_namelist, s2_namelist, g1_name, g2_name):
	len_S1 = len(s1_namelist)
	len_S2 = len(s2_namelist)

	S1_list = {}
	S2_list = {}
	with open(output_dir+AS+"_Output.csv",'w') as f:
		writer = csv.writer(f, dialect='excel',delimiter='\t')
		writer.writerow(['Chrom', 'Gene Name', 'Exon Start', 'Exon End', 'p-value', 'Ratio difference', 'Absolute Ratio difference', g1_name+':n1', g2_name+':n2',  g1_name+':N1', g2_name+':N2', 'Chrom region Long'])
		for k, sample1 in enumerate(s1_namelist):
			S1_reader = open(output_dir+sample1+"_"+AS+".csv", "rt", encoding='ascii')
			S1_read = csv.reader(S1_reader, delimiter="\t")
			S1_list[k] = list(S1_read)

		for p, sample2 in enumerate(s2_namelist):
			S2_reader = open(output_dir+sample2+"_"+AS+".csv", "rt", encoding='ascii')
			S2_read = csv.reader(S2_reader, delimiter="\t")
			S2_list[p] = list(S2_read)
			AS_file_length = len(S2_list[p])
		
		#print("AS- Length: ", AS_file_length)

		for i in range(1,AS_file_length):
			S1_n_total, S1_N_total, S2_n_total, S2_N_total = 0.0, 0.0, 0.0, 0.0
			for k in range(len_S1):
				full_list = S1_list[k]
				S1_n_total += float(full_list[i][10])
				S1_N_total += float(full_list[i][11])

				chrom = full_list[i][0]
				gene = full_list[i][1]
				start = full_list[i][2]
				end =  full_list[i][3]

				#print("S1 Printing: ", gene, start, end, float(full_list[i][10]), float(full_list[i][11]), S1_n_total, S1_N_total)

			n1 = S1_n_total/len_S1
			N1 = S1_N_total/len_S1

			for p in range(len_S2):
				full_list = S2_list[p]
				S2_n_total += float(full_list[i][10])
				S2_N_total += float(full_list[i][11])

				#print("S2 Printing: ", gene, start, end, float(full_list[i][10]), float(full_list[i][11]), S2_n_total, S2_N_total)

			n2 = S2_n_total/len_S2
			N2 = S2_N_total/len_S2

			N1 = N1 + n1
			N2 = N2 + n2
			if N1>0 and N2>0:
				ratio_diff = (n1/N1) - (n2/N2)
				abs_ratio_diff = abs(ratio_diff)
				#print("ratio_diff, abs_ratio_diff", ratio_diff, abs_ratio_diff)
			
				P0 = (n1+n2)/(N1+N2)
				n10 = N1 * P0
				n20 = N2 * P0
				exp=[n10, N1-n10, n20, N2-n20]
				if 0 not in exp:
					res = chisquare([n1, N1-n1, n2, N2-n2], f_exp=exp, ddof = 1)
					#print("p-value:",res[1])
					chrom_region_long = chrom+':'+gene+":"+start+"-"+end
					writer.writerow([chrom, gene, start, end, res[1], ratio_diff, abs_ratio_diff, n1, n2, N1-n1, N2-n2, chrom_region_long])
		
		f.close()
	
	for sample1 in s1_namelist:
		os.remove(output_dir+sample1+"_"+AS+".csv")
	for sample2 in s2_namelist:
		os.remove(output_dir+sample2+"_"+AS+".csv")
	
	return