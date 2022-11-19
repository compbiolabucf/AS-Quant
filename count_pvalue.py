import csv
from operator import attrgetter
import pandas as pd
import time
from scipy.stats import chisquare, ranksums
import sys
import os
import xlsxwriter
import numpy as np
from scipy import stats

def Count_pvalue_replicates(AS, output_dir, s1_namelist, s2_namelist, g1_name, g2_name):
	len_S1 = len(s1_namelist)
	len_S2 = len(s2_namelist)

	S1_list, S2_list = {}, {}
	writer_list = []
	output_columns = ['chrom', 'geneName', 'splicedExonStart', 'splicedExonEnd', 'P_value', 'ChromRegionLong']

	for k, sample1 in enumerate(s1_namelist):
		S1_df = pd.read_csv(os.path.join(output_dir, sample1+"_"+AS+".csv"), delimiter="\t")
		S1_list[k] = list(S1_df.values.tolist())

	for p, sample2 in enumerate(s2_namelist):
		S2_df = pd.read_csv(os.path.join(output_dir, sample2+"_"+AS+".csv"), delimiter="\t")
		S2_list[p] = list(S2_df.values.tolist())
		AS_file_length = len(S2_df)

	#### need to check this part, whether it can read the values
	for i in range(1, AS_file_length):
		X, Y = [], []
		N1s, N2s = [], []
		for k in range(len_S1):
			full_list = S1_list[k]
			chrom = full_list[i][0]
			gene = full_list[i][1]
			start = full_list[i][2]
			end =  full_list[i][3]

			n1 = float(full_list[i][10])
			N1 = float(full_list[i][11])
			N1s.append(N1)
			N1 = N1 + n1
			if N1 > 0:
				r1 = n1/N1
				X.append(r1)

		for p in range(len_S2):
			full_list = S2_list[p]
			n2 = float(full_list[i][10])
			N2 = float(full_list[i][11])
			N2s.append(N2)
			N1 = N2 + n2
			if N2 > 0:
				r2 = n2/N2
				Y.append(r2)

		##### added N filter: 11.14.2022
		##### both N should be >= 10
		if (min(np.mean(N1s), np.mean(N2s)) >= 10) and sum(np.array(X))!=0 and sum(np.array(Y))!=0:
			#stat, p_val = ranksums(X, Y)
			stat, p_val = stats.ttest_ind(X, Y)
			#print("P-val", p_val)
			chrom_region_long = chrom+':'+gene+":"+str(start)+"-"+str(end)
			writer_list.append((chrom, gene, start, end, p_val, chrom_region_long))
		
	
	df_out = pd.DataFrame(writer_list, columns = output_columns)
	df_out.to_csv(os.path.join(output_dir, AS+"_Output.csv"), sep='\t', index=False)
	
	for sample1 in s1_namelist:
		os.remove(os.path.join(output_dir, sample1+"_"+AS+".csv"))
	for sample2 in s2_namelist:
		os.remove(os.path.join(output_dir, sample2+"_"+AS+".csv"))

	return

def Count_pvalue(AS, output_dir, s1_namelist, s2_namelist, g1_name, g2_name):
	len_S1, len_S2 = len(s1_namelist), len(s2_namelist)
	S1_list, S2_list = {}, {}

	writer_list = []
	output_columns = ['chrom', 'geneName', 'splicedExonStart', 'splicedExonEnd', 'P_value', 'ratioDifference', 'absoluteRatioDifference', g1_name+':n1', g1_name+':N1', g2_name+':n2', g2_name+':N2', 'chromRegionLLong']
	
	for k, sample1 in enumerate(s1_namelist):
		S1_df = pd.read_csv(os.path.join(output_dir, sample1+"_"+AS+".csv"), delimiter="\t")
		S1_list[k] = list(S1_df.values.tolist())

	for p, sample2 in enumerate(s2_namelist):
		S2_df = pd.read_csv(os.path.join(output_dir, sample2+"_"+AS+".csv"), delimiter="\t")
		S2_list[p] = list(S2_df.values.tolist())
		AS_file_length = len(S2_df)

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

		n1 = S1_n_total/len_S1
		N1 = S1_N_total/len_S1

		for p in range(len_S2):
			full_list = S2_list[p]
			S2_n_total += float(full_list[i][10])
			S2_N_total += float(full_list[i][11])

		n2 = S2_n_total/len_S2
		N2 = S2_N_total/len_S2

		##### added N filter: 11.14.2022
		##### both N should be >= 10 and no negative values of n and N are allowed
		if (min(N1, N2) >= 10) and (n1 < N1) and (n2 < N2):
			N1 = N1 + n1
			N2 = N2 + n2
			if N1>0 and N2>0:
				ratio_diff = (n1/N1) - (n2/N2)
				abs_ratio_diff = abs(ratio_diff)
				
				P0 = (n1+n2)/(N1+N2)
				n10 = N1 * P0
				n20 = N2 * P0
				exp = [n10, N1-n10, n20, N2-n20]
				if 0 not in exp:
					res = chisquare([n1, N1-n1, n2, N2-n2], f_exp=exp, ddof = 1)
					chrom_region_long = chrom+':'+gene+":"+str(start)+"-"+str(end)
					writer_list.append((chrom, gene, start, end, res[1], ratio_diff, abs_ratio_diff, n1, N1-n1, n2, N2-n2, chrom_region_long))
	
	df_out = pd.DataFrame(writer_list, columns = output_columns)
	df_out.to_csv(os.path.join(output_dir, AS+"_Output.csv"), sep='\t', index=False)

	for sample1 in s1_namelist:
		os.remove(os.path.join(output_dir, sample1+"_"+AS+".csv"))
	for sample2 in s2_namelist:
		os.remove(os.path.join(output_dir, sample2+"_"+AS+".csv"))
	
	return