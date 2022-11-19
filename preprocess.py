import os
import time


def SamtoText(input_dir, current, bamfile_name, chromosomes):
	output_dir = os.path.join(input_dir, os.path.splitext(bamfile_name)[0])
	os.makedirs(output_dir, exist_ok=True)
	samtools_dir = os.path.join(current,"samtools")
	cmd1 = samtools_dir+" index "+os.path.join(input_dir, bamfile_name)		## make samtools index filename.bam.bai
	os.system(cmd1)
	print("bam index file generated.")

	for chrom in chromosomes:
		print("Chrom ", chrom, "... ...")
		tt = time.time()
		cmd2 = samtools_dir+" view -b "+os.path.join(input_dir, bamfile_name)+" "+chrom+" -o "+os.path.join(output_dir, chrom+".bam")
		cmd3 = samtools_dir+" pileup "+os.path.join(output_dir, chrom+".bam")+" | cut -f 2,4 > "+os.path.join(output_dir, chrom+".txt")    ### Need to use pileup, not mpileup
		command = cmd2+";"+cmd3
		os.system(command)
		os.system("rm "+os.path.join(output_dir, chrom+".bam"))
	print("Read coverage files generated.")
	return
