import os, sys
import time


def SamtoText(input_dir, current, bamfile_name, chromosomes):
	output_dir = os.path.join(input_dir, os.path.splitext(bamfile_name)[0])
	last_dir = os.path.basename(os.path.normpath(bamfile_name)).split('.bam')[0]
	os.makedirs(last_dir, exist_ok=True)

	samtools_dir = os.path.join(current, "samtools")
	try:
		cmd1 = samtools_dir+" index "+os.path.join(current, input_dir, bamfile_name)		## make samtools index filename.bam.bai
		os.system(cmd1)
		#print("bam index file generated.")
	except ValueError:
		print("Index file could not be generated")
		sys.exit()

	for chrom in chromosomes:
		#print("Chrom ", chrom, "... ...")
		tt = time.time()
		cmd2 = samtools_dir+" view -b "+os.path.join(current, input_dir, bamfile_name)+" "+chrom+" -o "+os.path.join(current, output_dir, chrom+".bam")
		cmd3 = samtools_dir+" pileup "+os.path.join(current, output_dir, chrom+".bam")+" | cut -f 2,4 > "+os.path.join(current, output_dir, chrom+".txt")    ### Need to use pileup, not mpileup
		command = cmd2+";"+cmd3
		#print(command)
		try:
			os.system(command)
			os.system("rm "+os.path.join(current, output_dir, chrom+".bam"))
		except ValueError:
			print("Read coverage file could not be generated")
			sys.exit()

	print("Read coverage files generated for", bamfile_name)
	return

