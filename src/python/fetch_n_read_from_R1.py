from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input file of read1")
ap.add_argument("-n", "--number", required=True, help="how many reads will be fetch")
ap.add_argument("-o1", "--output_txt", required=True, help="read name output file")
ap.add_argument("-o2", "--output_fq", required=True, help="fastq output file")

args = vars(ap.parse_args())

input_r2 = args["input"]
num = args["number"]
output_file_txt = args["output_txt"]
output_file_fq = args["output_fq"]

n_counts = 0
with gzopen(input_r2, "rt") as fi, open(output_file_fq, "w") as out_handle, open(output_file_txt, "w") as fo:
	for title, seq, qual in FastqGeneralIterator(fi):
		n_counts = n_counts + 1
		if n_counts <= int(num):
			read_1_name = title.split(" ")[0]
			read_2_name = read_1_name + '\n'
			out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
			fo.write(read_2_name)
		else:
			break