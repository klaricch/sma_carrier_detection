#!/usr/bin/env python
# this script merges the SMN coverage results at 6 positions into one file with all samples

import re
import os
import sys
import argparse
from os.path import isfile


p = argparse.ArgumentParser()
p.add_argument("-b", "--bams", required=True, help="file with one line per sample (tab delimited: absolute bam path and whether ice/agilent was used)")
p.add_argument("-o", "--output", required=True, help="path to output directory")

args = p.parse_args()
bam = args.bams
output_dir = args.output


OUT = open("{output_dir}/smn_doc_all.txt".format(**locals()), 'w')
OUT.write("sample\ttarget\tgene\tcoverage\n")


with open(bam, 'r') as IN:
	for line in IN:
		line = line.rstrip()
		items = line.split('\t')
		bam = items[0]
		sample = (os.path.splitext(bam)[0]).split("/")[-1] # fix line


		print sample
		cov = "{output_dir}/{sample}/smn_{sample}.sample_interval_summary".format(**locals())
		if isfile(cov):
			with open(cov, 'r') as IN:
				next(IN)
				for line in IN:
					line = line.rstrip()
					items = line.split('\t')
					target,coverage = items[0:2]

					if target == "5:70247724" or target == "5:70247773" or target == "5:70247921":
						gene = "SMN1"
					else:
						gene = "SMN2"

					OUT.write("{sample}\t{target}\t{gene}\t{coverage}\n".format(**locals()))
OUT.close()



