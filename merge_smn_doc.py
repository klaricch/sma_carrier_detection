#!/usr/bin/env python
# this script merges the SMN coverage reuslts at 6 postiions into one file with all samples

import re
import os
import sys
from os.path import isfile


OUT=open("smn_doc_all.txt", 'w')
OUT.write("sample\ttarget\tgene\tcoverage\n")

bam=sys.argv[1]
output_dir=sys.argv[2]

with open(bam, 'r') as IN:
	for line in IN:
		line=line.rstrip()
		items=line.split('\t')
		bam =items[0]
		sample = (os.path.splitext(bam)[0]).split("/")[-1] # fix line


		print sample
		cov="{output_dir}/{sample}/smn_{sample}.sample_interval_summary".format(**locals())
		if isfile(cov):
			with open(cov, 'r') as IN:
				next(IN)
				for line in IN:
					line=line.rstrip()
					items=line.split('\t')
					target,coverage=items[0:2]

					if target == "5:70247724" or target == "5:70247773" or target == "5:70247921":
						gene = "SMN1"
					else:
						gene = "SMN2"

					OUT.write("{sample}\t{target}\t{gene}\t{coverage}\n".format(**locals()))
OUT.close()



