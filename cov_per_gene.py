#!/usr/bin/env python
#this script calculates coverage for gencode target intervals

import sys
import os.path
import argparse
from subprocess import Popen, PIPE


p = argparse.ArgumentParser()
p.add_argument("-b", "--bams", required=True, help="file with one line per sample (tab delimited: absolute bam path and whether ice/agilent was used)")
p.add_argument("-r", "--reference", required=True, help="path to reference fasta")
p.add_argument("-i", "--intervals", required=True, help="specify if should run on ice or agilent")
p.add_argument("-o", "--output", required=True, help="path to output directory")
p.add_argument("-p", "--picard", required=True, help="path to picard installation")

args = p.parse_args()
bam=args.bams
reference=args.reference
target_intervals=args.intervals
output_dir=args.output
picard=args.picard

sample=os.path.splitext(os.path.basename(bam))[0]
bait_name=os.path.basename("{target_intervals}".format(**locals()))

cmd="""java -Xmx10g -jar {picard} CollectHsMetrics \
BAIT_INTERVALS={target_intervals} \
BAIT_SET_NAME={bait_name} \
TARGET_INTERVALS={target_intervals} \
PER_TARGET_COVERAGE={output_dir}/{sample}/gene_cov_{sample}.txt \
INPUT={bam} \
OUTPUT={output_dir}/{sample}/{sample}.hybrid_selection_metrics \
TMP_DIR={output_dir} \
REFERENCE_SEQUENCE={reference} \
MINIMUM_MAPPING_QUALITY=20 \
MINIMUM_BASE_QUALITY=20 \
CLIP_OVERLAPPING_READS=true \
NEAR_DISTANCE=250 \
COVERAGE_CAP=200 \
SAMPLE_SIZE=10000 \
VERBOSITY=INFO \
QUIET=false \
VALIDATION_STRINGENCY=STRICT \
COMPRESSION_LEVEL=5 \
MAX_RECORDS_IN_RAM=500000 \
CREATE_INDEX=false \
CREATE_MD5_FILE=false \
GA4GH_CLIENT_SECRETS=client_secrets.json""".format(**locals())
print cmd
result, err = Popen([cmd], stdout=PIPE, stderr=PIPE, shell=True).communicate()
print result
print err

print "DONE"