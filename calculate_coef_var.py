#!/usr/bin/env python
#this script calculates theta, di, ri, pi

import os
import re
import sys
import statistics
import argparse
import logging
from collections import defaultdict
from subprocess import Popen, PIPE
from operator import itemgetter
import numpy as np
import datetime
import pickle

logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

p = argparse.ArgumentParser()
p.add_argument("-o", "--output", required=True, help="path to output directory")
p.add_argument("-b", "--bams", required=True, help="file with one line per sample (tab delimited: absolute bam path and whether ice/agilent was used)")
p.add_argument("-i", "--intervals", required=True, help="specify if should run on ice or agilent")
p.add_argument("-d", "--datamash", required=True, help="path to datamash")

args = p.parse_args()
cov_directory = args.output
bam_files = args.bams
interval_of_interest = args.intervals
datamash = args.datamash


smn_coverages = "{cov_directory}/smn_doc_all.txt".format(**locals())

smn1_cov = {}
smn2_cov = {}


logging.info("Reading in coverage data...")
############################################
#       READ IN FILES
############################################

result, err = Popen(["""cat {smn_coverages} |{datamash} --sort --headers --group 1,3 mean 4 > smn12_doc_all.txt""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()

#get list of samples of interest
samples_of_interest = {}
with open(bam_files, 'r') as IN:
	for line in IN:
		line = line.rstrip()
		items = line.split('\t')
		bam,interval = items[0:2]
		if interval == interval_of_interest:
			sample = re.sub(".bam$","",os.path.basename(bam))
			samples_of_interest[sample] = 0

# read in the SMN1 and SMN2 mean/median gene coverage per sample
sample_list = {}
with open("{cov_directory}/smn12_doc_all.txt".format(**locals()), 'r') as IN:
	next(IN)
	for line in IN:
		line = line.rstrip()
		items = line.split('\t')
		sample,gene,cov = items[0:3]
		if gene == "SMN1" and sample in samples_of_interest:
			smn1_cov[sample] = cov
		elif gene == "SMN2" and sample in samples_of_interest:
			smn2_cov[sample] = cov
		sample_list[sample] = 0 

target_coverages = defaultdict(dict)
with open("{cov_directory}/smn_doc_all.txt".format(**locals()), 'r') as IN:
	for line in IN:
		line = line.rstrip()
		items = line.split('\t')
		sample,target,gene,coverage = items[0:4]
		info = gene + "_" + target
		target_coverages[sample][info] = coverage


## NEED TO ADD IN SUPPORT FOR FAILED SAMPLES
pre_files=[os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(cov_directory)) for f in fn]
files=[f for f in pre_files if os.path.isfile(f)  and os.path.basename(f).startswith("gene_cov") and f.split('/')[-2] in samples_of_interest]

logging.info("Calculating median coverage across all intervals...")
#get median coverage across all intervals per gene per sample
for file in files:
	filename = os.path.splitext(os.path.basename(file))[0]
	sample = re.sub("gene_cov_","",filename)
	cmd = """cat {file} | {datamash} --sort --headers --group 5 median 8 > {cov_directory}/{sample}/median_{filename}.txt""".format(**locals())
	result, err = Popen([cmd], stdout=PIPE, stderr=PIPE, shell=True).communicate()


#files of median(normalized_coverage)
#pre_files=[os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(cov_directory)) for f in fn]
medians = [f for f in pre_files if os.path.isfile(f)  and os.path.basename(f).startswith("median_gene_cov") and os.stat(f).st_size != 0 and f.split('/')[-2] in samples_of_interest]

sample_list = []


#get list of samples...remove this redundancy later with below chunk
for file in medians:
	match = re.search("median_gene_cov_([A-z\d\.-]+).txt", file)
	sample = match.group(1)
	sample_list.append(sample)


############################################
#       CHECK FOR VARIABILITY AT SMN SITES
############################################
logging.info("Checking variability in coverage at SMN sites...")
diff1_count = 0
diff2_count = 0
smn1_counts = {}
smn2_counts = {}
OUT = open("sample_smn_metrics.txt",'w')
for sample in sample_list:

	smn1_count = 0
	smn2_count = 0
	sample_info = target_coverages[sample]

	for region in sample_info:
		cov = float(target_coverages[sample][region])
		if re.search("SMN1", region):
			smn1_count += cov
		else:
			smn2_count += cov

	smn2_site1 = float(target_coverages[sample]["SMN2_5:69372304"])
	smn2_site2 = float(target_coverages[sample]["SMN2_5:69372353"])
	smn2_site3 = float(target_coverages[sample]["SMN2_5:69372501"])

	smn1_site1 = float(target_coverages[sample]["SMN1_5:70247724"])
	smn1_site2 = float(target_coverages[sample]["SMN1_5:70247773"])
	smn1_site3 = float(target_coverages[sample]["SMN1_5:70247921"])

	ratio_loci1 = smn1_site1/(smn1_site1 + smn2_site1)
	ratio_loci2 = smn1_site2/(smn1_site2 + smn2_site2)
	ratio_loci3 = smn1_site3/(smn1_site3 + smn2_site3)

	if abs(ratio_loci2 - ratio_loci1) > .10:
		final_ratio = ratio_loci2 # and output all the other specific info
		diff1_count += 1
		smn1_cov[sample] = smn1_site2 #average dictionary
		smn2_cov[sample] = smn2_site2 #average dictionary
		smn1_counts[sample] = smn1_site2 #raw count dictionary
		smn2_counts[sample] = smn2_site2 #raw count dictionary
		diff = "yes"
	elif abs(ratio_loci2 - ratio_loci3) > .10:
		final_ratio = ratio_loci2 # and output all the other speicfic info
		diff2_count += 1
		smn1_cov[sample] = smn1_site2 #average dictionary
		smn2_cov[sample] = smn2_site2 #average dictionary
		smn1_counts[sample] = smn1_site2 #raw count dictionary
		smn2_counts[sample] = smn2_site2 #raw count dictionary
		diff = "yes"
	else:
		# don't change pre-exisiting average coverage dictionary
		smn1_counts[sample] = smn1_count
		smn2_counts[sample] = smn2_count
		diff = "no"

	OUT.write("{sample}\t{smn1_site1}\t{smn1_site2}\t{smn1_site3}\t{smn2_site1}\t{smn2_site2}\t{smn2_site3}\t{diff}\n".format(**locals()))

OUT.close()

#get length value of 10% of samples
ten_percent = .10*len(medians)
coverages = defaultdict(list)
zkis = defaultdict(list)
sample_zkis = defaultdict(dict)


############################################
#       GET MEDIAN COVERAGE AND ZKI
############################################

logging.info("Calculating zki...")
#put medians for each gene into a list
for file in medians:
	match = re.search("median_gene_cov_([A-z\d\.-]+).txt", file)
	sample = match.group(1)
	with open(file,'r') as IN:
		next(IN)
		for line in IN:
			line = line.rstrip()
			items = line.split('\t')
			gene,median_cov = items[0:2]
			coverages[gene].append(float(median_cov))
			if sample != "IC2683758":

				sample_smn1 = smn1_cov[sample]
				sample_smn2 = smn2_cov[sample]

				if median_cov == "0":
					zkis[gene].append(0)
					sample_zkis[sample][gene]=0
				else:
					zki = (float(sample_smn1)+float(sample_smn2))/float(median_cov)
					zkis[gene].append(zki)
					sample_zkis[sample][gene] = zki

logging.info("Calculating MAD...")
def mad_sd(values_list):
	values_median = statistics.median(values_list)
	residuals = [abs(i-values_median) for i in values_list]
	mad = statistics.median(residuals)
	sd_convert = mad*1.4826
	return sd_convert

#for each gene, calculate median, mean, SD, and MAD, and zki mean/median
all_samples_median = {k:statistics.median(v) for k,v in coverages.iteritems()}
all_samples_MAD = {k:mad_sd(v) for k,v in coverages.iteritems()}
all_samples_median_zki = {k:statistics.median(v) for k,v in zkis.iteritems()}
all_samples_MAD_zki = {k:mad_sd(v) for k,v in zkis.iteritems()}


#PICKLE DUMP
with open("all_samples_median.txt", "wb") as fp: # Pickle for dev, remove later
		pickle.dump(all_samples_median, fp) 
with open("all_samples_MAD.txt", "wb") as fp: 
		pickle.dump(all_samples_MAD, fp) 
with open("all_samples_median_zki.txt", "wb") as fp: 
		pickle.dump(all_samples_median_zki, fp)
with open("all_samples_MAD_zki.txt", "wb") as fp: 
		pickle.dump(all_samples_MAD_zki, fp)


#PICKLE LOAD
with open("all_samples_median.txt", "rb") as fp: 
		all_samples_median = pickle.load(fp)
with open("all_samples_MAD.txt", "rb") as fp: 
		all_samples_MAD = pickle.load(fp)
with open("all_samples_median_zki.txt", "rb") as fp: 
		all_samples_median_zki = pickle.load(fp)
with open("all_samples_MAD_zki.txt", "rb") as fp: 
		all_samples_MAD_zki = pickle.load(fp)


logging.info("Calculating coefficient of variation...")
#for each gene, calculate SD/mean and MAD/median
coef_median = {k: float(all_samples_MAD[k])/float(all_samples_median[k]) if all_samples_median[k]!=0 else 0 for k in all_samples_MAD}
coef_median_zki = {k: float(all_samples_MAD_zki[k])/float(all_samples_median_zki[k]) if all_samples_median_zki[k]!=0 else 0 for k in all_samples_MAD_zki}




############################################
#       DECISION ON HOUSEKEEPING GENES
############################################
logging.info("Choosing housekeeping genes...")
#get 95th percentile value among gene medians
all_genes = [v for k,v in all_samples_median.iteritems()]
all_genes_95 = np.percentile(all_genes,95)
#number that defines bottom 5th percentile cutoff of median coverage
median_bottom_five_num = np.percentile(all_genes,5)
#get new list for all genes, only include values for samples with coverage less than the 5th percentile cutoff
median_bottom_five_per_gene = {k:filter(lambda a: a<median_bottom_five_num,v) for k,v in coverages.iteritems()}
#calculate how many samples for each gene in bottom 5th percentile coverage
median_samples_in_bottom = {k:len(v) for k,v in median_bottom_five_per_gene.iteritems()}

final_median_list = coef_median
final_median_list_zki = coef_median_zki
for gene,num in median_samples_in_bottom.items():
	if num > ten_percent:
		#delete gene from coef_median, need to check that deleted:
		del final_median_list[gene]
		del final_median_list_zki[gene]


#print list sorted by coef of variation
sorted_final_median_list = sorted(final_median_list.iteritems(), key=itemgetter(1), reverse=True)
sorted_final_median_list_zki = sorted(final_median_list_zki.iteritems(), key=itemgetter(1), reverse=True)

lowest_ten_median = sorted_final_median_list[-10:]
lowest_ten_median_zki = sorted_final_median_list_zki[-10:]
##### ADD PART TO GET ALL UNIQUE OUT OF THE TWO

house_keeping_genes = lowest_ten_median + lowest_ten_median_zki
house_keeping_genes = [i[0] for i in house_keeping_genes] # get just the gene and drop the value


OUT_MEDIAN_COV = open("rank_median_cov.txt",'w')
OUT_MEDIAN_ZKI = open("rank_median_zki.txt",'w')
OUT_MEDIAN_COV.write("gene\tcoef_var_cov\n")
OUT_MEDIAN_ZKI.write("gene\tcoef_var_zki\n")


for gene in sorted_final_median_list:
	OUT_MEDIAN_COV.write("{gene[0]}\t{gene[1]}\n".format(**locals()))

for gene in sorted_final_median_list_zki:
	OUT_MEDIAN_ZKI.write("{gene[0]}\t{gene[1]}\n".format(**locals()))


OUT_MEDIAN_COV.close()
OUT_MEDIAN_ZKI.close()

num_hk = len(house_keeping_genes)

pre_mean_zki_hk = {k:v for k,v in final_median_list_zki.iteritems() if k in house_keeping_genes}
mean_zki_hk = statistics.mean(pre_mean_zki_hk.values())


logging.info("Calculating variables...")#calculate pi
# Dbi= number of SMN1 reads
# rbi=number of SMN1 and SMN2 reads
# Di=theta Di
# pi=Di/ri
OUT_SAMPLE = open("{interval_of_interest}_sma_sample_stat.txt".format(**locals()), 'w')
OUT_SAMPLE.write("sample\tsmn1_reads\tsmn2_reads\ttheta\tdi\tri\tpi\n".format(**locals()))
for sample in sample_list:
	hk_sample_zki = []
	sum_zk = 0
	for gene in house_keeping_genes:
		zki_zk = (float(sample_zkis[sample][gene]))/float(mean_zki_hk)
		sum_zk += zki_zk
	theta = sum_zk/num_hk

	if theta > 1.0: #ceiling of 1 for theta
		theta = 1.0

	smn1_reads = smn1_counts[sample]
	smn2_reads = smn2_counts[sample]


	di = float(smn1_reads)*theta
	ri = float(smn1_reads) + float(smn2_reads)
	pi = di/ri
	OUT_SAMPLE.write("{sample}\t{smn1_reads}\t{smn2_reads}\t{theta}\t{di}\t{ri}\t{pi}\n".format(**locals()))
OUT_SAMPLE.close()
logging.info("DONE")
