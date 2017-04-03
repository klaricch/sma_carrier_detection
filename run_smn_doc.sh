#!/usr/bin/env bash
#this script calculates coverge per gene and at three SMN loci that distinguish SMN1 from SMN2

bam_list=$1 # file with one line per sample (tab delimited: absolute bam path and whether ice/agilent was used)
GATKjar=$2 # location of GATK jar installation
reference=$3 # path to reference hg37
sma_intervals=$4 # smn loci interval file
output_dir=$5 #name of output directory
picard=$6 #location of picard jar installation
scripts_dir=$7 #location of sma scripts


mkdir logs
> ${output_dir}/smn_doc_all.txt


ice_targets= #provide ice targets file
agilent_targets= # provide agilent target file

while read bam interval; do 
	sample=$(basename $bam .bam)
	echo $sample
	echo $bam
	qsub -q gsa -V -o ${output_dir}/logs/smn_doc_${sample}_o.txt -e ${output_dir}/logs/smn_doc_${sample}_e.txt ${scripts_dir}/smn_doc.sh $GATKjar $reference $sma_intervals $bam $output_dir

	if [[ $interval == "ice" ]];
		then \
			echo "ice"; \
			qsub -q gsa -V -l h_vmem=20g -o ${output_dir}/logs/cov_per_gene_${sample}_o.txt -e ${output_dir}/logs/cov_per_gene_${sample}_e.txt -S $(which python) ${scripts_dir}/cov_per_gene.py $bam $reference $ice_targets $output_dir $picard ; \

		elif [[ $interval == "agilent" ]];
			then \
				echo "agilent"; \
				qsub -q gsa -V -l h_vmem=20g -o ${output_dir}/logs/cov_per_gene_${sample}_o.txt -e ${output_dir}/logs/cov_per_gene_${sample}_e.txt -S $(which python) ${scripts_dir}/cov_per_gene.py $bam $reference $agilent_targets $output_dir $picard; \
		else
			echo "CUSTOM INTERVAL FILE"
	fi; \

done <$bam_list

echo "DONE"