#!/usr/bin/env bash
#this script calculates coverge per gene and at three SMN loci that distinguish SMN1 from SMN2

bam_list=$1 # file with one line per sample (tab delimited: absolute bam path and whether ice/agilent was used)
GATKjar=$2 # location of GATK jar installation
reference=$3 # path to reference hg37
sma_intervals=$4 # smn loci interval file
output_dir=$5 # name of output directory
picard=$6 # location of picard jar installation
scripts_dir=$7 # location of sma scripts


mkdir logs
echo ${output_dir}


ice_targets=${scripts_dir}/whole_exome_illumina_coding_v1.chr1_22.gencode.v19.targets.interval_list #provide ice targets file
agilent_targets=${scripts_dir}/whole_exome_agilent_1.1_refseq_plus_3_boosters.chr1_22.gencode.v19.targets.interval_list # provide agilent target file

while read bam interval; do 
	sample=$(basename $bam .bam)
	echo $sample
	echo $bam
	qsub -q gsa -V -o ${output_dir}/logs/smn_doc_${sample}_o.txt -e ${output_dir}/logs/smn_doc_${sample}_e.txt ${scripts_dir}/smn_doc.sh $bam $GATKjar $reference $sma_intervals $output_dir

	if [[ $interval == "ice" ]];
		then \
			echo "ice"; \
			qsub -q gsa -V -l h_vmem=20g -o ${output_dir}/logs/cov_per_gene_${sample}_o.txt -e ${output_dir}/logs/cov_per_gene_${sample}_e.txt -S $(which python) ${scripts_dir}/cov_per_gene.py -b $bam -r $reference -i $ice_targets -o $output_dir -p $picard; \

		elif [[ $interval == "agilent" ]];
			then \
				echo "agilent"; \
				qsub -q gsa -V -l h_vmem=20g -o ${output_dir}/logs/cov_per_gene_${sample}_o.txt -e ${output_dir}/logs/cov_per_gene_${sample}_e.txt -S $(which python) ${scripts_dir}/cov_per_gene.py -b $bam -r $reference -i $agilent_targets -o $output_dir -p $picard; \
		else
			echo "CUSTOM INTERVAL FILE"
	fi; \

done <$bam_list

echo "DONE"