#!/usr/bin/env bash
# this script uses GATK DepthOfCoverage to calculate read depth at 6 loci used to distinguish SMN1 and SMN2
#$ -cwd
#$ -o logs/
#$ -e logs/

bam=$1
GATKjar=$2
reference=$3
sma_intervals=$4
output_dir=$5



sample=$(basename $bam .bam)
mkdir ${output_dir}/${sample}


echo $GATKjar
echo $reference
echo $sma_intervals
echo $bam
echo $output_dir

#get coverage at SMN loci with GATk
java -jar $GATKjar \
   	-T DepthOfCoverage \
   	-R $reference \
   	-o ${output_dir}/${sample}/smn_${sample} \
   	-I $bam \
   	-L $sma_intervals

echo "DONE"
