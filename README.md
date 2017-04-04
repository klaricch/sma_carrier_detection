# SMA Carrier Detection

The scripts in this repository can be used to implement the methods described in Larson et al. 2015 (https://www.ncbi.nlm.nih.gov/pubmed/26510457) for detecting SMA carriers. This technique utilizes both carrier probabilities and coverage at SMN1 loci to investigate SMA carrier status. (in beta)


### run_smn_doc.sh
Calculate coverge per gene and at three SMN loci that distinguish SMN1 from SMN2
##### Input Files:
1) bam_list # file with one line per sample (tab delimited: absolute bam path and whether ice/agilent was used)
2) GATKjar # location of GATK jar installation
3) reference # path to reference hg37
4) sma_intervals # smn loci interval file
5) output_dir #name of output directory
6) picard #location of picard jar installation
7) scripts_dir #location of sma scripts

### merge_smn_doc.py
Merge SMN coverage results from all samples into one file
##### Input Files:
1) bam_list # file with one line per  sample (tab delimited: absolute bam path and whether ice/agilent was used)
2) output_dir #name of output directory (keep consistent with previous scripts)


### calculate_coef_var.py
Calculate theta, di, ri, pi
##### Input Files:
1) cov_directory  #name of output directory (keep consistent with previous scripts)
2) bam_files # file with one line per  sample (tab delimited: absolute bam path and whether ice/agilent was used)
3) interval_of_interest #specify if should run on ice or agilent

#### calculate_carrier_probability.R
Calculates the carrier probabilitiy and plots credible intervals


#### Requirements:
GATK
Picard
datamash

