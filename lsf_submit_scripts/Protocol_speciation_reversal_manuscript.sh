#!/bin/bash
#=========================================================================================
#=========================================================================================
#Submit script protocol for "Genomic consequences of speciation reversal"", David Frei, 2021
#=========================================================================================
#=========================================================================================


#=========================================================================================
#1 Trim polyg tails with fastp
#=========================================================================================
#Actually this is only useful in novaseq data. However, I ran it for all invididuals to 
#treat all the same.
#Nova
bsub -n 2 -W 24:00 -N -u "david.frei@eawag.ch" -J "fastp_n[1,3,4,5,6,7,8,9,11,12,121,122,123,126,127,128,131,132]%6" -R "rusage[mem=20000] span[hosts=1]" -o ./fastp_PROTOCOL_novaseq.txt < ./fastp_n.lsf

#Hiseq
bsub -n 2 -W 24:00 -N -u "david.frei@eawag.ch" -w "ended(fastp_nova)" -J "fastp_h[5,11,17,18,19,20,23,24,30,31,123446,123448,123470,123477,123458,123440]%6" -R "rusage[mem=20000] span[hosts=1]" -o ./fastp_PROTOCOL_hiseq.txt < ./fastp_h.lsf


#=========================================================================================
#2 Merging overlapping reads with Seqprep
#=========================================================================================
#First step: merge overlapping read pairs. This is especially important for the historical 
#samples, could theoretically be negelcted for recent samples, as fragments are longer and 
#there should be no overlapping reads
#novaseq
bsub -n 1 -W 24:00 -w "ended(fastp_h)" -J "SeqPrep_n[1,3,4,5,6,7,8,9,11,12,121,122,123,126,127,128,131,132]%6" -R "rusage[mem=22000] span[hosts=1]" -o ./SeqPrep_PROTOCOL_novaseq.txt < ./Seqprep_n.lsf
#hiseq
bsub -n 1 -W 24:00 -w "ended(SeqPrep_n)" -J "SeqPrep_h[5,11,17,18,19,20,23,24,30,31,123446,123448,123470,123477,123458,123440]%6" -R "rusage[mem=22000] span[hosts=1]" -o ./SeqPrep_PROTOCOL_hiseq.txt < ./Seqprep_h.lsf


#=========================================================================================
#3 Mapping to Rishi's whitefish reference with BWA
#=========================================================================================
#Second step: Align reads to the reference genome, in this case the Alpine whitefish 
#reference genome (C. steinmanni, De-Kayne et al. 2021).
#novaseq
bsub -n 10 -W 120:00 -w "ended(SeqPrep_h)" -J "BWA_n[1,3,4,5,6,7,8,9,11,12,121,122,123,126,127,128,131,132]%6" -R "rusage[mem=6500,scratch=10000] span[hosts=1]" -o ./BWA_PROTOCOL_novaseq.txt < ./BWA_n.lsf
#hiseq
bsub -n 10 -W 120:00 -w "ended(BWA_n)" -J "BWA_h[5,11,17,18,19,20,23,24,30,31,123446,123448,123470,123477,123458,123440]%6" -R "rusage[mem=6500,scratch=10000] span[hosts=1]" -o ./BWA_PROTOCOL_hiseq.txt < ./BWA_h.lsf


#=========================================================================================
#4 Processing the aligned bam files
#=========================================================================================
#novaseq
bsub -n 1 -W 24:00 -w "ended(BWA_h)" -J "Processing_n[1,3,4,5,6,7,8,9,11,12,121,122,123,126,127,128,131,132]%5" -R "rusage[mem=60000,scratch=100000] span[hosts=1]" -o ./Processing_novaseq_Protocol.txt < ./Processing_bam_novaseq.lsf
#hiseq
bsub -n 1 -W 24:00 -w "ended(Processing_n)" -J "Processing_h[5,11,17,18,19,20,23,24,30,31,123446,123448,123470,123477,123458,123440]%5" -R "rusage[mem=60000,scratch=100000] span[hosts=1]" -o ./Processing_hiseq_Protocol.txt < ./Processing_bam_hiseq.lsf


#=========================================================================================
#5 Merging bam files
#=========================================================================================
#Script that merges the bam files of merged reads and the one containing forward and 
#reverse reads. Additionally, the two samples that were sequenced on both hiseq and 
#novaseq have to be merged to, as Angsd only takes one bam per sample.
#novaseq:
bsub -n 2 -W 12:00 -w "ended(Processing_h)" -J "Merging_nova[1,3,4,5,6,7,8,9,11,12,121,122,123,126,127,128,131,132]%12" -R "rusage[mem=30000,scratch=30000] span[hosts=1]"  -o ./Merging_PROTOCOL.txt < ./Merging_bams_novaseq.lsf
#hiseq:
bsub -n 2 -W 12:00 -w "ended(Merging_nova)" -J "Merging_hiseq[5,11,17,18,19,20,23,24,30,31,123446,123448,123470,123477,123458,123440]%12" -R "rusage[mem=30000,scratch=30000] span[hosts=1]"  -o ./Merging_PROTOCOL.txt < ./Merging_bams_hiseq.lsf
#Merge hiseq and novaseq bam of DF5 and DF11
bsub -n 2 -W 12:00 -w "ended(Merging_hiseq)" -J "Merging_DF5_and_DF11[5,11]" -R "rusage[mem=30000,scratch=30000] span[hosts=1]"  -o ./Merge_DF11_and_DF5_PROTOCOL.txt < ./Merge_DF5_and_DF11.lsf


#=========================================================================================
#6 Genotype likelihoods in beagle format with Angsd for PCA
#=========================================================================================
#Calculation of genotypelikelihoods, inlcuding filtering, with Angsd
bsub -n 1 -W 12:00 -w "ended(Merging_hiseq)" -J "ANGSD_Per_Chromosome[1-20]%7" -R "rusage[mem=20000] span[hosts=1]" -o ./Angsd_protocol_1.txt < ./Angsd_per_chromosome.lsf
bsub -n 1 -W 12:00 -w "ended(ANGSD_Per_Chromosome)" -J "ANGSD_Per_Chromosome_2[21-40]%9" -R "rusage[mem=13000] span[hosts=1]" -o ./Angsd_protocol_2.txt < ./Angsd_per_chromosome.lsf


#=========================================================================================
#7 Merge single beagle files for each chromsome to a genome wide beagle/vcf file
#=========================================================================================
#for no missing data
bsub -n 1 -W 2:00 -w "ended(ANGSD_Per_Chromosome_2)" -J "Merging_single_beagle_min32" -o /dev/null < ./merging_single_beagle_min32.lsf
#for 2 individuals allowed to be missing (only names are different, otherwise the same script)
bsub -n 1 -W 2:00 -w "ended(Merging_single_beagle_min32)" -J "Merging_single_beagle_min30" -o /dev/null < ./merging_single_beagle_min30.lsf


#=========================================================================================
#8 PCA based on genotype likelihoods with PCAngsd
#=========================================================================================
bsub -n 1 -W 2:00 -w "ended(Merging_single_beagle_min30)" -J "PCAngsd" -R "rusage[mem=1000]" -o ./PCAngsd_LOGFILE.txt < ./PCAngsd.lsf


#=========================================================================================
#9:Also calculate genotype likelihoods of the polymorphic sites from step 5 in the salmon outgroup 
#=========================================================================================
bsub -n 1 -W 12:00 -w "ended(PCAngsd)" -J "ANGSD_Per_Chromosome_salmon[1-20]%7" -R "rusage[mem=20000] span[hosts=1]" -o ./Angsd_protocol_salmon_1.txt < ./Angsd_per_chromosome_salmon.lsf
bsub -n 1 -W 12:00 -w "ended(ANGSD_Per_Chromosome_salmon)" -J "ANGSD_Per_Chromosome_salmon_2[21-40]%9" -R "rusage[mem=13000] span[hosts=1]" -o ./Angsd_protocol_salmon_2.txt < ./Angsd_per_chromosome_salmon.lsf


#=========================================================================================
#10 Merge single beagle files again from step 9
#=========================================================================================
bsub -n 1 -W 2:00 -w "ended(ANGSD_Per_Chromosome_salmon_2)" -J "Merging_single_beagle_min30_salmon" -o /dev/null < ./merging_single_beagle_min30_salmon.lsf


#=========================================================================================
#11 D Stats in Angsd
#=========================================================================================
#Custom bash script that submits all needed pairwise comparisons of dstats
sh ./d_stats.sh


#=========================================================================================
#12 Phasing and calling genotype from genotype likelihoods in beagle
#=========================================================================================
bsub -n 1 -W 12:00 -w "ended(ABBA_BABA)" -J "BEAGLE_per_Chromosome_1[1-20]" -R "rusage[mem=2000]" -o /dev/null < ./Beagle_per_chromosome.lsf
bsub -n 1 -W 12:00 -w "ended(BEAGLE_Per_Chromosome_1)" -J "BEAGLE_Per_Chromosome2[21-40]" -R "rusage[mem=1500]" -o /dev/null < ./Beagle_per_chromosome.lsf



#=========================================================================================
#13 Merging phased vcf files inot a genome wide vcf 
#=========================================================================================
bsub -n 1 -W 2:00 -w "ended(BEAGLE_Per_Chromosome2)" -J "Merging_single_phased_vcf_min30_salmon" -o /dev/null < ./merging_single_phased_vcf_min_30.lsf
bsub -n 1 -W 2:00 -w "ended(Merging_single_phased_vcf_min30_salmon)" -J "Merging_single_phased_vcf_only_kilch" -o /dev/null < ./merging_single_phased_only_kilch.lsf




#=========================================================================================
#14 TWISST
#=========================================================================================
#Transform to .geno.gz format
bsub -n 1 -W 12:00 -J "Parse_VCF.py" -w "ended(Merging_single_phased_vcf_only_kilch)"  -R "rusage[mem=1000]" -o ./Parse_VCF_Protocol.txt < ./parse_vcf.lsf
#Run PhyML in windows along the genome
bsub -n 1 -W 48:00 -J "Sliding_PhyML" -w "ended(Parse_VCF.py)" -R "rusage[mem=12250]" -o ./Sliding_PhyML_Protocol.txt < ./Sliding_PhyML.lsf
#Actual TWISST run:
bsub -n 1 -W 48:00 -J "TWISST" -w "ended(Sliding_PhyML)" -R "rusage[mem=1000]" -o ./TWISST_Protocol.txt < ./twisst.lsf


#=========================================================================================
#15 Selscan
#=========================================================================================
bsub -n 1 -W 2:00 -w "ended(Sliding_PhyML)"-J "Selscan" -R "rusage[mem=1000]" -o ./Selscan_LOGFILE.txt < ./selscan.lsf


















