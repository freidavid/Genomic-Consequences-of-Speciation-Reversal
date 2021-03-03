#!/bin/bash
#Bash script to submit all possible d-stat comparisons, David Frei 2021
#No inputs needed
#example: sh /cluster/work/gdc/shared/p641/scripts/d_stats.sh

scr_path="/path/to/working/directory"
index="./WF_wtdbg2.chr.fasta.fai"
bam_path="/path/to/bam/files/bam"

	
#Limit stuff to only the chromosomes included in the global_intersect.txt file
rm ${scr_path}/chromosomes.txt
cat ${index} | cut -f1 | sed -n '1,3p' > ${scr_path}/chromosomes.txt
cat ${index} | cut -f1 | sed -n '5,6p' >> ${scr_path}/chromosomes.txt
cat ${index} | cut -f1 | sed -n '8,16p' >> ${scr_path}/chromosomes.txt
cat ${index} | cut -f1 | sed -n '18,21p' >> ${scr_path}/chromosomes.txt
cat ${index} | cut -f1 | sed -n '23,27p' >> ${scr_path}/chromosomes.txt
cat ${index} | cut -f1 | sed -n '29,31p' >> ${scr_path}/chromosomes.txt
cat ${index} | cut -f1 | sed -n '33,34p' >> ${scr_path}/chromosomes.txt
cat ${index} | cut -f1 | sed -n '36p' >> ${scr_path}/chromosomes.txt
cat ${index} | cut -f1 | sed -n '39,40p' >> ${scr_path}/chromosomes.txt





#get a file with only the sites we want to use (no missing data, file from step 7)
zcat ./min_32_ind.vcf.gz | cut -f 1-2 > ${scr_path}/global_intersect.txt
#Index intersect file
angsd sites index ${scr_path}/global_intersect.txt


#define outgroup
outgroup="salmon"


#Make bam file lists manually, separated by populations
#kilch
ls ${bam_path}/DF5_complete_nova_hiseq.bam > ${scr_path}/kilch_o.list
ls ${bam_path}/DF11_complete_nova_hiseq.bam >> ${scr_path}/kilch_o.list
ls ${bam_path}/DF1_complete_nova.bam >> ${scr_path}/kilch_o.list
ls ${bam_path}/DF3_complete_nova.bam >> ${scr_path}/kilch_o.list
ls ${bam_path}/DF4_complete_nova.bam >> ${scr_path}/kilch_o.list
ls ${bam_path}/DF6_complete_nova.bam >> ${scr_path}/kilch_o.list
ls ${bam_path}/DF7_complete_nova.bam >> ${scr_path}/kilch_o.list
ls ${bam_path}/DF8_complete_nova.bam >> ${scr_path}/kilch_o.list
ls ${bam_path}/DF9_complete_nova.bam >> ${scr_path}/kilch_o.list
ls ${bam_path}/DF12_complete_nova.bam >> ${scr_path}/kilch_o.list
ls ${bam_path}/DF20_complete_hiseq.bam >> ${scr_path}/kilch_o.list


#Sandfelchen
ls ${bam_path}/DF19_complete_hiseq.bam > ${scr_path}/sandfelchen_o.list
ls ${bam_path}/DF30_complete_hiseq.bam >> ${scr_path}/sandfelchen_o.list
ls ${bam_path}/DF31_complete_hiseq.bam >> ${scr_path}/sandfelchen_o.list
#Gangfisch
ls ${bam_path}/DF17_complete_hiseq.bam > ${scr_path}/gangfisch_o.list
ls ${bam_path}/DF18_complete_hiseq.bam >> ${scr_path}/gangfisch_o.list
#Blaufelchen
ls ${bam_path}/DF23_complete_hiseq.bam > ${scr_path}/blaufelchen_o.list
ls ${bam_path}/DF24_complete_hiseq.bam >> ${scr_path}/blaufelchen_o.list

#Sandfelchen
ls ${bam_path}/DF123477_hiseq.bam > ${scr_path}/sandfelchen_n.list
ls ${bam_path}/DF123440_hiseq.bam >> ${scr_path}/sandfelchen_n.list
ls ${bam_path}/DF126_complete_nova.bam >> ${scr_path}/sandfelchen_n.list
ls ${bam_path}/DF127_complete_nova.bam >> ${scr_path}/sandfelchen_n.list
ls ${bam_path}/DF128_complete_nova.bam >> ${scr_path}/sandfelchen_n.list
#Gangfisch
ls ${bam_path}/DF123470_complete_hiseq.bam > ${scr_path}/gangfisch_n.list
ls ${bam_path}/DF123458_complete_hiseq.bam >> ${scr_path}/gangfisch_n.list
ls ${bam_path}/DF132_complete_nova.bam >> ${scr_path}/gangfisch_n.list
#Blaufelchen
ls ${bam_path}/DF123446_complete_hiseq.bam > ${scr_path}/blaufelchen_n.list
ls ${bam_path}/DF123448_complete_hiseq.bam >> ${scr_path}/blaufelchen_n.list
ls ${bam_path}/DF121_complete_nova.bam >> ${scr_path}/blaufelchen_n.list
ls ${bam_path}/DF122_complete_nova.bam >> ${scr_path}/blaufelchen_n.list
ls ${bam_path}/DF123_complete_nova.bam >> ${scr_path}/blaufelchen_n.list
ls ${bam_path}/DF131_complete_nova.bam >> ${scr_path}/blaufelchen_n.list

#salmon outgroup
ls ${bam_path}/salmon_complete.bam > ${scr_path}/salmon.list



#Add Kilch comparisons by hand, as it's not symmetrical because there is no contemporary
#Kilch population
echo "kilch" > $scr_path/donor.list
echo "kilch" >> $scr_path/donor.list
echo "kilch" >> $scr_path/donor.list

echo "blaufelchen" > $scr_path/recipient.list
echo "gangfisch" >> $scr_path/recipient.list
echo "sandfelchen" >> $scr_path/recipient.list



#Create all pairwise combinations of populations
set -- "sandfelchen" "gangfisch" "blaufelchen"
for POP1; do
    #shift
    for POP2; do 
    if [ "$POP1" != "$POP2" ]; then
  	echo "$POP1" >> $scr_path/donor.list  
  	echo "$POP2" >> $scr_path/recipient.list 
  	
  	fi  
    done
done



#count how many jobs to submit
n=$(cat $scr_path/donor.list | wc -l)
seq $n > $HOME/jobindices
#Write all numbers from the $jobindices file into one comma separated string for the bsub
#command
index=$(awk -vORS=, '{ print $0 }' $HOME/jobindices | sed 's/,$/\n/')


bsub -n 1 -W 120:00 -N -u "david.frei@eawag.ch" -J "ABBA_BABA[${index}]" -R "rusage[mem=120000]" -o ./D-stat_LOGFILE.txt < ./d-stat._final.lsf

