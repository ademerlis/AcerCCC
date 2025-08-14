#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch4_AcerCCC/2_trimmed_reads/trimmed_trimgalore_NJ/trimmed_and_removedpolyA_fastqfiles"

cd "/scratch/projects/and_transcriptomics/Ch4_AcerCCC/2_trimmed_reads/trimmed_trimgalore_NJ/trimmed_and_removedpolyA_fastqfiles"

data=($(ls *.gz))

for samp in "${data[@]}" ; do \

#build script
echo "making bowtie2-align script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_bowtie2align_LocatelliShoguchi
#BSUB -e ${projdir}/bowtie2align_LocatelliShoguchi/logs/${samp}_bowtie2align_LocatelliShoguchi.err
#BSUB -o ${projdir}/bowtie2align_LocatelliShoguchi/logs/${samp}_bowtie2align_LocatelliShoguchi.out
#BSUB -q bigmem

cd \"/scratch/projects/and_transcriptomics/Ch4_AcerCCC/2_trimmed_reads/trimmed_trimgalore_NJ/trimmed_and_removedpolyA_fastqfiles\"

bowtie2 --local -U ${samp} -x /scratch/projects/and_transcriptomics/Ch4_AcerCCC/2_trimmed_reads/trimmed_trimgalore_NJ/trimmed_and_removedpolyA_fastqfiles/Locatelli_Shoguchi_concat --un ${samp}.unaligned -k 5 -S ${samp}.sam

" > ${projdir}/bowtie2align_LocatelliShoguchi/${samp}_bowtie2align_LocatelliShoguchi.job

bsub < ${projdir}/bowtie2align_LocatelliShoguchi/${samp}_bowtie2align_LocatelliShoguchi.job

done
