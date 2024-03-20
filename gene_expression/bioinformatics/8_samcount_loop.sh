#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch4_AcerCCC/3_bowtie2/alignment"

cd "/scratch/projects/and_transcriptomics/Ch4_AcerCCC/3_bowtie2/alignment"

data=($(ls *.sam))

for samp in "${data[@]}" ; do \

#build script
echo "making sam_counts script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_samcounts
#BSUB -e ${and}/Ch4_AcerCCC/3_bowtie2/alignment/logs/${samp}_samcounts.err
#BSUB -o ${and}/Ch4_AcerCCC/3_bowtie2/alignment/logs/${samp}_samcounts.out
#BSUB -W 12:00
#BSUB -n 8
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch4_AcerCCC/3_bowtie2/alignment\"

module load samtools/1.3

perl samcount.pl ${samp} /scratch/projects/and_transcriptomics/genomes/Host_concat_seq2iso.tab aligner=bowtie2 >${samp}.counts

" > ${and}/Ch4_AcerCCC/3_bowtie2/alignment/${samp}_samcounts.job

bsub < ${and}/Ch4_AcerCCC/3_bowtie2/alignment/${samp}_samcounts.job

done
