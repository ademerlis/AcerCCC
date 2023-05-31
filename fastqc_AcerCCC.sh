#!/bin/bash
#~/scratch/projects/and_transcriptomics/Allyson_CCC/scripts/fastqc/fastqc_AcerCCC.job
#/scratch/projects/and_transcriptomics/Allyson_CCC/scripts/fastqc/fastqc_AcerCCC.job
#purpose: quality checking of raw RNAseq reads using FASTQC on Pegasus compute node

#BSUB -J AcerCCC_fastqc
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o fastqc_AcerCCC.out
#BSUB -e fastqc_AcerCCC.err
#BSUB -n 8

and="/scratch/projects/and_transcriptomics/Allyson_CCC/" 

cd ${and}
for SAMP in *.fastq.gz

do

module load /share/Modules/general/java/1.8.0_60
module load /share/Modules/bio/fastqc/0.10.1 \
${and}/$SAMP \
--outdir ${and}/fastqc_results
done
