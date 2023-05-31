#!/bin/bash
#~/scripts/fastqc_AcerCCC.sh
#/scratch/projects/and_transcriptomics/Allyson_CCC/scripts/fastqc_AcerCCC.sh
#purpose: quality checking of raw RNAseq reads using FASTQC on Pegasus compute node

#BSUB -J AcerCCC_fastqc
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o fastqc_AcerCCC.out
#BSUB -e fastqc_AcerCCC.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics/Allyson_CCC"

cd ${and}
fastqc *.fastq.gz
--outdir ${and}/fastqc/

