#!/bin/bash
#BSUB -J AcerCCC_multiqc
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o multiqc_AcerCCC.out
#BSUB -e multiqc_AcerCCC.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

cd /scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/

multiqc .
