#!/bin/bash
#BSUB -J multiqc_aligned_updatedannotations
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o multiqc_STARalign_updatedannotations%J.out
#BSUB -e multiqc_STARalign_updatedannotations%J.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

cd /scratch/projects/and_transcriptomics/Allyson_CCC/aligned_updatedannotations

multiqc .
