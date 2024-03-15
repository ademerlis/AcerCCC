#!/bin/bash
#BSUB -J trim_polyA
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -o trim_polyA.out
#BSUB -e trim_polyA.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

for sample in ${and}/Allyson_CCC/trimmed/trimmed_fastq_files/*.gz ;

do \

${and}/programs/TrimGalore-0.6.10/trim_galore ${sample} \
--polyA ; \

done
