#!/bin/bash
#BSUB -J trim_CCC
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o trim_CCC%J.out
#BSUB -e trim_CCC%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="scratch/projects/and_transcriptomics"

for sample in ${and}/Allyson_CCC/fastq_files/*.gz ;

do \

${and}/programs/TrimGalore-0.6.10/trim_galore ${sample}
--gzip \
--fastqc \
--fastqc_args "--outdir ${and}/Acer_CCC/trimmed/" \
--illumina \
--cores 4 \
--three_prime_clip_R1 12 \
--nextseq 30 \
--length 20 \
-o ${and}/Acer_CCC/trimmed/ ; \

done

multiqc ${and}/Acer_CCC/trimmed/ \
--outdir ${and}/Acer_CCC/trimmed/
