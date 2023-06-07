#!/bin/bash
#BSUB -J trim_CCC
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -o trim_CCC.out
#BSUB -e trim_CCC.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

tworemainingsamples="2383"

echo "These are the samples to be processed:"
echo $tworemainingsamples

for sample in ${and}/Allyson_CCC/fastq_rawreads/${tworemainingsamples[*]}.fastq.gz ;

do \
cd ${and}/Allyson_CCC/fastq_rawreads/
${and}/programs/TrimGalore-0.6.10/trim_galore ${sample}
--illumina \
--cores 4 \
--three_prime_clip_R1 12 \
--nextseq 30 \
--length 20 \
--outdir ${and}/Allyson_CCC/trimmed/ ; \

done

cd ${and}/Allyson_CCC/trimmed/

mv *.fq.gz *.fastq.gz

module load fastqc/0.10.1

fastqc *.fastq.gz

multiqc ${and}/Allyson_CCC/trimmed/ \
--outdir ${and}/Allyson_CCC/trimmed/
