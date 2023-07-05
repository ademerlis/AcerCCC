#!/bin/bash
#BSUB -J featurecounts_trimmed
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o featurecounts%J.out
#BSUB -e featurecounts%J.err
#BSUB -n 8

and="/scratch/projects/and_transcriptomics"

module load subread/1.6.2

featureCounts -t gene \
-g ID \
-a ${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 \
-o ${and}/Allyson_CCC/subread_counts/AcerCCC.counts \
${and}/Allyson_CCC/aligned/*Aligned.sortedByCoord.out.bam

