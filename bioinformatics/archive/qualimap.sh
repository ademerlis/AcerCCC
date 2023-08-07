#!/bin/bash
#BSUB -J qualimap
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o qualimap%J.out
#BSUB -e qualimap%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_CCC/aligned"

data=($(ls *Aligned.sortedByCoord.out.bam))

module load java/7
module load R/4.1.0

for sample in ${data[@]} ;

do \
${and}/programs/qualimap_v2.3/qualimap bamqc \
-bam ${and}/Allyson_CCC/aligned/${sample} \
-outdir ${and}/Allyson_CCC/aligned \  
-gff ${and}/genomes/Acer/Acerv_assembly_v1.0.gff3 \  ; \

done
