#!/bin/bash
#BSUB -J stringtie_reestimateassembly
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie_reestimateassembly%J.out
#BSUB -e stringtie_reestimateassembly%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_CCC/aligned/"

data=($(ls *Aligned.sortedByCoord.out.bam))

for i in ${data[@]} ;

do \
${and}/programs/stringtie-2.2.1/stringtie -e -G /scratch/projects/and_transcriptomics/Allyson_CCC/aligned/stringtie_gtf_files/stringtie_acerv_merged.gtf -o ${i}_reestimated.merge.gtf ${i} ; \ 
done
