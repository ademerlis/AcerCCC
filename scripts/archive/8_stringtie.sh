#!/bin/bash
#BSUB -J stringtie_updatedannotations_take3
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o stringtie_updatedannotations_take3%J.out
#BSUB -e stringtie_updatedannotations_take3%J.err
#BSUB -u and128@miami.edu
#BSUB -N

module load python/3.8.7

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch4_AcerCCC/aligned_updatedannotations_take3"

data=($(ls *Aligned.sortedByCoord.out.bam))

for i in ${data[@]} ;

do \
/scratch/projects/and_transcriptomics/programs/stringtie-2.2.1/stringtie -p 8 -e -B -G /scratch/projects/and_transcriptomics/genomes/Acer/Acerv.GFFannotations.fixed_transcript_take3.gff3 -A ${i}.gene_abund.tab -o ${i}.gtf ${i} ; \
done
