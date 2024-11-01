#!/bin/bash
#BSUB -J star_align_trimmed_fixedannotations_take3
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o star_align%J.out
#BSUB -e star_align%J.err
#BSUB -u and128@miami.edu
#BSUB -N

# A soft clipping option is added to STAR to deal with any leftover polyA and adapter contamination 

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch4_AcerCCC/trimmed/trimmed_and_removedpolyA_fastqfiles/forSTAR"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
/scratch/projects/and_transcriptomics/programs/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
--genomeDir ${and}/genomes/Acer/Acer_STAR_index_gffannotations.fixed_take3/ \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 16 \
--readFilesCommand gunzip -c \
--readFilesIn ${sample} \
--quantMode TranscriptomeSAM GeneCounts \
--clip3pAdapterSeq AAAAAAAAAA \
--outReadsUnmapped Fastx \
--outFileNamePrefix ${and}/Ch4_AcerCCC/aligned_updatedannotations_take3/${sample} ; \

done
