#!/bin/bash
#BSUB -J star_align_rawfastq
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o star_align_rawfastq%J.out
#BSUB -e star_align_rawfastq%J.err
#BSUB -u and128@miami.edu
#BSUB -N

# try aligning raw reads to genome and see how alignment rates compare to trimmed reads
# A soft clipping option is added to STAR to deal with any polyA and adapter contamination

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Allyson_CCC/fastq_rawreads"

data=($(ls *.gz))

for sample in ${data[@]} ;

do \
/scratch/projects/and_transcriptomics/programs/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
--genomeDir ${and}/genomes/Acer/Acer_STAR_index/ \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 16 \
--readFilesCommand gunzip -c \
--readFilesIn ${sample} \
--quantMode TranscriptomeSAM GeneCounts \
--clip3pAdapterSeq AAAAAAAAAA \
--outReadsUnmapped Fastx \
--outFileNamePrefix ${and}/Allyson_CCC/aligned/aligned_rawfastq/${sample} ; \

done
