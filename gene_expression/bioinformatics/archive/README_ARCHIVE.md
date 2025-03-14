## Bioinformatics pipeline for *A.cervicornis* Coral City Camera vs. Emerald Reef Samples

Script written by: DeMerlis

Last updated: 20230808

Note: the file/folder hierarchies and names may have changed since uploading this, as I seem to keep editing the hierarchy (even while scripts are running and then they fail!!).

I followed the pipelines of [Dr. Natalia Andrade](https://github.com/China2302/SCTLD_RRC/tree/main/hpc) and [Jill Ashey](https://github.com/JillAshey/SedimentStress/blob/master/Bioinf/RNASeq_pipeline_FL.md?plain=1) for this analysis.

**Pipeline**: [FastQC](https://github.com/ademerlis/AcerCCC/tree/main/bioinformatics#1-fastqc-raw-reads) -> [TrimGalore (adapters + low-quality bp)](https://github.com/ademerlis/AcerCCC/tree/main/bioinformatics#2-trimgalore-part-1) -> [TrimGalore (polyA tail)](https://github.com/ademerlis/AcerCCC/tree/main/bioinformatics#3-trimgalore-part-2) -> [FastQC](https://github.com/ademerlis/AcerCCC/tree/main/bioinformatics#4-fastqc-trimmed-reads) -> [STAR](https://github.com/ademerlis/AcerCCC/tree/main/bioinformatics#7-star-index-genome) -> [DESeq2]()

Some notes: I had to edit the .gff3 file from the Acer genome because it was missing the "Parent_ID=" and "Transcript_ID=" flags that STAR specifically looks for during alignment.

## 1. FastQC Raw Reads

```{bash}
#!/bin/bash
#~/scripts/fastqc_AcerCCC.sh
#/scratch/projects/and_transcriptomics/Allyson_CCC/scripts/fastqc_AcerCCC.sh
#purpose: quality checking of raw RNAseq reads using FASTQC on Pegasus compute node

#BSUB -J AcerCCC_fastqc
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o fastqc_AcerCCC.out
#BSUB -e fastqc_AcerCCC.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

module load fastqc/0.10.1

and="/scratch/projects/and_transcriptomics/Allyson_CCC"

cd ${and}
fastqc *.fastq.gz
--outdir ${and}/fastqc/
```

## 2. TrimGalore part 1

I installed TrimGalore locally:

```{bash}
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
```

Then I ran this code, which first trims with specifications:

1. `--illumina` flag indicates to TrimGalore to use the correct quality score offset for Illumina-encoded quality scores (Phred+64)
2. `--cores 4`: This flag specifies the number of CPU cores to use for the trimming process.
3. `--three_prime_clip_R1 12`: This flag indicates that a 3'-end (end of read) clipping of 12 bases should be applied to the first read (R1) in a paired-end dataset. This means that the last 12 bases of the R1 read will be removed during the trimming process.
4. `--nextseq 30`:  specifies the Phred quality score threshold for trimming
5. `--length 20`: the minimum read length to retain after trimming. Any read that becomes shorter than 20 bases after trimming will be discarded.

Then, the second part of the code renames resulting trimmed files (.fq.gz) into .fastq.gz files so they can be read by FastQC. 

```{bash}
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

for sample in ${and}/Allyson_CCC/fastq_rawreads/*.gz ;

do \
${and}/programs/TrimGalore-0.6.10/trim_galore ${sample}
--illumina \
--cores 4 \
--three_prime_clip_R1 12 \
--nextseq 30 \
--length 20 \ ; \

done

cd ${and}/Allyson_CCC/scripts

for f in *.fq.gz;
do \
mv -v -- "$f" "${f%.fq.gz}.fastq.gz"; \

done
```

## 3. TrimGalore part 2

For some reason it doesn't work if you add the `--polyA` flag into the script above, so you have to run it twice. Idk why.

```{bash}
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
```

## 4. FastQC Trimmed Reads

```{bash}
#!/bin/bash
#BSUB -J fastqc_trimmed_polyA
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o fastqc_trimmed_polyA.out
#BSUB -e fastqc_trimmed_polyA.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics/Allyson_CCC"

module load fastqc/0.10.1

cd ${and}/trimmed/trimmed_and_removedpolyA_fastqfiles

for f in *.fq.gz;
do \
mv -v -- "$f" "${f%.fq.gz}.fastq.gz"; \

done

fastqc *.fastq.gz
--outdir ${and}/trimmed/
```

## 5. Download [*Acropora cervicornis*](https://usegalaxy.org/u/skitch/h/acervicornis-genome) genome files

Obtained from [Baums lab](http://baumslab.org/research/data/) with permission from Dr. Sheila Kitchen. Using Version v1.0_171209

Genome file: `Acerv_assembly_v1.0_171209.fasta`

GFF file: `Acerv_assembly_v1.0.gff3`

Protein file: `Acerv_assembly_v1.0.protein.fa`

Transcript file: `Acerv_assembly_v1.0.mRNA.fa`


## 6. Fix Acer .gff3 file

See [this blog post](https://github.com/ademerlis/ademerlis.github.io/blob/master/_posts/2023-06-27_stringtiecode_andredoingSTARwithupdatedannotationfile.md) for a detailed description as to why this needs to be done.

```{r}
# Title: A. cervicornis GFF adjustments
# Project: Sedimentation RNA-Seq / Mcap 2020
# Author: J. Ashey / A. Huffmyer -> Allyson DeMerlis
# Date: 06/28/2023

# Need to do some acerv gff adjustments so it can run properly in STAR. Here, I'll be adding transcript_id= to 'gene' column because STAR needs that label to run

#Load libraries
library(tidyverse)

#Load  gene gff
Acerv.gff <- read.csv(file="Downloads/Galaxy1-[Acerv_assembly_v1.0.gff3].gff3", header=FALSE, sep="\t", skip=1) 

#rename columns
colnames(Acerv.gff) <- c("scaffold", "Gene.Predict", "id", "gene.start","gene.stop", "pos1", "pos2","pos3", "gene")
head(Acerv.gff)

# Creating transcript id
Acerv.gff$transcript_id <- sub(";.*", "", Acerv.gff$gene)
Acerv.gff$transcript_id <- gsub("ID=", "", Acerv.gff$transcript_id) #remove ID= 
Acerv.gff$transcript_id <- gsub("Parent=", "", Acerv.gff$transcript_id) #remove Parent=
head(Acerv.gff)

# Create Parent ID 
Acerv.gff$parent_id <- sub(".*Parent=", "", Acerv.gff$gene)
Acerv.gff$parent_id <- sub(";.*", "", Acerv.gff$parent_id)
Acerv.gff$parent_id <- gsub("ID=", "", Acerv.gff$parent_id) #remove ID= 
head(Acerv.gff)

Acerv.gff <- Acerv.gff %>% 
  mutate(gene = ifelse(id != "gene", paste0(gene, ";transcript_id=", Acerv.gff$transcript_id, ";gene_id=", Acerv.gff$parent_id),  paste0(gene)))
head(Acerv.gff)

Acerv.gff<-Acerv.gff %>%
  select(!transcript_id)%>%
  select(!parent_id)

#save file
write.table(Acerv.gff, file="~/Downloads/Acerv.GFFannotations.fixed_transcript_take3.gff3", sep="\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
```

## 7. STAR Index Genome

First install STAR locally.

```{bash)
# Get latest STAR source from releases
#ran this all in login node 
wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
tar -xzf 2.7.10b.tar.gz
cd STAR-2.7.10b
# Compile
cd source
make STAR
```

Then run this code to index the genome.

The last three flags were added based on other people's codes:
1. `--sjdbOverhang 100`: the length of the splice junction overhang. A larger overhang can help improve the accuracy of alignment across splice junctions.
2. `--sjdbGTFtagExonParentTranscript Parent`: the "Parent" tag should be used to associate exons with their parent transcripts in the GTF file
3. `--genomeSAindexNbases 13`:  This flag determines the length of the "seed" used during the construction of the suffix array index. The suffix array is a data structure used by STAR for rapid string matching and alignment. The "seed" is a short sequence of nucleotides that is used to quickly identify potential alignment locations. The value "13" specifies that a 13-base seed will be used.

```{bash}
#!/bin/bash
#BSUB -J Acer_star_index_fixedannotations_take3
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -o /scratch/projects/and_transcriptomics/genomes/Acer/star_index%J.out
#BSUB -e /scratch/projects/and_transcriptomics/genomes/Acer/star_index%J.err
#BSUB -u and128@miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

/scratch/projects/and_transcriptomics/programs/STAR-2.7.10b/bin/Linux_x86_64_static/STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ${and}/genomes/Acer/Acer_STAR_index_gffannotations.fixed_take3 \
--genomeFastaFiles ${and}/genomes/Acer/Acerv_assembly_v1.0_171209.fasta \
--sjdbGTFfile ${and}/genomes/Acer/Acerv.GFFannotations.fixed_transcript_take3.gff3 \
--sjdbOverhang 100 \
--sjdbGTFtagExonParentTranscript Parent \
--genomeSAindexNbases 13
```

## 8. STAR Alignment

```{bash}
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
```
