## Bioinformatics pipeline for *A.cervicornis* Coral City Camera vs. Emerald Reef Samples

Script written by: DeMerlis

Last updated: 20240224

Following [Dr. Matz tag-based_RNAseq](https://github.com/z0on/tag-based_RNAseq) and [Dr. Studivan's annotated pipeline for tag-based RNAseq](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt). 

Also added my own steps (FastQC and multiqc). 

Note: I manually downloaded my sequence files from Illumina BaseSpace because I couldn't figure out how to use a script to do it.

**Pipeline**: FastQC -> [countreads.pl](https://github.com/z0on/tag-based_RNAseq/blob/master/countreads.pl) -> trimming: [tagseq_clipper.pl](https://github.com/z0on/tag-based_RNAseq/blob/master/tagseq_clipper.pl) + cutadapt -> [countreads_trim.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/countreads_trim.pl) -> download and format reference genome or transcriptome -> bowtie2 for index and alignment -> 

## 1. FastQC Raw Reads

```{bash}
module load fastqc/0.10.1
fastqc *.fastq.gz
--outdir ${and}/fastqc/

cd ${and}/fastqc/
multiqc .
```
This creates an html file that you can download from HPC and open in web browser. Download and view html files for this project [here](https://github.com/ademerlis/AcerCCC/tree/main/bioinformatics/multiqc_reports). 

## 2. Countreads.pl

```{bash}
#!/bin/bash
#BSUB -J countrawreads
#BSUB -q bigmem
#BSUB -P and_transcriptomics
#BSUB -n 8
#BSUB -W 120:00
#BSUB -o countrawreads.out
#BSUB -e countrawreads.err
#BSUB -u and128@miami.edu
#BSUB -N

#Purpose: counts the number of Illumina reads in a bunch of fastq files

#specify variables and paths

and="/scratch/projects/and_transcriptomics"

cd "/scratch/projects/and_transcriptomics/Ch4_AcerCCC/1_fastq_rawreads"

output_file="countreads_results.txt"

glob=".fastq.gz"
if [ ! -z "$1" ]; then
    glob="$1"
fi

fqs=(*$glob)
for f in "${fqs[@]}"; do
    gunzip -c "$f" > "temp.fastq"  # Decompress the file to a temporary file
    nrd=$(cat "temp.fastq" | wc -l)
    nrd=$((nrd / 4))
    echo "$f    $nrd"
    echo "$f    $nrd" >> "$output_file"  # Append the results to the output file
    rm "temp.fastq"  # Remove the temporary file
done

echo "Results have been saved to $output_file"
```

![Screen Shot 2024-02-24 at 1 35 07 PM](https://github.com/ademerlis/AcerCCC/assets/56000927/742688a5-b533-46b7-bd54-12b5df6f4fbc)


## 3. Trimming






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
