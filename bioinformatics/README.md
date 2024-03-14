## Bioinformatics pipeline for *A.cervicornis* Coral City Camera vs. Emerald Reef Samples

Script written by: DeMerlis

Last updated: 20240224

Following [Dr. Matz tag-based_RNAseq](https://github.com/z0on/tag-based_RNAseq), [Dr. Studivan's annotated pipeline for tag-based RNAseq](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt), and [Dr. Natalia Andrade's Lexogen QuantSeq trimming guidelines](https://github.com/China2302/SCTLD_RRC/tree/main/hpc). 

Also added my own steps (FastQC and multiqc). 

**Pipeline**: FastQC -> [countreads.pl](https://github.com/z0on/tag-based_RNAseq/blob/master/countreads.pl) -> TrimGalore -> [countreads_trim.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/countreads_trim.pl) -> FastQC -> download and format reference genome or transcriptome -> bowtie2 for index and alignment -> SAMtools for generating counts matrix -> DESeq2

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

Script 1: using TrimGalore to remove low-quality base pairs and adapters.

Specific flags:
--illumina \
--cores 4 \
--three_prime_clip_R1 12 \
--nextseq 30 \
--length 20 

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

cd "/scratch/projects/and_transcriptomics/Ch4_AcerCCC/1_fastq_rawreads"

for sample in ${and}/Ch4_AcerCCC/1_fastq_rawreads/*.gz ;

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

Script 2: Removing polyA tails

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

for sample in ${and}/Ch4_AcerCCC/2_trimmed_reads/trimmed_trimgalore_NJ/trimmed_fastq_files/*.gz ;

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

and="/scratch/projects/and_transcriptomics/Ch4_AcerCCC"

module load fastqc/0.10.1

cd ${and}/2_trimmed_reads/trimmed_trimgalore_NJ/trimmed_and_removedpolyA_fastqfiles

for f in *.fq.gz;
do \
mv -v -- "$f" "${f%.fq.gz}.fastq.gz"; \

done

fastqc *.fastq.gz
--outdir ${and}/trimmed/
```

MultiQC reports can be seen [here](https://github.com/ademerlis/AcerCCC/tree/main/bioinformatics/multiqc_reports). Need to download the html files first, then open them in browser.

## 5. Count trimmed reads

```{bash}
#!/bin/bash
#BSUB -J countreads_trim
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o countreads_trim.out
#BSUB -e countreads_trim.err
#BSUB -u and128@miami.edu
#BSUB -N

#Purpose: counts the number of Illumina reads in trimmed fastq files

#specify variables and paths

cd "/scratch/projects/and_transcriptomics/Ch4_AcerCCC/2_trimmed_reads/trimmed_trimgalore_NJ/trimmed_and_removedpolyA_fastqfiles"

output_file="countreads_results.txt"

# Default file pattern
glob="\.gz"

# Check if an argument is provided
if [ "$1" ]; then
    glob="$1"
fi

# Loop through files matching the pattern
for f in *$glob*; do
    
    # Count the number of lines in the file
    nrd=$(zcat "$f" | wc -l)

    # Divide the line count by 4
    nrd=$((nrd / 4))

    # Print the filename and the calculated number
    echo -e "$f\t$nrd" >> "$output_file"
done

echo "Results have been saved to $output_file"
```

![Screen Shot 2024-03-14 at 10 53 15 AM](https://github.com/ademerlis/AcerCCC/assets/56000927/260f79f8-70d0-4d17-9070-359e8ea0b966)


## 5. Download *Acropora cervicornis* genome files

Obtained most recent Baums lab genome. (Locatelli et al. *in prep*).

## 6. Build *Acropora cervicornis* + *Symbiodinium fitti* index with Bowtie2

```{bash}
#!/usr/bin/env bash
#BSUB -e bowtie2build_LocatelliShoguchiConcat.err
#BSUB -o bowtie2build_LocatelliShoguchiConcat.out
#BSUB -P and_transcriptomics
#BSUB -q bigmem

workdir="/scratch/projects/and_transcriptomics"

bowtie2-build ${workdir}/genomes/Acer/Locatelli_2023/Acer_Genome/Acropora_cervicornis.mrna-transcripts.fa, \
${workdir}/genomes/Symbiodinium/syma_transcriptome_37.fasta \
Locatelli_Shoguchi_concat
```

## 7. Run alignment with Bowtie2

**Note:** the .bt2 files from the bowtie-build step need to be in the same directory as the trimmed sequences, and cannot be in a subfolder. Otherwise, the alignment code won't be able to identify them. Also, the .bt2 files all need to have the same name before the extensions (i.e. "Locatelli_Shoguchi_concat" in my example), and that is what you write in the alignment script for it to find the files.

```{bash}
#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch4_AcerCCC/2_trimmed_reads/trimmed_trimgalore_NJ/trimmed_and_removedpolyA_fastqfiles"

cd "/scratch/projects/and_transcriptomics/Ch4_AcerCCC/2_trimmed_reads/trimmed_trimgalore_NJ/trimmed_and_removedpolyA_fastqfiles"

data=($(ls *.gz))

for samp in "${data[@]}" ; do \

#build script
echo "making bowtie2-align script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_bowtie2align_LocatelliShoguchi
#BSUB -e ${projdir}/bowtie2align_LocatelliShoguchi/logs/${samp}_bowtie2align_LocatelliShoguchi.err
#BSUB -o ${projdir}/bowtie2align_LocatelliShoguchi/logs/${samp}_bowtie2align_LocatelliShoguchi.out
#BSUB -W 12:00
#BSUB -q bigmem

cd \"/scratch/projects/and_transcriptomics/Ch4_AcerCCC/2_trimmed_reads/trimmed_trimgalore_NJ/trimmed_and_removedpolyA_fastqfiles\"

bowtie2 --local -U ${samp} -x Locatelli_Shoguchi_concat --un ${samp}.unaligned -k 5 -S ${samp}.sam

" > ${projdir}/bowtie2align_LocatelliShoguchi/${samp}_bowtie2align_LocatelliShoguchi.job

bsub < ${projdir}/bowtie2align_LocatelliShoguchi/${samp}_bowtie2align_LocatelliShoguchi.job

done
```


