## Bioinformatics pipeline for *A.cervicornis* Coral City Camera vs. Emerald Reef Samples

Script written by: DeMerlis

Last updated: 20240224

Following [Dr. Matz tag-based_RNAseq](https://github.com/z0on/tag-based_RNAseq), [Dr. Studivan's annotated pipeline for tag-based RNAseq](https://github.com/mstudiva/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt), and [Dr. Natalia Andrade's Lexogen QuantSeq trimming guidelines](https://github.com/China2302/SCTLD_RRC/tree/main/hpc). 

Also added my own steps (FastQC and multiqc). 

**Pipeline**: FastQC -> [countreads.pl](https://github.com/z0on/tag-based_RNAseq/blob/master/countreads.pl) -> TrimGalore -> [countreads_trim.pl](https://github.com/mstudiva/tag-based_RNAseq/blob/master/countreads_trim.pl) -> FastQC -> download and format reference genome or transcriptome -> bowtie2 for index and alignment -> SAMtools for generating counts matrix -> DESeq2

Programs I downloaded locally onto my HPC environment:
- multiQC (version 1.14)
- TrimGalore (version 0.6.10)
- Bowtie2 (version 2.5.2)

Programs I loaded in Pegasus environment that were already installed:
- fastQC (version 0.10.1)
- samtools (version 1.3)

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

- "--illumina": This option tells Trim Galore to remove Illumina universal adapters from the reads. Illumina sequencing platforms use a standard adapter sequence, and this flag ensures that any remnants of these adapters are trimmed off.
- "--cores 4": This flag sets the number of processing cores to use for the task. By specifying 4, you are instructing Trim Galore to use four cores, which can speed up the processing if your computer has multiple cores available.
- "--three_prime_clip_R1 12": This option instructs Trim Galore to trim 12 bases from the 3'-end of read 1. This is useful for removing any unwanted sequences or low-quality bases from the end of the reads.
- "--nextseq 30": This flag is used for NextSeq 500/550 data, which has two-color chemistry and can produce specific types of quality issues. The 30 value tells Trim Galore to trim bases at the 3'-end of each read that have a quality score of 30 or below. This is a more aggressive quality trimming approach suitable for NextSeq data to ensure high-quality output.
- "--length 20":  This option sets the minimum length of reads to keep after trimming. Reads that end up shorter than 20 bases after adapter and quality trimming will be discarded. This helps to ensure that only reads of sufficient length and quality are retained for subsequent analysis steps.

(thx ChatGPT)

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


## 5. Download *Acropora cervicornis* and *Symbiodinium* genome files

Obtained most recent Baums lab genome. (Locatelli et al. *in prep*).

Downloaded *Symbiodinium* clade A3 from [Shoguchi et al. 2021](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4857-9). 

[Link to download page](https://marinegenomics.oist.jp/symb/viewer/download?project_id=37)
"symatranscriptome_37.fasta.gz"

## 6. Build *Acropora cervicornis* + *Symbiodinium* index with Bowtie2

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

## 8. Count aligned reads and alignment rates

While Dr. Matz and Dr. Studivan have "countreads.pl" scripts that they used, I found that those codes did not result in accurate numbers for my aligned reads and rates. 

If you look at the .err files from the *_bowtie2align_LocatelliShoguchi.job files from the previous step, you will see this information:

![Screen Shot 2024-03-15 at 10 39 23 AM](https://github.com/ademerlis/AcerCCC/assets/56000927/27b6f885-6c52-467b-b137-9373f11318d1)

This information is the correct alignment rates and counts.

To extract this information from your .err files, run this script "alignment_extractinfo.sh" (be sure to run chmod +x to make it executable. Then, to run this script directly in your command line, run `./alignment_extractinfo.sh`)

```{bash}
#!/bin/bash

# Define the output file
output_file="summary.csv"

# Write the header to the output file
echo "Sample ID,Total Reads,Unpaired,Aligned 0 times,Aligned 1 time,Aligned >1 times,Overall Alignment Rate" > "$output_file"

# Loop through .err files in the current directory
for file in *.err
do
    # Extract the sample ID from the file name (assuming it's a number)
    sample_id=$(echo "$file" | grep -o '[0-9]\+')

    # Extract the necessary lines and format the output
    awk -v id="$sample_id" '/reads; of these:/ { total_reads = $1 }
         /were unpaired; of these:/ { unpaired = $1 }
         /aligned 0 times/ { aligned_0 = $1 }
         /aligned exactly 1 time/ { aligned_1 = $1 }
         /aligned >1 times/ { aligned_more_1 = $1 }
         /overall alignment rate/ { overall_rate = $1 }
         END {
            print id "," total_reads "," unpaired "," aligned_0 "," aligned_1 "," aligned_more_1 "," overall_rate
         }' "$file" >> "$output_file"
done

# Output the result
echo "Extraction complete. Data saved in $output_file"
```

Then I copied this information into an excel spreadsheet for my records.

#### A bit about Bowtie2 and aligning to a contatenated Host-Symbiont genome (for my own sanity's sake):

- "By default, Bowtie 2 performs end-to-end read alignment. That is, it searches for alignments involving all of the read characters. This is also called an "untrimmed" or "unclipped" alignment. When the **--local option** is specified, Bowtie 2 performs local read alignment. In this mode, Bowtie 2 might "trim" or "clip" some read characters from one or both ends of the alignment if doing so maximizes the alignment score."
- "For an alignment to be considered "valid" (i.e. "good enough") by Bowtie 2, it must have an alignment score no less than the minimum score threshold."
- "In local alignment mode, the default minimum score threshold is 20 + 8.0 * ln(L), where L is the read length."
- "When we say that a read has multiple alignments, we mean that it has multiple alignments that are valid and distinct from one another."
- "By default, Bowtie 2 searches for distinct, valid alignments for each read. When it finds a valid alignment, it generally will continue to look for alignments that are nearly as good or better. It will eventually stop looking, either because it exceeded a limit placed on search effort (see -D and -R) or because it already knows all it needs to know to report an alignment. Information from the best alignments are used to estimate mapping quality (the MAPQ SAM field) and to set SAM optional fields, such as AS:i and XS:i. Bowtie 2 does not guarantee that the alignment reported is the best possible in terms of alignment score."
- "In -k mode, Bowtie 2 searches for up to N distinct, valid alignments for each read, where N equals the integer specified with the -k parameter. That is, if -k 2 is specified, Bowtie 2 will search for at most 2 distinct alignments. It reports all alignments found, in descending order by alignment score."
- "This mode [-k mode] can be effective and fast in situations where the user cares more about whether a read aligns (or aligns a certain number of times) than where exactly it originated."

[Source](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-aligner)

So, when it comes to determining whether a read aligned to the Host genome or the Symbiont genome, this is actually not "filtered out" in the bowtie2 alignment step. Since -k is set to 5 in the below code, this means that a read can align to the concatenated Host-Symbiont genome up to 5 times successfully (which means that a read could align to both the Host and the Symbiont genome). When looking at percent alignment rates, the rates that are usually reported are actually the combined rates of aligning 1 time and aligning multiple times (see screenshot below). This is because the assumption is that aligning multiple times does not equal duplicated reads -- de-duplication should have occurred in the trimming step (if necessary -- depends on the type of sequencing performed). 

The filtering of whether a read is a true Host or Symbiont read, and the discarding of reads that align to both "equally well" (and thus should not be counted as a read) actually occurs in the next step of the process using SAMtools. When generating the read counts per gene, you first create a concatenated file that has the sequence ID corresponding gene/isogroup ID for each species. Then, when you run [samcount.pl](https://github.com/z0on/tag-based_RNAseq/blob/master/samcount.pl), there are some defaults that are specified in the code itself:

- **Dup.reads:** You can specify whether you want to keep or remove exact sequence-duplicate reads mapping to the same position in the reference genome/transcriptome. By default, the setting is “keep” (as duplicates are supposed to have been removed during the trimming stage).
- **Aligner:** “samcount.pl” by default assumes that you used Bowtie2 in -k mode for your alignment step. 
- **Mult.iso:** By default, if reads map to multiple isogroups, then the read is disregarded. You can change this setting to “random” if you want an isogroup to be randomly selected and used to assign a count to.

So, based on the "mult.iso" default, if a read maps to multiple isogroups (for example both a Host isogroup AND a symbiont isogroup), then that read is discarded. This is the step where we determine if a read is truly host or truly symbiont.

## 9. Generating read counts per gene

First, you need to make a "two-column tab-delimited table transcriptome_seq2gene.tab giving correspondence between entries in the transcriptome fasta file and genes. In de novo transcriptomes, several fasta contigs may correspond to the same gene (e.g., splice isoforms, or assembly uncertainties)." (from [Dr. Matz](https://github.com/z0on/tag-based_RNAseq/blob/master/tagSeq_processing_README.txt))

To generate these files for *A.cervicornis* and *S.fitt*, run this code on the fasta files:

```{bash}
# making seq2iso.tab files
grep ">" Acropora_cervicornis.mrna-transcripts.fa | perl -pe 's/>FUN_(\d+)(\S+)\s.+/FUN_$1$2\t FUN_$1/'>Acervicornis_seq2iso.tab
grep ">" syma_transcriptome_37.fasta | perl -pe 's/>comp(\d+)(\S+)\s.+/comp$1$2\t comp$1/'>Symbiodinium_seq2iso.tab

# create combo file

cat Acer/Locatelli_2023/Acer_Genome/Acervicornis_seq2iso.tab Symbiodinium/Symbiodinium_seq2iso.tab > Host_concat_seq2iso.tab
```

Next, download [samcount.pl](https://github.com/z0on/tag-based_RNAseq/blob/master/samcount.pl).

Then, run the following script to create .counts files:

```{bash}
#! /usr/bin/env bash

#define variables for directories and files
and="/scratch/projects/and_transcriptomics"
project="and_transcriptomics"
projdir="/scratch/projects/and_transcriptomics/Ch4_AcerCCC/3_bowtie2/alignment"

cd "/scratch/projects/and_transcriptomics/Ch4_AcerCCC/3_bowtie2/alignment"

data=($(ls *.sam))

for samp in "${data[@]}" ; do \

#build script
echo "making sam_counts script for ${samp}..."
echo "
#! /usr/bin/env bash
#BSUB -P ${project}
#BSUB -J ${samp}_samcounts
#BSUB -e ${and}/Ch4_AcerCCC/3_bowtie2/alignment/logs/${samp}_samcounts.err
#BSUB -o ${and}/Ch4_AcerCCC/3_bowtie2/alignment/logs/${samp}_samcounts.out
#BSUB -W 12:00
#BSUB -q general

cd \"/scratch/projects/and_transcriptomics/Ch4_AcerCCC/3_bowtie2/alignment\"

module load samtools/1.3

perl samcount.pl ${samp} /scratch/projects/and_transcriptomics/genomes/Host_concat_seq2iso.tab aligner=bowtie2 >${samp}.counts

" > ${and}/Ch4_AcerCCC/3_bowtie2/alignment/${samp}_samcounts.job

bsub < ${and}/Ch4_AcerCCC/3_bowtie2/alignment/${samp}_samcounts.job

done
```

Then, run these lines of code directly in the command line:

```{bash}
perl expression_compiler.pl *.counts > allcounts.txt

# let's remove those annoying chains of extensions from sample names
cat allcounts.txt | perl -pe 's/\_trimmed\.fastq\.gz\.sam\.counts//g'> counts.txt

#and rename FUN -> Acropora and comp -> symbiodinium
sed -i 's/FUN/Acropora/g' counts.txt
sed -i 's/comp/Symbiodinium/g' counts.txt
```

Then, use scp to move the counts.txt file to your local drive.
