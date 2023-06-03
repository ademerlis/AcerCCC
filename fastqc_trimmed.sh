#!/bin/bash
#~/scripts/fastqc_trimmed.sh
#/scratch/projects/and_transcriptomics/Allyson_CCC/scripts/fastqc_trimmed.sh
#purpose: quality checking of trimmed reads using FASTQC on Pegasus compute node

#BSUB -J AcerCCC_fastqc_trimmed
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -o fastqc_trimmed.out
#BSUB -e fastqc_trimmed.err
#BSUB -n 8
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics/Allyson_CCC"

samples=${and}/Allyson_CCC/trimmed/samples.txt

# trimming the files
for SAMP in `cat ${samples}`
do

echo '#!/bin/bash' > /scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/${SAMP}_trimming.job
echo '#BSUB -q general' >> /scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/${SAMP}_trimming.job
echo '#BSUB -J '$SAMP'' >> /scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/${SAMP}_trimming.job
echo '#BSUB -o /scratch/projects/transcriptomics/allysondemerlis/scripts/'$SAMP'_error_trimming.txt' >> /scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/${SAMP}_trimming.job
echo '#BSUB -e /scratch/projects/transcriptomics/allysondemerlis/scripts/'$SAMP'_output_trimming.txt' >> /scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/${SAMP}_trimming.job

echo 'module load fastqc/0.10.1' >> /scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/${SAMP}_trimming.job

echo 'cd ${and}/trimmed/' >> /scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/${SAMP}_trimming.job

echo 'fastqc *.fastq.gz
--outdir ${and}/trimmed/' >> /scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/${SAMP}_trimming.job

done

