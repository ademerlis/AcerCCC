#!/bin/bash
#BSUB -J trim_CCC_parallel_trimgalore
#BSUB -q general
#BSUB -P and_transcriptomics
#BSUB -n 16
#BSUB -W 120:00
#BSUB -o trim_CCC_parallel_trimgalore%J.out
#BSUB -e trim_CCC_parallel_trimgalore%J.err
#BSUB -u allyson.demerlis@earth.miami.edu
#BSUB -N

and="/scratch/projects/and_transcriptomics"

for sample in ${and}/Allyson_CCC/fastq_files/*.gz ;

do
echo '#!/bin/bash' > /scratch/projects/and_transcriptomics/Allyson_CCC/scripts/${sample}_trimming.job
echo '#BSUB -q general' >> /scratch/projects/and_transcriptomics/Allyson_CCC/scripts/${sample}_trimming.job
echo '#BSUB -J '$sample'' >> /scratch/projects/and_transcriptomics/Allyson_CCC/scripts/${sample}_trimming.job
echo '#BSUB -o /scratch/projects/and_transcriptomics/Allyson_CCC/scripts/'$sample'_error_trimming.txt' >> /scratch/projects/and_transcriptomics/Allyson_CCC/scripts/${sample}_trimming.job
echo '#BSUB -e /scratch/projects/and_transcriptomics/Allyson_CCC/scripts/'$sample'_output_trimming.txt' >>
/scratch/projects/and_transcriptomics/Allyson_CCC/scripts/${sample}_trimming.job
echo '/scratch/projects/and_transcriptomics/programs/TrimGalore-0.6.10/trim_galore ${sample} \
--gzip \
--fastqc \
--fastqc_args "--outdir /scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/" \
--illumina \
--cores 4 \
--three_prime_clip_R1 12 \
--nextseq 30 \
--length 20 \
-o /scratch/projects/and_transcriptomics/Allyson_CCC/trimmed/ ; \' >>
/scratch/projects/and_transcriptomics/Allyson_CCC/scripts/${sample}_trimming.job

done
