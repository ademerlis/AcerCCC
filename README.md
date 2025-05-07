This repository contains data and analysis scripts for the manuscript:

## 
#### **Authors:** 
#### **Journal:** 

-----

### Description:
These repository contains all data and code used to study the physiological impact of translocation of _Acropora cervicorns_ to an urbanized environment in the Port of Miami.

### Contents:

#### Environmental Data:
* In the subfolder **raw_data**, 

#### Physiology:
* The complete data file that contains all metadata and buoyant weight mesaurements is the file **metadata.csv**.
* Physiology data analysis is broken up into three subfolders: **Calcification**, **R_intensity**, and **Photosynthetic Efficiency**.
* In the subfolder **Calcification**, **calcification.Rmd** has all the code needed for Supplementary Figure 1 and Table S3 statistics.
* In the subfolder **R_intensity**, **colorscores.Rmd** has all the code needed for Supplementary Figure 2 and Table S4 statistics.
* The subfolder **Photosynthetic Efficiency** has several R files that are important for analysis.
  - **1_importtidyPAMdata.Rmd** includes the function to import raw IPAM data into R, and creates a format file that matches coral fragment metadata to the area of interest (AOI) and YII (photosynthetic efficiency) values to the IPAM image metadata files. The raw files for this can be found in the **ipam_data** subfolder.
 
#### Gene Expression:
* The raw sequence .fastq files will be made publicly available on the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) under BioProject PRJNA1196005, upon publication.
* The file **readcounts_rawtrimmedaligned_summary.csv** has a description of the raw, trimmed, and aligned reads for *A. cervicornis* and *P. clivosa*.
* The file **treatment_metadata.csv** has the sample metadata, including date sampled, treatment group, and genotype.
* For each species, there is a corresponding subfolder.
* In **Acervicornis**:
  - **1_bioinformatics** subfolder contains the scripts needed to process raw 3' RNA-Seq sequences on the UM HPC, [Pegasus](https://acs-docs.readthedocs.io/pegasus/README.html), which uses an LSF resource manager. The bioinformatics pipeline outlined in this folder includes: FastQC > Cutadapt > bowtie2 > samtools
  - **2_DESeq2_host.Rmd** R markdown file for analyzing *A. cervicornis* host differential gene expression as a result of the variable temperature treatment. This has the code needed for Figures 3 and 4, as well as Supplementary Tables 8, 9, and 11.
  - **2_DESeq2_symbiont.Rmd** R markdown file for analyzing *A. cervicornis* symbiont differential gene expression as a result of the variable temperature treatment. This has the code needed for Figures 3 and 4, as well as Supplementary Tables 8, 9, and 12.
  - **2_outlier_detection** subfolder for the *ArrayQualityMetrics* outlier detection method used for removing samples based on gene expression patterns.
  - **3_GO-MWU** subfolder contains the R markdown files and functions necessary for running the specific method of Gene Ontology (GO) enrichment analysis used in this study. This has the code needed for Figures 5 and 6, as well as Supplementary Tables 10, 15, 16, and 17.
  - **results_csv** subfolder contains all the results files for the *A. cervicornis* host and symbiont differential gene expression.
 
#### Figures:
* This folder contains all the manuscript figures.
</br>

#### Notes
* This README.md was formatted following [Dr. Ana Palacio's GitHub repository](https://github.com/anampc/Acer_NH4_disease/tree/master).
