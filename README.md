This repository contains data and analysis scripts for the manuscript:

## Characterization of gene expression patterns in *Acropora cervicornis* following transplantation to a marginal environment
#### **Authors:** Allyson DeMerlis1,2,3, Colin Foord4, Dalton J. Hesley1, Richard F. Karp1,2,3, Lys M. Isma1, Natalia Andrade-Rodriguez1, Kevin Wong1, Michael S. Studivan2,3, Andrew C. Baker1, Diego Lirman1, Nikki Traylor-Knowles1, Ian C. Enochs3
1. Rosenstiel School for Marine, Atmospheric, and Earth Science, University of Miami, Miami, FL, USA
2. University of Miami Cooperative Institute for Marine and Atmospheric Studies, Miami, FL, USA
3. Atlantic Oceanographic and Meteorological Laboratory, Ocean Chemistry and Ecosystems Division, U.S. National Oceanic and Atmospheric Administration, Miami, FL, USA
4. Coral Morphologic, Miami, FL, USA

#### **Journal:** 

-----

### Description:
These repository contains all data and code used to study the physiological impact of translocation of _Acropora cervicornis_ to an urbanized environment in the Port of Miami.

### Contents:

#### environmental_data:
* This folder contains 3 Rmarkdown files:
    * **1_tidying_env_data.Rmd**: tidying temperature data
    * **2_thermalvariability.Rmd**: calculating and plotting thermal variability
    * **3_MMMDHW.Rmd**: obtaining the mean monthly maximum (MMM) data from NOAA Coral Watch.
* In the subfolder **raw_data**, there are .csv files for the temperature data at the Key Biscayne (KB) Nursery (KBNursery_X_.csv), for the [CCC site](https://github.com/ademerlis/AcerCCC/blob/main/environmental_data/raw_data/Oct2020_Jan2022_1901102_urban-2020-10_(0)_Temperature.xlsx), and the [NOAA Coral Watch 5km climatology](https://github.com/ademerlis/AcerCCC/blob/main/environmental_data/raw_data/ct5km_climatology_v3.1.nc).
* In the subfolder **results_csv**, there are tables and plots for the temperature data comparisons between sites.

#### gene_expression:
* The raw sequence .fastq files will be made publicly available on the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) under BioProject PRJNA1333470, upon publication.
* The file **readcounts_rawtrimmedaligned_summary.csv** has a description of the raw, trimmed, and aligned reads for *A. cervicornis*.
* The file **sample_metadata.csv** has the sample metadata, including date sampled, genotype, and location.
* Subfolders:
    * **1_bioinformatics** subfolder contains the scripts needed to process raw 3' RNA-Seq sequences on the UM HPC, [Pegasus](https://acs-docs.readthedocs.io/pegasus/README.html), which uses an LSF resource manager. The bioinformatics pipeline outlined in this folder includes: FastQC > Cutadapt > bowtie2 > samtools
    * **2_OutlierDetection_Acer** subfolder contains results of *ArrayQualityMetrics* for filtering outliers before DESeq2 analysis for the _A. cervicornis_ host genes.
    * **2_OutlierDetection_Sym** subfolder contains results of *ArrayQualityMetrics* for filtering outliers before DESeq2 analysis for the symbiont genes (_S. fitti_).
    * **3_WGCNA** subfolder contains results following the weighted gene co-expression network analysis (WGCNA)
* R code:
    * **1_Acer_deseq2.R**: R file for analyzing *A. cervicornis* host differential gene expression.
    * **1_Sym_deseq2.R**: R file for analyzing *S. fitti* symbiont differential gene expression.
    * **2_Acer_PCA.Rmd** For making the principal components analysis (PCA) plots.
    * **3_VolcanoPlots.Rmd** Code for making volcano plots.
    * **4_specific_gene_expression.Rmd** Code for making boxplots of specific gene expression patterns.
    * **5_wgcna_Acer.R** Code for running weighted gene co-expression network analysis (WGCNA).
    * **6_wgcna_GO.Rmd** Code for running TopGO on significant WGCNA modules.
 
#### Figures:
* This folder contains all the manuscript figures.
</br>

#### Notes
* This README.md was formatted following [Dr. Ana Palacio's GitHub repository](https://github.com/anampc/Acer_NH4_disease/tree/master).
