#### PACKAGES ####

# run these once, then comment out
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.15")
# BiocManager::install("DESeq2",dependencies=T)
#BiocManager::install("arrayQualityMetrics",dependencies=T)  # requires Xquartz, xquartz.org
# BiocManager::install("BiocParallel")

# install.packages("pheatmap")
# install.packages("VennDiagram")
# install.packages("gplots")
# install.packages("vegan")
# install.packages("plotrix")
# install.packages("ape")
# install.packages("ggplot2")
# install.packages("rgl")
# install.packages("adegenet")


#### DATA IMPORT ####
# assembling data, running outlier detection, and fitting models
# (skip this section if you don't need to remake models)

library(DESeq2)
library(arrayQualityMetrics)
library(tidyverse)

#read in counts
counts = read.delim("results/counts.txt")

#select only Acropora genes, remove Symbiodinium
Acer_counts <- counts %>% filter(!grepl("^Symbiodinium[0-9]+$", X))
Sym_counts <- counts %>% filter(!grepl("^Acropora_[0-9]+$", X))

column_to_rownames(Acer_counts, var ="X") -> Acer_counts

Acer_counts %>% 
  rename_with(~sub("^X(.+)_trimmed$", "Acer_\\1", .)) %>% 
  select(!X.1) -> Acer_counts
 
# how many genes we have total?
nrow(Acer_counts) #24,453
ncol(Acer_counts) #20 samples

# Based on alignment rates and number of reads (and the one "staghybrid" sample), I will be removing the following samples from downstream analysis: 
#Acer_1087, Acer_1088, Acer_1096, Acer_2264, Acer_2383

Acer_counts %>% 
  select(!c("Acer_1087", "Acer_1088", "Acer_1096", "Acer_2264", "Acer_2383")) -> Acer_counts

nrow(Acer_counts) #24,453
ncol(Acer_counts) #15 samples

# filtering out low-count genes
keep <- rowSums(Acer_counts) >= 10
countData <- Acer_counts[keep,]
nrow(countData) #21814
ncol(countData) #15
#write.csv(countData, file = "Acer_countdata.csv")

# for WCGNA: removing all genes with counts of <10 in more than 90 % of samples
counts4wgcna = Acer_counts[apply(Acer_counts,1,function(x) sum(x<10))<ncol(Acer_counts)*0.9,]
nrow(counts4wgcna) #18192
ncol(counts4wgcna) #15
#write.csv(counts4wgcna, file="Acer_counts4wgcna.csv")

# importing a design .csv file
design = readxl::read_xlsx("sample_metadata.xlsx")

design %>% 
   filter(!Sample_ID %in% c("1087", "1088", "1096", "2264", "2383")) -> design

design %>% 
  mutate(Sample_ID = paste("Acer_", Sample_ID, sep = "")) -> design

column_to_rownames(design, var="Sample_ID") -> design
design$Genotype <- as.factor(design$Genotype)
design$Genotype <- factor(gsub("-", "_", design$Genotype)) #DESeq2 does not like hyphens in factor names
design$Location <- as.factor(design$Location)
design$Genotype <- factor(gsub("'", "", design$Genotype)) #DESeq2 does not like hyphens in factor names
design$Genotype <- factor(gsub(" ", "", design$Genotype)) #DESeq2 does not like hyphens in factor names

str(design)

#### FULL MODEL DESIGN (Genotype + Location) and OUTLIERS ####

#when making dds formula, it is CRITICAL that you put the right order of variables in the design. The design indicates how to model the samples, 
#(here: design = ~batch + condition), 
#that we want to measure the effect of the condition, controlling for batch differences. The two factor variables batch and condition should be columns of coldata.

# make big dataframe including all factors and interaction, getting normalized data for outlier detection
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Location)
#note: I will not be doing an interaction term of Genotype:Location because I do not have enough biological replicates of each genotype per location

# for large datasets, rlog may take too much time, especially for an unfiltered dataframe
# vsd is much faster and still works for outlier detection
Vsd=varianceStabilizingTransformation(dds)

colnames(Vsd)

library(Biobase)
e=ExpressionSet(assay(Vsd), AnnotatedDataFrame(as.data.frame(colData(Vsd))))

# running outlier detection
arrayQualityMetrics(e,intgroup=c("Location", "Genotype"),force=T) 
# open the directory "arrayQualityMetrics report for e" in your working directory and open index.html
# Array metadata and outlier detection overview gives a report of all samples, and which are likely outliers according to the 3 methods tested.
#I typically remove the samples that violate *1 (distance between arrays).
# Figure 2 shows a bar plot of array-to-array distances and an outlier detection threshold based on your samples. 
#Samples above the threshold are considered outliers
# under Figure 3: Principal Components Analyses, look for any points far away from the rest of the sample cluster
# use the array number for removal in the following section

# 3 outliers detected: 5,6,8

# if there were outliers:
outs=c(5,6,8) #these numbers were taken from the index.html report from arrayQualityMetrics Figure 2 "Outlier detection"
countData=countData[,-outs]
Vsd=Vsd[,-outs]
counts4wgcna=counts4wgcna[,-outs]
design=design[-outs,]

# remaking model with outliers removed from dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Location)

# save all these dataframes as an Rdata package so you don't need to rerun each time
save(dds,design,countData,Vsd,counts4wgcna,file="initial_fullddsdesigncountsVsdcountsWGCNA.RData")

# generating normalized variance-stabilized data for PCoA, heatmaps, etc
vsd=assay(Vsd)
colnames(vsd)
# takes the sample IDs and factor levels from the design to create new column names for the dataframe
snames=paste(colnames(countData),design[,2], design[,3],sep="_")
snames #i.e. 

# renames the column names
colnames(vsd)=snames

save(vsd,design,file="vsd.RData")

# more reduced stabilized dataset for WGCNA
wg = DESeqDataSetFromMatrix(countData=counts4wgcna, colData=design, design=~ Genotype + Location)
vsd.wg=assay(varianceStabilizingTransformation(wg), blind=FALSE) #blind=TRUE is the default, and it is a fully unsupervised transformation. However, the creator of DESeq2,
#Michael Love, recommends using blind=FALSE for downstream analyses because when transforming data, the full use of the design information should be made. If many genes have
#large differences in counts due to experimental design, then blind=FALSE will account for that.

head(vsd.wg)
colnames(vsd.wg)=snames
colnames(vsd.wg)
colnames(vsd)
save(vsd.wg,design,file="data4wgcna.RData")

#### DESEQ ####

# with multi-factor, multi-level design
load("Rdata_files/initial_fullddsdesigncountsVsdcountsWGCNA.RData")
library(DESeq2)
library(BiocParallel)

# Running full model for contrast statements
#dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~ Genotype + Location)
rownames(design)
colnames(countData)
dds=DESeq(dds, parallel=TRUE)

# saving all models
save(dds,file="realModels_Acer.RData")

#### DEGs and CONTRASTS ####

load("Rdata_files/realModels_Acer.RData")
library(DESeq2)
library(tidyverse)

read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                          sep = "\t",
                          quote="", fill=FALSE) %>%
  rename(gene = V1,
         annot = V2) -> iso2geneName

resultsNames(dds)
#primary ones of interest: "Location_nursery_vs_CCC"....

# Location_nursery_vs_CCC
Location_nursery_vs_CCC=results(dds,contrast=c("Location","nursery","CCC"))
summary(Location_nursery_vs_CCC, alpha = 0.05)

Location_CCC_vs_nursery= results(dds,contrast=c("Location","CCC","nursery"))
summary(Location_CCC_vs_nursery, alpha = 0.05)

degs_Location_nursery_vs_CCC=row.names(Location_nursery_vs_CCC)[Location_nursery_vs_CCC$padj<0.05 & !(is.na(Location_nursery_vs_CCC$padj))]
length(degs_Location_nursery_vs_CCC) #828

save(Location_nursery_vs_CCC, degs_Location_nursery_vs_CCC, file="pvals.RData")

Location_nursery_vs_CCC %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  filter(padj<0.05) %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  drop_na(baseMean) %>% 
  select(gene, annot, baseMean:padj) %>% 
  write_csv("Location_nursery_vs_CCC_annotDGEs.csv")


#### GO/KOG EXPORT ####

load("RData_files/realModels_Acer.RData")

# fold change (fc) can only be used for binary factors, such as control/treatment, or specific contrasts comparing two factor levels
# log p value (lpv) is for multi-level factors, including binary factors

# Untreated vs Initial
# log2 fold changes:
source=Location_nursery_vs_CCC[!is.na(Location_nursery_vs_CCC$padj),]
Location_nursery_vs_CCC.fc=data.frame("gene"=row.names(source))
Location_nursery_vs_CCC.fc$lfc=source[,"log2FoldChange"]
head(Location_nursery_vs_CCC.fc)
write.csv(Location_nursery_vs_CCC.fc,file="nursery_vs_CCC_fc.csv",row.names=F,quote=F)
save(Location_nursery_vs_CCC.fc,file="Rdata_files/nursery_vs_CCC_fc.RData")

# signed log FDR-adjusted p-values: -log(p-adj)* direction:
nursery_vs_CCC.p=data.frame("gene"=row.names(source))
nursery_vs_CCC.p$lpv=-log(source[,"padj"],10)
nursery_vs_CCC.p$lpv[source$stat<0]=nursery_vs_CCC.p$lpv[source$stat<0]*-1
head(nursery_vs_CCC.p)
write.csv(nursery_vs_CCC.p,file="nursery_vs_CCC_lpv.csv",row.names=F,quote=F)
save(nursery_vs_CCC.p,file="Rdata_files/nursery_vs_CCC_lpv.RData")


#### ANNOTATING DGES ####
load("RData_files/realModels_Acer.RData")

library(tidyverse)

#Untreated vs. Initial
Location_nursery_vs_CCC %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  mutate(lpv = -log(padj, base = 10)) %>%
  mutate(lpv = if_else(stat < 0, lpv * -1, lpv)) %>% 
  filter(abs(lpv) >= 1.3) %>% 
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = c("gene" = "gene")) %>% write_csv("Location_nursery_vs_CCC_annotatedDGEs.csv") %>% str()
#830 genes 


#### PCOA and PERMANOVA ####

# heatmap and hierarchical clustering:
load("Rdata_files/vsd.RData")
library(pheatmap)
# similarity among samples
pdf(file="heatmap_fullmodel.pdf", width=15, height=15)
head(vsd)
pheatmap(cor(vsd))
dev.off()

# Principal coordinates analysis
library(vegan)
# library(rgl)
library(ape)

conditions=design
colnames(vsd)
rownames(conditions)
#while the names aren't equal, they are in the same order, so when doing the PCoA and using conditions it still should work

# creating a PCoA eigenvalue matrix
dds.pcoa=pcoa(dist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors
# copy this table for % variation explained by each axis (Relative_eig column)
dds.pcoa$values
#axis 1 = 39.05%
#axis 2 = 11.88%

# how many good PC's do we have? Compared to random ("broken stick") model
# plotting PCoA eigenvalues 
pdf(file="PCoA_Manhattan.pdf", width=6, height=6)
plot(dds.pcoa$values$Relative_eig)
points(dds.pcoa$values$Broken_stick,col="red",pch=3)
dev.off()
# the number of black points above the line of red crosses (random model) corresponds to the number of good PC's
#there is 1 "good PC" based on this figure

# plotting PCoA by treatment and Genotype (axes 1 and 2)
par(mfrow=c(1,2))
plot(scores[,1], scores[,2],col=c("red","blue")[as.numeric(as.factor(conditions$Location))],pch=c(15,17,25)[as.numeric((as.factor(conditions$Genotype)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Location")
ordispider(scores, conditions$Location, label=F, col=c("red","blue"))
legend("topright", legend=c("CCC", "Nursery"), fill = c("red","blue"), bty="n")
legend("topleft", legend=c("Cheetos_B", "MiamiBeach_C", "SunnyIsles_E"), pch=c(15,17,25), bty="n")
plot(scores[,1], scores[,2],col=c("orange","lightblue", "pink")[as.numeric(as.factor(conditions$Genotype))],pch=c(15,17,25)[as.numeric((as.factor(conditions$Location)))], xlab="Coordinate 1", ylab="Coordinate 2", main="Genotype")
ordispider(scores, conditions$Genotype, label=F, col=c("orange","darkblue", "magenta"))
legend("topleft", legend=c("Cheetos_B", "MiamiBeach_C", "SunnyIsles_E"), fill = c("orange","darkblue", "magenta"), bty="n")
legend("topright", legend=c("CCC", "Nursery"), pch=c(15,17,25), bty="n")
dev.off()
#manually save as pdf

# plotting PCoA by location
plot(scores[,1], scores[,2],col=c("red","blue")[as.numeric(as.factor(conditions$Location))], xlab="Coordinate 1", ylab="Coordinate 2", main="Location")
ordispider(scores, conditions$Location, label=F, col=c("red","blue"))
legend("topleft", legend=c("CCC", "Nursery"), fill = c("red","blue"), bty="n")
dev.off()
#manually save as pdf


# formal analysis of variance in distance matricies: 
ad=adonis2(t(vsd)~Genotype + Location,data=design,method="manhattan",permutations=1e6)
ad
summary(ad)
as.data.frame(ad) %>% write_csv("permanova_results.csv")

