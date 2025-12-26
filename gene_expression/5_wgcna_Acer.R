#### PACKAGES ####

# installing WGCNA:
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.15")
# BiocManager::install(c("AnnotationDbi", "impute", "GO.db", "preprocessCore"))
# install.packages("flashClust")
# install.packages("WGCNA",dependencies=TRUE)
# repos="http://cran.us.r-project.org"
# run these above commands once, then comment out

# note: all WGCNA tutorials have moved to dropbox: https://www.dropbox.com/scl/fo/4vqfiysan6rlurfo2pbnk/h?rlkey=thqg8wlpdn4spu3ihjuc1kmlu&e=1&dl=0

# always run these before running any of the following script chunks
library(tidyverse)
library(WGCNA)
library(flashClust)
library(ape)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()
#BiocManager::install(c("GO.db", "preprocessCore", "impute") ) #WGCNA dependencies

#### DATA IMPORT and TRAITS ####

# importing data generated from DESeq2 script
lnames=load("RData_files/data4wgcna.RData")
lnames # "vsd.wg"  "design" # log-transformed variance-stabilized gene expression, and table or experimental conditions
datt=t(vsd.wg)
ncol(datt) #18192
nrow(datt) #12

head(design)
str(design)
colnames(vsd.wg)

all.equal(colnames(vsd.wg), rownames(design)) #TRUE

#change Location to be binary (0 = FALSE, 1 = TRUE)
CCC = as.numeric(design$Location=="CCC")
Nursery = as.numeric(design$Location == "nursery")

traits <- data.frame(cbind(CCC, Nursery))

#### OUTLIER DETECTION ####

# identifies outlier genes
gsg = goodSamplesGenes(datt, verbose = 3);
gsg$allOK #if TRUE, no outlier genes
#TRUE!

# calculates mean expression per array, then the number of missing values per array
meanExpressionByArray=apply(datt,1,mean, na.rm=T)
NumberMissingByArray=apply( is.na(data.frame(datt)),1, sum)
NumberMissingByArray
# keep samples with missing values under 500
# in this case, all samples OK

# plots mean expression across all samples
quartz()
barplot(meanExpressionByArray,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:12), cex.names = 0.7)
# look for any obvious deviations in expression across samples

# sample dendrogram and trait heat map showing outliers
A=adjacency(t(datt),type="signed")                 #SELECT SIGNED OR UNSIGNED HERE
#I am going with signed because if you pick unsigned, it mixes negatively and positiviely correlated nodes together and 
#direction of correlation does matter for downstream analysis. 
#(full explanation here: https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/)

# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(traits,signed=TRUE))
dimnames(traitColors)[[2]]=paste(names(traits))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
quartz()
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")
# the resulting plot shows a sample dendrogram and the spread of your traits across the cluster
# outlier samples will show as red in the outlierC row

# Remove outlying samples from expression and trait data
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datt=datt[!remove.samples,]
traits=traits[!remove.samples,] #1 sample removed

str(datt) #11
str(traits) #11
  
write.csv(traits, file="traits.csv")

save(datt,traits,file="wgcnaData.RData")


#### SOFT THRESHOLDS ####

library(tidyverse)
library(WGCNA)
library(flashClust)
library(ape)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

load("wgcnaData.RData")

#datt=t(vsd.wg)
str(datt) #11

# Try different betas ("soft threshold") - power factor for calling connections between genes
powers = c(seq(from = 2, to=26, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(datt, powerVector = powers, verbose = 8,networkType="signed")

# Plot the results:
# Run from the line below to dev.off()
sizeGrWindow(9, 5)
quartz()
pdf("soft_threshold_signed.pdf",height=4, width=8)

par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#### MAKING MODULES ####

# take a look at the threshold plots produced above, and the output table from the pickSoftThreshold command
# so Michael's code recommends this: "pick the power that corresponds with a SFT.R.sq value above 0.90."

#However, in doing some googling (because none of my SFTs in this dataset exceed 0.9), I came across this 
#helpful tutorial: https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#45_Format_normalized_data_for_WGCNA
# it says this: WGCNA’s authors recommend using a power that has an signed R2
# above 0.80, otherwise they warn your results may be too noisy to be meaningful.
# If you have multiple power values with signed R2
# above 0.80, then picking the one at an inflection point, in other words where the R2
# values seem to have reached their saturation (Zhang and Horvath 2005). 
#You want to a power that gives you a big enough R2 but is not excessively large.

# they also provide this code to visualize the sft cut-off better

sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

#based on this graph above, I'm going to pick #14 ! 
#14 led to NA for the MEgrey module, so i'm going to try to run 15 instead and see if that helps

# ok after all of this back and forth, it seems like the "lower" sft numbers like 14 and 15 
# result in only 2 genes for the grey module, which ends up causing problems for the 
# correlation of the MEs. Doesn't matter if i use the traits_n_temp or just traits.
#so i think i will just stick to 26 because that gives the grey module 7 genes and seems to work
# for the correlation plot. 

# run from the line below to the save command
s.th=26 # re-specify according to previous section
adjacency = adjacency(datt, power = s.th,type="signed");
TOM = TOMsimilarity(adjacency,TOMType="signed");
rm(adjacency) #for memory space
dissTOM = 1-TOM
rm(TOM)
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average")

plot(geneTree, xlab="", sub="", main="Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 100; 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Calculate eigengenes
MEList = moduleEigengenes(datt, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = flashClust(as.dist(MEDiss), method = "average");

save(dynamicMods,dynamicColors,MEs,METree,geneTree,file="1stPassModules.RData")


#### MERGING MODULES ####

mm=load('RData_files/1stPassModules.RData')
mm
lnames=load('RData_files/wgcnaData.RData')
traits
head(datt)

quartz()

MEDissThres = 0.6 # in the first pass, set this to 0 - no merging (we want to see the module-traits heatmap first, then decide which modules are telling us the same story and better be merged)
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")  # on 2nd pass: does this cut height meet your merging goals? If not, reset MEDissThres and replot

# Call an automatic merging function
merge = mergeCloseModules(datt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# plotting the fabulous ridiculogram
quartz()
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = FALSE, guideHang = 0.05,lwd=0.3)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Calculate dissimilarity of module eigengenes
quartz()

MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

# how many genes in each module?
table(moduleColors)
# Save module colors and labels for use in subsequent parts
save(MEs, geneTree, moduleLabels, moduleColors, file = "networkdata_signed.RData")


#### MODULE CORRELATIONS ####
# plotting correlations with traits:
load(file = "RData_files/networkdata_signed.RData")
load(file = "RData_files/wgcnaData.RData");
traits <- read_csv("4_WGCNA/traits.csv")

# Define numbers of genes and samples
nGenes = ncol(datt); #18192
nSamples = nrow(datt); #11
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# module-trait correlations

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(traits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)

table(moduleColors) # gives numbers of genes in each module

# shows only significant correlations
quartz()
library(RColorBrewer)
modLabels=sub("ME","",names(MEs))

ps=signif(moduleTraitPvalue,1)
cors=signif(moduleTraitCor,2)
textMatrix = cors;
# paste(cors, "\n(",ps, ")", sep = "");
textMatrix[ps>0.05]="-"
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               ySymbols = modLabels,
               yLabels = modLabels,
               colorLabels = FALSE,
               colors = colorRampPalette(c("blue","lightblue","white","coral","red"))(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-0.7,0.7),
               main = paste("A.cervicornis Module-Trait correlations"))

# module size barplot
quartz()
labelShift=750 # increase to move module size labels to the right
par(mar = c(6, 8.5, 3, 3));
mct=table(moduleColors)
mct[modLabels]
x=barplot(mct[rev(modLabels)],horiz=T,las=1,xlim=c(0,16000),col=rev(modLabels))
text(mct[rev(modLabels)]+labelShift,y=x,mct[rev(modLabels)],cex=0.9) 

# If it was first pass with no module merging, this is where you examine your heatmap and dendrogram of module eigengenes 
#to see where you would like to set cut height (MEDissThres parameter) 
#in the previous section to merge modules that are telling the same story for your trait data 
# A good way to do it is to find a group of similar modules in the heat map and then see at which tree height they connect in the dendrogram

#### GO BACK AND MERGE ####
#done!

#### MODULE MEMBERSHIP SCATTERPLOTS ####

# scatterplots of gene significance (correlation-based) vs kME
load(file = "RData_files/networkdata_signed.RData")
load(file = "RData_files/wgcnaData.RData");
traits
table(moduleColors)

# run for each of these statements individually
whichTrait="CCC"

nGenes = ncol(datt);
nSamples = nrow(datt);
selTrait = as.data.frame(traits[,whichTrait]);
names(selTrait) = whichTrait
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(signedKME(datt, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");

rownames_to_column(geneTraitSignificance, var="Locus") -> geneTraitSignificance
rownames_to_column(geneModuleMembership, var="Locus") -> geneModuleMembership

merged_GS_MM <- inner_join(geneTraitSignificance, geneModuleMembership, by = "Locus")

merged_GS_MM %>% 
  tidyr::gather(module, module_membership, MMdarkturquoise:MMgrey) %>%
  filter(module== "MMdarkmagenta" | module == "MMmediumpurple3") %>% 
  ggplot(aes(module_membership, GS.CCC, color = module)) +
  geom_point(size = 2, alpha = 1/100) +
  geom_smooth(color="black", method = lm, se = T, fill = "grey") +
  theme_classic() +
  ylim(c(-1,1)) +
  xlim(c(-1,1)) +
  facet_wrap(~module, scales = "free") +
  scale_color_manual(values = c("darkmagenta", "mediumpurple")) +
  theme(legend.position = "none") +
  labs(x="Module Membership",
       y="Gene Significance for Coral City Camera")
ggsave("MM_GS_CCC.png", width = 8, height = 5)


# selecting specific modules to plot (change depending on which trait you're looking at)
moduleCols=c("darkmagenta","mediumpurple3")
quartz()
# set par to be big enough for all significant module correlations, 
#then run the next whichTrait and moduleCols statements above and repeat from the 'for' loop
par(mfrow=c(1,6)) 

counter=0
# shows correlations for all modules
for(module in modNames[1:length(modNames)]){
counter=counter+1}

quartz()
	par(mfrow=c(3,3))

# shows correlations for significant modules only as specified above
for (module in moduleCols) {
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
column = match(module, modNames);
moduleGenes = moduleColors==module;
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste(module,"module membership"),
ylab = paste("GS for", whichTrait),
col = "grey50",mgp=c(2.3,1,0))
}
	
load("Rdata_files/initial_fullddsdesigncountsVsdcounts.RData")

design %>% 
  rownames_to_column(var="Sample_ID") %>% 
  filter(!Sample_ID == "Acer_1099") -> design

#plotting eigengene expression
 MEs %>%
	  tibble::rownames_to_column("accession_code") %>%
	  dplyr::inner_join(design %>%
	                      dplyr::select(Sample_ID, Location),
	                    by = c("accession_code" = "Sample_ID")) %>% 
   select(accession_code, MEdarkmagenta, MEmediumpurple3, Location) %>% 
   pivot_longer(MEdarkmagenta:MEmediumpurple3, names_to = "module", values_to = "collective_expression_level") %>% 
	ggplot(.,aes(x = Location, y = collective_expression_level,fill = Location)) +
	  geom_boxplot(width = 0.2) +
   facet_wrap(~module, ncol = 1) +
	  theme_classic() +
   scale_fill_manual(values = c("darkblue", "orange"))
#ggsave("ME_expression_location.pdf")
 
 # stats
 MEs %>%
   tibble::rownames_to_column("accession_code") %>%
   dplyr::inner_join(design %>%
                       dplyr::select(Sample_ID, Location),
                     by = c("accession_code" = "Sample_ID")) %>% 
   select(accession_code, MEdarkmagenta, MEmediumpurple3, Location) -> ME_stats
 
 # Darkmagenta module
 t_test_dark <- t.test(ME_stats$MEdarkmagenta ~ ME_stats$Location)
 
 # Mediumpurple3 module  
 t_test_purple <- t.test(ME_stats$MEmediumpurple3 ~ ME_stats$Location)
 

 design %>% 
   dplyr::rename(site = Location) -> design
 
#  temp_data %>% 
#    mutate(site = case_when(site == "Nursery" ~ "nursery",
#                            site == "CCC" ~ "CCC")) -> temp_data
#  
# full_join(design, temp_data, by = "site") -> sample_data_temp
# 
# MEs %>%
#   tibble::rownames_to_column("accession_code") %>%
#   dplyr::inner_join(sample_data_temp %>%
#                       dplyr::select(Sample_ID, site, Max:DHW),
#                     by = c("accession_code" = "Sample_ID")) %>% 
#   select(accession_code, MEbrown2, MEmagenta3, site:DHW) %>% 
#   pivot_longer(MEbrown2:MEmagenta3, names_to = "module", values_to = "collective_expression_level") %>% 
#   ggplot(.,aes(x = Seasonal, y = collective_expression_level,fill = site)) +
#   geom_boxplot(width = 0.2) +
#   facet_wrap(~module) +
#   theme_classic() +
#   scale_fill_manual(values = c("orange", "darkblue"))

#### EIGENGENE SANITY CHECK ####

# eigengene-heatmap plot (sanity check - is the whole module driven by just one crazy sample?)
# note: this part does not make much sense for unsigned modules
load(file = "Rdata_files/networkdata_signed.RData")
load(file = "Rdata_files/wgcnaData.RData");

# run for each of these statements individually

#significant modules
which.module="darkmagenta" 
which.module="mediumpurple3"

#non-significant modules
which.module="darkturquoise"
which.module="red"
which.module="orangered4"
which.module="plum1"
which.module="black"
which.module="brown"
which.module="darkolivegreen"
which.module="grey"


datME=MEs
datExpr=datt
quartz()
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
nrgcols=30,rlabels=F,rcols=which.module,
main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
ylab="eigengene expression",xlab="sample")

length(datExpr[1,moduleColors==which.module ]) # darkmagenta = 658 genes

as.data.frame(datExpr[1,moduleColors==which.module ])  %>% 
  rownames_to_column(var="gene") %>%
  write_csv("darkmagenta_genelist.csv")

length(datExpr[1,moduleColors==which.module ]) # mediumpurple3 = 751 genes

as.data.frame(datExpr[1,moduleColors==which.module ])  %>% 
  rownames_to_column(var="gene") %>% 
  write_csv("mediumpurple3_genelist.csv")

length(datExpr[1,moduleColors==which.module ]) # darkturquoise = 323 genes 

as.data.frame(datExpr[1,moduleColors==which.module ])  %>% 
  rownames_to_column(var="gene") %>% 
  write_csv("darkturquoise_genelist.csv")

length(datExpr[1,moduleColors==which.module ]) # red = 4165 genes 

as.data.frame(datExpr[1,moduleColors==which.module ])  %>% 
  rownames_to_column(var="gene") %>% 
  write_csv("red_genelist.csv")

length(datExpr[1,moduleColors==which.module ]) # orangered4 = 400 genes 

as.data.frame(datExpr[1,moduleColors==which.module ])  %>% 
  rownames_to_column(var="gene") %>% 
  write_csv("orangered4_genelist.csv")

length(datExpr[1,moduleColors==which.module ]) # plum1 = 184 genes 

as.data.frame(datExpr[1,moduleColors==which.module ])  %>% 
  rownames_to_column(var="gene") %>% 
  write_csv("plum1_genelist.csv")

length(datExpr[1,moduleColors==which.module ]) # black = 1444 genes 

as.data.frame(datExpr[1,moduleColors==which.module ])  %>% 
  rownames_to_column(var="gene") %>% 
  write_csv("black_genelist.csv")

length(datExpr[1,moduleColors==which.module ]) # brown = 7312 genes 

as.data.frame(datExpr[1,moduleColors==which.module ])  %>% 
  rownames_to_column(var="gene") %>% 
  write_csv("brown_genelist.csv")

length(datExpr[1,moduleColors==which.module ]) # darkolivegreen = 2759 genes 

as.data.frame(datExpr[1,moduleColors==which.module ])  %>% 
  rownames_to_column(var="gene") %>% 
  write_csv("darkolivegreen_genelist.csv")

length(datExpr[1,moduleColors==which.module ]) # grey = 196 genes 

as.data.frame(datExpr[1,moduleColors==which.module ])  %>% 
  rownames_to_column(var="gene") %>% 
  write_csv("grey_genelist.csv")


# If individual samples appear to be driving expression of significant modules, they are likely outliers

#### GO/KOG EXPORT ####

# saving selected modules for GO and KOG analysis (two-parts: Fisher test, MWU test within-module)
library(WGCNA)
load(file = "RData_files/networkdata_signed.RData") # moduleColors, MEs
load(file = "RData_files/wgcnaData.RData") # vsd table
load(file = "RData_files/data4wgcna.RData") # vsd table

# calculating module memberships for all genes for all modules
allkME =as.data.frame(signedKME(datt, MEs)) 
names(allkME)=gsub("kME","",names(allkME))

# run for each of these statements individually
which.module="darkmagenta"
which.module="mediumpurple3"


# Saving data for Fisher-MWU combo test (GO_MWU)
inModuleBinary=as.numeric(moduleColors==which.module)
combo=data.frame("gene"=row.names(t(datt)),"Fish_kME"=allkME[,which.module]*inModuleBinary)
write.csv(combo,file=paste(which.module,".csv",sep=""),row.names=F,quote=F)


#### HEATMAPS ####

# plotting heatmap for named top-kME genes
library(WGCNA)
library(pheatmap)

load(file = "RData_files/networkdata_signed.RData")
load(file = "RData_files/data4wgcna.RData") 
load(file = "RData_files/wgcnaData.RData");

allkME =signedKME(datt, MEs)
gg=read.delim(file="bioinformatics/Acervicornis_iso2geneName.tab",sep="\t")

which.module="darkmagenta"
which.module="mediumpurple3"
design %>% 
  rownames_to_column(var = "Sample_ID") -> design

design %>% 
  filter(!Sample_ID == "Acer_1099") -> design

snames=paste(rownames(datt),design[,3],sep="_")
rownames(datt)
vsd.wg = t(datt)

top=30 # number of named top-kME genes to plot

datME=MEs
datExpr=datt
modcol=paste("kME",which.module,sep="")
sorted=vsd.wg[order(allkME[,modcol],decreasing=T),]
# selection top N names genes, attaching gene names
gnames=c();counts=0;hubs=c()
for(i in 1:length(sorted[,1])) {
	if (row.names(sorted)[i] %in% gg[,1]) { 
		counts=counts+1
		gn=gg[gg[,1]==row.names(sorted)[i],2]
		gn=paste(gn,row.names(sorted)[i],sep=".")
		if (gn %in% gnames) {
			gn=paste(gn,counts,sep=".")
		}
		gnames=append(gnames,gn) 
		hubs=data.frame(rbind(hubs,sorted[i,]))
		if (counts==top) {break}
	}
} 
row.names(hubs)=gnames

str(gnames)
str(hubs)

colnames(hubs) = snames

rownames(hubs)

categories <- c("CCC", "nursery")

# Extract and sort columns for each category
category_columns <- lapply(categories, function(cat) {
  matching_columns <- grep(cat, names(hubs), value = TRUE)
  hubs[, matching_columns, drop = FALSE]
})

# Bind the columns back together in the desired order
reordered_df <- do.call(cbind, category_columns)

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#contrasting2 = colorRampPalette(rev(c("chocolate1","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
#contrasting3 = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)

#pdf(file="heatmap_top30_turquoise.pdf")
pheatmap(as.matrix(reordered_df),scale="row",col=contrasting,border_color=NA,treeheight_col=0,cex=0.9,cluster_rows = F, cluster_cols = F) 
dev.off()

#### HUB GENES ####

library(WGCNA)
library(tidyverse)
load(file = "networkdata_signed.RData")
load(file = "data4wgcna.RData") 
load(file = "wgcnaData.RData");
allkME =as.data.frame(signedKME(datt, MEs))

hubgenes <- chooseTopHubInEachModule(datt, moduleColors, omitColors = "grey", 
                                     power = 2, 
                                     type = "signed")
hubgenes <-data.frame(hubgenes)
hubgenes <- tibble::rownames_to_column(hubgenes, "module")
hubgenes

hubgenes %>%
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") -> hubgenes
hubgenes

write.csv(hubgenes, file="hubgenes.csv")


