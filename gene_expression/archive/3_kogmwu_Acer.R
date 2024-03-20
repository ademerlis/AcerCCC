#### packages ####

#install.packages("KOGMWU")
library(KOGMWU)
library(tidyverse)

#### loading KOG and gene annotations #### 
gene2kog=read.table("bioinformatics/Acervicornis_iso2kogClass.tab",sep="\t", fill=T) #iso2kogClass.tab not iso2kogClass1.tab because that file has an "error" when you try to view it using the terminal
head(gene2kog)

gene2kog %>% 
  rename(gene = V1, KOG = V2) -> gene2kog_table

read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                        sep = "\t",
                        quote="", fill=FALSE) %>%
  rename(gene = V1,
         annot = V2) -> iso2geneName

#### pairwise treatments (fc) ####

load("RData_files/nursery_vs_CCC_fc.RData") 


fc.Location_nursery_vs_CCC=kog.mwu(Location_nursery_vs_CCC.fc,gene2kog)
#write_csv(fc.Location_nursery_vs_CCC, "KOG_pvalues_fc_nurseryvsCCC.csv")

full_join(Location_nursery_vs_CCC.fc, gene2kog_table, by="gene") %>% 
  drop_na() %>% 
  filter(!KOG == "" & !KOG == "Function Unknown") %>% 
  rename(term = KOG) %>% 
  full_join(., fc.Location_nursery_vs_CCC, by = "term") %>% 
  full_join(., iso2geneName, by = "gene") %>% 
  select(gene, annot, lfc, term, nseqs:padj) %>% 
  rename(KOG = term) %>% 
  drop_na() %>% 
  write_csv("KOGterms_allgenes_pvalues_fc_nurseryvsCCC.csv")

fc.Location_nursery_vs_CCC %>% 
  filter(!term == "" & !term == "Function Unknown") -> fc.Location_nursery_vs_CCC


#each fc table has 23 KOG terms. I think this is a setting within the kog.mwu function

# compiling a table of delta-ranks to compare these results:
ktable=makeDeltaRanksTable(list("nursery_vs_CCC"=fc.Location_nursery_vs_CCC))

library(RColorBrewer)
color = colorRampPalette(rev(c(brewer.pal(n = 7, name ="RdBu"),"royalblue","darkblue")))(100)

# Making a heatmap with hierarchical clustering trees: 
#pdf(file="KOG_Acer_host_fc.pdf", width=7, height=8)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation",color=color, cellwidth=15, cellheight=15) 
while (!is.null(dev.list()))  dev.off()
#needed to manually save this as a PDF

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#scatterplots between pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval) 

