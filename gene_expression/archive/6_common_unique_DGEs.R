#### PACKAGES ####

library(tidyverse)

CCC_vs_nursery_annotatedDGEs <- read_csv("results/Location_CCC_vs_nursery_annotDGEs_padj05.csv")

CCC_vs_nursery_annotatedDGEs %>% 
  column_to_rownames(var="gene") -> CCC_vs_nursery_annotatedDGEs

degs_CCCnursery <- rownames(CCC_vs_nursery_annotatedDGEs)

Acer_stresshardening_TvU_annotDGEs <- read_csv("../../Ch2_temperaturevariability2023/gene_expression/Acervicornis/results_csv/TreatedvUntreated_annotDGEs.csv")

Acer_stresshardening_TvU_annotDGEs %>% 
  column_to_rownames(var="gene") -> Acer_stresshardening_TvU_annotDGEs

degs_Acerstresshardening <- rownames(Acer_stresshardening_TvU_annotDGEs)

pairwise=list("Acer Stress-Hardening: Treated vs. Untreated"=degs_Acerstresshardening, 
              "CCC vs. Nursery"=degs_CCCnursery)

library(ggvenn)

ggvenn(pairwise) + 
  scale_fill_manual(values = c("#ca0020", "orange"))
#ggsave("venndiagram_ch1vsch2.pdf")

#common genes

right_join(Acer_stresshardening_TvU_annotDGEs, CCC_vs_nursery_annotatedDGEs)

Acer_stresshardening_TvU_annotDGEs %>% 
  rownames_to_column(var = "gene") -> Acer_stresshardening_TvU_annotDGEs

CCC_vs_nursery_annotatedDGEs %>% 
  rownames_to_column(var = "gene") -> CCC_vs_nursery_annotatedDGEs

Acer_stresshardening_TvU_annotDGEs %>% 
  select(gene, annot, log2FoldChange, padj) %>% 
  rename(log2FoldChange_SH = log2FoldChange, padj_SH = padj) -> Acer_stresshardening_TvU_annotDGEs

CCC_vs_nursery_annotatedDGEs %>% 
  select(gene, annot, log2FoldChange, padj) %>% 
  rename(log2FoldChange_CCCN = log2FoldChange, padj_CCCN = padj) -> CCC_vs_nursery_annotatedDGEs

commongenes<-inner_join(Acer_stresshardening_TvU_annotDGEs, CCC_vs_nursery_annotatedDGEs)
#write_csv("common_genes.csv")


negative_rows_both <- sum(commongenes$log2FoldChange_SH < 0 & commongenes$log2FoldChange_CCCN < 0)
# 8 genes where both are downregulated

positive_rows_both <- sum(commongenes$log2FoldChange_SH > 0 & commongenes$log2FoldChange_CCCN > 0)
#19 genes where both are upregulated

negative_rows_SH <- sum(commongenes$log2FoldChange_SH < 0)
# 44 genes
positive_rows_SH <- sum(commongenes$log2FoldChange_SH > 0)
# 30 genes

negative_rows_CCCN <- sum(commongenes$log2FoldChange_CCCN < 0)
# 19 genes
positive_rows_CCCN <- sum(commongenes$log2FoldChange_CCCN > 0)
# 55 genes

common_lfc<-readxl::read_xlsx("common_genes_LFC.xlsx")

ggplot(data = common_lfc, aes(x = Experiment, y = genes, fill= Regulation))+
  geom_col()+
  theme_classic()
ggsave("common_lfc_barplot.pdf")
  
# what about genes unique to each treatment?

Acer_stresshardening_TvU_annotDGEs %>% 
  select(gene, annot) -> Acer_stresshardening_TvU_annotDGEs

CCC_vs_nursery_annotatedDGEs %>% 
  select(gene, annot) -> CCC_vs_nursery_annotatedDGEs

unique_CCCN <- anti_join(CCC_vs_nursery_annotatedDGEs, Acer_stresshardening_TvU_annotDGEs)
write_csv(unique_CCCN, "unique_CCCvsNursery.csv")
unique_SH <- anti_join(Acer_stresshardening_TvU_annotDGEs, CCC_vs_nursery_annotatedDGEs)
write_csv(unique_SH, "unique_stresshardening.csv")
