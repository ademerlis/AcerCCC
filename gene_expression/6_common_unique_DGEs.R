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

str(Acer_stresshardening_TvU_annotDGEs)

Acer_stresshardening_TvU_annotDGEs %>% 
  rename(baseMean_SH = baseMean)

CCC_vs_nursery_annotatedDGEs %>% 
  rownames_to_column(var = "gene") -> CCC_vs_nursery_annotatedDGEs
  
# what about genes unique to each treatment?

unique_Untreated_vs_Initial <- anti_join(Untreated_vs_Initial_sig, Treated_vs_Initial_sig)
unique_Untreated_vs_Initial <- anti_join(unique_Untreated_vs_Initial, Treated_vs_Untreated_sig)
unique_Untreated_vs_Initial %>% 
  as.data.frame() %>% 
 left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                     sep = "\t",
                     quote="", fill=FALSE) %>%
            mutate(gene = V1,
                   annot = V2) %>%
            dplyr::select(-V1, -V2), by = "gene") -> unique_Untreated_vs_Initial_annotated

str(unique_Untreated_vs_Initial_annotated) #2435

load("RData_files/realModels_Acer.RData")

Treatment_Untreated_vs_Initial=results(dds,contrast=c("Treatment","Untreated","Initial"))

Treatment_Untreated_vs_Initial %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  dplyr::filter(padj<0.05) %>% 
  right_join(., unique_Untreated_vs_Initial_annotated, by = "gene") %>%
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") %>% view()
  write_csv("unique_Untreated_vs_Initial_annotated_KOG.csv")

write_csv(unique_Untreated_vs_Initial_annotated, "unique_Untreated_vs_Initial_annotated.csv")

unique_Treated_vs_Initial <- anti_join(Treated_vs_Initial_sig, Untreated_vs_Initial_sig)
unique_Treated_vs_Initial <- anti_join(unique_Treated_vs_Initial, Treated_vs_Untreated_sig)
unique_Treated_vs_Initial %>% 
  as.data.frame() %>% 
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") -> unique_Treated_vs_Initial_annotated

str(unique_Treated_vs_Initial_annotated) #1429

Treatment_Treated_vs_Initial=results(dds,contrast=c("Treatment","Treated","Initial"))

Treatment_Treated_vs_Initial %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  dplyr::filter(padj<0.05) %>% 
  right_join(., unique_Treated_vs_Initial_annotated, by = "gene") %>%
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") %>% view()
  write_csv("unique_Treated_vs_Initial_annotated_KOG.csv")

write_csv(unique_Treated_vs_Initial_annotated, "unique_Treated_vs_Initial_annotated.csv")

unique_Treated_vs_Untreated <- anti_join(Treated_vs_Untreated_sig, Untreated_vs_Initial_sig)
unique_Treated_vs_Untreated <- anti_join(unique_Treated_vs_Untreated, Treated_vs_Initial_sig)
unique_Treated_vs_Untreated %>% 
  as.data.frame() %>% 
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2geneName.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     annot = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") -> unique_Treated_vs_Untreated_annotated

str(unique_Treated_vs_Untreated_annotated) #116 

write_csv(unique_Treated_vs_Untreated_annotated, "unique_Treated_vs_Untreated_annotated.csv")

Treatment_Treated_vs_Untreated=results(dds,contrast=c("Treatment","Treated","Untreated"))

Treatment_Treated_vs_Untreated %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  drop_na(padj) %>% 
  dplyr::filter(padj<0.05) %>% 
  right_join(., unique_Treated_vs_Untreated_annotated, by = "gene") %>%
  left_join(read.table(file = "bioinformatics/Acervicornis_iso2kogClass.tab",
                       sep = "\t",
                       quote="", fill=FALSE) %>%
              mutate(gene = V1,
                     KOG = V2) %>%
              dplyr::select(-V1, -V2), by = "gene") %>% 
  write_csv("unique_Treated_vs_Untreated_annotated_KOG.csv")
