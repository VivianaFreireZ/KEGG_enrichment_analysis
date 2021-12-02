
## Author: Viviana Freire-Zapata
## Date: 12/12/21
## Topic: KEGG enrichment analysis 


## Required libraries
library(tidyverse)
library(clusterProfiler)
library(RColorBrewer)
library(here)

## This function will set up your working directory

here()

## Loading data

metaG <- read_csv( here('input', 'metaG-raw-copy-numbers.csv')) 

metaT_DEG <- read_csv(here('input', 'pre-drought-vs-drought.csv')) %>% 
  drop_na(padj)

metaB <- read_csv(here('input' , 'KEGG-table-corrected.csv'))

metadata <- read_csv(here('input', 'metadata-metaG-metaT-metaB.csv'))

metadata_metaB <- read_csv(here('input', 'metadata_metaB.csv'))

## Creating vector with Human pathways to be deleted

delete <- c("Pathways of neurodegeneration - multiple diseases","Amyotrophic lateral sclerosis", 
            "Alzheimer disease", "Huntington disease", "Parkinson disease", "Pathways of neurodegeneration - multiple diseases",
            "Amyotrophic lateral sclerosis", "Alzheimer disease", "Huntington disease", 
            "Pathways in cancer", "Prion disease", "Thermogenesis", "Central carbon metabolism in cancer",
            "Shigellosis", "Human papillomavirus infection", "Diabetic cardiomyopathy", "Coronavirus disease - COVID-19",
            "Chemical carcinogenesis - reactive oxygen species", "Taurine and hypotaurine metabolism")

## KEGG ENRICHMENT ANALYSIS 

###### METAGENOME ######

## Organizing data for analysis

## Transforming data to presence/absence form

for(i in 3:36){
  metaG[metaG[,i] > 0, i] <- 1 
}

## Treatment: Drought and time
## Calculating number of samples per treatment in which KO's were present 

metadata_G_drought <- metadata %>% 
  filter(Dataset == "metaG") %>% 
  filter(Condition == "Drought")

metaG_drought <- metaG %>% 
  select(Feature, all_of(metadata_G_drought$SampleID)) %>% 
  mutate(Time0 = select(., contains("time0")) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(Time6h = select(., contains("time6")) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(Time48h = select(., contains("drought_48")) %>% rowSums(na.rm = TRUE))
  
## Creating tables per time
## Filtering KO's that were present in all samples per time

metaG_drought_0 <- metaG_drought %>% 
  select(c(Feature, Time0)) %>% 
  filter(Time0 == 6) %>% 
  mutate(Time = "t0", Condition = "drought") %>% 
  select(-Time0)

metaG_drought_6 <- metaG_drought %>% 
  select(c(Feature, Time6h)) %>% 
  filter(Time6h == 7) %>% 
  mutate(Time = "t6", Condition = "drought") %>% 
  select(-Time6h)

metaG_drought_48 <- metaG_drought %>% 
  select(c(Feature, Time48h)) %>% 
  filter(Time48h == 5) %>% 
  mutate(Time = "t48", Condition = "drought") %>% 
  select(-Time48h)

##Joining tables 

drought_ko <- rbind(metaG_drought_0, metaG_drought_6, metaG_drought_48)

## Treatment: Pre_Drought and time
## Calculating number of samples per treatment in which KO's were present 

metadata_G_pre <- metadata %>% 
  filter(Dataset == "metaG") %>% 
  filter(Condition == "PreDrought")

metaG_pre <- metaG %>% 
  select(Feature, all_of(metadata_G_pre$SampleID)) %>% 
  mutate(Time0 = select(., contains("time0")) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(Time6h = select(., contains("time6")) %>% rowSums(na.rm = TRUE)) %>% 
  mutate(Time48h = select(., contains("drought_48")) %>% rowSums(na.rm = TRUE))

## Creating tables per time
## Filtering KO's that were present in all samples per time

metaG_pre_0 <- metaG_pre %>% 
  select(c(Feature, Time0)) %>% 
  filter(Time0 == 5) %>% 
  mutate(Time = "t0", Condition = "predrought") %>% 
  select(-Time0)

metaG_pre_6 <- metaG_pre %>% 
  select(c(Feature, Time6h)) %>% 
  filter(Time6h == 5) %>% 
  mutate(Time = "t6", Condition = "predrought") %>% 
  select(-Time6h)

metaG_pre_48 <- metaG_pre %>% 
  select(c(Feature, Time48h)) %>% 
  filter(Time48h == 6) %>% 
  mutate(Time = "t48", Condition = "predrought") %>% 
  select(-Time48h)

##Joining tables 

pre_drought_ko <- rbind(metaG_pre_0, metaG_pre_6, metaG_pre_48)

## Joining drought and predrought ko tables

metaG_ko <- rbind(drought_ko, pre_drought_ko)

metaG_ko$Feature <- str_remove(metaG_ko$Feature, "KO:")


## ENRICHMENT ANALYSIS

metagenome_enrich <- compareCluster(Feature~Condition+Time, data = metaG_ko, fun = enrichKEGG,
                                    organism = "ko", minGSSize = 10)

metagenome_enrich@compareClusterResult$Cluster <- factor(metagenome_enrich@compareClusterResult$Cluster, 
                                                           levels = c("drought.t0", "drought.t6",
                                                                      "drought.t48", "predrought.t0", 
                                                                      "predrought.t6","predrought.t48"))
## Plotting results

metaG_plot <- dotplot(metagenome_enrich, showCategory = 20, 
                   title = "Metagenome",
                   font.size = 12) +
  facet_wrap(~Condition, scales = 'free_x')

metaG_plot

ggsave(filename = 'output/KEGG_enrichment_metagenome_per_time.png', metaG_plot, height = 10, width = 10)

## Saving enrichment table results

genome_table_results <- metagenome_enrich@compareClusterResult

genome_table_results$GeneRatio <- paste0(' ', genome_table_results$GeneRatio)

genome_table_results$BgRatio <- paste0(' ', genome_table_results$GeneRatio)

write_csv(genome_table_results, 'output/metagenome_enrichment_result_table.csv')


## Treatment: Drought
## Calculating number of samples per treatment in which KO's were present 

metaG_drought_only <- metaG %>% 
  select(Feature, all_of(metadata_G_drought$SampleID)) %>% 
  mutate(Drought = select(., contains("_drought_")) %>% rowSums(na.rm = TRUE))

## Creating table
## Filtering KO's that were present in all samples per condition

metaG_drought_ko <- metaG_drought_only %>% 
  select(c(Feature, Drought)) %>% 
  filter(Drought == 18) %>% # you can change this parameter to be less restrictive
  mutate(Condition = "drought") %>% 
  select(-Drought)

## Treatment: Pre_Drought 
## Calculating number of samples per treatment in which KO's were present 

metaG_pre_only <- metaG %>% 
  select(Feature, all_of(metadata_G_pre$SampleID)) %>% 
  mutate(Predrought = select(., contains("_pre_")) %>% rowSums(na.rm = TRUE))

## Creating table
## Filtering KO's that were present in all samples per condition

metaG_pre_ko <- metaG_pre_only %>% 
  select(c(Feature, Predrought)) %>% 
  filter(Predrought == 16) %>% 
  mutate(Condition = "predrought") %>% 
  select(-Predrought)

## Joining drought and predrought ko tables

metaG_ko_condition <- rbind(metaG_drought_ko, metaG_pre_ko)

metaG_ko_condition$Feature <- str_remove(metaG_ko_condition$Feature, "KO:")


## ENRICHMENT ANALYSIS

metagenome_enrich_condition <- compareCluster(Feature~Condition, data = metaG_ko_condition, fun = enrichKEGG,
                                    organism = "ko", minGSSize = 5)

## Plotting results

metaG_plot_condition <- dotplot(metagenome_enrich_condition, showCategory = 20, 
                      title = "Metagenome_condition",
                      font.size = 12) +
  facet_wrap(~Condition, scales = 'free_x')

metaG_plot_condition

ggsave(filename = 'output/KEGG_enrichment_metagenome_per_condition.png', metaG_plot, height = 10, width = 10)

## Saving enrichment table results

condition_results <- metagenome_enrich_condition@compareClusterResult

condition_results$GeneRatio <- paste0(' ', condition_results$GeneRatio)

condition_results$BgRatio <- paste0(' ', condition_results$GeneRatio)

write_csv(genome_table_results, 'output/metagenome_enrichment_result_table_condition.csv')


###### METATRASNCRIPTOME ######

## Drought vs Pre_drought

DEG_significant <- metaT_DEG %>% 
  filter(padj < 0.05) %>% ## setting a padj < 0.05 for significant DEG, you can change this if you want 
  mutate(DE = ifelse(log2FoldChange > 0, "upregulated", "downregulated")) %>% 
  rename(Feature = ...1) %>% 
  select(Feature, DE)

DEG_significant$Feature <- str_remove(DEG_significant$Feature, "KO:")

## ENRICHMENT ANALYSIS 

metaT_enrich <- compareCluster(Feature~DE, data = DEG_significant, fun = enrichKEGG,
                                    organism = "ko", minGSSize = 10)

## Filtering out Human KEGG pathways

metaT_enrich@compareClusterResult <- metaT_enrich@compareClusterResult %>% 
  filter(!(Description %in% delete))


metaT_plot <- dotplot(metaT_enrich, showCategory = 20, 
                         title = "Metatranscriptome",
                         font.size = 12)+
  facet_wrap(~DE, scales = 'free_x')

metaT_plot

ggsave(filename = 'output/KEGG_enrichment_metatranscriptome.png', metaT_plot, height = 12, width = 10)

trans_table_results <- metaT_enrich@compareClusterResult

trans_table_results$GeneRatio <- paste0(' ', trans_table_results$GeneRatio)

trans_table_results$BgRatio <- paste0(' ', trans_table_results$GeneRatio)

write_csv(trans_table_results, 'output/metatranscriptome_enrichment_result_table.csv')


###### METABOLOME ######

## Organizing data 

df <- read_csv(here('input', 'Report_processed_MolecFormulas.csv')) %>% 
  pivot_longer(all_of(metadata_metaB$SampleID), names_to = 'SampleID', values_to = 'NormIntensity') %>% 
  filter(NormIntensity != 0)

df_kegg <- left_join(df, metaB, by = c('Mass', 'MolecularFormula'))

df_kegg_meta <- left_join(df_kegg, metadata_metaB, by = 'SampleID')

df_kegg_meta_l <- df_kegg_meta %>% 
  pivot_wider(names_from = Condition, values_from = Condition, names_prefix = 'GRP_') %>% 
  select(Mass, contains('KEGG_'), contains('GRP_')) %>% 
  group_by(Mass) %>% 
  fill(contains('GRP_'), .direction = 'downup')  %>% 
  distinct() %>% 
  unite(Presence, contains('GRP_'), sep = ', ', na.rm = TRUE)


## Creating KEGG universe for metabolite analysis 
## Loading .csv file 

compound_df <- read.csv("KEGG_compound_db.csv")

path2id <- compound_df %>% 
  select(KEGG_pathway, KEGG_id) %>% 
  separate_rows(KEGG_pathway, sep = ';') %>% 
  filter(!is.na(KEGG_pathway))

## Metabolite vector of interest

## Drought

metabolite_drought <- df_kegg_meta_l %>%
  ungroup() %>% 
  filter(Presence == "Drought") %>% #selecting KEGG compounds unique in drought
  select(KEGG_id) %>%
  filter(!is.na(KEGG_id)) %>% 
  mutate(Condition = "drought")

## Pre_Drought

metabolite_pre <- df_kegg_meta_l %>%
  ungroup() %>% 
  filter(Presence == "PreDrought") %>%#selecting KEGG compounds unique in predrought
  select(KEGG_id) %>%
  filter(!is.na(KEGG_id)) %>% 
  mutate(Condition = "predrought")

# Joining tables

metabolite_kegg <- rbind(metabolite_drought, metabolite_pre)

## Enrichment analysis

metabolite_enrich <- compareCluster(KEGG_id~Condition, data = metabolite_kegg, fun = enricher,
                                    TERM2GENE = path2id, minGSSize = 5)


metaB_plot <- dotplot(metabolite_enrich, showCategory = 20, 
                   title = "Metabolite",
                   font.size = 12)
metaB_plot

ggsave(filename = 'output/KEGG_metabolite_erichment.png', metaB_plot)

metabolite_table_result <- metabolite_enrich@compareClusterResult

metabolite_table_result$GeneRatio <- paste0(' ', metabolite_table_result$GeneRatio)
metabolite_table_result$BgRatio <- paste0(' ', metabolite_table_result$BgRatio)

write_csv(metabolite_table_result, 'output/metabolite_enrichment_result.csv')



## JOINING metaG and metaB results in one graph

## Analysis by condition

gene <- metagenome_enrich_condition@compareClusterResult %>% 
  mutate(dataset = "metaG")

metabo <- metabolite_enrich@compareClusterResult %>% 
  mutate(dataset = "metaB")

join <- rbind(gene, metabo) %>% 
  separate(GeneRatio, c("numerator", "denominator"), sep = "/") %>% 
  filter(p.adjust < 0.05) %>% 
  mutate(numerator = as.numeric(numerator), 
         denominator = as.numeric(denominator),
         GeneRatio = numerator/ denominator) %>% 
  group_by(dataset) %>%
  slice_max(order_by = c(Count, GeneRatio), n = 30) 


join_plot <- ggplot(join, aes( x= GeneRatio, y = Description, color = dataset, size = Count))+
  geom_point()+
  theme_bw()+
  facet_wrap(~Condition)
join_plot

ggsave(filename = 'output/KEGG_join_plot.png', join_plot, width = 10, height = 6)


