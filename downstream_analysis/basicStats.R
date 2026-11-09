#basis stats

library(remotes)
library(readr)
library(Biostrings)
library(treeio)
library(tidyverse)
library(UpSetR)
library(stringr)
library(patchwork)
load("data/A2TEA_finished_potatofix2.RData")


##########
conserved_HOGS <- HOG_level_list$all_species_overview %>% filter(if_all(ends_with("_total"), ~ . > 0)) %>% pull(HOG)

#test <-HOG_DE.a2tea %>% filter(species=="Solanum_pennellii") %>% filter(HOG!="singleton")
genes_per_species <- HOG_DE.a2tea %>% #replace_na(list(species = "Solanum_tuberosum")) %>% 
  group_by(species) %>% 
  summarise(genes_total = n()) 

genes_per_species_in_hogs <- HOG_DE.a2tea %>% replace_na(list(species = "Solanum_tuberosum")) %>% filter(HOG!="singleton") %>% 
  group_by(species) %>% 
  summarise(genes_in_HOGs = n()) 

genes_per_species_in_conserved_hogs <- HOG_DE.a2tea %>% replace_na(list(species = "Solanum_tuberosum")) %>% filter(HOG!="singleton") %>% filter(HOG %in% conserved_HOGS) %>% 
  group_by(species) %>% 
  summarise(genes_in_conserved_HOGs = n()) 

HOG_DE.a2tea %>%
  group_by(species) %>%
  summarise(n = n(),
            in_HOG    = sum(HOG != "singleton"),
            singleton = sum(HOG == "singleton"),
            other     = n - in_HOG - singleton)

#hier noch interessant: gesamtzahl (nicht nur signifikant -> wie hoch? muss ich noch gene ausschließen? Das gleiche gilt auch für die konservierten Genefamilien)
DE_per_species <- HOG_DE.a2tea %>% 
  filter(significant=="yes") %>% 
  #replace_na(list(species = "Solanum_tuberosum")) %>% #does not matter, because none of the NAs are DEGs 
  mutate(
    regulation = case_when(
      log2FoldChange > 0 ~ "up",
      log2FoldChange < 0 ~ "down",
    ),
    .groups = "drop"
  ) %>% 
  group_by(species) %>%
  summarise(
    total_significant = n(),
    upregulated = sum(regulation == "up"),
    downregulated = sum(regulation == "down"),
  )

DE_per_species_inHogs <- HOG_DE.a2tea %>% 
  filter(significant=="yes") %>% filter(HOG!="singleton") %>% 
  #replace_na(list(species = "Solanum_tuberosum")) %>% #does not matter, because none of the NAs are DEGs 
  mutate(
    regulation = case_when(
      log2FoldChange > 0 ~ "up",
      log2FoldChange < 0 ~ "down",
    ),
    .groups = "drop"
  ) %>% 
  group_by(species) %>%
  summarise(
    significant_in_HOGs = n(),
    upregulated = sum(regulation == "up"),
    downregulated = sum(regulation == "down"),
  )

DE_per_species_conserved <- HOG_DE.a2tea %>% 
  filter(significant=="yes") %>% 
  filter(HOG %in% conserved_HOGS) %>% 
  #replace_na(list(species = "Solanum_tuberosum")) %>% #does not matter, because none of the NAs are DEGs 
  mutate(
    regulation = case_when(
      log2FoldChange > 0 ~ "up",
      log2FoldChange < 0 ~ "down",
    ),
    .groups = "drop"
  ) %>% 
  group_by(species) %>%
  summarise(
    significant_in_conserved_HOGs = n(),
    upregulated = sum(regulation == "up"),
    downregulated = sum(regulation == "down"),
  )


######Visualisierung:








genes_per_species_per_hog <- HOG_DE.a2tea %>%
  filter(HOG %in% conserved_HOGS) %>% 
  replace_na(list(species = "Solanum_tuberosum")) %>% 
  group_by(species, HOG) %>%
  summarise(gene_count = n(), .groups = 'drop') #%>% filter(gene_count>1)



mean_genes_per_species <- genes_per_species_per_hog %>%
  group_by(species) %>%
  summarise(mean_genes = mean(gene_count),
            median_genes = median(gene_count))

# for conserved set
summary_stats <- genes_per_species_per_hog %>%
  group_by(species) %>%
  summarise(
    total_hogs = n(),
    total_genes = sum(gene_count),
    single_copy = sum(gene_count == 1),
    two_genes = sum(gene_count == 2),
    three_genes = sum(gene_count == 3),
    multi_copy = sum(gene_count > 3),
    mean_genes = mean(gene_count),
    median_genes = median(gene_count),
    max_genes = max(gene_count)
  )


DE_per_species_per_hog <- HOG_DE.a2tea %>%
  filter(HOG %in% conserved_HOGS) %>% 
  replace_na(list(species = "Solanum_tuberosum")) %>% 
  filter(significant=="yes") %>% 
  group_by(species, HOG) %>%
  summarise(gene_count = n(), .groups = 'drop') #%>% filter(gene_count>1)


# for conserved set
summary_stats_DE <- DE_per_species_per_hog %>%
  group_by(species) %>%
  summarise(
    total_hogs = n(),
    total_genes = sum(gene_count),
    single_copy = sum(gene_count == 1),
    two_genes = sum(gene_count == 2),
    three_genes = sum(gene_count == 3),
    multi_copy = sum(gene_count > 3),
    mean_genes = mean(gene_count),
    median_genes = median(gene_count),
    max_genes = max(gene_count)
  )


#########
#plotting
#not really useful yet (log transform? use only counts >2)
genes_per_species_per_hog_over1 <- HOG_DE.a2tea %>%
  filter(HOG %in% conserved_HOGS) %>% 
  replace_na(list(species = "Solanum_tuberosum")) %>% 
  group_by(species, HOG) %>%
  summarise(gene_count = n(), .groups = 'drop') %>% filter(gene_count>1)

ggplot(genes_per_species_per_hog_over1, aes(x = gene_count)) +
  geom_histogram( fill = "steelblue", color = "black", alpha = 0.7) +
  #geom_vline(data = mean_genes_per_species, 
  #          aes(xintercept = mean_genes), 
  #         color = "red", linetype = "dashed", linewidth = 1) +
  facet_wrap(~ species, ncol = 5) +
  labs(title = "Distribution of Genes per HOG by Species",
       x = "Number of Genes per HOG",
       y = "Frequency") +
  theme(strip.text = element_text(face = "bold"))



