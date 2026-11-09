library(remotes)
library(readr)
library(Biostrings)
library(treeio)
library(tidyverse)
library(UpSetR)
library(stringr)
#load("data/A2TEA_finished_8_29.RData") #old
load("data/A2TEA_finished_potatofix2.RData")

workset <- HOG_level_list$all_species_overview %>%
  dplyr::select(HOG, ends_with("_total"))

# remove '_total' from column names
#str_remove(colnames(workset), "_total")
workset <- workset %>% setNames(names(workset) %>% stringr::str_replace("_total",""))
str(workset)
# a vector with all species names for currently selected hypothesis
species_vec <- names(workset %>% dplyr::select(-HOG))
head(workset)

workset_1 <- workset %>% mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0))) %>% data.frame()
# using map_at from purrr
# we perform the operation for all columns in species_vec
# filter/only keep HOG "rows" per species if species possesses at least 1 gene inside the HOG
# with purrr::list_modify("HOG" = NULL) the untouched pure HOG column (now a list) is removed
# genotypeList <- workset %>% purrr::map_at(species_vec, 
#                                           ~filter(workset, .!= 0) %>% 
#                                             dplyr::select(HOG) %>% 
#                                             as_vector() %>%
#                                            unname() )%>% purrr::list_modify("HOG" = NULL)


species_list <- list("Solanum_pennellii" , "Solanum_tuberosum","Lycium_barbarum", "Solanum_lycopersicum","Capsicum_annuum")
#upset plot with intersection of genes in 
pdf("plots/8_29_plots/UpsetRHOGs.pdf", 15, 8, onefile = FALSE)
a <- upset(workset_1, queries = list(list(query = intersects, params = list("Solanum_pennellii" , "Solanum_tuberosum","Lycium_barbarum", "Solanum_lycopersicum","Capsicum_annuum"), color = "orange", active = T)),
           order.by = "freq",
           empty.intersections = "on",
           sets.bar.color ="cornflowerblue",
           main.bar.color = "cornflowerblue",
           matrix.color= "cornflowerblue",
           
           text.scale = 1.7,
           #number.angles = 45
)

a

dev.off()


                                                                                                               
#################
##upsetr for DE sets

#conserved set

conserved_OGs <- HOG_level_list$all_species_overview %>%
  #at least 1 gene from every species thus "conserved"
  filter_at(vars(c(ends_with("_total"))), all_vars(. > 0)) %>%
  pull(HOG)

workset <- HOG_level_list$all_species_overview %>% 
  filter(HOG %in% conserved_OGs) %>% 
  dplyr::select(HOG, ends_with("_sigDE"))



# remove '_sigDE' from column names
workset <- workset %>% setNames(names(workset) %>% stringr::str_replace("_sigDE",""))
workset_1 <- workset %>%
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0))) %>%
  data.frame()
species_vec <- names(workset %>% dplyr::select(-HOG))




a <- upset(workset_1,
      sets = species_vec,
      order.by = "freq",
      empty.intersections = "on",
      set_size.show = TRUE,
      set_size.numbers_size = 5,
      set_size.scale_max = 7000,      # headroom past the ~6000 Capsicum bar
      sets.bar.color = "cornflowerblue",
      main.bar.color = "cornflowerblue",
      matrix.color = "cornflowerblue",
      text.scale = 1.7,
      )
a



#####query-list for sets that include both tolerant (for coloring)
tolerant <- c("Solanum_tuberosum", "Solanum_pennellii")
others <- setdiff(species_vec, tolerant)

# every subset of the non-tolerant species, each paired with BOTH tolerant species
other_combos <- unlist(
  lapply(0:length(others), function(k) combn(others, k, simplify = FALSE)),
  recursive = FALSE
)

de_queries <- lapply(other_combos, function(extra) {
  list(query  = intersects,
       params = as.list(c(tolerant, extra)),
       color  = "#3B7DB4",   # the figure's blue; swap for cornflowerblue if you like
       active = TRUE)
})
###############

pdf("figures/UpsetRconservedHOGs.pdf", 15, 8, onefile = FALSE)
a <-upset(workset_1,
      sets = species_vec,
      order.by = "freq",
      empty.intersections = "on",
      set_size.show = TRUE,
      set_size.numbers_size = 5,
      set_size.scale_max = 7000,
      queries = de_queries,
      main.bar.color = "gray23",      # non-highlighted bars dark
      matrix.color   = "gray23",
      sets.bar.color = "gray23",
      text.scale = 1.7)
a
dev.off()

########################################
#conserved, expanded in 2tol

expanded_HOG_list <- HOG_level_list$hypothesis_4 %>% 
  filter(expansion=="yes") %>% 
  pull(HOG)

workset <- HOG_level_list$all_species_overview %>% 
  filter(HOG %in% conserved_OGs) %>% 
  filter(HOG %in% expanded_HOG_list) %>% 
  dplyr::select(HOG, ends_with("_sigDE"))



# remove '_sigDE' from column names
workset <- workset %>% setNames(names(workset) %>% stringr::str_replace("_sigDE",""))
workset_1 <- workset %>%
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0))) %>%
  data.frame()

pdf("figures/UpsetR2tolexpandedHOGs.pdf", 15, 8, onefile = FALSE)
a <- upset(workset_1,
           sets = species_vec,
           order.by = "freq",
           empty.intersections = "on",
           set_size.show = TRUE,
           set_size.numbers_size = 5,
           set_size.scale_max = 70,
           queries = de_queries,
           main.bar.color = "gray23",      # non-highlighted bars dark
           matrix.color   = "gray23",
           sets.bar.color = "gray23",
           text.scale = 1.7)
a

dev.off()





#conserved, expanded in 3tol

expanded_HOG_list <- HOG_level_list$hypothesis_7 %>% 
  filter(expansion=="yes") %>% 
  pull(HOG)

workset <- HOG_level_list$all_species_overview %>% 
  filter(HOG %in% conserved_OGs) %>% 
  filter(HOG %in% expanded_HOG_list) %>% 
  dplyr::select(HOG, ends_with("_sigDE"))



# remove '_total' from column names

workset <- workset %>% setNames(names(workset) %>% stringr::str_replace("_sigDE",""))
workset_1 <- workset %>%
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0))) %>%
  data.frame()

a <- upset(workset_1,
      sets = species_vec,
      order.by = "freq",
      empty.intersections = "on",
      set_size.show = TRUE,
      set_size.numbers_size = 5,
      set_size.scale_max = 70,
      queries = de_queries,
      main.bar.color = "gray23",      # non-highlighted bars dark
      matrix.color   = "gray23",
      sets.bar.color = "gray23",
      text.scale = 1.7)
a

#same set, but in color :D
tolerant <- c("Solanum_tuberosum", "Solanum_pennellii")
others <- setdiff(species_vec, tolerant)

# every subset of the non-tolerant species, each paired with BOTH tolerant species
other_combos <- unlist(
  lapply(0:length(others), function(k) combn(others, k, simplify = FALSE)),
  recursive = FALSE
)

de_queries <- lapply(other_combos, function(extra) {
  list(query  = intersects,
       params = as.list(c(tolerant, extra)),
       color  = "#3B7DB4",   # the figure's blue; swap for cornflowerblue if you like
       active = TRUE)
})

pdf("figures/UpsetR3tolexpandedHOGs.pdf", 15, 8, onefile = FALSE)
a <- upset(workset_1,
      sets = species_vec,
      order.by = "freq",
      empty.intersections = "on",
      set_size.show = TRUE,
      set_size.numbers_size = 5,
      set_size.scale_max = 70,
      queries = de_queries,
      main.bar.color = "gray23",      # non-highlighted bars dark
      matrix.color   = "gray23",
      sets.bar.color = "gray23",
      text.scale = 1.7)
a
dev.off()

####UpsetR Plot for Regulation categories per species


conserved_HOGS <- HOG_level_list$all_species_overview %>% filter(if_all(ends_with("_total"), ~ . > 0)) %>% pull(HOG)


#background - with or without at least 2 DEgenes in HOG
HOGS2filter4 <- HOG_level_list$hypothesis_4 %>% filter(HOG %in% conserved_HOGS) %>% filter(total_sigDE>1) #%>% filter(expansion=="yes")

#filter HOG_DE based on the HOGS we just identified, only use significant
filtered_df <- HOG_DE.a2tea %>%
  filter(HOG %in% HOGS2filter4$HOG) %>% filter(significant=="yes") #%>%  filter(species %in% c("Solanum_pennellii", "Solanum_tuberosum"))
hog_number_background <- length(unique(filtered_df$HOG))

#group hogs and summarise direction of regulation
hog_regulation_per_species <- filtered_df %>%
  #filter(!is.na(log2FoldChange)) %>%  # Exclude NAs in log2FoldChange
  group_by(HOG, species) %>%
  mutate( category_of_regulation_per_species= case_when(
    all(log2FoldChange > 0) ~ "up",
    all(log2FoldChange < 0) ~ "down",
    any(log2FoldChange > 0) & any(log2FoldChange < 0) ~ "up_and_down",
  )
  )


hog_regulation_per_species$category_of_regulation_per_species_with_name <-
  paste(hog_regulation_per_species$species, hog_regulation_per_species$category_of_regulation_per_species, sep=" ")


long_df <- hog_regulation_per_species %>%
  group_by(HOG, category_of_regulation_per_species_with_name) %>%
  summarise(gene_count = n(), .groups = 'drop')

wide_df <- long_df %>%
  pivot_wider(names_from = category_of_regulation_per_species_with_name, values_from = gene_count, values_fill = list(gene_count = 0))


workset <- wide_df

# a vector with all species names for currently selected hypothesis
species_vec <- names(workset %>% dplyr::select(-HOG))

# using map_at from purrr
# we perform the operation for all columns in species_vec
# filter/only keep HOG "rows" per species if species possesses at least 1 gene inside the HOG
# with purrr::list_modify("HOG" = NULL) the untouched pure HOG column (now a list) is removed
genotypeList <- workset %>% purrr::map_at(species_vec, 
                                          ~filter(workset, .!= 0) %>% 
                                            dplyr::select(HOG) %>% 
                                            as_vector() %>%
                                            unname() )%>% purrr::list_modify("HOG" = NULL)

# return genotypeList
#genotypeList

pdf("figures/Intersection_of_HOGs_with_regulation_categories_per_species_conservedHOGS_atleast2DEG.pdf", 16, 8, onefile = FALSE)

a <- upset(fromList(genotypeList), 
           order.by = "freq",
           nsets=12,
           #empty.intersections = "off",
           sets.bar.color ="cornflowerblue",
           main.bar.color = "cornflowerblue",
           matrix.color= "cornflowerblue",
           text.scale = 1.7,
           )

a
dev.off()

#complexupsetr might be easier for coloring (?)

library(ComplexUpset)
library(tidyverse)

sets <- setdiff(colnames(wide_df), "HOG")
memb <- wide_df %>% mutate(across(all_of(sets), ~ .x > 0))

cat_colors <- c(up = "#2C7FB8", down = "#D95F02", up_and_down = "#66A61E")

# category per set (up_and_down first, since it also ends in "down")
category <- case_when(
  grepl("up_and_down$", sets) ~ "up_and_down",
  grepl(" up$",   sets)       ~ "up",
  grepl(" down$",  sets)      ~ "down"
)
set_meta <- data.frame(set = sets, Regulation = category)

pdf("figures/Intersection_of_HOGs_with_regulation_categories_per_species_conservedHOGS_atleast2DEG_CompexUpsetR.pdf", 18, 8, onefile = FALSE)

a <- upset(
  memb,
  intersect         = sets,
  encode_sets       = FALSE,
  n_intersections   = 40,
  keep_empty_groups = TRUE,
  name              = "",
  sort_intersections_by = "cardinality",
  sort_sets         = "descending",
  set_sizes = (
    upset_set_size()
    + geom_text(aes(label = after_stat(count)), stat = "count", hjust = 1.1, size = 3)
    + expand_limits(y = max(sapply(memb[sets], sum)) * 1.1)
  ),
  stripes = upset_stripes(
    geom    = geom_segment(size = 7, alpha = 0.35),
    mapping = aes(color = Regulation),
    colors  = cat_colors,
    data    = set_meta
  ),
  base_annotations = list(
    "Intersection size" = (
      intersection_size(
        bar_number_threshold = 1,
        text = list(size = 2.5)          # smaller count labels (default is ~3.9)
      )
      + scale_y_continuous(expand = expansion(mult = c(0, 0.12)))
    )
  )
    )
a
dev.off()

###############################
#upsetR for 2tol-expanded set (at least 2DEG per HOG), intersection of regulation per species

#conserved_HOGS <- HOG_level_list$all_species_overview %>% filter(if_all(ends_with("_total"), ~ . > 0)) %>% pull(HOG)


#background - with or without at least 2 DEgenes in HOG
HOGS2filter4 <- HOG_level_list$hypothesis_4 %>% filter(HOG %in% conserved_HOGS) %>% filter(total_sigDE>1) %>% filter(expansion=="yes")

#filter HOG_DE based on the HOGS we just identified, only use significant
filtered_df <- HOG_DE.a2tea %>%
  filter(HOG %in% HOGS2filter4$HOG) %>% filter(significant=="yes") #%>%  filter(species %in% c("Solanum_pennellii", "Solanum_tuberosum"))
hog_number_background <- length(unique(filtered_df$HOG))

#group hogs and summarise direction of regulation
hog_regulation_per_species <- filtered_df %>%
  #filter(!is.na(log2FoldChange)) %>%  # Exclude NAs in log2FoldChange
  group_by(HOG, species) %>%
  mutate( category_of_regulation_per_species= case_when(
    all(log2FoldChange > 0) ~ "up",
    all(log2FoldChange < 0) ~ "down",
    any(log2FoldChange > 0) & any(log2FoldChange < 0) ~ "up_and_down",
  )
  )


hog_regulation_per_species$category_of_regulation_per_species_with_name <-
  paste(hog_regulation_per_species$species, hog_regulation_per_species$category_of_regulation_per_species, sep=" ")


long_df <- hog_regulation_per_species %>%
  group_by(HOG, category_of_regulation_per_species_with_name) %>%
  summarise(gene_count = n(), .groups = 'drop')

wide_df <- long_df %>%
  pivot_wider(names_from = category_of_regulation_per_species_with_name, values_from = gene_count, values_fill = list(gene_count = 0))


workset <- wide_df

# a vector with all species names for currently selected hypothesis
species_vec <- names(workset %>% dplyr::select(-HOG))

# using map_at from purrr
# we perform the operation for all columns in species_vec
# filter/only keep HOG "rows" per species if species possesses at least 1 gene inside the HOG
# with purrr::list_modify("HOG" = NULL) the untouched pure HOG column (now a list) is removed
genotypeList <- workset %>% purrr::map_at(species_vec, 
                                          ~filter(workset, .!= 0) %>% 
                                            dplyr::select(HOG) %>% 
                                            as_vector() %>%
                                            unname() )%>% purrr::list_modify("HOG" = NULL)

# return genotypeList
#genotypeList

pdf("figures/Intersection_of_HOGs_with_regulation_categories_per_species_2tolexpHOGS_atleast2DEG.pdf", 16, 8, onefile = FALSE)

a <- upset(fromList(genotypeList), 
           order.by = "freq",
           nsets=12,
           #empty.intersections = "off",
           sets.bar.color ="cornflowerblue",
           main.bar.color = "cornflowerblue",
           matrix.color= "cornflowerblue",
           text.scale = 1.7,
)

a
dev.off()

#complexupsetr might be easier for coloring (?)

#library(ComplexUpset)
#library(tidyverse)

sets <- setdiff(colnames(wide_df), "HOG")
memb <- wide_df %>% mutate(across(all_of(sets), ~ .x > 0))

cat_colors <- c(up = "#2C7FB8", down = "#D95F02", up_and_down = "#66A61E")

# category per set (up_and_down first, since it also ends in "down")
category <- case_when(
  grepl("up_and_down$", sets) ~ "up_and_down",
  grepl(" up$",   sets)       ~ "up",
  grepl(" down$",  sets)      ~ "down"
)
set_meta <- data.frame(set = sets, Regulation = category)

pdf("figures/Intersection_of_HOGs_with_regulation_categories_per_species_2tolexpHOGS_atleast2DEG_CompexUpsetR.pdf", 18, 8, onefile = FALSE)

a <- upset(
  memb,
  intersect         = sets,
  encode_sets       = FALSE,
  n_intersections   = 40,
  keep_empty_groups = TRUE,
  name              = "",
  sort_intersections_by = "cardinality",
  sort_sets         = "descending",
  set_sizes = (
    upset_set_size()
    + geom_text(aes(label = after_stat(count)), stat = "count", hjust = 1.1, size = 3)
    + expand_limits(y = max(sapply(memb[sets], sum)) * 1.1)
  ),
  stripes = upset_stripes(
    geom    = geom_segment(size = 7, alpha = 0.35),
    mapping = aes(color = Regulation),
    colors  = cat_colors,
    data    = set_meta
  ),
  base_annotations = list(
    "Intersection size" = (
      intersection_size(
        bar_number_threshold = 1,
        text = list(size = 2.5)          # smaller count labels (default is ~3.9)
      )
      + scale_y_continuous(expand = expansion(mult = c(0, 0.12)))
    )
  )
)
a
dev.off()


