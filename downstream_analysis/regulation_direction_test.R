
# Direction of regulation in expanded HOGs — per-species enrichment test

suppressPackageStartupMessages(library(tidyverse))


## 0. LOAD DATA  &  SETTINGS

load("data/A2TEA_finished_potatofix2.RData")

outdir <- "data/GOEnrichment_potatofix"


## 1. BACKGROUND  &  EXPANDED SET

## Background = Conserved Set (HOGs present in every species).
background <- HOG_level_list$all_species_overview %>%
  filter(if_all(ends_with("_total"), ~ . > 0)) %>%
  pull(HOG)

## Expanded in both tolerant Solanum species (hypothesis_1 = S. pennellii,
## hypothesis_2 = S. tuberosum), restricted to the Conserved Set. Same
## definition as in GO_enrichment_topGO_withGoji.R (164 HOGs).
expanded_penn <- HOG_level_list$hypothesis_1 %>%
  filter(HOG %in% background, expansion == "yes") %>% pull(HOG)
expanded_tub  <- HOG_level_list$hypothesis_2 %>%
  filter(HOG %in% background, expansion == "yes") %>% pull(HOG)
expanded_both <- intersect(expanded_penn, expanded_tub)


## 2. REGULATION CATEGORY PER HOG x SPECIES

## hypothesis_7 is used only as the HOG table here (it carries DE data from all
## four species with expression data); it is NOT filtered by expansion, so this
## is every conserved HOG. A former filter total_sigDE > 1 (>= 2 DEGs pooled
## across species, as in the reference paper) was removed: the category is
## assigned per species, where that filter does not prevent single-DEG
## categories (see section 3).
HOGS2filter4 <- HOG_level_list$hypothesis_7 %>% filter(HOG %in% background)

## significant DEGs of these HOGs
filtered_df <- HOG_DE.a2tea %>%
  filter(HOG %in% HOGS2filter4$HOG) %>% filter(significant == "yes")

## per species: all DEGs up / all down / both
regulation_long <- filtered_df %>%
  group_by(HOG, species) %>%
  summarise(
    regulation = case_when(
      all(log2FoldChange > 0) ~ "up",
      all(log2FoldChange < 0) ~ "down",
      any(log2FoldChange > 0) & any(log2FoldChange < 0) ~ "up_and_down"
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(regulation))

print(count(regulation_long, species, regulation))


## 3. NUMBER OF DEGs PER HOG x SPECIES

## Most HOG x species combinations contain a single DEG, so the category mostly
## reflects the direction of one gene and "up_and_down" is rarely possible.
filtered_df %>%
  count(HOG, species) %>%
  count(n, name = "pairs") %>%
  mutate(pct = round(100 * pairs / sum(pairs), 1))


## 4. Test function
## For one species and one regulation category:
##   foreground = HOGs categorised in that species AND expanded
##   background = all HOGs categorised in that species
## one-sided Fisher's exact test for over-representation of the category in
## the foreground. The comparison stays within one species and therefore within
## one cold experiment.

test_regulation <- function(reg_tbl, expanded_hogs, sp, category) {

  bg <- reg_tbl %>% dplyr::filter(species == sp)
  fg <- bg %>% dplyr::filter(HOG %in% expanded_hogs)

  # 2x2: rows = in foreground / only in background, cols = category / other
  a <- sum(fg$regulation == category)                 # expanded, category
  b <- sum(fg$regulation != category)                 # expanded, other
  c <- sum(bg$regulation == category) - a             # not expanded, category
  d <- sum(bg$regulation != category) - b             # not expanded, other

  if (a + b == 0) {
    return(tibble(species = sp, category = category, n_expanded = 0,
                  pct_expanded = NA_real_, n_background = a + c,
                  pct_background = NA_real_, p = NA_real_))
  }

  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
                    alternative = "greater")

  tibble(
    species        = sp,
    category       = category,
    n_expanded     = a,
    n_expanded_tot = a + b,
    pct_expanded   = round(100 * a / (a + b), 1),
    n_background   = a + c,
    n_bg_tot       = a + b + c + d,
    pct_background = round(100 * (a + c) / (a + b + c + d), 1),
    p              = signif(ft$p.value, 3)
  )
}


## 5. Run for every species x category

species_list  <- sort(unique(regulation_long$species))
category_list <- c("up", "down", "up_and_down")

reg_results <- tidyr::crossing(sp = species_list, category = category_list) %>%
  purrr::pmap_dfr(function(sp, category)
    test_regulation(regulation_long, expanded_both, sp, category))

print(reg_results, n = Inf)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(reg_results, row.names = FALSE,
          file = file.path(outdir, "Table4_regulation_direction.csv"))


## 6. Power check
## Smallest proportion of "up" in the foreground that would reach p <= 0.05,
## given the foreground size and the background rate of that species.

power_check <- reg_results %>%
  dplyr::filter(category == "up") %>%
  dplyr::select(species, n_expanded_tot, pct_background) %>%
  mutate(
    min_detectable_pct = map2_dbl(n_expanded_tot, pct_background, function(n, bgpct) {
      if (is.na(n) || n == 0) return(NA_real_)
      p0 <- bgpct / 100
      k <- 0:n
      ok <- which(stats::pbinom(k - 1, n, p0, lower.tail = FALSE) <= 0.05)
      if (!length(ok)) return(NA_real_)
      round(100 * (k[min(ok)]) / n, 1)
    })
  )

print(power_check)
