# ============================================================
# Build the 76-HOG candidate table for the appendix
# (Supplementary Table tab:candidate-hogs)
#
# Input : data/A2TEA_finished_potatofix2.RData
# Output: tables/candidate_rows.tex  ->  \input{} into the longtable
# ============================================================

library(tidyverse)

load("data/A2TEA_finished_potatofix2.RData")

# ---- 1. Conserved Set (14,996 HOGs) ------------------------
conserved <- HOG_level_list$all_species_overview %>%
  filter(if_all(ends_with("_total"), ~ . > 0)) %>%
  pull(HOG)

# ---- 2. Primary candidate set: hypothesis 4, DE in >=1 tolerant
candidates <- HOG_level_list$hypothesis_4 %>%
  filter(HOG %in% conserved, expansion == "yes") %>%
  filter(if_any(c(Solanum_pennellii_sigDE, Solanum_tuberosum_sigDE), ~ . > 0)) %>%
  pull(HOG)

stopifnot(length(candidates) == 76)   # sanity check against the thesis

# ---- 3. Copy numbers per species ---------------------------
copies <- HOG_level_list$all_species_overview %>%
  filter(HOG %in% candidates) %>%
  select(HOG, ends_with("_total"), ends_with("_sigDE"))

# ---- 4. AHRD annotation: most frequent description per HOG --
# SFA is a named list of per-species tibbles; column names contain
# dashes, hence the backticks.
ann <- bind_rows(SFA) %>%
  filter(HOG %in% candidates) %>%
  count(HOG, `Human-Readable-Description`, sort = TRUE) %>%
  group_by(HOG) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(HOG, annotation = `Human-Readable-Description`)

# ---- 5. DE evidence with direction -------------------------
# padj < 0.1 = DEG_FDR in config.yaml
de_dir <- HOG_DE.a2tea %>%
  filter(HOG %in% candidates,
         species %in% c("Solanum_pennellii", "Solanum_tuberosum"),
         padj < 0.1) %>%
  group_by(HOG, species) %>%
  summarise(dir = case_when(all(log2FoldChange > 0) ~ "up",
                            all(log2FoldChange < 0) ~ "down",
                            TRUE                    ~ "both"),
            .groups = "drop") %>%
  mutate(sp = if_else(species == "Solanum_pennellii", "penn", "tub")) %>%
  group_by(HOG) %>%
  summarise(de = paste(paste(sp, dir), collapse = ", "), .groups = "drop")

# ---- 6. Also expanded in L. barbarum (hypothesis 7)? --------
h7 <- HOG_level_list$hypothesis_7 %>%
  filter(expansion == "yes") %>%
  pull(HOG)

# ---- 7. The six families examined in Section 3.7 ------------
examined <- c("N0.HOG0001280", "N0.HOG0006681", "N0.HOG0006736",
              "N0.HOG0005673", "N0.HOG0000291", "N0.HOG0003220")

# ---- 8. Assemble and write LaTeX rows -----------------------
esc <- function(x) str_replace_all(x, "([#$%&_])", "\\\\\\1")

tab <- copies %>%
  left_join(ann,    by = "HOG") %>%
  left_join(de_dir, by = "HOG") %>%
  arrange(HOG) %>%
  mutate(
    mark = paste0(if_else(HOG %in% examined, "$^{*}$",       ""),
                  if_else(HOG %in% h7,       "$^{\\dagger}$", "")),
    row  = sprintf("%s%s & %s & %d & %d & %d & %d & %d & %s \\\\",
                   esc(HOG), mark,
                   esc(str_trunc(replace_na(annotation, "--"), 55)),
                   Solanum_pennellii_total,  Solanum_tuberosum_total,
                   Lycium_barbarum_total,    Solanum_lycopersicum_total,
                   Capsicum_annuum_total,
                   replace_na(de, "--"))
  )

writeLines(tab$row, "tables/candidate_rows.tex")

# ---- 9. Bonus: the red [X] in Discussion 4.6 ----------------
# candidate HOGs where the expansion call is marginal, i.e. at least
# one tolerant species has EXACTLY twice the copies of the larger
# cold-sensitive count (a single background copy flips the call)
marginal <- copies %>%
  mutate(bg = pmax(Solanum_lycopersicum_total, Capsicum_annuum_total),
         marginal = Solanum_pennellii_total == 2 * bg |
                    Solanum_tuberosum_total == 2 * bg)
cat("Marginal expansion calls [X]:", sum(marginal$marginal), "of 76\n")
