#Expansion analysis

library(dplyr)
library(purrr)
library(tibble)
library(stringr)
load("data/A2TEA_finished_potatofix2.RData")   # HOG_level_list, hypotheses

species <- c("Solanum_pennellii", "Solanum_tuberosum", "Lycium_barbarum",
             "Solanum_lycopersicum", "Capsicum_annuum")

conserved_set <- HOG_level_list$all_species_overview %>%
  filter(if_all(ends_with("_total"), ~ . > 0))

M <- nrow(conserved_set)   # 14996, the conserved background

## Generalized expansion filter: min(expanded_in) >= factor * max(compared_to)
expanded_hogs <- function(data, expanded_in, compared_to, factor = 2) {
  stopifnot(length(intersect(expanded_in, compared_to)) == 0)
  e_min <- Reduce(pmin, data[paste0(expanded_in, "_total")])
  c_max <- Reduce(pmax, data[paste0(compared_to, "_total")])
  data[e_min >= factor * c_max, ]
}

## All 2-vs-2 (30) and 3-vs-2 (10) combinations
make_combos <- function(species, n_exp, n_cmp) {
  combn(species, n_exp, simplify = FALSE) |>
    lapply(function(E) {
      rest <- setdiff(species, E)
      combn(rest, n_cmp, simplify = FALSE) |>
        lapply(function(C) list(expanded_in = E, compared_to = C))
    }) |>
    unlist(recursive = FALSE)
}

combos <- c(make_combos(species, 2, 2),
            make_combos(species, 3, 2))   # 40

abbr <- c(Solanum_pennellii = "Spen", Solanum_tuberosum = "Stub",
          Lycium_barbarum   = "Lbar", Solanum_lycopersicum = "Slyc",
          Capsicum_annuum   = "Cann")
lbl <- function(v) paste(abbr[v], collapse = "+")

## Count table 
expansion_summary <- map_dfr(combos, function(cmb) {
  n <- expanded_hogs(conserved_set, cmb$expanded_in, cmb$compared_to) |> nrow()
  tibble(expanded_in = lbl(cmb$expanded_in), compared_to = lbl(cmb$compared_to),
         n_exp = length(cmb$expanded_in), n_cmp = length(cmb$compared_to),
         n_hogs = n)
}) |> arrange(desc(n_hogs))

## Exact hypergeometric test for the 2-set combinations
two_set <- map_dfr(keep(combos, ~ length(.x$expanded_in) == 2), function(cmb) {
  marg <- lapply(cmb$expanded_in,
                 \(e) expanded_hogs(conserved_set, e, cmb$compared_to) |> pull(HOG))
  K <- length(marg[[1]]); n <- length(marg[[2]])
  k <- length(intersect(marg[[1]], marg[[2]]))
  tibble(expanded_in = lbl(cmb$expanded_in), compared_to = lbl(cmb$compared_to),
         n_exp = 2, observed = k, expected = K * n / M, fold = k / (K * n / M),
         p = phyper(k - 1, K, M - K, n, lower.tail = FALSE))
})

## Randomisation test only for the 3-set combinations, 1e6 resamples
n_perm <- 100  # this is just a debugging number 1mio was used for my thesis
set.seed(1)
three_set <- map_dfr(keep(combos, ~ length(.x$expanded_in) == 3), function(cmb) {
  marg  <- lapply(cmb$expanded_in,
                  \(e) expanded_hogs(conserved_set, e, cmb$compared_to) |> pull(HOG))
  sizes <- lengths(marg)
  obs   <- length(Reduce(intersect, marg))
  perm  <- replicate(n_perm, {
    draws <- lapply(sizes, \(s) sample.int(M, s))
    length(Reduce(intersect, draws))
  })
  tibble(expanded_in = lbl(cmb$expanded_in), compared_to = lbl(cmb$compared_to),
         n_exp = 3, observed = obs, expected = mean(perm), fold = obs / mean(perm),
         p = (sum(perm >= obs) + 1) / (n_perm + 1))
})

overlap_results <- bind_rows(two_set, three_set) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(desc(fold))
## ---- single-species expanded sets (hypotheses 1, 2, 3, 8, 9) ----
background <- conserved_set$HOG
expanded_of <- function(h) {
  HOG_level_list[[h]] %>% filter(HOG %in% background, expansion == "yes") %>% pull(HOG)
}
sets <- list(pennellii    = expanded_of("hypothesis_1"),
             potato       = expanded_of("hypothesis_2"),
             goji         = expanded_of("hypothesis_3"),
             lycopersicum = expanded_of("hypothesis_8"),
             capsicum     = expanded_of("hypothesis_9"))


## ---- pairwise overlap: OR/CI from Fisher, one-sided p from phyper ----
pair_test <- function(a, b) {
  A <- sets[[a]]; B <- sets[[b]]
  K <- length(A); n <- length(B); k <- length(intersect(A, B))
  ft <- fisher.test(matrix(c(k, K - k, n - k, M - K - n + k), nrow = 2))  # two-sided: OR + 95% CI
  tibble(set_A = a, set_B = b, K = K, n = n, observed = k,
         expected = K * n / M, fold = k / (K * n / M),
         odds_ratio = unname(ft$estimate),
         ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
         p_enrich  = phyper(k - 1, K, M - K, n, lower.tail = FALSE),
         p_deplete = phyper(k, K, M - K, n))
}


## shared expansions as fraction of each species-specific set
tibble(pair    = c("Spen∩Stub", "Spen∩Stub", "Slyc∩Cann", "Slyc∩Cann"),
       species = c("pennellii", "potato", "lycopersicum", "capsicum"),
       shared  = c(rep(length(intersect(sets$pennellii, sets$potato)), 2),
                   rep(length(intersect(sets$lycopersicum, sets$capsicum)), 2))) %>%
  mutate(set_size = lengths(sets)[species],
         pct      = round(100 * shared / set_size, 1))

## the three tolerant pairs, BH over the three
tol_pairs <- map_dfr(list(c("pennellii", "potato"), c("pennellii", "goji"),
                          c("potato", "goji")), ~ pair_test(.x[1], .x[2])) %>%
  mutate(p_adj = p.adjust(p_enrich, method = "BH"))

## sensitive pair (73 HOGs) and the four sensitive-vs-tolerant controls
control <- map_dfr(list(c("lycopersicum", "capsicum"),
                        c("lycopersicum", "pennellii"), c("lycopersicum", "potato"),
                        c("capsicum", "pennellii"),     c("capsicum", "potato")),
                   ~ pair_test(.x[1], .x[2]))

overlap_results %>%                                   # cann expanded vs not, 2-set and 3-set
  mutate(cann_expanded = str_detect(expanded_in, "Cann")) %>%
  group_by(n_exp, cann_expanded) %>%
  summarise(median_fold = median(fold), n = n(), .groups = "drop")

overlap_results %>%                                   # Capsicum anywhere in the combo
  filter(n_exp == 2) %>%
  mutate(has_capsicum = str_detect(paste(expanded_in, compared_to), "Cann")) %>%
  group_by(has_capsicum) %>% summarise(median_fold = median(fold), n = n())

ov <- HOG_level_list$all_species_overview
colSums(ov[grepl("_total$", names(ov))])              # total gene counts per species
summary(ov$Capsicum_annuum_total)
summary(ov$Solanum_lycopersicum_total)


#visualising all 2vs2 pairs:
library(dplyr); library(stringr); library(purrr); library(ggplot2)

## tolerance status of each species
tol <- c(Spen = "T", Stub = "T", Lbar = "T", Slyc = "S", Cann = "S")

tol_pattern <- function(exp_pair, cmp_pair) {
  e  <- paste(sort(tol[str_split(exp_pair, "\\+")[[1]]]), collapse = "")
  cm <- paste(sort(tol[str_split(cmp_pair, "\\+")[[1]]]), collapse = "")
  if (e == "TT" && cm == "SS")      "Tolerant vs sensitive"
  else if (e == "SS" && cm == "TT") "Sensitive vs tolerant"
  else                              "Mixed (crosses tolerance)"
}

viz2 <- overlap_results %>%
  filter(n_exp == 2) %>%
  mutate(tol_class = map2_chr(expanded_in, compared_to, tol_pattern))

p_heat <- ggplot(viz2, aes(compared_to, expanded_in)) +
  geom_tile(aes(fill = fold), colour = "white", linewidth = 0.4) +
  ## outline only the tolerance-coherent combinations
  geom_tile(data = ~ filter(.x, tol_class != "Mixed (crosses tolerance)"),
            fill = NA, aes(colour = tol_class), linewidth = 1.3) +
  geom_text(aes(label = observed), size = 3) +
  scale_fill_viridis_c(trans = "log10", option = "magma", direction = -1,
                       name = "Fold") +
  scale_colour_manual(values = c("Tolerant vs sensitive" = "#1b9e77",
                                 "Sensitive vs tolerant" = "#7570b3"),
                      name = "Tolerance-coherent") +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/coexpansion_heatmap.pdf", p_heat, width = 8, height = 7)