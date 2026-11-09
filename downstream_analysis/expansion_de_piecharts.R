
# Nested "background vs. expanded" DE pie charts



suppressMessages({ library(stringr); library(ggplot2); library(patchwork); library(ggtext)})


RDATA_IN           <- "data/A2TEA_finished_potatofix2.RData"

RESTRICT_CONSERVED <- TRUE
SELECTORS          <- c("any_all", "any_expanded", "both_expanded")  # used only if PANELS = NULL
PANELS <- NULL


## ---------- function for creating the plots (x-axis title suppressed) -------------
create_custom_plot <- function(background_set, background_subset, interest_set, interest_subset,
                               axis_name = "interest (inside) vs. background (outside)") {
  p <- ggplot() +
    geom_bar(aes(x=log10(length(background_set))/2,
                 y=c(length(background_subset), (length(background_set)-length(background_subset)))),
             width = log10(length(background_set)),
             position="fill", stat="identity",
             fill=c("#F39C12", "#428bca"), color="white", alpha=0.5)+
    geom_bar(aes(x=log10(length(interest_set))/2,
                 y=c(length(interest_subset), length(interest_set)-length(interest_subset))),
             width = log10(length(interest_set)),
             position="fill", stat="identity",
             fill=c("#F39C12", "#428bca"), color="white", alpha=0.6)+
    coord_polar(theta = "y", start=0)+
    scale_y_continuous(name = NULL, breaks = NULL, expand = c(0,0))+
    scale_x_continuous(name = NULL, breaks = seq(0,0,0),
                       expand = expansion(mult = c(0, 0.08)))+
    labs(caption = axis_name)+
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill="transparent", colour=NA),
          panel.grid  = element_blank(),
          plot.caption = ggtext::element_markdown(size = 14, hjust = 0.5,
                                                  margin = margin(t = 2)))+
    geom_text(aes(x=2.7, y=0.1,
                  label=paste(length(background_subset), "/", length(background_set), sep = " ")),
              size=6) +
    geom_text(aes(x=1.5, y=0.1,
                  label=paste(length(interest_subset), "/", length(interest_set), sep = " ")),
              size=6)
  return(p)
}
## ---------------------------------------------------------------------------

abbrev <- function(s) sub("^([A-Za-z])[a-z]+_", "\\1. ", s)

load(RDATA_IN)
ov <- HOG_level_list[["all_species_overview"]]

background <- HOG_level_list$all_species_overview %>%
  dplyr::filter(if_all(ends_with("_total"), ~ . > 0)) %>% pull(HOG)

total_cols <- grep("_total$", names(ov), value = TRUE)
if (RESTRICT_CONSERVED) {
  ov <- ov[apply(as.matrix(ov[, total_cols]) > 0, 1, all), ]
  message(sprintf("Conserved-set restriction: %d HOGs kept.", nrow(ov)))
}
sig_species <- sub("_sigDE$", "", grep("_sigDE$", names(ov), value = TRUE))
message("Species WITH expression (DE) data: ", paste(sig_species, collapse = ", "))

deg_vec <- function(sp) { c <- paste0(sp,"_sigDE"); if (!c %in% names(ov)) NULL else ov[[c]] > 0 }
group_deg <- function(sps, mode) {            # mode = "any" | "all"
  mats <- lapply(sps, deg_vec)
  if (any(vapply(mats, is.null, logical(1)))) return(NULL)
  M <- do.call(cbind, mats)
  if (mode == "all") rowSums(M) == length(sps) else rowSums(M) >= 1
}
exp_hogs <- function(h) {
  tb <- HOG_level_list[[paste0("hypothesis_", h)]]
  if (is.null(tb) || !"expansion" %in% names(tb)) character(0) else tb$HOG[tb$expansion == "yes"]
}
compute_panel <- function(h, deg_species, deg_mode) {
  degm <- group_deg(deg_species, deg_mode); if (is.null(degm)) return(NULL)
  is_exp <- ov$HOG %in% exp_hogs(h)
  list(bg_total=nrow(ov), bg_deg=sum(degm), exp_total=sum(is_exp), exp_deg=sum(degm & is_exp))
}

message("\nHypotheses found:")
print(hypotheses[, intersect(c("hypothesis","name","expanded_in","compared_to"), names(hypotheses))])

if (is.null(PANELS)) {
  PANELS <- list(); L <- 0
  for (h in sort(unique(as.integer(hypotheses$hypothesis)))) {
    E    <- str_trim(str_split(hypotheses$expanded_in[hypotheses$hypothesis == h], ";")[[1]])
    E_de <- E[E %in% sig_species]
    if (length(E_de) == 0) { message(sprintf("hyp %d: no expanded DE species - skipped", h)); next }
    exp_desc <- paste(abbrev(E_de), collapse = " + ")
    for (sel in SELECTORS) {
      if (sel == "both_expanded" && length(E_de) < 2) next
      spec <- switch(sel, any_all=sig_species, any_expanded=E_de, both_expanded=E_de)
      mode <- if (sel == "both_expanded") "all" else "any"
      tag  <- switch(sel, any_all="DE: any species", any_expanded="DE: any expanded",
                     both_expanded="DE: all expanded")
      L <- L + 1
      PANELS[[L]] <- list(letter=letters[L], label=paste0(exp_desc, " — ", tag),
                          hypothesis=h, deg_species=spec, deg_mode=mode, italic=FALSE)
    }
  }
}


render_figure <- function(PANELS, out_pdf) {
  rows <- list(); keep <- list()
  for (p in PANELS) {
    r <- compute_panel(p$hypothesis, p$deg_species, p$deg_mode)
    if (is.null(r)) { message(sprintf("panel %s skipped (selector species lacks DE)", p$letter)); next }
    if (r$exp_total == 0) { message(sprintf("panel %s skipped (no expanded HOGs)", p$letter)); next }
    p$nums <- r; keep[[length(keep)+1]] <- p
    
    ## one-sided hypergeometric test for over-representation of DE HOGs
    ## in the expanded set, against the  background set
    p_over <- phyper(r$exp_deg - 1, r$bg_deg, r$bg_total - r$bg_deg,
                     r$exp_total, lower.tail = FALSE)
    F
    rows[[length(rows)+1]] <- data.frame(
      panel       = p$letter,
      hyp         = p$hypothesis,
      deg_species = paste(p$deg_species, collapse = "|"),
      mode        = p$deg_mode,
      background  = sprintf("%d/%d", r$bg_deg, r$bg_total),
      bg_pct      = round(100 * r$bg_deg / r$bg_total, 1),
      expanded    = sprintf("%d/%d", r$exp_deg, r$exp_total),
      exp_pct     = round(100 * r$exp_deg / r$exp_total, 1),
      expected    = round(r$exp_total * r$bg_deg / r$bg_total, 1),
      p_over      = signif(p_over, 3)
    )
  }
  PANELS <- keep
  stats_tbl <- do.call(rbind, rows)
  message("\n== ", out_pdf, " ==")
  print(stats_tbl, row.names = FALSE)
  dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
  write.csv(stats_tbl, sub("\\.pdf$", "_stats.csv", out_pdf), row.names = FALSE)
  
  plots <- lapply(PANELS, function(p) {
    n <- p$nums
    create_custom_plot(
      background_set    = seq_len(n$bg_total),
      background_subset = seq_len(n$bg_deg),
      interest_set      = seq_len(n$exp_total),
      interest_subset   = seq_len(n$exp_deg),
      axis_name         = p$label)
  })
  ncol <- min(4, length(plots)); nrow <- ceiling(length(plots) / ncol)
  combined <- wrap_plots(plots, ncol = ncol) +
    plot_annotation(tag_levels = "a",
                    #caption = "Outer ring = background HOGs; inner disc = expanded subset (radius \u221d log10 total). Orange = HOG with \u22651 DEG, blue = without."
    )
  ggsave(out_pdf, combined, width = ncol * 2.9 + 0.4, height = nrow * 3.2 + 0.5,
         device = cairo_pdf, bg = "white")
  message("Wrote: ", normalizePath(out_pdf))
}



# ---- Figure 1: tolerant expansions ----
panels_tolerant <- list(
  list(letter="a", label="*S. pennellii*", hypothesis=1,
       deg_species="Solanum_pennellii", deg_mode="any"),
  list(letter="b", label="*S. tuberosum*", hypothesis=2,
       deg_species="Solanum_tuberosum", deg_mode="any"),
  list(letter="c", label="*S. pennellii* and *S. tuberosum*", hypothesis=4,
       deg_species=c("Solanum_pennellii","Solanum_tuberosum"), deg_mode="any"),
  list(letter="d", label="*S. pennellii* and *S. tuberosum*", hypothesis=4,
       deg_species=c("Solanum_pennellii","Solanum_tuberosum"), deg_mode="all"),
  list(letter="e", label="all tolerant", hypothesis=7,
       deg_species=c("Solanum_pennellii","Solanum_tuberosum"), deg_mode="any"),
  list(letter="f", label="all tolerant", hypothesis=7,
       deg_species=c("Solanum_pennellii","Solanum_tuberosum"), deg_mode="all")
)

render_figure(panels_tolerant, "figures/piecharts_tolerant_expandedDE.pdf")

# ---- Figure 2: sensitive expansions ----
panels_sensitive <- list(
  list(letter="a", label="*S. lycopersicum*", hypothesis=8,
       deg_species="Solanum_lycopersicum", deg_mode="any"),
  list(letter="b", label="*C. annuum*", hypothesis=9,
       deg_species="Capsicum_annuum", deg_mode="any"),
  list(letter="c", label="*S. lycopersicum* and *C. annuum*", hypothesis=10,
       deg_species=c("Solanum_lycopersicum","Capsicum_annuum"), deg_mode="any"),
  list(letter="d", label="*S. lycopersicum* and *C. annuum*", hypothesis=10,
       deg_species=c("Solanum_lycopersicum","Capsicum_annuum"), deg_mode="all")
)
render_figure(panels_sensitive, "figures/piecharts_sensitive_expandedDE.pdf")

# ---- Figure 3: tolerant expansions - DE from any of all the 4 species with DE-data ----
panels_tolerant_any <- list(
  list(letter="a", label="*S. pennellii*", hypothesis=1,
       deg_species=sig_species, deg_mode="any"),
  list(letter="b", label="*S. tuberosum*", hypothesis=2,
       deg_species=sig_species, deg_mode="any"),
  list(letter="c", label="*S. pennellii* and *S. tuberosum*", hypothesis=4,
       deg_species=sig_species, deg_mode="any"),
  list(letter="d", label="*S. pennellii* and *S. tuberosum*", hypothesis=4,
       deg_species=sig_species, deg_mode="all"),
  list(letter="e", label="all tolerant", hypothesis=7,
       deg_species=sig_species, deg_mode="any"),
  list(letter="f", label="all tolerant", hypothesis=7,
       deg_species=sig_species, deg_mode="all")
)

render_figure(panels_tolerant_any, "figures/piecharts_tolerant_anyDE.pdf")

# ---- Figure 2: sensitive expansions - DE from any of all the 4 species with DE-data----
panels_sensitive_any <- list(
  list(letter="a", label="*S. lycopersicum*", hypothesis=8,
       deg_species=sig_species, deg_mode="any"),
  list(letter="b", label="*C. annuum*", hypothesis=9,
       deg_species=sig_species, deg_mode="any"),
  list(letter="c", label="*S. lycopersicum* and *C. annuum*", hypothesis=10,
       deg_species=sig_species, deg_mode="any"),
  list(letter="d", label="*S. lycopersicum* and *C. annuum*", hypothesis=10,
       deg_species=sig_species, deg_mode="all")
)
render_figure(panels_sensitive_any, "figures/piecharts_sensitive_anyDE.pdf")

############################
#df with nr. of expansions per hyp:

library(dplyr)
library(purrr)

summary_df <- paste0("hypothesis_", 1:11) %>%
  set_names(1:11) %>%
  map_dfr(
    ~ {
      df <- HOG_level_list[[.x]]
      tibble(
        n_expanded  = sum(df$HOG %in% background & df$expansion == "yes")
      )
    },
    .id = "hypothesis"
  ) %>%
  mutate(hypothesis = as.numeric(hypothesis))

result <- hypotheses %>%
  dplyr::select(hypothesis, expanded_in, compared_to) %>%
  left_join(summary_df, by = "hypothesis")

library(kableExtra)

fmt_species <- function(x) {
  x %>%
    str_split(";") %>%
    map_chr(~ paste(sprintf("\\textit{%s}", str_replace_all(.x, "_", " ")),
                    collapse = ", "))
}

result_fmt <- result %>%
  mutate(
    expanded_in = fmt_species(expanded_in),
    compared_to = fmt_species(compared_to)
  )

tbl <- result_fmt %>%
  kbl(
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    col.names = c("Hyp.", "Expanded in", "Compared to", "Expanded"),
    caption   = "Number of expanded HOGs per hypothesis.",
    label     = "hog-counts",
    align     = "cllr"
  ) %>%
  kable_styling(latex_options = "hold_position")

writeLines(tbl, "tables/hog_counts.tex")



library(dplyr); library(purrr); library(tibble)

M <- length(background)

expanded_of <- function(h) {
  HOG_level_list[[h]] %>%
    dplyr::filter(HOG %in% background, expansion == "yes") %>%
    pull(HOG)
}

sets <- list(
  pennellii = expanded_of("hypothesis_1"),
  potato    = expanded_of("hypothesis_2"),
  goji      = expanded_of("hypothesis_3")
)

overlap_test <- function(A, B, M) {
  K <- length(A); n <- length(B); k <- length(intersect(A, B))
  expected <- K * n / M
  tibble(
    K = K, n = n, observed = k, expected = expected,
    fold = k / expected,
    p_enrichment = phyper(k - 1, K, M - K, n, lower.tail = FALSE), # P(X >= k)
    p_depletion  = phyper(k,     K, M - K, n)                       # P(X <= k)
  )
}

overlap_test(sets$pennellii, sets$potato, M)
overlap_test(sets$pennellii, sets$goji, M)
overlap_test(sets$potato, sets$goji, M)


pairs <- list(
  c("pennellii", "potato"),
  c("pennellii", "goji"),
  c("potato",    "goji")
)

fisher_pair <- function(a, b, M) {
  A <- sets[[a]]; B <- sets[[b]]
  K <- length(A); n <- length(B); k <- length(intersect(A, B))
  # 2x2: both / A-only / B-only / neither
  ft <- fisher.test(matrix(c(k, K - k, n - k, M - K - n + k), nrow = 2))
  tibble(
    set_A = a, set_B = b,
    K = K, n = n, observed = k, expected = K * n / M,
    odds_ratio = unname(ft$estimate),
    ci_low  = ft$conf.int[1],
    ci_high = ft$conf.int[2],
    p_value = ft$p.value
  )
}

fisher_all <- map_dfr(pairs, ~ fisher_pair(.x[1], .x[2], M)) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

fisher_all


set.seed(1)
n_perm <- 10000
grp    <- c("pennellii", "potato", "goji")

sizes <- lengths(sets[grp])
obs   <- length(Reduce(intersect, sets[grp]))

perm <- replicate(n_perm, {
  draws <- lapply(sizes, function(s) sample.int(M, s))
  length(Reduce(intersect, draws))
})

p_enrich  <- (sum(perm >= obs) + 1) / (n_perm + 1)
p_deplete <- (sum(perm <= obs) + 1) / (n_perm + 1)
c(observed = obs, expected = mean(perm), p_enrich = p_enrich, p_deplete = p_deplete)


###########
#sensitive control

## Sensitive-species expanded sets, restricted to the conserved background
sets$lycopersicum <- expanded_of("hypothesis_8")
sets$capsicum     <- expanded_of("hypothesis_9")



overlap_test(sets$lycopersicum, sets$capsicum, M)
fisher_pair("lycopersicum", "capsicum", M)

phyper(164 - 1, 270, M - 270, 897, lower.tail = FALSE)

set.seed(1)
n_perm <- 10000

perm_overlap_test <- function(members, M, n_perm = 10000) {
  sizes <- lengths(sets[members])
  obs   <- length(Reduce(intersect, sets[members]))
  perm  <- replicate(n_perm, {
    draws <- lapply(sizes, function(s) sample.int(M, s))
    length(Reduce(intersect, draws))
  })
  tibble(
    sets      = paste(members, collapse = " ∩ "),
    observed  = obs,
    expected  = mean(perm),
    fold      = obs / mean(perm),
    p_enrich  = (sum(perm >= obs) + 1) / (n_perm + 1),
    p_deplete = (sum(perm <= obs) + 1) / (n_perm + 1)
  )
}

pairs3 <- list(
  c("pennellii", "potato"),
  c("pennellii", "goji"),
  c("potato",    "goji")
)

perm_pairs <- purrr::map_dfr(pairs3, ~ perm_overlap_test(.x, M)) %>%
  mutate(p_adj = p.adjust(p_enrich, method = "BH"))

perm_pairs

set.seed(1)

perm_sensitive <- perm_overlap_test(c("lycopersicum", "capsicum"), M)
perm_sensitive

perm_overlap_test(c("lycopersicum", "pennellii"), M)   # expect observed = 0
perm_overlap_test(c("capsicum",     "potato"),    M)   # expect observed = 0
