
# GO term enrichment analysis (topGO) — semi-automated over HOG subsets

suppressPackageStartupMessages({
  library(topGO)
  library(tidyverse)
})


## 0. LOAD DATA  &  SETTINGS

## A2TEA results — load the potato-fix version so tuberosum DE is correct.
load("data/A2TEA_finished_potatofix2.RData")

outdir     <- "data/GOEnrichment_potatofix"   # created if missing
ontologies <- c("BP", "MF", "CC")        # which ontologies to test

P_CUTOFF        <- 0.05   # one-sided Fisher p <= this               (paper)
MIN_SIGNIFICANT <- 1      # keep terms with > this many sig. HOGs, i.e. >= 2
                          # ("annotated at least twice in the set")   (paper)
NODE_SIZE       <- 1      # min. HOGs annotated to a GO term. topGO default = 1
                          # (your original run); 5-10 prunes the smallest terms.

WRITE_REVIGO <- TRUE      # also write 2-column (GO.ID, p) files for REVIGO


## 1. Full HOG -> GO term map  (the annotation; built once)
## GO annotation per HOG = UNION of the GO IDs of every gene in the orthogroup,
## across all species (bind_rows(SFA)). Same cleaning as your original: collapse
## per HOG, strip the "NA" placeholders, drop HOGs left with no real GO term.
## (Assumes the column holds GO IDs, so "NA" only ever appears as a placeholder.)

hog2go_tbl <- bind_rows(SFA) %>%
  dplyr::select(HOG, `Gene-Ontology-Term`) %>%
  group_by(HOG) %>%
  summarise(`Gene-Ontology-Term` = paste(unique(`Gene-Ontology-Term`),
                                          collapse = ", "),
            .groups = "drop") %>%
  mutate(`Gene-Ontology-Term` = str_remove_all(`Gene-Ontology-Term`,
                                               "NA, |, NA")) %>%
  filter(!str_detect(`Gene-Ontology-Term`, "NA"))

# named list:  HOG  ->  character vector of GO IDs
hog2go_map <- setNames(hog2go_tbl$`Gene-Ontology-Term`, hog2go_tbl$HOG) %>%
  strsplit(split = ", ")


## 2. BACKGROUND  &  HOG SETS

overview <- HOG_level_list$all_species_overview

## Background distribution — the Conserved Set: HOGs present in every species
## (all `_total` counts > 0). Used as the background for every test (paper).
background <- conserved_HOGS <- overview %>%
  dplyr::filter(if_all(ends_with("_total"), ~ . > 0)) %>%
  pull(HOG)

## Expanded HOGs, restricted to the Conserved Set. hypothesis_1 = S. pennellii,
## hypothesis_2 = S. tuberosum (expanded_both could equally be pulled from
## hypothesis_4 — logically the same result).
expanded_penn <- HOG_level_list$hypothesis_1 %>%
  dplyr::filter(HOG %in% background, expansion == "yes") %>% pull(HOG)
expanded_tub  <- HOG_level_list$hypothesis_2 %>%
  dplyr::filter(HOG %in% background, expansion == "yes") %>% pull(HOG)
expanded_both <- intersect(expanded_penn, expanded_tub)

## Expanded in ALL THREE tolerant species (incl. Lycium barbarum), hypothesis_7.
## There is no expression data for L. barbarum, so this set can only ever be
## filtered by penn/tub (or sensitive-species) DE — the Lycium contribution is
## the expansion evidence alone.
expanded_all_tol <- HOG_level_list$hypothesis_7 %>%
  dplyr::filter(HOG %in% background, expansion == "yes") %>% pull(HOG)

## DE definition. A HOG counts as DE in a species when that species contributes
## >= 1 significant DEG to the HOG, i.e. its `_sigDE` count is > 0 in the
## overview — the same definition used for the DEG pie charts.
de_penn <- overview %>% dplyr::filter(Solanum_pennellii_sigDE > 0) %>% pull(HOG)
de_tub  <- overview %>% dplyr::filter(Solanum_tuberosum_sigDE  > 0) %>% pull(HOG)

de_any_tol  <- union(de_penn, de_tub)      # DE in EITHER tolerant species (penn OR tub)
de_both_tol <- intersect(de_penn, de_tub)  # DE in BOTH tolerant species (penn AND tub)

## DE in ANY species with expression data (penn, tub, lycopersicum, capsicum).
## Rationale: DE is used here as evidence that a gene family responds to cold at
## all. A HOG expanded in the tolerant species that is DE only in a sensitive
## species is still cold-responsive, so it stays a candidate under this wider
## definition. Picks up every `_sigDE` column present, so it adapts if the set
## of species with expression data changes.
de_any_species <- overview %>%
  dplyr::filter(if_any(ends_with("_sigDE"), ~ . > 0)) %>% pull(HOG)

## Species-exclusive DE within the tolerant pair. The penn and tub cold
## experiments differ in severity and duration (short severe near-freezing vs.
## prolonged mild), so the two halves of the union may reflect different
## physiology. Testing them separately makes that visible instead of confounded.
de_penn_only <- setdiff(de_penn, de_tub)
de_tub_only  <- setdiff(de_tub,  de_penn)

# Note: `_sigDE` counts carry no direction. For up/down-specific sets, build
# de_penn / de_tub from the per-species regulation categories instead.


## 3. JOB TABLE  — each entry is one set of interest (foreground)

## The background is always the Conserved Set, so each job is just the
## foreground HOG set. Any foreground HOG that is not in the Conserved Set (or
## carries no GO annotation) is necessarily excluded from its test — the
## function reports how many of each foreground actually enter the analysis.

jobs <- list(
  # --- primary candidate sets: expanded in BOTH tolerant species -------------
  expanded_penn_tub           = expanded_both,
  # DE in either tolerant species — the direct analogue of the reference
  # paper's 93 candidate HOGs
  expanded_penn_tub_DEanyTol  = intersect(expanded_both, de_any_tol),
  # DE in ANY species with expression data (incl. the sensitive references)
  expanded_penn_tub_DEany     = intersect(expanded_both, de_any_species),
  # DE in both tolerant species (focus subset; may be very small)
  expanded_penn_tub_DEbothTol = intersect(expanded_both, de_both_tol),

  # --- the two halves of the tolerant union, tested separately --------------
  # Splits the candidate set by which cold experiment contributed the DE
  expanded_penn_tub_DEpennOnly = intersect(expanded_both, de_penn_only),
  expanded_penn_tub_DEtubOnly  = intersect(expanded_both, de_tub_only),

  # --- expanded in all three tolerant species (incl. L. barbarum) -----------
  # Expansion evidence from three species, DE evidence still only from penn/tub
  expanded_all_tol             = expanded_all_tol,
  expanded_all_tol_DEanyTol    = intersect(expanded_all_tol, de_any_tol),
  expanded_all_tol_DEany       = intersect(expanded_all_tol, de_any_species),

  # --- backgrounds for comparison: conserved HOGs by DE status --------------
  conserved_DEanyTol          = de_any_tol,
  conserved_DEbothTol         = de_both_tol,
  conserved_DEany             = de_any_species
)


## 3b. SET SIZES  &  DE PROVENANCE  (supplementary tables)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## How many HOGs are in each foreground, and how many survive the two filters
## every test applies (must be in the Conserved Set, must carry a GO term).
set_sizes <- tibble::tibble(
  subset   = names(jobs),
  n_total  = purrr::map_int(jobs, ~ length(unique(.x))),
  n_in_bg  = purrr::map_int(jobs, ~ length(intersect(unique(.x), background))),
  n_tested = purrr::map_int(jobs, ~ length(intersect(intersect(unique(.x), background),
                                                     names(hog2go_map))))
)
print(set_sizes, n = Inf)
write.csv(set_sizes, row.names = FALSE,
          file = file.path(outdir, "GOenrichment_set_sizes.csv"))

## Which cold experiment puts each candidate HOG into the tolerant DE union.
## The penn and tub experiments differ in DEG yield, so this shows whether the
## union is balanced or effectively driven by one of the two.
cand <- intersect(expanded_both, de_any_tol)
de_provenance <- tibble::tibble(
  source = c("S. pennellii only", "S. tuberosum only", "both tolerant species"),
  n_HOGs = c(length(intersect(cand, de_penn_only)),
             length(intersect(cand, de_tub_only)),
             length(intersect(cand, de_both_tol)))
) %>% mutate(percent = round(100 * n_HOGs / sum(n_HOGs), 1))
print(de_provenance)
write.csv(de_provenance, row.names = FALSE,
          file = file.path(outdir, "GOenrichment_DE_provenance.csv"))


## 4. Enrichment function  (one foreground set, one ontology)

run_go_enrichment <- function(gene2go, foreground, universe = names(gene2go),
                              ontology = "BP", node_size = NODE_SIZE,
                              p_cutoff = P_CUTOFF,
                              min_significant = MIN_SIGNIFICANT) {

  # restrict the annotation to the background distribution
  gene2go <- gene2go[names(gene2go) %in% universe]
  if (length(gene2go) == 0) {
    warning("Empty background — skipping."); return(NULL)
  }

  # membership vector over the background: 1 = foreground, 0 = rest
  in_genes <- factor(as.integer(names(gene2go) %in% foreground))
  names(in_genes) <- names(gene2go)

  n_fg <- sum(in_genes == 1)
  message(sprintf("    %d of %d foreground HOGs enter the test (in background & annotated)",
                  n_fg, length(unique(foreground))))
  if (n_fg == 0) {
    warning("No foreground HOGs in the background — skipping."); return(NULL)
  }

  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = in_genes,
                nodeSize = node_size,
                annot    = annFUN.gene2GO,
                gene2GO  = gene2go)

  # one-sided Fisher's exact test (topGO tests over-representation)
  result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

  # note: "classicFisher" is only the column label; the p-values are from the
  # weight01 algorithm run above
  res <- GenTable(GOdata,
                  classicFisher = result,
                  orderBy  = "classicFisher",          # must match the arg name
                  ranksOf  = "classicFisher",
                  topNodes = length(usedGO(GOdata)),   # all tested terms, no cut
                  numChar  = 1000)

  # topGO prints very small p-values as "< 1e-30"; strip non-numeric chars
  # before coercing, otherwise those (most significant!) terms become NA.
  res$classicFisher <- as.numeric(gsub("[^0-9.eE+-]", "", res$classicFisher))

  # over-represented: p <= 0.05 AND term annotated to >= 2 HOGs of interest
  res[res$Significant > min_significant & res$classicFisher <= p_cutoff, ]
}


## 5. Run everything  &  write CSVs

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if (WRITE_REVIGO) dir.create(file.path(outdir, "revigo"),
                             showWarnings = FALSE, recursive = TRUE)

run_grid <- tidyr::crossing(job_name = names(jobs), ontology = ontologies)

all_results <- purrr::pmap(run_grid, function(job_name, ontology) {

  message(sprintf("  enrichment:  %-24s  %s", job_name, ontology))

  res <- run_go_enrichment(
    gene2go    = hog2go_map,
    foreground = jobs[[job_name]],
    universe   = background,          # Conserved Set, always
    ontology   = ontology
  )
  if (is.null(res) || nrow(res) == 0) return(NULL)

  # main per-subset / per-ontology table
  write.csv(res, row.names = FALSE,
            file = file.path(outdir,
                             sprintf("GOenrichment_%s_%s.csv", job_name, ontology)))

  # optional REVIGO input (GO ID + p-value)
  if (WRITE_REVIGO) {
    write.csv(res[, c("GO.ID", "classicFisher")], row.names = FALSE,
              file = file.path(outdir, "revigo",
                               sprintf("REVIGO_%s_%s.csv", job_name, ontology)))
  }

  tibble::as_tibble(res) %>%
    mutate(subset = job_name, ontology = ontology, .before = 1)
}) %>%
  bind_rows()

# one combined master table across all subsets/ontologies (single supp. file)
write.csv(all_results, row.names = FALSE,
          file = file.path(outdir, "GOenrichment_ALL_combined.csv"))

message("\nDone. ", nrow(all_results), " enriched GO terms across ",
        length(jobs), " subsets written to '", outdir, "'.")


## 6. PLOTTING & TABLES  — settings

indir   <- "data/GOEnrichment_potatofix"                  # output of sections 0-5
figdir  <- file.path(indir, "figures")
tabdir  <- file.path(indir, "tables")

## The primary candidate set (matches the reference paper's definition:
## expanded in both tolerant species, DE in at least one tolerant species).
PRIMARY   <- "expanded_penn_tub_DEanyTol"
## Background shown alongside the candidate sets in Figure 2.
BACKGROUND <- "conserved_DEanyTol"

## Subsets shown in Figure 2, in plotting order (left to right).
CANDIDATE_SUBSETS <- c("expanded_penn_tub",
                       "expanded_all_tol",
                       "expanded_penn_tub_DEanyTol",
                       "expanded_all_tol_DEanyTol",
                       "expanded_penn_tub_DEpennOnly",
                       "expanded_penn_tub_DEtubOnly")

## Short display labels — keep axis text readable.
SUBSET_LABELS <- c(
  expanded_penn_tub            = "Exp. 2 tol.",
  expanded_all_tol             = "Exp. 3 tol.",
  expanded_penn_tub_DEanyTol   = "Exp. 2 tol.\n+ DE",
  expanded_all_tol_DEanyTol    = "Exp. 3 tol.\n+ DE",
  expanded_penn_tub_DEpennOnly = "+ DE penn\nonly",
  expanded_penn_tub_DEtubOnly  = "+ DE tub\nonly",
  conserved_DEanyTol           = "Background\n(cons. + DE)"
)

TERM_WRAP <- 45   # wrap long GO term labels at this many characters

dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
dir.create(tabdir, showWarnings = FALSE, recursive = TRUE)


## 7. Load the enrichment output

res <- read_csv(file.path(indir, "GOenrichment_ALL_combined.csv"),
                show_col_types = FALSE)

sizes <- read_csv(file.path(indir, "GOenrichment_set_sizes.csv"),
                  show_col_types = FALSE)

prov_file <- file.path(indir, "GOenrichment_DE_provenance.csv")
prov <- if (file.exists(prov_file)) read_csv(prov_file, show_col_types = FALSE) else NULL

## A subset that was tested but yielded nothing has no rows in `res`. Recover
## those from `sizes` so "tested, no enrichment" can be reported explicitly
## rather than silently missing from every figure and table.
tested_no_hits <- setdiff(sizes$subset, unique(res$subset))
if (length(tested_no_hits)) {
  message("Tested but no significant enrichment: ",
          paste(tested_no_hits, collapse = ", "))
}

## Helper: write a table as CSV and (if knitr is installed) as LaTeX.
write_table <- function(x, name, caption = NULL, digits = 4) {
  write.csv(x, row.names = FALSE, file = file.path(tabdir, paste0(name, ".csv")))
  if (requireNamespace("knitr", quietly = TRUE)) {
    tex <- knitr::kable(x, format = "latex", booktabs = TRUE,
                        digits = digits, caption = caption, linesep = "")
    writeLines(as.character(tex), file.path(tabdir, paste0(name, ".tex")))
  }
  invisible(x)
}


## 8. TABLE 2 — set sizes and DE rates
## Documents the nested series and shows what each filter costs. The DE rate
## column makes explicit whether a stricter expansion filter also selects for
## cold responsiveness (in this data it does not).

expansion_sets <- c("expanded_penn_tub", "expanded_all_tol")

tab_sizes <- sizes %>%
  mutate(
    enriched_terms = map_int(subset, ~ sum(res$subset == .x)),
    note = if_else(subset %in% tested_no_hits,
                   "tested, no significant enrichment", "")
  ) %>%
  dplyr::select(subset, n_total, n_in_bg, n_tested, enriched_terms, note)

write_table(tab_sizes, "Table2_set_sizes",
            caption = "HOG numbers per analysed subset. n_in_bg: within the Conserved Set; n_tested: additionally carrying at least one GO annotation.")

## DE rate within each expansion set (how many expanded HOGs are cold-responsive)
de_rate <- tibble(
  expansion_set = expansion_sets,
  n_expanded    = sizes$n_total[match(expansion_sets, sizes$subset)],
  n_DE_tolerant = sizes$n_total[match(paste0(expansion_sets, "_DEanyTol"), sizes$subset)]
) %>%
  mutate(percent_DE = round(100 * n_DE_tolerant / n_expanded, 1))

write_table(de_rate, "Table2b_DE_rates",
            caption = "Proportion of expanded HOGs differentially expressed in at least one tolerant species.")
print(de_rate)


## 9. TABLE 1 — the primary candidate set

tab_primary <- res %>%
  dplyr::filter(subset == PRIMARY) %>%
  arrange(ontology, classicFisher) %>%
  dplyr::select(Ontology = ontology, `GO ID` = GO.ID, Term,
                Annotated, Significant, Expected, p = classicFisher)

write_table(tab_primary, "Table1_primary_candidate_set",
            caption = paste0("GO terms over-represented in the primary candidate set (",
                             PRIMARY, ")."))
message("Table 1: ", nrow(tab_primary), " terms")

## Which candidate terms are absent from every conserved DE background — the
## central contrast with the reference study. Reported as a column, so the
## claim in the text is backed by the table itself. Note: "background" here is
## ANY of the three conserved_* subsets, whereas Figure 2 shows only BACKGROUND.
background_subsets <- grep("^conserved_", unique(res$subset), value = TRUE)
bg_terms <- res %>% dplyr::filter(subset %in% background_subsets) %>% pull(GO.ID) %>% unique()

tab_primary_flagged <- tab_primary %>%
  mutate(in_background = if_else(`GO ID` %in% bg_terms, "yes", "no"))

write_table(tab_primary_flagged, "Table1b_primary_vs_background",
            caption = "Primary candidate set terms, flagged by whether they are also enriched in any conserved DE background.")
message("  of these, absent from all DE backgrounds: ",
        sum(tab_primary_flagged$in_background == "no"))


## 10. FIGURE 1 — set sizes across the nested series

fig1_data <- tab_sizes %>%
  dplyr::filter(subset %in% c(CANDIDATE_SUBSETS, "expanded_penn_tub_DEbothTol")) %>%
  mutate(subset = fct_reorder(subset, n_total))

fig1 <- ggplot(fig1_data, aes(x = n_total, y = subset)) +
  geom_col(fill = "grey75", width = 0.7) +
  geom_text(aes(label = n_total), hjust = -0.25, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "HOGs", y = NULL,
       title = "Nested HOG subsets") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(figdir, "Figure1_set_sizes.pdf"), fig1,
       width = 6, height = 3.4)

## DE provenance — which cold experiment contributed the DE evidence
if (!is.null(prov)) {
  fig1b <- ggplot(prov, aes(x = n_HOGs, y = fct_reorder(source, n_HOGs))) +
    geom_col(fill = "grey55", width = 0.6) +
    geom_text(aes(label = paste0(n_HOGs, " (", percent, "%)")),
              hjust = -0.15, size = 3.2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = "candidate HOGs", y = NULL,
         title = "Source of the DE evidence") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank())
  ggsave(file.path(figdir, "Figure1b_DE_provenance.pdf"), fig1b,
         width = 5.5, height = 2.4)
}


## 11. FIGURE 2 — enriched terms across candidate subsets vs. background
## Terms are those enriched in ANY candidate subset; the background column
## shows which of them are also enriched in the conserved DE set. The visual
## claim is that most candidate terms have no background counterpart.

cand_terms <- res %>%
  dplyr::filter(subset %in% CANDIDATE_SUBSETS) %>%
  pull(GO.ID) %>% unique()

fig2_data <- res %>%
  dplyr::filter(subset %in% c(CANDIDATE_SUBSETS, BACKGROUND),
                GO.ID %in% cand_terms) %>%
  mutate(
    subset    = factor(subset, levels = c(CANDIDATE_SUBSETS, BACKGROUND)),
    ontology  = factor(ontology, levels = c("BP", "MF", "CC")),
    neglog_p  = -log10(classicFisher),
    Term_wrap = str_wrap(Term, TERM_WRAP)
  )

## order terms within each ontology by their best p-value across subsets
term_order <- fig2_data %>%
  group_by(ontology, Term_wrap) %>%
  summarise(best = min(classicFisher), .groups = "drop") %>%
  arrange(ontology, desc(best))

fig2_data <- fig2_data %>%
  mutate(Term_wrap = factor(Term_wrap, levels = term_order$Term_wrap))

fig2 <- ggplot(fig2_data, aes(x = subset, y = Term_wrap)) +
  geom_point(aes(size = Significant, colour = neglog_p)) +
  facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = SUBSET_LABELS) +
  scale_colour_viridis_c(option = "C", direction = -1,
                         name = expression(-log[10](p))) +
  scale_size_continuous(name = "HOGs in set", range = c(1.5, 6),
                        trans = "sqrt") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(size = 8),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

ggsave(file.path(figdir, "Figure2_enrichment_dotplot.pdf"), fig2,
       width = 8.5, height = 11)

message("Figure 2: ", n_distinct(fig2_data$GO.ID), " terms across ",
        n_distinct(fig2_data$subset), " columns")

## Note: `Significant` is on very different scales between candidate subsets
## (single digits) and the background (hundreds), so the size legend spans a
## wide range. sqrt scaling keeps the small sets visible.


## 12. TABLE 3 — candidate HOGs with functional annotation
## With several terms resting on only 2-4 HOGs, the named gene families carry
## more weight than the p-values. Needs SFA, hog2go_map and jobs from
## sections 0-3 in the session.
##
## SFA is a named list of per-species tibbles with columns
##   HOG, Protein-Accession, Blast-Hit-Accession,
##   Human-Readable-Description, Gene-Ontology-Term
## i.e. one row per GENE. Table 3 collapses that to one row per HOG x GO term
## (for the thesis), Table 3b keeps the gene-level detail (for the appendix).

if (exists("SFA") && exists("hog2go_map") && exists("jobs")) {

  sfa_all <- bind_rows(SFA, .id = "species")

  focal_terms    <- tab_primary$`GO ID`
  candidate_hogs <- jobs[[PRIMARY]]

  ## which candidate HOGs carry which of the enriched terms
  hog_terms <- tibble(HOG = names(hog2go_map), GO = hog2go_map) %>%
    dplyr::filter(HOG %in% candidate_hogs) %>%
    unnest(GO) %>%
    dplyr::filter(GO %in% focal_terms)

  term_info <- res %>%
    dplyr::filter(subset == PRIMARY) %>%
    dplyr::select(GO = GO.ID, Term, ontology)

  ## Table 3b: gene-level detail
  tab_hogs_genes <- hog_terms %>%
    left_join(sfa_all, by = "HOG", relationship = "many-to-many") %>%
    left_join(term_info, by = "GO") %>%
    arrange(ontology, Term, HOG, species) %>%
    dplyr::select(Ontology = ontology, Term, `GO ID` = GO, HOG, Species = species,
                  Gene = `Protein-Accession`,
                  Annotation = `Human-Readable-Description`)

  write_table(tab_hogs_genes, "Table3b_candidate_HOG_genes",
              caption = "Genes of the candidate HOGs underlying each enriched GO term, with AHRD annotation.")

  ## Table 3: one row per HOG x term, annotations collapsed.
  ## Uninformative AHRD strings are dropped before collapsing so the column
  ## shows real functions; if a HOG has only such hits it is kept and marked.
  uninformative <- regex("uncharacterized|unknown protein|hypothetical|LOW QUALITY",
                         ignore_case = TRUE)

  hog_annot <- sfa_all %>%
    dplyr::filter(HOG %in% candidate_hogs) %>%
    group_by(HOG) %>%
    summarise(
      n_genes    = n(),
      n_species  = n_distinct(species),
      Annotation = {
        a <- unique(na.omit(`Human-Readable-Description`))
        a <- a[!str_detect(a, uninformative)]
        if (length(a) == 0) "(uncharacterized)" else paste(a, collapse = "; ")
      },
      .groups = "drop"
    )

  tab_hogs <- hog_terms %>%
    left_join(term_info, by = "GO") %>%
    left_join(hog_annot, by = "HOG") %>%
    arrange(ontology, Term, HOG) %>%
    dplyr::select(Ontology = ontology, Term, `GO ID` = GO, HOG,
                  n_genes, n_species, Annotation)

  write_table(tab_hogs, "Table3_candidate_HOGs_annotation",
              caption = "HOGs of the primary candidate set underlying each enriched GO term, with AHRD functional annotation.")

  message("Table 3: ", n_distinct(tab_hogs$HOG), " HOGs across ",
          n_distinct(tab_hogs$`GO ID`), " terms (",
          nrow(tab_hogs_genes), " genes in Table 3b)")

  ## Table 3c: inverted view, one row per HOG.
  ## Shows how many of the enriched terms each HOG contributes to. Terms that
  ## are annotation-redundant (e.g. translation / ribosome / RNP complex on the
  ## same ribosomal-protein HOG) collapse onto one row here, so this is what
  ## tells you how many independent gene families are behind a signature.
  tab_hogs_wide <- tab_hogs %>%
    group_by(HOG, n_genes, n_species, Annotation) %>%
    summarise(n_terms = n_distinct(`GO ID`),
              Terms   = paste(sort(unique(Term)), collapse = "; "),
              .groups = "drop") %>%
    arrange(desc(n_terms), HOG) %>%
    dplyr::select(HOG, n_genes, n_species, n_terms, Terms, Annotation)

  write_table(tab_hogs_wide, "Table3c_candidate_HOGs_by_HOG",
              caption = "Candidate HOGs contributing to the enriched GO terms, with the number of terms each carries.")

} else {
  message("Table 3 skipped: SFA / hog2go_map / jobs not in the session ",
          "(run sections 0-3 first).")
}
