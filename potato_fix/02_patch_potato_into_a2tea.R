#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Splice corrected Solanum_tuberosum DE results into an existing A2TEA
# final .RData
#
# Touches ONLY DE-derived objects/columns:
#   - combined_AllSpeciesDEResultsDataFrames   (potato rows)
#   - HOG_DE.a2tea                             (rebuilt via the same join)
#   - HOG_level_list[[*]]  Solanum_tuberosum_sigDE, total_sigDE,
#                          norm_sum_*_sigDE, oaes_hyper, tea_value
#   - all_species_overview  Solanum_tuberosum_sigDE
# Leaves untouched: HYPOTHESES.a2tea, A2TEA.fa.seqs, SFA, trees, and every
# *_total , expansion , cafe_pvalue column.
#
# It writes to RDATA_OUT (a NEW file) and never overwrites RDATA_IN.
# ---------------------------------------------------------------------------

suppressMessages({
  library(DESeq2); library(dplyr); library(tibble); library(tidyr); library(stringr)
})

## ---- EDIT THESE -----------------------------------------------------------
RDATA_IN      <- "data/A2TEA_finished_8_29.RData"                # final tea output
RDATA_OUT     <- "data/A2TEA_finished_potatofix2.RData"       # new output (won't overwrite)
DEA_POT       <- "data/dea_Solanum_tuberosum"                    # corrected potato dea object (RDS)
SPECIES       <- "Solanum_tuberosum"
DEG_FDR       <- 0.1                                        # must match config.yaml["DEG_FDR"]
RECOMPUTE_TEA <- FALSE   # TRUE = also recompute norm_sum_*/oaes_hyper/tea_value (needs survcomp)
## ---------------------------------------------------------------------------

load(RDATA_IN)                 # brings HOG_DE.a2tea, HOG_level_list, hypotheses, etc.
dea_pot <- readRDS(DEA_POT)

## species that actually have DE data (the ones with *_sigDE columns)
species_list <- combined_AllSpeciesDEResultsDataFrames %>%
  filter(!is.na(baseMean)) %>% pull(species) %>% unique()

## ---- GATE: do potato DE gene IDs match the HOG table? ---------------------
pot_res   <- as.data.frame(results(dea_pot))
pot_genes <- rownames(pot_res)
overlap   <- length(intersect(pot_genes, HOG_file_long$gene))
overlap_table <- 
old_hog   <- HOG_DE.a2tea %>% filter(species == SPECIES, HOG != "singleton") %>% nrow()
cat(sprintf("Potato DE genes: %d | overlap with HOG map: %d | old potato-in-HOG rows: %d\n",
            length(pot_genes), overlap, old_hog))
if (overlap == 0)
  stop("No overlap between potato DE gene IDs and the HOG table. The peptide/HOG IDs ",
       "and the featureCounts gene_ids don't match -> harmonise IDs before patching.")
if (overlap < 0.5 * old_hog)
  warning("Overlap is much smaller than before - check ID formats before trusting output.")

## ---- 1. rebuild potato per-gene rows (exactly as the script does) ---------
pot_df <- pot_res
pot_df <- add_column(pot_df, gene = rownames(pot_df), .before = "baseMean")
pot_df <- add_column(pot_df, species = SPECIES, .before = "gene")
rownames(pot_df) <- NULL

combined_AllSpeciesDEResultsDataFrames <- combined_AllSpeciesDEResultsDataFrames %>%
  filter(species != SPECIES) %>%
  bind_rows(pot_df)

## ---- 2. rebuild HOG_DE.a2tea via the same join + significance rule --------
HOG_DE.a2tea <- full_join(combined_AllSpeciesDEResultsDataFrames, HOG_file_long,
                          by = "gene") %>%
  replace_na(list(HOG = "singleton"))
HOG_DE.a2tea$significant <- ifelse(!is.na(HOG_DE.a2tea$padj) & HOG_DE.a2tea$padj < DEG_FDR,
                                   "yes", "no")
HOG_DE.a2tea <- HOG_DE.a2tea %>% relocate(significant, .after = padj)

## ---- 3. new potato sig-DE count per HOG -----------------------------------
pot_sig <- HOG_DE.a2tea %>%
  filter(species == SPECIES, significant == "yes", HOG != "singleton") %>%
  count(HOG, name = "n")
cat(sprintf("Potato significant genes (padj < %.2g) in HOGs: %d across %d HOGs\n",
            DEG_FDR, sum(pot_sig$n), nrow(pot_sig)))

col <- paste0(SPECIES, "_sigDE")

## helper: set potato's sigDE column in a table to the new per-HOG counts (0-filled)
set_pot <- function(tbl) {
  if (!col %in% names(tbl)) return(tbl)
  tbl[[col]] <- pot_sig$n[match(tbl$HOG, pot_sig$HOG)]
  tbl[[col]][is.na(tbl[[col]])] <- 0
  tbl
}

## ---- 4. update each hypothesis table + the overview -----------------------
for (nm in names(HOG_level_list)) {
  tbl <- HOG_level_list[[nm]]

  if (nm == "all_species_overview") {
    HOG_level_list[[nm]] <- set_pot(tbl)          # only its potato sigDE column
    next
  }

  i <- as.integer(sub("hypothesis_", "", nm))
  hrow <- hypotheses[hypotheses$hypothesis == i, ]
  h_exp <- str_split(hrow$expanded_in,  ";")[[1]]
  h_cmp <- str_split(hrow$compared_to, ";")[[1]]
  h_species <- intersect(c(h_exp, h_cmp), species_list)

  if (!(SPECIES %in% h_species) || !(col %in% names(tbl))) next  # potato not in this hyp

  tbl <- set_pot(tbl)
  sig_cols <- paste0(h_species, "_sigDE")
  tbl$total_sigDE <- rowSums(tbl[, sig_cols, drop = FALSE], na.rm = TRUE)

  ## faithful recomputation of norm/oaes/tea (same formulas as the script)
  if (RECOMPUTE_TEA &&
      all(h_exp %in% species_list) && all(h_cmp %in% species_list) &&
      all(c("norm_sum_expanded_sigDE", "oaes_hyper", "cafe_pvalue") %in% names(tbl))) {
    if (!requireNamespace("survcomp", quietly = TRUE))
      stop("RECOMPUTE_TEA=TRUE needs the 'survcomp' package (bioconductor-survcomp).")
    ne <- rowSums(tbl[, paste0(h_exp, "_sigDE"), drop = FALSE], na.rm = TRUE) / length(h_exp)
    nc <- rowSums(tbl[, paste0(h_cmp, "_sigDE"), drop = FALSE], na.rm = TRUE) / length(h_cmp)
    tbl$norm_sum_expanded_sigDE <- ne
    tbl$norm_sum_compared_sigDE <- nc
    oa <- phyper(ne - 1, sum(ne), sum(nc), ne + nc, lower.tail = FALSE)
    oa[oa == 1] <- NA_real_
    tbl$oaes_hyper <- oa
    tbl$tea_value <- vapply(seq_len(nrow(tbl)), function(r) {
      vals <- c(tbl$oaes_hyper[r], tbl$cafe_pvalue[r])
      if (all(is.na(vals))) return(NA_real_)
      survcomp::combine.test(vals, na.rm = TRUE)
    }, numeric(1))
  } else if (!RECOMPUTE_TEA && "tea_value" %in% names(tbl)) {
    warning(nm, ": sigDE updated but norm_sum_*/oaes_hyper/tea_value left STALE ",
            "(RECOMPUTE_TEA=FALSE).")
  }

  HOG_level_list[[nm]] <- tbl
}

## ---- 5. save to a NEW file (originals preserved) --------------------------
save(hypotheses, HYPOTHESES.a2tea, A2TEA.fa.seqs, HOG_DE.a2tea, HOG_file_long,
     combined_AllSpeciesDEResultsDataFrames, HOG_level_list, SFA, all_speciesTree,
     file = RDATA_OUT, compress = "xz", compression_level = 9)

cat("\nWrote patched objects to:", RDATA_OUT, "\n")
cat("Potato rows in HOG_DE.a2tea now:",
    HOG_DE.a2tea %>% filter(species == SPECIES) %>% nrow(),
    "| in HOGs:", HOG_DE.a2tea %>% filter(species == SPECIES, HOG != "singleton") %>% nrow(), "\n")



####some sanity checks

new <- new.env(); load("data/A2TEA_finished_potatofix2.RData", envir = new)
old <- new.env(); load("data/A2TEA_finished_8_29.RData",      envir = old)   # <-- envir = old

# potato now present, significant count should be ~4,220 at padj < 0.1:
nrow(subset(new$HOG_DE.a2tea, species == "Solanum_tuberosum"))
nrow(subset(new$HOG_DE.a2tea, species == "Solanum_tuberosum" & significant == "yes"))

# confirm the OLD file really had the broken numbers (sanity check):
nrow(subset(old$HOG_DE.a2tea, species == "Solanum_tuberosum"))          # expect ~19,982

# expansion side untouched (expect TRUE):
identical(old$HOG_level_list$all_species_overview$Solanum_tuberosum_total,
          new$HOG_level_list$all_species_overview$Solanum_tuberosum_total)

# other species' DE untouched (expect TRUE):
identical(subset(old$HOG_DE.a2tea, species == "Capsicum_annuum"),
          subset(new$HOG_DE.a2tea, species == "Capsicum_annuum"))

col <- paste0(SPECIES, "_sigDE")

ov_new <- new$HOG_level_list[["all_species_overview"]]
ov_old <- old$HOG_level_list[["all_species_overview"]]
stopifnot(col %in% names(ov_new), col %in% names(ov_old))

## align old vs new potato sigDE by HOG
m <- merge(data.frame(HOG = ov_old$HOG, old = ov_old[[col]]),
           data.frame(HOG = ov_new$HOG, new = ov_new[[col]]),
           by = "HOG", all = TRUE)
m$old[is.na(m$old)] <- 0; m$new[is.na(m$new)] <- 0
m$delta   <- m$new - m$old
m$changed <- m$delta != 0
m$status  <- with(m, ifelse(old == 0 & new > 0, "gained",
                            ifelse(old > 0 & new == 0, "lost",
                                   ifelse(delta > 0, "increased",
                                          ifelse(delta < 0, "decreased", "unchanged")))))

cat("=== Global potato DE change (per HOG, count of significant potato genes) ===\n")
print(table(m$status))
cat(sprintf("\nHOGs changed in total: %d of %d\n",
            sum(m$changed), nrow(m)))
cat("  (gained/increased = recovered genes now significant;\n",
    "   decreased/lost can occur because adding ~13k genes re-estimates\n",
    "   size factors and independent filtering, shifting some padj values.)\n\n")

exp_hogs <- function(env, h) {
  tb <- env$HOG_level_list[[paste0("hypothesis_", h)]]
  if (is.null(tb) || !"expansion" %in% names(tb)) character(0) else tb$HOG[tb$expansion == "yes"]
}
if (is.null(HYPS)) HYPS <- sort(unique(as.integer(new$hypotheses$hypothesis)))

cat("=== Per-hypothesis: changed HOGs vs. the expanded set of interest ===\n")
report <- list()
for (h in HYPS) {
  eh  <- exp_hogs(new, h)
  sub <- m[m$HOG %in% eh, ]
  changed_bg_only <- sum(m$changed) - sum(sub$changed)   # changed but NOT expanded here
  report[[length(report)+1]] <- data.frame(
    hypothesis      = h,
    expanded_HOGs   = length(eh),
    expanded_changed= sum(sub$changed),
    exp_gained      = sum(sub$status == "gained"),
    exp_with_DE_old = sum(sub$old > 0),
    exp_with_DE_new = sum(sub$new > 0),
    changed_not_here= changed_bg_only)
}
print(do.call(rbind, report), row.names = FALSE)
cat("\nReading it:\n",
    " expanded_changed / expanded_HOGs = how much of THIS expanded set the fix moved\n",
    " exp_with_DE_old -> exp_with_DE_new = expanded HOGs carrying a potato DEG, before vs after\n",
    " changed_not_here = changed HOGs that fall OUTSIDE this expanded set (background)\n\n")



########
# quick sanity check that DE is calculated the same as in the workflow (apply same workflow as above to dea_solanum_pennellii)  

suppressMessages({ library(DESeq2); library(dplyr) })

DEA_PEN <- "data/dea_Solanum_pennellii"                 # the object you downloaded
RDATA   <- "data/A2TEA_finished_8_29.RData"        # ORIGINAL (pre-patch) object
SPECIES <- "Solanum_pennellii"
DEG_FDR <- 0.1
## --------------------

dea <- readRDS(DEA_PEN)

# our significance rule, exactly as run_potato_deseq2.R / A2TEA apply it
res <- as.data.frame(results(dea))
ours <- tibble(gene = rownames(res), padj = res$padj) %>%
  mutate(significant_ours = ifelse(!is.na(padj) & padj < DEG_FDR, "yes", "no"))

# significance already stored in the pipeline's object
e <- new.env(); load(RDATA, envir = e)
stored <- e$HOG_DE.a2tea %>%
  filter(species == SPECIES) %>%
  distinct(gene, significant)                      # HOG_DE has one sig call per gene

cmp <- inner_join(ours, stored, by = "gene")
cat(sprintf("genes compared: %d\n", nrow(cmp)))
cat(sprintf("matches: %d   mismatches: %d\n",
            sum(cmp$significant_ours == cmp$significant),
            sum(cmp$significant_ours != cmp$significant)))

# show any disagreements (expect none, or only NA-padj edge cases)
cmp %>% filter(significant_ours != significant) %>% arrange(padj) %>% head(20)

# no mismatches, i.e. significance calculated same as in original pipeline! Everything correct
