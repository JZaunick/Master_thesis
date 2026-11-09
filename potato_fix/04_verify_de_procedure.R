#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Verification: is the manual DESeq2 rerun identical to the A2TEA procedure?
# Recomputes significance for a species that was NOT modified (S. pennellii)
# with the same rule used for the potato rerun (results(dea), padj < DEG_FDR)
# and compares it gene by gene with the calls stored in the original A2TEA
# object. Zero mismatches = identical procedure.
# ---------------------------------------------------------------------------
suppressMessages({ library(DESeq2); library(dplyr) })

## ---- edit these ----
DEA_PEN <- "dea_Solanum_pennellii"            # dea object produced by the workflow
RDATA   <- "data/A2TEA_finished_original.RData" # ORIGINAL (pre-fix) A2TEA object
SPECIES <- "Solanum_pennellii"
DEG_FDR <- 0.1                                 # config.yaml: DEG_FDR
## --------------------

res  <- as.data.frame(results(readRDS(DEA_PEN)))
ours <- tibble(gene = rownames(res), padj = res$padj) %>%
  mutate(significant_ours = ifelse(!is.na(padj) & padj < DEG_FDR, "yes", "no"))

e <- new.env(); load(RDATA, envir = e)
stored <- e$HOG_DE.a2tea %>%
  dplyr::filter(species == SPECIES) %>%
  dplyr::distinct(gene, significant)

cmp <- inner_join(ours, stored, by = "gene")
cat(sprintf("genes compared: %d | mismatches: %d\n",
            nrow(cmp), sum(cmp$significant_ours != cmp$significant)))
