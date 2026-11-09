# ============================================================
# Microsynteny validation, step 2: shared flanking HOGs and
# orthologous loci (R)
#
# Input : data/A2TEA_finished_potatofix2.RData   (A2TEA workflow objects)
#         data/flank_*.tsv                        (from 01_extract_flanks.sh)
# Output: printed tables; optionally results/*.tsv (see end of script)
#
# Objects used from the RData file:
#   HOG_DE.a2tea        per-gene table: species, gene, HOG, baseMean,
#                       log2FoldChange, padj
#   A2TEA.fa.seqs       protein sequences of all genes (named list)
#   HYPOTHESES.a2tea    expanded_OGs per hypothesis, incl. MSA member names
#
# Logic:
#   1. Every flanking gene is mapped to its HOG via HOG_DE.a2tea.
#   2. Two loci that share many flanking HOGs (out of at most 30) are
#      orthologous/same lineage; loci of different lineages share few.
#   3. The flanking HOGs of a copy that is present in one species are looked
#      up in the species that lacks the copy. The tomato genes carrying those
#      HOGs delimit the syntenic interval that is then searched with tblastn
#      (validation protocol, step 3).
#   4. A protein-length screen over the examined families flags members whose
#      length deviates strongly from the family median (fused or truncated
#      gene models).
# ============================================================

library(dplyr); library(readr); library(purrr); library(tidyr); library(tibble)

load("data/A2TEA_finished_potatofix2.RData")

# ---- helpers ------------------------------------------------------------
# GFF gene IDs carry source-specific prefixes ("gene-" RefSeq, "gene:"
# Ensembl); the workflow gene table mostly does not.
norm   <- function(x) sub("^gene[:-]", "", x)
hogmap <- HOG_DE.a2tea %>% distinct(gene, HOG) %>%
  mutate(key = norm(gene)) %>% distinct(key, HOG)

# Read all flank_<prefix>_*.tsv files of one family and attach HOGs.
read_flanks <- function(prefix, dir = "data") {
  files <- list.files(dir, pattern = paste0("^flank_", prefix, ".*\\.tsv$"),
                      full.names = TRUE)
  stopifnot(length(files) > 0)
  map_dfr(files, ~ read_tsv(.x, col_names = c("offset", "chr", "start", "strand", "id"),
                            col_types = "icicc") %>%
            mutate(locus = gsub(paste0("^flank_", prefix, "|\\.tsv$"), "", basename(.x)))) %>%
    mutate(key = norm(id)) %>%
    left_join(hogmap, by = "key")
}

# Pairwise number of shared flanking HOGs (target gene itself excluded).
shared_matrix <- function(fl) {
  m  <- fl %>% filter(offset != 0, grepl("^N0\\.HOG", HOG)) %>% distinct(locus, HOG)
  sh <- inner_join(m, m, by = "HOG", relationship = "many-to-many") %>%
    filter(locus.x < locus.y)
  list(pairs  = sh,
       matrix = sh %>% count(locus.x, locus.y, name = "shared") %>%
         pivot_wider(names_from = locus.y, values_from = shared, values_fill = 0))
}

# Genes of `target_species` that belong to the HOGs flanking `locus`
# (|offset| <= window). Their chromosomal order in the target species
# delimits the syntenic interval.
syntenic_genes <- function(fl, locus, target_species, window = 5) {
  fl %>% filter(locus == !!locus, abs(offset) <= window, offset != 0,
                grepl("^N0\\.HOG", HOG)) %>%
    select(offset, id, HOG) %>%
    left_join(HOG_DE.a2tea %>% filter(species %in% target_species) %>%
                distinct(species, gene, HOG), by = "HOG") %>%
    arrange(offset)
}

# ======================================================================
# A. SOB FIVE-LIKE (SOFL) family (N0.HOG0006736), chromosome 12 locus L12 -- loci L12_Stub, L12_Spen, L12_Lbar, L12_Cann
# ======================================================================
fl_l12 <- read_flanks("L12_") %>% mutate(locus = paste0("L12_", locus))

# how many of the 30 flanking genes could be assigned to a HOG
fl_l12 %>% group_by(locus) %>% summarise(genes = n(), with_HOG = sum(!is.na(HOG)))

sh_l12 <- shared_matrix(fl_l12)
sh_l12$matrix

# which flanking HOGs the L. barbarum and S. tuberosum loci share, side by side
pair <- c("L12_Lbar", "L12_Stub")
fl_l12 %>%
  filter(locus %in% pair,
         HOG %in% (sh_l12$pairs %>% filter(locus.x %in% pair, locus.y %in% pair) %>% pull(HOG))) %>%
  select(locus, offset, strand, id, HOG) %>% arrange(HOG, locus)

# tomato genes belonging to the HOGs that flank the potato/S. pennellii copies
l12_hogs <- fl_l12 %>% filter(locus %in% c("L12_Stub", "L12_Spen"), offset != 0,
                              grepl("^N0\\.HOG", HOG)) %>% distinct(HOG)
HOG_DE.a2tea %>% filter(species == "Solanum_lycopersicum", HOG %in% l12_hogs$HOG) %>%
  distinct(gene, HOG) %>% arrange(gene)

# the five genes directly flanking the potato copy and their tomato orthologues
syntenic_genes(fl_l12, "L12_Stub", "Solanum_lycopersicum")

# ======================================================================
# B. ARR9 family (N0.HOG0006681) -- lineage A: Stub, Spen, Lbar, Slyc, Cann;
#    lineage B: Stub, Spen, Lbar
# ======================================================================
fl_arr <- read_flanks("ARR_")

shared_matrix(fl_arr)$matrix          # A vs B lineages share few flanking HOGs

# tomato and pepper orthologues of the genes flanking the potato lineage-B copy
syntenic_genes(fl_arr, "B_Stub", c("Solanum_lycopersicum", "Capsicum_annuum")) %>%
  print(n = Inf)

# The interval between the tomato orthologues of the immediate neighbours
# (Solyc10g079690 / Solyc10g079710) was then searched with tblastn
# (validation protocol, step 3). The missing lineage-B copy turned out to
# be concealed in the fused gene model Solyc10g079710.

# properties of the fused tomato model and of the C. annuum member (T459_ = pepper)
HOG_DE.a2tea %>% filter(gene %in% c("T459_26538", "Solyc10g079710.2")) %>%
  select(species, gene, HOG, baseMean, log2FoldChange, padj)

# ======================================================================
# C. SOFL family (N0.HOG0006736), chromosome 4 locus L04 --
#    loci L04_Stub, L04_Spen, L04_Lbar, L04_Slyc (same steps as for L12)
# ======================================================================
fl_l04 <- read_flanks("L04_") %>% mutate(locus = paste0("L04_", locus))

fl_l04 %>% group_by(locus) %>% summarise(genes = n(), with_HOG = sum(!is.na(HOG)))
shared_matrix(fl_l04)$matrix
syntenic_genes(fl_l04, "L04_Stub", "Solanum_lycopersicum")

# ======================================================================
# D. XTH family (N0.HOG0003220) -- potato copy on chromosome 7
#    tomato orthologues of the genes flanking the potato copy
# ======================================================================
fl_xth <- read_flanks("XTH_")
syntenic_genes(fl_xth, "Stub7", "Solanum_lycopersicum") %>% print(n = Inf)

# ======================================================================
# E. Protein-length screen of the examined families
# ======================================================================
len <- sapply(A2TEA.fa.seqs, function(s) nchar(paste(unlist(s), collapse = "")))

len[c("T459_26538", "T459_26557", "Solyc10g079710.2")]   # ARR9: fused tomato model vs. C. annuum members

# all members of one family, sorted by length
ids <- names(HYPOTHESES.a2tea$hypothesis_4@expanded_OGs[["N0.HOG0005673"]]@add_OG_analysis[[1]]@msa)
sort(len[ids])

# members deviating > 1.5x or < 0.6x from the family median, over the six
# examined candidate families
candidate_hogs <- c("N0.HOG0001280",   # PUB22/23
                    "N0.HOG0006681",   # ARR9
                    "N0.HOG0006736",   # SOFL
                    "N0.HOG0005673",   # PME
                    "N0.HOG0000291",   # UGT
                    "N0.HOG0003220")   # XTH

eogs <- HYPOTHESES.a2tea$hypothesis_4@expanded_OGs[candidate_hogs]

lenflags <- imap_dfr(eogs, function(eog, hog) {
  ids <- names(eog@add_OG_analysis[[1]]@msa)
  l   <- len[ids]
  tibble(HOG = hog, gene = ids, length = l, ratio = round(l / median(l), 2))
}) %>% filter(ratio > 1.5 | ratio < 0.6) %>% arrange(HOG, ratio)

lenflags
lenflags %>% distinct(HOG)

# ---- optional: write results ----------------------------------------------
# dir.create("results", showWarnings = FALSE)
# write_tsv(sh_l12$matrix,                 "results/shared_flanking_HOGs_L12.tsv")
# write_tsv(shared_matrix(fl_arr)$matrix,  "results/shared_flanking_HOGs_ARR9.tsv")
# write_tsv(lenflags,                      "results/length_flags.tsv")
