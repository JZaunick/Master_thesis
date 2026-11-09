#analysis flanking genes


library(dplyr); library(readr); library(purrr); library(tidyr)
load("data/A2TEA_finished_potatofix2.RData")


dir   <- "data"
files <- list.files(dir, pattern = "^flank_.*\\.tsv$", full.names = TRUE)
norm   <- function(x) sub("^gene[:-]", "", x)   # GFF IDs carry "gene-"/"gene:" prefixes; HOG table mostly doesn't
hogmap <- HOG_DE.a2tea %>% distinct(gene, HOG) %>% mutate(key = norm(gene)) %>% distinct(key, HOG)

fl <- map_dfr(files, ~
                read_tsv(.x, col_names = c("offset","chr","start","strand","id"),
                         col_types = "icicc") %>%
                mutate(locus = gsub("^flank_|\\.tsv$", "", basename(.x)))) %>%
  mutate(key = norm(id)) %>%
  left_join(hogmap, by = "key")

fl %>% group_by(locus) %>% summarise(genes = n(), with_HOG = sum(!is.na(HOG)))

m  <- fl %>% filter(offset != 0, !is.na(HOG)) %>% distinct(locus, HOG)
sh <- inner_join(m, m, by = "HOG", relationship = "many-to-many") %>% filter(locus.x < locus.y)

sh %>% count(locus.x, locus.y, name = "shared") %>%
  pivot_wider(names_from = locus.y, values_from = shared, values_fill = 0)

pair <- c("L12_Lbar", "L12_Stub")
fl %>% filter(locus %in% pair, HOG %in% (sh %>% filter(locus.x %in% pair, locus.y %in% pair) %>% pull(HOG))) %>%
  select(locus, offset, strand, id, HOG) %>% arrange(HOG, locus)


l12_hogs <- fl %>% filter(locus %in% c("L12_Stub","L12_Spen"), offset != 0, grepl("^N0\\.HOG", HOG)) %>% distinct(HOG)
HOG_DE.a2tea %>% filter(species == "Solanum_lycopersicum", HOG %in% l12_hogs$HOG) %>% distinct(gene, HOG) %>% arrange(gene)

#Map potato genes directly flanking potato SOFL copy to their tomato orthologues:
fl %>% filter(locus == "L12_Stub", abs(offset) <= 5, offset != 0, grepl("^N0\\.HOG", HOG)) %>%
  select(offset, id, HOG) %>%
  left_join(HOG_DE.a2tea %>% filter(species == "Solanum_lycopersicum") %>% distinct(gene, HOG), by = "HOG") %>%
  arrange(offset)

# 




###################################


ids <- names(HYPOTHESES.a2tea$hypothesis_4@expanded_OGs[["N0.HOG0005673"]]@add_OG_analysis[[1]]@msa)
sort(sapply(A2TEA.fa.seqs[ids], function(s) nchar(paste(unlist(s), collapse = ""))))


len <- sapply(A2TEA.fa.seqs, function(s) nchar(paste(unlist(s), collapse = "")))
eogs <- HYPOTHESES.a2tea$hypothesis_4@expanded_OGs[candidate_hogs]

lenflags <- purrr::imap_dfr(eogs, function(eog, hog) {
  ids <- names(eog@add_OG_analysis[[1]]@msa)
  l <- len[ids]
  tibble::tibble(HOG = hog, gene = ids, length = l, ratio = round(l / median(l), 2))
}) %>% dplyr::filter(ratio > 1.5 | ratio < 0.6) %>% dplyr::arrange(HOG, ratio)

lenflags %>%  select(HOG) %>% unique()



##############################
#check for ARR9 HOG

files <- list.files("data", pattern = "^flank_ARR_.*\\.tsv$", full.names = TRUE)
fl <- map_dfr(files, ~ read_tsv(.x, col_names = c("offset","chr","start","strand","id"), col_types = "icicc") %>%
                mutate(locus = gsub("^flank_ARR_|\\.tsv$", "", basename(.x)))) %>%
  mutate(key = norm(id)) %>% left_join(hogmap, by = "key")

m  <- fl %>% filter(offset != 0, grepl("^N0\\.HOG", HOG)) %>% distinct(locus, HOG)
inner_join(m, m, by = "HOG", relationship = "many-to-many") %>% filter(locus.x < locus.y) %>%
  count(locus.x, locus.y, name = "shared") %>% pivot_wider(names_from = locus.y, values_from = shared, values_fill = 0)


fl %>% filter(locus == "B_Stub", abs(offset) <= 5, offset != 0, grepl("^N0\\.HOG", HOG)) %>%
  select(offset, id, HOG) %>%
  left_join(HOG_DE.a2tea %>% filter(species %in% c("Solanum_lycopersicum","Capsicum_annuum")) %>%
              distinct(species, gene, HOG), by = "HOG") %>% arrange(offset) %>% print(n = Inf)


HOG_DE.a2tea %>% filter(gene %in% c("T459_26538", "Solyc10g079710.2")) %>%
  select(species, gene, HOG, baseMean, log2FoldChange, padj)

len[c("T459_26538", "T459_26557", "Solyc10g079710.2")]   # len from the length screen; recreate if needed:
len <- sapply(A2TEA.fa.seqs, function(s) nchar(paste(unlist(s), collapse = "")))


HOG_DE.a2tea %>% filter(gene %in% c("T459_26538", "Solyc10g079710.2")) %>%
  select(species, gene, HOG, baseMean, log2FoldChange, padj)
