#!/usr/bin/env Rscript
# Manual, Snakemake-free reproduction of A2TEA's gen_fasta_deseq2.R for potato.
# Reads from the (read-only is fine) home dir, writes everything into WORKDIR so it
# works even when your home quota is full. Run inside the pipeline's DESeq2 conda env:
#
#   WORKDIR=/scratch/$USER/potato_fix Rscript run_potato_deseq2.R
#
# or just edit the paths below.

suppressMessages(library(DESeq2))

## ---- settings ---------------------------------------------------------------
SPECIES   <- "Solanum_tuberosum"                       # exact string in gen_fasta_samples.csv
WORKDIR   <- Sys.getenv("WORKDIR", ".")                # writable dir for all outputs
COUNTFILE <- file.path(WORKDIR, "counts_potato.txt")   # featureCounts output from previous step
SAMPLES   <- "R/gen_fasta_samples.csv"                 # A2TEA sample table (tab-delimited); copy down if local
REF_LEVEL <- "control"                                 # reference level of 'condition', if present
OUTRDS    <- file.path(WORKDIR, "dea_Solanum_tuberosum")            # the dea object (A2TEA's format/name)
OUTCSV    <- file.path(WORKDIR, "results_Solanum_tuberosum.csv")    # readable results, for you
## -----------------------------------------------------------------------------

## counts: raw featureCounts output has 1 comment line, then the header
fc <- read.delim(COUNTFILE, skip = 1, header = TRUE,
                 check.names = FALSE, row.names = "Geneid")
counts <- as.matrix(fc[, -(1:5)])                      # drop Chr,Start,End,Strand,Length
colnames(counts) <- sub("_Aligned.sortedByCoord.out.bam$", "",
                        basename(colnames(counts)))     # recover sample_unit keys

## sample table, subset to this species
st <- read.delim(SAMPLES)
st <- st[st$species == SPECIES, ]
if (nrow(st) == 0)
  stop("No rows for species '", SPECIES, "'. Present: ",
       paste(unique(read.delim(SAMPLES)$species), collapse = ", "))
st$su <- paste0(st$sample, "_", st$unit)
rownames(st) <- st$su

## align condition to count columns BY NAME (order-independent)
missing <- setdiff(colnames(counts), rownames(st))
if (length(missing))
  stop("Count columns absent from sample table: ", paste(missing, collapse = ", "))
coldata <- data.frame(condition = st[colnames(counts), "condition"],
                      row.names = colnames(counts))
coldata$condition <- factor(coldata$condition)
if (REF_LEVEL %in% levels(coldata$condition))
  coldata$condition <- relevel(coldata$condition, ref = REF_LEVEL)

cat("Samples:\n"); print(coldata)
cat("\nGenes in count matrix:", nrow(counts), "\n\n")

## DESeq2: same design + dispersion fallback as A2TEA
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
tryCatch({
  dea <<- DESeq(dds)
}, error = function(b) {
  dea  <- estimateSizeFactors(dds)
  dea  <- estimateDispersionsGeneEst(dea)
  dispersions(dea) <- mcols(dea)$dispGeneEst
  dea <<- nbinomWaldTest(dea)
})

dir.create(dirname(OUTRDS), recursive = TRUE, showWarnings = FALSE)
saveRDS(dea, file = OUTRDS)
res <- results(dea)
write.csv(as.data.frame(res), OUTCSV)

cat("Saved dea object ->", OUTRDS, "\n")
cat("Saved results     ->", OUTCSV, "\n\n")
cat("Genes tested:", nrow(dea), "  |  padj < 0.1:", sum(res$padj < 0.1, na.rm = TRUE), "\n\n")
summary(res)
