#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Potato fix, steps 0-2: diagnose the annotation, re-convert it with AGAT,
# and re-count reads with featureCounts against the existing STAR alignments.
#
# Run from the A2TEA.Workflow working directory (where star/ lives).
# Requires: gffread 0.12.1 (diagnosis only), AGAT, subread/featureCounts.
# Adjust the variables below; no paths outside these need changing.
# ---------------------------------------------------------------------------
set -euo pipefail

GFF3="Stuberosum_686_v6.1.gene.gff3"          # Phytozome DM v6.1 annotation (input)
GTF_GFFREAD="Solanum_tuberosum_gffread.gtf"   # standard A2TEA conversion (diagnosis only)
GTF_AGAT="Solanum_tuberosum_agat.gtf"         # corrected conversion (used downstream)
BAM_DIR="star/Solanum_tuberosum"              # existing STAR alignments (not regenerated)
OUTDIR="potato_fix_out"                       # writable output directory
THREADS=4

mkdir -p "$OUTDIR"

# --- Step 0: diagnosis -------------------------------------------------------
# The DM v6.1 GFF3 has CDS and UTR features but no 'exon' features.
echo "== Feature types in the source annotation =="
grep -v '^#' "$GFF3" | cut -f3 | sort | uniq -c

# Standard conversion as performed by the A2TEA workflow (rule annotation_standardization).
gffread "$GFF3" -T -o "$OUTDIR/$GTF_GFFREAD"

echo "== Genes with transcript lines after gffread (expected 32917) =="
grep -P '\ttranscript\t' "$OUTDIR/$GTF_GFFREAD" | grep -oP 'gene_id "\K[^"]+' | sort -u | wc -l
echo "== Genes with exon lines after gffread (expected 19982) =="
grep -P '\texon\t' "$OUTDIR/$GTF_GFFREAD" | grep -oP 'gene_id "\K[^"]+' | sort -u | wc -l
# featureCounts counts with -t exon, so only the second number is quantifiable.
# The missing 12,935 genes are those whose transcripts carry CDS but no UTR features.

# --- Step 1: re-convert with AGAT ------------------------------------------
agat_convert_sp_gff2gtf.pl --gff "$GFF3" -o "$OUTDIR/$GTF_AGAT"

echo "== Genes with exon lines after AGAT (expected 32917) =="
grep -P '\texon\t' "$OUTDIR/$GTF_AGAT" | grep -oP 'gene_id "\K[^"]+' | sort -u | wc -l
# Note: AGAT labels transcript lines 'mRNA' instead of 'transcript'; this does not
# affect featureCounts, which reads only the exon lines.

# --- Step 2: re-count reads (A2TEA flags + paired-end) -----------------------
# Same parameters as the workflow rule (-O -M -t exon -g gene_id), run in the
# workflow's own subread environment (workflow/envs/subread.yaml, subread ==2.0).
# The reads are paired-end (fq2 present in config/samples.tsv), hence -p.
# In subread < 2.0.2, -p alone counts fragments; from 2.0.2 on, --countReadPairs
# is needed as well. The flag is added only if the installed version needs it.
PAIR_FLAGS="-p"
if featureCounts -v 2>&1 | grep -qE 'v2\.0\.([2-9]|[1-9][0-9])|v2\.[1-9]'; then
  PAIR_FLAGS="-p --countReadPairs"
fi
echo "featureCounts version: $(featureCounts -v 2>&1 | grep -oE 'v[0-9.]+')  -> using: $PAIR_FLAGS"

featureCounts -T "$THREADS" $PAIR_FLAGS -O -M -t exon -g gene_id \
  -a "$OUTDIR/$GTF_AGAT" \
  -o "$OUTDIR/counts_potato.txt" \
  "$BAM_DIR"/*_Aligned.sortedByCoord.out.bam

echo "== Genes in count matrix (expected 32917) =="
tail -n +3 "$OUTDIR/counts_potato.txt" | wc -l

echo "Done. Next: Rscript run_potato_deseq2.R (with WORKDIR=$OUTDIR)"
