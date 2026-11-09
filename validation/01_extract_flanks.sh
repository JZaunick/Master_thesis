#!/usr/bin/env bash
# ============================================================
# Microsynteny validation, step 1: flanking genes of each family copy
#
# For every copy of a candidate family (and for the position where a copy
# is suspected to be missing), extract the 15 protein-coding genes on either
# side from the genome annotation. One TSV per locus is written to data/:
#
#   offset  chr  start  strand  gene_id
#
# offset 0 is the target gene itself, negative offsets lie upstream
# (lower coordinates), positive offsets downstream. Only genes on the same
# chromosome/scaffold as the target are reported, so a locus near a
# scaffold end yields fewer than 31 rows.
#
# The TSV files are the input of 02_microsynteny_analysis.R.
# ============================================================
set -euo pipefail

# ---- annotations (adjust to your layout) --------------------------------
FS=~/experiments/exp_solcoldres_woTob_maniso/A2TEA.Workflow/FS

SLYC_GFF=$FS/Solanum_lycopersicum.SL3.0.60.gff3.gz          # Ensembl Plants, IDs "gene:Solyc..N"
SPEN_GFF=$FS/GCF_001406875.1_SPENNV200_genomic.gff.gz       # RefSeq,         IDs "gene-LOC..."
LBAR_GFF=$FS/GCF_019175385.1_ASM1917538v2_genomic.gff.gz    # RefSeq,         IDs "gene-LOC..."
CANN_GFF=$FS/Capsicum_annuum.ASM51225v2.60.gff3.gz          # Ensembl Plants
STUB_GFF=$FS/Stuberosum_686_v6.1.gene.gff3.gz              # Phytozome,      IDs "Soltu.DM..."

# ---- flank(): $1 = GFF(.gz), $2 = exact "ID=...;" string, $3 = output -----
# Protein-coding filter: GFFs without a biotype attribute (Spud DB) pass all
# genes; RefSeq/Ensembl GFFs keep only biotype=protein_coding.
# The trailing ";" in $2 prevents partial matches (Solyc10g079595.1 vs .10),
# so the exact ID including any version suffix must be given. Look it up
# with:  zcat $SLYC_GFF | grep -m1 "Solyc10g079595"
flank () {
  zcat "$1" \
  | awk -F'\t' '$3=="gene" && ($9 !~ /biotype=/ || $9 ~ /biotype=protein_coding/)' \
  | sort -t$'\t' -k1,1 -k4,4n \
  | awk -F'\t' -v t="$2" -v n=15 '
      { c[NR]=$1; s[NR]=$4; st[NR]=$7
        id=$9; sub(/^ID=/,"",id); sub(/;.*/,"",id); g[NR]=id
        if (index($9,t)) h=NR }
      END {
        if (!h) { print "TARGET NOT FOUND: " t > "/dev/stderr"; exit 1 }
        for (i=h-n; i<=h+n; i++)
          if (i>0 && c[i]==c[h]) print i-h"\t"c[i]"\t"s[i]"\t"st[i]"\t"g[i]
      }' > "$3"
  echo "$(wc -l < "$3") genes -> $3"
}

mkdir -p data

# ---- loci ---------------------------------------------------------------
# File names encode <family>_<lineage/copy>_<species>; the R script derives
# the locus label from the file name (flank_<label>.tsv).

# The target gene of every existing file is its offset-0 row:
#   for f in data/flank_*.tsv; do printf "%s\t" "$f"; awk -F'\t' '$1==0{print $5}' "$f"; done

# N0.HOG0006681 -- type-A response regulator (ARR9), Section 3.7.2
# lineage A
flank $STUB_GFF "ID=Soltu.DM.10G027810;"     data/flank_ARR_A_Stub.tsv
flank $SPEN_GFF "ID=<ID>;"                   data/flank_ARR_A_Spen.tsv
flank $LBAR_GFF "ID=<ID>;"                   data/flank_ARR_A_Lbar.tsv
flank $SLYC_GFF "ID=gene:Solyc10g079595.1;"  data/flank_ARR_A_Slyc.tsv
flank $CANN_GFF "ID=<ID>;"                   data/flank_ARR_A_Cann.tsv
# lineage B (tomato copy concealed in the fused model Solyc10g079710)
flank $STUB_GFF "ID=Soltu.DM.10G027680;"     data/flank_ARR_B_Stub.tsv
flank $SPEN_GFF "ID=gene-LOC107032291;"      data/flank_ARR_B_Spen.tsv
flank $LBAR_GFF "ID=<ID>;"                   data/flank_ARR_B_Lbar.tsv

# <family/HOG>, locus L04
flank $STUB_GFF "ID=<ID>;"                   data/flank_L04_Stub.tsv
flank $SPEN_GFF "ID=<ID>;"                   data/flank_L04_Spen.tsv
flank $LBAR_GFF "ID=<ID>;"                   data/flank_L04_Lbar.tsv
flank $SLYC_GFF "ID=<ID>;"                   data/flank_L04_Slyc.tsv

# N0.HOG0000291 -- SOFL, locus L12
flank $STUB_GFF "ID=<ID>;"                   data/flank_L12_Stub.tsv
flank $SPEN_GFF "ID=<ID>;"                   data/flank_L12_Spen.tsv
flank $LBAR_GFF "ID=<ID>;"                   data/flank_L12_Lbar.tsv
flank $CANN_GFF "ID=<ID>;"                   data/flank_L12_Cann.tsv

# N0.HOG0003220 -- XTH, potato copy on chromosome 7
flank $STUB_GFF "ID=<ID>;"                   data/flank_XTH_Stub7.tsv
