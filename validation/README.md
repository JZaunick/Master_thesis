# Validation of candidate gene families

This protocol documents the checks used to validate the candidate HOGs examined in
detail in the thesis (Methods 2.3.4, Table `tab:validation`). Each check answers one
question:

| Step | Question | Tools |
|---|---|---|
| 0 | Which genes belong to the HOG, and are their lengths consistent? | `grep`, `samtools faidx` |
| 1 | Does an overlong gene model contain several family members? | `blastp` |
| 2 | Where is the orthologous locus in a species that seems to lack a copy? | `awk`, R |
| 3 | Is there an unannotated or concealed copy at that locus? | `samtools faidx`, `tblastn` |
| 4 | Do RNA-seq reads support the annotated gene model? | STAR `SJ.out.tab`, IGV |

The worked example throughout is the type-A response regulator family
**N0.HOG0006681** (ARR9). Its lineage-B copy is present in *S. tuberosum*
(`Soltu.DM.10G027680`) and *S. pennellii* (`LOC107032291`) but missing from the
*S. lycopersicum* proteome. It turned out to be concealed in the fused gene model
`Solyc10g079710` (chromosome 10).

Software: BLAST+ 2.17.0, samtools 1.24, R 4.x with tidyverse, IGV.

Scripts in this folder:

- `01_extract_flanks.sh` – step 2a for all loci examined (flanking genes of each copy).
- `02_microsynteny_analysis.R` – steps 2b/2c (shared flanking HOGs, orthologous loci in the
  species lacking a copy) and the protein-length screen of step 0.

The remaining steps are run interactively with the commands given below.

---

## Setup

Adjust the paths to your own layout. Gene and protein identifiers must match the
FASTA headers of the longest-isoform proteomes; check them with `grep ">"` before
running the commands below.

```bash
FS=~/experiments/exp_solcoldres_woTob_maniso/A2TEA.Workflow/FS     # proteomes and annotations
N0=path/to/Phylogenetic_Hierarchical_Orthogroups/N0.tsv            # OrthoFinder HOG table

STUB_PEP=$FS/Stuberosum_longestiso.pep.fa
SPEN_PEP=$FS/SPENNV200_longestiso.pep.fa
SLYC_PEP=$FS/Solanum_lycopersicum_longestiso.pep.fa

SLYC_GFF=$FS/Solanum_lycopersicum.SL3.0.60.gff3.gz
SPEN_GFF=$FS/GCF_001406875.1_SPENNV200_genomic.gff.gz
LBAR_GFF=$FS/GCF_019175385.1_ASM1917538v2_genomic.gff.gz
CANN_GFF=$FS/Capsicum_annuum.ASM51225v2.60.gff3.gz
STUB_GFF=$FS/Stuberosum_686_v6.1.gene.gff3.gz

# index the proteomes once (creates *.fai with sequence lengths)
for f in $STUB_PEP $SPEN_PEP $SLYC_PEP; do samtools faidx $f; done
```

---

## Step 0 – Members and protein lengths of a HOG

List the members of a HOG. In `N0.tsv`, columns 1–3 hold the HOG, the OG and the gene
tree parent clade, and the remaining columns hold the genes per species:

```bash
grep -w "N0.HOG0006681" $N0 | cut -f4- | tr '\t,' '\n\n' | sed 's/^ *//' | grep -v '^$' > hog_members.txt
```

Protein lengths come from column 2 of the `.fai` index:

```bash
cat $STUB_PEP.fai $SPEN_PEP.fai $SLYC_PEP.fai | grep -F -f hog_members.txt | cut -f1,2
```

Alternatively, compute the lengths in R from the workflow object:

```r
len <- sapply(A2TEA.fa.seqs, function(s) nchar(paste(unlist(s), collapse = "")))
```

A model about twice or three times the family length suggests a fusion (step 1).
A model about half the family length suggests a split or truncated model.

Check whether a gene was assigned to any HOG. A count of 0 means it is unassigned:

```bash
grep -c "T459_26538" $N0
```

---

## Step 1 – Testing an overlong model for fused genes

Align one family member of typical length against the overlong protein. Several
full-length hits to consecutive, non-overlapping segments of the subject indicate
fused genes.

```bash
samtools faidx $STUB_PEP Soltu.DM.10G027680.1 > query.faa      # typical member
samtools faidx $SLYC_PEP Solyc10g079710.1.1   > subject.faa    # overlong model

blastp -query query.faa -subject subject.faa -evalue 1e-10 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue"
```

Read the `sstart`/`send` columns. In the PME family, for example, the *S. tuberosum*
protein hit `Solyc07g064190` at positions 1–540, 538–1062 and 1063–1605, which shows
three fused copies.

---

## Step 2 – Locating the orthologous locus by microsynteny

### 2a. Flanking genes of each copy

For each family member, extract the 15 protein-coding genes on either side from its
genome annotation. The function writes one row per gene: offset from the target,
chromosome, start, strand and gene ID.

```bash
flank () {  # $1 = GFF, $2 = exact "ID=...;" string of the target gene, $3 = output file
zcat "$1" | awk -F'\t' '$3=="gene" && ($9 !~ /biotype=/ || $9 ~ /biotype=protein_coding/)' \
 | sort -t$'\t' -k1,1 -k4,4n \
 | awk -F'\t' -v t="$2" -v n=15 '{c[NR]=$1; s[NR]=$4; st[NR]=$7; id=$9; sub(/^ID=/,"",id); sub(/;.*/,"",id); g[NR]=id; if(index($9,t)) h=NR}
   END{ if(!h){print "TARGET NOT FOUND: " t > "/dev/stderr"; exit}
        for(i=h-n;i<=h+n;i++) if(i>0 && c[i]==c[h]) print i-h"\t"c[i]"\t"s[i]"\t"st[i]"\t"g[i] }' > "$3"
}

mkdir -p data
flank $STUB_GFF "ID=Soltu.DM.10G027680;"   data/flank_ARR_B_Stub.tsv
flank $STUB_GFF "ID=Soltu.DM.10G027810;"   data/flank_ARR_A_Stub.tsv
flank $SPEN_GFF "ID=gene-LOC107032291;"    data/flank_ARR_B_Spen.tsv
flank $SLYC_GFF "ID=gene:Solyc10g079595.1;" data/flank_ARR_A_Slyc.tsv
# ... one call per copy and species
```

The ID format differs between sources:
- RefSeq: prefix `gene-`.
- Ensembl: prefix `gene:`, and the tomato IDs carry a version suffix (`.1`, `.2`, `.3`).
- Potato annotation: no prefix.

Including the trailing `;` prevents partial matches, so the version must be written
out. Look the ID up first, for example with
`zcat $SLYC_GFF | grep -m1 "Solyc10g079595"`.

### 2b. Shared flanking HOGs between loci (R)

```r
library(tidyverse)
load("data/A2TEA_finished_potatofix2.RData")

norm   <- function(x) sub("^gene[:-]", "", x)
hogmap <- HOG_DE.a2tea %>% distinct(gene, HOG) %>%
  mutate(key = norm(gene)) %>% distinct(key, HOG)

files <- list.files("data", pattern = "^flank_.*\\.tsv$", full.names = TRUE)
fl <- map_dfr(files, ~ read_tsv(.x, col_names = c("offset", "chr", "start", "strand", "gene"),
                                col_types = "icicc") %>%
                mutate(locus = sub("\\.tsv$", "", basename(.x)))) %>%
  filter(offset != 0) %>%
  mutate(key = norm(gene)) %>%
  left_join(hogmap, by = "key")

sets   <- split(fl$HOG, fl$locus) %>% map(~ unique(na.omit(.x)))
shared <- outer(names(sets), names(sets),
                Vectorize(function(a, b) length(intersect(sets[[a]], sets[[b]]))))
dimnames(shared) <- list(names(sets), names(sets))
shared   # number of shared flanking HOGs (out of at most 30) per pair of loci
```

Copies of the same lineage share many flanking HOGs. Copies of different lineages
share few.

### 2c. Finding the syntenic interval in the species lacking a copy (R)

Map all genes of the target species to their HOGs. Then look up where the flanking
HOGs of a reference copy lie in that species.

```r
gff_genes <- function(gff) {
  read_tsv(gff, comment = "#", col_names = FALSE, col_types = cols(.default = "c")) %>%
    filter(X3 == "gene") %>%
    transmute(chr = X1, start = as.integer(X4), end = as.integer(X5), strand = X7,
              key = norm(str_match(X9, "ID=([^;]+)")[, 2]))
}

slyc <- gff_genes("FS/Solanum_lycopersicum.SL3.0.60.gff3.gz") %>%   # FS = folder from the setup
  inner_join(hogmap, by = "key")

fl %>% filter(locus == "flank_ARR_B_Stub") %>%
  select(offset, HOG) %>%
  inner_join(slyc, by = "HOG") %>%
  arrange(offset)
```

A run of target genes in conserved order on one chromosome identifies the locus. The
interval between the genes corresponding to the immediate neighbours of the reference
copy is searched in step 3. In the example, these were `Solyc10g079690` and
`Solyc10g079710` on chromosome 10.

---

## Step 3 – Searching the interval with tblastn

`samtools faidx` needs an uncompressed or bgzip-compressed genome:

```bash
zcat Solanum_lycopersicum.SL3.0.dna.toplevel.fa.gz > SL3.0.fa
samtools faidx SL3.0.fa

# interval around the syntenic position (here generously ±30 kb)
samtools faidx SL3.0.fa 10:61300000-61370000 > interval.fa

# queries: the family members of the species that carry the copy
samtools faidx $STUB_PEP Soltu.DM.10G027680.1 >  query_B.faa
samtools faidx $SPEN_PEP LOC107032291         >> query_B.faa

tblastn -query query_B.faa -subject interval.fa -evalue 1e-5 \
  -outfmt "6 qseqid pident length qstart qend sstart send sframe evalue sseq"
```

How to read the output:

- **Coverage:** do the hits together span the query from position 1 to its end? A
  complete copy is covered in several segments, one per exon.
- **Identity:** orthologues between the *Solanum* species are typically above 90 %.
- **Stop codons:** an asterisk (`*`) in `sseq` marks an internal stop codon. No stop
  codons means an intact reading frame; stop codons indicate a pseudogene.
- **Genomic position:** interval start + `sstart` − 1.

Compare the hit coordinates with the annotated gene models:

```bash
zcat $SLYC_GFF | awk -F'\t' -v c=10 -v s=61332360 -v e=61334115 \
  '$1==c && ($3=="gene" || $3=="mRNA") && $4<=e && $5>=s {print $3, $4, $5, $7, $9}'
```

Possible outcomes:

- No overlapping model: an unannotated copy.
- The hit lies inside a model that extends well beyond it: a fused model. Test that
  model with step 1.
- Hits only within the exons of the known copy: no additional copy. This is a valid
  negative result.

---

## Step 4 – RNA-seq support for the gene model

### 4a. Splice junctions from STAR (scriptable)

`SJ.out.tab` columns: chromosome, intron start, intron end, strand (1 = +, 2 = −),
motif, annotated (1/0), uniquely mapping reads, multi-mapping reads, maximum overhang.

```bash
for f in path/to/star/*SJ.out.tab; do
  echo "== $f"
  awk -v c=10 -v s=61325000 -v e=61345000 '$1==c && $2>=s && $3<=e' "$f"
done
```

Check two things:

- **Junctions inside the family copy:** these show that the copy is transcribed and
  spliced on its own.
- **The annotated intron joining the copy to the next gene:** no reads, or only one
  or two, means the fused model is not supported.

### 4b. Visual inspection in IGV

1. Load the genome (`SL3.0.fa`), the annotation (GFF3) and the sorted, indexed BAM
   files of all samples of the species (`samtools index sample.bam` if no `.bai`
   exists).
2. Navigate to the interval (e.g. `10:61,325,000-61,345,000`).
3. Compare the coverage track with the exons of the annotated model. Then open a
   Sashimi plot (right-click on a BAM track) to see junction read counts.
4. Save a screenshot (`File > Save PNG Image`) for the record.

Genes with little expression in leaves cannot be confirmed this way. In that case,
rely on steps 1 and 3.

---

## Helpers

Chromosome names for RefSeq assemblies, where sequences are named by accession:

```bash
zcat $SPEN_GFF | awk -F'\t' '$3=="region" && match($9, /chromosome=[^;]+/) \
  {print $1, substr($9, RSTART+11, RLENGTH-11)}'
```

All genes in a window, for example to describe a tandem array:

```bash
zcat $SLYC_GFF | awk -F'\t' '$3=="gene" && $1=="1" && $4>1560000 && $5<1600000 \
  {print $4, $5, $7, $9}' | sort -n
```
