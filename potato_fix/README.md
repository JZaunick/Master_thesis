# Potato fix: recovering *S. tuberosum* genes lost in GFF3-to-GTF conversion

This folder documents how the differential-expression (DE) results for
*Solanum tuberosum* were corrected after the A2TEA.Workflow run, and contains
everything needed to reproduce the correction.

## The problem

Before read counting, A2TEA converts every annotation to GTF with
`gffread <gff3> -T` (rule `annotation_standardization`). The Phytozome
*S. tuberosum* DM v6.1 annotation (`Stuberosum_686_v6.1.gene.gff3`,
Pham et al. 2020) describes transcripts through `CDS` and UTR features and has
**no `exon` features**. gffread rebuilds `exon` lines only for transcripts that
carry at least one UTR feature; CDS-only transcripts are written without them.
featureCounts in the workflow counts with `-t exon`, so these genes received no
reads and were silently missing from the count matrix and all DE results.

| | genes |
|---|---|
| Genes in the DM v6.1 annotation | 32,917 |
| Genes with `exon` lines after gffread (quantifiable) | 19,982 |
| Genes lost (CDS-only transcripts) | 12,935 (39.3 %) |
| Genes with `exon` lines after AGAT | 32,917 |

The other four species are not affected: their annotations contain explicit
`exon` features, which gffread passes through unchanged.

The expansion analysis is not affected either, because it runs on the
longest-isoform proteomes, not on the GTF. Only the DE results of
*S. tuberosum* and everything derived from them in the final A2TEA object
(`HOG_DE.a2tea`, the `Solanum_tuberosum_sigDE` and `total_sigDE` columns of
`HOG_level_list`) needed correction.

## The fix, step by step

Run in this order. Steps 1–3 reuse the existing STAR alignments; no
realignment is needed, since the annotation is only used at the counting step.

| Step | Script | Input | Output |
|---|---|---|---|
| 0–2 | `01_preprocess_annotation_and_counts.sh` | DM v6.1 GFF3, STAR BAMs | AGAT GTF, `counts_potato.txt` |
| 3 | `run_potato_deseq2.R` | `counts_potato.txt`, `R/gen_fasta_samples.csv` | `dea_Solanum_tuberosum`, results table |
| 4 | `patch_potato_into_a2tea.R` | original A2TEA RData, `dea_Solanum_tuberosum` | `A2TEA_finished_potatofix.RData` |
| 5 (optional) | `04_verify_de_procedure.R` | `dea_Solanum_pennellii`, original RData | match report |

### Steps 0–2: diagnosis, AGAT conversion, read counting

`01_preprocess_annotation_and_counts.sh` does three things and prints the
expected counts after each:

0. **Diagnosis.** Lists the feature types in the GFF3 (no `exon`), runs the
   standard gffread conversion and counts genes with transcript lines (32,917)
   versus genes with exon lines (19,982).
1. **Conversion with AGAT.** `agat_convert_sp_gff2gtf.pl` rebuilds the full
   gene → mRNA → exon hierarchy, giving exon lines for all 32,917 genes.
   AGAT labels transcript lines `mRNA` rather than `transcript`; featureCounts
   only reads the exon lines, so this makes no difference.
2. **Read counting** with the workflow's parameters
   (`-O -M -t exon -g gene_id`) plus `-p` for the paired-end libraries
   (`fq2` set in `config/samples.tsv`), run in the workflow's own subread
   environment (`workflow/envs/subread.yaml`, which installs subread 2.0.0).
   In subread versions before 2.0.2, `-p` alone counts fragments; from 2.0.2
   on, `--countReadPairs` is also required. The original fix was run with
   subread 2.0.0 and `-p` only; the script checks the installed version and
   adds `--countReadPairs` only when a newer subread needs it.

### Step 3: DESeq2

`run_potato_deseq2.R` reproduces `workflow/scripts/gen_fasta_deseq2.R`
without Snakemake: same design (`~ condition`), same `DESeq()` call with the
same dispersion fallback. It was run in the workflow's own DESeq2 environment
(`workflow/envs/deseq2_tximport.yaml`; DESeq2 1.34.0, R 4.1.3). Sample-to-condition assignment is matched by name.
Run it from the A2TEA working directory (it reads `R/gen_fasta_samples.csv`
by relative path) with `WORKDIR` pointing at the output directory of step 2.

Result: 32,917 genes tested (29,473 with non-zero counts), 4,220 DEGs at
padj < 0.1 (1,924 up, 2,296 down), compared with 3,541 (1,645 up, 1,896 down)
before the fix.

### Step 4: patching the final A2TEA object

`patch_potato_into_a2tea.R` replaces the *S. tuberosum* rows of the per-gene
DE tables with the corrected results, rebuilds `HOG_DE.a2tea` with the same
join and significance rule as `workflow/scripts/final_tea_computation.R`
(`padj < DEG_FDR`), and recomputes `Solanum_tuberosum_sigDE`, `total_sigDE`
and the derived columns in `HOG_level_list`. Expansion columns (`*_total`,
`expansion`, `cafe_pvalue`) are not touched. It writes a new RData file and
never overwrites the original.

Before writing, it checks that the potato gene IDs match the HOG table
(30,584 of 32,917 genes are assigned to HOGs; the rest are singletons).
`RECOMPUTE_TEA = TRUE` also recomputes `oaes_hyper` and `tea_value` and
requires the `survcomp` package.

### Step 5 (optional): verifying the procedure

`04_verify_de_procedure.R` recomputes significance for *S. pennellii*, which
was not modified, with the same rule and compares it with the calls stored
by the workflow. Zero mismatches confirm that the manual rerun follows the
workflow's procedure exactly.

## Expected effect on the final object

| | before | after |
|---|---|---|
| *S. tuberosum* genes in `HOG_DE.a2tea` | 19,982 | 32,917 |
| ... assigned to HOGs | 19,449 | 30,584 |
| *S. tuberosum* DEGs (padj < 0.1) | 3,541 | 4,220 |

At the HOG level, 580 HOGs gained a significant potato gene and 80 gained
additional ones; 151 lost their only significant potato gene and 2 lost one
of several (26,689 unchanged). The losses are borderline genes that crossed
the 0.1 threshold when the multiple-testing correction was recomputed on the
larger gene set (816 potato genes now have 0.1 < padj < 0.15).

## Software

| Tool | Version | Environment | Used for |
|---|---|---|---|
| gffread | 0.12.1 | own conda env, same pin as `workflow/envs/gffread.yaml` | diagnosis |
| AGAT | 1.7.0 | own conda env (`conda create -n agat -c bioconda -c conda-forge agat=1.7.0`) | GFF3 → GTF conversion |
| subread / featureCounts | 2.0.0 | workflow env `workflow/envs/subread.yaml` | read counting (`-p` only; `--countReadPairs` does not exist in 2.0.0) |
| STAR | 2.7.8a | workflow env `workflow/envs/star.yaml` | alignments (original run) |
| R / DESeq2 | 4.1.3 / 1.34.0 | workflow env `workflow/envs/deseq2_tximport.yaml` | DE analysis |
| A2TEA.Workflow | commit 4a00657 (groupschoof/A2TEA.Workflow) | – | original workflow run |

## Not included

The STAR BAM files, the count matrices and the A2TEA RData objects are too
large for this repository. The BAMs and the original RData are produced by
the A2TEA.Workflow; the remaining files are regenerated by the scripts above.

## References

- Pham GM et al. (2020) Construction of a chromosome-scale long-read reference
  genome assembly for potato. *GigaScience* 9: giaa100.
- Pertea G, Pertea M (2020) GFF Utilities: GffRead and GffCompare.
  *F1000Research* 9: 304.
- Dainat J. AGAT: Another Gff Analysis Toolkit. Zenodo,
  doi:10.5281/zenodo.3552717.
- Liao Y, Smyth GK, Shi W (2014) featureCounts. *Bioinformatics* 30: 923–930.
- Love MI, Huber W, Anders S (2014) DESeq2. *Genome Biology* 15: 550.
