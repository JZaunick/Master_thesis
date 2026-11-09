# Cold adaptation in Solanaceae: analysis code

Code accompanying the master's thesis

> Zaunick, J. (2026). *Integrating differential expression with gene family evolution to
> prioritize candidate genes for cold adaptation in Solanaceae.* Master's thesis,
> University of Bonn.

Gene families (hierarchical orthogroups, HOGs) of three cold-tolerant species
(*Solanum pennellii*, *S. tuberosum*, *Lycium barbarum*) and two cold-sensitive species
(*S. lycopersicum*, *Capsicum annuum*) were inferred and tested for expansion with the
A2TEA workflow, combined with cold-stress RNA-seq data, and the resulting candidate
families were validated at the genome level.

## Repository structure

```
config/                 A2TEA configuration: config.yaml, hypotheses.tsv, species.tsv,
                        samples.tsv (with SRA accessions), primer_adapter/adapters.fa
isoform_filtering/      longest-isoform proteomes used as A2TEA input
potato_fix/             corrected S. tuberosum annotation, DE rerun, patching the A2TEA output
downstream_analysis/    statistics, expansion tests, candidate set, GO enrichment (R)
validation/             validation protocol and microsynteny analysis of candidate families
data/                   GO enrichment results and derived tables, flank files
```

## Order of the analysis

| Step | Folder / script | Thesis section |
|---|---|---|
| 1 | `isoform_filtering/` – one protein per gene (see its README) | 2.2.1 |
| 2 | A2TEA run with `config/` | 2.2.1 |
| 3 | `potato_fix/` – scripts 01–04 in order (see its README) → `A2TEA_finished_potatofix2.RData` | 2.2.2 |
| 4 | `downstream_analysis/basicStats.R`, `DE_gene_overview.R`, `newick_trees.R`, `upsetR_plots.R` | 3.1–3.3, 3.5, 3.6 |
| 5 | `downstream_analysis/Expansion_analysis_complete.R`, `expansion_de_piecharts.R` | 3.4, 3.5 |
| 6 | `downstream_analysis/candidate_table.R`, `GO_enrichment_topGO_withGoji.R` (enrichment, then tables and figures in `data/GOEnrichment_potatofix/`), `regulation_direction_test.R` (Table4_regulation_direction.csv) | 3.6 |
| 7 | `validation/` – protocol (README), `01_extract_flanks.sh`, `02_microsynteny_analysis.R` | 2.3.4, 3.7 |

The R scripts in steps 4–7 read `data/A2TEA_finished_potatofix2.RData` and are run from the
repository root.

## Data

Raw RNA-seq reads are available from the NCBI Sequence Read Archive under the accessions
listed in `config/samples.tsv`. Genome assemblies and annotations: RefSeq
GCF_001406875.1 (*S. pennellii*), GCF_019175385.1 (*L. barbarum*), Ensembl Plants release 60
(*S. lycopersicum* SL3.0, *C. annuum* ASM51225v2), Phytozome (*S. tuberosum* DM v6.1).

The A2TEA output object (`A2TEA_finished_potatofix2.RData`) is not included because of its
size; it is available from the author on request.

Included in `data/`:

- `GOEnrichment_potatofix/` – GO enrichment results, tables and figures
  (output of `downstream_analysis/GO_enrichment_topGO_withGoji.R`) and
  `Table4_regulation_direction.csv` (output of `downstream_analysis/regulation_direction_test.R`)
- `flank_*.tsv` – flanking genes of each examined copy (output of
  `validation/01_extract_flanks.sh`, input of `validation/02_microsynteny_analysis.R`)

## Software

A2TEA, OrthoFinder, STAR, featureCounts, DESeq2, gffread 0.12.7, R with tidyverse, topGO,
UpSetR, ComplexUpset, ape/ggtree; BLAST+ 2.17.0, samtools 1.24 and IGV for validation.

## Licence

GNU General Public License v3.0. `isoform_filtering/custom_longest_isoforms2.py` is a
modified version of `primary_transcript.py` from
[OrthoFinder](https://github.com/davidemms/OrthoFinder) (Emms & Kelly 2019), which is
distributed under the same licence.
