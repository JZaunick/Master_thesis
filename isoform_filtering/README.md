# Input preparation: longest-isoform proteomes

A2TEA was run with `auto_isoform_filtering: "NO"` (config.yaml), so one representative
protein per gene was selected beforehand. For every species, the longest protein isoform
per gene was kept, using `custom_longest_isoforms2.py`.

`custom_longest_isoforms2.py` is a modified version of `primary_transcript.py` from
OrthoFinder (Emms & Kelly 2019; https://github.com/davidemms/OrthoFinder), distributed
under the GNU GPL v3. Modifications: gene-name extraction by attribute key (`by_key`) or by
user-defined keywords (`custom`), optional prefix/suffix for gene names, and optional removal
of an isoform suffix pattern (regex) before grouping (`--trim_suffix_pattern`). Among isoforms
of equal length, the first one in the file is kept.

The gene identifiers written to the FASTA headers must match the gene IDs in each species'
annotation (config.yaml: "the names between pep. fasta and either genes or transcripts in
the species specific GTF has to match"), because they connect orthology to expression.

Software: gffread 0.12.7, Python 3 (standard library only).

### *Solanum pennellii* (RefSeq GCF_001406875.1)

Protein sequences derived from the annotation with gffread; gene ID from the `gene=`
attribute, prefixed with `gene-` to match the gene IDs of the RefSeq GFF.

```bash
gffread GCF_001406875.1_SPENNV200_genomic.gff -g GCF_001406875.1_SPENNV200_genomic.fna \
  -FSy SPENNV200_proteins_with_all_infos.fasta

python3 custom_longest_isoforms2.py SPENNV200_proteins_with_all_infos.fasta SPENNV200_longestiso.pep.fa \
  --gene_name_function by_key --key gene --prefix "gene-"
```

### *Solanum tuberosum* (Phytozome, DM v6.1)

Filtered from the published protein FASTA; gene ID from the `locus=` field, with the
suffix `.v6.1` added to match the gene IDs of the Phytozome GFF3.

```bash
zcat Stuberosum_686_v6.1.protein.fa.gz > Stuberosum_686_v6.1.protein.fa

python3 custom_longest_isoforms2.py Stuberosum_686_v6.1.protein.fa Stuberosum_longestiso.pep.fa \
  --gene_name_function custom --keyword locus --suffix .v6.1
# (--keyword is accepted by argparse as an abbreviation of --keywords)
```

Cross-check 1 – the same selection starting from the annotation. Here the isoforms share no
common ID (`.1`, `.2`, … per isoform), so the isoform suffix is removed before grouping:

```bash
gffread Stuberosum_686_v6.1.gene.gff3 -g Stuberosum_686_v6.1.fa -FSy Stub_proteins_with_all_infos.fasta

python3 custom_longest_isoforms2.py Stub_proteins_with_all_infos.fasta stub_longestisofromgff.pep.fa \
  --gene_name_function by_key --key Name --trim_suffix_pattern '\.[0-9]+' --suffix .v6.1
```

Both routes give 32,917 proteins.

Cross-check 2 – the Phytozome annotation flags each gene's longest transcript with
`longest=1`, which gffread `-F` carries into the headers:

```bash
grep -c 'longest=1' Stub_proteins_with_all_infos.fasta    # 32917
```

The count equals the number of selected proteins.

### *Lycium barbarum* (RefSeq GCF_019175385.1)

gffread cannot parse `?` in the strand column, so these were set to `.` first. Gene ID from
the `gene=` attribute.

```bash
sed 's/\t?\t/\t.\t/g' GCF_019175385.1_ASM1917538v2_genomic.gff > Lbarb_fixed.gff

gffread Lbarb_fixed.gff -g GCF_019175385.1_ASM1917538v2_genomic.fna -Sy Lbarb_proteins_with_all_infos.fasta

python3 custom_longest_isoforms2.py Lbarb_proteins_with_all_infos.fasta Lbar_longestiso.pep.fasta \
  --gene_name_function by_key --key gene
```

The output was used in A2TEA as `Lbarbarum_longestiso.pep.fa` (species.tsv).

### *Solanum lycopersicum* and *Capsicum annuum* (Ensembl Plants, release 60)

Filtered directly from the Ensembl protein FASTA files; gene ID from the `gene:` field of
the Ensembl headers (e.g. `Solyc01g005000.3`).

```bash
zcat Solanum_lycopersicum.SL3.0.pep.all.fa.gz > Solanum_lycopersicum.SL3.0.pep.all.fa
python3 custom_longest_isoforms2.py Solanum_lycopersicum.SL3.0.pep.all.fa Solanum_lycopersicum_longestiso.pep.fa \
  --gene_name_function ensembl

zcat Capsicum_annuum.ASM51225v2.pep.all.fa.gz > Capsicum_annuum.ASM51225v2.pep.all.fa
python3 custom_longest_isoforms2.py Capsicum_annuum.ASM51225v2.pep.all.fa Capsicum_annuum_longestiso.pep.fa \
  --gene_name_function ensembl
```

## Verifying that the deposited script reproduces the proteomes used

Rerun the commands above into a scratch directory and compare with the files A2TEA used:

```bash
md5sum FS/*_longestiso.pep.fa* scratch/*_longestiso.pep.fa*
```

## Checking the headers against the annotations

```bash
grep -m3 ">" *_longestiso.pep.fa*
grep -c ">" *_longestiso.pep.fa*
```
