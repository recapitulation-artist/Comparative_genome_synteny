# Cd Uptake Synteny

Three reproducible entries generate synteny figures. Paths are relative to the
repository root; original source files are preserved.

| Entry | Inputs | Output |
| --- | --- | --- |
| `plot-main` | Supplied sequence index, curated annotations, BLAST hits, reference tree and species metadata | `results/main/` |
| `plot-curated` | Raw FASTA/GFF, curated gene annotations and genome metadata | `results/curated/` |
| `plot-auto` | Raw FASTA/GFF, genome metadata and target gene name | `results/auto/` |

## Setup

Linux x86_64, Python >=3.9, Make and Slurm are required. Locked versions:
BEDTools 2.31.0, BLAST+ 2.17.0+, MAFFT 7.526, stable IQ-TREE 3.1.4 and R 4.6.0.

```sh
# Restore R if a matching environment is unavailable.
conda create --name cd_uptake_synteny --file environment/conda_linux_64.lock
conda activate cd_uptake_synteny
make setup
```

Submission recognizes the existing local R environment. Executable overrides:
`RSCRIPT`, `BEDTOOLS`, `BLASTN`, `MAKEBLASTDB`, `MAFFT`, `IQTREE`.
Bundled data require no downloads.

## Run an entry

```sh
# Main: plot supplied files directly; no raw genomes or reconstruction required.
make submit TARGET=plot-main

# Curated: extract, BLAST, align, infer a tree and plot using manual annotations.
make submit TARGET=plot-curated

# Automatic: discover every target copy and select its coding-gene neighborhood.
make submit TARGET=plot-auto TARGET_GENE=Nramp NEIGHBORS=10

# Build and verify all entries, or verify one entry.
make submit TARGET=verify
make submit TARGET=verify-auto TARGET_GENE=Nramp NEIGHBORS=10

# Replay using this entry's recorded initial random seed.
make submit TARGET=verify-curated MODE=replay
make submit TARGET=verify-auto MODE=replay
```

`verify-main` and `verify-curated` also build and verify individual entries.
Submission uses node2, one task and 96 CPUs; logs are in `results/slurm/`.
BLAST, MAFFT and IQ-TREE use 96 threads. MAFFT and IQ-TREE retain default analysis
settings. Actual seeds, models, commands, versions and hashes are recorded.
Matching tree caches are reused; incompatible caches stop with an error.
Use a fresh repository copy for a different target or neighborhood configuration.
Replay outputs are kept under the entry's `phylogeny_replay/` directory.

## Inputs and parameters

Main reads `data/sequences.csv`, `data/genes.tsv`,
`data/reference_blast_hits.csv`, `data/phylogeny.nwk` and `data/genomes.tsv`.

Curated and automatic accept `GENOME_DIR` and `GENOME_METADATA`; curated also
accepts `GENE_ANNOTATIONS`:

```sh
make submit TARGET=plot-auto TARGET_GENE=Nramp NEIGHBORS=10 \
  GENOME_DIR=data/genomes GENOME_METADATA=data/genomes.tsv
make submit TARGET=plot-curated GENE_ANNOTATIONS=data/genes.tsv
```

Each assembly directory contains `genome.fna.gz` and `annotation.gff.gz`.
Metadata requires `assembly_accession` and `organism_name`; curated additionally
requires unique `seq_id` values matching its annotation table. Curated derives
bounds from the manual annotations, preserving the historical one-base coordinate
convention and reproducing the reference FASTA without reading it.

Automatic matches `mntH` or products containing `Nramp` for the default target;
other names match gene names or product tokens case-insensitively. Every target
copy is retained with up to `NEIGHBORS` coding genes on each side. RNA and
pseudogenes are excluded. Contig ends truncate the neighborhood without
compensation. At least four regions are required for tree inference.

`prepare-curated` / `prepare-auto`, `blast-curated` / `blast-auto` and
`phylogeny-curated` / `phylogeny-auto` run individual stages. Make targets run
directly inside the allocation; preparation can also run locally. `make run`
builds all figures. Earlier plotting target names remain aliases.

## Figures

All entries accept independent gene and homology-link palettes. Defaults are
`GENE_PALETTE=main-muted` and `LINK_PALETTE=blue`.

| `GENE_PALETTE` | 42-color palette |
| --- | --- |
| `main-muted` | Historical main colors at 0.70 saturation, extended with 16 darker tones |
| `automatic` | Existing seven-hue, six-tone automatic palette |
| `nature` | [NPG-inspired](https://nanx.me/ggsci/reference/pal_npg.html), extended with lighter tones |
| `cell` | Project-defined soft Cell-inspired colors |
| `science` | [AAAS-inspired](https://nanx.me/ggsci/reference/pal_aaas.html), extended with lighter tones |

`main-muted` retains the historical main category mapping. Other main palettes
rank its 26 ordinary categories by count descending, then name descending in
byte order. Curated and automatic rank the same way and use the first 42 colors;
overflow and unknown annotations become gray `Other`. Automatic uses gene names,
falling back to products. Targets remain white with black dots; main `Unknown`
remains black. Genus circles are 2.2 mm. Journal labels describe color inspiration.

`LINK_PALETTE` accepts `blue`, `green`, `purple`, `orange` (three pale tones), or
`original` (historical grays). `LINK_COLORS` overrides that preset with three
distinct comma-separated `#RRGGBB` colors, ordered by BLAST identity:
**>=90%, 80–<90%, 70–<80%**.

```sh
make submit TARGET=plot-main GENE_PALETTE=main-muted LINK_PALETTE=blue
make submit TARGET=plot-curated GENE_PALETTE=nature LINK_PALETTE=green
make submit TARGET=plot-auto GENE_PALETTE=automatic LINK_PALETTE=original
make submit TARGET=plot-auto GENE_PALETTE=science \
  LINK_COLORS="#AECDE3,#D0E1EE,#EBF3F8"
```

To redraw and verify existing results without extraction, BLAST or tree inference,
run R directly (activate the R environment first):

```sh
for entry in main curated auto; do
  Rscript --vanilla scripts/plot_synteny.R --entry "$entry" \
    --gene-palette main-muted --link-palette blue
  Rscript --vanilla scripts/verify_results.R --entry "$entry" \
    --gene-palette main-muted --link-palette blue
done
```

Both R commands also accept `--link-colors "#AECDE3,#D0E1EE,#EBF3F8"`,
`--target-gene NAME`, and `--phylogeny-dir results/curated/phylogeny_replay`
(or the corresponding automatic directory) when replay inputs are selected.
Use the same palette options for plotting and verification.

Each output directory contains `synteny.pdf`, `synteny.png`, `caption.txt`,
`gene_colors.tsv`, `link_colors.tsv`, selected palettes in `figure_style.tsv`,
and plotting audit tables. Curated and automatic also contain
`legend.pdf` with full category names, extracted sequences, normalized annotations,
BED intervals, BLAST results, alignments, trees, species tables and provenance.
PDFs are vector graphics with embedded fonts; PNG previews use 600 dpi.
