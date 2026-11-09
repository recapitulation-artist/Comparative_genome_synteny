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

Main retains its category mapping with saturation multiplied by 0.70.
Other entries rank categories by count descending, then name descending in byte
order, using 42 colors (seven hues, six tones each); overflow becomes gray `Other`.
Automatic uses gene names, falling back to products; unknown annotations also
become `Other`. Targets remain white with black dots. Genus circles are 2.2 mm.

Each output directory contains `synteny.pdf`, `synteny.png`, `caption.txt`,
`gene_colors.tsv` and plotting audit tables. Curated and automatic also contain
`legend.pdf` with full category names, extracted sequences, normalized annotations,
BED intervals, BLAST results, alignments, trees, species tables and provenance.
PDFs are vector graphics with embedded fonts; PNG previews use 600 dpi.
