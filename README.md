# Enterobacteriaceae Cd Uptake Comparative Genome Synteny

Compare cadmium-uptake-associated gene neighborhoods in Enterobacteriaceae.
The workflow extracts genomic regions with BEDTools, compares nucleotide
homology with BLAST, aligns sequences with MAFFT, infers phylogenies with stable
IQ-TREE, and renders species-annotated synteny figures. Three entries reproduce
the supplied figure, reconstruct with curated annotations, or discover targets
from original genome annotations without manual curation.

All analysis inputs are bundled. Complete original FASTA/GFF files for 29
assemblies are losslessly compressed (about 53 MiB); their decompressed bytes
match the original files. Unrelated reports and downloads are omitted. No step
reads neighboring data directories. This repository starts from `7bce487`.

## Installation and deployment

Supported: **Linux x86_64 / Windows WSL2 x86_64**, glibc >=2.28. Use a Linux
filesystem on WSL2. Software and package requirements are listed separately in
[DEPENDENCIES.md](DEPENDENCIES.md); exact locks and fonts are in `environment/`.

After obtaining the repository, enter its root and run:

```bash
# Install Miniforge/Conda first if conda is unavailable.
conda create --yes --prefix .tools/r --file environment/conda_linux_64.lock
conda activate "$PWD/.tools/r"
make setup
make check
make verify-main
```

The host needs Python >=3.9 and GNU coreutils. The R lock provides Make,
compilers, curl and graphics libraries. `make setup` restores remaining tools and
the five required R packages; `setup-r` installs packages into `.r_library/`
after the locked R runtime has been created. Software installation requires
internet access; analysis does not download genome data or install system-wide.

Matching existing tools can be selected explicitly or through PATH. Resolution
is explicit override, repository `.tools`, then PATH. Personal environment names
and home-directory paths are never searched. After moving to another machine or
path, restore Conda prefixes and compiled R libraries there instead of copying
them. Transfer all tracked files, including data, locks and fonts.

## Three entries and required files

| Entry | Required inputs | Output |
| --- | --- | --- |
| `plot-main` | Supplied index, manual annotations, historical BLAST, reference tree, species metadata | `results/main/` |
| `plot-curated` | Original FASTA/GFF, manual annotations, species metadata | `results/curated/` |
| `plot-auto` | Original FASTA/GFF, species metadata, target name | `results/auto/` |

### Main: reproduce the supplied figure

Prepare or verify these files:

| File | Key checks |
| --- | --- |
| `data/sequences.csv` | Unique `seq_id`, positive `length`; preserve the supplied sequence index |
| `data/genes.tsv` | `seq_id`, `start`, `end`, `strand`, `note`, `organisms`; target rows use `note=Uptake-Gene` |
| `data/reference_blast_hits.csv` | Both sequence IDs exist in the index; coordinates refer to those regions |
| `data/phylogeny.nwk` | Valid Newick; exactly one tip per indexed sequence |
| `data/genomes.tsv` | Sequence IDs map to assembly accessions and species names |

```bash
make check-data
make plot-main
make verify-main
```

Main reads these fixed filenames. It needs R and basic system commands, and does
not read raw genomes or run BLAST, MAFFT or IQ-TREE. Genome-path overrides below
apply to the reconstructed entries.

### Curated: reconstruct from genomes and manual annotations

Required: `data/genomes/`, `data/genomes.tsv`, and `data/genes.tsv`. The last is
a byte-identical standardized copy of `29gene2.txt`. This entry does not read the
historical extracted FASTA, supplied sequence index, historical BLAST or reference
tree. Region bounds are derived from the manual annotations.

Metadata requires `seq_id`, `assembly_accession`, and `organism_name`; `seq_id`
must be unique, while assemblies and species may repeat.
Annotations require `seq_id`, `start`, `end`, `strand`, `type`, `feat_id`,
`geom_id`, `note`, and `organisms`. Each selected sequence needs a target row and
feature IDs matching the original GFF. Keep the supplied annotation convention:
its historical one-base start offset is preserved for exact FASTA reproduction.

```bash
make plot-curated
make verify-curated
# Optional exact historical comparison.
cmp data/reference_sequences.fasta results/curated/extracted_sequences.fasta
```

### Automatic: discover targets without manual annotations

Only original FASTA/GFF, species metadata and a target query are required.
`data/genes.tsv`, historical fragments, historical BLAST and the reference tree
are never read. Standardized input layout:

```text
data/genomes/
  GCF_003470965.1/
    genome.fna.gz
    annotation.gff.gz
  ...
data/genomes.tsv
```

Uncompressed `genome.fna` and `annotation.gff` are also accepted. FASTA contig
names must match GFF sequence IDs. GFF must provide parent genes, CDS features,
their `ID`/`Parent` relationships, and gene or product descriptions. The automatic
metadata table minimally needs these tab-separated columns:

```tsv
assembly_accession	organism_name
GCF_003470965.1	Escherichia coli
```

Index every assembly directory; extra metadata columns are preserved. Curated
metadata additionally needs `seq_id` matching the manual annotation table.

```bash
make plot-auto TARGET_GENE=Nramp NEIGHBORS=10
make verify-auto TARGET_GENE=Nramp NEIGHBORS=10
```

`Nramp` matches `mntH` gene names or products containing `Nramp`. Other targets
match gene names or whole product tokens case-insensitively. This searches
existing GFF annotations; it does not perform de novo annotation or protein
homology discovery. Inspect `target_matches.tsv` for matches and exclusions.

Every eligible target copy is retained. `NEIGHBORS=10` means up to ten coding
parent genes on each genomic side, not 10 kb. CDS fragments count through their
parent gene; RNAs and pseudogenes are excluded. Contig ends truncate the window
without compensation. At least four regions are required for tree inference.

## Local and Slurm execution

The default total CPU budget is **56**. Slurm is optional:

```bash
make run THREADS=56                  # All three figures, locally
make verify THREADS=56               # Build and validate all three
make plot-auto THREADS=56 TARGET_GENE=Nramp NEIGHBORS=10
```

BLAST, MAFFT and IQ-TREE each receive the requested count. Preparation uses up to
`min(THREADS, assembly count)` workers with independent single-thread BEDTools
processes. Results retain stable input order. CPU-heavy stages run sequentially,
even under `make -j`, to avoid multiplying the CPU budget. Plotting, short checks
and required data dependencies retain necessary serial work. Actual CPU
utilization depends on the tool and data size.

For a smaller machine, explicitly choose a suitable budget, for example:

```bash
make plot-auto THREADS=28 TARGET_GENE=Nramp
```

Slurm submission requests one node, one task and `THREADS` CPUs. Node and
partition default to the scheduler's choice:

```bash
make submit TARGET=plot-auto THREADS=56 TARGET_GENE=Nramp
make submit TARGET=verify-curated THREADS=56 \
  SLURM_PARTITION=compute SLURM_NODE=nodeA SLURM_TIME=04:00:00 SLURM_MEMORY=32G
# Inspect arguments without submitting.
python3 scripts/submit_workflow.py --dry-run --target plot-auto --threads 56
make -s submit TARGET=plot-auto SUBMIT_DRY_RUN=1
```

Replace the example partition/node with cluster-specific values. The login node
need not have 56 CPUs; allocation checks run on the compute node. Direct Make
targets also work inside an existing suitable allocation. Logs go to
`results/slurm/`.
Compute nodes must see the same repository path, installed environments and
explicit executable paths as the submission host, normally through shared storage.

## Optional parameters

Each table gives defaults, applicability and a separate example. Analysis input
paths must resolve inside this repository; copy custom data into `data/` first.
Tool executables may live outside it. Run commands from the root and quote paths
and colors. Use a fresh working copy when changing datasets or target settings.

### Resources and inputs

| Parameter | Default and meaning | Example |
| --- | --- | --- |
| `THREADS` | `56`; positive CPU budget, reconstructed entries and Slurm | `make plot-curated THREADS=56` |
| `TARGET_GENE` | `Nramp`; automatic query and reconstructed-figure target label | `make plot-auto TARGET_GENE=mntH` |
| `NEIGHBORS` | `10`; nonnegative coding genes per side, automatic only | `make plot-auto NEIGHBORS=5` |
| `GENOME_DIR` | `data/genomes`; curated/automatic FASTA/GFF directory | `make plot-auto GENOME_DIR=data/custom/genomes` |
| `GENOME_METADATA` | `data/genomes.tsv`; curated/automatic species index | `make plot-auto GENOME_METADATA=data/custom/genomes.tsv` |
| `GENE_ANNOTATIONS` | `data/genes.tsv`; curated only | `make plot-curated GENE_ANNOTATIONS=data/custom/genes.tsv` |
| `MODE` | `default` or `replay`; reconstructed-tree mode | `make verify-auto MODE=replay` |
| `TARGET` | `verify`; target selected by `make submit` | `make submit TARGET=plot-main` |
| `SUBMIT_DRY_RUN` | `0`; set `1` to print submission JSON without writes or submission | `make -s submit TARGET=plot-auto SUBMIT_DRY_RUN=1` |

Combined examples, after preparing the indicated custom files:

```bash
make plot-auto THREADS=56 TARGET_GENE=Nramp NEIGHBORS=10 \
  GENOME_DIR=data/custom/genomes GENOME_METADATA=data/custom/genomes.tsv
make plot-curated THREADS=56 GENOME_DIR=data/custom/genomes \
  GENOME_METADATA=data/custom/genomes.tsv GENE_ANNOTATIONS=data/custom/genes.tsv
```

### Gene and homology-link palettes

| `GENE_PALETTE` | Meaning |
| --- | --- |
| `main-muted` (default) | Historical main saturation multiplied by 0.70, extended to 42 colors |
| `automatic` | Seven hue families, six tones each |
| `nature` | NPG-inspired colors extended with lighter tones |
| `cell` | Project-defined soft Cell-inspired colors |
| `science` | AAAS-inspired colors extended with lighter tones |

Journal names describe inspiration. `main-muted` keeps the historical main
mapping. Other main palettes rank ordinary categories by count descending, then
name descending in C byte order. Curated/automatic rank the same way and use the
first 42 colors; overflow/unknown annotations become gray `Other`. Automatic
labels prefer gene names, falling back to products. Targets remain white with
black dots; main `Unknown` remains black. Genus circles are 2.2 mm.

`LINK_PALETTE` accepts `blue` (default), `green`, `purple`, `orange` (three pale
tones), or `original` (historical grays). `LINK_COLORS` overrides the preset with
three distinct comma-separated `#RRGGBB` values ordered by BLAST identity:
**>=90%, 80–<90%, 70–<80%**.

```bash
make plot-main GENE_PALETTE=main-muted LINK_PALETTE=blue
make plot-curated GENE_PALETTE=automatic LINK_PALETTE=green
make plot-auto GENE_PALETTE=nature LINK_PALETTE=purple
make plot-main GENE_PALETTE=cell LINK_PALETTE=orange
make plot-auto GENE_PALETTE=science LINK_PALETTE=original
make plot-auto LINK_COLORS="#AECDE3,#D0E1EE,#EBF3F8"
```

### Executables and scheduler settings

| Parameter | Default / meaning | Example |
| --- | --- | --- |
| `PYTHON` | `python3` on PATH | `make check PYTHON=/opt/python/bin/python3` |
| `RSCRIPT` | Repository R, then PATH | `make plot-main RSCRIPT=/opt/R/bin/Rscript` |
| `BEDTOOLS` | Repository binary, then PATH | `make prepare-auto BEDTOOLS=/opt/bedtools/bin/bedtools` |
| `BLASTN` | Repository BLAST, then PATH | `make blast-auto BLASTN=/opt/blast/bin/blastn` |
| `MAKEBLASTDB` | Repository BLAST, then PATH | `make blast-auto MAKEBLASTDB=/opt/blast/bin/makeblastdb` |
| `MAFFT` | Repository MAFFT, then PATH | `make phylogeny-auto MAFFT=/opt/mafft/bin/mafft` |
| `IQTREE` | Repository IQ-TREE, then PATH | `make phylogeny-auto IQTREE=/opt/iqtree/bin/iqtree3` |
| `SLURM_PARTITION` | Unset; scheduler default | `make submit TARGET=plot-auto SLURM_PARTITION=compute` |
| `SLURM_NODE` | Unset; scheduler chooses | `make submit TARGET=plot-auto SLURM_NODE=nodeA` |
| `SLURM_TIME` | Unset; scheduler default | `make submit TARGET=plot-auto SLURM_TIME=04:00:00` |
| `SLURM_MEMORY` | Unset; scheduler default | `make submit TARGET=plot-auto SLURM_MEMORY=32G` |

Versions must match [DEPENDENCIES.md](DEPENDENCIES.md). Supply both BLAST
executables for a separate installation. Overrides apply to setup, checks, local
execution and submission.

## Stages, replay and redraw

Stages: `prepare-curated`, `blast-curated`, `phylogeny-curated`, and corresponding
`*-auto` targets. Prerequisites are built automatically. `verify-main`,
`verify-curated` and `verify-auto` validate an entry. `check-data` checks input
hashes; `check` additionally checks all tools.

```bash
make prepare-auto THREADS=56 TARGET_GENE=Nramp NEIGHBORS=10
make blast-auto THREADS=56 TARGET_GENE=Nramp NEIGHBORS=10
make phylogeny-auto THREADS=56 TARGET_GENE=Nramp NEIGHBORS=10
# Record an initial run before replaying its seed.
make verify-auto THREADS=56
make verify-auto THREADS=56 MODE=replay
```

MAFFT/IQ-TREE keep default analysis settings apart from resource flags. Initial
random seeds, selected models and commands are recorded. Replay needs that
entry's completed `phylogeny/replay_manifest.json`; results are written to
`phylogeny_replay/`. Compatible tree caches are reused. Different inputs,
scripts, tool builds or threads invalidate caches; preserve old results in a
separate working copy before changing settings. Different hardware/parallelism
may yield different inferred trees; byte-identical trees are not guaranteed.

Redraw existing results without extraction, BLAST or inference:

```bash
Rscript --vanilla scripts/plot_synteny.R --entry auto --target-gene Nramp \
  --gene-palette main-muted --link-palette blue
Rscript --vanilla scripts/verify_results.R --entry auto --target-gene Nramp \
  --gene-palette main-muted --link-palette blue
```

Both commands accept `--entry main|curated|auto`, `--target-gene NAME`,
`--gene-palette NAME`, `--link-palette NAME`, and `--link-colors HEX,HEX,HEX`.
`--input-dir` selects existing reconstructed tables; `--phylogeny-dir` selects
their tree, including replay outputs:

```bash
Rscript --vanilla scripts/plot_synteny.R --entry auto \
  --input-dir results/auto --phylogeny-dir results/auto/phylogeny_replay \
  --gene-palette nature --link-colors "#AECDE3,#D0E1EE,#EBF3F8"
```

Use identical options for plotting and verification. These R options do not
change extraction or inference.

## Outputs and notes

Each entry exports `synteny.pdf`, `synteny.png`, `caption.txt`,
`gene_colors.tsv`, `link_colors.tsv`, `figure_style.tsv` and audit tables.
PDFs are vector graphics with embedded fonts; PNG previews use 600 dpi.
Reconstructed entries additionally export `legend.pdf`, extracted FASTA,
normalized annotations, BED regions, BLAST, alignments, trees, species mappings
and execution provenance.

- Species labels come from assembly metadata, not taxonomic inference.
- GFF is 1-based inclusive; BED is 0-based half-open. Automatic local coordinates
  are 1-based inclusive; curated retains its historical offset.
- No automatic target gives a diagnostic table and an error. Old scientific
  outputs do not represent that new query.
- Insufficient CPUs, missing source features, mismatched versions, missing replay
  manifests and incompatible caches produce errors. Correct the input/runtime or
  use a fresh working copy for different data/settings.
- Results and installed environments are Git-ignored. Share needed results
  separately; tracked inputs, code and locks can rebuild all three entries.
