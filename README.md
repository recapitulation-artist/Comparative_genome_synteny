# Cd Uptake Synteny

Extract 29 uptake-gene neighborhoods from the bundled full genome assemblies
and GFF annotations with BEDTools, reproduce the supplied reference FASTA exactly,
and render the original phylogeny and synteny figure from the prepared plot inputs.
All paths resolve within this repository; original project files are untouched.

## Run

Extraction requires Linux x86_64, Python >=3.9, Make, curl and sha256sum.
The first command downloads the official BEDTools 2.31.0 static binary into
`.tools/`, checking its locked SHA256 before execution. No system installation
is performed. Subsequent extraction runs work offline.

```sh
make setup-bedtools
make extract
```

To restore the R environment and run both extraction and plotting:

```sh
conda create --name cd_uptake_synteny --file environment/conda_linux_64.lock
conda activate cd_uptake_synteny
make setup
make verify
```

The Conda lock requires Linux x86_64 with glibc >=2.28 and <3.0. It supplies
R 4.6.0, compilers, Cairo and fonts. `make setup-r` installs 23 exact CRAN package
versions into `.r_library/`. Initial environment restoration requires internet
access. Python >=3.9 must also be available; extraction uses only its standard
library and records the actual Python version.

With an existing matching R environment, use
`make verify RSCRIPT=/path/to/Rscript`. `make plot` renders only the figure;
`make check` checks all input hashes and software locks. `PYTHON` and `BEDTOOLS`
can also be overridden; the BEDTools artifact must match the binary lock.

## Extraction method

`data/sequences.csv` defines region order and lengths. The 29 `Uptake-Gene` rows
in `data/genes.tsv` identify the assembly and target protein. Each target must
match exactly one GFF CDS by `protein_id`. All 579 legacy annotations are also
matched independently: 476 CDS by protein accession and 103 RNA genes by locus
tag. Their contigs, strands and coordinate offsets must agree.

The supplied reference has a historical coordinate offset on both strands.
To reproduce it, the workflow uses:

```text
BED_start = GFF_target_start_1based - legacy_target_start + 1
BED_end   = BED_start + region_length
```

These boundaries are independently checked against every mapped legacy feature.
Each reference region begins one nucleotide after the lowest genomic feature
start and ends at the greatest feature end. Thus its first annotated feature
is clipped by one nucleotide. This behavior is retained for exact reproduction;
new coordinate tables use standard conventions and record clipping explicitly.
All 29 regions retain the forward genomic orientation, including those with a
reverse-strand target. No fixed flank distance is applied.

[BEDTools getfasta](https://bedtools.readthedocs.io/en/stable/content/tools/getfasta.html)
extracts each genome-derived interval. The reference FASTA is read only after
extraction for validation.
[BEDTools intersect](https://bedtools.readthedocs.io/en/stable/content/tools/intersect.html)
recovers every overlapping GFF gene or pseudogene on either strand, including
RNA genes. Neighbor tables retain full genomic coordinates and clipped local
coordinates. Source annotations identify the targets as Nramp divalent-metal
transporters; `Uptake-Gene` is the supplied analysis label.

## Files and outputs

```text
config/       Plot colors and expected plot metrics
data/         Prepared plot inputs, reference FASTA, genome index and checksums
data/genomes/ Full genome FASTA and GFF files, losslessly compressed
environment/  BEDTools binary lock, Conda runtime lock and R package lock
scripts/      Setup, extraction, plotting and verification
results/      Generated outputs, ignored by Git
```

Input provenance and coordinate conventions are in [`data/README.md`](data/README.md).
Only the two necessary genome files are included per assembly; download logs and
redundant genome formats are omitted. Assembly status and biological notes are
retained in `data/genomes.tsv`, including the suppressed assembly GCF_003460245.1.

| Output in `results/` | Content |
| --- | --- |
| `extracted_29_sequences.fasta` | Genome-derived reproduction of the original FASTA |
| `regions.tsv` | Region boundaries, target positions and source annotations |
| `neighbors.tsv` | All overlapping genes, source coordinates and clipping flags |
| `regions/*.bed` | Reusable region and gene BED intervals per assembly |
| `extraction_validation.tsv` | Per-record reference comparisons and sequence hashes |
| `extraction_metrics.tsv`, `extraction_versions.tsv` | Extraction counts, FASTA hash and runtime versions |
| `synteny.pdf`, `synteny.png` | Vector figure and 4,200 by 4,500 pixel raster figure |
| `sequence_order.csv`, `plotted_blast_hits.csv`, `metrics.tsv` | Plot order, displayed alignments and plot counts |
| `session_info.txt` | Actual R package environment |

The plot uses the supplied historical annotation, BLAST and tree inputs. It
ladderizes the topology without branch lengths or an inferred biological root.
Only hits between adjacent tracks with identity >=70% are drawn. Gene arrows
follow the supplied strands; red stars mark the uptake genes. Alignment, BLAST
execution and tree inference are outside this workflow.

## Verification

`make extract` checks input hashes, the locked BEDTools binary, all annotation
mappings, region boundaries and reference equality. It requires identical FASTA
headers, order, lengths, nucleotide sequences and complete file bytes.
`make verify` additionally checks the R package lock, plot coordinates and
expected counts: 29 sequences, 579 prepared annotations, 29 uptake markers and
606 overlapping genomic genes, plus 135 displayed BLAST hits. PDF and PNG
formats and PNG dimensions are checked.
Fonts and PDF timestamps may vary across machines.

Verified with the matching local environments and a separate export of
Git-tracked files; all 64 extraction outputs were byte-identical. Deliberate
coordinate, strand, duplicate-annotation and reference errors were rejected.
Fresh network restoration of the R environment has not been tested.
