# Cd Uptake Synteny

Extract 29 uptake-gene neighborhoods from the bundled full genome assemblies
and GFF annotations with BEDTools, reproduce the supplied reference FASTA exactly,
run local BLAST comparisons, and render the phylogeny and synteny figure.
All paths resolve within this repository; original project files are untouched.

## Run

Extraction requires Linux x86_64, Python >=3.9, Make, curl and sha256sum.
The first command downloads the official BEDTools 2.31.0 static binary into
`.tools/`, checking its locked SHA256 before execution. No system installation
is performed. Subsequent extraction runs work offline.

```sh
make setup-bedtools
make setup-blast
make blast
```

`make blast` runs extraction first, builds a database from the 29 regions with
local `makeblastdb`, runs local `blastn`, and writes `data/blast_hits.csv`.
BLAST+ 2.17.0+ is required. `make setup-blast` uses a matching installation on
PATH or restores all 51 locked packages into `.tools/blast/` with Conda.

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
`make verify RSCRIPT=/path/to/Rscript`. `make extract` generates only the FASTA;
`make plot` runs extraction and BLAST before rendering the figure;
`make check` checks all input hashes and software locks. `PYTHON` and `BEDTOOLS`
can also be overridden; the BEDTools artifact must match the binary lock.
Use `BLASTN=/path/to/blastn MAKEBLASTDB=/path/to/makeblastdb` to select explicit
executables. Otherwise the workflow checks `.tools/blast/bin/` and then PATH.

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

## BLAST method

`config/blast.json` specifies all search parameters. The pipeline builds one
nucleotide database containing all 29 extracted regions and searches those same
regions with `blastn -task megablast`. It uses a 28-base word, reward 1, penalty
-2, gap-open and gap-extension costs 0, DUST `20 64 1`, soft masking, both
strands, E-value threshold 10, and one thread. It reports all HSPs and retains
only pairs with query rank less than subject rank in `data/sequences.csv`.
Self hits and reciprocal comparisons are removed; HSP order and directed
coordinates are retained. The 4,513 raw hits produce 2,189 rows across all
406 sequence pairs. No additional length or identity filter is applied here.

The original table is preserved as `data/reference_blast_hits.csv`. All generated
rows match it in order and in all 11 columns other than E-value. The 994 nonzero
E-values differ; the original BLAST command and database statistics were not
supplied. [E-values depend on database search space](https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html).
The workflow retains the E-values calculated by the current 29-sequence search
and reports the differences explicitly. It does not copy historical statistics
into the generated table. The figure uses identities and coordinates, so this
difference does not alter its links.

## Files and outputs

```text
config/       BLAST parameters, plot colors and expected metrics
data/         Source inputs, references and generated blast_hits.csv
data/genomes/ Full genome FASTA and GFF files, losslessly compressed
environment/  BEDTools binary lock, BLAST/R Conda locks and R package lock
scripts/      Setup, extraction, BLAST, plotting and verification
results/      Generated outputs, ignored by Git
```

Input provenance and coordinate conventions are in [`data/README.md`](data/README.md).
`data/blast_hits.csv` is generated and ignored by Git; the reference table and
all raw inputs are tracked. The generated table can be removed and rebuilt.
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
| `blast/` | Raw BLAST output, sequence IDs, parameters, versions and reference comparison |
| `synteny.pdf`, `synteny.png` | Vector figure and 4,200 by 4,500 pixel raster figure |
| `sequence_order.csv`, `plotted_blast_hits.csv`, `metrics.tsv` | Plot order, displayed alignments and plot counts |
| `session_info.txt` | Actual R package environment |

The plot uses the supplied historical annotation and tree with generated BLAST hits. It
ladderizes the topology without branch lengths or an inferred biological root.
Only hits between adjacent tracks with identity >=70% are drawn. Gene arrows
follow the supplied strands; red stars mark the uptake genes. Multiple-sequence
alignment and tree inference are outside this workflow.

## Verification

`make extract` checks input hashes, the locked BEDTools binary, all annotation
mappings, region boundaries and reference equality. It requires identical FASTA
headers, order, lengths, nucleotide sequences and complete file bytes.
`make verify` additionally checks the R package lock, plot coordinates and
expected counts: 29 sequences, 579 prepared annotations, 29 uptake markers and
606 overlapping genomic genes, plus 135 displayed BLAST hits. PDF and PNG
formats and PNG dimensions are checked.
Fonts and PDF timestamps may vary across machines.

BLAST verification requires matching identifiers, HSP order, all 11 non-E-value
fields and all 406 pairs, and records every E-value difference. The input hash
check excludes generated `data/blast_hits.csv`. Network restoration of the R
and BLAST environments has not been tested.

The complete workflow passed in the local environments and a clean export of
Git-tracked files without generated inputs. Rebuilt FASTA, BLAST outputs and
PNG matched the local results byte-for-byte.
