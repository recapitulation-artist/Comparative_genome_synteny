# Cd Uptake Synteny

A reproducible workflow from 29 bundled genome assemblies to uptake-gene
neighborhoods, local BLAST comparisons, multiple-sequence alignments, region
phylogenies and two synteny figures. Original project files remain unchanged.
Analysis paths resolve relative to the repository root.

## Two plotting entry points

Both figures use **the same manually curated `data/genes.tsv`**, a byte-identical
copy of `29gene2.txt` (579 annotations, including 29 uptake targets). They also
use the same extracted regions and newly calculated BLAST hits.

| Make target | Tree | Output directory |
| --- | --- | --- |
| `plot-reference` | Supplied `data/phylogeny.nwk` | `results/reference/` |
| `plot-inferred` | Newly inferred default-parameter tree | `results/inferred/` |

Tree topology determines row order and therefore adjacent BLAST links. The
plot does not substitute automatically recovered GFF neighbors for the curated
annotations. Species labels come from `data/genomes.tsv`; unique sequence IDs
remain visible. Colors retain the eight recorded genera.

## Figure layout

Both entry points share `config/figure_style.tsv` and the same renderer. The
layout follows the [Nature figure specifications](https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/)
and [sizing guidance](https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/):
183 x 170 mm, 5.2–6.5 pt black sans-serif lettering, 0.25–0.5 pt strokes,
editable vector PDF and a 600 dpi PNG preview. Cairo may round the PDF media
box down by less than 0.4 mm. Liberation Sans 2.1.5, an Arial-metric-compatible
sans-serif family, is bundled under `environment/fonts/` with its SIL Open Font
License and SHA256 checksums. The regular, italic and bold faces come from the
[official release](https://github.com/liberationfonts/liberation-fonts/releases/tag/2.1.5).
`config/fonts.conf` selects these files on every node. The verifier inspects
the exported Cairo PDF for all three embedded TrueType faces and Unicode
mappings; silent font substitution fails verification. When Poppler's
`pdfinfo` and `pdffonts` are available, it also independently checks the
page dimensions, single-page export and font embedding.

The left tree is a cladogram: row order and topology are preserved, branch
lengths and node supports are not displayed, and no biological root is inferred.
The 22 mm tree area is separated from italic full species names and upright
unique sequence IDs. Tree and gene rows share the same vertical coordinates.

Genus markers use the Wong eight-colour palette and BLAST identity uses three
grey levels. Gene-category colours retain the original annotation mapping;
28 categories should not be interpreted using colour alone. Arrow keys and
labels identify each category, and a small black dot marks the target Nramp
transporter. The coordinate axis reports position within each extracted region
in kb. Each entry also writes `caption.txt` and `figure_style.tsv` alongside its
figures, with the style configuration included in the input-hash manifest.

## Setup and execution

The workflow requires Linux x86_64, Python >=3.9, Make and the locked tools.
It uses BEDTools 2.31.0, BLAST+ 2.17.0+, MAFFT 7.526, IQ-TREE 3.1.4 and R 4.6.0.
The IQ-TREE installation is the local stable release, not a beta build.

```sh
# Restore R when a matching installation is not available.
conda create --name cd_uptake_synteny --file environment/conda_linux_64.lock
conda activate cd_uptake_synteny
make setup

# Submit both plot entry points and all verification to node2.
make submit MODE=default

# Submit either plot entry point separately.
make submit TARGET=plot-reference
make submit TARGET=plot-inferred MODE=default

# Rebuild using the recorded first-run random seeds.
make submit MODE=replay
```

`make submit` requests the Slurm `normal` partition, node2, one task and 96 CPUs,
waits for completion, and writes `results/slurm/<job_id>.log`. Parallel tool
stages are configured for 96 threads; the two tree analyses run sequentially.
Some operations are inherently serial or have fewer available parallel tasks.
The workflow checks the node and allocation before calculation. Node2 SSH may
be denied before a Slurm allocation is active; submission does not require SSH.

Inside that allocation, `make plot-reference`, `make plot-inferred` and
`make verify` are the corresponding direct commands. `make run` generates both
figures; `make plot` remains an alias for `plot-reference`. `make extract` can
be used separately without a cluster allocation.

`make setup-bedtools` verifies the official static binary against its SHA256.
`make setup-blast` and `make setup-phylogeny` use matching local installations
or restore exact Conda package locks under ignored `.tools/`. The phylogeny
lock contains only the required IQ-TREE/MAFFT runtime dependencies.
`make setup-r` installs the 23 locked R packages into `.r_library/`.
Initial restoration needs internet access; bundled analysis data need none.

The submission helper recognizes the existing local R environment. To select
other matching tools, set `RSCRIPT`, `BEDTOOLS`, `BLASTN`, `MAKEBLASTDB`, `MAFFT`
and `IQTREE`. Example: `make submit RSCRIPT=/path/to/Rscript`.
`make check RSCRIPT=/path/to/Rscript` checks data hashes and software versions
without submitting a job. Environment restoration on a separate machine must
meet the virtual requirements recorded in each lock.

## Extraction and BLAST

`data/sequences.csv` defines the order and lengths of the 29 regions.
The `Uptake-Gene` rows in the curated annotations identify each assembly and
target protein. All 579 historical features are independently matched to GFF:
476 CDS by protein accession and 103 RNA genes by locus tag.

The reference has a historical one-base offset, retained for exact reproduction:

```text
BED_start = GFF_target_start_1based - legacy_target_start + 1
BED_end   = BED_start + region_length
```

All matched features must agree on these boundaries. Each region starts one
base after its lowest mapped feature start and ends at the greatest mapped
feature end. All regions retain forward genomic orientation; no fixed flank
size is imposed. BEDTools `getfasta -nameOnly` performs extraction, and
`intersect` recovers all overlapping source genes and pseudogenes separately.
The recovered neighbors are an audit output, not the plotting annotations.

The resulting 29 sequences contain 567,757 bases and reproduce the reference
FASTA bytes exactly. Source target products are annotated Nramp divalent-metal
transporters; `Uptake-Gene` is the supplied analysis label.

`config/blast.json` specifies the search parameters. Local `makeblastdb` builds
one database from the extracted regions and `blastn -task megablast` searches
those same regions with 96 threads. Other settings are word size 28, reward 1,
penalty -2, gap-open/extension 0, DUST `20 64 1`, soft masking, both strands,
E-value threshold 10 and `max_target_seqs=500`. Output retains only pairs with
query rank less than subject rank, preserving directed endpoints and HSP order.

The 4,513 raw HSPs yield 2,189 rows across all 406 unordered pairs. All 11 fields
other than E-value match the historical table. The 994 nonzero E-values differ:
original database/search statistics were not supplied. Current values remain
unchanged, with detailed differences reported. Plotting uses identities and
coordinates, so these statistical differences do not affect its links.

## Alignment, trees and species

The workflow converts full extracted FASTA headers to stable legacy `seq_id`
values through `results/regions.tsv`, never by truncating strings. Repeated
species names are metadata and do not replace unique tip identifiers.

Two branches are generated under `results/phylogeny/`:

- `default/`: original forward-genomic regions, used by `plot-inferred`.
- `oriented/`: a diagnostic control reverse-complementing the 12 regions whose
  GFF target strand is negative, so all targets have a common direction.

Both use MAFFT's default FFT-NS-2 alignment and IQ-TREE's default model selection
and tree search:

```sh
mafft --thread 96 --threadtb 96 input.fasta > aligned.fasta
iqtree3 -s aligned.fasta -T 96 --prefix tree
```

The thread options change resources only. No trimming, automatic direction
adjustment, specified substitution model, outgroup or support replicates are
added. The initial `MODE=default` run does not specify a seed. Its actual seeds,
model choices, versions, commands, alignment statistics and hashes are recorded.
The oriented tree is a diagnostic result, not a third plotting entry point.

`MODE=replay` uses seeds and expected hashes in
`config/phylogeny_replay.json` and writes trees under `results/phylogeny_replay/`,
leaving the initial trees available for comparison. Its inferred synteny figure
is written to the same `results/inferred/` entry-point directory. Completed
matching runs are cached; an interrupted matching run can resume from its
IQ-TREE checkpoint. To recompute independently,
use a clean export or a new `--output-root` with `scripts/run_phylogeny.py`.
Existing incompatible outputs cause an error instead of silent replacement.
`make record-phylogeny` records the completed initial runs in the replay manifest;
it refuses to replace a different recorded baseline. Tree byte differences during
replay are reported alongside topology and branch-length comparisons.

Species and strains are inherited from bundled assembly metadata: 13 species
and eight genera. Assembly status and notes are retained, including the
suppressed assembly GCF_003460245.1. Generated region trees are not described as
species trees. Annotated Newick files retain distinct sequence IDs, and PDF/PNG
phylograms show branch lengths and species/strain labels.

The original alignment, software settings, seed and run logs were not supplied.
The reference has numeric node labels, while default inference does not request
support tests. Exact historical file reproduction is therefore assessed rather
than assumed. Reports compare unrooted splits/RF distance, shared splits,
406 pairwise patristic distances and branch lengths. Raw historical labels >=95
are reported without assigning an unknown support-test method.
The reference's longest internal branch separates the 17 positive-target-strand
regions from the 12 negative-strand regions. The diagnostic control explicitly
checks whether this partition persists after direction normalization.

The initial default tree shares 24 of 26 nontrivial splits with the reference
(RF = 4; normalized RF = 0.0769); it is not an exact reconstruction. All 18
reference splits carrying raw numeric labels >=95 are retained. The two absent
reference splits carry labels 64 and 86, with unknown support-test provenance.
The oriented tree shares only five splits with either unnormalized tree
(RF = 42) and removes the exact 17-positive/12-negative partition. These
comparisons indicate strong orientation sensitivity, without establishing a
species phylogeny or proving the cause of the historical topology. The two
alignments have 41,061 and 35,844 columns, with fitted models TVM+F+G4 and
GTR+F+I+R5, respectively; branch-length comparisons remain descriptive.

## Outputs and verification

| Location under `results/` | Content |
| --- | --- |
| `extracted_29_sequences.fasta` | Genome-derived, byte-identical reference reproduction |
| `regions.tsv`, `neighbors.tsv`, `regions/` | Source coordinates, target mapping and BED intervals |
| `extraction_*.tsv` | Mapping/reference validation, counts and versions |
| `blast/` | Raw HSPs, commands, versions and historical comparison |
| `reference/`, `inferred/` | Synteny PDF/PNG, caption, figure style, plot order, links, metrics and input hashes |
| `phylogeny/species.tsv` | Explicit ID, species, strain, assembly and strand mapping |
| `phylogeny/{default,oriented}/` | Inputs, alignments, IQ-TREE outputs, annotated trees and run manifests |
| `phylogeny/comparison.tsv`, `phylogeny/summary.md` | Historical-tree and direction-control comparisons |
| `phylogeny_replay/` | Seeded replay trees and comparisons with the recorded initial trees |
| `validation/clean_replay.json` | Observed clean-export hashes, topology comparisons and Slurm status |

`data/blast_hits.csv` is generated and ignored by Git. Raw inputs and references
are tracked and covered by `data/checksums.sha256`; generated results and logs
are ignored. Only full assembly FASTA and GFF are bundled per genome.

Verification checks exact extraction bytes and all annotation mappings; BLAST
geometry/order; alignment length, ID sets and ungapped-input equality; tree tip
sets and species labels; plotting coordinates, annotation counts and image
formats. The reference figure retains 135 links. The inferred figure derives
its own adjacent links from the inferred topology. Both figures retain all 579
manual annotations and 29 uptake markers.

The full workflow passed on node2 with 96 CPUs, then passed again from a clean
Git-index export containing no generated analysis inputs. The extracted FASTA,
BLAST CSV, both alignments and the then-current synteny PNGs were byte-identical.
The publication layout was subsequently rerendered and verified for both entry
points, retaining the same row order, annotation counts and links. Both new
PNGs were also reproduced byte-for-byte from clean copies of the plotting
inputs, using the bundled fonts from an unrelated working directory. Both seeded
tree replays retained all 26 initial splits (RF = 0). Their 406-distance RMSEs
were 4.31e-9 (default) and 6.38e-6 (oriented); Newick bytes differed, so exact
numerical tree reproduction is not claimed. Full evidence is retained under
`results/validation/` and `results/phylogeny_replay/`.

The 15-package phylogeny runtime lock was also restored successfully into a fresh
isolated prefix; both tools passed version and dynamic-library checks. PDF
timestamps, fonts and floating-point results may vary across environments;
versions, input hashes and topology/branch comparisons make these differences
visible.
