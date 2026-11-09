# Inputs

| File | Content |
| --- | --- |
| `sequences.csv` | Supplied sequence identifiers, order and region lengths |
| `genes.tsv` | Byte-identical copy of curated `29gene2.txt` |
| `reference_blast_hits.csv` | Supplied historical BLAST hits |
| `phylogeny.nwk` | Supplied reference tree |
| `reference_sequences.fasta` | Historical FASTA for optional comparison |
| `genomes.tsv` | Assembly accessions and species metadata |
| `genomes/GCF_*/genome.fna.gz` | Full raw assembly FASTA |
| `genomes/GCF_*/annotation.gff.gz` | Original assembly GFF |

`checksums.sha256` covers bundled inputs. Genome gzip files decompress to the
original source bytes. Generated files are written under `results/`.

Main uses the supplied tables and tree. Curated uses only raw genomes, genome
metadata and `genes.tsv`. Automatic uses only raw genomes, genome metadata and
the requested target name. Neither reconstructed entry reads historical BLAST,
reference tree, reference FASTA or the supplied sequence index.

GFF coordinates are 1-based inclusive; BED intervals are 0-based half-open.
Curated coordinates retain the original one-base start offset for exact FASTA
reproduction. Automatic local coordinates are true 1-based inclusive.
