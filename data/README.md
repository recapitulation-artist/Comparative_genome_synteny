# Inputs

The four inputs are byte-identical copies with standardized filenames. SHA256
checksums are recorded in `checksums.sha256`; original files remain untouched.

| File | Original filename | Content |
| --- | --- | --- |
| `sequences.csv` | `29seqs2.csv` | 29 sequence identifiers and lengths |
| `genes.tsv` | `29gene2.txt` | 579 gene annotations |
| `blast_hits.csv` | `29uptake_avad3.csv` | 2,189 BLAST hits |
| `phylogeny.nwk` | `aligned_29_sequences_high.fasta.treefile` | Supplied 29-tip Newick tree |

Required columns:

- Sequences: `seq_id`, `length`.
- Genes: `seq_id`, `start`, `end`, `strand`, `note`, `organisms` (genus).
- BLAST hits: `seq_id`, `seq_id2`, `pident`, `start`, `end`, `start2`, `end2`.

Coordinates are 1-based and inclusive. Reverse BLAST endpoints retain their
direction. Identifiers must match across inputs; each sequence has one
`Uptake-Gene` annotation. The repository starts from these prepared inputs;
sequence alignment, BLAST execution and tree inference are outside its scope.
