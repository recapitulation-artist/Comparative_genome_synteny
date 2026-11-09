# Inputs

All inputs are bundled. `checksums.sha256` covers every data file except itself
and this documentation. Original project files remain unchanged.

| File | Source | Content |
| --- | --- | --- |
| `sequences.csv` | `29seqs2.csv` | 29 sequence identifiers, order and region lengths |
| `genes.tsv` | `29gene2.txt` | 579 historical annotations and 29 uptake targets |
| `blast_hits.csv` | `29uptake_avad3.csv` | 2,189 BLAST hits |
| `phylogeny.nwk` | `aligned_29_sequences_high.fasta.treefile` | Supplied 29-tip tree |
| `reference_sequences.fasta` | `extracted_29_sequences.fasta` | Validation reference, never used to derive extraction coordinates |
| `genomes.tsv` | `Enterobacteriaceae_genomes/genomes.tsv` | Assembly accessions, sequence identifiers and biological metadata |
| `genomes/GCF_*/genome.fna.gz` | Assembly `*_genomic.fna` | Full assembly FASTA |
| `genomes/GCF_*/annotation.gff.gz` | Assembly `*_genomic.gff` | Full source GFF annotation |

The six top-level data files are byte-identical source copies. Genome inputs
are deterministic gzip copies whose decompressed bytes match their sources.
No downloads are needed for the data. `genomes.tsv` preserves source assembly
status, annotation dates, completeness notes and contig counts.

## Schemas

- Sequences: `seq_id`, `length`; row order defines output FASTA order.
- Genes: `seq_id`, `start`, `end`, `strand`, `type`, `feat_id`, `geom_id`,
  `note`, `organisms`. `geom_id` is the assembly accession; `feat_id` is a
  protein accession for CDS or a locus tag for RNA genes. Exactly one row per
  sequence has `note=Uptake-Gene`.
- BLAST: `seq_id`, `seq_id2`, `pident`, `start`, `end`, `start2`, `end2`.
  Endpoints retain their direction, including reverse alignments.
- Genomes: `assembly_accession` maps each legacy sequence to its bundled
  FASTA/GFF directory. The FASTA output identifier uses the full GFF contig
  accession followed by `_` and the target protein accession. Legacy `seq_id`
  values can be truncated contig prefixes.

## Coordinates

The historical `genes.tsv` starts are one nucleotide larger than true local
1-based starts in the reference FASTA; its ends equal true local inclusive ends.
This applies to all 579 records on both strands. A historical `start=1`
therefore corresponds to an annotation extending one nucleotide before the
reference region. These original bytes are retained for plotting.

Extraction anchors to source GFF coordinates, which are 1-based inclusive, and
converts to 0-based, half-open BED intervals. All mapped legacy annotations must
independently yield the same BED start:

```text
BED_start = GFF_start - legacy_start + 1 = GFF_end - legacy_end
BED_end   = BED_start + sequences.length
```

The generated `regions.tsv` and `neighbors.tsv` explicitly label coordinate
systems. Neighbor genomic spans remain complete; local spans are clipped to the
extracted region and clipping flags identify truncated features. Neighbor rows
come from all overlapping source GFF genes and pseudogenes and may differ in
count from the curated historical plot annotations.
