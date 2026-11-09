"""Scientific selection checks independent of the repository's historical fixtures."""

from pathlib import Path
import csv
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from extract_regions import Feature, annotation_indices
from prepare_entry import coding_gene, load_inputs, match_target, select_neighbors


def gene(identifier, start, contig="contig", strand="+", **attrs):
    return Feature(contig, "gene", start, start + 8, strand,
                   dict(ID=identifier, gene_biotype="protein_coding", **attrs))


def cds(parent, part=0, **attrs):
    return Feature(parent.contig, "CDS", parent.start + part, parent.end, parent.strand,
                   dict(ID=f"cds-{parent.identifier}-{part}", Parent=parent.identifier, **attrs))


class NeighborSelectionTests(unittest.TestCase):
    def test_counts_parent_genes_once_and_keeps_genomic_order_on_negative_strand(self):
        genes = [gene(f"g{i:02}", 20 * i + 1, strand="-") for i in range(25)]
        features = genes + [cds(parent) for parent in genes] + [cds(genes[12], 2)]
        _, children = annotation_indices(features)
        selected, before, after = select_neighbors(genes, genes[12], children, 10)
        self.assertEqual([parent.identifier for parent in selected], [f"g{i:02}" for i in range(2, 23)])
        self.assertEqual((before, after), (12, 12))
        self.assertEqual(sum(parent is genes[12] for parent in selected), 1)

    def test_contig_boundary_does_not_rebalance_or_cross_contigs(self):
        genes = [gene(f"a{i}", 20 * i + 1) for i in range(8)]
        genes += [gene(f"b{i}", 20 * i + 1, contig="other") for i in range(20)]
        _, children = annotation_indices(genes + [cds(parent) for parent in genes])
        selected, before, after = select_neighbors(genes, genes[1], children, 3)
        self.assertEqual([parent.identifier for parent in selected], ["a0", "a1", "a2", "a3", "a4"])
        self.assertEqual((before, after), (1, 6))

    def test_excludes_rna_and_pseudogenes_even_when_they_have_cds(self):
        valid, rna, pseudo, flagged, child_flagged = [gene(name, i * 20 + 1) for i, name
                                                    in enumerate(["valid", "rna", "pseudo", "flag", "child"]) ]
        rna.attrs["gene_biotype"] = "tRNA"
        pseudo.kind = "pseudogene"
        flagged.attrs["pseudo"] = "true"
        features = [valid, rna, pseudo, flagged, child_flagged]
        features += [cds(parent) for parent in features]
        features[-1].attrs["pseudo"] = "true"
        _, children = annotation_indices(features)
        self.assertEqual([coding_gene(parent, children) for parent in features[:5]],
                         [True, False, False, False, False])
        self.assertFalse(coding_gene(gene("without-cds", 200), children))

    def test_zero_neighbors_keeps_target_only(self):
        target = gene("target", 10)
        _, children = annotation_indices([target, cds(target)])
        selected, before, after = select_neighbors([target], target, children, 0)
        self.assertEqual(selected, [target])
        self.assertEqual((before, after), (0, 0))

    def test_nramp_alias_and_product_match_are_case_insensitive(self):
        named = gene("named", 1, gene="MNTH")
        product = gene("product", 20)
        _, children = annotation_indices([named, product, cds(named),
                                          cds(product, product="NRAMP family divalent metal transporter")])
        self.assertTrue(match_target(named, children, "nramp"))
        self.assertTrue(match_target(product, children, "Nramp"))
        self.assertTrue(match_target(named, children, "mntH"))
        self.assertFalse(match_target(product, children, "Nr"))

    def test_multiple_targets_remain_distinct_parent_loci(self):
        targets = [gene("one", 1, gene="mntH"), gene("two", 100, gene="mntH")]
        _, children = annotation_indices(targets + [cds(parent) for parent in targets])
        matched = [parent for parent in targets if match_target(parent, children, "Nramp")]
        self.assertEqual(len(matched), 2)
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            select_neighbors(targets, targets[0], children, -1)


class InputIsolationTests(unittest.TestCase):
    def create_inputs(self, directory, rows):
        genome_dir = directory / "genomes"
        assembly_dir = genome_dir / "assembly"
        assembly_dir.mkdir(parents=True)
        (assembly_dir / "genome.fna").write_text(">contig\n" + "A" * 100 + "\n")
        (assembly_dir / "annotation.gff").write_text(
            "##gff-version 3\n"
            "contig\tRefSeq\tgene\t1\t9\t.\t+\t.\tID=gene-one;gene=unrelated;gene_biotype=protein_coding\n"
            "contig\tRefSeq\tCDS\t1\t9\t.\t+\t0\tID=cds-one;Parent=gene-one;protein_id=protein-one;product=unrelated\n")
        metadata = directory / "metadata.tsv"
        with metadata.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        return SimpleNamespace(entry="auto", genome_dir=genome_dir, genome_metadata=metadata,
                               annotations=directory / "missing-manual.tsv")

    def test_auto_metadata_does_not_require_legacy_sequence_ids_or_manual_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            args = self.create_inputs(Path(temporary), [dict(assembly_accession="assembly", organism_name="Species name")])
            metadata, manual, _, _, hashes = load_inputs(args)
            self.assertIsNone(manual)
            self.assertNotIn("seq_id", metadata[0])
            self.assertFalse(any("missing-manual" in row["file"] for row in hashes))

    def test_conflicting_metadata_for_same_assembly_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            args = self.create_inputs(Path(temporary), [dict(assembly_accession="assembly", organism_name="Species one"),
                                                      dict(assembly_accession="assembly", organism_name="Species two")])
            with self.assertRaisesRegex(ValueError, "Conflicting metadata"):
                load_inputs(args)

    def test_no_target_preserves_previous_outputs_and_reports_their_status(self):
        root = Path(__file__).resolve().parents[1]
        binary = root / ".tools/bedtools"
        if not binary.is_file():
            self.skipTest("Locked BEDTools is not installed")
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            args = self.create_inputs(directory, [dict(assembly_accession="assembly", organism_name="Species name")])
            output = directory / "output"
            (output / "regions").mkdir(parents=True)
            previous_fasta, previous_bed = output / "extracted_sequences.fasta", output / "regions/old.bed"
            previous_fasta.write_text("previous FASTA\n")
            previous_bed.write_text("previous BED\n")
            result = subprocess.run([sys.executable, str(root / "scripts/prepare_entry.py"), "--entry", "auto",
                                     "--genome-dir", str(args.genome_dir), "--genome-metadata", str(args.genome_metadata),
                                     "--annotations", str(args.annotations), "--bedtools", str(binary),
                                     "--output", str(output)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("previous scientific outputs", result.stderr)
            self.assertTrue((output / "target_matches.tsv").is_file())
            self.assertEqual(previous_fasta.read_text(), "previous FASTA\n")
            self.assertEqual(previous_bed.read_text(), "previous BED\n")

    def test_output_cannot_overwrite_supplied_metadata(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            args = self.create_inputs(directory, [dict(assembly_accession="assembly", organism_name="Species name")])
            metadata = directory / "genomes.tsv"
            args.genome_metadata.rename(metadata)
            original = metadata.read_bytes()
            result = subprocess.run([sys.executable, str(root / "scripts/prepare_entry.py"), "--entry", "auto",
                                     "--genome-dir", str(args.genome_dir), "--genome-metadata", str(metadata),
                                     "--output", str(directory)], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("overwrite", result.stderr)
            self.assertEqual(metadata.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
