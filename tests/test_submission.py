"""Command and environment checks; no scientific tools or Slurm jobs are run."""

import contextlib
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import submit_workflow as submit


class SubmissionTests(unittest.TestCase):
    def args(self, *values):
        return submit.argument_parser().parse_args(values)

    def test_default_allocation_is_generic_and_56_cpus(self):
        command = submit.build_command(self.args())
        self.assertIn("--cpus-per-task=56", command)
        self.assertIn("--nodes=1", command)
        self.assertIn("--ntasks=1", command)
        self.assertIn("THREADS=56", command)
        self.assertIn("--wait", command)
        self.assertFalse(any(item.startswith(("--partition=", "--nodelist=", "--mem=", "--time="))
                             for item in command))
        self.assertFalse(any("node2" in item or "96" == item for item in command))

    def test_future_allocation_does_not_validate_login_affinity(self):
        # Submission can originate on a small login allocation; execution validates CPUs later.
        with patch.dict(os.environ, {"SLURM_CPUS_PER_TASK": "unrelated-login-allocation"}):
            self.assertIn("--cpus-per-task=56", submit.build_command(self.args("--threads", "56")))

    def test_optional_scheduler_settings_are_forwarded(self):
        command = submit.build_command(self.args("--partition", "compute", "--node", "worker7",
                                                  "--time", "02:00:00", "--memory", "32G"))
        for item in ("--partition=compute", "--nodelist=worker7", "--time=02:00:00", "--mem=32G"):
            self.assertIn(item, command)

    def test_thread_count_must_be_positive_integer(self):
        for value in ("0", "-56", "56.0", "many", True):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "positive integer"):
                submit.positive_threads(value)

    def test_data_tools_and_palette_values_remain_single_arguments(self):
        args = self.args("--target", "plot-auto", "--mode", "replay",
                         "--target-gene", "transporter's family", "--neighbors", "7",
                         "--genome-dir", "data/new genomes", "--genome-metadata", "data/new metadata.tsv",
                         "--annotations", "data/new annotations.tsv", "--rscript", "/opt/R env/bin/Rscript",
                         "--bedtools", "/opt/tools path/bedtools", "--blastn", "/opt/BLAST path/blastn",
                         "--makeblastdb", "/opt/BLAST path/makeblastdb", "--mafft", "/opt/MAFFT path/mafft",
                         "--iqtree", "/opt/IQ tree/iqtree3", "--gene-palette", "cell",
                         "--link-palette", "purple", "--link-colors", "#123456,#abcdef,#fedcba")
        command = submit.build_command(args, root=ROOT / "path with spaces")
        for item in ("MODE=replay", "TARGET_GENE=transporter's family", "NEIGHBORS=7",
                     "GENOME_DIR=data/new genomes", "GENOME_METADATA=data/new metadata.tsv",
                     "GENE_ANNOTATIONS=data/new annotations.tsv", "RSCRIPT=/opt/R env/bin/Rscript",
                     "BEDTOOLS=/opt/tools path/bedtools", "BLASTN=/opt/BLAST path/blastn",
                     "MAKEBLASTDB=/opt/BLAST path/makeblastdb", "MAFFT=/opt/MAFFT path/mafft",
                     "IQTREE=/opt/IQ tree/iqtree3", "GENE_PALETTE=cell", "LINK_PALETTE=purple",
                     "LINK_COLORS=#123456,#abcdef,#fedcba"):
            self.assertIn(item, command)

    def test_literal_dollars_are_protected_from_make_expansion(self):
        self.assertEqual(submit.make_assignment("TARGET_GENE", "name $(literal) $value"),
                         "TARGET_GENE=name $$(literal) $$value")

    def test_invalid_targets_and_palettes_are_rejected(self):
        for values in (("--target", "plot-reference"), ("--gene-palette", "unknown"),
                       ("--link-palette", "unknown")):
            with self.subTest(values=values), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                self.args(*values)

    def test_invalid_neighbor_target_and_custom_colors_are_rejected(self):
        for values in (("--neighbors", "-1"), ("--target-gene", " "),
                       ("--link-colors", "#123456,#123456,#abcdef"),
                       ("--link-colors", "#123456,#abcdef")):
            with self.subTest(values=values), self.assertRaises(ValueError):
                submit.build_command(self.args(*values))

    def test_inherited_sbatch_settings_are_removed(self):
        with patch.dict(os.environ, {"SBATCH_PARTITION": "old", "SBATCH_NODELIST": "node2",
                                     "SBATCH_CPUS_PER_TASK": "96", "KEEP_ME": "yes"}, clear=True):
            self.assertEqual(submit.submission_environment(), {"KEEP_ME": "yes"})

    def test_dry_run_is_read_only(self):
        with patch.object(sys, "argv", ["submit_workflow.py", "--dry-run", "--threads", "56"]), \
                patch.object(Path, "mkdir") as mkdir, patch.object(submit.subprocess, "run") as run, \
                contextlib.redirect_stdout(io.StringIO()) as output:
            submit.main()
        mkdir.assert_not_called()
        run.assert_not_called()
        payload = json.loads(output.getvalue())
        self.assertEqual(payload["threads"], 56)
        self.assertIn("THREADS=56", payload["argv"])


class MakeInterfaceTests(unittest.TestCase):
    def dry_commands(self, *args):
        result = subprocess.run(["make", "--no-print-directory", "-sn", *args], cwd=ROOT,
                                check=True, capture_output=True, text=True)
        return [shlex.split(line) for line in result.stdout.splitlines() if line]

    def test_both_reconstructed_entries_propagate_default_budget(self):
        commands = self.dry_commands("run")
        compute = [command for command in commands if len(command) > 1 and command[1]
                   in ("scripts/prepare_entry.py", "scripts/run_blast.py", "scripts/infer_entry.py")
                   and "--check-tools" not in command]
        self.assertEqual(len(compute), 6)
        self.assertTrue(all("--threads=56" in command for command in compute))

    def test_quoted_parameters_survive_make_template_expansion(self):
        name = "target's transporter"
        directory = "data/input's genomes"
        commands = self.dry_commands("plot-auto", "TARGET_GENE=" + name, "GENOME_DIR=" + directory,
                                     "RSCRIPT=/opt/R env/bin/Rscript", "GENE_PALETTE=science",
                                     "LINK_COLORS=#123456,#abcdef,#fedcba")
        prepared = next(command for command in commands if "scripts/prepare_entry.py" in command)
        plotted = next(command for command in commands if "scripts/plot_synteny.R" in command)
        self.assertIn("--target-gene=" + name, prepared)
        self.assertIn("--genome-dir=" + directory, prepared)
        self.assertEqual(plotted[0], "/opt/R env/bin/Rscript")
        self.assertIn("--target-gene=" + name, plotted)
        self.assertIn("--gene-palette=science", plotted)
        self.assertIn("--link-colors=#123456,#abcdef,#fedcba", plotted)

    def test_setup_passes_explicit_tools(self):
        commands = self.dry_commands("setup", "BEDTOOLS=/opt/tool space/bedtools",
                                     "BLASTN=/opt/blastn", "MAKEBLASTDB=/opt/makeblastdb",
                                     "MAFFT=/opt/mafft", "IQTREE=/opt/iqtree3")
        self.assertIn("--bedtools=/opt/tool space/bedtools", commands[0])
        self.assertIn("--blastn=/opt/blastn", commands[1])
        self.assertIn("--makeblastdb=/opt/makeblastdb", commands[1])
        self.assertIn("--mafft=/opt/mafft", commands[2])
        self.assertIn("--iqtree=/opt/iqtree3", commands[2])

    def test_main_has_no_genome_tools_or_slurm_requirement(self):
        commands = self.dry_commands("plot-main")
        self.assertEqual(len(commands), 2)
        self.assertTrue(all(command[0].endswith("Rscript") for command in commands))


if __name__ == "__main__":
    unittest.main()
