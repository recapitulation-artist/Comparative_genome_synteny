"""CPU allocation checks without launching analysis at alternate thread counts."""

import os
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import runtime


class RuntimeTests(unittest.TestCase):
    def test_local_default_56_has_no_scheduler_requirement(self):
        with patch.dict(os.environ, {}, clear=True), patch.object(runtime, "available_cpus", return_value=64):
            self.assertEqual(runtime.validate_threads(runtime.DEFAULT_THREADS), 56)
            self.assertEqual(runtime.execution_context(56)["mode"], "local")

    def test_insufficient_affinity_does_not_silently_reduce_request(self):
        with patch.dict(os.environ, {}, clear=True), patch.object(runtime, "available_cpus", return_value=32):
            with self.assertRaisesRegex(ValueError, "THREADS=56.*32 CPUs"):
                runtime.validate_threads(56)

    def test_slurm_allocation_can_exceed_request_without_specific_node(self):
        env = dict(SLURM_JOB_ID="123", SLURM_CPUS_PER_TASK="64", SLURM_NTASKS="1")
        with patch.dict(os.environ, env, clear=True), patch.object(runtime, "available_cpus", return_value=64):
            self.assertEqual(runtime.validate_threads(56), 56)
            context = runtime.execution_context(56)
            self.assertEqual(context["mode"], "slurm")
            self.assertEqual(context["cpus_per_task"], 64)

    def test_slurm_cpu_shortage_and_multiple_tasks_are_rejected(self):
        with patch.object(runtime, "available_cpus", return_value=64):
            with patch.dict(os.environ, dict(SLURM_JOB_ID="123", SLURM_CPUS_PER_TASK="32"), clear=True):
                with self.assertRaisesRegex(ValueError, "at least 56 CPUs"):
                    runtime.validate_threads(56)
            with patch.dict(os.environ, dict(SLURM_JOB_ID="123", SLURM_CPUS_PER_TASK="56", SLURM_NTASKS="2"), clear=True):
                with self.assertRaisesRegex(ValueError, "one Slurm task"):
                    runtime.validate_threads(56)

    def test_native_budget_does_not_multiply_nested_library_threads(self):
        env = runtime.thread_environment(56)
        self.assertEqual(env["OMP_NUM_THREADS"], "56")
        self.assertEqual(env["OMP_DYNAMIC"], "FALSE")
        for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            self.assertEqual(env[key], "1")


if __name__ == "__main__":
    unittest.main()
