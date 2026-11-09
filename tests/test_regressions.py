import json
import io
from contextlib import redirect_stdout
from pathlib import Path
import tempfile
import unittest

import numpy as np

from polypangf.core import InputError, Runner, project_lock
from polypangf.evolution import merge_intervals, covered
from polypangf.family import parse_groups
from polypangf.inputs import fasta
from polypangf.pan import bh_adjust, classify
from polypangf.pan import pan
from polypangf.core import write_tsv
from polypangf.core import read_tsv
from polypangf.cli import main


class StatisticalRegressions(unittest.TestCase):
    def test_bh_preserves_zero_hit_test_count(self):
        np.testing.assert_allclose(bh_adjust([0.01, 0.04, 1, 1]), [0.04, 0.08, 1, 1])
        self.assertEqual(len(bh_adjust([])), 0)

    def test_coverage_edge_cases(self):
        self.assertEqual(classify(13, 13), "Core")
        self.assertEqual(classify(12, 13), "Softcore")
        self.assertEqual(classify(11, 13), "Shell")
        self.assertEqual(classify(1, 13), "Specific")
        self.assertEqual(classify(1, 2), "Specific")

    def test_nested_te_intervals_count_once(self):
        intervals = merge_intervals([(0, 10), (2, 4), (8, 12), (12, 15)])
        self.assertEqual(intervals, [(0, 15)])
        self.assertEqual(covered(intervals, 5, 20), 10)
        with self.assertRaises(InputError):
            merge_intervals([(5, 5)])


class InputRegressions(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_duplicate_fasta_and_internal_stop_rejected(self):
        path = self.root / "bad.fa"
        path.write_text(">a\nMAA\n>a\nMAA\n")
        with self.assertRaises(InputError):
            fasta(path, protein=True)

    def test_masking_is_explicit_and_audited(self):
        path = self.root / "proteins.fa"
        path.write_text(">a\nMA*A\n>b\nMAGA\n")
        audit = []
        self.assertEqual(
            fasta(path, protein=True, invalid_policy="mask", audit=audit)["a"], "MAXA"
        )
        self.assertEqual(audit[0]["action"], "mask")
        self.assertEqual(set(fasta(path, protein=True, invalid_policy="skip")), {"b"})
        path.write_text(">a\nMA*A\n")
        with self.assertRaises(InputError):
            fasta(path, protein=True)
        path.write_text(">a\nMA*A\n>a\nMAGA\n")
        with self.assertRaises(InputError):
            fasta(path, protein=True, invalid_policy="skip")

    def test_duplicate_tsv_columns_are_rejected(self):
        path = self.root / "bad.tsv"
        path.write_text("gene_id\tgene_id\na\tb\n")
        with self.assertRaises(InputError):
            read_tsv(path)

    def test_dataset_scoped_original_ids(self):
        rows = [
            {"uid": "A__1", "dataset": "A", "gene_id": "g1", "protein_id": "p1"},
            {"uid": "B__1", "dataset": "B", "gene_id": "g1", "protein_id": "p1"},
        ]
        path = self.root / "groups.tsv"
        path.write_text("Orthogroup\tA\tB\nOG1\tp1\tp1\n")
        groups, membership, _ = parse_groups(path, rows)
        self.assertEqual(set(groups["OG1"]), {"A__1", "B__1"})
        path.write_text("Orthogroup\tA\tB\nOG1\tp1\tp1\nOG2\tp1\t\n")
        with self.assertRaises(InputError):
            parse_groups(path, rows)

    def test_unknown_members_require_explicit_scope(self):
        rows = [{"uid": "A__1", "dataset": "A", "gene_id": "g1", "protein_id": "p1"}]
        path = self.root / "groups.tsv"
        path.write_text("Orthogroup\tA\nOG1\tnonfamily\n")
        with self.assertRaises(InputError):
            parse_groups(path, rows)
        self.assertEqual(parse_groups(path, rows, allow_nonfamily=True)[2], 1)
        path.write_text("Orthogroup\tTYPO\nOG1\tnonfamily\n")
        with self.assertRaises(InputError):
            parse_groups(path, rows, allow_nonfamily=True)


class TransactionRegressions(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.input = self.root / "input.txt"
        self.input.write_text("input")
        self.cfg = {"output": str(self.root / "results")}
        Path(self.cfg["output"]).mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def build(self, work):
        (work / "value.txt").write_text("complete")
        return {"ok": True}

    def test_resume_verifies_outputs_and_inputs(self):
        runner = Runner(self.cfg, resume=True)
        dest = runner.run("stage", self.build, [self.input], {})
        state = (dest / "state.json").read_text()
        runner.run(
            "stage", lambda _: self.fail("cache was not reused"), [self.input], {}
        )
        self.assertEqual((dest / "state.json").read_text(), state)
        (dest / "value.txt").write_text("damaged")
        runner.run("stage", self.build, [self.input], {})
        self.assertEqual((dest / "value.txt").read_text(), "complete")
        self.input.write_text("changed")
        runner.run("stage", self.build, [self.input], {})
        self.assertNotEqual((dest / "state.json").read_text(), state)

    def test_failed_recompute_preserves_complete_outputs(self):
        runner = Runner(self.cfg, resume=True)
        dest = runner.run("stage", self.build, [self.input], {})

        def fail(work):
            (work / "partial.txt").write_text("incomplete")
            raise InputError("expected failure")

        self.input.write_text("changed")
        with self.assertRaises(InputError):
            runner.run("stage", fail, [self.input], {})
        self.assertEqual((dest / "value.txt").read_text(), "complete")
        self.assertTrue((self.root / "results/failures/stage.json").exists())

    def test_corrupt_state_can_be_rebuilt(self):
        runner = Runner(self.cfg, resume=True)
        dest = runner.run("stage", self.build, [self.input], {})
        (dest / "state.json").write_text("{invalid JSON")
        runner.run("stage", self.build, [self.input], {})
        self.assertEqual(
            json.loads((dest / "state.json").read_text())["status"], "complete"
        )

    def test_input_change_during_run_is_not_published(self):
        runner = Runner(self.cfg, resume=True)
        dest = runner.run("stage", self.build, [self.input], {})
        self.input.write_text("new input")

        def mutate(work):
            self.input.write_text("changed while running")
            return self.build(work)

        with self.assertRaises(InputError):
            runner.run("stage", mutate, [self.input], {})
        self.assertEqual((dest / "value.txt").read_text(), "complete")

    def test_concurrent_writes_are_rejected(self):
        with project_lock(self.cfg["output"]):
            with self.assertRaises(InputError):
                with project_lock(self.cfg["output"]):
                    pass


class WorkflowRegressions(unittest.TestCase):
    def test_single_dataset_identification_needs_no_orthology_config(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "proteins.fa").write_text(">short\nMAGA\n>long\nMAGAAAAAAA\n")
            write_tsv(
                root / "models.tsv",
                [
                    {"protein_id": "short", "gene_id": "g1"},
                    {"protein_id": "long", "gene_id": "g1"},
                ],
                ["protein_id", "gene_id"],
            )
            write_tsv(
                root / "datasets.tsv",
                [
                    {
                        "dataset": "ONE",
                        "species": "one",
                        "ploidy": 2,
                        "represented_sets": 1,
                        "proteins": "proteins.fa",
                        "model_map": "models.tsv",
                    }
                ],
                [
                    "dataset",
                    "species",
                    "ploidy",
                    "represented_sets",
                    "proteins",
                    "model_map",
                ],
            )
            write_tsv(
                root / "candidates.tsv",
                [{"dataset": "ONE", "protein_id": "short"}],
                ["dataset", "protein_id"],
            )
            (root / "config.json").write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "manifest": "datasets.tsv",
                        "output": "results",
                        "identification": {"candidates": "candidates.tsv"},
                    }
                )
            )
            with redirect_stdout(io.StringIO()):
                self.assertEqual(main(["identify", "-c", str(root / "config.json")]), 0)
            self.assertEqual(
                read_tsv(root / "results/identify/candidates.tsv")[0]["protein_id"],
                "short",
            )

    def test_same_species_assemblies_do_not_add_presence_votes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "results"
            (output / "identify").mkdir(parents=True)
            (output / "orthology").mkdir()
            dest = output / "pan"
            dest.mkdir()
            (root / "proteins.fa").write_text(">p\nMAGA\n")
            samples = [
                {
                    "dataset": dataset,
                    "species": species,
                    "ploidy": 2,
                    "represented_sets": 1,
                    "proteins": "proteins.fa",
                }
                for dataset, species in [("A1", "A"), ("A2", "A"), ("B", "B")]
            ]
            write_tsv(root / "datasets.tsv", samples, list(samples[0]))
            rows = [
                {
                    "uid": dataset + "__1",
                    "dataset": dataset,
                    "species": "A",
                    "gene_id": "g1",
                    "role": "analysis",
                    "ploidy": "2",
                    "represented_sets": "1",
                }
                for dataset in ["A1", "A2"]
            ]
            write_tsv(output / "identify/candidates.tsv", rows, list(rows[0]))
            write_tsv(
                output / "orthology/membership.tsv",
                [{"uid": row["uid"], "group": "OG1"} for row in rows],
                ["uid", "group"],
            )
            summary = pan(
                {"manifest": str(root / "datasets.tsv"), "output": str(output)}, dest
            )
            self.assertEqual(summary["units"], 2)
            self.assertEqual(summary["classes"], {"Specific": 1})
            self.assertEqual(int(read_tsv(dest / "groups.tsv")[0]["present_units"]), 1)


if __name__ == "__main__":
    unittest.main()
