import itertools
import json
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from timestampseeds.cli import parser, run
from timestampseeds.core import CURVES, MAX_KEYS_PER_SEED, MAX_SKIP_DRAWS, attempt_at, candidate, scalar


class SequenceTests(unittest.TestCase):
    def test_native_python_continuing_stream(self):
        rng = random.Random(42)
        for _ in range(7):
            rng.getrandbits(32)
        order = CURVES["p256"][1]
        for index in range(5):
            value = rng.getrandbits(256)
            self.assertEqual(scalar("python-mt19937", 42, "python-getrandbits-reject", order, 7, index),
                             value if 1 <= value < order else None)

    def test_estimate_order_and_dedup(self):
        result = run(parser().parse_args(["estimate", "--start", "1", "--end", "3", "--rng", "lrand48", "msvc-rand",
                                         "--skip-draws", "0", "4", "4", "--keys-per-seed", "3"]))
        self.assertEqual(result["attempts"], 24)
        expected = [(t, r, "concat-be-reject", "p256", s, k)
                    for t, r, s, k in itertools.product([1, 2], ["lrand48", "msvc-rand"], [0, 4], range(3))]
        self.assertEqual([attempt_at(result["config"], i) for i in range(24)], expected)

    def test_invalid_bounds(self):
        for option, value in [("--skip-draws", -1), ("--skip-draws", MAX_SKIP_DRAWS + 1),
                              ("--keys-per-seed", 0), ("--keys-per-seed", MAX_KEYS_PER_SEED + 1)]:
            with self.assertRaises(ValueError):
                run(parser().parse_args(["estimate", option, str(value)]))

    def test_resume_mid_sequence_matches_direct_candidates(self):
        with tempfile.TemporaryDirectory() as directory:
            out, cp = (Path(directory) / name for name in ("out", "cp"))
            run(parser().parse_args(["generate", "--start", "42", "--end", "43", "--rng", "drand48-state",
                                    "--skip-draws", "0", "7", "--keys-per-seed", "4", "--max-candidates", "2",
                                    "--output", str(out), "--checkpoint", str(cp)]))
            run(parser().parse_args(["resume", str(cp)]))
            expected = [candidate(42, "seconds", "drand48-state", "concat-be-reject", "p256", skip, index)
                        for skip, index in itertools.product([0, 7], range(4))]
            self.assertEqual([json.loads(line) for line in out.read_text().splitlines()], [x for x in expected if x])

    def test_rejected_blocks_advance_sequence(self):
        # Blocks encode 0 (reject), 1 (accept), 2 (reject at order). Reusing or
        # retrying a block would change the accepted record's key_index.
        blocks = iter([(0, 4)] * 8 + [(0, 4)] * 7 + [(1, 4)] + [(0, 4)] * 7 + [(2, 4)])
        with (tempfile.TemporaryDirectory() as directory,
              patch.dict(CURVES, {"p256": (CURVES["p256"][0], 2)}),
              patch("timestampseeds.cli.stream_at", return_value=blocks) as initialize):
            out = Path(directory) / "out"
            result = run(parser().parse_args(["generate", "--start", "1", "--end", "2", "--keys-per-seed", "3", "--output", str(out)]))
            self.assertEqual(result["attempted_this_run"], 3)
            self.assertEqual(result["rejected_this_run"], 2)
            self.assertEqual(result["emitted_total"], 1)
            self.assertEqual(json.loads(out.read_text())["key_index"], 1)
            initialize.assert_called_once()

    def test_schema_one_resume_keeps_default_output_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            out, cp, full = (Path(directory) / name for name in ("out", "cp", "full"))
            common = ["generate", "--start", "1", "--end", "4"]
            run(parser().parse_args(common + ["--output", str(out), "--checkpoint", str(cp), "--max-candidates", "1"]))
            data = json.loads(cp.read_text())
            data["schema_version"] = 1
            del data["config"]["skip_draws"], data["config"]["keys_per_seed"]
            cp.write_text(json.dumps(data))
            run(parser().parse_args(["resume", str(cp)]))
            run(parser().parse_args(common + ["--output", str(full)]))
            self.assertEqual(out.read_bytes(), full.read_bytes())
            self.assertEqual(json.loads(cp.read_text())["schema_version"], 2)
