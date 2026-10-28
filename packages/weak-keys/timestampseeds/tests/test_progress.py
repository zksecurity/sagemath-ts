import io
import json
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest.mock import patch

from timestampseeds.cli import parser, run


class ProgressTests(unittest.TestCase):
    def test_progress_and_summary_count_all_attempts(self):
        with tempfile.TemporaryDirectory() as directory:
            out, match = (Path(directory) / name for name in ("out", "match"))
            match.write_text("a" * 64)
            errors = io.StringIO()
            records = [None, {"spki_sha256": "b" * 64}, {"spki_sha256": "a" * 64}]
            with (redirect_stderr(errors), patch("timestampseeds.cli.candidate", side_effect=records),
                  patch("timestampseeds.cli.time.monotonic", side_effect=[0, 1, 2, 3, 4])):
                result = run(parser().parse_args(["generate", "--start", "1", "--end", "10", "--max-candidates", "3",
                                                 "--match", str(match), "--output", str(out), "--progress-every", "1"]))
            progress = [json.loads(line) for line in errors.getvalue().splitlines()]
            self.assertEqual([x["attempted_this_run"] for x in progress], [1, 2, 3])
            self.assertEqual([x["remaining_this_run"] for x in progress], [2, 1, 0])
            self.assertEqual(progress[0]["eta_seconds_this_run"], 2)
            self.assertTrue(all(x["event"] == "progress" for x in progress))
            self.assertEqual(result["rejected_this_run"], 1)
            self.assertEqual(result["nonmatching_this_run"], 1)
            self.assertEqual(result["emitted_total"], 1)
            self.assertEqual(result["attempts_per_second"], 0.75)
            self.assertFalse(result["complete"])
            self.assertEqual(json.loads(out.read_text()), records[-1])

    def test_disabled_progress_and_invalid_intervals(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / "out"
            errors = io.StringIO()
            with redirect_stderr(errors):
                run(parser().parse_args(["generate", "--start", "1", "--end", "2", "--output", str(out), "--progress-every", "0"]))
            self.assertEqual(errors.getvalue(), "")
            for value in ["nan", "inf", "-1"]:
                with self.assertRaisesRegex(ValueError, "Progress interval"):
                    run(parser().parse_args(["generate", "--progress-every", value]))
