import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from timestampseeds.cli import checkpoint_lock


class CommandLineTests(unittest.TestCase):
    def invoke(self, *args):
        return subprocess.run([sys.executable, "-m", "timestampseeds", *args], capture_output=True, text=True)

    def test_generate_resume_and_progress_are_separate_from_output(self):
        with tempfile.TemporaryDirectory() as directory:
            out, cp, full = (Path(directory) / name for name in ("out", "cp", "full"))
            options = ["--start", "1700000000", "--end", "1700000002", "--rng", "java-random",
                       "--profile", "java-nextbytes-reject", "--curve", "p256", "secp256k1",
                       "--skip-draws", "0", "1", "8", "--keys-per-seed", "4"]
            estimate = self.invoke("estimate", *options)
            self.assertEqual(estimate.returncode, 0, estimate.stderr)
            self.assertEqual(json.loads(estimate.stdout)["attempts"], 48)
            first = self.invoke("generate", *options, "--output", str(out), "--checkpoint", str(cp),
                                "--max-candidates", "5", "--progress-every", "0.000001")
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertEqual(first.stdout, "")
            reports = [json.loads(line) for line in first.stderr.splitlines()]
            self.assertTrue(any(report.get("event") == "progress" for report in reports))
            self.assertEqual(reports[-1]["next_index"], 5)
            resumed = self.invoke("resume", str(cp), "--progress-every", "0")
            self.assertEqual(resumed.returncode, 0, resumed.stderr)
            self.assertTrue(json.loads(resumed.stderr)["complete"])
            complete = self.invoke("generate", *options, "--output", str(full), "--progress-every", "0")
            self.assertEqual(complete.returncode, 0, complete.stderr)
            self.assertEqual(out.read_bytes(), full.read_bytes())
            records = [json.loads(line) for line in out.read_text().splitlines()]
            self.assertEqual(len(records), 48)
            self.assertTrue(all("event" not in record for record in records))

    def test_separate_process_cannot_claim_locked_checkpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            cp, out = (Path(directory) / name for name in ("cp", "out"))
            with checkpoint_lock(cp):
                result = self.invoke("generate", "--start", "1", "--end", "2", "--output", str(out), "--checkpoint", str(cp))
            self.assertEqual(result.returncode, 2)
            self.assertIn("Checkpoint is in use", result.stderr)
            self.assertNotIn("Traceback", result.stderr)
            self.assertFalse(out.exists())
