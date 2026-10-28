import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import timestampseeds.cli as cli
from timestampseeds.core import parse_time


def run(*args):
    return cli.run(cli.parser().parse_args(args))


class IntegrityTests(unittest.TestCase):
    def test_exact_timestamp_fractions(self):
        for stamp, unit in [
            ("1970-01-01T00:00:01.0000001Z", "seconds"),
            ("1970-01-01T00:00:00.0010001Z", "milliseconds"),
            ("1970-01-01T00:00:01,0000001Z", "seconds"),
            ("1970-01-01T00:00:01+00:00:00.0000001", "seconds"),
            ("1970-01-01T00:00:01+00:00:00.5", "seconds"),
        ]:
            with self.subTest(stamp=stamp), self.assertRaises(ValueError):
                parse_time(stamp, unit)
        self.assertEqual(parse_time("1970-01-01T00:00:01.0000000Z", "seconds"), 1)
        self.assertEqual(parse_time("1970-01-01T00:00:00.001000000Z", "milliseconds"), 1)

    def test_invalid_checkpoint_never_changes_output(self):
        with tempfile.TemporaryDirectory() as directory:
            out, cp = (Path(directory) / name for name in ("out", "cp"))
            run("generate", "--start", "1", "--end", "4", "--max-candidates", "1",
                "--output", str(out), "--checkpoint", str(cp))
            saved = json.loads(cp.read_text())
            with out.open("ab") as stream:
                stream.write(b"uncommitted tail\n")
            original = out.read_bytes()
            bad_values = {"emitted": ["bad", True, -1, 2], "next_index": [True, -1, 4],
                          "output_offset": [True, -1], "output_sha256": ["bad", None],
                          "output": ["relative", None], "config": [None, []]}
            for field, values in bad_values.items():
                for value in values:
                    data = copy.deepcopy(saved)
                    data[field] = value
                    cp.write_text(json.dumps(data))
                    with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                        run("resume", str(cp))
                    self.assertEqual(out.read_bytes(), original)
            for value in ["bad", {}, ["not-a-hash"]]:
                data = copy.deepcopy(saved)
                data["config"]["match_fingerprints"] = value
                cp.write_text(json.dumps(data))
                with self.assertRaises(ValueError):
                    run("resume", str(cp))
                self.assertEqual(out.read_bytes(), original)
            for value in [[], None, {}, {"schema_version": True}]:
                cp.write_text(json.dumps(value))
                with self.assertRaises(ValueError):
                    run("resume", str(cp))
                self.assertEqual(out.read_bytes(), original)

    def test_shared_checkpoint_is_locked_before_first_save(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second, cp = (Path(directory) / name for name in ("first", "second", "cp"))
            real_save = cli.atomic_json
            def competing_writer(path, data):
                with self.assertRaisesRegex(ValueError, "Checkpoint is in use"):
                    run("generate", "--start", "20", "--end", "22",
                        "--output", str(second), "--checkpoint", str(cp))
                self.assertFalse(second.exists())
                real_save(path, data)
            with patch.object(cli, "atomic_json", side_effect=competing_writer):
                run("generate", "--start", "1", "--end", "3",
                    "--output", str(first), "--checkpoint", str(cp))
            self.assertEqual(json.loads(cp.read_text())["output"], str(first.resolve()))
            self.assertTrue(run("resume", str(cp))["complete"])

    def test_output_cannot_alias_checkpoint_lock(self):
        with tempfile.TemporaryDirectory() as directory:
            cp = Path(directory) / "cp"
            with self.assertRaisesRegex(ValueError, "must differ"):
                run("generate", "--start", "1", "--end", "2",
                    "--output", str(cp) + ".lock", "--checkpoint", str(cp))
