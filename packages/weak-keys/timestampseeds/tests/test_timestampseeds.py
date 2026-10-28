import argparse
import base64
import ctypes
import hashlib
import itertools
import json
import platform
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

from cryptography.hazmat.primitives import serialization

from timestampseeds.cli import parser, run, targets
from timestampseeds.core import CURVES, candidate, combination_at, parse_time, scalar
from timestampseeds.rngs import words


class RngTests(unittest.TestCase):
    def test_glibc_seed_one_vector(self):
        self.assertEqual([x for x, _ in itertools.islice(words("glibc-rand", 1), 5)],
                         [1804289383, 846930886, 1681692777, 1714636915, 1957747793])

    def test_ansi_explicit_vector(self):
        self.assertEqual([x for x, _ in itertools.islice(words("ansi-lcg31", 1), 3)],
                         [1103527590, 377401575, 662824084])

    def test_python_not_cpp_mt_seed(self):
        self.assertEqual(next(words("python-mt19937", 1)), (577090037, 4))

    @unittest.skipUnless(platform.libc_ver()[0] == "glibc", "glibc reference only")
    def test_linux_libc_vectors(self):
        libc = ctypes.CDLL(None)
        libc.srand48.argtypes = [ctypes.c_long]
        libc.lrand48.restype = ctypes.c_long
        libc.drand48.restype = ctypes.c_double
        for seed in [0, 1, 42, 1700000000, 0x80000000, 0xFFFFFFFF, 0x100000001]:
            libc.srand48(seed)
            self.assertEqual([libc.lrand48() for _ in range(20)],
                             [v for v, _ in itertools.islice(words("lrand48", seed), 20)])
            libc.srand48(seed)
            self.assertEqual([libc.drand48() for _ in range(20)],
                             [v / (1 << 48) for v, _ in itertools.islice(words("drand48-state", seed), 20)])
            libc.srand(ctypes.c_uint(seed))
            self.assertEqual([libc.rand() for _ in range(20)],
                             [v for v, _ in itertools.islice(words("glibc-rand", seed), 20)])

    def test_srand48_seed_truncation(self):
        self.assertEqual(list(itertools.islice(words("lrand48", 1), 8)),
                         list(itertools.islice(words("lrand48", 2**32 + 1), 8)))


class CoreTests(unittest.TestCase):
    def test_utc_units_and_timezones(self):
        self.assertEqual(parse_time("1970-01-01T01:00:01+01:00", "seconds"), 1)
        self.assertEqual(parse_time("1970-01-01T00:00:00.001Z", "milliseconds"), 1)
        self.assertEqual(parse_time("1700000000", "seconds"), 1700000000)
        for stamp in ["2020-01-01T00:00:00", "1970-01-01T00:00:00.001Z"]:
            with self.assertRaises(ValueError):
                parse_time(stamp, "seconds")

    def test_serialization_profile(self):
        order = CURVES["p256"][1]
        expected = b"".join(v.to_bytes(w, "big") for v, w in itertools.islice(words("lrand48", 42), 8))
        self.assertEqual(scalar("lrand48", 42, "concat-be-reject", order), int.from_bytes(expected, "big"))
        little = b"".join(v.to_bytes(w, "little") for v, w in itertools.islice(words("lrand48", 42), 8))
        self.assertEqual(scalar("lrand48", 42, "concat-le-reject", order),
                         int.from_bytes(little, "big") if 1 <= int.from_bytes(little, "big") < order else None)
        self.assertEqual(scalar("lrand48", 42, "concat-be-mod", 7), int.from_bytes(expected, "big") % 6 + 1)
        self.assertIsNone(scalar("lrand48", 42, "concat-be-reject", 7))

    def test_public_only_and_reparse(self):
        for curve in CURVES:
            record = candidate(1700000000, "seconds", "lrand48", "concat-be-reject", curve)
            der = base64.b64decode(record["spki_der_base64"])
            key = serialization.load_der_public_key(der)
            self.assertEqual(hashlib.sha256(der).hexdigest(), record["spki_sha256"])
            self.assertEqual(len(bytes.fromhex(record["sec1_compressed_hex"])), 33)
            self.assertEqual(set(record), {"timestamp", "unit", "rng", "profile", "curve", "spki_sha256", "spki_der_base64", "sec1_compressed_hex"})
            self.assertEqual(key.key_size, 256)


class CliTests(unittest.TestCase):
    def parse(self, *args):
        return parser().parse_args(args)

    def test_estimate_end_exclusive_step_and_dedup(self):
        result = run(self.parse("estimate", "--start", "1", "--end", "6", "--step", "2", "--rng", "lrand48", "lrand48", "--curve", "p256", "secp256k1"))
        self.assertEqual(result["attempts"], 6)
        self.assertEqual(result["timestamps"], 3)
        self.assertEqual(combination_at(result["config"], 5), (5, "lrand48", "concat-be-reject", "secp256k1"))

    def test_invalid_ranges(self):
        for extra in [("--start", "-1"), ("--step", "0"), ("--start", "5", "--end", "4")]:
            with self.assertRaises(ValueError):
                run(self.parse("estimate", *extra))

    def test_checkpoint_resume_exact_bytes_and_crash_tail(self):
        with tempfile.TemporaryDirectory() as d:
            out, cp, full = [str(Path(d) / n) for n in ("keys.jsonl", "checkpoint.json", "full.jsonl")]
            first = run(self.parse("generate", "--start", "1", "--end", "9", "--max-candidates", "3", "--output", out, "--checkpoint", cp, "--checkpoint-every", "2"))
            self.assertEqual(first["next_index"], 3)
            with open(out, "ab") as f:
                f.write(b"uncommitted crash tail\n")
            second = run(self.parse("resume", cp, "--max-candidates", "2"))
            self.assertEqual(second["next_index"], 5)
            final = run(self.parse("resume", cp))
            self.assertTrue(final["complete"])
            run(self.parse("generate", "--start", "1", "--end", "9", "--output", full))
            self.assertEqual(Path(out).read_bytes(), Path(full).read_bytes())
            self.assertEqual(len(Path(out).read_text().splitlines()), 8)

    def test_checkpoint_detects_modified_prefix(self):
        with tempfile.TemporaryDirectory() as d:
            out, cp = str(Path(d) / "out"), str(Path(d) / "cp")
            run(self.parse("generate", "--start", "1", "--end", "3", "--output", out, "--checkpoint", cp, "--max-candidates", "1"))
            data = Path(out).read_bytes()
            Path(out).write_bytes(b"X" + data[1:])
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                run(self.parse("resume", cp))

    def test_offline_fingerprint_match(self):
        with tempfile.TemporaryDirectory() as d:
            match, out = Path(d) / "match", Path(d) / "out"
            desired = candidate(2, "seconds", "lrand48", "concat-be-reject", "p256")["spki_sha256"]
            match.write_text(desired + "\n")
            result = run(self.parse("generate", "--start", "1", "--end", "5", "--match", str(match), "--output", str(out)))
            self.assertEqual(result["emitted_total"], 1)
            self.assertEqual(json.loads(out.read_text())["timestamp"], 2)
            match.write_text("badfingerprint")
            with self.assertRaises(ValueError):
                targets(str(match))

    def test_interrupt_between_write_and_digest_retains_last_checkpoint(self):
        import timestampseeds.cli as cli
        real_sha256 = hashlib.sha256
        class InterruptingDigest:
            def __init__(self):
                self.inner = real_sha256()
                self.calls = 0
            def update(self, block):
                self.calls += 1
                if self.calls == 2:
                    raise KeyboardInterrupt()
                self.inner.update(block)
            def hexdigest(self):
                return self.inner.hexdigest()
        with tempfile.TemporaryDirectory() as d:
            out, cp = str(Path(d) / "out"), str(Path(d) / "cp")
            with patch.object(cli, "candidate", side_effect=lambda stamp, *a: {"timestamp": stamp}), patch.object(cli.hashlib, "sha256", side_effect=InterruptingDigest):
                with self.assertRaises(KeyboardInterrupt):
                    run(self.parse("generate", "--start", "1", "--end", "4", "--output", out, "--checkpoint", cp, "--checkpoint-every", "1"))
            saved = json.loads(Path(cp).read_text())
            self.assertEqual(saved["next_index"], 1)
            self.assertGreater(Path(out).stat().st_size, saved["output_offset"])
            with patch.object(cli, "candidate", side_effect=lambda stamp, *a: {"timestamp": stamp}):
                run(self.parse("resume", cp))
            self.assertEqual([json.loads(x)["timestamp"] for x in Path(out).read_text().splitlines()], [1, 2, 3])

    def test_checkpoint_changed_before_lock_does_not_truncate(self):
        import timestampseeds.cli as cli
        with tempfile.TemporaryDirectory() as d:
            out, cp = str(Path(d) / "out"), str(Path(d) / "cp")
            run(self.parse("generate", "--start", "1", "--end", "4", "--output", out, "--checkpoint", cp, "--max-candidates", "1"))
            original = Path(out).read_bytes()
            real_flock = cli.fcntl.flock
            def advance_before_lock(stream, flags):
                if isinstance(stream.name, str) and Path(stream.name).resolve() == Path(out).resolve():
                    advanced = json.loads(Path(cp).read_text())
                    advanced["emitted"] += 1
                    Path(cp).write_text(json.dumps(advanced))
                return real_flock(stream, flags)
            with patch.object(cli.fcntl, "flock", side_effect=advance_before_lock):
                with self.assertRaisesRegex(ValueError, "Checkpoint changed"):
                    run(self.parse("resume", cp))
            self.assertEqual(Path(out).read_bytes(), original)

    def test_existing_file_refused(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "out"
            out.write_text("do not overwrite")
            with self.assertRaises(FileExistsError):
                run(self.parse("generate", "--start", "1", "--end", "2", "--output", str(out)))
            self.assertEqual(out.read_text(), "do not overwrite")


if __name__ == "__main__":
    unittest.main()
