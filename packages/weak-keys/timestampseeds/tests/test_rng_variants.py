import itertools
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from timestampseeds.rngs import words


REFERENCE = Path(__file__).parent / "reference"
SEEDS = [0, 1, 42, 1700000000, 0x80000000, 0xFFFFFFFF, 0x100000001, 2**48 + 42]


class VariantTests(unittest.TestCase):
    def test_known_vectors(self):
        # Java Random(0); Microsoft srand(1); standard MT19937 default seed 5489.
        for rng, seed, expected in [
            ("java-random", 0, [3139482720, 3571011896, 1033096058, 2604232894, 2737687030]),
            ("msvc-rand", 1, [41, 18467, 6334, 26500, 19169]),
            ("cpp-mt19937", 5489, [3499211612, 581869302, 3890346734, 3586334585, 545404204]),
        ]:
            with self.subTest(rng=rng):
                self.assertEqual([v for v, _ in itertools.islice(words(rng, seed), 5)], expected)

    def test_seed_widths(self):
        for rng, bits in [("java-random", 48), ("msvc-rand", 32), ("cpp-mt19937", 32)]:
            self.assertEqual(list(itertools.islice(words(rng, 42), 20)),
                             list(itertools.islice(words(rng, 42 + 2**bits), 20)))

    @unittest.skipUnless(shutil.which("c++"), "C++ reference compiler unavailable")
    def test_cpp_native_reference_across_twists(self):
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / "mt-reference"
            subprocess.run(["c++", "-std=c++11", str(REFERENCE / "mt19937.cpp"), "-o", str(executable)], check=True)
            for seed in SEEDS:
                expected = subprocess.check_output([str(executable), str(seed), "1300"], text=True)
                self.assertEqual([v for v, _ in itertools.islice(words("cpp-mt19937", seed), 1300)],
                                 [int(line) for line in expected.splitlines()])

    @unittest.skipUnless(shutil.which("javac"), "Java reference compiler unavailable")
    def test_java_native_reference(self):
        # macOS may expose a javac launcher with no JDK behind it.
        probe = subprocess.run(["javac", "-version"], capture_output=True)
        if probe.returncode:
            self.skipTest("Java reference compiler unavailable")
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run(["javac", "-d", directory, str(REFERENCE / "RandomReference.java")], check=True)
            for seed in SEEDS:
                expected = subprocess.check_output(["java", "-cp", directory, "RandomReference", str(seed), "100"], text=True)
                self.assertEqual([v for v, _ in itertools.islice(words("java-random", seed), 100)],
                                 [int(line) for line in expected.splitlines()])

