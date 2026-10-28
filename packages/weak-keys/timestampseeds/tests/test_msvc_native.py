"""Native Windows oracle; the Unix CLI itself is intentionally not imported."""
import ctypes
import itertools
import sys
import unittest

from timestampseeds.rngs import words


@unittest.skipUnless(sys.platform == "win32", "Microsoft UCRT reference requires Windows")
class MicrosoftNativeTests(unittest.TestCase):
    def test_ucrt_srand_rand(self):
        crt = ctypes.CDLL("ucrtbase")
        crt.srand.argtypes = [ctypes.c_uint]
        crt.rand.restype = ctypes.c_int
        for seed in [0, 1, 1792, 1700000000, 0x80000000, 0xFFFFFFFF, 0x100000001]:
            crt.srand(seed)
            self.assertEqual([crt.rand() for _ in range(100)],
                             [v for v, _ in itertools.islice(words("msvc-rand", seed), 100)])
