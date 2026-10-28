import random
import unittest
from unittest.mock import patch

from timestampseeds.core import CURVES, scalar, scalar_from_words
from timestampseeds.cli import parser, run


class ProfileTests(unittest.TestCase):
    def test_native_python_apis(self):
        order = CURVES["p256"][1]
        for seed in [0, 1, 1700000000, 2**32 + 1, 2**64 + 42]:
            for profile, expected in [
                ("python-getrandbits-reject", random.Random(seed).getrandbits(256)),
                ("python-randbytes-reject", int.from_bytes(random.Random(seed).randbytes(32), "big")),
            ]:
                with self.subTest(seed=seed, profile=profile):
                    self.assertEqual(scalar("python-mt19937", seed, profile, order),
                                     expected if 1 <= expected < order else None)

    def test_word_order_and_scalar_order_are_distinct(self):
        material = bytes(range(1, 33))
        stream = [(int.from_bytes(material[i:i+4], "big"), 4) for i in range(0, 32, 4)]
        order = 2**256
        self.assertEqual(scalar_from_words(iter(stream), "buffer-le-reject", order), int.from_bytes(material, "little"))
        self.assertNotEqual(scalar_from_words(iter(stream), "concat-le-reject", order), int.from_bytes(material, "little"))

    def test_low_byte_consumes_one_word_per_byte(self):
        stream = iter((0x123400 + i, 4) for i in range(33))
        self.assertEqual(scalar_from_words(stream, "low-byte-reject", 2**256), int.from_bytes(bytes(range(32)), "big"))
        self.assertEqual(next(stream), (0x123420, 4))

    def test_six_byte_words_discard_unused_tail(self):
        stream = [(int.from_bytes(bytes(range(i, i+6)), "big"), 6) for i in range(1, 43, 6)]
        for profile, expected in [
            ("concat-be-reject", bytes(range(1, 33))),
            ("concat-le-reject", b"".join(bytes(range(i, i+6))[::-1] for i in range(1, 37, 6))[:32]),
        ]:
            iterator = iter(stream)
            self.assertEqual(scalar_from_words(iterator, profile, 2**256), int.from_bytes(expected, "big"))
            self.assertEqual(next(iterator), stream[6])

    def test_rejection_boundaries_and_modulo(self):
        for _, order in CURVES.values():
            for value in [0, 1, order - 1, order, order + 1, 2**256 - 1]:
                material = value.to_bytes(32, "big")
                stream = [(int.from_bytes(material[i:i+4], "big"), 4) for i in range(0, 32, 4)]
                with patch("timestampseeds.core.words", return_value=iter(stream)):
                    self.assertEqual(scalar("lrand48", 1, "concat-be-reject", order),
                                     value if 1 <= value < order else None)
                self.assertEqual(scalar_from_words(iter(stream), "concat-be-mod", order), value % (order - 1) + 1)

    def test_native_profiles_reject_incompatible_rng(self):
        for profile in ["python-getrandbits-reject", "python-randbytes-reject", "java-nextbytes-reject"]:
            with self.assertRaisesRegex(ValueError, "requires --rng"):
                run(parser().parse_args(["estimate", "--profile", profile]))
