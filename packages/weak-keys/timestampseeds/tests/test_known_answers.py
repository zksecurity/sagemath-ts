import json
import unittest
from pathlib import Path

from timestampseeds.core import candidate


class KnownAnswerTests(unittest.TestCase):
    def test_independent_seed_to_public_key_fixtures(self):
        fixture = json.loads((Path(__file__).parent / "fixtures" / "public_keys.json").read_text())
        for case in fixture["cases"]:
            with self.subTest(rng=case["rng"], profile=case["profile"], curve=case["curve"]):
                actual = candidate(case["seed"], "seconds", case["rng"], case["profile"], case["curve"],
                                   case["skip_draws"], case["key_index"])
                if case["public"] is None:
                    self.assertIsNone(actual)
                else:
                    self.assertIsNotNone(actual)
                    for field, expected in case["public"].items():
                        self.assertEqual(actual[field], expected)
