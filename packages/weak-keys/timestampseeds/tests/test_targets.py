import base64
import json
import tempfile
import unittest
from datetime import datetime, timedelta, timezone
from pathlib import Path

from cryptography import x509
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import ec, ed25519
from cryptography.x509.oid import NameOID

from timestampseeds.cli import parser, run
from timestampseeds.core import CURVES, candidate, scalar
from timestampseeds.targets import fingerprint, fingerprints


class TargetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.key = ec.derive_private_key(scalar("lrand48", 2, "concat-be-reject", CURVES["p256"][1]), ec.SECP256R1())
        cls.public = cls.key.public_key()
        name = x509.Name([x509.NameAttribute(NameOID.COMMON_NAME, "synthetic test")])
        now = datetime(2024, 1, 1, tzinfo=timezone.utc)
        cls.cert = (x509.CertificateBuilder().subject_name(name).issuer_name(name).public_key(cls.public)
                    .serial_number(1).not_valid_before(now).not_valid_after(now + timedelta(days=1))
                    .sign(cls.key, hashes.SHA256()))
        cls.expected = candidate(2, "seconds", "lrand48", "concat-be-reject", "p256")["spki_sha256"]

    def test_pem_der_keys_and_certificates(self):
        for encoding in [serialization.Encoding.PEM, serialization.Encoding.DER]:
            for data in [self.public.public_bytes(encoding, serialization.PublicFormat.SubjectPublicKeyInfo),
                         self.cert.public_bytes(encoding)]:
                self.assertEqual(fingerprints(data), {self.expected})

    def test_pem_bundle_and_duplicate_fingerprints(self):
        key = self.public.public_bytes(serialization.Encoding.PEM, serialization.PublicFormat.SubjectPublicKeyInfo)
        cert = self.cert.public_bytes(serialization.Encoding.PEM)
        self.assertEqual(fingerprints(key + cert + key), {self.expected})
        self.assertEqual(fingerprints(("# comment\n\n" + self.expected.upper() + "\n" + self.expected).encode()), {self.expected})
        self.assertEqual(fingerprints(b"# no targets\n"), set())

    def test_compressed_spki_normalization(self):
        point = self.public.public_bytes(serialization.Encoding.X962, serialization.PublicFormat.CompressedPoint)
        compressed_der = bytes.fromhex("3039301306072a8648ce3d020106082a8648ce3d030107032200") + point
        self.assertEqual(fingerprints(compressed_der), {self.expected})

    def test_reject_private_unsupported_and_malformed_inputs(self):
        private_pem = self.key.private_bytes(serialization.Encoding.PEM, serialization.PrivateFormat.PKCS8, serialization.NoEncryption())
        private_der = self.key.private_bytes(serialization.Encoding.DER, serialization.PrivateFormat.PKCS8, serialization.NoEncryption())
        public = self.public.public_bytes(serialization.Encoding.PEM, serialization.PublicFormat.SubjectPublicKeyInfo)
        der = self.cert.public_bytes(serialization.Encoding.DER)
        unsupported = ed25519.Ed25519PrivateKey.generate().public_key().public_bytes(serialization.Encoding.PEM, serialization.PublicFormat.SubjectPublicKeyInfo)
        for data in [private_pem, private_der, public + private_pem, unsupported, der + b"tail",
                     public + b"garbage", b"-----BEGIN PUBLIC KEY-----\nbroken", b"not a fingerprint"]:
            with self.subTest(data=data[:20]), self.assertRaises(ValueError):
                fingerprints(data)

    def test_matching_resume_needs_no_original_target_files(self):
        with tempfile.TemporaryDirectory() as directory:
            key, cert, out, cp = (Path(directory) / name for name in ("key", "cert", "out", "cp"))
            key.write_bytes(self.public.public_bytes(serialization.Encoding.DER, serialization.PublicFormat.SubjectPublicKeyInfo))
            cert.write_bytes(self.cert.public_bytes(serialization.Encoding.PEM))
            run(parser().parse_args(["generate", "--start", "1", "--end", "4", "--max-candidates", "1",
                                    "--match", str(key), "--match", str(cert), "--output", str(out), "--checkpoint", str(cp)]))
            key.unlink()
            cert.unlink()
            result = run(parser().parse_args(["resume", str(cp)]))
            self.assertEqual(result["emitted_total"], 1)
            self.assertEqual(json.loads(out.read_text())["timestamp"], 2)
