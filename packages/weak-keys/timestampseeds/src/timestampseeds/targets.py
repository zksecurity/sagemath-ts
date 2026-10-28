"""Normalize local public inputs to SHA256 of canonical, uncompressed SPKI DER."""
import hashlib
import re
from pathlib import Path

from cryptography import x509
from cryptography.exceptions import UnsupportedAlgorithm
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import ec


PEM_BLOCK = re.compile(rb"-----BEGIN ([A-Z0-9 ]+)-----\s+.*?-----END \1-----", re.DOTALL)


def valid_fingerprint(value) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def fingerprint(public) -> str:
    if not isinstance(public, ec.EllipticCurvePublicKey) or public.curve.name not in ("secp256r1", "secp256k1"):
        raise ValueError("Targets must use P-256 or secp256k1 public keys")
    # Reconstruct from coordinates to normalize compressed input SPKIs as well.
    public = public.public_numbers().public_key()
    der = public.public_bytes(serialization.Encoding.DER, serialization.PublicFormat.SubjectPublicKeyInfo)
    return hashlib.sha256(der).hexdigest()


def fingerprints(data: bytes) -> set[str]:
    if b"-----BEGIN " in data:
        result, end = set(), 0
        for block in PEM_BLOCK.finditer(data):
            if data[end:block.start()].strip():
                raise ValueError("Unexpected content outside PEM blocks")
            label, payload = block[1], block[0]
            if label == b"PUBLIC KEY":
                public = serialization.load_pem_public_key(payload)
            elif label == b"CERTIFICATE":
                public = x509.load_pem_x509_certificate(payload).public_key()
            else:
                raise ValueError("Only PUBLIC KEY and CERTIFICATE PEM blocks are accepted")
            result.add(fingerprint(public))
            end = block.end()
        if not result or data[end:].strip():
            raise ValueError("Malformed PEM input")
        return result
    # DER keys and certificates are both ASN.1 SEQUENCEs. Reject trailing DER data.
    if data.startswith(b"\x30") and not _looks_like_text(data):
        _validate_der_length(data)
        try:
            public = serialization.load_der_public_key(data)
        except ValueError:
            public = x509.load_der_x509_certificate(data).public_key()
        return {fingerprint(public)}
    result = set()
    for line in data.decode("ascii").splitlines():
        line = line.strip().lower()
        if not line or line.startswith("#"):
            continue
        if not valid_fingerprint(line):
            raise ValueError("Expected SHA256 SPKI fingerprints, PEM public keys/certificates, or a DER public key/certificate")
        result.add(line)
    return result


def _looks_like_text(data: bytes) -> bool:
    return all(byte in (9, 10, 13) or 32 <= byte < 127 for byte in data)


def _validate_der_length(data: bytes) -> None:
    if len(data) < 2:
        raise ValueError("Truncated DER input")
    length, header = data[1], 2
    if length & 128:
        count = length & 127
        if count == 0 or count > 4 or len(data) < 2 + count:
            raise ValueError("Invalid DER length")
        length, header = int.from_bytes(data[2:2+count], "big"), 2 + count
    if header + length != len(data):
        raise ValueError("DER input must contain exactly one object")


def targets(paths: str | list[str] | None) -> set[str] | None:
    if paths is None:
        return None
    if isinstance(paths, str):
        paths = [paths]
    result = set()
    for path in paths:
        try:
            result.update(fingerprints(Path(path).read_bytes()))
        except (ValueError, UnsupportedAlgorithm) as exc:
            raise ValueError(f"Invalid target file {path}: {exc}") from exc
    return result
