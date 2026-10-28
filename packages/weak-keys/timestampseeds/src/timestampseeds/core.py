import base64
import hashlib
import re
from collections import deque
from datetime import datetime, timezone
from itertools import islice

from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import ec

from .rngs import WIDTHS, words

PROFILES = ("concat-be-reject", "concat-le-reject", "concat-be-mod", "low-byte-reject",
            "buffer-le-reject", "python-getrandbits-reject", "python-randbytes-reject",
            "java-nextbytes-reject")
PROFILE_RNG = {"python-getrandbits-reject": "python-mt19937",
               "python-randbytes-reject": "python-mt19937", "java-nextbytes-reject": "java-random"}
MAX_SKIP_DRAWS = 1_000_000
MAX_KEYS_PER_SEED = 1_000
CURVES = {
    "p256": (ec.SECP256R1, int("FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551", 16)),
    "secp256k1": (ec.SECP256K1, int("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16)),
}


def parse_time(value: str, unit: str) -> int:
    """Integer epoch units or timezone-explicit ISO8601. Endpoints are exact."""
    if unit not in ("seconds", "milliseconds"):
        raise ValueError("Unknown timestamp unit")
    try:
        return int(value)
    except ValueError:
        # datetime truncates fractions after six digits. Never silently change a seed.
        for fraction in re.findall(r"[.,](\d+)", value):
            if any(digit != "0" for digit in fraction[6:]):
                raise ValueError(f"Timestamp is not aligned to {unit}")
        # datetime also ignores fractional offsets whose integral part is zero.
        if re.search(r"[+-]\d{2}:?\d{2}(?::?\d{2})?[.,]", value):
            raise ValueError("Fractional timezone offsets are not supported")
        dt = datetime.fromisoformat(value.replace("Z", "+00:00"))
        if dt.tzinfo is None or dt.utcoffset() is None:
            raise ValueError("ISO timestamps require Z or an explicit timezone offset")
        delta = dt.astimezone(timezone.utc) - datetime(1970, 1, 1, tzinfo=timezone.utc)
        micros = (delta.days * 86400 + delta.seconds) * 1_000_000 + delta.microseconds
        divisor = 1_000_000 if unit == "seconds" else 1000
        if micros % divisor:
            raise ValueError(f"Timestamp is not aligned to {unit}")
        return micros // divisor


def validate_profile(rng: str, profile: str) -> None:
    if rng not in WIDTHS:
        raise ValueError(f"Unknown RNG: {rng}")
    if profile not in PROFILES:
        raise ValueError("Unknown derivation profile")
    if profile in PROFILE_RNG and rng != PROFILE_RNG[profile]:
        raise ValueError(f"{profile} requires --rng {PROFILE_RNG[profile]}")


def draws_per_key(rng: str, profile: str) -> int:
    validate_profile(rng, profile)
    return 32 if profile == "low-byte-reject" else (32 + WIDTHS[rng] - 1) // WIDTHS[rng]


def scalar_from_words(stream, profile: str, order: int) -> int | None:
    if profile not in PROFILES or type(order) is not int or order < 2:
        raise ValueError("Invalid derivation profile or subgroup order")
    # Fixed width per output (including leading zeros), concatenate FIRST 32 bytes.
    little_words = {"concat-le-reject", "python-getrandbits-reject", "python-randbytes-reject", "java-nextbytes-reject"}
    endian = "little" if profile in little_words else "big"
    material = bytearray()
    while len(material) < 32:
        value, width = next(stream)
        if profile == "low-byte-reject":
            material.append(value & 255)
        else:
            material.extend(value.to_bytes(width, endian))
    scalar_endian = "little" if profile in ("buffer-le-reject", "python-getrandbits-reject") else "big"
    value = int.from_bytes(material[:32], scalar_endian)
    if profile.endswith("-mod"):
        return value % (order - 1) + 1
    # Reject this entire candidate, do NOT retry the RNG or silently reduce.
    return value if 1 <= value < order else None


def stream_at(rng: str, seed: int, profile: str, skip_draws: int = 0, key_index: int = 0):
    if type(skip_draws) is not int or not 0 <= skip_draws <= MAX_SKIP_DRAWS:
        raise ValueError(f"Skipped draws must be between 0 and {MAX_SKIP_DRAWS}")
    if type(key_index) is not int or not 0 <= key_index < MAX_KEYS_PER_SEED:
        raise ValueError(f"Key index must be between 0 and {MAX_KEYS_PER_SEED - 1}")
    count = skip_draws + key_index * draws_per_key(rng, profile)
    stream = words(rng, seed)
    deque(islice(stream, count), maxlen=0)
    return stream


def scalar(rng: str, seed: int, profile: str, order: int, skip_draws: int = 0, key_index: int = 0) -> int | None:
    return scalar_from_words(stream_at(rng, seed, profile, skip_draws, key_index), profile, order)


def candidate(timestamp: int, unit: str, rng: str, profile: str, curve: str,
              skip_draws: int = 0, key_index: int = 0, stream=None) -> dict | None:
    curve_type, order = CURVES[curve]
    validate_profile(rng, profile)
    value = (scalar(rng, timestamp, profile, order, skip_draws, key_index) if stream is None
             else scalar_from_words(stream, profile, order))
    if value is None:
        return None
    public = ec.derive_private_key(value, curve_type()).public_key()
    der = public.public_bytes(serialization.Encoding.DER, serialization.PublicFormat.SubjectPublicKeyInfo)
    point = public.public_bytes(serialization.Encoding.X962, serialization.PublicFormat.CompressedPoint)
    # Deliberately never serialize the scalar or any private-key object.
    record = {
        "timestamp": timestamp, "unit": unit, "rng": rng, "profile": profile,
        "curve": curve, "spki_sha256": hashlib.sha256(der).hexdigest(),
        "spki_der_base64": base64.b64encode(der).decode("ascii"),
        "sec1_compressed_hex": point.hex(),
    }
    # Preserve the original default record format, including schema-1 resume bytes.
    if skip_draws or key_index:
        record.update(skip_draws=skip_draws, key_index=key_index)
    return record


def combination_at(config: dict, index: int) -> tuple[int, str, str, str]:
    """Stable order: timestamp, RNG, profile, curve. No range materialization."""
    ncurves, nprofiles, nrngs = len(config["curves"]), len(config["profiles"]), len(config["rngs"])
    curve = config["curves"][index % ncurves]
    index //= ncurves
    profile = config["profiles"][index % nprofiles]
    index //= nprofiles
    rng = config["rngs"][index % nrngs]
    index //= nrngs
    return config["start"] + index * config["step"], rng, profile, curve


def total_candidates(config: dict) -> int:
    return (timestamp_count(config) * len(config["rngs"]) * len(config["profiles"]) * len(config["curves"])
            * len(config.get("skip_draws", [0])) * config.get("keys_per_seed", 1))


def timestamp_count(config: dict) -> int:
    return (config["end"] - config["start"] + config["step"] - 1) // config["step"]


def attempt_at(config: dict, index: int) -> tuple[int, str, str, str, int, int]:
    """Order: timestamp, RNG, profile, curve, skipped draws, successive key index."""
    key_index = index % config.get("keys_per_seed", 1)
    index //= config.get("keys_per_seed", 1)
    skips = config.get("skip_draws", [0])
    skip = skips[index % len(skips)]
    index //= len(skips)
    return (*combination_at(config, index), skip, key_index)
