"""Regenerate public known-answer fixtures with independent reference programs.

Requires Linux glibc, Python 3.11+, c++ and a JDK. Does not import timestampseeds
or cryptography. Public points use a separate, simple affine curve implementation.
Run from any directory; writes JSON to stdout, never serializes private scalars.
"""
import ctypes
import hashlib
import json
import platform
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).parent
CURVES = {
    "p256": (
        int("ffffffff00000001000000000000000000000000ffffffffffffffffffffffff", 16), -3,
        int("ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551", 16),
        (int("6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296", 16),
         int("4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5", 16)),
        "3059301306072a8648ce3d020106082a8648ce3d03010703420004"),
    "secp256k1": (
        2**256 - 2**32 - 977, 0,
        int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16),
        (int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16),
         int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16)),
        "3056301006072a8648ce3d020106052b8104000a03420004"),
}
RNG_NAMES = ["lrand48", "drand48-state", "glibc-rand", "ansi-lcg31", "python-mt19937",
             "java-random", "msvc-rand", "cpp-mt19937"]


def multiply(value, prime, a, point):
    def add(left, right):
        if left is None:
            return right
        if right is None:
            return left
        x, y = left
        u, v = right
        if x == u and (y + v) % prime == 0:
            return None
        slope = ((3*x*x + a) * pow(2*y, -1, prime) if left == right else
                 (v-y) * pow(u-x, -1, prime)) % prime
        result_x = (slope*slope - x - u) % prime
        return result_x, (slope*(x-result_x) - y) % prime
    result = None
    while value:
        if value & 1:
            result = add(result, point)
        point = add(point, point)
        value >>= 1
    return result


class References:
    def __init__(self, directory):
        self.directory = directory
        self.libc = ctypes.CDLL(None)
        self.libc.srand48.argtypes = [ctypes.c_long]
        self.libc.lrand48.restype = ctypes.c_long
        self.libc.drand48.restype = ctypes.c_double
        self.libc.srand.argtypes = [ctypes.c_uint]
        for name in ["mt19937", "lcg_examples"]:
            subprocess.run(["c++", "-std=c++11", str(ROOT / (name + ".cpp")), "-o", str(Path(directory) / name)], check=True)
        subprocess.run(["javac", "-d", directory, str(ROOT / "RandomReference.java")], check=True)

    def java(self, seed, count, *extra):
        return subprocess.check_output(["java", "-cp", self.directory, "RandomReference", str(seed), str(count), *map(str, extra)], text=True)

    def words(self, name, seed, count):
        if name in ("lrand48", "drand48-state"):
            self.libc.srand48(seed)
            return [self.libc.lrand48() if name == "lrand48" else int(self.libc.drand48() * 2**48) for _ in range(count)]
        if name == "glibc-rand":
            self.libc.srand(seed)
            return [self.libc.rand() for _ in range(count)]
        if name == "python-mt19937":
            generator = random.Random(seed)
            return [generator.getrandbits(32) for _ in range(count)]
        if name == "java-random":
            return list(map(int, self.java(seed, count).splitlines()))
        args = ([str(Path(self.directory) / "mt19937")] if name == "cpp-mt19937" else
                [str(Path(self.directory) / "lcg_examples"), name])
        return list(map(int, subprocess.check_output([*args, str(seed), str(count)], text=True).splitlines()))

    def value(self, rng, seed, profile, skip, index):
        if profile.startswith("python-"):
            generator = random.Random(seed)
            for _ in range(skip):
                generator.getrandbits(32)
            for _ in range(index + 1):
                value = (generator.getrandbits(256) if profile == "python-getrandbits-reject" else
                         int.from_bytes(generator.randbytes(32), "big"))
            return value
        if profile == "java-nextbytes-reject":
            material = bytes.fromhex(self.java(seed, 32 * (index + 1), "bytes", skip))[-32:]
            return int.from_bytes(material, "big")
        width = 6 if rng == "drand48-state" else 4
        count = 32 if profile == "low-byte-reject" else (32 + width - 1) // width
        values = self.words(rng, seed, skip + (index + 1)*count)[skip + index*count:]
        material = (bytes(value & 255 for value in values) if profile == "low-byte-reject" else
                    b"".join(value.to_bytes(width, "little" if profile == "concat-le-reject" else "big") for value in values)[:32])
        return int.from_bytes(material, "little" if profile == "buffer-le-reject" else "big")


def build():
    if platform.libc_ver()[0] != "glibc":
        raise SystemExit("Fixture generation requires glibc Linux")
    cases = [(rng, 42, "concat-be-reject", 0, 0) for rng in RNG_NAMES]
    cases += [(rng, 1700000000, "low-byte-reject", 3, 2) for rng in RNG_NAMES]
    cases += [("drand48-state", 0xffffffff, "concat-le-reject", 7, 1),
              ("cpp-mt19937", 42, "buffer-le-reject", 1, 1),
              ("cpp-mt19937", 2, "concat-be-mod", 0, 0),
              ("python-mt19937", 2**64 + 42, "python-getrandbits-reject", 5, 3),
              ("python-mt19937", 2**64 + 42, "python-randbytes-reject", 5, 3),
              ("java-random", 2**48 + 42, "java-nextbytes-reject", 1, 2)]
    fixtures = []
    with tempfile.TemporaryDirectory() as directory:
        reference = References(directory)
        for rng, seed, profile, skip, index in cases:
            value = reference.value(rng, seed, profile, skip, index)
            for curve, (prime, a, order, generator, prefix) in CURVES.items():
                scalar = value % (order - 1) + 1 if profile == "concat-be-mod" else value
                public = None
                if 1 <= scalar < order:
                    x, y = multiply(scalar, prime, a, generator)
                    encoded_x, encoded_y = x.to_bytes(32, "big"), y.to_bytes(32, "big")
                    public = {"sec1_compressed_hex": (bytes([2 + y % 2]) + encoded_x).hex(),
                              "spki_sha256": hashlib.sha256(bytes.fromhex(prefix) + encoded_x + encoded_y).hexdigest()}
                fixtures.append({"rng": rng, "seed": seed, "profile": profile, "curve": curve,
                                 "skip_draws": skip, "key_index": index, "public": public})
    return {"reference": "tests/reference/build_fixtures.py; native APIs plus independent affine EC", "cases": fixtures}


if __name__ == "__main__":
    print(json.dumps(build(), indent=2) + "\n", end="")
