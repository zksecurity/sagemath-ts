# timestampseeds

Offline generation of **public-key candidates** for explicit timestamp-seeded RNG hypotheses. It never connects to a target and never outputs private scalars, private-key files, or credentials. It can compare generated public fingerprints against a local public-key corpus.

**A timestamp is not a universal way to predict a key.** An actual generator's seed handling, RNG algorithm, output consumption, byte order, scalar mapping, and curve must match a profile here. Failure to match proves nothing about a different generator. This does not emulate modern OS CSPRNGs, OpenSSL key generation, Debian's historical OpenSSL flaw, or every program calling `srand48`.

## Install / test

Python 3.11+ on Linux or macOS (uses `fcntl.flock` for output-file locking;
Windows is not supported). Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then:

```sh
uv sync --locked
uv run --locked python -m unittest discover -s tests -v
uv run --locked timestampseeds --help
```

The checkout defaults to Python 3.13 via `.python-version`. `uv` creates `.venv`
and installs the package in editable mode. `uv.lock` records exact runtime
dependencies, including `cryptography==50.0.2`; keep it in version control.
Tests use the standard-library `unittest` runner, so no test dependencies are needed.
Native reference tests use `c++` and a JDK when installed. The libc reference
test runs only on glibc Linux. Missing native toolchains are reported as skips;
the checked-in public-key fixtures run everywhere. CI tests Python 3.11 and 3.13
on Linux/macOS, regenerates fixtures on Linux, and separately checks the RNG
against the native Microsoft CRT on Windows (the CLI remains Unix-only).

Use `uv run timestampseeds ...` as below, or activate `.venv` to use the bare
`timestampseeds` command. `uv run python -m timestampseeds ...` is also supported.
After intentionally changing dependencies in `pyproject.toml`, run `uv lock`
and `uv sync --locked`.

## Usage

First estimate the number of attempts; no curve multiplication occurs:

```
uv run timestampseeds estimate
# Default window: 2000-01-01 UTC through the current UTC second, end-exclusive.

uv run timestampseeds estimate --start 2024-01-01T00:00:00Z --end 2025-01-01T00:00:00Z \
  --rng lrand48 glibc-rand python-mt19937 --curve p256 secp256k1
```

Generate a narrower range. The default per-invocation cap is **10,000 attempts**, including rejected scalars and nonmatches. The broad default window is never traversed without that cap. Set an explicit `--max-candidates` to change it:

```
uv run timestampseeds generate --start 2024-03-01T12:00:00Z --end 2024-03-01T12:05:00Z \
  --rng lrand48 glibc-rand ansi-lcg31 python-mt19937 \
  --profile concat-be-reject concat-le-reject --curve p256 secp256k1 \
  --output candidates.jsonl --checkpoint progress.json

uv run timestampseeds resume progress.json --max-candidates 10000
```

Seconds/milliseconds are explicit seed hypotheses, not just display formats:

```
uv run timestampseeds generate --unit milliseconds --start 1700000000000 --end 1700000001000 \
  --step 10 --rng lrand48 --output candidates-ms.jsonl
```

All ranges are `[start, end)`. `--step` is a positive integer in selected units. ISO timestamps require a timezone (`Z` or offset) and must align exactly to the selected unit, including digits beyond microsecond precision. Fractional timezone offsets are rejected. Seeds are nonnegative integers. Integer endpoints mean epoch values in the selected unit, **not** years. Duplicate RNG/profile/curve/skipped-draw arguments are removed.

By default, JSONL goes to stdout and a summary goes to stderr. File outputs must be new and are created mode 0600, as are checkpoint files. Checkpoints require file output, save resolved ranges/algorithm choices and selected matching fingerprints, and record a SHA256 of committed public output. Resume validates the full checkpoint before changing output, verifies that prefix and truncates any uncommitted crash tail. Checkpoints are atomic and file output is flushed/fsynced first. Default checkpoint interval: 100 attempts. Ctrl-C retains the last consistent checkpoint rather than saving possibly partial work; resume truncates any later output tail and recomputes up to that checkpoint interval. If interrupted before the first checkpoint exists, restart with new output/checkpoint paths. `--max-candidates` is a fresh cap per resume. Do not delete/modify committed output. Unix advisory locks protect both output and checkpoint ownership. The persistent `<checkpoint>.lock` file must not be removed while writers are running; its presence alone does not mean a writer is active. No scalar or RNG state is checkpointed. Schema-1 checkpoints remain resumable and are upgraded to schema 2 on the next save.

### Prior draws and successive keys

```sh
uv run timestampseeds generate --start 1700000000 --end 1700000010 \
  --rng java-random --profile java-nextbytes-reject \
  --skip-draws 0 1 8 --keys-per-seed 4 \
  --output java-candidates.jsonl --checkpoint java-checkpoint.json
```

Each `--skip-draws` value is a separate hypothesis, from 0 to 1,000,000 inclusive.
It discards that many RNG **output words**, after the algorithm's normal seed
initialization/warm-up. Java counts `nextInt()` calls, Python counts
`getrandbits(32)` calls, and the other RNGs count the outputs in the table below.
These are not counts of bytes or arbitrary API calls such as `nextDouble()`.

`--keys-per-seed` is 1..1,000 (default 1). A fresh stream starts for each
timestamp/RNG/profile/curve/skip combination, then continues through that many
32-byte candidate blocks. A rejected block still consumes draws and a key index;
there is no rejection retry. The implementation reuses the stream within a run
and replays only the needed prefix when resuming mid-sequence. For six-byte
words, unused bytes at the end of each candidate block are discarded.

### Progress and throughput

Generation and resume emit JSON progress events to stderr every five seconds
at attempt boundaries. Set `--progress-every 1` for more frequent reports or
`--progress-every 0` to disable them. Candidate JSONL is unaffected.
Reports include attempt count, emitted count, rejected/nonmatching attempts,
elapsed seconds, and attempts/second. Progress ETA is for the **current capped
invocation**, not the entire search. The final summary always includes measured
throughput; rejection/nonmatch counts and elapsed time reset on resume.

### Match public fingerprints offline

`targets.txt` contains one SHA256 SPKI-DER fingerprint (64 hex characters) per line; blank lines/comments are allowed:

```
uv run timestampseeds generate --start 1700000000 --end 1700000100 \
  --rng lrand48 glibc-rand --match targets.txt --output matches.jsonl
```

Only matching public records are emitted, but **all attempts count toward the cap**. No matches is not an exhaustive weak-key test. Matches identify an explicit seed/RNG/profile hypothesis; keep affected-key associations confidential for human verification/disclosure. The tool makes no live verification attempts.

`--match` also accepts a PEM bundle containing `PUBLIC KEY` and/or `CERTIFICATE`
blocks, or one DER public key/certificate. Repeat it to combine files:

```sh
uv run timestampseeds generate --start 1700000000 --end 1700000100 \
  --match public-keys.pem --match certificate.der --match targets.txt \
  --output matches.jsonl --checkpoint match-checkpoint.json
```

Keys must use P-256 or secp256k1. Private keys, unsupported key types, malformed
bundles, and trailing DER objects are rejected. Certificates are only containers
for public keys; signatures, trust, and expiry are not checked. Keys are
normalized to uncompressed SPKI DER before hashing, so compressed SPKI inputs
match the same public point. Text fingerprints must already use that canonical
encoding; certificate fingerprints are different. Selected hashes are saved in
the checkpoint, so resume does not need the source target files. Empty fingerprint
files intentionally match nothing. Duplicate keys are deduplicated across files.

## Exact RNG variants

| Name | Definition / seed semantics | Serialized output |
|---|---|---|
| `lrand48` | POSIX `srand48`: state `(low32(seed) << 16) | 0x330e`, then `state = (0x5deece66d * state + 0xb) mod 2^48` | Each updated state's top 31 bits in **4 bytes**, including leading zeros |
| `drand48-state` | Same 48-bit state as above | Full updated state in **6 bytes**; exactly the numerator of `drand48() * 2^48`. This is **not** serialization of IEEE754 double bytes |
| `glibc-rand` | glibc TYPE_3 additive `random`/`rand`, degree 31, separation 3, Park–Miller signed initialization, 310 warm-up iterations; seed low32, zero becomes one | 31-bit result in **4 bytes** |
| `ansi-lcg31` | Explicit common example `state = (1103515245 * state + 12345) mod 2^31`; seed low31 | Updated state in **4 bytes**. This is **not** a portable definition of C `rand()` |
| `python-mt19937` | Python `random.Random(integer_seed).getrandbits(32)`; full nonnegative seed integer, Python integer seed expansion | Each result in **4 bytes**; **not** C++ `std::mt19937(seed)` |
| `java-random` | `java.util.Random(long)`: initial state `(seed ^ 0x5deece66d) mod 2^48`, same 48-bit LCG step as `srand48`, different seeding and extraction | `nextInt()` bit pattern as an unsigned **4-byte** word (top 32 updated state bits) |
| `msvc-rand` | Microsoft CRT: seed low32; `state = (214013 * state + 2531011) mod 2^32` | `(state >> 16) & 0x7fff` in **4 bytes**, including leading zeros |
| `cpp-mt19937` | C++ `std::mt19937(integer_seed)` with low32 seed and standard integer initialization; no `seed_seq` or distribution adapter | Full unsigned result in **4 bytes** |

Linux glibc `lrand48`, `drand48`, and `rand` outputs are cross-tested with libc using seeds including `0`, high-bit uint32 values, and values above `2^32`. libc tests skip on non-glibc systems. Fixed glibc/ANSI/Python vectors are always tested.

Reference definitions: [POSIX drand48/srand48](https://pubs.opengroup.org/onlinepubs/9699919799/functions/drand48.html), [glibc random implementation](https://sourceware.org/git/?p=glibc.git;a=blob;f=stdlib/random_r.c), [Python random documentation](https://docs.python.org/3/library/random.html), [Java Random specification](https://docs.oracle.com/en/java/javase/22/docs/api/java.base/java/util/Random.html), [Microsoft rand](https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/rand), [C++ MT specification](https://eel.is/c++draft/rand.eng.mers). This package's explicit profiles, rather than a generic RNG name, define each candidate.

## Exact derivation profiles

Each profile specifies exactly how to consume a 32-byte candidate block. Unless
stated otherwise, the full buffer is an unsigned **big-endian scalar integer**.
All `-reject` profiles accept only `1 <= d < curve_order`:

- `concat-be-reject`: each word serialized big-endian; accept only `1 <= d < curve_order`; otherwise reject that candidate without further draws.
- `concat-le-reject`: each word serialized little-endian, then whole 32-byte buffer interpreted big-endian; same rejection rule. This is **not** whole-buffer little-endian interpretation.
- `concat-be-mod`: each word serialized big-endian; map with `d = value % (curve_order - 1) + 1`.
- `low-byte-reject`: consume 32 words, keeping only `word & 255` from each. This models `rand() % 256` for these nonnegative outputs.
- `buffer-le-reject`: serialize words big-endian, take the first 32 bytes, then interpret the **whole buffer little-endian**.
- `python-getrandbits-reject`: Python `Random(seed).getrandbits(256)` as the scalar integer; eight successive 32-bit words, first word least significant. Requires `python-mt19937`.
- `python-randbytes-reject`: Python `Random(seed).randbytes(32)` interpreted big-endian. Requires `python-mt19937`.
- `java-nextbytes-reject`: Java `new Random(seed).nextBytes(new byte[32])` interpreted big-endian. Requires `java-random`.

`python-randbytes-reject` and `java-nextbytes-reject` are named API equivalents
of `concat-le-reject` for their respective RNGs at this fixed 32-byte length.
They remain separate provenance hypotheses; identical records are not collapsed.
Incompatible RNG/profile combinations fail validation, rather than silently
skipping work. Use separate invocations for native profiles with different RNGs.

These are deliberate generation hypotheses, **not assertions that keygen libraries
use these mappings**. There is no automatic entropy mixing, hashing, retry,
floating-point conversion, or hidden modulus reduction. Except for low-byte
extraction, `drand48-state` consumes six 6-byte words and truncates the sixth to
its first two serialized bytes. Rejection and modulo profiles are not interchangeable.

Supported curves: NIST P-256 and secp256k1. Their subgroup orders are explicit constants. Cryptography derives a transient private-key object only to compute the public point; private material is never serialized or logged. Python/OpenSSL memory is not guaranteed zeroized.

## Output / scaling

Each public JSONL record has timestamp, unit, RNG, profile, curve, compressed SEC1 point, SPKI DER (base64), and SHA256 SPKI fingerprint. When skipped draws or key index are nonzero, `skip_draws` and zero-based `key_index` are also included. If absent, both are zero; the default record format stays backward compatible. Records contain enough provenance to reproduce the hypothesis. **Seed metadata is sensitive when associated with a real weak public key** even though it is not a private-key field: do not publish such matches indiscriminately.

Stable attempt order: timestamp → RNG → profile → curve → skipped-draw choice → key index. Work and memory are streamed, except optional target files loaded one at a time and their matching fingerprint set kept in memory. Timestamp count × RNGs × profiles × curves × skipped-draw choices × keys per seed determines attempts. Seconds versus milliseconds changes the range size by up to 1000×. Different seeds can produce the same public key (e.g. 32-bit seed truncation or glibc 0/1); records intentionally preserve distinct provenance rather than keeping an unbounded global dedup set. No claim to exhaustive seed or RNG coverage.

## Independent reference fixtures

`tests/fixtures/public_keys.json` contains public points and SPKI hashes generated
without importing this package or `cryptography`. The generator uses glibc,
`std::mt19937`, `java.util.Random`, Python's native byte APIs, standalone C++ LCG
examples, and a separate affine elliptic-curve implementation. It covers all RNGs
and profiles on both curves, including prior consumption and successive keys.
Boundary tests separately cover zero/order rejection and six-byte truncation.
Native Windows CI checks the Microsoft LCG against `ucrtbase`.

To regenerate on glibc Linux with Python, `c++`, and a JDK installed:

```sh
uv run python tests/reference/build_fixtures.py > /tmp/public_keys.json
diff -u tests/fixtures/public_keys.json /tmp/public_keys.json
```

No real-world key corpus has been matched as part of repository tests; all positive matches are synthetic. The implementation must be evaluated on known generator-specific weak-key samples before treating it as a broad detection tool.
