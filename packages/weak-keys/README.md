# @zksecurity/weak-keys

Offline weak-key research tools in one workspace, separate from the reusable math in
`@zksecurity/cheatsheets` and SageMath-compatible packages. **Research only; not an
assurance that an unmatched key is secure.** No network acquisition or operator probes.

## Layout and migration

- `src/`: TypeScript public-only screening API and local JSON CLI.
- `timestampseeds/`: the complete Python timestamp engine, CLI, tests, independent
  native references, public synthetic fixtures and locked dependencies.
- `MIGRATION.json`: source repository, branch, exact commit and SHA-256 of every
  imported file. The initial import is byte-for-byte from
  [`zksecurity/timestampseeds`](https://github.com/zksecurity/timestampseeds), feature
  commit `6e68fdc`, not its older `main`. These are migration provenance hashes, not
  a claim that the engine cannot evolve here.

The Python implementation is deliberately retained rather than silently replacing
it with a new TS RNG approximation. This preserves seed handling, draw consumption,
Java/MSVC/C++ references, derivation profiles, skipped draws/successive keys, public-key
normalization, UTC windows, locking and resumable checkpoints. TypeScript does not
invoke Python implicitly; the two commands have separate, explicit runtimes.

The original timestamp repository is not deleted or modified by this migration.
Its imported tree had no separate LICENSE or project license declaration; this
owner-authorized consolidation is covered by this repository's GPL-2.0-or-later
license. No third-party authorship or license text has been invented.

## TypeScript screening

From the repository root (`bun install --frozen-lockfile --ignore-scripts` first):

```sh
bun run weak-keys --help
printf '[{"type":"rsa","n":"27","e":"5"}]' | bun run weak-keys --input -
```

Or use the library inside this checkout:

```ts
import { screenPublicKeys } from '@zksecurity/weak-keys';

const report = screenPublicKeys([{ type: 'rsa', n: 27n, e: 5n }]);
// report.records[0].checks: prime-modulus no-match, perfect-power matched;
// Fermat and EC checks not-run. No factors, private exponents or input values.
```

Records are RSA `{type:'rsa', n, e}` with bigint or decimal/`0x` integer strings, or
EC `{type:'ec', curve, point}` with `curve` in `p256`, `p384`, `secp256k1` and a SEC1
hex public point. JavaScript numbers are deliberately not accepted for RSA integers.
Unsupported curves are not silently substituted. At most 10,000 records per API
report; the CLI additionally limits local input to 1 MiB. No target discovery.

### Implemented checks

| Check ID | Default? | Meaning and scope |
|---|---|---|
| `rsa-prime-modulus` | Yes | Reuses existing PARI-port BPSW screen. A match is **probable**, not a primality certificate; rejection is composite. Does not test hidden RSA factors. |
| `rsa-perfect-power` | Yes | Exact general perfect powers via existing `Z_isanypower`, not only squares or prime bases. |
| `rsa-fermat` | No | Reuses existing Fermat search; a match requires a verified proper factorization internally. Factors are never returned. |
| `ec-small-scalar` | No | Reuses existing EC helper for a declared +/-small-scalar interval, verifies valid public points, returns no recovered scalar/pattern. |

Examples with explicitly bounded extra work:

```sh
bun run weak-keys --input /absolute/path/to/local-public-keys.json \
  --checks rsa-fermat --fermat-iterations 256 --max-fermat-steps 8192 --max-keys 32
bun run weak-keys --input /absolute/path/to/local-ec-keys.json \
  --checks ec-small-scalar --ec-table-bits 16 --max-ec-keys 8
```

Defaults/hard maxima: `maxKeys=32/256`, `maxModulusBits=4096/8192`,
`fermatIterations=256/4096` per key, `maxFermatSteps=8192/65536` per batch,
`ecTableBits=8/16`, `maxEcKeys=8/32`. Zero key/Fermat/EC-key budgets do not run
searches. EC table bits must be 1..16. Limits bound counts and input sizes, **not
interruptible wall-clock time**. Run inside your own time/resource limit if needed.
The helper may reuse an existing larger EC table; matches are filtered to the
explicitly declared +/-2^bits interval so process history cannot broaden results.
No all-pairs EC comparisons, Pollard p-1, kangaroo or complete attack catalogue runs
are enabled by this adapter.

### Reports and interpretation

Only input indices, fixed check/status/confidence enums and configured numeric
search scopes are returned. No identifiers, moduli, points, fingerprints, factors,
private exponents, scalars, pattern words, recovery messages or arbitrary metadata.
Do not forward the underlying cheatsheet recovery reports as public screening output.

| Status | Meaning |
|---|---|
| `matched` | The declared check found a condition; consult `confidence`. |
| `no-match` | Completed the declared check/range without a match; not a universal safety claim. |
| `not-run` | Check was not selected. |
| `budget-exhausted` | No work could start for that record/check under remaining budget. |
| `malformed` | Required fields, integer/exponent constraints or EC encoding/point are invalid. |
| `unsupported` | Unknown key/curve type, wrong check family or RSA input beyond the bit limit. |

A bounded negative is not proof against other ranges/profiles. Reports intentionally
distinguish exact matches, probable-prime matches and bounded search negatives.
Even a numeric finding count can be sensitive when tied to a source: publication is
still the caller's decision. Repeated certificates/keys do not themselves establish
weakness; this adapter does not classify duplicate moduli as vulnerabilities.

## Timestamp-seeded hypotheses (Python)

Install Python 3.11+ and `uv`, then from the repository root:

```sh
bun run weak-keys:timestampseeds estimate --start 1 --end 3
```

Or enter the explicit Python subproject:

```sh
cd packages/weak-keys/timestampseeds
uv sync --locked
uv run --locked timestampseeds --help
uv run --locked timestampseeds estimate --start 1 --end 3
```

See [the complete engine README](timestampseeds/README.md) for generate/match/resume,
profiles, skipped draws, successive keys, PEM/DER targets, exact reference vectors,
checkpoint crash recovery and limitations. The default interpreter pin is Python
3.13; Python 3.11 is also supported (`uv ... --python 3.11`, or `UV_PYTHON=3.11`).
The CLI remains Unix-only (`fcntl.flock`); a separate native CRT reference test runs
on Windows. Native Java/C++/glibc tests report unavailable toolchains as skips;
CI supplies the required native toolchains and regenerates the independent fixtures.

Timestamp candidates are **public outputs only**, but a matching seed/profile
association can make a private value reproducible. Keep candidates, target material,
matches, logs and checkpoints confidential and outside git, preferably under the
ignored `runtime/` directory or outside the checkout. The root convenience script uses `uv run --project` and preserves the caller
working directory: relative target/output/checkpoint paths stay relative to the
repository root. Direct invocations from inside the Python subproject remain
relative to that directory. Absolute paths are unambiguous in either case.

## Tests and build

```sh
# Repository root
bun run test:weak-keys
bun run test:weak-keys:python
bun run --cwd packages/weak-keys build
bun tests/weak-keys-typecheck.ts
bun x biome check packages/weak-keys/src tests/weak-keys-typecheck.ts
# Inside timestampseeds/: uv build; uv run --locked python tests/reference/build_fixtures.py
```

TS tests use synthetic inputs only, a live small-integer Python oracle with fresh
replay seeds, and Node/OpenSSL-generated synthetic EC points. Missing reference
runtimes fail those comparisons; generated passing transcripts are not committed.
The copied Python tests/native reference assets are preserved unchanged.

The current dependency graph has existing TypeScript errors. `typecheck` retains
normal `tsc --noEmit` semantics and does **not** claim to pass. The separate scoped
regression command compares live diagnostics to the existing cheatsheets graph,
fails on new errors/tool failures/global compiler errors, and explicitly reports
inherited dependency errors; it is not a full-repository clean-typecheck claim.

No collected registry data or real-key findings are included. Future collection
adapters, signature corpus work and additional attacks should remain opt-in,
source-provenance-aware and independently budgeted rather than widening defaults.
