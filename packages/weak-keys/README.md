# @zksecurity/weak-keys

Offline weak-key research tools in one workspace, separate from the reusable math in
`@zksecurity/cheatsheets` and SageMath-compatible packages. **Research only; not an
assurance that an unmatched key is secure.** No network acquisition or operator probes.

## Layout and migration

- `src/`: TypeScript screening API/CLI and the timestamp engine in `src/timestampseeds/`.
- `timestampseeds/`: complete usage guide, original public synthetic fixtures and
  independent native reference sources with TypeScript drivers.
- `MIGRATION.json`: upstream source commit and historical import hashes. The
  initial import from feature commit `6e68fdc` is preserved in Git history;
  the executable implementation and tests are now TypeScript.

All eight RNGs and profiles retain seed handling, draw consumption, skipped draws,
successive keys, public-key normalization, exact windows and resumable checkpoints.
The existing TypeScript CPython RNG and EC arithmetic are reused. No Python
interpreter, Python source, cryptography package or uv project is required.

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
Unsupported curves are not silently substituted. At most 256 records eligible for work per API
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

## Timestamp-seeded hypotheses (TypeScript)

Install Bun, then run from the repository root:

```sh
bun run weak-keys:timestampseeds estimate --start 1 --end 3
bun run weak-keys:timestampseeds --help
```

See [the complete engine README](timestampseeds/README.md) for generate/match/resume,
profiles, skipped draws, successive keys, PEM/DER targets and checkpoint recovery.
The CLI keeps Linux/macOS advisory locks through system libc via Bun FFI. Locks are
released by the OS after crashes; schema-1/2 checkpoints remain compatible.
The native Microsoft CRT comparison runs separately on Windows.

Timestamp candidates are **public outputs only**, but a matching seed/profile
association can make a private value reproducible. Keep candidates, target material,
matches, logs and checkpoints confidential and outside git, preferably under the
ignored `runtime/` directory or outside the checkout. The root convenience script invokes the TypeScript CLI directly and preserves the
caller working directory. Relative target/output/checkpoint paths remain relative
to that directory.

## Tests and build

```sh
# Repository root
bun run test:weak-keys
bun run test:weak-keys:native
bun run --cwd packages/weak-keys build
bun tests/weak-keys-typecheck.ts
bun x biome check packages/weak-keys/src tests/weak-keys-typecheck.ts
# Native comparisons require c++, a JDK and (on Linux) glibc.
```

Tests use synthetic inputs only, a live TypeScript small-integer oracle, native
C++/Java RNG output, independently computed affine EC points, and fresh replay
seeds. The original 44 public fixtures remain regression controls. The CPython
seed profile uses a standalone C++ integer-seeding adapter with the standard MT
engine; it does not invoke Python. Linux compares all RNGs and regenerates the
original fixture cases live; macOS tests seven RNGs (glibc rand is Linux-only).
Missing tools fail the selected native suite, rather than skipping comparisons.
Set `WEAK_KEYS_SEED` to replay properties; `CXX`, `JAVA` and `JAVAC` may select
native toolchain executables. Generated passing output is never committed.

The current dependency graph has existing TypeScript errors. `typecheck` retains
normal `tsc --noEmit` semantics and does **not** claim to pass. The separate scoped
regression command compares live diagnostics to the existing cheatsheets graph,
fails on new errors/tool failures/global compiler errors, and explicitly reports
inherited dependency errors; it is not a full-repository clean-typecheck claim.

No collected registry data or real-key findings are included. Future collection
adapters, signature corpus work and additional attacks should remain opt-in,
source-provenance-aware and independently budgeted rather than widening defaults.
