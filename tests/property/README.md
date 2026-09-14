# Live comparative property tests

Generate inputs once, execute the same inputs against the original implementation
and TypeScript, and compare exact results, exception classes, messages and native
state. Successful runs do not write transcripts or expected-output snapshots.

```sh
bun run test:property -- --case arith
bun run test:property -- --case arith --seed 123 --runs 100
bun run test:property -- --cases arith,finite_fields --seed 123
bun run test:property:fast
bun run test:property:slow
```

The default seed comes from the operating system and is printed before execution.
`--runs` is the minimum number of fresh trials per generator (default 25); a
recipe suite also visits each constructor variant at least once. Explicit
regressions and retained regression seeds run as well. The same seed, generator
version, case definitions and run count reproduce the same inputs; selecting other
areas does not change an area's generated inputs.

## References are required

Live comparisons require SageMath (`sage` on PATH). Native adapters compile and
execute the bundled PARI, NTL, FLINT and other original sources where the relevant
Sage entry point delegates to them. The source/build adapters are in `python/` and
`native/`; compiled references are cached in the operating system's temporary
folder. They are never replaced by the TypeScript implementation as an oracle.

Some adapters also require a native compiler or CPython 3.12 (`PYTHON312` may name
its executable). Missing references, dispatch entries or result rows fail the run.
There is no cached-output fallback: `--generate` and `--typescript-only` have been
retired. Native comparison tests in package directories also require the original
runtime. CI should install these references before running comparative tests.

## Adding a property

An area has three files:

- `cases/<area>.cases.json`: compact generator definitions and focused regressions.
- `python/areas/<area>.py`: exports `FUNCTIONS` and optional `FORMATTERS`.
- `typescript/areas/<area>.ts`: exports matching `functions` and optional `formatters`.

Both area modules receive concrete arguments. The comparison harness generates
those arguments once with `seeded.ts`, avoiding differences between two input
implementations. Large integers are encoded as decimal strings and decoded to
native integers / BigInt without a floating-point conversion.

```json
{
  "module": "arith",
  "cases": [
    {
      "function": "gcd",
      "argGenerators": ["randomBigint(-10000, 10000)", "randomBigint(-10000, 10000)"]
    },
    {
      "function": "xgcd",
      "rows": [[42, 0, 0]]
    }
  ]
}
```

Supported generator specifications:

| Specification | Domain |
| --- | --- |
| `randomBigint(min, max)` | Inclusive arbitrary-precision integer interval |
| `randomPrime(min, max)` | Primes in an interval within 0..1,000,000; empty domains fail |
| `randomList(generator, length)` | Fixed-length list generated from the inner specification |
| `fixedValue(value)` | Constant scalar or integer list used alongside random arguments |

Use a small explicit `rows` entry for an entirely fixed regression, with the shape
`[seed, arg1, ...]`. Add `seeds` to a generator only when retaining useful replay
seeds. Constrained cases can retain `replaySeeds: [{"seed":456,"runs":2,
"reason":"bug description"}]` to reproduce the complete generation plan for that
function alongside fresh trials. Keep its recipe and case position stable, or
review/migrate the replay when changing the generator version. These seeds are
regression tests even while a recorded port mismatch remains unfixed.
New passing samples and generated expected results must not be committed.
Document each bug example's purpose in its area module or the audit notes.

```sh
bun tests/property/normalize-cases.ts arith
bun tests/property/normalize-cases.ts --check
```

Area discovery is automatic. Both dispatch tables must have the requested
function. Result formatting belongs in the area module; exception comparisons
must preserve the exact class and complete message. Intentional native-version
adaptations require source evidence and documentation in `DEVIATIONS.md`.

## Failures and regression replay

On a mismatch, the runner prints failing comparisons and a replay command:

```sh
bun tests/property/compare.ts --replay /path/printed/by/the/run/arith.json
```

This temporary failure artifact contains the failing inputs, run seed, generator
version, case-source hash, project revision and SageMath version. It contains no
saved expected outputs. The replay executes both implementations again. Minimize
the counterexample, fix the port and retain the small input or stable generator
seed as a permanent comparative regression. Temporary artifacts must be promoted
into the suite when fixing the bug; they are not permanent tests on their own.

The standalone fast-check properties use fresh default seeds. Their regression
reporter preserves the structured minimized counterexample, seed, shrink path
and fast-check version, and propagates failure. This replaces parsing nested
counterexamples out of human-readable exception text.

## Constrained recipes and native tests

The old multi-million-row sweeps have been retired with the user's explicit
acceptance that unclassified historical bug inputs may be lost. Each area now
stores small constructor domains and generates fresh inputs. This reduces stored
input coverage; it does not establish full behavioral equivalence.

A `recipes` entry is a list of constructor variants. Each variant has one recipe
per argument: `{"integer":[min,max]}`, `{"bytes":32}`,
`{"list":{"length":4,"range":[-8,8]}}`,
`{"scale":{"values":[1,0,0,1],"range":[1,7]}}`,
`{"text":"expression"}` / `{"text":"identifier"}`, or `{"constant":value}`.
List recipes may set `nonzeroLast` to keep a divisor’s leading coefficient nonzero.
Ranges are inclusive; large integers use decimal strings. Positive whole-matrix
scaling preserves the original shape and Gram definiteness. Constants retain
valid field moduli, constructor tags and packed native command structure. For
finite domains, runs select among the declared variants; they do not invent
unsupported constructor combinations. Expand domains by reading the original
wrapper's argument contract, not by blindly randomizing every numeric slot.

97 bulk native output files have been removed. Package tests call
`native-live.mjs`, which generates inputs using `native-suites.json`, executes the
original implementation and checks exact results. Former snapshot filenames are
stable suite IDs only. Named handwritten checks use retained regression inputs;
the two former positional consumers now resolve stable regression IDs.

```sh
# A fresh OS seed is printed by default.
bun test packages/ntl-ts/src/ZZ_random_stream.test.ts
# Replay the same native generated inputs:
SAGEMATH_TEST_SEED=123 SAGEMATH_TEST_RUNS=8 bun test packages/ntl-ts/src/ZZ_random_stream.test.ts
```

Native runs default to 8 trials per profile, plus all constructor variants and
named regressions. Replay requires the same generator version, definitions and
run count. Original/test subprocesses have a 120-second default deadline, settable
with `SAGEMATH_TEST_TIMEOUT_MS`; timeout/crash inputs are saved for area replay.

Large integer regressions can use `shifted` list recipes (`coefficient`, `shift`,
`offset`) to reconstruct exact powers and offsets without decimal blobs.

Small standalone resource-limit controls remain separate from the bulk replay
corpora. `tests/audit/storage-manifest.json` documents the earlier gzip migration;
its filenames and hashes are historical records, not a list of current fixtures.

Use temporary directories for exploratory output. Keep audit handoffs in TODO.md
and `tests/audit/notes/`; do not accumulate successful transcripts, duplicate input
corpora, or generated LCOV reports in the repository.
