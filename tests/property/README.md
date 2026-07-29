# Property Testing Framework

This directory contains cross-language property tests that ensure our TypeScript implementations produce **identical outputs** to their reference implementations:

1. **sagemath-ts vs SageMath** - Python/SageMath comparison
2. **parigp-ts vs PARI/GP** - GP script comparison

## Architecture

```
tests/property/
├── cases/                          # Test case definitions (JSON), one file per area
│   ├── arith.cases.json
│   ├── finite_fields.cases.json
│   └── ...
├── python/                         # Python/SageMath test runner
│   ├── runner.py                   # Area-agnostic driver (discovers areas/)
│   └── areas/                      # One module per area
│       ├── __init__.py
│       ├── _helpers.py             # Helpers shared by >1 area (not an area)
│       ├── arith.py
│       └── ...
├── typescript/                     # TypeScript test runner
│   ├── runner.ts                   # Area-agnostic driver (discovers areas/)
│   └── areas/                      # One module per area
│       ├── _helpers.ts             # Helpers shared by >1 area (not an area)
│       ├── arith.ts
│       └── ...
├── pari/                           # PARI/GP scripts (for parigp-ts)
│   └── elliptic.gp
├── transcripts/                    # Generated output (gitignored)
│   ├── python/                     # SageMath results
│   └── typescript/                 # TypeScript results
├── compare.ts                      # Main comparison harness
├── run-comparison.sh               # PARI/GP comparison script
└── README.md                       # This file
```

### Why areas are separate files

`runner.py` and `runner.ts` contain **no area-specific code**. They glob their
`areas/` directory at run time and load the one area a suite asks for. That means:

- Adding an area is **three new files and zero edits to shared files**, so many
  agents can add areas in parallel without merge conflicts.
- Areas are imported **lazily and one at a time**, so a broken or half-finished
  area module can never break anybody else's area. (Verified: corrupting one area
  module leaves every other area's transcript byte-identical.)

An "area" is just a namespace for property tests. Its name is used in three
places and **must be spelled identically in all three**:

| | |
|---|---|
| `cases/<area>.cases.json` | the `"module"` field inside it must also be `<area>` |
| `python/areas/<area>.py` | exports `FUNCTIONS` |
| `typescript/areas/<area>.ts` | exports `functions` |

## How It Works

### SageMath vs sagemath-ts Comparison

```
┌─────────────────────┐    ┌─────────────────────┐
│   Test Cases JSON   │    │   Test Cases JSON   │
│  (cases/*.json)     │    │  (cases/*.json)     │
│         ↓           │    │         ↓           │
│  python/runner.py   │    │  typescript/runner.ts│
│  (runs in SageMath) │    │  (runs in Bun)      │
│         ↓           │    │         ↓           │
│  python/areas/<a>.py│    │  typescript/areas/<a>.ts │
│         ↓           │    │         ↓           │
│  JSON results       │    │  JSON results       │
└─────────────────────┘    └─────────────────────┘
            ↓                       ↓
            └───────────┬───────────┘
                        ↓
              ┌─────────────────┐
              │  compare.ts     │
              │  (diff check)   │
              │       ↓         │
              │  PASS or FAIL   │
              └─────────────────┘
```

---

## Adding a New Area

Three new files. **Do not edit `runner.py`, `runner.ts`, `compare.ts`, or any
existing area** — that is the whole point of this layout.

Pick an area name matching the SageMath module you are porting (e.g.
`groups_generic`, `padics`, `quadratic_forms`). Check `cases/` first so you do
not collide with an area another agent is adding.

### 1. `cases/<area>.cases.json`

```json
{
  "module": "quadratic_forms",
  "cases": [
    {
      "function": "qf_discriminant",
      "seeds": [42, 123, 999],
      "argGenerators": ["randomBigint(1, 1000)", "randomBigint(1, 1000)"]
    },
    {
      "function": "qf_is_definite",
      "seeds": [1],
      "argGenerators": ["fixedValue(17)"]
    }
  ]
}
```

The top-level `"module"` **must** equal the file's base name, and both must equal
the area module base names below.

### 2. `python/areas/<area>.py` — the SageMath oracle

```python
"""SageMath side of the ``quadratic_forms`` property-test area.

Cases: tests/property/cases/quadratic_forms.cases.json
"""

from sage.all import *


def qf_discriminant(a, b):
    Q = QuadraticForm(ZZ, 2, [a, b, 1])
    return Q.disc()


FUNCTIONS = {
    'qf_discriminant': qf_discriminant,
}
```

Required export: `FUNCTIONS`, a dict of `name -> callable`. The callable receives
the generated arguments positionally.

Optional export: `FORMATTERS`, a dict of `name -> (result) -> str`, for results
the generic `format_result()` in `runner.py` cannot render (see
[Result Formatting](#result-formatting-rules)).

### 3. `typescript/areas/<area>.ts` — the port under test

```typescript
/**
 * sagemath-ts side of the `quadratic_forms` property-test area.
 *
 * Cases: tests/property/cases/quadratic_forms.cases.json
 * SageMath counterpart: tests/property/python/areas/quadratic_forms.py
 */

import { QuadraticForm } from '../../../../packages/sagemath-ts/src/quadratic_forms/index.js';

export const functions = {
  qf_discriminant: (a: bigint, b: bigint) => new QuadraticForm(2, [a, b, 1n]).disc(),
};
```

Required export: `functions`, an object of `name -> callable`, with **exactly the
same keys** as the Python `FUNCTIONS`.

Optional export: `formatters`, an object of `name -> (result: unknown) => string`.

Note the import depth: area modules are one level deeper than the old runner, so
package imports start with `../../../../packages/sagemath-ts/src/...`.

### 4. Run it

```bash
bun run test:property -- --case quadratic_forms   # your area only
bun run test:property                             # everything
```

### Rules for area authors

- **Never edit another area's module.** If you need a helper that already exists
  in someone else's area, copy it or lift it into `_helpers.*` — but lifting into
  `_helpers.*` touches a shared file, so prefer keeping helpers area-private.
- **Files starting with `_` are not areas.** Discovery skips them; that is how
  `_helpers.py` / `_helpers.ts` and `__init__.py` stay out of the area list.
- **Do not add special cases to `formatResult` / `format_result`.** Those live in
  the runners and would make them a merge hotspot. Either return an
  already-formatted string from your area function, or export a `formatters` /
  `FORMATTERS` entry.
- **The oracle is SageMath, not your intuition.** Every expected value must come
  from actually running `sage`. Never "fix" a mismatch by changing the expected
  value to what the port produces.
- A case whose function is missing from **both** area modules produces
  `Unknown function: <area>.<name>` on both sides, and `compare.ts` scores
  "both errored" as a **pass**. Adding a case JSON without both area modules
  therefore gives you a vacuously green area — see
  [Vacuous passes](#vacuous-passes-both-sides-errored) below.

---

## Test Case Format

Test cases are defined in JSON files in the `cases/` directory:

```json
{
  "module": "arith",
  "cases": [
    {
      "function": "gcd",
      "seeds": [42, 123, 999],
      "argGenerators": ["randomBigint(1, 1000)", "randomBigint(1, 1000)"]
    }
  ]
}
```

### Field Definitions

| Field | Description |
|-------|-------------|
| `module` | Area name; must match `cases/<area>.cases.json`, `python/areas/<area>.py`, `typescript/areas/<area>.ts` |
| `function` | Key in that area's `FUNCTIONS` / `functions` dispatch table |
| `seeds` | Array of random seeds; the case runs once per seed (so N seeds = N tests) |
| `argGenerators` | Argument generator specifications, one per positional argument |

### Argument Generators

| Generator | Description | Example |
|-----------|-------------|---------|
| `randomBigint(min, max)` | Random integer in range [min, max] | `randomBigint(1, 1000)` |
| `randomPrime(min, max)` | Random prime in range [min, max] | `randomPrime(2, 1000)` |
| `fixedValue(value)` | Fixed/constant value; also `fixedValue([1, 2, 3])` for a list | `fixedValue(17)` |
| `randomList(gen, len)` | List of random values | `randomList(randomBigint(1, 100), 5)` |

Both runners drive these from the **same Mersenne Twister**, seeded identically,
so `randomBigint`/`randomPrime` produce the same values in Python and TypeScript.

## Running Tests

### Run All Comparison Tests

```bash
bun run test:property
```

This runs both SageMath and TypeScript with the same test cases and compares results.

### Run Specific Area

```bash
bun run test:property -- --case arith
```

### Generate SageMath Results Only

```bash
bun run test:property:generate
```

Useful for generating expected results when SageMath is available.

### Run TypeScript Tests Only

```bash
bun run test:property:typescript
```

Compares TypeScript results against saved SageMath transcripts (from `transcripts/python/`).

### Verbose Output

```bash
bun run test:property:verbose
```

Shows all test results, not just failures.

### PARI/GP Comparison (for parigp-ts)

```bash
bun run test:property:pari           # Run all PARI tests
bun run test:property:pari:generate  # Generate PARI output only
bun run test:property:pari:compare   # Run TypeScript comparison only
```

#### Elliptic Curve Test Categories (75 tests)

The elliptic curve property tests cover:

| Category | Tests | Functions Tested |
|----------|-------|------------------|
| ellinit | 3 | `ellinit` with various Weierstrass forms |
| ellcard | 8 | Cardinality computation (various primes) |
| Point operations | 7 | `elladd`, `ellmul`, `ellneg` |
| ellorder | 3 | Point order computation |
| Identity tests | 5 | P+O=P, P+(-P)=O, etc. |
| Curve invariants | 2 | Discriminant, j-invariant |
| ellisoncurve | 2 | Valid/invalid point detection |
| ellgroup | 3 | Group structure [d1] or [d1,d2] |
| Associativity | 2 | (P+Q)+R = P+(Q+R), -nP = -(nP) |
| elllog | 5+ | Discrete logarithm (Pohlig-Hellman) |
| elldivpol | 7 | Division polynomials psi_n |
| elltatepairing | 4 | Tate pairing |
| ellweilpairing | 3 | Weil pairing |
| Known curves | 3 | secp256k1-like, P-256-like |
| Random curves | 3 | Deterministic random parameters |
| Edge cases | 4 | Small primes p=5,7,11, singular curves |
| Large fields | 2 | F_4999 |
| Bilinearity | 2 | e(2P,Q) = e(P,Q)^2 |

## Output Format

Both runners produce JSON output with this structure:

```json
[
  {
    "function": "gcd",
    "args": ["12", "8"],
    "result": "4",
    "error": null,
    "seed": 42
  }
]
```

### Result Formatting Rules

Results must match **exactly** between SageMath and TypeScript:

| Type | Python Format | TypeScript Format |
|------|---------------|-------------------|
| Integer | `str(n)` | `n.toString()` |
| Boolean | `True`/`False` | `'True'`/`'False'` |
| List | `[1, 2, 3]` | `[1, 2, 3]` |
| Tuple | `(a, b, c)` | `(a, b, c)` |
| Factorization | `2^2 * 3` | `2^2 * 3` |
| `None` / `null` | `null` | `null` |

If your area's natural return value does not land on the same string on both
sides, the cheapest fix is to **format it yourself inside the area function** and
return a plain string (this is what `elliptic_curves` does for points). The next
cheapest is a `FORMATTERS` / `formatters` entry. Do not touch the runners.

### Vacuous passes (both sides errored)

`compare.ts` scores a case as **passed** whenever *both* runners raised — it only
checks that the two implementations agree on failing, not on *why*. This is
intentional (it lets a case pin "SageMath also rejects this input"), but it means
a case whose area module is missing on both sides passes without testing
anything.

Three areas are in exactly that state today — `cases/*.cases.json` exists but
neither `python/areas/` nor `typescript/areas/` has a module, so all their cases
error with `Unknown module: <area>` on both sides and are scored as passes:

| Area | Cases scored as vacuous passes |
|------|-------------------------------|
| `arith_extended` | 26 |
| `lwe` | 15 |
| `matrix` | 16 |

That is 57 of the 433 currently-green tests. Whoever implements those areas
should expect the number of *real* assertions to go up, not the pass count.
To check whether an area is real, look at its transcript:

```bash
bun run test:property -- --case lwe
python3 -c "import json;d=json.load(open('tests/property/transcripts/python/lwe.json'));print(set(r['error'] for r in d))"
```

## Debugging Failures

When tests fail:

1. **Check the diff** - Compare `transcripts/python/<area>.json` with `transcripts/typescript/<area>.json`

2. **Verify result formatting** - Often mismatches are due to string representation differences

3. **Check for edge cases** - SageMath may handle edge cases differently

4. **Run with verbose mode**:
   ```bash
   bun run test:property:verbose -- --case arith
   ```

5. **Run one runner by hand** - both read a cases JSON on stdin and write results JSON on stdout:
   ```bash
   sage tests/property/python/runner.py      < tests/property/cases/arith.cases.json
   bun run tests/property/typescript/runner.ts < tests/property/cases/arith.cases.json
   ```

### Manual Verification

You can manually verify results in SageMath:

```python
sage: gcd(12, 8)
4
sage: factor(100)
2^2 * 5^2
sage: is_prime(17)
True
```

And in TypeScript:

```typescript
import { gcd, factor, is_prime } from '@sagemath-ts/sagemath-ts/arith';

console.log(gcd(12n, 8n));        // 4n
console.log(factor(100n));        // [[2n, 2n], [5n, 2n]]
console.log(is_prime(17n));       // true
```

## Requirements

- **SageMath** - For Python/SageMath comparison tests
  - macOS: `brew install sage`
  - Ubuntu: `sudo apt-get install sagemath`

- **PARI/GP** - For parigp-ts comparison tests
  - macOS: `brew install pari`
  - Ubuntu: `sudo apt-get install pari-gp`

- **Bun** - For running TypeScript tests
  - https://bun.sh/

## CI Integration

In CI environments without SageMath, use `--typescript-only` mode with pre-generated transcripts:

```yaml
steps:
  - name: Run property tests
    run: bun run test:property:typescript
```

Store generated transcripts in the repository for CI comparison.
