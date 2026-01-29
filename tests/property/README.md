# Property Testing Framework

This directory contains cross-language property tests that ensure our TypeScript implementations produce **identical outputs** to their reference implementations:

1. **sagemath-ts vs SageMath** - Python/SageMath comparison
2. **parigp-ts vs PARI/GP** - GP script comparison

## Architecture

```
tests/property/
├── cases/                      # Test case definitions (JSON)
│   ├── arith.cases.json        # Arithmetic function test cases
│   ├── finite_fields.cases.json
│   ├── polynomials.cases.json
│   └── elliptic_curves.cases.json
├── python/                     # Python/SageMath test runner
│   └── runner.py               # Executes tests in SageMath
├── typescript/                 # TypeScript test runner
│   └── runner.ts               # Executes tests in sagemath-ts
├── pari/                       # PARI/GP scripts (for parigp-ts)
│   └── elliptic.gp
├── transcripts/                # Generated output (gitignored)
│   ├── python/                 # SageMath results
│   └── typescript/             # TypeScript results
├── compare.ts                  # Main comparison harness
├── run-comparison.sh           # PARI/GP comparison script
└── README.md                   # This file
```

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
    },
    {
      "function": "is_prime",
      "seeds": [1],
      "argGenerators": ["fixedValue(17)"]
    }
  ]
}
```

### Field Definitions

| Field | Description |
|-------|-------------|
| `module` | Module name (arith, finite_fields, polynomials, elliptic_curves) |
| `function` | Function name to test |
| `seeds` | Array of random seeds for reproducibility |
| `argGenerators` | Argument generator specifications |

### Argument Generators

| Generator | Description | Example |
|-----------|-------------|---------|
| `randomBigint(min, max)` | Random integer in range [min, max] | `randomBigint(1, 1000)` |
| `randomPrime(min, max)` | Random prime in range [min, max] | `randomPrime(2, 1000)` |
| `fixedValue(value)` | Fixed/constant value | `fixedValue(17)` |
| `randomList(gen, len)` | List of random values | `randomList(randomBigint(1, 100), 5)` |

## Running Tests

### Run All Comparison Tests

```bash
bun run test:property
```

This runs both SageMath and TypeScript with the same test cases and compares results.

### Run Specific Module

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

## Adding New Tests

### 1. Add Test Cases

Edit or create a file in `cases/`:

```json
{
  "module": "arith",
  "cases": [
    {
      "function": "new_function",
      "seeds": [42, 123],
      "argGenerators": ["randomBigint(1, 100)"]
    }
  ]
}
```

### 2. Add Function to Runners

**Python runner** (`python/runner.py`):
```python
function_map = {
    'arith': {
        # ... existing functions ...
        'new_function': new_function,
    },
}
```

**TypeScript runner** (`typescript/runner.ts`):
```typescript
const functionMap = {
  arith: {
    // ... existing functions ...
    new_function: (n: bigint) => new_function(n),
  },
};
```

### 3. Run Tests

```bash
bun run test:property -- --case arith
```

## Seeded Random Number Generation

### Matching RNG Between Languages

For deterministic tests, both implementations use seeded random number generators.

**Python (SageMath):**
```python
set_random_seed(SEED)
ZZ.random_element(lower, upper)
```

**TypeScript:**
```typescript
const random = new SeededRandom(SEED);
random.integer(lower, upper);
```

**Note:** The current TypeScript SeededRandom uses a simple LCG algorithm. For exact RNG compatibility with SageMath's Mersenne Twister, additional work may be needed. For now, tests use `fixedValue()` generators where exact reproducibility is critical.

## Debugging Failures

When tests fail:

1. **Check the diff** - Compare `transcripts/python/<module>.json` with `transcripts/typescript/<module>.json`

2. **Verify result formatting** - Often mismatches are due to string representation differences

3. **Check for edge cases** - SageMath may handle edge cases differently

4. **Run with verbose mode**:
   ```bash
   bun run test:property:verbose -- --case arith
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
