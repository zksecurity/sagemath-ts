# sagemath-ts LLM Reference

TypeScript port of SageMath's number-theory and cryptography modules. Pure TypeScript,
no native dependencies, ~3.5 MB.

**Warning**: Educational only. Not constant-time, not audited. Do not use for production
cryptography.

Every example below is executed by `tests/llm-doc.test.ts`. If an example here disagrees
with the code, the test fails — trust this file over your memory of the API, and run
`bun test tests/llm-doc.test.ts` if you suspect drift.

---

## Getting the library

**sagemath-ts is not published to npm.** `npm install sagemath-ts` and
`bun add sagemath-ts` both 404. The only way to get it is a git checkout:

```bash
git clone https://github.com/ZkSecurity/sagemath-ts
cd sagemath-ts
bun install          # required: creates the node_modules links the packages resolve through
```

`./scripts/clone-references.sh` is only needed to work *on* the port (it fetches the
SageMath/PARI/FLINT sources for reference). It is not needed to *use* the library.

## Importing

Which specifier works depends on where your file lives. Pick the row that matches you:

| Your file lives... | Import as | Notes |
|---|---|---|
| Anywhere inside the checkout | `from 'sagemath-ts'` | Works at the repo root, in `packages/*`, in `tests/`, `tutorial/`, `playground/`. **Prefer this.** |
| Inside `packages/sagemath-ts/` | `from '@sagemath-ts/sagemath-ts'` | Also resolves, but only from that one directory. Does **not** resolve from the repo root or from sibling packages. |
| **Outside** the checkout (vendored copy, read-only Docker mount, separate project) | `from '<path-to-checkout>/packages/sagemath-ts/src/index.ts'` | Bare specifiers cannot resolve from outside the tree; a direct path to the entry point can. Use an absolute path, or a relative one from your file. |

`@sagemath-ts/sagemath-ts` is the real `name` in `packages/sagemath-ts/package.json`;
`sagemath-ts` is the repo-root workspace alias. Both point at the same code. The bare
alias is the one that works from every location, so reach for it first.

Out-of-tree example — a script at `/workspace/run.ts` against a checkout mounted
read-only at `/tools/sagemath-ts`:

```typescript
import { factor, GF, EllipticCurve } from '/tools/sagemath-ts/packages/sagemath-ts/src/index.ts';
```

The mount may be read-only; the library never writes to its own tree. It does need
`node_modules` to already exist inside the checkout (from `bun install`), because
`sagemath-ts` depends on the sibling `parigp-ts`, `flint-ts` and `ntl-ts` packages.

Subpath entry points also work, from any of the locations above:

```typescript
import { gcd } from 'sagemath-ts/arith';
import { GF } from 'sagemath-ts/rings/finite_rings';
```

Available subpaths: `arith`, `crypto`, `crypto/lwe`, `rings`, `rings/finite_rings`,
`rings/polynomial`, `rings/function_field`, `matrix`, `algebras`, `algebras/quatalg`,
`quadratic_forms`, `schemes/elliptic_curves`, `schemes/hyperelliptic_curves`,
`misc/randstate`, `stats`, `stats/distributions`.

### The package root is a curated subset — most of the library is not on it

`import { … } from 'sagemath-ts'` exposes **159 names**. The subpaths expose far more,
and a missing name fails loudly at import time (`SyntaxError: Export named 'nth_prime'
not found`). If a function you expect is not at the root, it is almost certainly on its
subpath:

| Subpath | Exports | Not re-exported at the root |
|---|---:|---:|
| `sagemath-ts/rings` | 257 | 241 |
| `sagemath-ts/matrix` | 239 | 226 |
| `sagemath-ts/arith` | 100 | 58 |
| `sagemath-ts/schemes/elliptic_curves` | 99 | 86 |
| `sagemath-ts/rings/polynomial` | 82 | **all 82** |
| `sagemath-ts/rings/finite_rings` | 43 | 35 |
| `sagemath-ts/crypto` | 35 | 22 |
| `sagemath-ts/stats` | 13 | 7 |

Common names that are **not** at the root:

```typescript
import { nth_prime, binomial, factorial, fibonacci, kronecker, hilbert_symbol,
         bernoulli, primitive_root, CRT, trial_division, random_prime } from 'sagemath-ts/arith';
import { PolynomialRing, Polynomial, NumberField } from 'sagemath-ts/rings';
```

The full list of `arith` names missing from the root: `CRT`, `algdep`,
`algebraic_dependency`, `bernoulli`, `binomial`, `binomial_coefficients`,
`carmichael_lambda`, `continuant`, `coprime_part`, `dedekind_psi`, `dedekind_sum`,
`differences`, `eratosthenes`, `factorial`, `falling_factorial`, `fibonacci`,
`four_squares`, `fundamental_discriminant`, `gauss_sum`, `get_gcd`, `get_inverse_mod`,
`hilbert_conductor`, `hilbert_conductor_inverse`, `hilbert_symbol`, `integer_ceil`,
`integer_floor`, `integer_trunc`, `is_power_of_two`, `is_pseudoprime`,
`is_pseudoprime_power`, `is_strong_probable_prime`, `kronecker`, `lucas_number`,
`mqrr_rational_reconstruction`, `multinomial`, `multinomial_coefficients`,
`next_prime_power`, `next_probable_prime`, `nth_prime`, `odd_part`,
`previous_prime_power`, `prime_powers`, `prime_to_m_part`, `primes`, `primes_first_n`,
`primitive_root`, `quadratic_residues`, `random_prime`, `rising_factorial`,
`smooth_part`, `sort_complex_numbers_for_display`, `squarefree_divisors`, `subfactorial`,
`sum_of_k_squares`, `three_squares`, `trial_division`, `two_squares`, `xlcm`.

When in doubt, import from the subpath — everything at the root is also on its subpath,
so the subpath import always works.

## Runtime

The library ships as TypeScript source, not compiled JavaScript. Bun (`bun run x.ts`)
runs it directly. Node needs `>=22` plus a TypeScript loader, or run
`bun run build` in `packages/sagemath-ts` first.

---

## Four rules that cause most first-attempt failures

**1. JavaScript `number` is rejected for integer arguments.** `IntegerLike` is
`bigint | Integer` — deliberately *not* `number`, because IEEE-754 silently loses
precision past 2^53. Use bigint literals everywhere.

```typescript
gcd(12n, 8n)   // 4n
gcd(12, 8)     // TypeError: JavaScript numbers are not accepted due to precision
               // loss risk; use bigint literals (e.g., 123n)
```

This rule covers the arbitrary-precision free functions. Three places *do* take `number`:

- genuinely non-integer parameters (`sigma`, `delta`, `eta`) and array indices / matrix
  dimensions (`nrows`, `ncols`, `get(i, j)`);
- ring constructors — `GF(7)`, `Zmod(7)`, `Mod(10, 7)`;
- ring *element* arithmetic — `F.__call__(3n).add(4)` works, because the value is already
  reduced into a bounded ring.

Note also that `IntegerLike` is not applied uniformly: about 37 exported `arith` functions
declare a bare `bigint` parameter instead (`nth_prime`, `primitive_root`, `random_prime`,
`two_squares`, `eratosthenes`, `next_prime_power`, …). They never call `toBigInt`, so
passing an `Integer` is a type error that can also misbehave at runtime —
`nth_prime(new Integer(5n))` throws `ValueError: nth prime not found` rather than
returning `11n`. **Pass bigints to those, not `Integer` wrappers.**

**2. Strings are not coerced.** `toBigInt('123')` throws. Use `BigInt('123')` first.

**3. Rings are not callable — use `.__call__(x)`.** TypeScript classes cannot be made
callable, so SageMath's `F(3)` becomes `F.__call__(3n)`.

```typescript
const F = GF(13n);
F(5n)              // TypeError: F is not a function
F.__call__(5n)     // correct
```

**4. Method names follow SageMath, not JS shorthand.** `determinant()` not `det()`;
`.call()` not `.sample()` on LWE oracles.

---

## Arithmetic

```typescript
import { gcd, lcm, xgcd, factor, is_prime, power_mod, inverse_mod, euler_phi, crt, CRT_list } from 'sagemath-ts';

gcd(12n, 8n)                     // 4n
lcm(12n, 8n)                     // 24n
xgcd(15n, 6n)                    // [3n, 1n, -2n]  — an array [g, s, t] with s*15 + t*6 = g
factor(60n)                      // [[2n, 2n], [3n, 1n], [5n, 1n]]  — both entries are bigint
is_prime(97n)                    // true
power_mod(2n, 100n, 1000000007n) // 976371285n
inverse_mod(3n, 7n)              // 5n
euler_phi(60n)                   // 16n
crt(2n, 3n, 5n, 7n)              // 17n  — four scalars, not two arrays
CRT_list([2n, 3n], [5n, 7n])     // 17n  — the list form
```

## Primes

```typescript
import { next_prime, previous_prime, prime_range, prime_factors, is_prime_power } from 'sagemath-ts';

next_prime(10n)                  // 11n
previous_prime(10n)              // 7n
prime_range(10n, 30n)            // [11n, 13n, 17n, 19n, 23n, 29n]
prime_factors(60n)               // [2n, 3n, 5n]
is_prime_power(8n)               // true          — a boolean by default
is_prime_power(8n, true)         // [2n, 3n]      — pass get_data for [base, exponent]
```

## Divisibility

```typescript
import { divisors, sigma, valuation, is_square, isqrt, squarefree_part } from 'sagemath-ts';

divisors(12n)                    // [1n, 2n, 3n, 4n, 6n, 12n]
sigma(12n, 0n)                   // 6n   — number of divisors (note the bigint 0n)
sigma(12n, 1n)                   // 28n  — sum of divisors
valuation(12n, 2n)               // 2n
is_square(16n)                   // true
is_square(15n, true)             // [false, null]
isqrt(17n)                       // 4n
squarefree_part(12n)             // 3n
```

## Factorization: what it actually does

`factor()` delegates to `@sagemath-ts/parigp-ts`, a port of PARI's `ifactor1.c`. It is
**not** trial division. The full PARI cascade is implemented: trial division, then
Shanks' SQUFOF, Pollard–Brent rho, Lenstra–Montgomery ECM, and MPQS.

Measured wall-clock on a balanced semiprime (Apple Silicon, Bun 1.3):

| Input size | Time |
|---|---|
| 64-bit | ~8 ms |
| 128-bit | ~56 ms |
| 144-bit | ~0.6 s |
| 160-bit | ~1.5 s |
| 176-bit | ~3.7 s |
| 200-bit | ~32 s |

So: comfortable through ~160 bits, painful past ~200, hopeless at RSA sizes — the same
shape as PARI itself, just slower by the JS constant factor. Reach for an external tool
only above that range, not for ordinary composites.

```typescript
factor(12345678901234567890n)
// [[2n,1n],[3n,2n],[5n,1n],[101n,1n],[3541n,1n],[3607n,1n],[3803n,1n],[27961n,1n]]
```

`factor(0n)` throws `ArithmeticError`, matching SageMath.

## Rings and fields

```typescript
import { ZZ, QQ, Integer, Rational, Zmod, GF, Mod } from 'sagemath-ts';

// Integers. ZZ is a singleton ring object, not a constructor.
ZZ.__call__(42n)                 // 42n (a bigint — ZZ's elements are plain bigints)
new Integer(42n).factor()        // [[2n,1n],[3n,1n],[7n,1n]]  — the wrapper class

// Rationals
new Rational(3n, 4n).add(new Rational(1n, 2n))   // 5/4

// Modular integers (Z/nZ)
const Z7 = Zmod(7n);
Z7.__call__(3n).mul(Z7.__call__(5n))             // 1  (15 mod 7)
Mod(10n, 7n)                                     // 3  — shorthand for a single element

// Prime fields
const F13 = GF(13n);
F13.__call__(5n).pow(12n)                        // 1  (Fermat)
```

`GF` returns a `FiniteFieldPrime`. Its methods: `__call__`, `zero`, `one`, `gen`,
`cardinality`, `elements`, `is_field`, `random_element`, `multiplicative_generator`,
`primitive_element`, `primitiveRoot`.

## Matrices

There are two families. **Use `IntegerMatrix` for integer linear algebra** — the generic
`Matrix<R>` needs entries that are ring-element objects, and `ZZ`'s elements are raw
bigints without `.mul`, so `matrix(ZZ, ...)` builds but cannot multiply. Over a field it
is fine: `matrix(GF(7n), [[1n, 2n], [3n, 4n]])` multiplies correctly, and the factory
coerces raw bigints for you.

```typescript
import { IntegerMatrixFromEntries, identity_integer_matrix, zero_integer_matrix } from 'sagemath-ts';

const A = IntegerMatrixFromEntries([[1n, 2n], [3n, 4n]]);
const B = IntegerMatrixFromEntries([[5n, 6n], [7n, 8n]]);

A.add(B)
A.mul(B)                         // [[19, 22], [43, 50]]
A.determinant()                  // -2  (an Integer — note: determinant(), not det())
A.transpose()
A.rank()                         // 2
A.hermite_form()
A.smith_form()
A.elementary_divisors()
A.right_kernel_matrix()
identity_integer_matrix(3)
zero_integer_matrix(2, 3)
```

`IntegerMatrix` has **no** `inverse()`. For LLL and lattice work see below.

## Elliptic curves

```typescript
import { EllipticCurve, GF } from 'sagemath-ts';

const F = GF(101n);
const E = EllipticCurve(F, [1n, 2n]);   // y^2 = x^3 + x + 2

E.order()                        // 100n
E.discriminant()                 // 26
E.j_invariant()                  // 4

const P = E.point(1n, 2n);       // bigints accepted directly over a finite field
P.order()                        // 4n
P.add(P)
P.mul(4n)                        // (0 : 1 : 0), the point at infinity

E.is_on_curve(F.__call__(1n), F.__call__(2n))   // true — takes two field elements, not an array
E.random_point()
E.lift_x(F.__call__(4n))         // throws ValueError when no point has that x
```

Point methods: `add`, `sub`, `neg`, `double`, `mul`, `order`, `has_order`, `isZero`,
`weil_pairing`, `tate_pairing`, `ate_pairing`.

Curve order and point order delegate to PARI (`ellcard`/`ellorder`), so they use SEA and
BSGS rather than enumeration.

## Pairings

```typescript
import { weil_pairing, tate_pairing, ate_pairing, embedding_degree } from 'sagemath-ts';

weil_pairing(P, Q, n)            // n: bigint
tate_pairing(P, Q, n, k)         // k: embedding degree
ate_pairing(P, Q, n, k, t)       // t: trace of Frobenius
embedding_degree(E, 5n)          // smallest k with n | q^k - 1
```

## Lattices and LLL

Both take an `IntegerMatrix` or a plain `bigint[][]` — **not** an array of `vector()`s.

```typescript
import { IntegerLattice, lllReduce, IntegerMatrixFromEntries } from 'sagemath-ts';

const L = IntegerLattice([[1n, 0n, 0n], [0n, 1n, 0n], [0n, 0n, 1n]]);

lllReduce(IntegerMatrixFromEntries([[1n, 2n, 3n], [4n, 5n, 6n], [7n, 8n, 10n]]))
// [[0, 0, 1], [-1, 1, 0], [2, 1, 0]]
```

`lllReduce(basis, { delta, eta })` defaults to `delta = 0.99`, `eta = 0.501`.

## LWE

These mirror `sage.crypto.lwe`: they are **oracles**, sampled with `.call()`. There is no
`keygen`/`encrypt`/`decrypt` — `Regev` here is Sage's LWE parameter set, not a full
encryption scheme.

```typescript
import { LWE, Regev, LindnerPeikert, RingLWE, DiscreteGaussianInteger } from 'sagemath-ts';

const D = DiscreteGaussianInteger(3.2);   // factory function; sigma is a plain number
D.call()                                  // one integer sample

const lwe = new LWE(32n, 40961n, D);      // positional: (n, q, D, secret_dist?, m?)
const [a, c] = lwe.call();                // a: LWEVector, c: IntegerMod
lwe.samples(10n)                          // ten [a, c] pairs

const regev = new Regev(32n);             // positional: (n, secret_dist?, m?)
const [a2, c2] = regev.call();
```

`secret_dist` is `'uniform'` (default), `'noise'`, or `[lb, ub]`.

## Discrete Gaussian

```typescript
import { DiscreteGaussianInteger, DiscreteGaussianDistributionIntegerSampler } from 'sagemath-ts';

DiscreteGaussianInteger(3.2).call()                             // sigma = 3.2, centered at 0
DiscreteGaussianInteger(3.2, 5n, 6n)                            // (sigma, c, tau, algorithm?)
new DiscreteGaussianDistributionIntegerSampler({ sigma: 3.2 })  // the class takes an options object
```

## Coding theory

`length` must divide `q - 1`.

```typescript
import { createClassicalReedSolomonCode, GF } from 'sagemath-ts';

const F = GF(929n);
const RS = createClassicalReedSolomonCode(F, 928n, 4n);   // n and k are IntegerLike (bigint)
RS.minimum_distance()            // 925 = n - k + 1
```

## Groups

```typescript
import { bsgs, pohlig_hellman, order_from_multiple, GF } from 'sagemath-ts';

const F = GF(101n);
const g = F.__call__(2n);
bsgs(g, g.pow(37n), [0n, 100n])  // 37n
```

`bsgs(a, b, bounds, operation?, identity?, inverse?, op?)` is **positional**, not an
options object. `operation` defaults to `'*'`; use `'+'` for additive groups. For
`'other'` you must supply all three of `identity`, `inverse` and `op`.

## Type coercion

```typescript
import { toBigInt, toRational, type IntegerLike, type RationalLike } from 'sagemath-ts';

toBigInt(42n)             // 42n
toBigInt(new Integer(42n))// 42n
toBigInt('123')           // TypeError — convert with BigInt('123') first
toBigInt(42)              // TypeError — numbers are rejected on purpose
toRational(3n)            // 3  — takes bigint | Integer | Rational, not '3/4'
```

`toSafeNumber` (bigint to number, `RangeError` past 2^53) exists in
`src/types/coercion.ts` but is **not** re-exported from the package root, and the
`exports` map blocks unlisted subpaths — reach it by file path, or just use
`Number(x)` where the precision loss is intentional.

## Namespaced imports

```typescript
import { arith, crypto, coding, groups, modules, schemes, stats } from 'sagemath-ts';

arith.gcd(12n, 8n)
crypto.LWE
coding.ReedSolomonCode
groups.bsgs
modules.lllReduce
schemes.EllipticCurve
stats.DiscreteGaussianInteger
```

## Related packages

Same checkout, same import rules (`@sagemath-ts/…` from inside `packages/*`, or a direct
path to `packages/<pkg>/src/index.ts` from outside):

- `@sagemath-ts/parigp-ts` — PARI/GP algorithms (factorization, elliptic curves, SEA)
- `@sagemath-ts/flint-ts` — FLINT port
- `@sagemath-ts/ntl-ts` — NTL port
- `@zksecurity/cheatsheets` — named curve parameters

## Further reading

- `README.md` — setup, layout, playground
- `DESIGN.md` — type mappings and architectural conventions
- `DEVIATIONS.md` — every documented behavioral difference from SageMath
- `SCOPE.md` — which modules are implemented
- `tutorial/` — 77 runnable lessons, browsable via `bun playground`
