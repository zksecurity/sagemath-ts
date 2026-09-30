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

`valuation(n: IntegerLike, p: IntegerLike)` returns `bigint | 'Infinity'` and delegates
to Integer. Zero returns `'Infinity'` before base validation. Its generic overload forwards
to an object's `valuation(p)` method, preserving that method's result type; Rational values
are supported. Objects with `_integer_(parent: typeof ZZ): IntegerLike` also work when
no valuation method exists. An `AttributeError` from lookup or invocation triggers Integer coercion.
Rational `val_unit`, `padic_valuation`, `numerator_valuation` and `denominator_valuation`
also accept IntegerLike bases. `val_unit` returns `[bigint | 'Infinity', Rational]`;
component and padic aliases retain their existing return types.

```typescript
import { valuation } from 'sagemath-ts';
valuation(0n, 0n)                              // 'Infinity'
valuation({_integer_: () => 40n}, 2n)          // 3n
valuation(new Rational(5n, 9n), 3n)            // -2n
new Rational(-4n, 17n).val_unit(new Integer(2n)) // [2n, Rational(-1n, 17n)]
new Integer((1n << 262144n) - 1n).nbits()       // 262144n
```

Integer `bit_length` and its `nbits` alias read magnitude-size metadata cached when the
immutable wrapper is constructed; repeated queries do not scan or shift the integer.

`sort_complex_numbers_for_display` accepts an array of `{ re: number; im: number }`
values, or tuple records whose first entry has that shape. It preserves each item and
its metadata, returning a sorted copy for nonempty input and the same array for empty
input. Only an exactly zero imaginary part counts as real. Non-real keys round their
real parts to 34 binary bits for display; the ordering also follows Python for NaNs.
This input representation corresponds to Sage's binary64 complex domain.

```typescript
import { sort_complex_numbers_for_display } from 'sagemath-ts/arith';
const nonreal = [{ re: -1, im: 1e-15 }, 'nonreal'] as const;
const real = [{ re: 2, im: 0 }, 'real'] as const;
sort_complex_numbers_for_display([nonreal, real]) // [real, nonreal], same records
```

`mqrr_rational_reconstruction(u: IntegerLike, m: IntegerLike, T: IntegerLike)` returns
`[bigint, bigint] | null`. It preserves the original Python body's exact rational division,
including cases where that body returns no reconstruction. Its source comment about
implicit floor division is stale; this is distinct from `rational_reconstruction`.

```typescript
import { mqrr_rational_reconstruction } from 'sagemath-ts/arith';
mqrr_rational_reconstruction(new Integer(0n), 100n, 50n) // [0n, 1n]
mqrr_rational_reconstruction(15n, 31n, 2n)              // [15n, 1n]
mqrr_rational_reconstruction(16n, 31n, 2n)              // null
```

Square decompositions accept `IntegerLike`: `two_squares(n)`, `three_squares(n)` and
`four_squares(n)` return bigint tuples of length two, three and four. Impossible
representations raise `ValueError` with Sage's message. `sum_of_k_squares(k, n)` returns
`bigint[]`; `k` accepts `IntegerLike | number`, truncating finite numeric counts as
Python `int` does. Nonfinite numeric counts raise conversion errors. Prefer bigint counts.
The native `two_squares_pyx`, `three_squares_pyx`, `four_squares_pyx` and
`is_sum_of_two_squares_pyx` are exported by `sagemath-ts/rings`; their IntegerLike inputs
must lie in `[0, 2^32)`, otherwise they raise `OverflowError`.

```typescript
import { two_squares, three_squares, sum_of_k_squares } from 'sagemath-ts/arith';
import { is_sum_of_two_squares_pyx } from 'sagemath-ts/rings';
two_squares(new Integer(106n))      // [5n, 9n]
three_squares(new Integer(0n))      // [0n, 0n, 0n]
sum_of_k_squares(2n, 9634n)         // [15n, 97n]
sum_of_k_squares(2.9, 9634n)        // [15n, 97n]
is_sum_of_two_squares_pyx(21n)      // false
two_squares(21n)                   // ValueError: 21 is not a sum of 2 squares
```

Prime iteration follows Sage's range conventions. `primes(start?: IntegerLike,
stop?: IntegerLike | number, proof?: boolean)` accepts `Infinity` as the only numeric
stop; finite bounds use `bigint` or `Integer`. A single argument is the exclusive upper
bound, and no arguments produce an empty iterator. Draws are lazy for finite bounds too.
`differences(lis: IntegerLike[], n?: IntegerLike): bigint[]` requires positive order;
orders at least the list length produce an empty list. `subfactorial(n: IntegerLike)`
returns a `bigint`.

```typescript
import { primes, differences, subfactorial } from 'sagemath-ts/arith';
[...primes(13n)]                         // [2n, 3n, 5n, 7n, 11n]
[...primes()]                            // []
primes(10n, Infinity).next().value        // 11n
differences([1n, 4n, 9n], 3n)            // []
subfactorial(8n)                         // 14833n
```

`dedekind_sum(p: IntegerLike, q: IntegerLike, algorithm?: 'default' | 'flint' | 'pari')`
returns a reduced `{numerator: bigint, denominator: bigint}`. The default delegates to
FLINT, whose result is zero for `q <= 2`; the PARI option preserves PARI's signed-denominator
behavior and raises `PariError` for `q = 0`. Unknown algorithms raise `ValueError`.
The dependencies export `fmpq_dedekind_sum(h: bigint, k: bigint)` (FLINT) and
`sumdedekind(h: bigint, k: bigint)` (PARI), each returning `[bigint, bigint]`.

```typescript
import { dedekind_sum } from 'sagemath-ts/arith';
dedekind_sum(5n, -7n)          // {numerator: 0n, denominator: 1n}
dedekind_sum(5n, -7n, 'pari')  // {numerator: 1n, denominator: 14n}
```

`integer_floor`, `integer_ceil` and `integer_trunc` return `bigint`. They accept numeric
inputs, Integer/Rational/RealNumber objects and prime modular/field elements. Floor and
ceiling also parse strings or byte arrays as Python floats; truncation rejects those before
conversion. Exact numeric objects use their floor/ceiling methods. NaN/infinity errors
therefore depend on whether the input is a JavaScript number or a real-field element.

```typescript
import { integer_floor, integer_ceil, integer_trunc } from 'sagemath-ts/arith';
integer_floor(new Rational(-5n, 2n))  // -3n
integer_ceil('2.5')                  // 3n
integer_trunc(new Rational(-5n, 2n))  // -2n
```

When in doubt, import from the subpath — everything at the root is also on its subpath,
so the subpath import always works.

Two names are bound to *different* implementations depending on where you import them
from, so pick deliberately: `rational_reconstruction` (root vs `sagemath-ts/matrix`) and
`order_from_multiple` (root vs `sagemath-ts/schemes/elliptic_curves`).

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

Some exported `arith` functions still declare bare `bigint` parameters
(`primitive_root`, `random_prime`, …).
Use bigints for those. `nth_prime(n: IntegerLike): bigint` supports Integer wrappers and
delegates to the PARI prime-table/successor implementation.
`primes_first_n(n: IntegerLike | number): bigint[]` keeps its numeric-count API;
positive fractional counts are truncated after rejecting negative inputs.

```typescript
import { nth_prime, primes_first_n } from 'sagemath-ts/arith';
nth_prime(new Integer(5n))       // 11n
nth_prime(10000000n)             // 179424673n
primes_first_n(2.5)              // [2n, 3n]
primes_first_n(0.5)              // []
```

The dependency package also exports `prime(n: bigint): bigint`, corresponding to PARI's
`prime(long)`. Its input follows the signed 64-bit C-long range used by the oracle.

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

`Integer.factorial(): Integer` uses the native GMP factorial algorithm. Negative inputs
raise `ValueError('factorial only defined for nonnegative integers')`; values at least
`2^64` raise `OverflowError('argument too large for factorial')`. Physical BigInt allocation
limits depend on the runtime. Free `factorial(n)` and `factorial(n, 'gmp')` delegate here;
unknown selectors raise ValueError after negative-input validation.

`factorial(n, 'pari')` returns a PARI real record `{s, e, m, p}`: its represented value
is `s * m * 2^(e + 1 - p)`. The sign and precision are numbers; the mantissa and
exponent are bigints. This backend rounds to 64 bits and uses gamma for large inputs.
It does not materialize a huge exact factorial. Inputs at least `2^63` raise
`OverflowError('Python int too large to convert to C long')`; results beyond PARI's
exponent range raise PariError. Negative inputs retain the free factorial ValueError.

```typescript
const f = factorial(5n, 'pari');
[f.s, f.e, f.m, f.p]              // [1, 6n, 17293822569102704640n, 64]
factorial(1n << 54n, 'pari').e     // 946788236117799971n
```

At the dependency boundary, `@sagemath-ts/parigp-ts` exports
`mpfactr(n: bigint, precisionBits = 64): MpReal<bigint>`. Precision must be a positive
multiple of 64 bits. `MpReal<E = number>` retains number exponents by default for
existing quadratic-form kernels; factorial and exponential results use bigint
exponents. Full cypari2 Gen methods and global precision settings are not exposed.
The same package exports `divru(x: MpReal, n: number | bigint): MpReal` for unsigned
word divisors (`1..2^64-1`). Use bigint above the safe integer range.
Real zeros have `m = 0n`; `p` records allocated precision and may be nonzero.
`itor(0n, 128)` returns `{s: 0, e: -128, m: 0n, p: 128}`.
`rtor` retains the requested allocation and clamps zero accuracy to `min(e, -p)`.
`mulrr(x, x)` selects PARI's square kernel; equal but distinct operand records can
round differently because the multiplication and square thresholds differ.
`sqrtr(x: MpReal): MpReal | MpComplex` returns a complex result for negative input;
`MpComplex` has `re` and `im` components, each an exact bigint or an `MpReal`.
`sqrtr_abs(x)` returns the real square root of the magnitude. The low-level
`sqrti(n: bigint): bigint` also reads the integer magnitude, including negative inputs.
`gcvtoi(x: MpReal): [bigint, number]` returns the truncated integer and native error
exponent. Exact cancellation retains finite real accuracy: `gcvtoi(itor(1n, 64))`
returns `[1n, -63]`.

```typescript
import { itor, sqrtr, sqrti } from '@sagemath-ts/parigp-ts';
sqrtr(itor(-4n, 64))
// { re: 0n, im: { s: 1, e: 1, p: 64, m: 9223372036854775808n } }
sqrti(-9n) // 3n
```


```typescript
import { divru, real_1 } from '@sagemath-ts/parigp-ts';
divru(real_1(64), (1n << 64n) - 1n)
// { s: 1, e: -64, p: 64, m: 9223372036854775809n }
```


```typescript
new Integer(25n).factorial().value // 15511210043330985984000000n
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
ZZ.__call__()                    // 0n; null also maps to zero
ZZ.__call__('010')               // 10n; Integer uses Sage's PEP 3127 prefix rules
ZZ.__call__('1_000')             // 1000n
new Integer('ff', 16n).value    // 255n; optional base for strings
new Integer([1n, 0n, 1n], 2n).value // 5n; least-significant digit first
new Integer(new Rational(3n)).value // 3n; nonintegral rationals raise TypeError
// Integer and ZZ also accept booleans, Integer wrappers, _integer_ hooks and integer lifts.
// ZZ.__call__(rational, explicitBase) reproduces Sage's Q_to_Z NotImplementedError.
new Integer(42n).factor()        // [[2n,1n],[3n,1n],[7n,1n]]  — the wrapper class
new Integer(2n).pow(-3n)         // Rational(1/8); pow returns Integer | Rational
new Integer(0n).exact_log(2n)    // '-Infinity'; exact_log returns bigint | '-Infinity'
new Integer(0n).log(2n)          // '-Infinity'; log returns Integer | '-Infinity'
new Integer(0n).log()            // '-Infinity'; natural log also handles huge integers
new Integer(5n).__invert__()     // Rational(1/5); __invert__ always returns Rational
new Integer(0n).valuation(2n)    // 'Infinity'; valuation and ord return bigint | 'Infinity'
new Integer(-1n).popcount()      // 'Infinity'; popcount and hamming_weight use the same union
new Integer(5n).mod(-7n)         // Integer(5); ideal reduction uses |modulus|
new Integer(5n).__mod__(-7n)     // Integer(-2); % follows the divisor's sign
new Integer(5n).mod(0n)          // Integer(5); reduction by the zero ideal
new Integer(-1n).catalan_number() // Rational(-1/2); returns Integer | Rational
new Integer(-179n).class_number() // 5n; delegates to the PARI quadratic backend

// Rationals
new Rational().toString()                      // '0'
Rational.from(new Integer(7n)).toString()       // '7'; accepts IntegerLike
Rational.from().toString()                     // '0'; null and false also give zero
QQ.__call__().toString()                       // '0'
QQ.__call__(null).toString()                   // '0'; QQ.__call__([null]) raises TypeError
QQ.__call__(true).toString()                   // '1'
new Rational(3n, 4n).add(new Rational(1n, 2n))   // 5/4
Rational.from(1 / 3).toString()                 // '1/3'; Sage's simplest rational
Rational.fromString('010').toString()          // '8'; GMP base-zero integer syntax
// Rational.fromString('1.5') raises TypeError; use Rational.from(1.5).
// new Rational(1n, 0n) raises ValueError.
new Rational(1n, 2n ** 1074n).toNumber()        // 5e-324; correctly rounded binary64
new Rational(0n).ord(2n)                       // 'Infinity'; bigint | 'Infinity'
new Rational(0n).gamma().toString()            // 'Infinity'; Rational | UnsignedInfinityElement
new Rational(250n).roundToRational(-2).toString() // '200'; ties to even
QQ.cardinality()                              // 'Infinity'; inherited set cardinality
[...QQ.some_elements()].length                 // 100, in Sage's exact order
[...QQ.primes_of_bounded_norm_iter(new Rational(5n, 2n))] // [2n, 3n]; bound is rounded up
// QQ.gen(n: unknown = 0) tests equality to zero; QQ.gen(0n) returns one.
// QQ.gen(1) raises IndexError; QQ.__call__([new Rational(1n, 2n), 3n]) raises TypeError.
// QQ.random_element(numBound?: IntegerLike, denBound?: IntegerLike): Rational
// With no bounds it uses Sage's unbounded distribution; numerator bound 0 means 2.
new Integer(2n).nth_root(2147483647n, true)     // [Integer(1), false]
// Integer.nth_root and Rational.nth_root/is_nth_power require a signed 32-bit exponent.
// An out-of-range exponent raises OverflowError, as in Sage's Cython interface.

// Modular integers (Z/nZ)
Zmod() === ZZ                                 // true; default order is zero
Zmod(0n) === ZZ                               // true; quotient by zero is ZZ
Zmod(-7n).modulus                             // 7n
Mod(10n, 0n)                                  // 10n; preserves the original value/type
Mod(10n, 12n, Zmod(7n)).value                  // 3n; an explicit parent takes precedence
Mod(new Rational(1n, 2n), 12n, GF(7n)).is_square() // true; explicit prime parent
// Factory aliases: Integers and IntegerModRingFactory have the same signature as Zmod.
// Zmod accepts an optional IntegerLike, number, boolean or Rational order.
// Literal nonzero orders return IntegerModRing; zero/default returns IntegerRing (ZZ).
// A dynamic order returns IntegerModRing | IntegerRing: narrow before finite-ring methods.
// Mod(value, modulus, parent?) returns the original value for zero, otherwise an element
// of the explicit modular/prime parent or the factory ring. Dynamic moduli widen its type.
// new IntegerModRing(order) is the direct positive-order constructor; it never returns ZZ.
// Repeated Zmod calls use Sage's signed-order weak cache. Zmod(7) and Zmod(-7) are distinct.
const Z7 = Zmod(7n);
Z7.__call__(3n).mul(Z7.__call__(5n))             // 1  (15 mod 7)
Mod(10n, 7n)                                     // 3  — shorthand for a single element
Zmod(12n).__call__(new Rational(1n, 5n)).toString() // '5'; explicit rational conversion
Zmod(12n).__call__().toString()                   // '0'; null also means zero
Mod(3n, 12n).add(Mod(2n, 30n)).modulus           // 6n; quotient-ring pushout
Mod(3n, 12n).add(GF(2n).one()).toString()        // '0'; result belongs to GF(2)
// IntegerMod add/sub/mul/div accept IntegerLike, numbers, booleans and modular elements.
// Prime-field operands return that prime-element type when a canonical map exists.
// Compatible extension operands return an extension element; this also applies to prime elements.
Mod(2n, 7n).mul('ab')                           // 'abab'; Python sequence repetition
Mod(2n, 7n).mul([3n, 5n])                       // [3n, 5n, 3n, 5n]
// Native modular/prime mul(string): string and mul<T>(readonly T[]): T[] preserve item identity.
// Extension elements reject sequence multiplication, including zero and constant elements.
// eq(other: unknown) returns false when there is no compatible arithmetic parent.
// pow(exponent: IntegerLike | number | Rational | boolean | string | null);
// new IntegerMod(value: unknown, parent).
// IntegerModRing.__call__(value?: unknown), constructor(order: IntegerLike | number).
// random_element(bound?: IntegerLike | number | null) uses Sage's seeded Python stream;
// with a bound it reduces a uniformly drawn integer in [-bound, bound].

// Finite fields — prime and extension
const F13 = GF(13n);
F13.__call__(5n).pow(12n)                        // 1  (Fermat)
F13.__call__(5n).pow(new Integer(-2n)).toString() // '12'; exponent coercion preserves Sage errors
// Both prime-element pow methods share IntegerMod's exponent types and native/GMP dispatch.
// Non-integral Rational/number exponents raise TypeError; null is converted to exponent zero.
F13.__call__(new Integer(20n)).toString()          // '7'
F13.__call__(new Rational(1n, 2n)).toString()     // '7'; explicit rational conversion
F13.__call__().toString()                        // '0'; null and undefined also mean zero
F13.__call__(1n).eq(new Integer(14n))            // true
// add/sub/mul/div/eq also accept Integer wrappers and booleans.
// Implicit arithmetic with Rational raises TypeError; convert it with __call__ first.
const F4 = GF(4n);                               // extension field GF(2^2)
const F8 = GF(8n, 'b');                          // GF(2^3), generator named 'b'
F8.gen().frobenius(-1)                          // inverse Frobenius (powers reduce modulo degree)
// Extension.fromInteger(n: bigint) requires 0 <= n < cardinality(); otherwise ValueError.
// PrimeField.__call__(x?: unknown) and FiniteFieldExtension.__call__(x?: unknown)
// accept Integer, Rational, boolean, integer strings, null/undefined and existing elements.
// Extension coefficient arrays and polynomial inputs are converted coefficient by coefficient.
GF(9n).__call__('a^2 + 1').toString(); // 'a + 2'
GF(9n).__call__('1/a*a').toString(); // '1'
// Strings must first represent a polynomial: 'a^-1' is rejected; use gen().pow(-1n) for inversion.
// Incompatible extension parents raise TypeError.
Mod(2n, 12n).add(GF(9n).gen()).toString(); // 'a + 2'; promoted to GF(9)
GF(3n).one().eq(GF(9n).one()); // true; canonical prime-subfield embedding
// Prime/extension arithmetic accepts compatible elements from either prime implementation and Zmod.
// Nonintegral numbers, null, strings and arrays have no canonical arithmetic conversion;
// equality with them returns false, while native prime/modular multiplication repeats sequences.
// Prime-field sqrt() follows Sage extend=True and may return a quadratic-extension element.
GF(7n).__call__(3n).sqrt().toString(); // 'sqrt3'
GF(7n).__call__(3n).sqrt().pow(2n).toString(); // '3'
GF(7n).__call__(2n).sqrt({ all: true }).map(String); // ['3', '4']
GF(7n).__call__(3n).sqrt({ extend: false, all: true }); // []
// sqrt({ extend: false }) returns an element in the original field; nonsquares raise ValueError.
// random_element() uses Sage's Python random stream and the bundled Sage version's draw order.
```

`GF` and `FiniteField` alias `GFExtended(q: IntegerLike | number,
nameOrOptions: string | FiniteFieldOptions = 'a', options?: FiniteFieldOptions)`.
They return `PrimeField` for a prime order and `FiniteFieldExtension` for an extension.
`FiniteFieldOptions` is exported from `sagemath-ts/rings/finite_rings` and accepts
`name?: string` and `modulus?: readonly (IntegerLike | number)[] |
Polynomial<PrimeFieldElement> | FiniteFieldModulusAlgorithm | null`.

```typescript
GF(9n, { name: 'b', modulus: [1n, 0n, 1n] }).gen().pow(2n).toString(); // '2'
GF(7n, 'a', { modulus: [-2n, 1n] }).gen().toString(); // '2'
GF(9n, { name: 'b', modulus: 'conway' }).gen().toString(); // 'b'
```

Factory coefficient lists include the leading term. The existing `GFpn` coefficient
lists omit the monic leading 1. Unknown options such as `check` raise TypeError.
The separate prime-only `GF(order, { check })` in `finite_field_constructor.ts`
is reachable by importing that file directly. Broader Sage factory keywords remain unsupported.
The direct `PrimeField(p: unknown, options?: { modulus?:
Polynomial<PrimeFieldElement> | FiniteFieldModulusAlgorithm })` constructor also accepts
custom linear moduli; import that class from `rings/finite_rings/finite_field_extension`.

Both prime-element implementations accept `sqrt(options?: { extend?: boolean; all?: boolean })`.
The default `extend: true, all: false` returns a prime-field or extension element. Explicit
`extend: false, all: false` narrows the return type to the original element type, and `all: true`
returns an array of base-field roots. As in Sage, the unimplemented extension/all branch raises
NotImplementedError, while Sage's small-modulus search can return an empty array first.

`PrimeField` methods: `__call__`, `zero`, `one`, `gen`, `cardinality`, `elements`,
`is_field`, `random_element`, `multiplicative_generator`, `primitive_element`.

Note that `PrimeField` is not assignable to `FiniteFieldPrime`, the type `EllipticCurve`
declares for its base field — the two prime-field classes have diverged. This is invisible
at runtime (`EllipticCurve(GF(101n), …)` works), but it is a type error wherever the code
is actually typechecked.

`PolynomialRing<C>.__call__(x?: unknown)` validates inputs through the coefficient ring.
Omitted/null inputs return zero, and coefficient arrays are coerced before normalization.
Calling a ring on its own polynomial preserves identity; an element of the exact polynomial
base ring becomes a constant. `gen(n: unknown = 0)` requires equality to zero and rejects
list indices as unhashable. `PolynomialRingBase` exposes these same signatures.

```typescript
import { PolynomialRing } from 'sagemath-ts/rings';
const P = new PolynomialRing(GF(7n), 'x');
P.__call__([9n, 3n, 7n]).toString()             // '3*x + 2'
P.__call__().toString()                        // '0'
P.gen(0n).toString()                           // 'x'
QQ.__call__(GF(7n).__call__(-2n)).toString()    // '5'
```

Both prime-element classes and IntegerMod expose `_rational_(): Rational`, lifting the
canonical residue into QQ. Extension parents accept modular coefficients when their
characteristic divides the source modulus.

Over prime fields, strings support `+`, `-`, `*`, `/`, `^`, `**`, parentheses, numeric
factorials, comparisons, and implicit multiplication before a name (`2a`). Fractions must
cancel to a polynomial. General parser function calls, sequences, matrices and symbolic promotions remain
incomplete; unknown functions are rejected. GF(2) and large-prime polynomial ring displays include Sage's backend label.

The internal `misc/parser.ts` module has no package subpath export. It exports `Tokenizer`
(`next`, `peek`, `last`, `backtrack`, `reset`, `last_token_string`, `test`), `token_to_str`,
token constants, `LookupNameMaker<T>` (`__call__`, `set_names`), and `Parser<T>`.
Parser takes `(make_int, make_float, make_var, make_function = {},
implicit_multiplication = true, operations)`; `operations` supplies required
`binary(op, left, right)` and `unary(op, value)` callbacks. Name constructors accept a
callback or dictionary. Arithmetic methods are `parse`, `parse_expression`, `p_eqn`, `p_expr`,
`p_term`, `p_factor`, `p_power`, `p_atom`, `_variable_constructor`, `_callable_constructor`.
Sequence/matrix entry points and calls to supplied functions remain explicit stubs.

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

Short-model point methods: `add`, `sub`, `neg`, `double`, `mul`, `order`, `has_order`, `isZero`,
`weil_pairing`, `tate_pairing`, `ate_pairing`.
`order(options?: {algorithm?: 'pari' | 'generic_small' | 'hybrid'})` and
`additive_order(options?)` share Sage's algorithm selection and point cache.
`_compute_order(algorithm?)` performs the uncached calculation for a nonzero point.
The default/PARI path also populates the curve's `_order: bigint | null` cache.
Short-model `mul(n: IntegerLike | number)` coerces through ZZ, delegates to PARI
and propagates a known point order. `pari_curve()`, `__pari__()` and `toPari()`
return the same cached short PARI record (implicit a1=a2=a3=0). PARI ellmul
converts these records through native ellinit invariants; characteristic three
uses the finite-field model and FlxqE scalar backend.

Short-model curve and point orders delegate to the PARI port (`ellcard`/`ellorder`).
Cardinality uses exhaustive counting for tiny primes, CM formulas when applicable,
then Shanks below `expi(p) = 56` and SEA at or above that threshold.
The middle-range counter follows PARI's word Shanks search with exact integer bounds.

General Weierstrass isomorphisms use distinct polynomial roots in every
characteristic, including extension fields and QQ. Tuple order is the native
polynomial-root order; `WeierstrassIsomorphism(E, null, F)` selects its first tuple.

```typescript
import { GFpn } from 'sagemath-ts/rings/finite_rings';
import { EllipticCurveGeneric, _isomorphisms } from 'sagemath-ts/schemes/elliptic_curves';
const K = GFpn(2n, 2, [1, 1], 'a');
const C = new EllipticCurveGeneric(K, [K.zero(), K.zero(), K.one(), K.zero(), K.zero()]);
const automorphisms = [..._isomorphisms(C, C)];
automorphisms.length; // 24
automorphisms[0].map(String); // ['1', '0', '0', '0']
```

`C.isomorphism_to(F)` returns the first generator tuple. `C.isomorphisms(F)`
and `C.automorphisms()` return sorted tuple lists, putting identity and negation
first when present. `C.is_isomorphic(F)` checks base-parent equality and stops
as soon as it finds an isomorphism.

```typescript
C.automorphisms().slice(0, 2).map(t => t.map(String));
// [['1', '0', '0', '0'], ['1', '0', '0', '1']]
C.is_isomorphic(C); // true
C.isomorphism_to(C).map(String); // ['1', '0', '0', '0']
```

`WeierstrassIsomorphism._comparison_impl(left, right, op)` returns a boolean
for morphism operands; its unknown-operand overload returns `boolean | null`.
Use `eq`, `ne`, `lt`, `le`, `gt`, `ge` or `==`, `!=`, `<`, `<=`, `>`, `>=`.
Unknown operands return `null` (Python's `NotImplemented`). Ordering morphisms
with unequal domains or codomains raises TypeError over QQ/finite fields.

The source-internal export `_same_base_ring(left: FieldRing, right: FieldRing):
boolean` lives in `schemes/elliptic_curves/types.ts`. It recognizes equivalent
prime implementations (including the chosen generator for custom linear
moduli) and compares named extension definitions, including the
modulus; it is used by the curve APIs above.

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

Core coercion and factorization contracts: prime and extension field `add`, `sub`, `mul`,
`div`, and `eq` accept integer scalars (`number | bigint`); extension elements also accept
prime-subfield elements. Arithmetic rejects incompatible field parents, while equality returns
false. Explicit `__call__` conversion remains separate. `GFpn(p, n, coefficients)` checks
custom moduli for irreducibility (coefficients omit the leading 1). Polynomial arithmetic
requires compatible parents; explicit `R.__call__(f)` can rename a variable. Both `factor()`
and `squarefree_decomposition()` preserve a nontrivial unit as a degree-zero factor.

`AttributeError` and `OverflowError` are exported from `sagemath-ts` alongside the other Python-style errors.
`sigma(n, k)` preserves Sage’s `AttributeError` for negative `k` when `|n| > 1`.

The PARI dependency also exports `Fp_order(a: bigint, N: bigint, p: bigint)` from
`@sagemath-ts/parigp-ts`. It computes multiplicative order given a positive known multiple
`N`; inputs must satisfy `a^N = 1 (mod p)`. For example, `Fp_order(2n, 12n, 13n)` returns `12n`.

The PARI package exports `Z_isanypower(x: bigint): [number, bigint]`: maximal perfect-power
exponent and base, with `[0, x]` when no power is found (including `-1`, `0`, `1`).
For negative inputs, the maximal exponent is odd and the base is negative.
`Z_isanypower(676n)` returns `[2, 26n]`; `Z_isanypower(37n)` returns `[0, 37n]`;
`Z_isanypower(-64n)` returns `[3, -4n]`; `Z_isanypower(-16n)` returns `[0, -16n]`.


PARI's root export `znorder(a: bigint, modulus: bigint, order?: bigint): bigint` computes
multiplicative order using prime-power decomposition. The optional order is a positive
multiple of the true order. Unlike PARI's GEN argument, the residue and modulus are separate
BigInts; inputs must have nonzero modulus and an invertible residue.

Internal finite-ring adapters in `rings/finite_rings/integer_mod.ts` are exported for sharing
between the element classes: `canonicalFiniteOperands(left: FiniteArithmeticElement,
right: unknown, operation: string)` returns a pair in a common finite parent;
`finiteArithmeticEquals(left, right): boolean` handles incompatible equality, and
`repeatFiniteSequence<T>(count: bigint, sequence: string | readonly T[]): string | T[]`
implements native sequence repetition. `FiniteArithmeticElement` is the union of modular,
both prime-element classes and extension elements. These are a bounded implementation detail,
not Sage's general coercion API.


Direct `PrimeField(p: unknown, options?)` and `FiniteFieldPrime(p: unknown, check = true)`
coerce their order through Sage Integer conversion before validation. `check:false` on
the legacy class skips primality checking but still requires a positive order; it permits
order one. Arithmetic on these unchecked nonprime parents is still under audit.

All four finite-ring parents expose `gen(n: unknown = 0)`. Field parents accept a false
Python scalar/list/field-element index and otherwise raise `IndexError('only one generator')`;
quotient rings require equality to zero and raise `IndexError('n must be 0')`. For example,
`GF(7n).gen(false)` returns one, while `Zmod(7n).gen(null)` raises IndexError.

Modern prime and extension elements expose `_integer_(parent?: unknown): bigint`, which
ignores the optional parent as Sage does. Integer constructors call this hook automatically;
extension elements outside the prime subfield raise ValueError. The internal
`checkFiniteGeneratorIndex(index: unknown, field: boolean): void` adapter is exported from
`rings/finite_rings/integer_mod.ts` for the four parent classes.

```typescript
ZZ.__call__(GF(7n).__call__(5n)); // 5n
new Integer(GF(9n).__call__(2n)).value; // 2n
GF(7n).gen(false).toString(); // '1'
Zmod(7n).gen(new Rational(0n)).toString(); // '1'
// ZZ.__call__(GF(9n).gen()) raises ValueError: element is not in the prime field.
```


Extension `pow(exponent: unknown)` validates Sage's exponent comparison and conversion order.
Integer, integral Rational and prime-subfield finite-element exponents are supported;
nonconstant extension exponents raise ValueError. Unlike prime-element pow, extension pow
rejects string/list/null exponents during comparison. Zero to a negative nonintegral exponent
raises ZeroDivisionError before the integer-conversion check. Extension arithmetic delegates
to the PARI-port polynomial/quotient kernels and preserves the original field parent.

```typescript
GF(9n).gen().pow(new Rational(-1n)).toString(); // 'a + 2'
// GF(9n).zero().pow(new Rational(-1n, 2n)) raises ZeroDivisionError.
// GF(9n).gen().pow(new Rational(-1n, 2n)) raises TypeError.
```

`@sagemath-ts/parigp-ts` exports `FpXQ_pow(a, exponent: bigint, T, p: bigint)` and the
identical `FpXQ_powBig` alias. Inputs/results are ascending bigint coefficient arrays,
with zero represented by `[]`. Negative exponents invert first; a zero exponent returns
`[1n]` before quotient reduction. `FpXQ_inv(a, T, p)` is the original PARI signature;
the existing `FpXQ_inv(a, T, q, p)` form lifts to the prime-power modulus q=p^e.

```typescript
import { FpXQ_pow, FpXQ_inv } from '@sagemath-ts/parigp-ts';
FpXQ_pow([0n, 1n], -1n, [1n, 0n, 1n], 3n); // [0n, 2n]
FpXQ_inv([1n, 1n], [1n, 0n, 1n], 9n, 3n); // [5n, 4n]
```

The FLINT root exports dense, constant-first integer polynomial kernels:
`_fmpz_poly_gcd(a: readonly bigint[], b: readonly bigint[]): bigint[]`, and the same
signature for `_fmpz_poly_gcd_subresultant` and `_fmpz_poly_gcd_modular`.
`_fmpz_poly_gcd_heuristic(a, b): bigint[] | null` returns `null` when its heuristic fails.
For example, `_fmpz_poly_gcd([-6n, 0n, 6n], [9n, -9n])` returns `[-3n, 3n]`.
Zero is `[]`; inputs are not mutated. `_nmod_poly_gcd(a, b, p: bigint): bigint[]`
requires a prime modulus and returns a remainder that need not be monic.
These array kernels are implemented; the older mutable `fmpz_poly`/`nmod_poly` APIs
remain stubbed. Sage's integer-polynomial GCD delegates to the array backend.

`Polynomial.monic()` preserves a monic integer polynomial; otherwise it returns a polynomial
over QQ, even when the leading coefficient is -1. The Integer-coefficient overload of `monic`
and `_monic` returns an integer-or-rational polynomial (using the existing `RingElement`
intersection types); inspect `result.parent.base_ring` when the distinction matters.
`is_monic()` compares the leading coefficient with the base ring's one and is false for zero.
Zero normalization raises the backend's inversion error; it does not generally return zero.

FLINT exports `_nmod_poly_make_monic(a: readonly bigint[], n: bigint): bigint[]` for nonzero
polynomials with a unit leading coefficient modulo `n >= 2`. It does not mutate its input.
`_nmod_poly_make_monic([1n, 2n, 3n], 5n)` returns `[2n, 4n, 1n]`.


`Polynomial.getCoeff(n: unknown)`, `truncate(n: unknown)` and `reverse(degree?: unknown)`
preserve Sage's distinct index conversions. `getCoeff` requires an exact index; `truncate`
uses a signed 64-bit C-long conversion, including truncation of fractional values. Reversal
uses backend-specific unsigned conversion; over QQ it converts `degree + 1`, so degree -1
returns zero. `shift(n: number | bigint)` retains the polynomial return overload used by
finite numeric callers; its `unknown` overload returns `Polynomial<C> | null`. A nonzero
QQ/ZZ polynomial shifted by NaN returns null, so callers allowing nonfinite numbers must
handle that result. Zero-input shortcuts and result identity depend on the backend.

```typescript
const indexExample = new PolynomialRing(QQ, 'x').__call__([1n, 2n, 3n]);
indexExample.getCoeff(new Integer(1n)).toString(); // '2'
indexExample.truncate(new Rational(3n, 2n)).coeffs.map(String); // ['1']
indexExample.reverse(new Rational(3n, 2n)).coeffs.map(String); // ['2', '1']
indexExample.shift(new Rational(-3n, 2n))?.coeffs.map(String); // ['2', '3']
```

`Polynomial.add`, `sub` and `mul` accept a polynomial with a different coefficient type.
The existing same-type overload returns `Polynomial<C>`; the mixed overload returns
`Polynomial<C | D>`. Inspect `result.parent.base_ring` for the actual common coefficient
ring. The implemented ZZ/QQ/finite/nested parents are coerced before coefficient arithmetic,
including zero operands. Incompatible characteristics or variable names raise TypeError.
`eq(other: unknown): boolean` compares compatible polynomial and integer/rational/finite
scalar operands after coercion; coefficient reduction can lower the degree before equality
is tested. A coefficient polynomial embeds as a constant in its outer polynomial ring.
Nonintegral numbers compare after binary64 conversion over ZZ/QQ polynomial bases;
integral-valued JavaScript numbers represent Sage integers.

```typescript
const primePolynomial = new PolynomialRing(GF(7n), 'x').__call__([1n, 2n]);
const extensionPolynomial = new PolynomialRing(GF(49n, 'a'), 'x').__call__([3n, 1n]);
primePolynomial.add(extensionPolynomial).coeffs.map(String); // ['4', '3']; result is over GF(49)
new PolynomialRing(QQ, 'x').__call__([1n, -2n, -1n]).toString(); // '-x^2 - 2*x + 1'
```

`Polynomial.quo_rem` and `mod` also accept mixed coefficient types, returning polynomials
with coefficient type `C | D` in the common parent. Over ZZ, a nonintegral quotient
coefficient can leave a remainder of degree at least the divisor's degree. Generic
polynomial coefficient rings still require each quotient coefficient to lie in that ring.

FLINT exports `_fmpz_poly_divrem(a, b, exact?)` for dense constant-first bigint arrays.
The default/false flag returns `[quotient, remainder]`. With `exact: true`, it returns
`null` if a required leading quotient is not integral; a lower-degree nonzero remainder
is allowed. A dynamic boolean flag has the nullable return type. The nonzero divisor is
required, inputs are not mutated, and trailing zeros are normalized.
`_nmod_poly_divrem(a, b, n: bigint)` returns the corresponding reduced arrays modulo
`n >= 2`, with a nonzero divisor whose leading coefficient is a unit. These kernels
share the FLINT-port division algorithms used by polynomial GCD.

```typescript
import { _fmpz_poly_divrem, _nmod_poly_divrem } from '@sagemath-ts/flint-ts';
_fmpz_poly_divrem([-7n, -5n], [3n]); // [[-3n, -2n], [2n, 1n]]
_fmpz_poly_divrem([1n, 3n], [0n, 2n], true); // null
_nmod_poly_divrem([1n, 0n, 1n], [1n, 1n], 7n); // [[6n, 1n], [2n]]
```

`NTLError` is exported with the other Sage error classes and extends `RuntimeError`.
It preserves NTL's division-by-zero error for polynomials over large modular parents.

`Polynomial.gcd` accepts mixed coefficient types with the same `Polynomial<C | D>`
return convention as arithmetic. The result's parent is the canonical common parent.
Native finite backends return the original nonmonic operand when the other input is zero;
QQ normalizes and allocates a result. Composite nonunit failures preserve Sage's errors.

FLINT exports `_fmpq_poly_gcd(a, b)` for dense constant-first integer numerator arrays.
Input rational common denominators do not affect a monic GCD. The result is
`[numeratorCoefficients: bigint[], denominator: bigint]` with positive denominator;
zero is `[[], 1n]`. Inputs are normalized without mutation.

```typescript
import { _fmpq_poly_gcd } from '@sagemath-ts/flint-ts';
_fmpq_poly_gcd([-6n, -9n], [4n, 6n]); // [[2n, 3n], 3n]
```

`Polynomial.xgcd` accepts mixed polynomial coefficient types. Field results are three
`Polynomial<C | D>` values in the common parent. For `Polynomial<Integer & RingElement>`
over ZZ, each component is `Polynomial<Integer & RingElement> | Integer`: zero-input
shortcuts retain the nonzero polynomial and return Integer cofactors; two nonzero
constants return three Integers. Other ZZ inputs return three polynomials, using FLINT's
length-ordered resultant or a denominator-cleared QQ triple when it vanishes. Thus the
first component can differ in sign from `this.resultant(other)` after input reordering.
Native finite zero shortcuts retain nonmonic operands; over QQ both-zero inputs return
three distinct zero polynomials. Over a large prime using NTL, they return `(0, 1, 0)`.

FLINT exports `_fmpz_poly_xgcd(a, b): [bigint, bigint[], bigint[]]`, returning
`[r, s, t]` with `s*a + t*b = r`. A zero resultant returns `[0n, [], []]`.
For two nonzero constant arrays, the right constant must be a unit; otherwise RangeError
replaces the native reconstruction's undefined behavior. Sage handles constants separately.
`_fmpq_poly_xgcd(a, denA, b, denB)` returns three `DenseRationalPolynomial` pairs
`[numeratorCoefficients: bigint[], denominator: bigint]`, with monic GCD and positive
output denominators. Input denominators must be positive. Both kernels normalize dense
constant-first bigint arrays without mutation.

```typescript
import { _fmpz_poly_xgcd, _fmpq_poly_xgcd } from '@sagemath-ts/flint-ts';
_fmpz_poly_xgcd([2n, 1n], [3n, 0n, 0n, 1n]); // [5n, [4n, -2n, 1n], [-1n]]
_fmpq_poly_xgcd([2n, 2n], 6n, [], 1n); // [[[1n, 1n], 1n], [[3n], 1n], [[], 1n]]
```

FLINT exports `_nmod_poly_xgcd(a, b, n)` for dense constant-first bigint arrays modulo
`n >= 2`. It returns `[g, s, t]` with monic `g = s*a + t*b`, normalizes inputs without
mutation, and raises RangeError when a required inverse does not exist. Its zero-input
normalization follows FLINT directly, before Sage's finite-template zero shortcuts.

```typescript
import { _nmod_poly_xgcd } from '@sagemath-ts/flint-ts';
_nmod_poly_xgcd([2n, 2n], [], 7n); // [[1n, 1n], [4n], []]
```

FLINT exports `_fmpz_poly_resultant(a, b): bigint` for dense integer arrays and
`_nmod_poly_resultant(a, b, n): bigint` modulo `n >= 2`. Arrays are constant-first and
normalized without mutation. Modular inversions may raise RangeError over a composite
modulus. The modular kernel preserves native pointer aliasing: passing the same nonzero
constant array twice returns zero; passing two separate constant arrays returns one.
Sage polynomial resultants now delegate to these kernels for ZZ and word-sized prime fields.
QQ uses `_fmpq_poly_resultant(a, denA, b, denB): [bigint, bigint]`, which requires positive
input denominators and returns a reduced numerator and positive denominator.

```typescript
import { _fmpz_poly_resultant, _nmod_poly_resultant } from '@sagemath-ts/flint-ts';
_fmpz_poly_resultant([2n, 4n], [-6n, -9n]); // -6n
_nmod_poly_resultant([1n], [1n], 7n); // 1n
```


`Polynomial.sylvester_matrix(other, variable?)` returns coefficient rows as a dense array.
Mixed polynomial inputs return `(C | D)[][]`; scalar inputs return `RingElement[][]`.
Sage's canonical-parent, zero-input and optional-variable error order is preserved: a
constant left input can raise IndexError when parent conversion is needed. The optional
variable is converted, but does not otherwise change a univariate matrix. Default small
PARI-field matrix backend failures are preserved too (see DEVIATIONS.md).

Explicit polynomial construction from a foreign finite-field element can expand its
prime-field polynomial representative. A QQ coefficient list uses coefficient conversion
instead, and may reject that same element.

```typescript
import { QQ, Rational, GF } from 'sagemath-ts';
import { PolynomialRing } from 'sagemath-ts/rings';
const sylvesterRing = new PolynomialRing(QQ, 'x');
sylvesterRing.__call__([1n, 0n, 1n]).sylvester_matrix(new Rational(2n)).map(row => row.map(String));
// [['2', '0'], ['0', '2']]
sylvesterRing.__call__(GF(49n, 'a').gen()).coeffs.map(String); // ['0', '1']
```


`Polynomial.resultant(other, options?: { proof?: unknown })` accepts a polynomial or
scalar in a common canonical parent. Same-coefficient polynomial operands return `C`,
mixed operands return `C | D`, and scalar operands return `RingElement`. Only the integer
backend accepts `proof` (the native implementation ignores its value); supplying it after
promotion to another backend raises TypeError. Composite word-modular rings use Sage's
matrix determinant path, including its zero-polynomial error.

`RDF`, `RealDoubleField_class` and `RealDoubleElement` are exported from `sagemath-ts`
and `sagemath-ts/rings`. This is the numeric coefficient subset: `__call__` accepts numbers,
bigints, Integers, Rationals, booleans and existing RDF elements; elements expose `value`,
`parent`, `add`, `sub`, `mul`, `div`, `neg`, `pow(IntegerLike)`, `inv`, `eq`, `isZero` and `toString`. The field exposes
`zero`, `one`, `is_field` and `characteristic`. Broader Sage RDF methods are not yet ported.
Nonintegral number operands promote ZZ/QQ polynomial resultants to RDF coefficients.
General binary64 real resultants use native PARI Sylvester elimination. The SciPy
finite-overflow fallback and nonconstant nested real coefficients remain unimplemented
(see DEVIATIONS.md).

```typescript
import { QQ, RDF } from 'sagemath-ts';
import { PolynomialRing } from 'sagemath-ts/rings';
const resultantRing = new PolynomialRing(QQ, 'x');
resultantRing.__call__([1n, 0n, 1n]).resultant(0.5).toString(); // '0.25'
RDF.__call__(1n << 1200n).toString(); // '+infinity'
```

PARI exports `resultant(a: readonly number[], b: readonly number[]): number` for
constant-first binary64 coefficient arrays. It retains inexact trailing zero coefficients,
as PARI `Polrev` does; an all-zero array still represents a zero polynomial. Sage removes
trailing zeros before calling this dependency. It uses native real conversion,
exponent-based pivots and whole-word elimination arithmetic; NaN, infinity and conversion
overflow raise the exported `PariError`. Native `rtodbl` flushes the smallest subnormals to positive zero.

```typescript
import { _fmpq_poly_resultant } from '@sagemath-ts/flint-ts';
import { PariError, resultant } from '@sagemath-ts/parigp-ts';
_fmpq_poly_resultant([1n, 1n], 2n, [1n, -1n], 3n); // [1n, 3n]
resultant([1, 2, 3, 4], [0.5]); // 0.125
try { resultant([Infinity, 1], [1, 1]); } catch (error) { error instanceof PariError; } // true

```


```typescript
import { PariError, resultant } from '@sagemath-ts/parigp-ts';
resultant([1, 0, 1], [-2, 0, 0, 1]); // 5
resultant([1, 0], [1, 1]); // -1; PARI retains the inexact leading zero
```


`Polynomial.derivative(...args)` accepts variables and nonnegative integer repetition
counts. A single list argument is already expanded: its integer entries are treated as
variables, not counts. Zero counts and an empty single list return the original polynomial.
`diff` and `differentiate` are aliases. `_derivative(variable?)` differentiates once;
`gradient()` returns the one-element list of partial derivatives. Nested polynomial
variables recurse into coefficients. Generic zero derivatives return the same polynomial;
QQ delegates to FLINT and allocates a new result.

```typescript
import { QQ, Rational } from 'sagemath-ts';
import { PolynomialRing } from 'sagemath-ts/rings';
const derivativeRing = new PolynomialRing(QQ, 'x');
const derivativeExample = derivativeRing.__call__([-1n, 0n, new Rational(1n, 2n), 0n, -1n]);
derivativeExample.derivative(derivativeRing.gen(), 2n).coeffs.map(String); // ['1', '0', '-12']
derivativeExample.diff([]) === derivativeExample; // true
```

`sagemath-ts/misc/derivative` exports `derivative_parse(args: readonly unknown[]): unknown[]`
and `multi_derivative(F, args)` for objects implementing `_derivative(variable?)`. Counts
are BigInt, Integer, integral number or boolean values and must fit a signed C int before
the nonnegative check. Other values are variables. `derivative_parse([list])` returns that
same list. Symbolic differentiation and general non-polynomial coefficient hooks remain
outside the supported polynomial subset.

FLINT exports `_fmpz_poly_derivative(a: readonly bigint[]): bigint[]` and
`_fmpq_poly_derivative(a: readonly bigint[], den: bigint): [bigint[], bigint]`.
Arrays are constant-first, normalized without mutation; the rational input denominator
must be positive and the returned numerator/denominator is canonical.

```typescript
import { derivative_parse } from 'sagemath-ts/misc/derivative';
import { _fmpq_poly_derivative } from '@sagemath-ts/flint-ts';
derivative_parse(['x', 2n, 2n]); // ['x', 'x', null, null]
_fmpq_poly_derivative([1n, 0n, 3n], 6n); // [[0n, 1n], 1n]
```

`Polynomial.pseudo_quo_rem(other: unknown)` returns `[quotient, remainder]`. Its quotient
can be a `Polynomial` or a `FractionElement` (including `FpTElement`); coefficient parents can change during
arithmetic. Degree gaps use negative powers of the divisor's leading coefficient. Scalar
operands must supply Sage's zero-test protocol (BigInt and Integer are supported); the
original operand is checked before conversion. Dividing a zero polynomial by a nonzero
constant polynomial can return a fraction-field zero quotient.

`PolynomialRing.fraction_field()` and `FractionField(R)` (exported from
`sagemath-ts/rings/fraction_field`) return a cached polynomial fraction-field parent.
Over fields this selects `FractionField_1poly_field`; odd word primes below 46,341
select its `FpT` subclass. Direct `FractionField_generic(R)` keeps the generic protocol.
The supported fraction arithmetic covers univariate polynomials over ZZ, QQ and finite
fields, including nonconstant denominators. Polynomial/fraction conversion follows the
source partial sections: generic fractions require a unit denominator; native FpT
sections can normalize their input in place. Direct scalar conversions cover the compared ZZ, QQ and finite-field parents;
general nested coefficient backends remain under audit.

Parents provide `__call__(numerator, denominator = 1n)`, cached `zero()`/`one()`, `ring()`,
`base_ring()`, `is_field()`, `is_exact()`, `characteristic()`, `gen(i = 0)` and `ngens()`.
Elements provide `numerator()`, `denominator()`, `add`, `sub`, `mul`, `div`, `neg`, `inv`,
`pow(number | bigint)`, `eq`, `isZero`, `parent` and `toString`. Direct constructors accept
`{ coerce?: boolean, reduce?: boolean }` as their fourth argument. Generic elements and
`FractionFieldElement_1poly_field` expose `reduce()`; `FpTElement` has `numer()`/`denom()`
aliases and no `reduce()` method. `FpT` and `FpTElement` are exported from
`sagemath-ts/rings/fraction_field_FpT`. Unreduced native fractions retain Sage's structural
equality and zero test; their numerator/denominator accessors return fresh polynomials.

```typescript
import { QQ } from 'sagemath-ts';
import { PolynomialRing } from 'sagemath-ts/rings';
import { FractionField_generic } from 'sagemath-ts/rings/fraction_field';
const pseudoRing = new PolynomialRing(QQ, 'x');
pseudoRing.__call__(1n).pseudo_quo_rem(pseudoRing.__call__([2n, 0n, 0n, 2n])).map(String); // ['0', '1/4']
const constantFractions = new FractionField_generic(pseudoRing);
constantFractions.__call__([1n, 2n], 3n).mul(constantFractions.__call__(2n, 3n)).toString(); // '4/9*x + 2/9'
```

```typescript
import { FractionField } from 'sagemath-ts/rings/fraction_field';
const fractions = FractionField(pseudoRing);
fractions.__call__([1n, 2n, 1n], [1n, 1n]).toString(); // 'x + 1'
fractions.__call__([1n, 1n]).inv().toString(); // '1/(x + 1)'
```

`PolynomialRing.has_coerce_map_from(other: CoefficientRing<RingElement>)` reports the
supported canonical parent maps. Polynomial variables must agree unless the source embeds
as a coefficient. An exact generic fraction section with denominator one retains numerator
identity.

```typescript
const sectionRing = new PolynomialRing(QQ, 'x');
const sectionField = sectionRing.fraction_field();
sectionRing.has_coerce_map_from(new PolynomialRing(QQ, 'x')); // true
sectionRing.has_coerce_map_from(new PolynomialRing(QQ, 'y')); // false
const sectionValue = sectionField.__call__([1n, 2n, 1n], [1n, 1n]);
sectionRing.__call__(sectionValue).toString(); // 'x + 1'
sectionRing.__call__(sectionValue) === sectionValue.numerator(); // true
```

Generic fraction elements provide `_conversion<T>(R: { __call__(x: unknown): T }): T`,
`_integer_(R: IntegerRing): bigint` and `_rational_(): Rational`. These hooks can reduce
the element in place. FpT uses native partial sections and has no such public hooks.
Polynomials provide `_scalar_conversion<T>(R: { __call__(x: unknown): T }): T`,
`_integer_(R: IntegerRing): bigint` and `_rational_(): Rational`. Direct hooks and parent
constructors can use different section maps and report different errors for nonconstants.
The GF2 parent exposes `__call__(x?: unknown): GF2Element` with runtime checking;
integer strings reduce exactly, and null/undefined or an omitted input give zero.
GF2 elements also provide `_integer_(parent?: unknown): bigint` and `_rational_(): Rational`.

```typescript
import { FractionFieldElement } from 'sagemath-ts/rings/fraction_field_element';
const scalarRing = new PolynomialRing(QQ, 'x');
QQ.__call__(scalarRing.__call__(new Rational(3n, 2n))).toString(); // '3/2'
scalarRing.__call__(3n)._integer_(ZZ); // 3n
scalarRing.__call__(3n)._rational_().toString(); // '3'
scalarRing.__call__(3n)._scalar_conversion(QQ).toString(); // '3'
const scalarFraction = new FractionFieldElement(new FractionField_generic(scalarRing), [1n, 1n], [2n, 2n]);
scalarFraction._rational_().toString(); // '1/2'
scalarFraction._conversion(QQ).toString(); // '1/2'
```

FLINT exports `_nmod_poly_add(a,b,p)`, `_nmod_poly_sub(a,b,p)`, `_nmod_poly_mul(a,b,p)`
and `_nmod_poly_pow(a,e,p)` on constant-first BigInt arrays, with zero represented by `[]`.
The modulus satisfies `1 <= p < 2^64`; powers accept `0 <= e < 2^64`. Inputs are preserved
and coefficients are reduced. Explicit multiplication kernels `_nmod_poly_mul_classical`,
`_nmod_poly_mul_KS`, `_nmod_poly_mul_KS2` and `_nmod_poly_mul_KS4` share the product signature.

```typescript
import { _nmod_poly_mul, _nmod_poly_pow } from '@sagemath-ts/flint-ts';
_nmod_poly_mul([1n, 1n], [6n, 1n], 7n); // [6n, 0n, 1n]
_nmod_poly_pow([1n, 1n], 2n, 7n); // [1n, 2n, 1n]
```

### Polynomial roots

`Polynomial<C>.roots()` and `roots({ multiplicities: true })` return
`Array<[C, number]>` for field and integer/rational coefficients.
`roots({ multiplicities: false })` returns `C[]`;
a dynamic boolean option returns their union. Only base-ring roots are supported;
ring overrides and algorithm options are not exposed. For finite fields, the false
option first computes `gcd(f, x^q-x)` and factors that smaller squarefree polynomial.
Its order can differ from the multiplicity-preserving result. The `division_points`
helper uses this distinct-root route; its returned division polynomials expose
`roots()` and `roots({ multiplicities: false })`.

```ts
import { PolynomialRing } from 'sagemath-ts/rings';
import { GF } from 'sagemath-ts';
const rootRing = new PolynomialRing(GF(7n), 'x');
const rootX = rootRing.gen();
const rootPolynomial = rootX.pow(2n).mul(rootX.sub(rootRing.one()));
rootPolynomial.roots().map(([r, m]) => [String(r), m]); // [['1', 1], ['0', 2]]
rootPolynomial.roots({ multiplicities: false }).map(String); // ['0', '1']
```

### Modular-ring polynomial roots

For `Polynomial<IntegerMod & RingElement>`, default `roots()` returns
`Array<[FiniteFieldElement, number]>`, where the element class is from
`finite_field_prime.ts` and the parent is the ring's cached `field()`.
`roots({ multiplicities: false })` returns `IntegerMod[]` in the input ring;
a dynamic boolean returns the union. Composite rings reject multiplicities,
including constant polynomials. Distinct zero roots list every residue class.
Composite distinct roots use factorization of the modulus, Hensel lifting and CRT.

`IntegerModRing.field(): FiniteFieldPrime` is cached and raises
`ValueError('self must be a field')` for a nonprime modulus.
`factored_order(): Array<[bigint, bigint]>` caches the factorization array.
`_roots_univariate_polynomial(f, options?)` exposes the ring hook, returning
`Array<[FiniteFieldElement, number]> | IntegerMod[]`; options are
`{ ring?: IntegerModRing | null; multiplicities?: boolean; algorithm?: unknown }`.
The hook accepts only its own ring or null, and ignores algorithm as Sage does.
`IntegerModRing._lift_residue_field_root(p, e, f, fprime, root): IntegerMod[]`
is static; p/e accept IntegerLike, p must be prime and e positive, f/fprime are
polynomials over Zmod(p^e), and root must be a root in Zmod(p).

```ts
import { Zmod } from 'sagemath-ts';
import { PolynomialRing } from 'sagemath-ts/rings';
const K = Zmod(7n), R = new PolynomialRing(K, 'x');
String(K.field()); // 'Finite Field of size 7'
K.field() === K.field(); // true
K.factored_order(); // [[7n, 1n]]
const f = R.__call__([0n, -1n, 1n]);
f.roots().map(([r,m]) => [String(r),m,r.parent === K.field()]);
// [['0', 1, true], ['1', 1, true]]
f.roots({ multiplicities: false }).map(r => [String(r),r.parent === K]);
// [['0', true], ['1', true]]
const S = new PolynomialRing(Zmod(8n), 'x');
S.__call__([-1n,0n,1n]).roots({ multiplicities: false }).map(String);
// ['1', '5', '3', '7']
```

**Bundled upstream bug:** the noninvertible linear branch recurses without dividing
coefficients by their common factor. The compatibility port preserves that result;
it is not a guarantee that every returned value is a mathematical root:

```ts
const T = new PolynomialRing(Zmod(4n), 'x');
T.__call__([2n,2n]).roots({ multiplicities: false }).map(String); // ['0', '2']
// These are not roots of 2*x + 2 modulo 4; the actual roots are 1 and 3.
```

### Polynomial integer powers

`Polynomial<C>.pow(n: unknown, modulus?: unknown)` returns `Polynomial<C> | FractionElement<C>`.
A nonnegative numeric literal retains the narrower `Polynomial<C>` overload. Dynamic or
negative exponents require narrowing before calling polynomial-only methods. This corrects
an earlier return type that could not represent negative powers (a TypeScript breaking change).
Integer powers delegate to the coefficient backend; zero-negative-power behavior and signed
word bounds follow that backend. Over ZZ and QQ, exact rational exponents use polynomial
roots when they exist; unsupported roots raise the original error. An optional modulus
follows the coefficient backend; ZZ/QQ reject it with NotImplementedError.

```typescript
const powerRing = new PolynomialRing(QQ, 'x');
const powerBase = powerRing.gen().add(powerRing.one());
powerBase.pow(3n).toString()  // 'x^3 + 3*x^2 + 3*x + 1'
powerBase.pow(-2n).toString() // '1/(x^2 + 2*x + 1)'
const powerExponent: bigint = -2n;
const powerResult = powerBase.pow(powerExponent);
if (!(powerResult instanceof Polynomial)) powerResult.denominator().toString(); // 'x^2 + 2*x + 1'
```

The FLINT root exports `_fmpz_poly_pow(a: readonly bigint[], e: bigint): bigint[]`, its
`_fmpz_poly_pow_small`, `_fmpz_poly_pow_binomial`, `_fmpz_poly_pow_multinomial`
and `_fmpz_poly_pow_binexp` variants, and
`_fmpq_poly_pow(a: readonly bigint[], denominator: bigint, e: bigint): [bigint[], bigint]`.
Exponents are unsigned words; `_small` accepts 0–4 and `_binomial` requires a normalized length-two
input. Rational power preserves the raw numerator/common-denominator
representation: `_fmpq_poly_pow([2n], 2n, 3n)` is `[[8n], 8n]`.
`_fmpz_poly_pow([1n, 1n], 3n)` is `[1n, 3n, 3n, 1n]`.

The NTL root exports `GF2X_power(a: GF2X, e: bigint): GF2X` (also `GF2X.power`),
`ZZ_pX_power(a: readonly bigint[], e: bigint, p: bigint): bigint[]` and
`ZZ_pEX_power(a: readonly (readonly bigint[])[], e: bigint, f: readonly bigint[], p: bigint): bigint[][]`.
Exponents are nonnegative signed words. Extension modulus `f` is monic; inner arrays are
coefficients in ascending powers of its generator. The array adapters do not implement the
still-stubbed stateful `ZZ_pX` class.
`GF2X_power(new GF2X([1, 1]), 3n).rep()` is `15n`;
`ZZ_pX_power([1n, 1n], 3n, 7n)` is `[1n, 3n, 3n, 1n]`;
`ZZ_pEX_power([[0n, 1n]], 2n, [1n, 1n, 1n], 7n)` is `[[6n, 6n]]`.


### Polynomial roots and truncated series

`Polynomial<C>` exposes `nth_root(n: unknown): Polynomial<C>`,
`_nth_root_series(n: unknown, prec: unknown, start?: unknown): Polynomial<C>`,
`inverse_series_trunc(prec: unknown): Polynomial<C>`,
`power_trunc(n: unknown, prec: unknown): Polynomial<C>`,
`_power_trunc(n: unknown, prec: unknown): Polynomial<C>`,
`_mul_trunc_(right: Polynomial<C> | null, n: unknown): Polynomial<C>` and
`multiplication_trunc(other: unknown, n: unknown): Polynomial<RingElement>`.
The public multiplication method resolves common coefficient parents. The underscore method
expects polynomial operands in the same backend. Precision and exponent coercions follow
Sage's individual methods; they are not interchangeable integer parsers. The direct module
`rings/polynomial/polynomial_element` additionally exports
`generic_power_trunc<C>(p: Polynomial<C>, n: bigint, prec: number): Polynomial<C>`.

```typescript
const seriesRing = new PolynomialRing(QQ, 'x');
const seriesBase = seriesRing.__call__([1n, 1n]);
seriesRing.__call__([1n, 2n, 1n]).pow(new Rational(1n, 2n)).toString(); // 'x + 1'
seriesRing.__call__([1n, 2n, 1n]).nth_root(2n).toString(); // 'x + 1'
seriesBase._nth_root_series(2n, 3n).toString(); // '-1/8*x^2 + 1/2*x + 1'
seriesBase.inverse_series_trunc(4n).toString(); // '-x^3 + x^2 - x + 1'
seriesBase.power_trunc(5n, 3n).toString(); // '10*x^2 + 5*x + 1'
seriesBase.multiplication_trunc(seriesBase, 2n).toString(); // '2*x + 1'
```

FLINT exports `_fmpz_poly_mullow(a, b, n)`, `_fmpq_poly_mullow(a, denA, b, denB, n)`,
`_nmod_poly_mullow(a, b, n, p)`, `_fmpz_poly_pow_trunc(a, e, n)`,
`_nmod_poly_pow_trunc(a, e, n, p)`, `_fmpz_poly_inv_series(a, n)` and
`_fmpq_poly_inv_series_newton(a, den, n)`. Polynomial inputs are readonly bigint arrays;
precisions are numbers; exponents, moduli and denominators are bigints. Rational results
are `[bigint[], bigint]`; other results are `bigint[]`. Inverses require positive precision
and an invertible constant (±1 over ZZ). For example,
`_fmpq_poly_inv_series_newton([2n, 1n], 1n, 3)` returns `[[4n, -2n, 1n], 8n]`.

NTL exports `ZZ_pEX_InvTrunc(a: readonly (readonly bigint[])[], n: number,
f: readonly bigint[], p: bigint): bigint[][]`. It uses the same extension representation
as `ZZ_pEX_power`, accepts precision zero, and requires an invertible constant at positive
precision. `ZZ_pEX_InvTrunc([[1n], [1n]], 4, [1n, 1n, 1n], 7n)` returns
`[[1n], [6n], [1n], [6n]]`.


Truncated multiplication uses canonical coercion for scalar operands. A JavaScript
nonintegral number promotes a QQ polynomial to RDF; a string is not a canonical scalar.
The internal `_mul_trunc_` also accepts `null`, matching the compared native `None` case.
For a QQ polynomial `f = R.__call__([1n, 2n, 1n])`,
`f.multiplication_trunc(0.5, 3n).toString()` is `'0.5*x^2 + x + 0.5'` over RDF, and
`f._mul_trunc_(null, 3n).toString()` is `'0'`. Real-coefficient products preserve the
source's ordered sums, term shortcuts and aliased-square operations, including signed zeros.

### Full exact polynomial products

`Polynomial.mul(other)` delegates integer and rational products to FLINT, modular
products to FLINT/NTL, and extension products to NTL. Generic exact coefficients use
SageMath's Karatsuba thresholds, including nested polynomial and fraction fields.

FLINT exports `_fmpz_poly_mul(a, b): bigint[]` and
`_fmpq_poly_mul(a, denA, b, denB): [bigint[], bigint]`, with readonly bigint coefficient
arrays and positive bigint denominators. Full rational multiplication preserves native
cross-content cancellation; it does not canonicalize arbitrary input pairs globally.
For example, `_fmpq_poly_mul([2n, 2n], 2n, [1n, 1n], 1n)` returns
`[[2n, 4n, 2n], 2n]`. Passing the same array twice selects the native square behavior.

NTL exports `ZZ_pX_mul(a, b, p): bigint[]` and
`ZZ_pEX_mul(a, b, f, p): bigint[][]`. The first takes readonly bigint coefficient arrays;
the second takes readonly arrays of readonly bigint coefficient arrays over the monic
extension modulus `f`. Both take a bigint characteristic `p` and normalize inputs.
`ZZ_pX_mul([1n, 2n], [3n, 4n], 7n)` returns `[3n, 3n, 1n]`.
`ZZ_pEX_mul([[1n], [0n, 1n]], [[1n], [0n, 1n]], [1n, 0n, 1n], 7n)`
returns `[[1n], [0n, 2n], [6n]]`.


### Modular polynomial powers and RDF integer powers

A nonnegative numeric literal with a same-coefficient polynomial modulus retains the
polynomial return type; dynamic arguments return the polynomial/fraction union. Backend
rules can ignore a modulus for generic constants or negative powers. Negative native
powers can return a fraction-field reciprocal rather than a modular inverse.

```typescript
import { RDF, GF } from 'sagemath-ts';
import { PolynomialRing } from 'sagemath-ts/rings';
const modularRing = new PolynomialRing(GF(7n), 'x');
const mx = modularRing.gen();
mx.add(modularRing.one()).pow(5n, mx.pow(2n).add(modularRing.one())).toString(); // '3*x + 3'
RDF.__call__(2).pow(-3n).value; // 0.125
RDF.__call__(4).inv().value; // 0.25
```

FLINT exports `_nmod_poly_inv_series_newton(a, n, p)`,
`_nmod_poly_powmod_ui_binexp(a, e, f, p)`,
`_nmod_poly_powmod_fmpz_binexp_preinv(a, e, f, finv, p)` and
`_nmod_poly_powmod_x_fmpz_preinv(e, f, finv, p)`. Coefficient arrays are readonly bigint
arrays in increasing degree order; results are normalized bigint arrays. `e` is bigint,
nonnegative (unsigned 64-bit for the ui function); `n` is a positive safe integer number.
`finv` is the inverse series of reversed `f`, to at least `f.length` terms. The modulus
`p` and every required leading/constant inverse must be valid for native modular arithmetic.
The internal module `packages/flint-ts/src/nmod_poly/powmod_binexp_preinv.ts` also exports
`_nmod_poly_preinv_remainder(a, f, finv, p)` for already-reduced products of degree below
`2*deg(f)`; it is not a package-root export. With `a = [1n, 2n, 1n]`, `f = [1n, 0n, 1n]`,
`finv = [1n, 0n, 6n]`, and `p = 7n`, it returns `[0n, 2n]`.

```typescript
import * as FM from '@sagemath-ts/flint-ts';
const modulus = [1n, 0n, 1n];
const inverse = FM._nmod_poly_inv_series_newton(modulus.slice().reverse(), 3, 7n); // [1n, 0n, 6n]
FM._nmod_poly_powmod_ui_binexp([1n, 1n], 5n, modulus, 7n); // [3n, 3n]
FM._nmod_poly_powmod_fmpz_binexp_preinv([1n, 1n], 5n, modulus, inverse, 7n); // [3n, 3n]
FM._nmod_poly_powmod_x_fmpz_preinv(5n, modulus, inverse, 7n); // [0n, 1n]
```

NTL exports `ZZ_pEX_PowerMod(a, e, m, f, p)`, `ZZ_pEX_PowerXMod(e, m, f, p)` and
`ZZ_pEX_XGCD(a, b, f, p)`. Outer polynomials are readonly arrays of readonly bigint
coefficient arrays; `f` is the monic coefficient-field modulus and `p` its prime
characteristic. Powers accept signed bigint exponents, require positive-degree `m`,
and require `deg(a) < deg(m)`. XGCD returns `[g, s, t]` with monic `g = s*a + t*b`.
`ZZ_pEX_PowerXMod(5n, [[1n], [], [1n]], [1n, 0n, 1n], 7n)` returns `[[], [1n]]`.

`@sagemath-ts/gsl-ts` exports `gsl_pow_int(x: number, n: number): number`,
`gsl_sf_log(x: number): number` and `gsl_sf_exp(x: number): number`.
The integer exponent must fit signed 32 bits. `gsl_pow_int(2, -3)` is `0.125`,
`gsl_sf_log(1)` is `0`, and `gsl_sf_exp(0)` is `1`. This is only a three-function
numeric subset; see DEVIATIONS.md for native error-handler and libm limitations.

### Polynomial evaluation and composition

`f.evaluate(x)` implements Sage's `f(x)` dispatch; `f.compose(g)` uses the same dispatch
for polynomial substitution. Inputs include compatible scalars, polynomials, `IntegerMatrix`, `Matrix_modn_dense` and
`Matrix<RingElement>` (generic matrices), and `Matrix_mod2_dense`. Same-coefficient
composition returns `Polynomial<C>`; a different coefficient type returns `Polynomial<RingElement>`. `null`,
omitted input and an empty argument list substitute the generator; a list unwraps one
level. Extra arguments evaluate nested polynomial coefficients. Incompatible inputs can
raise even for a zero or constant polynomial. Use bigint/Integer for integers: JavaScript
numbers mean Python floats here and can produce a number or an RDF element. A same-coefficient
input or bigint/Integer/boolean retains its coefficient return type; dynamic unknown inputs return `RingElement | number | IntegerMatrix | Matrix_modn_dense | Matrix_mod2_dense | Matrix<RingElement>`.
Matrix inputs return that matrix union because coefficient promotion can change the concrete class.

`new Polynomial(coefficients, parent, is_gen?)` accepts an optional boolean generator flag
(default false). QQ and binary backends initialize the indeterminate when it is true; other
array constructors retain their supplied coefficients. Prefer `R.gen()` for a generator: it is cached, and `f.is_gen(): boolean | bigint`
reflects the original backend's distinguished-flag or native indeterminate test. Large
composite-modulus NTL polynomials return `0n`/`1n`, matching that native integer predicate.

```typescript
import { GF, RDF } from 'sagemath-ts';
import { PolynomialRing, CompiledPolynomialFunction } from 'sagemath-ts/rings/polynomial';
const evaluationRing = new PolynomialRing(GF(7n), 'x');
const ex = evaluationRing.gen();
const ef = ex.pow(2n).add(evaluationRing.one());
ef.evaluate(3n).toString(); // '3'
ef.compose(ex.add(evaluationRing.one())).toString(); // 'x^2 + 2*x + 2'
ex === evaluationRing.gen(); // true
ex.is_gen(); // true
const compiled = new CompiledPolynomialFunction([1, 0, 0, 0, 1].map(v => RDF.__call__(v)));
compiled.eval(RDF.__call__(2)).value; // 17
compiled.__call__(RDF.__call__(2)).value; // 17
compiled.toString(); // 'CompiledPolynomialFunction((a4*((x)^2)^2+a0))'
```

`CompiledPolynomialFunction<C>` accepts readonly coefficient elements and an optional
algorithm string (default `'binary'`). Empty coefficients raise IndexError. Other algorithms
raise the original errors when a nontrivial gap requires them. `eval(x: C): C` and
`__call__(x: C): C` execute the graph. An internal overload `eval(x: unknown, arithmetic)`
accepts `multiply(a: unknown, b: unknown): unknown` and `add(a: unknown, b: unknown): unknown`
callbacks for heterogeneous scalar coercion. Direct constant-list compilation retains the
original class behavior; ordinary constant polynomials should use `evaluate`.

The FLINT package exports the following array adapters. Coefficients are readonly bigint
arrays in increasing degree order; all denominators are positive bigints. Evaluation pairs
retain native **unreduced** denominators; reconstruct a Rational to canonicalize them.

| API | Result |
|---|---|
| `_fmpz_poly_evaluate_fmpz(a, x)` | bigint |
| `_fmpz_poly_evaluate_fmpq(a, num, den)` | `[bigint, bigint]` |
| `_fmpq_poly_evaluate_fmpz(a, polyDen, x)` | `[bigint, bigint]` |
| `_fmpq_poly_evaluate_fmpq(a, polyDen, num, den)` | `[bigint, bigint]` |
| `_nmod_poly_evaluate_nmod(a, x, p)` | bigint |
| `_fmpz_poly_compose(a, b)` | bigint[] |
| `_fmpq_poly_compose(a, denA, b, denB)` | `[bigint[], bigint]` |
| `_nmod_poly_compose(a, b, p)` | bigint[] |
| `_fmpz_poly_taylor_shift(a, c)` | bigint[] representing `a(x+c)` |

For example, integer evaluation of `[1n, 2n, 3n]` at `2n` returns `17n`; rational evaluation
at `1n/2n` returns `[11n, 4n]`. With polynomial denominator `3n`, the corresponding QQ
adapters return `[17n, 3n]` and `[11n, 12n]`. Composing `[1n, 0n, 1n]` with `[1n, 1n]`
returns `[2n, 2n, 1n]` (also modulo `7n`). Giving those arrays denominators `2n` and `3n`
returns `[[10n, 2n, 1n], 18n]`. Taylor shift of `[1n, 0n, 1n]` by `2n` returns `[5n, 4n, 1n]`.

NTL exports `ZZ_pX_evaluate(a, x, p): bigint` and
`ZZ_pEX_eval(a, x, f, p): bigint[]`. The latter takes readonly arrays of readonly bigint
coefficient arrays, a readonly bigint point vector `x`, monic extension modulus `f`, and
prime characteristic `p`. `ZZ_pX_evaluate([1n, 2n, 3n], 2n, 7n)` returns `3n`;
`ZZ_pEX_eval([[1n], [0n, 1n]], [2n], [1n, 0n, 1n], 7n)` returns `[1n, 2n]`.
The source module `ZZ_pEX.ts` exports this function under the native reserved name `eval`.

### Polynomial matrix evaluation and integer matrix conversion

Matrix evaluation preserves scalar actions, coefficient promotion and sparse compiled powers.
A scalar constant becomes a diagonal matrix. Rectangular inputs require Sage's integer-zero
coercion exception; nonzero scalar addition or incompatible matrix products still raise.

```typescript
import { QQ, IntegerMatrixFromEntries } from 'sagemath-ts';
import { PolynomialRing } from 'sagemath-ts/rings/polynomial';
import { integer_to_real_double_dense } from 'sagemath-ts/matrix';
const matrixPolynomial = new PolynomialRing(QQ, 'x').__call__([1n, 0n, 1n]);
const matrixPoint = IntegerMatrixFromEntries([[1n, 2n], [3n, 4n]]);
const matrixValue = matrixPolynomial.evaluate(matrixPoint);
matrixValue.get(0, 0).toString(); // '8'
matrixValue.get(1, 1).toString(); // '23'
const realMatrix = integer_to_real_double_dense(IntegerMatrixFromEntries([[9007199254740995n]]));
realMatrix.get(0, 0).value; // 9007199254740994
```

`integer_to_real_double_dense(a: IntegerMatrix): Matrix<RealDoubleElement>` mirrors
`sage.matrix.change_ring` and delegates each entry to FLINT. The FLINT package exports
`fmpz_get_d(value: bigint): number`: it truncates towards zero to 53 significant bits and
returns signed infinity on overflow. This is different from scalar RDF integer conversion.
`fmpz_get_d(9007199254740995n)` returns `9007199254740994`.

### Binary matrix multiplication and evaluation

`Matrix_mod2_dense.mul(B)` delegates to M4RI. Explicit original methods are
`_multiply_classical(B)`, `_multiply_m4rm(B, k: IntegerLike)` and
`_multiply_strassen(B, cutoff: IntegerLike)`. Use `0n` for automatic tuning.
M4RM clamps nonzero k to 2–8; Strassen normalizes its cutoff to 64-bit word boundaries.
Native small-cutoff faults raise RuntimeError. Argument lists and mixed scalar/matrix
inputs to polynomial `evaluate` retain the possible matrix result in the return type.

```typescript
import { Matrix_mod2_dense } from 'sagemath-ts/matrix';
import { GF } from 'sagemath-ts';
import { PolynomialRing } from 'sagemath-ts/rings/polynomial';
const binaryPoint = new Matrix_mod2_dense(2, 2, [[1, 1], [0, 1]]);
const binarySquare = binaryPoint.mul(binaryPoint);
binarySquare.row(0); // [1, 0]
binarySquare.row(1); // [0, 1]
binaryPoint._multiply_m4rm(binaryPoint, 0n).row(0); // [1, 0]
binaryPoint._multiply_strassen(binaryPoint, 64n).row(1); // [0, 1]
const binaryPolynomial = new PolynomialRing(GF(2n), 'x').__call__([1n, 1n]);
const binaryValue = binaryPolynomial.evaluate([binaryPoint]);
(binaryValue as Matrix_mod2_dense).row(0); // [0, 1]
```

The dependency `@sagemath-ts/m4ri-ts` exports `mzd_t` with readonly `nrows`, `ncols`,
and packed bigint `rows`. Its array adapters are `mzd_init(nrows, ncols, rows?)`,
`mzd_add(A, B)`, `mzd_init_window(A, lowRow, lowCol, highRow, highCol)`,
`mzd_mul_naive(A, B)`, `mzd_mul_m4rm(A, B, k = 0)`, `mzd_mul(A, B, cutoff = 0)` and
`mzd_make_table(A, startRow, k)`. A window starts at a column divisible by 64 and owns
its copied rows. Table output has Gray-ordered `rows` and inverse `lookup` arrays.

For `A = mzd_init(2, 2, [3n, 2n])`, all three product adapters return rows `[1n, 2n]`
when squaring A. `mzd_add(A, A).rows` is `[0n, 0n]`.
`mzd_make_table(A, 0, 2)` returns rows `[0n, 3n, 1n, 2n]` and lookup `[0, 1, 3, 2]`.

### Binary matrix indices and slices

`Matrix_mod2_dense.get(row, col)`, `set(row, col, value)`, `row(index, from_list?)` and
`submatrix(row = 0, col = 0, nrows = -1, ncols = -1)` accept `IntegerLike`, `Rational`,
boolean and the existing integer-valued number indices. Reads and writes support negative
entry indices. Any negative slice size means the remaining extent; slice starts must be
nonnegative. Rows and entries raise `IndexError` outside their bounds; invalid slice
extents raise the original `TypeError`. C-index conversion overflow is also preserved.

Rational indices follow the original conversion distinction: `get` truncates them, while
`set`, `row` and `submatrix` require denominator one. `set` still accepts number/boolean
entry values. Row arrays are copies; native vector/cache identity is not exposed.

```typescript
import { Integer, Rational } from 'sagemath-ts';
import { Matrix_mod2_dense } from 'sagemath-ts/matrix';
const indexed = new Matrix_mod2_dense(2, 3, [[0, 1, 0], [1, 0, 1]]);
indexed.get(-1n, new Integer(-1n)); // 1
indexed.row(new Rational(1n)); // [1, 0, 1]
indexed.submatrix(0n, 1n, -2n, -1n).row(1n); // [0, 1]
indexed.get(new Rational(1n, 2n), 1n); // 1 (row index truncates to zero)
indexed.set(-1n, -1n, 0);
indexed.get(1n, 2n); // 0
```

The M4RI dependency also exports `mzd_submatrix(A, startRow, startCol, endRow, endCol)`.
Unlike `mzd_init_window`, this copied slice accepts unaligned column starts. For
`A = mzd_init(1, 5, [22n])`, `mzd_submatrix(A, 0, 1, 1, 4).rows` is `[3n]`.

### Exact binary density and decimal real literals

`Matrix_mod2_dense.density(approx?: false)` returns `Rational | bigint` (`0n` for empty
matrices). `density(true)` returns `RealLiteral` using original M4RI sampling. This return
type changed in 4.0.0. Approximate density of a positive-row, zero-column matrix raises
`RuntimeError`, translating the original native fault.

```typescript
import { Matrix_mod2_dense } from 'sagemath-ts/matrix';
import { create_RealNumber } from 'sagemath-ts/rings';
const densityMatrix = new Matrix_mod2_dense(1, 3, [[1, 0, 0]]);
densityMatrix.density().toString(); // '1/3'
densityMatrix.density(true).toString(); // '0.333333333333333'
const exactLiteral = create_RealNumber('9007199254740993');
exactLiteral.exact_rational(); // [9007199254740993n, 1n]
exactLiteral.toNumber(); // 9007199254740992
exactLiteral.numerical_approx(128).exact_rational(); // [9007199254740993n, 1n]
```

The `sagemath-ts/rings` exports include `RealLiteral(parent: RealField, text: string, base = 10)`
and `create_RealNumber(value, {base?, pad?, rnd?, min_prec?})`. Factory defaults are base 10,
pad 0, rounding `'RNDN'`, minimum precision 53. Literals retain readonly `literal` and `base`.
`numerical_approx(prec?, digits?, algorithm?)` reparses the text and returns `RealNumber`;
precision takes precedence over digits, the default is 53 bits, and algorithm is ignored.
Exact conversion supports decimal nearest-even rounding and precision through 4096 bits.
Input significant digits and decimal scale magnitude are bounded by 4096. Inherited real
arithmetic still has the existing binary64 limitations; this is not an arbitrary-precision
arithmetic replacement.

The dependency `@sagemath-ts/mpfr-ts` exports `mpfr_init2(precision)`,
`mpfr_set_str(value, text, base = 10, mode = 'RNDN')`, `mpfr_get_d(value, mode = 'RNDN')`,
`mpfr_get_str(base, digits, value, mode = 'RNDN')`, and types `mpfr_t`, `mpfr_rnd_t`.
The string setter mutates its destination and returns 0 for complete parsing, −1 otherwise;
nonempty partial parses can still update the value. String output is a digits/exponent
pair with at most 4096 digits; zero requests enough digits to round-trip. Other bases,
rounding modes and larger domains fail explicitly. `m4ri-ts` additionally exports
`mzd_density(A, resolution = 0)`, where zero chooses the native automatic stride.

### Binary elimination, cached forms and mutability

`Matrix_mod2_dense.rank(algorithm = 'ple')` accepts `'ple'` or `'m4ri'` and caches the rank.
`echelonize(algorithm = 'heuristic', cutoff?, reduced = true, options?: {k?: IntegerLike})`
mutates the matrix and accepts `'heuristic'`, `'m4ri'`, `'pluq'` or `'classical'`. Classical
always computes reduced form. The k option applies only to M4RI; cutoff is retained for
signature compatibility. Empty matrices return `this`; other successful calls return
`undefined`. The return type is `void | this` as of 5.0.0.

`echelon_form(algorithm = 'default', cutoff?, reduced = true, options?: {k?: IntegerLike})`
returns and caches an immutable form. Cached results bypass later algorithm validation,
and a cached unreduced form stays unreduced. Copies are mutable. The original mutability
methods are `set_immutable()`, `is_immutable()`, `is_mutable()` and `_clear_cache()`.
Mutation checks invalidate computed properties; the original `doubly_lexical_ordering(true)`
retains its cache and rejects immutable inputs with `TypeError`.

```typescript
import { Matrix_mod2_dense } from 'sagemath-ts/matrix';
const eliminationMatrix = new Matrix_mod2_dense(2, 2, [[1, 1], [0, 1]]);
const cachedForm = eliminationMatrix.echelon_form('m4ri', 0, false, { k: 1n });
cachedForm.row(0); // [1, 1]
cachedForm.is_immutable(); // true
eliminationMatrix.echelon_form('bogus') === cachedForm; // true
eliminationMatrix.rank('m4ri'); // 2
eliminationMatrix.rank('bogus'); // 2 (cached)
eliminationMatrix.set(0, 0, 0);
eliminationMatrix.rank(); // 1
cachedForm.copy().is_mutable(); // true
```

The M4RI dependency exports `mzd_echelonize(A, full = true)`,
`mzd_echelonize_m4ri(A, full = true, k = 0)` and `mzd_echelonize_pluq(A, full = true)`.
Each returns `{matrix, rank}` with an owned result. `mzd_ple(A, cutoff = 0)` and
`mzd_pluq(A, cutoff = 0)` return `{matrix, rank, P, Q}`; P/Q are native transposition lists.
`mzd_trsm_upper_left(A, B, cutoff = 0)` and `mzd_trsm_lower_left(A, B, cutoff = 0)` return
owned solutions for unit triangular systems. Matrix adapter arguments are immutable values;
none of these dependency calls modifies A or B.

### Standalone binary factorizations

`ple(A, algorithm = 'standard', param = 0)` and
`pluq(A, algorithm = 'standard', param = 0)` are exported from `sagemath-ts/matrix`.
Both return `[LU, P, Q]`, where LU is a mutable copy and P/Q are native transposition
lists. PLE accepts `'standard'`, `'russian'`, `'naive'`; PLUQ accepts `'standard'`,
`'mmpf'`, `'naive'`. Every algorithm delegates to its M4RI implementation.

The parameter accepts `IntegerLike`, the existing number inputs, `Rational` and boolean.
Numbers and rationals truncate toward zero before signed C-int bounds checking. Conversion
precedes algorithm validation. Standard uses a multiplication cutoff, Russian/MMPF uses k,
and naive ignores the converted parameter. The supported explicit Russian table range is
1–16; zero chooses automatic tuning. Nonempty matrices reject other table sizes with `RangeError` in the dependency.
Standard factorization of a positive-row, zero-column matrix raises `RuntimeError`, safely
translating the original native fault.

```typescript
import { Matrix_mod2_dense, ple, pluq } from 'sagemath-ts/matrix';
import { Rational } from 'sagemath-ts';
const factorMatrix = new Matrix_mod2_dense(4, 4, [
  [0, 1, 0, 1], [0, 1, 1, 1], [0, 0, 0, 1], [0, 1, 1, 0]
]);
ple(factorMatrix, 'russian', new Rational(3n, 2n))[0].row(0); // [1, 0, 0, 1]
pluq(factorMatrix, 'naive', 0n)[0].row(0); // [1, 0, 1, 0]
```

The M4RI dependency also exports `_mzd_ple_naive(A)`, `_mzd_pluq_naive(A)`,
`_mzd_ple_russian(A, k = 0)` and `_mzd_pluq_russian(A, k = 0)`. They return owned
`{matrix, rank, P, Q}` results. Standard PLE/PLUQ preserve the native zero-column fault
as `RuntimeError`; explicit naive/Russian empty cases return their native zero rank.


### Binary swaps, permutations and inverse

`swap_rows(i, j)` and `swap_columns(i, j)` accept bigint, Integer, integral Rational,
boolean and integer-valued number indices. Both indices convert before immutability
and bounds checks. Negative swap indices are invalid.
`permute_rows(images)` and `permute_columns(images)` accept 1-based number image lists.
Their degree need not equal the matrix size: omitted points and outside fixed points
are unchanged, while moved outside points raise IndexError, possibly after earlier cycles.
`inverse()` requires a square nonsingular matrix and returns a mutable copy through M4RI.
The dependency exports `mzd_inv_m4ri(A, k = 0)`; its native k is ignored and direct callers
must check nonsingularity themselves.

```typescript
import { Matrix_mod2_dense } from 'sagemath-ts/matrix';
const swapExample = new Matrix_mod2_dense(3, 3, [[1, 1, 0], [0, 1, 1], [0, 0, 1]]);
swapExample.swap_rows(0n, 1n);
swapExample.permute_rows([2, 1]);
swapExample.permute_columns([1, 2, 3, 4]);
swapExample.inverse().row(0); // [1, 1, 1]
```


### Binary linear systems and kernel options

`Matrix_mod2_dense.solve_right(B, check = true)` solves `A * X = B` for binary matrix B.
It supports square, rectangular and singular A. The row counts must match; a checked
inconsistent system raises ValueError. `check = false` skips consistency checking.
Results have `A.ncols` rows and `B.ncols` columns.

`right_kernel_matrix({basis, algorithm, proof}?)` returns kernel vectors as rows.
Basis choices are `'default'`/`'echelon'`, `'pivot'` and `'computed'`.
Algorithms are `'default'`/`'pluq'` and `'generic'`. The original other algorithm names,
`'LLL'` basis and non-null proof flags raise their field-specific errors over GF(2).
A zero-row input returns an immutable identity kernel.

```typescript
import { Matrix_mod2_dense } from 'sagemath-ts/matrix';
const system = new Matrix_mod2_dense(2, 3, [[1, 0, 1], [0, 1, 1]]);
const rhs = new Matrix_mod2_dense(2, 1, [[1], [0]]);
system.solve_right(rhs).list(); // [1, 0, 0]
system.right_kernel_matrix({basis: 'pivot', algorithm: 'pluq'}).row(0); // [1, 1, 1]
```

The M4RI dependency exports `mzd_solve_left(A, B, cutoff = 0, check = true)` returning
`{matrix, rhs, status}` and `mzd_kernel_left_pluq(A, cutoff = 0)` returning
`{matrix, kernel}`. These contain owned copies of native mutated buffers; status is
zero on success and -1 on inconsistency, and a trivial kernel is null.


### Binary columns and native basic operations

`Matrix_mod2_dense.columns(copy = true)` returns a shallow outer-array copy by default;
its column vectors are immutable and shared with the cache. `columns(false)` returns the
cached outer array itself. Copy a column with `.slice()` before changing its entries.
Matrix mutations invalidate the cache. `add`, `sub`, `transpose` and `augment` return
mutable matrices and delegate their bit operations to M4RI.

```typescript
import { Matrix_mod2_dense } from 'sagemath-ts/matrix';
const columnExample = new Matrix_mod2_dense(2, 2, [[1, 0], [1, 1]]);
const cachedColumns = columnExample.columns(false);
columnExample.columns()[0] === cachedColumns[0]; // true
const editableColumn = cachedColumns[0].slice();
editableColumn[0] = 0;
columnExample.columns()[0]; // [1, 1]
columnExample.transpose().row(0); // [1, 1]
```

The M4RI dependency exports `mzd_transpose(A)` and `mzd_concat(A, B)`, returning owned
packed matrices. `mzd_concat` requires matching row counts. Sage safely handles empty
transposes; direct native rectangular-empty transposes raise the original abort exception.


### Binary formatting and subdivisions

`str(mapping?, zero?, plus_one?, minus_one?)` accepts a record keyed by 0/1 or a callback
from a number bit to a string. Zero/one keywords override and mutate record mappings;
callbacks take precedence over those keywords. Minus-one is ignored over GF(2). Entries
align to a shared width, counting Unicode code points.

`subdivide(row_lines?, col_lines?)` accepts individual integer-coercible values or arrays.
Omitted/null arguments clear that axis. Lists are sorted and duplicate lines remain.
`subdivisions()` and `get_subdivisions()` return copied `[bigint[], bigint[]]` pairs.
Copies, negation and transpose preserve subdivisions. `augment(B, true)` inserts the
boundary and combines compatible subdivisions; arithmetic and extracted slices start fresh.

```typescript
import { Matrix_mod2_dense } from 'sagemath-ts/matrix';
const formatted = new Matrix_mod2_dense(2, 3, [[0, 1, 0], [1, 0, 1]]);
const mapping = {0: 'zero', 1: 'x'};
formatted.str(mapping, '.', 'one'); // '[  . one   .]\n[one   . one]'
mapping; // {0: '.', 1: 'one'}
formatted.subdivide([1n], [1n]);
formatted.augment(formatted, true).subdivisions(); // [[1n], [1n, 3n, 4n]]
formatted.transpose().get_subdivisions(); // [[1n], [1n]]
```


### Binary construction and cached rows

`new Matrix_mod2_dense(nrows, ncols, entries?)` accepts number, bigint, Integer, Rational
or boolean dimensions. Dimensions truncate as positional MatrixSpace integer inputs do.
Entries may be omitted/null, a scalar, a flat array or a nested row array. Nonzero scalars
require a square matrix; array lengths must match exactly. Entries and `set(i, j, value)`
convert through GF(2), including exact large integers and rationals with odd denominators.

`row(i)` returns a new mutable array. `row(i, true)` returns the cached immutable row;
repeated calls share that row until a matrix mutation clears the cache.

```typescript
import { Matrix_mod2_dense } from 'sagemath-ts/matrix';
import { Rational } from 'sagemath-ts';
const constructed = new Matrix_mod2_dense(2n, 2n, [1n, 0n, 0n, 1n]);
constructed.set(0, 1, new Rational(5n, 3n));
constructed.list(); // [1, 1, 0, 1]
const cachedRow = constructed.row(0, true);
constructed.row(0, true) === cachedRow; // true
constructed.row(0) === cachedRow; // false
constructed.set(0, 0, 0n);
constructed.row(0, true) === cachedRow; // false
```


### Seeded binary matrices

`randomize(density = 1, nonzero = false)` uses the shared Sage random state. Full density
with `nonzero = false` gives uniform random bits; lower densities update selected positions
with replacement and leave other entries unchanged. With `nonzero = true`, selected entries
become one. Empty matrices and nonpositive densities return without changing cached state.

`zero_matrix_gf2`, `identity_matrix_gf2` and `random_matrix_gf2` accept the constructor's
number/IntegerLike/Rational/boolean dimensions. Omitting random-factory density gives
uniform bits; supplying density selects nonzero entries. `matrix_gf2_from_entries` accepts
nested arrays of constructor-compatible coefficients and validates consistent row lengths.

```typescript
import { random_matrix_gf2 } from 'sagemath-ts/matrix';
import { set_random_seed } from 'sagemath-ts/misc/randstate';
set_random_seed(42n);
const firstRandom = random_matrix_gf2(2n, 65n);
set_random_seed(42n);
firstRandom.eq(random_matrix_gf2(2n, 65n)); // true
random_matrix_gf2(2, 2, 1).list(); // [1, 1, 1, 1]
```


### Binary density inputs and PNG data errors

Density also accepts bigint, Integer, Rational, RealNumber/RealLiteral, boolean, prime/modular
field elements, numeric strings and Uint8Array text bytes. Conversion follows Python float
rules. Oversized bigint raises OverflowError; Sage Integer/Rational can convert to infinity,
which density clamps to one. `randomize(null)` raises TypeError on nonempty matrices, while
null factory density selects its default behavior. Empty matrices skip density conversion.

```typescript
import { Matrix_mod2_dense, random_matrix_gf2, to_png_data } from 'sagemath-ts/matrix';
const densityInput = new Matrix_mod2_dense(1, 2, [[1, 0]]);
densityInput.randomize('0');
densityInput.list(); // [1, 0]
random_matrix_gf2(1, 2, new TextEncoder().encode('1')).list(); // [1, 1]
// densityInput.randomize(null) throws TypeError.
// to_png_data(new Matrix_mod2_dense(0, 0)) throws TypeError:
// cannot write image with dimensions 0 x 0
```


### Lattice Gaussian defaults and deferred centers

`new DiscreteGaussianDistributionLatticeSampler(basis, options?)` and
`DiscreteGaussianLattice(basis, sigma = 1, c?, tau = 6n)` default sigma to 1.
`options.c`/`set_c` accept a vector, `0`, `0n`, or `null`; omitted constructor centers
mean zero. Explicit null defers precomputation: `c(): Rational[] | null` and
`cNumeric(): number[] | null`. Center changes preserve the original sampler cache behavior.
Use `sampleExact()` for Rational samples; `call()`/`sample()` require integral results.
The low-level `_call(): Rational[] | 0` returns zero on a basis with no rows.

```ts
import { DiscreteGaussianDistributionLatticeSampler as DGL,
         DiscreteGaussianLattice } from 'sagemath-ts/stats';
const gaussian = new DGL([[1]]);
gaussian.sigma();                 // 1
DiscreteGaussianLattice([[1]]).sigma(); // 1
gaussian.set_c(null);
gaussian.c();                     // null
gaussian.cNumeric();              // null
new DGL([])._call();              // 0
new DGL([]).sampleExact();        // []
```


### Gaussian parameter strings and MPFR double conversion

Sampler `repr()` follows Sage's numeric formatting, including nearest-even halfway
rounding and negative zero. With `@sagemath-ts/mpfr-ts` installed as a direct dependency, it exports
`mpfr_set_d(destination, value, rounding = 'RNDN'): number`: it mutates an `mpfr_t` and
returns zero for exact conversion, negative for rounding down, positive for rounding up.
The port supports nearest-even mode and precisions 1..4096.

```ts
import { DiscreteGaussianDistributionIntegerSampler as DGI,
         DiscreteGaussianDistributionLatticeSampler as DGL } from 'sagemath-ts/stats';
import { mpfr_init2, mpfr_set_d, mpfr_get_d } from '@sagemath-ts/mpfr-ts';
new DGI({ sigma: 1/128, c: -0 }).repr();
// Discrete Gaussian sampler over the Integers with sigma = 0.007812 and c = -0.000000
new DGL([[2]], { sigma: 1e6 }).repr().split(', c=')[0];
// Discrete Gaussian sampler with Gaussian parameter σ = 1.00000000000000e6
const rounded = mpfr_init2(2);
mpfr_set_d(rounded, 1.25); // -1
mpfr_get_d(rounded);      // 1
```


### Integer Gaussian centers and tail cutoffs

Integer centers are rounded to the sampler's 53-bit real field. Samples and exposed
bounds retain exact integer offsets around that rounded center. A fractional numeric
`tau` raises `TypeError`, as in Sage's final integer conversion.

```ts
import { DiscreteGaussianDistributionIntegerSampler as DGI } from 'sagemath-ts/stats';
const shifted = new DGI({ sigma: 2, c: 9007199254740993n, tau: 1 });
shifted.c;          // 9007199254740992
shifted.lowerBound; // 9007199254740990n
shifted.upperBound; // 9007199254740994n
```


### Native real numeric construction

Nearest-even real construction rounds doubles and integers to the field precision, parses
strings through MPFR, and retains the exact converted value for formatting and rational
extraction. Conversion supports precisions 1..4096; ordinary arithmetic retains its
existing binary64 limitations. `mpfr_set_z(destination, value: bigint, rounding = 'RNDN')`
returns zero for exact conversion, negative for rounding down, positive for rounding up.

```ts
import { RealField } from 'sagemath-ts/rings';
import { mpfr_init2, mpfr_set_z, mpfr_get_d } from '@sagemath-ts/mpfr-ts';
new RealField(2).__call__(1.25).exact_rational(); // [1n, 1n]
new RealField(100).__call__(9007199254740993n).exact_rational();
// [9007199254740993n, 1n]
const integer = mpfr_init2(2);
mpfr_set_z(integer, 7n); // 1
mpfr_get_d(integer);    // 8
```


### Native real predicates and integer rounding

Native-backed reals retain finite/zero/sign distinctions beyond the binary64 range.
`floor`, `ceil`, `round` (ties away from zero) and `trunc` return exact BigInts; infinity
and NaN raise Sage's `ValueError` messages.

The MPFR dependency also exports `mpfr_nan_p`, `mpfr_inf_p`, `mpfr_number_p`,
`mpfr_integer_p`, `mpfr_sgn` and `mpfr_cmp_si(value, integer: bigint)`.
`mpfr_rint(destination, source, mode = 'RNDN')` and its two-argument wrappers
`mpfr_roundeven`, `mpfr_round`, `mpfr_trunc`, `mpfr_ceil`, `mpfr_floor` support in-place calls.
`mpfr_get_z(source, mode = 'RNDN')` returns `[integer, inexactStatus]`. For these operations,
status is zero when exact; its sign gives the rounding direction, with magnitude 2 for
fractional input and 1 for an integer losing destination precision. Integer rounding
supports `RNDN`, `RNDZ`, `RNDU`, `RNDD`, `RNDA`; `mpfr_rint` also accepts `RNDNA` for ties away.

```ts
import { RealField } from 'sagemath-ts/rings';
import { mpfr_init2, mpfr_set_d, mpfr_floor, mpfr_get_z } from '@sagemath-ts/mpfr-ts';
const R = new RealField(128);
R.__call__('1e400').is_infinity();          // false
R.__call__('-1e-400').sign();               // -1
R.__call__('9007199254740993.25').floor();  // 9007199254740993n
const source = mpfr_init2(53), rounded = mpfr_init2(2);
mpfr_set_d(source, 1.5);
mpfr_floor(rounded, source); // -2
mpfr_get_z(source);          // [2n, 2]
```


### Native fractional parts and real comparisons

`frac()` retains the native fractional value and signed zero. `cmp` and `equals` compare
two real fields at their lower precision; literals reparse their text during coercion.
Comparisons with a JavaScript number follow Python-float coercion: real parents above
53 bits compare in binary64, while lower-precision parents receive the float.

The MPFR dependency exports `mpfr_set(destination, source, mode = 'RNDN')`,
`mpfr_frac(destination, source, mode = 'RNDN')` and `mpfr_cmp(left, right)`. Copy/fraction
support nearest-even mode and return native inexact status. Primitive comparison returns
−1, 0 or 1, with 0 for NaN; real `equals` handles NaN separately.

```ts
import { RealField } from 'sagemath-ts/rings';
import { mpfr_init2, mpfr_set_str, mpfr_frac, mpfr_get_d } from '@sagemath-ts/mpfr-ts';
const R = new RealField(128);
R.__call__('9007199254740993.25').frac().exact_rational(); // [1n, 4n]
const nearOne = R.__call__('1.00000000000000000001');
nearOne.equals(R.__call__(1)); // false
nearOne.equals(1);             // true, after Python-float coercion
const source = mpfr_init2(17), fraction = mpfr_init2(1);
mpfr_set_str(source, '1.75');
mpfr_frac(fraction, source); // 2, the native halfway status
mpfr_get_d(fraction);        // 1
```


### Copying Gaussian parameters

`sampler.withOptions(overrides)` constructs a sampler with selected parameter overrides.
Undefined values retain the original parameter; supplied values use constructor validation.
The copy has independent sampling caches. `DiscreteGaussianOptions.c` and `.tau` accept
`IntegerLike | number`; numeric tail cutoffs must be integral. This also allows mixed
integer/numeric options in the constructor and `withOptions`.

```ts
import { DiscreteGaussianDistributionIntegerSampler as DGI } from 'sagemath-ts/stats';
const original = new DGI({ sigma: 2 });
original.withOptions({ c: 0.5 }).c;               // 0.5
original.withOptions({ c: 9007199254740993n }).c; // 9007199254740992
original.c;                                      // 0
```

### Algebraic dependencies and native real arithmetic

`algebraic_dependency(z, degree, options?)` (alias `algdep`) accepts a number,
IntegerLike, Rational or RealNumber and an IntegerLike degree. It returns ascending
bigint coefficients; a supplied height bound may return `null` for exact/real inputs.
Numbers follow Python-float/RDF dispatch to PARI: precision hints are ignored, and
nonzero height bounds/proofs are unsupported. RealNumber inputs preserve their stored
precision, apply `known_bits`, `use_bits`, `known_digits`, `use_digits` and `proof`
options, and select the closest irreducible factor. All precision/height options are bigint.
`RealNumber.algebraic_dependency(n)` and `.algdep(n)` retain their number degree parameter.

```typescript
import { algebraic_dependency } from 'sagemath-ts/arith';
import { RealField } from 'sagemath-ts/rings';
const R = new RealField(128);
algebraic_dependency(R.__call__('1.5'), 2n); // [-3n, 2n]
algebraic_dependency(5n, 0n, { height_bound: 5n }); // null
algebraic_dependency(Math.sqrt(2), 3n); // [-2n, 0n, 1n]
```

The MPFR dependency exports `mpfr_add(destination, left, right, mode = 'RNDN')`
and `mpfr_mul(destination, left, right, mode = 'RNDN')`. They mutate the destination
at its precision, support operand aliasing, and return native inexact status. Only
nearest-even mode is implemented; other modes raise NotImplementedError. These
primitives do not change the existing limits of general RealNumber arithmetic.

```typescript
import { mpfr_init2, mpfr_set_str, mpfr_add, mpfr_mul, mpfr_get_d } from '@sagemath-ts/mpfr-ts';
const a = mpfr_init2(100), b = mpfr_init2(100), result = mpfr_init2(2);
mpfr_set_str(a, '1.5'); mpfr_set_str(b, '1.5');
mpfr_mul(result, a, b); // -1 (2.25 rounds down to 2 at precision 2)
mpfr_get_d(result); // 2
mpfr_add(result, a, b); // 0 (exact)
mpfr_get_d(result); // 3
```

The PARI package also exports `algdep(z: number, degree: bigint): bigint[]`, returning
its raw ascending polynomial coefficients (before Sage's irreducible-factor selection).
It preserves 64-bit C-long argument limits and native conversion/degree errors.


Prime/power traversal and ranges accept IntegerLike inputs: `next_prime_power`,
`previous_prime_power`, `next_probable_prime`, `is_power_of_two`,
`is_pseudoprime_power`, `prime_powers(start, stop?)` and `eratosthenes(n)`.
`is_pseudoprime_power(n, true)` returns `[base, exponent]`, or `[n, 0n]` on failure.
The free `previous_prime` function raises `ValueError('no previous prime')` below 3;
`Integer.previous_prime()` retains its original `ValueError('no prime less than 2')`.

```typescript
import { Integer } from 'sagemath-ts';
import { next_prime_power, is_pseudoprime_power, prime_powers, eratosthenes } from 'sagemath-ts/arith';
next_prime_power(new Integer(255n)); // 256n
is_pseudoprime_power(new Integer(1000003n ** 6n), true); // [1000003n, 6n]
prime_powers(new Integer(8n), new Integer(17n)); // [8n, 9n, 11n, 13n, 16n]
eratosthenes(new Integer(11n)); // [2n, 3n, 5n, 7n, 11n]
```

The PARI dependency exports `nextprime(n: bigint): bigint` and
`precprime(n: bigint): bigint` with inclusive bounds, plus
`isprimepower(n: bigint): [bigint, number] | null`. The latter uses its existing
BPSW backend and returns a number exponent. `precprime` returns zero below two;
`nextprime` returns two for inputs below two.

`hilbert_symbol(a: RationalLike, b: RationalLike, p: IntegerLike, algorithm = 'pari')`
accepts integer or rational coefficients and a prime place (`-1n` for the real place).
Algorithms are `'pari'`, `'direct'` and `'all'`; the last compares both computations.
Invalid algorithms raise ValueError even for zero coefficients. `hilbert_conductor(a, b)`
and `hilbert_conductor_inverse(d)` accept IntegerLike inputs. The PARI dependency exports
`hilbert(a: bigint, b: bigint, p = 0n): number` for integer coefficients, where zero denotes
the real place and the result is -1, 0 or 1.

```typescript
import { Integer, Rational } from 'sagemath-ts';
import { hilbert_symbol, hilbert_conductor, hilbert_conductor_inverse } from 'sagemath-ts/arith';
hilbert_symbol(new Rational(1n, 2n), new Rational(1n, 3n), 2n, 'all'); // -1n
hilbert_conductor(new Integer(-3n), new Integer(-17n)); // 17n
hilbert_conductor_inverse(new Integer(30n)); // [-3n, -10n]
```

The internal PARI module `src/gen2.ts` exports
`Z_pvalrem(n: bigint, p: bigint): [number, bigint]` for nonzero n and p>1, returning the
valuation and signed unit. Invalid preconditions raise RangeError. It is used by the
Hilbert backend and is not re-exported from the package root.

The same module exports `cmp_universal(x: PariFfelt, y: PariFfelt): number` for
reduced finite-field elements at one native variable index. It returns -1, 0 or 1
in PARI's representation order. Provide coefficient arrays and `definingPoly` for
extensions; a scalar value without a defining polynomial denotes the prime-field
modulus X. This is the finite-field subset, not a comparator for arbitrary GENs.

```ts
import { cmp_universal } from '@sagemath-ts/parigp-ts/src/gen2.js';
import { PariType } from '@sagemath-ts/parigp-ts/src/types.js';
const ff = { type: PariType.t_FFELT as const, p: 7n, degree: 2,
  definingPoly: [3n, 6n, 1n] };
cmp_universal({ ...ff, value: [1n, 5n] }, { ...ff, value: [2n, 1n] }); // -1
cmp_universal({ ...ff, value: [2n, 1n] }, { ...ff, value: [1n, 5n] }); // 1
cmp_universal({ ...ff, value: [0n] }, { ...ff, value: [] }); // 0
```


`quadratic_residues`, `fundamental_discriminant`, `squarefree_divisors`, `odd_part`
and both arguments of `prime_to_m_part` accept IntegerLike. `continuant(v, n?)`
accepts IntegerLike entries and an optional IntegerLike order, returning bigint.
Orders above the list length are truncated; negative orders return the first entry
(and raise IndexError for an empty list). `fundamental_discriminant(0n)` is zero;
`quadratic_residues(0n)` raises ZeroDivisionError. Squarefree divisors retain Sage's
subset order and support short lazy prefixes even with more than 31 prime factors.

```typescript
import { Integer } from 'sagemath-ts';
import { continuant, quadratic_residues, fundamental_discriminant, squarefree_divisors, odd_part, prime_to_m_part } from 'sagemath-ts/arith';
continuant([new Integer(7n), 2n], -1n); // 7n
quadratic_residues(new Integer(-8n)); // [0n, 1n, 4n]
fundamental_discriminant(0n); // 0n
[...squarefree_divisors(30n)]; // [1n, 2n, 3n, 6n, 5n, 10n, 15n, 30n]
odd_part(new Integer(0n)); // 0n
prime_to_m_part(new Integer(-240n), new Integer(-3n)); // -80n
```

`get_gcd(order: IntegerLike)` and `get_inverse_mod(order: IntegerLike)` return callbacks
accepting two IntegerLike arguments and returning bigint. They select Sage's 32-bit
implementation for `order <= 46340`, its 64-bit implementation through 2147483647, and
the generic function above that. Bounded callbacks enforce native argument widths and
preserve their distinct inverse errors and negative-modulus conventions. See
DEVIATIONS.md for native undefined arithmetic and the Cython -1 error sentinel.

```ts
import { Integer } from 'sagemath-ts';
import { get_gcd, get_inverse_mod, gcd } from 'sagemath-ts/arith';
get_gcd(new Integer(4000n))(new Integer(-18n), 24n); // 6n
get_inverse_mod(6000n)(2n, -3n); // 2n
get_inverse_mod(600000n)(2n, -3n); // -4n
get_gcd(4000000000n) === gcd; // true
```

`crt` and its `CRT` alias also accept `(residues: IntegerLike[], moduli: IntegerLike[])`.
`CRT_list(residues, moduli)` returns bigint and combines pairs in a balanced tree;
a singleton returns the original integer value without reducing it. The one-argument
`CRT_list(modularElements)` form accepts IntegerMod and both prime-field element classes,
returns a modular element, and preserves the input object for a singleton. Empty input
returns zero modulo one. Input arrays are not mutated.

`CRT_basis(moduli: IntegerLike[], require_coprime_moduli = true)` returns a bigint array,
or `[basis, coprime]` with the option false, except that an empty input always returns `[]`.
Negative moduli retain Sage's signed remainder. The bundled source can retain extra
prefix coefficients when entering its non-coprime fallback; those are preserved.
`CRT_vectors(X: IntegerLike[][], moduli: IntegerLike[])` uses that same basis and preserves
Sage's empty-input, ragged-vector and inconsistent-system errors.

```ts
import { Mod } from 'sagemath-ts';
import { crt, CRT_list, CRT_basis } from 'sagemath-ts/arith';
crt([2n, 3n], [5n, 7n]); // 17n
CRT_list([10n], [3n]); // 10n
const residue = Mod(2n, 3n);
CRT_list([residue]) === residue; // true
CRT_list([residue, Mod(3n, 5n)]).modulus; // 15n
CRT_basis([-3n]); // [-2n]
CRT_basis([7n, 6n, 10n], false); // [[120n, 120n, -140n, 21n], false]
```

`ProductTree` and `prod_with_derivative` are exported from `sagemath-ts/rings`.
`ProductTree(leaves: Iterable<IntegerLike>)` provides `root(): bigint`,
`leaves(): readonly bigint[]`, `remainders(x: IntegerLike): bigint[]`,
`interpolation(xs: readonly IntegerLike[]): bigint`, `__len__()` and iteration.
Its public `layers` contains frozen arrays; interpolation caches CRT bases across calls.
Empty-tree root/interpolation raises AssertionError, and singleton interpolation returns
the supplied value without reduction. These signatures cover the integer domain.

`smooth_part(x: IntegerLike, base: Iterable<IntegerLike> | ProductTree)` returns sorted
`[factor, exponent]` bigint pairs. `coprime_part(x, base)` returns the integer floor quotient
by that factorization. Use pairwise coprime nonunit bases when mathematical divisibility
is required: overlapping bases retain Sage's surprising quotient behavior. Inputs that
would make Sage's extraction loop run forever raise the documented NotImplementedError.

```ts
import { ProductTree, prod_with_derivative } from 'sagemath-ts/rings';
import { smooth_part, coprime_part } from 'sagemath-ts/arith';
const tree = new ProductTree([2n, 3n, 5n]);
tree.root(); // 30n
tree.remainders(17n); // [1n, 2n, 2n]
tree.interpolation([1n, 2n, 2n]); // 17n
smooth_part(240n, tree); // [[2n, 4n], [3n, 1n], [5n, 1n]]
coprime_part(-126n, [6n, 14n]); // -2n
prod_with_derivative([[2n, 1n], [3n, 1n]]); // [6n, 5n]
```

`xlcm(m: IntegerLike, n: IntegerLike)` returns `[multiple, m1, n1]` as bigints and
uses Sage's factorization-free gcd reduction. Signed inputs retain the sign of
`m*n/gcd(m,n)`; a zero first argument raises ZeroDivisionError, while `(m, 0)` with
nonzero m returns zero with a signed unit first factor. `dedekind_psi(N: IntegerLike)`
accepts negative integers and preserves their sign; zero raises the original factorization
error. `carmichael_lambda(n: IntegerLike)` requires a positive integer and returns bigint.

```ts
import { Integer } from 'sagemath-ts';
import { xlcm, dedekind_psi, carmichael_lambda } from 'sagemath-ts/arith';
xlcm(new Integer(-12n), 18n); // [-36n, -4n, 9n]
xlcm(12n, 0n); // [0n, 1n, 0n]
dedekind_psi(new Integer(-6n)); // -12n
carmichael_lambda(new Integer(16n)); // 4n
```

The PARI package exports `eulerphi(n: bigint): bigint` and `numdiv(n: bigint): bigint`.
They accept signed integers; native `eulerphi(0n)` returns `2n`, while `numdiv(0n)` raises
PariError. The Sage `euler_phi` and `number_of_divisors` wrappers delegate to these functions
and preserve Sage's distinct zero checks. PARI factorization and quadratic roots share the
native valuation/unit implementation; large prime powers use its divide-and-conquer path.

```ts
import { eulerphi, numdiv } from '@sagemath-ts/parigp-ts';
eulerphi(-12n); // 4n
eulerphi(0n); // 2n
numdiv(-12n); // 6n
```

The Buchmann submodule's `dbltor(d: number, p?: number): Real` delegates to native
conversion. Omit `p` to retain native allocation (zero has `p: 0, e: -1023`); an
explicit positive working precision applies native `rtor` after rounding the requested
bits up to a multiple of 64. `rtodbl(x: Real | MpReal<bigint>): number` uses the same native conversion
as the shared kernel, including its underflow and overflow rules.


The Buchmann real constructors `itor(n: bigint, p?: number)`, `real_0(p?: number)`
and `real_1(p?: number)` use a default of 64 bits and round requested working bits up
to whole 64-bit words. `setprec(x: Real, p: number)` performs immutable native `rtor`
conversion, including zero-accuracy clamping. `real_0_bit(e: number)` now has the
native one-argument signature and returns a canonical zero with `p: 0`.
`mulur(n: number | bigint, x: Real)` and `divru(x: Real, n: number | bigint)` accept
BigInt to retain the full unsigned-word domain. Real operands use complete native words.

`truncr(x: MpReal): bigint` rejects a nonzero input with `x.e >= x.p` using `PariError`
because its precision cannot determine the integer part. `gcvtoi` still returns the
stored integer approximation with its error exponent in that case.

```typescript
import { gcvtoi, itor, truncr } from '@sagemath-ts/parigp-ts';
const lowPrecisionInteger = itor(1n << 64n, 64);
gcvtoi(lowPrecisionInteger); // [18446744073709551616n, 1]
try { truncr(lowPrecisionInteger); } catch (error) { (error as Error).message; }
// 'precision too low in truncr (precision loss in truncation)'
```


The public `mpexp(x: Real): MpReal<bigint>` (Buchmann `expr`) now delegates to the
native exponential kernel, including zero accuracy and precision growth. Its exponent
is always BigInt. Public `rtodbl` accepts that result directly. Buchmann `sqrtr` now
returns `MpReal | MpComplex`; logarithms and logarithm constants share the native
series/AGM and binary-splitting kernels with qfb.

```typescript
import { itor, mpexp, rtodbl } from '@sagemath-ts/parigp-ts';
mpexp(itor(0n, 64)); // {s: 1, e: 0n, m: 9223372036854775808n, p: 64}
rtodbl(mpexp(itor(0n, 64))); // 1
```


Internal dependency modules expose the native split helpers: `trans1.ts` exports
`abpq_init(n: number): Abpq`, `abpq_sum(n1: number, n2: number, A: Abpq): AbpqResult`
and `constlog2(p: number): MpReal`. `Abpq` has BigInt arrays `a`, `b`, `p`, `q`;
initialize indices `0..n` before summing a nonempty interval `[n1,n2)`.
`AbpqResult` contains BigInts `P`, `Q`, `B`, `T`. The constant cache may retain a
larger precision; `mplog2(p)` converts it to the requested precision.
`trans2.ts` exports `atanhuu(u: bigint, v: bigint, p: number): MpReal` for unsigned
words `0 < u < v` and a positive whole-word bit precision.

The internal `trans1.ts` helper `exp1r_abs(x: MpReal): MpReal<bigint>` computes
`exp(abs(x)) - 1` and retains the full native exponent, including overflow errors.
Its input must be a nonzero normalized real. As with `mpexp`, the result's exponent
is BigInt even for small results; range-reduced internal callers convert only known
small exponents when entering the number-exponent arithmetic kernel.


Rational function fields are available from `sagemath-ts/rings/function_field`.
`FunctionField(k, name)` reuses the same parent for a constant-field object and name;
`['x']` and `'x'` share a key. `K.change_variable_name(name)` returns
`[L, fromL, toL]`, with plain callable isomorphisms in the indicated directions.
For an unchanged name, both maps are the same identity function. The current
structural constant-field interface requires the following cast for `QQ`:

```typescript
import { QQ } from 'sagemath-ts';
import { FunctionField, type ConstantField, type ConstantFieldElement } from 'sagemath-ts/rings/function_field';
const K = FunctionField(QQ as unknown as ConstantField<ConstantFieldElement>, 'x');
const [L, fromL, toL] = K.change_variable_name('y');
toL(K.gen()).toString(); // 'y'
fromL(L.gen()).toString(); // 'x'
L.change_variable_name('x')[0] === K; // true
```

`K.gen(n)` and `K.maximal_order().gen(n)` raise `IndexError` for nonzero indices.
Function-field element comparison follows Sage's selected fraction backend, including
FpT's ordering for odd prime constant fields below 46341 and numerical rational
coefficient ordering. Square-root lists contain one root in characteristic two.

Function-field `valuation` converts the entire polynomial argument, rejects nonintegral
fractions and nonzero constants, and validates the argument even for a zero element.
`K.zero().factor()` and factorization/divisor operations on its finite ideal raise
`ArithmeticError('factorization of 0 is not defined')`, following the bundled Sage source.


`K.field(): FractionField_generic<C>` returns the actual `Frac(k[x])` parent;
`f.element(): FractionElement<C>` returns the underlying fraction object. Their
identities and the fraction's normalization state are preserved. `K.__call__` accepts
fraction elements and rational expressions as strings, using the fraction parent's
coercion rules. `new FunctionFieldElement_rational(K, fraction)` wraps that fraction;
its existing polynomial numerator/denominator constructor remains available.

```typescript
const underlying = K.field();
const reciprocal = K.__call__('1/x');
reciprocal.element().parent === underlying; // true
K.__call__(reciprocal.element()).element() === reciprocal.element(); // true
underlying.__call__('x==x', '1/x').toString(); // 'x'
reciprocal.pow(1n) === reciprocal; // true
```


`K.maximal_order_infinite().basis(): [bigint]` returns `[1n]`, matching the bundled
Sage integer tuple. The finite order's basis consists of function-field elements.
Order construction preserves syntax and division errors; modular inversion preserves
Sage's empty `AssertionError` for nonintegral arguments or ideals. Inversion modulo
the zero ideal succeeds for nonzero constants.

```typescript
const finiteOrder = K.maximal_order();
K.maximal_order_infinite().basis(); // [1n]
K.__call__(2n).inverse_mod(finiteOrder.ideal(K.zero())).toString(); // '1/2'
const ideal = finiteOrder.ideal(K.gen());
finiteOrder.ideal_monoid() === finiteOrder.ideal_monoid(); // true
finiteOrder.ideal_monoid().__call__(ideal) === ideal; // true
```


FLINT exports `n_sqrtmod(a, p)`, `n_is_square(a)`, `n_jacobi_unsigned(a, p)`,
`n_jacobi(a, p)`, `_n_jacobi_unsigned(a, p, parity)`, `n_preinvert_limb(p)` and
`n_powmod2_ui_preinv(a, exponent, p, reciprocal)`. Word values are BigInt; the
Jacobi parity argument is a machine-sized number. Square tests return Boolean,
Jacobi symbols return numbers, and the other results are BigInt.

`n_sqrtmod` expects a reduced residue and a prime modulus fitting an unsigned
64-bit word; zero denotes either the zero root or failure. Jacobi denominators
must be positive odd unsigned words; `n_jacobi` accepts a signed 64-bit numerator.
Modular power expects unsigned words, nonzero modulus and its correct reciprocal.
The port uses exact BigInt residues internally instead of consuming that reciprocal.

```typescript
import { n_sqrtmod, n_jacobi, n_is_square, n_preinvert_limb, n_powmod2_ui_preinv } from '@sagemath-ts/flint-ts';
n_sqrtmod(100n, 101n); // 10n
n_jacobi(-1n, 7n); // -1
n_is_square(4294967295n ** 2n); // true
n_powmod2_ui_preinv(100n, 3n, 101n, n_preinvert_limb(101n)); // 100n
```


FLINT also exports `_nmod_poly_sqrt(f, p): bigint[] | null`,
`_nmod_poly_sqrt_series(f, n, p)`, `_nmod_poly_mulhigh(a, b, start, p)` and
`_nmod_poly_mulhigh_classical(a, b, start, p)`. Series helpers are
`_gr_poly_sqrt_series_basecase(f, n, p)`, `_gr_poly_sqrt_series_newton(f, n, cutoff, p)`,
`_gr_poly_rsqrt_series_basecase(f, n, p)` and `_gr_poly_inv_series_basecase(f, n, p)`.
They return coefficient arrays; lengths and cutoffs are numbers. GR helpers currently
support modular coefficients only. Root series require a nonzero square constant in
an odd prime field; inverse series require an invertible constant. Full polynomial
roots additionally support characteristic two and return `null` for a nonsquare.

NTL exports `ZZX_GCD(a, b, state?: PolynomialProductState): bigint[]` and
`ZZX_SquareFreeDecomp(f, state?: PolynomialProductState): Array<[bigint[], number]>`.
Squarefree input must be primitive with positive leading coefficient;
Sage's integer polynomial method removes signed content before calling it.

```typescript
import { _nmod_poly_sqrt, _gr_poly_sqrt_series_newton } from '@sagemath-ts/flint-ts';
import { ZZX_GCD, ZZX_SquareFreeDecomp } from '@sagemath-ts/ntl-ts';
_nmod_poly_sqrt([4n, 8n, 4n], 17n); // [2n, 2n]
_gr_poly_sqrt_series_newton([1n, 2n, 1n], 8, 2, 17n); // [1n, 1n]
ZZX_GCD([-2n, 0n, 2n], [2n, -4n, 2n]); // [-2n, 2n]
ZZX_SquareFreeDecomp([1n, 4n, 4n]); // [[[1n, 2n], 2]]
```

`Polynomial.is_square()` returns Boolean; `is_square(true)` returns
`[boolean, Polynomial | null]`. Generic `FractionFieldElement` has the same optional
root protocol, returning a fraction root. Both generic fractions and FpT expose
`sqrt(extend = true, all = false)`, returning one fraction or an array. Generic
fractions additionally accept `name: string | null = null`; named nonsquare extensions
are currently unimplemented. FpT exposes `_sqrt_or_None()` and accepts no root argument
on `is_square()`. Function-field `sqrt(all = false)` delegates and wraps each fraction.

```typescript
import { FractionFieldElement } from 'sagemath-ts/rings/fraction_field_element';
const F = K.field(), R = K._ring;
const raw = new FractionFieldElement(F, R.gen(), R.gen(), { reduce: false });
raw.is_square(); // true
raw.sqrt().toString(); // '1'
R.__call__([1n, 2n, 1n]).is_square(true)[1]?.toString(); // 'x + 1'
const zero = K.zero();
zero.sqrt() === zero; // false
F.gen().sqrt(false, true); // [] over QQ; FpT instead raises ValueError
```


`Polynomial.squarefree_decomposition()` preserves each backend's unit and factor order.
Over QQ the first factor can be nonmonic; over word prime fields the nonconstant factors
are monic. A constant pair carries a nontrivial unit, including `[[zero, 1]]` for ZZ zero.
Generic zero polynomials raise ValueError, while word modular zero polynomials raise
ArithmeticError. Composite word moduli raise NotImplementedError for nonzero inputs.

FLINT exports `nmod_poly_factor_squarefree(f, p): Array<[bigint[], number]>` for prime
unsigned-word moduli. It omits the leading-coefficient unit, which its caller retains.

```typescript
import { nmod_poly_factor_squarefree } from '@sagemath-ts/flint-ts';
nmod_poly_factor_squarefree([0n, 1n, 0n, 0n, 0n, 0n, 1n], 5n);
// [[[0n, 1n], 1], [[1n, 1n], 5]]
const Qx = new PolynomialRing(QQ, 'x');
Qx.__call__([0n, 2n]).squarefree_decomposition().map(([f, e]) => [f.toString(), e]);
// [['2*x', 1]]
```

`GF2Element.is_square(): boolean` and `sqrt({extend?, all?})` support both binary
coefficients. Square roots return a coefficient, or a singleton array with `all:true`.
The parent caches its canonical entries; conversion of an existing GF2Element preserves
that input, and negation preserves a directly constructed zero.

```typescript
import { GF2, GF2Element } from 'sagemath-ts/rings/finite_rings';
const bit = GF2.__call__(1n);
GF2.__call__('9007199254740993') === GF2.one(); // true
GF2.__call__(null) === GF2.zero(); // true
GF2.__call__() === GF2.zero(); // true
bit.is_square(); // true
bit.sqrt() === bit; // true
bit.sqrt({all: true}); // [bit]
bit.add(1n) === GF2.zero(); // true
const rawZero = new GF2Element(0n, GF2);
rawZero.neg() === rawZero; // true
rawZero.sqrt() === GF2.zero(); // true
```


PARI exports `FpX_factor_squarefree(f: bigint[], p: bigint): bigint[][]`. Inputs use
little-endian reduced coefficients and a prime modulus. Entry `i` is the squarefree
component of multiplicity `i+1`, with `[1n]` placeholders. Components need not be
irreducible. Word-prime constants return `[]`; constants above the 64-bit word limit
return `[f]`. Zero raises (large-prime zero is rejected defensively).

```typescript
import { FpX_factor_squarefree } from '@sagemath-ts/parigp-ts';
FpX_factor_squarefree([2n, 3n, 1n], 101n); // [[2n, 3n, 1n]]
FpX_factor_squarefree([1n, 2n, 1n], 101n); // [[1n], [1n, 1n]]
FpX_factor_squarefree([2n, 2n], 101n); // [[2n, 2n]]
FpX_factor_squarefree([2n], 101n); // []
FpX_factor_squarefree([2n], (1n << 127n) - 1n); // [[2n]]
```


`Polynomial.denominator(): bigint | RingElement` returns the coefficient denominator. QQ uses
a positive bigint; generic coefficient rings return their native denominator or one.
`Polynomial.numerator()` clears that denominator. For QQ its result belongs to a
cached ZZ polynomial parent with Integer object coefficients; other backends retain
the polynomial parent. `Polynomial.lcm(other)` returns Sage's signed integer result
over ZZ and a monic result over fields.

```typescript
const T = new PolynomialRing(QQ, 't');
const f = T.__call__([QQ.__call__(1n).div(QQ.__call__(3n)), QQ.__call__(1n).div(QQ.__call__(2n))]);
f.denominator(); // 6n
f.numerator().toString(); // '3*t + 2'
f.numerator().parent.base_ring.toString(); // 'Integer Ring'
f.numerator().parent === f.numerator().parent; // true
f.lcm(T.__call__([1n, 1n])).toString(); // 't^2 + 5/3*t + 2/3'
```

FLINT exports `_fmpz_poly_lcm(a,b): bigint[]` (positive leading coefficient) and
`_fmpq_poly_lcm(a,b): [bigint[],bigint]` (monic rational result). Dense coefficient
arrays have no trailing zero entries. Rational input denominators cancel in the LCM.
`fmpq_poly_get_numerator(a,den)` copies canonical integer numerator storage, and
`fmpq_poly_get_denominator(a,den)` returns its positive denominator.

```typescript
import { _fmpz_poly_lcm, _fmpq_poly_lcm, fmpq_poly_get_numerator, fmpq_poly_get_denominator } from '@sagemath-ts/flint-ts';
_fmpz_poly_lcm([-2n, 2n], [3n, 3n]); // [-6n, 0n, 6n]
_fmpq_poly_lcm([1n, 2n], [1n, 2n]); // [[1n, 2n], 2n]
fmpq_poly_get_numerator([-1n, 0n, 3n], 2n); // [-1n, 0n, 3n]
fmpq_poly_get_denominator([-1n, 0n, 3n], 2n); // 2n
```

Place sets and valuation rings preserve native parent identity. Residue triples and
callable maps are cached by name. Public residue maps coerce into their domains and
reject poles with Sage's map-conversion TypeError. The direct `_residue_field()` maps
retain their lower-level protocol. Converting an existing function-field element
into its own field or a containing valuation ring retains that element.


`FunctionFieldDivisor._format(formatter, mul, cr)` takes a formatter for places and
bigint multiplicities, followed by multiplication and line-separator strings.
`_repr_(split = true)` uses native place formatting. Explicit zero multiplicities
remain in constructed divisors; addition and scalar multiplication discard zeros.
Divisor groups preserve identity, cache their zero, and accept existing divisors or
places. `function_space()` caches its `[dimension, from_V, to_V]` triple. Invalid
coordinates preserve Sage's pivot-list ValueError.

```typescript
import { DivisorGroup } from 'sagemath-ts/rings/function_field';
const DG = K.divisor_group();
const P = K.maximal_order().ideal(K.gen()).place();
new DivisorGroup(K) === DG; // true
DG.__call__(P).toString(); // 'Place (x)'
const D = P.divisor(0n);
D._format(v => `<${v}>`, ' @ ', ' / '); // '<0> @ <Place (x)>'
D.add(DG.zero()).support().length; // 0
DG.zero().function_space() === DG.zero().function_space(); // true
```

For supported QQ and prime-field vector domains, wrong-length coordinate arrays
raise the native map-conversion TypeError. Generic dense vectors (including prime
moduli at least 2147483647) expand `[]` to the zero vector. A zero-dimensional
function-space lifting map raises `TypeError('Cannot convert int to
sage.structure.element.Element')`, matching the original empty-sum boundary.


`FunctionFieldIdeal.parent()` returns its unique `IdealMonoid`; `IdealMonoid.one()`
returns the cached unit ideal. Ideal `pow(1)` retains the original ideal (including
a raw generator), and `pow(0)` returns that monoid identity. Finite and infinite
ideal `contains(value)` accept values coercible into the function field. Membership
in a zero ideal raises the native division-by-zero error, including for zero.

```typescript
const O = K.maximal_order(), I = O.ideal(K.gen());
I.parent() === O.ideal_monoid(); // true
I.pow(1n) === I; // true
I.pow(0n) === I.parent().one(); // true
O.ideal(K.one()).contains(1n); // true
```

`PolynomialRing.polynomials(options)` and `.monics(options)` return lazy iterators.
Pass exactly one of `of_degree?: number | null` and `max_degree?: number | null`.
Coefficients follow the base ring's iteration order, with the constant coefficient
varying fastest. Public factories reject infinite base rings with an empty
`NotImplementedError`; direct `_polys_degree(n)`, `_polys_max(n)`,
`_monics_degree(n)` and `_monics_max(n)` retain Sage's lazy internal protocol.
The internal `src/misc/mrange.ts` dependency `_xmrange_iter<T,U = T[]>(inputs:
Iterable<T>[], typ?: (row: T[]) => U): Generator<U>` restarts input iterators
without materializing them. Its converter is invoked without arguments on an empty
product; the default converter returns an empty array and otherwise copies rows.

```ts
import { PolynomialRing } from 'sagemath-ts/rings/polynomial';
import { GF } from 'sagemath-ts/rings/finite_rings';
const E = new PolynomialRing(GF(3n), 'x');
[...E.monics({ max_degree: 1 })].map(String) // ['1', 'x', 'x + 1', 'x + 2']
[...E.polynomials({ max_degree: -1 })].map(String) // ['0']
```

`get_place(1)` and `_places_finite(1)` now obtain their first place without allocating
all constant-field elements. Degree zero yields no finite places; `get_place(0)`
raises `AssertionError('there is a bug around')`. Negative `monomial(n)` behavior
follows the selected polynomial backend, except unsafe FLINT calls use the
`RangeError` documented in **Negative Polynomial Monomial Degrees** in DEVIATIONS.md.

`MPolynomial.toString()` uses balanced prime-field coefficients through characteristic
2147483647 (Singular's backend), and the coefficient's usual display above that
boundary. For example, over GF(5), `3*x` prints `-2*x`. Hyperelliptic constructors
reject low-degree models when their homogenization has a denominator, preserving
the original projective-coordinate-ring TypeError.

`ConstantField.cardinality?(): bigint | 'Infinity'` includes the inherited QQ method.
The exported `constant_field_cardinality(k): bigint` helper retains its finite-field
contract and rejects an infinite field with `TypeError`.

Function-field element `.matrix()` returns an immutable nested array. Matrix entries,
`.trace()` and `.norm()` run the native multiplication/summation steps and can
normalize a raw fraction without changing the original element. Writes raise Sage's
immutable-matrix `ValueError`; make a mutable copy with `M.map(row => row.slice())`.
Fresh calls return fresh matrices, as in Sage's Cython method lookup behavior.

```ts
import { QQ } from 'sagemath-ts/rings';
import { FunctionField, type ConstantField, type ConstantFieldElement } from 'sagemath-ts/rings/function_field';
const F = FunctionField(QQ as unknown as ConstantField<ConstantFieldElement>, 'x');
const M = F.gen().matrix();
Object.isFrozen(M) && Object.isFrozen(M[0]) // true
const C = M.map(row => row.slice());
C[0]![0] = F.one();
String(M[0]![0]) // 'x'; the original matrix is unchanged
```

`is_FunctionFieldElement(x)` and `is_FunctionField(x)` also check native-style
parent/category protocols. Missing element parents raise `AttributeError`; category
lookup catches `AttributeError`, while other category failures propagate.

```ts
import { FunctionField, is_FunctionFieldElement } from 'sagemath-ts/rings/function_field';
import { QQ } from 'sagemath-ts/rings';
const K = FunctionField(QQ as never, 'x');
is_FunctionFieldElement({ parent: () => K }) // true
is_FunctionFieldElement(1n) // false
```


`FunctionFieldElement_rational.matrix(base?: unknown)` accepts its own parent or
`null`/omitted base. On a fresh parent, any other hashable base raises
`ValueError('base must be the rational function field itself')`. Sage's cached
free-module result ignores that base after a successful matrix, trace or norm call
on any element of the same parent. Lists, dictionaries and sets still raise the
native unhashable-type `TypeError`, even after warming. Matrices themselves remain
fresh immutable arrays. This specialization does not expose `free_module()`.


`change_variable_name(name)` returns `[field, from_field, to_field]`; both maps accept
`unknown` and convert inputs into their source field. Identity maps also convert
scalars and retain existing source-field elements. Nested singleton lists follow
the fraction-field constructor.

`PolynomialRing.__call__` accepts `Map` exponent dictionaries and numeric-key Records.
For example, over QQ, `{ 2: 3n, 0: -2n }` constructs `3*x^2 - 2`. Backend-specific
coefficient and exponent rules are preserved; tuple keys are accepted by ZZ and
generic dense backends, while QQ and word-prime templates reject them.

### FLINT deterministic word random state

`@sagemath-ts/flint-ts` exports `flint_rand_init()`, `flint_rand_set_seed(state, a, b)`,
`flint_rand_get_seed(state)` and the no-op `flint_rand_clear(state)`. The state type
is `flint_rand_t`; seeds and word results are BigInts. `n_randlimb`, `n_randint`,
`n_urandint`, `n_randbits`, `n_randtest_bits`, `n_randtest` and `n_randtest_not_zero`
use that mutable state. Bit-count arguments are numbers from 0 through 64.

```ts
import { flint_rand_init, n_randlimb } from '@sagemath-ts/flint-ts';
n_randlimb(flint_rand_init()) // 15737102703946861599n
```

### Native FLINT polynomial factor backend

From `@sagemath-ts/flint-ts` (coefficients in ascending degree):

```ts
const [unit, factors] = nmod_poly_factor([4n, 0n, 0n, 0n, 1n], 5n);
// unit === 1n; factors: [[4n,1n], [2n,1n], [3n,1n], [1n,1n]], each exponent 1.
```

`P = readonly bigint[]`, `M = readonly P[]`, `Factors = Array<[bigint[], number]>`.
All `p` arguments are word primes. Mutable random-state inputs use `flint_rand_t`.

| Export | Signature / result |
|---|---|
| `nmod_poly_factor`, `nmod_poly_factor_with_cantor_zassenhaus`, `nmod_poly_factor_with_kaltofen_shoup` | `(f: P, p: bigint): [bigint, Factors]`; native insertion order |
| `nmod_poly_factor_cantor_zassenhaus`, `nmod_poly_factor_kaltofen_shoup` | `(f: P, p: bigint): Factors`; positive-degree input |
| `nmod_poly_factor_distinct_deg` | `(f: P, p: bigint): Factors`; squarefree input; pairs contain common irreducible degree |
| `nmod_poly_factor_equal_deg` | `(f: P, d: number, p: bigint): Factors`; squarefree, all irreducible factors degree `d` |
| `nmod_poly_factor_equal_deg_prob` | `(state, f: P, d: number, p: bigint): bigint[] \| null`; degree at least two |
| `nmod_poly_deflation` | `(f: P, p: bigint): number` |
| `nmod_poly_deflate`, `nmod_poly_inflate` | `(f: P, stride: number, p: bigint): bigint[]`; positive deflation / nonnegative inflation |
| `nmod_poly_remove` | `(f: P, divisor: P, p: bigint): [bigint[], number]`; remaining polynomial and multiplicity |
| `_nmod_vec_rand`, `_nmod_vec_randtest`, `nmod_poly_randtest` | `(state, len: number, p: bigint): bigint[]`; polynomial variant removes trailing zeros |
| `nmod_mat_mul` | `(a: M, b: M, p: bigint, ncols?: number): bigint[][]`; `ncols` retains empty shapes |
| `nmod_poly_precompute_matrix` | `(inner: P, f: P, finv: P, p: bigint): bigint[][]` |
| `_nmod_poly_reduce_matrix_mod_poly` | `(rows: M, f: P, p: bigint): bigint[][]` |
| `nmod_poly_compose_mod_brent_kung_precomp_preinv` | `(outer: P, rows: M, f: P, finv: P, p: bigint): bigint[]` |
| `nmod_poly_compose_mod_brent_kung_vec_preinv` | `(polys: M, count: number, inner: P, f: P, finv: P, p: bigint): bigint[][]` |
| `_nmod_poly_mod_matrix_rows_evaluate` | `(rows: M, h: P, f: P, finv: P, p: bigint): bigint[]` |

`finv` is the power-series inverse of the reversed modulus through its length,
computed with `_nmod_poly_inv_series_newton`. Composition requires outer degree
smaller than modulus degree; the matrix-reduction helper is used when reducing a
precomputed matrix to a smaller positive-degree modulus. See `DEVIATIONS.md` for
native preconditions and safe error boundaries. Sage's odd word-prime polynomial
`factor()` delegates to this backend and sorts the returned factors.

### PARI native random state

From `@sagemath-ts/parigp-ts`:

```ts
setrand(1n);
const saved = getrand();
pari_rand(); // 13282407956253574712n
setrand(saved); // restores the same next output
```

| Export | Signature / result |
|---|---|
| `pari_init_rand` | `(): void`; reset to native seed 1 |
| `setrand`, `getrand` | `(seed: bigint): void`, `(): bigint`; positive word seed or saved native state |
| `pari_rand` | `(): bigint`; unsigned 64-bit output |
| `random_bits` | `(bits: number): bigint`; widths 1–64, signed native result at width 64 |
| `random_Fl`, `randomi` | `(limit: bigint): bigint`; positive word / arbitrary positive integer limit, exclusive |
| `random_F2x` | `(length: number): bigint`; packed binary coefficients |
| `random_zv` | `(length: number): bigint[]`; signed 64-bit words |
| `random_Flx`, `random_FpX` | `(length: number, p: bigint): bigint[]`; ascending coefficients modulo positive word / arbitrary positive modulus |

All calls share the module's native state. Array lengths are nonnegative safe
integers; polynomial outputs omit trailing zeros. Empty samples consume no state.
See `DEVIATIONS.md` for guards outside native C's defined input contract.

### PARI binary polynomial and matrix kernels

From `@sagemath-ts/parigp-ts`, a nonnegative bigint encodes coefficient `i` in bit `i`.

```ts
F2x_mul(7n, 3n); // 9n: (x^2+x+1)(x+1) = x^3+1
F2x_factor(63n); // [[3n, 1], [7n, 2]]
const columns = [3n, 3n];
F2m_ker_sp(columns, 2); // [3n]; columns becomes [3n, 1n]
```

| Export | Signature / result |
|---|---|
| `F2x_degree` | `(f: bigint): number`; zero has degree -1 |
| `F2x_add`, `F2x_mul`, `F2x_gcd`, `F2x_rem` | `(a: bigint, b: bigint): bigint` |
| `F2x_divrem` | `(a: bigint, b: bigint): [bigint, bigint]`; quotient/remainder |
| `F2x_sqr`, `F2x_sqrt`, `F2x_deriv`, `F2x_Frobenius` | `(f: bigint): bigint`; sqrt requires a square |
| `F2x_valrem` | `(f: bigint): [bigint, bigint]`; valuation/remainder, zero valuation `(1n << 63n)-1n` |
| `F2xq_mul`, `F2xq_sqr` | `(a, b, modulus: bigint): bigint`, `(a, modulus: bigint): bigint` |
| `F2xq_powers` | `(a: bigint, lastPower: number, modulus: bigint): bigint[]`; includes power zero and original power one |
| `F2x_matFrobenius` | `(modulus: bigint): bigint[]`; packed columns, positive-degree modulus |
| `F2m_ker_sp` | `(columns: bigint[], rows: number, deplin = 0): bigint[] \| bigint \| null`; mutates columns |
| `F2m_ker` | `(columns: readonly bigint[], rows: number): bigint[]`; preserves input |
| `F2x_factor` | `(f: bigint): Array<[bigint, number]>`; sorted factor/multiplicity pairs |
| `F2x_factor_squarefree` | `(f: bigint): bigint[]`; component index is multiplicity minus one, placeholders are `1n` |
| `F2x_ddf` | `(squarefree: bigint): Array<[bigint, number]>`; factor/common-irreducible-degree pairs |

For nonzero `deplin`, the kernel call returns the first dependency vector or null;
empty input returns the same empty array. Direct `F2x_factor(0n)` returns `[[0n,1]]`.
Sage's `Polynomial.factor()` rejects zero, delegates binary fields to PARI, then
sorts by degree, multiplicity and polynomial comparison. Native preconditions and
safe guards are documented in `DEVIATIONS.md`.


PARI quotient helper exports accept ascending bigint coefficient arrays:
`FpXQ_powers(x, n: number, T, p: bigint)` returns `[1, x, x², ..., xⁿ]`,
with the initial x copied raw for counts at most two or a generic-size modulus.
For counts greater than two and 0 < abs(p) < 2^64, native word conversion reduces
that initial x first. `FpXQ_autpow(x, n, T, p)` performs binary
composition. `FpXQ_autpowers(x, n, T, p)` has an unused array at slot zero,
then unreduced identity and x, followed by successive compositions. Counts must
be nonnegative safe integers.

```typescript
import { FpXQ_powers, FpXQ_autpowers } from '@sagemath-ts/parigp-ts';
FpXQ_powers([1n, 1n], 2, [1n], 3n); // [[1n], [1n, 1n], []]
FpXQ_autpowers([1n, 1n], 1, [1n], 3n); // [[], [0n, 1n], [1n, 1n]]
```


PARI matrix products use columns with unused slot zero in both dimensions:
`ZM_mul(A, B)`, `Flm_mul(A, B, p: bigint)` and `FpM_mul(A, B, p: bigint)`.
Inputs must be rectangular and dimension-compatible. Flm coefficients are reduced
modulo a positive 64-bit p; FpM accepts integer coefficients and any positive p.
`F2m_mul(A: readonly bigint[], B: readonly bigint[], rows: number): bigint[]`
uses packed binary columns without dummy slots, with A's row count supplied.

```typescript
import { ZM_mul, FpM_mul, F2m_mul } from '@sagemath-ts/parigp-ts';
const A = [[], [0n, 1n, 3n], [0n, 2n, 4n]];
ZM_mul(A, A); // [[], [0n, 7n, 15n], [0n, 10n, 22n]]
FpM_mul(A, A, 3n); // [[], [0n, 1n, 0n], [0n, 1n, 1n]]
F2m_mul([3n, 2n], [3n, 1n], 2); // [1n, 3n]
```


PARI composition exports include `brent_kung_optpow(d: number, n: number, m: number)`,
`Flxq_powers(x, n: number, T, p: bigint)`, `Flx_Flxq_eval(Q, x, T, p)`,
`Flx_FlxqV_eval(Q, powers, T, p)` and `FpX_FpXQV_eval(Q, powers, T, p)`.
Polynomials are ascending bigint arrays; power tables are zero-indexed. Word
inputs use reduced coefficients. Precomputed tables must be large enough, and
word powers copied into a matrix must have degree smaller than the modulus.
`FpXQ_auttrace([phi, a], n: bigint, T, p)` powers the automorphism/additive-trace
pair for a nonzero unsigned 64-bit exponent. Both automorphism vectors canonicalize
trailing zero storage without reducing copied coefficients. Trace and successive
automorphism tables prepare native reciprocals even for counts one and zero;
`FpXQ_autpow` retains its count-zero/count-one remainder shortcuts.

```typescript
import { brent_kung_optpow, Flxq_powers, FpXQ_auttrace } from '@sagemath-ts/parigp-ts';
brent_kung_optpow(10, 2, 1); // 5
Flxq_powers([1n, 1n], 3, [1n, 0n, 1n], 5n); // [[1n], [1n, 1n], [0n, 2n], [3n, 2n]]
FpXQ_auttrace([[0n, 1n], [1n, 1n]], 3n, [1n, 0n, 1n], 5n); // [[0n, 1n], [3n, 3n]]
```


PARI also exports `ZX_sqr(a: readonly bigint[]): bigint[]`,
`Flx_mul(a: readonly bigint[], b: readonly bigint[], p: bigint): bigint[]`,
`Flx_sqr(a: readonly bigint[], p: bigint): bigint[]` and
`FpX_sqr(a: bigint[], p: bigint): bigint[]`. Direct Flx inputs must have reduced
coefficients and a positive 64-bit modulus. FpX inputs may have signed coefficients.

```typescript
import { ZX_sqr, Flx_mul, Flx_sqr, FpX_sqr } from '@sagemath-ts/parigp-ts';
ZX_sqr([1n, -1n, 1n]); // [1n, -2n, 3n, -2n, 1n]
Flx_mul([1n, 2n], [3n, 4n], 5n); // [3n, 0n, 3n]
Flx_sqr([1n, 2n], 5n); // [1n, 4n, 4n]
FpX_sqr([-1n, 1n], 3n); // [1n, 1n, 1n]
```

PARI division exports `FpX_divrem(a, T, p)` / `Flx_divrem(a, T, p)` returning
`[bigint[], bigint[]]`, and `FpX_rem(a, T, p)` / `Flx_rem(a, T, p)` returning
`bigint[]`. Here `a` and `T` are `bigint[]`, and `p` is `bigint`.
`FpX_invBarrett(T: bigint[], p: bigint): bigint[]` and
`Flx_invBarrett(T: bigint[], p: bigint): bigint[]` return the reciprocal-series
inverse through degree(T)-2. Direct Flx inputs follow the reduced-word contract.

```typescript
import { FpX_divrem, FpX_rem, Flx_divrem, Flx_rem, FpX_invBarrett, Flx_invBarrett } from '@sagemath-ts/parigp-ts';
FpX_divrem([1n, 2n, 3n], [1n, 1n], 5n); // [[4n, 3n], [2n]]
FpX_rem([1n, 2n, 3n], [1n, 1n], 5n); // [2n]
Flx_divrem([1n, 2n, 3n], [1n, 1n], 5n); // [[4n, 3n], [2n]]
Flx_rem([1n, 2n], [2n], 8n); // []
FpX_invBarrett([1n, 2n, 3n, 1n], 5n); // [1n, 2n]
Flx_invBarrett([1n, 2n, 3n, 1n], 5n); // [1n, 2n]
```

PARI `FpX_gcd(a, b, p)` and `Flx_gcd(a, b, p)` return the **unscaled** native gcd.
`FpX_extgcd(a, b, p)` / `Flx_extgcd(a, b, p)` return `[g, u, v]` with
`u*a + v*b = g`; g is not automatically monic. All polynomial arguments/results
are ascending `bigint[]`, and p is `bigint`.

`FpX_halfgcd(a, b, p)` / `Flx_halfgcd(a, b, p)` return a zero-indexed 2x2 **row**
matrix of polynomials. `FpX_halfgcd_all` / `Flx_halfgcd_all` additionally return the
transformed pair as `[M, a1, b1]`, where `[a1,b1]^t = M*[a,b]^t`.

```typescript
import { FpX_gcd, FpX_extgcd, Flx_gcd, Flx_extgcd, FpX_halfgcd, FpX_halfgcd_all, Flx_halfgcd, Flx_halfgcd_all, gen_pow_i } from '@sagemath-ts/parigp-ts';
const a = [1n, 2n, 3n], b = [1n, 1n];
FpX_gcd(a, b, 5n); // [2n]
Flx_gcd(a, b, 5n); // [2n]
FpX_extgcd(a, b, 5n); // [[2n], [1n], [1n, 2n]]
Flx_extgcd(a, b, 5n); // [[2n], [1n], [1n, 2n]]
FpX_halfgcd(a, b, 5n); // [[[], [1n]], [[1n], [1n, 2n]]]
Flx_halfgcd(a, b, 5n); // [[[], [1n]], [[1n], [1n, 2n]]]
FpX_halfgcd_all(a, b, 5n); // [[[[], [1n]], [[1n], [1n, 2n]]], [1n, 1n], [2n]]
Flx_halfgcd_all(a, b, 5n); // [[[[], [1n]], [[1n], [1n, 2n]]], [1n, 1n], [2n]]
gen_pow_i(3n, 3n, x => x*x, (x,y) => x*y); // 27n
```

`gen_pow_i<T>(x: T, n: bigint, square: (x: T) => T, multiply: (x: T, y: T) => T): T`
uses the magnitude of a nonzero exponent; callers invert x themselves for negative
powers. Quotient-power APIs perform that inversion automatically.

The exported `PolynomialMatrix` type is `[[bigint[], bigint[]], [bigint[], bigint[]]]`;
`HalfGcdResult` is `[PolynomialMatrix, bigint[], bigint[]]`.


PARI quotient minimal polynomials use Shoup's randomized algorithm. Both functions
consume the package's global PARI random state; use `setrand` for reproducibility.
`FpXQ_minpoly` accepts integer coefficient arrays and delegates word primes;
`Flxq_minpoly` requires reduced coefficients and a positive word prime. The
polynomial modulus must have positive degree. Word input degree must be smaller
than the polynomial modulus degree.

```ts
import { FpXQ_minpoly, Flxq_minpoly, setrand } from '@sagemath-ts/parigp-ts';
setrand(1n);
FpXQ_minpoly([0n, 1n], [1n, 0n, 1n], 7n); // [1n, 0n, 1n]
setrand(1n);
Flxq_minpoly([0n, 1n], [1n, 0n, 1n], 7n); // [1n, 0n, 1n]
```

PARI distinct-degree factorization groups squarefree factors by their common
irreducible degree. `FpX_ddf` returns the existing `Map<number, bigint[]>`;
`Flx_ddf` returns `[bigint[], number][]` (`PolynomialDegreeFactor[]`). Word
factors retain PARI's scale and can be nonmonic. `Flx_nbfact_by_degree` returns
`{D, nb}`, with unused `D[0] = 0`; `FpX_nbfact` returns the total factor count.
These kernels require a prime modulus; direct Flx input coefficients are reduced.

```ts
import { FpX_ddf, Flx_ddf, FpX_nbfact, Flx_nbfact_by_degree } from '@sagemath-ts/parigp-ts';
const f = [5n, 2n, 6n, 0n, 0n, 2n, 6n, 3n, 1n];
FpX_ddf(f, 7n).get(1); // [2n, 4n, 5n]
Flx_ddf(f, 7n)[0]; // [[2n, 4n, 5n], 1]
FpX_nbfact(f, 7n); // 4
Flx_nbfact_by_degree(f, 7n); // {D: [0, 2, 1, 0, 1, 0, 0, 0, 0], nb: 4}
```

PARI `Fp_sqrt(a: bigint, p: bigint): bigint | null` returns the native smallest
square root or `null` for a nonresidue. The modulus is assumed prime.
`Fp_sqrt_i(a: bigint, y: bigint | null, p: bigint)` has the same return type;
`y` is a generator of the multiplicative group's 2-Sylow subgroup, or `null` to
search deterministically. The word route ignores `y`.
`Fl_sqrt(a: bigint, p: bigint): bigint | null` requires `0 <= a < p < 2^64`.
`gen_pow_fold<T>(x: T, n: bigint, square: (x: T) => T, multiplySquare: (x: T) => T): T`
uses nonzero `|n|` and a fused callback computing `x * value^2`.

```ts
import { Fp_sqrt, Fp_sqrt_i, Fl_sqrt, gen_pow_fold } from '@sagemath-ts/parigp-ts';
Fp_sqrt(2n, 17n); // 6n
Fp_sqrt(3n, 17n); // null
Fp_sqrt_i(2n, null, 17n); // 6n
Fl_sqrt(2n, 17n); // 6n
gen_pow_fold(3n, -5n, x => x*x, x => 3n*x*x); // 243n
```


PARI full polynomial factors use ascending bigint coefficients:
`FpX_factor(f: bigint[], p: bigint): PolynomialFactor[]` and
`Flx_factor(f: bigint[], p: bigint): PolynomialFactor[]`, where
`PolynomialFactor = [bigint[], number]`. They return monic irreducible factors with
multiplicities, discard the leading unit and consume the global PARI random state.
Both require p >= 2; direct Flx requires p < 2^64 and reduced coefficients.
`FpX_normalize(f: bigint[], p: bigint): bigint[]` and
`Flx_normalize(f: bigint[], p: bigint): bigint[]` make the polynomial monic.
Generic normalization preserves zero; direct Flx normalization/factorization of
zero raises the native inverse error for p > 1.

```ts
import { FpX_factor, Flx_factor, FpX_normalize, Flx_normalize } from '@sagemath-ts/parigp-ts';
FpX_factor([1n, 0n, 0n, 0n, 1n], 5n); // [[[2n, 0n, 1n], 1], [[3n, 0n, 1n], 1]]
Flx_factor([2n, 4n, 2n], 5n); // [[[1n, 1n], 2]]
FpX_factor([], 7n); // [[[], 1]]
FpX_normalize([2n, 4n, 2n], 5n); // [1n, 2n, 1n]
Flx_normalize([2n, 4n, 2n], 5n); // [1n, 2n, 1n]
```

Sage `Polynomial.factor()` also delegates large-prime coefficient fields to PARI,
then restores the leading unit and Sage factor order. Extension-field factorization
remains a separate implementation under audit.


PARI `FpX_roots(f: bigint[], p: bigint): bigint[]` and
`Flx_roots(f: bigint[], p: bigint): bigint[]` return distinct roots. Generic inputs
are reduced before checking their degree; a zero polynomial raises `PariError`.
Both root entries require p >= 2; direct Flx entries require a positive 64-bit
modulus and reduced coefficients.

`FpX_nbroots(f, p)` / `Flx_nbroots(f, p)` return a number (zero polynomial: -1).
`FpX_is_totally_split(f, p)` / `Flx_is_totally_split(f, p)` return a boolean and
retain native raw-degree shortcuts. Total splitting means distinct linear factors:
a repeated linear factor fails the predicate. The generic zero predicate is false
for primes below 2^64 and raises an inverse error for larger primes.

Root order follows PARI exactly: the general word-field path sorts signed 64-bit
words, so roots >= 2^63 precede smaller roots. Its degree-at-most-two shortcut
sorts unsigned values; the larger-prime path sorts ordinary nonnegative integers.

```ts
import { FpX_roots, Flx_roots, Flx_nbroots, FpX_is_totally_split, Flx_is_totally_split } from '@sagemath-ts/parigp-ts';
FpX_roots([1n, 7n], 7n); // []
Flx_roots([6n, 0n, 1n], 7n); // [1n, 6n]
Flx_nbroots([1n, 2n, 1n], 7n); // 1
Flx_is_totally_split([1n, 2n, 1n], 7n); // false
FpX_is_totally_split([7n], 7n); // true: native raw constant-degree shortcut
const p = (1n << 64n) - 59n;
FpX_roots([2n, p - 1n, p - 2n, 1n], p); // [p - 1n, 1n, 2n]
```


PARI extension-polynomial arithmetic preserves native coefficient storage:
`ExtensionCoefficient = bigint | bigint[]` and
`ExtensionPolynomial = ExtensionCoefficient[]`. A bigint is a native integer
coefficient; an array is a polynomial coefficient in the inner variable `y`.
All polynomial arrays use ascending powers. Direct word coefficients are always
arrays (`bigint[][]`); binary coefficients are packed bit polynomials (`bigint[]`).

| Operation | Generic / word signatures | Binary signature |
|---|---|---|
| Product | `FpXQX_mul(x,y,T,p)` / `FlxqX_mul(x,y,T,p)` | `F2xqX_mul(x,y,T)` |
| Square | `FpXQX_sqr(x,T,p)` / `FlxqX_sqr(x,T,p)` | `F2xqX_sqr(x,T)` |
| Reduce coefficients | `FpXQX_red(x,T,p)` / `FlxqX_red(x,T,p)` | `F2xqX_red(x,T)` |
| Make monic | `FpXQX_normalize(x,T,p)` / `FlxqX_normalize(x,T,p)` | `F2xqX_normalize(x,T)` |

Each returns its input representation. Generic and word `T` are `bigint[]`,
`p` is a bigint >= 2, and binary `T` is packed bigint. The extension modulus
must have positive degree. Direct word inputs require p < 2^64 and coefficients
in [0,p). Generic/word Kronecker products and squares require inner coefficient degree below
deg(T); mixed generic scalar products retain PARI's transpose shortcut. Binary
products preserve native XOR packing even for overlapping raw coefficient blocks;
each shifted coefficient must fit the native 64-bit word allocation. Binary squares, reduction and normalization allow unreduced inner inputs.
These are native field primitives: no primality/irreducibility validation is added.

Native generic squaring of all-integer coefficients returns an unreduced integer
polynomial. Generic normalization with leading one preserves lower coefficients;
word and binary normalization still invert and reduce. Use `_red` explicitly when
canonical coefficient reduction is required.

```ts
import { FpXQX_sqr, FpXQX_red, FpXQX_normalize, FlxqX_mul, FlxqX_normalize,
  F2xqX_sqr, F2xqX_normalize } from '@sagemath-ts/parigp-ts';
const T = [2n, 0n, 1n];
FpXQX_sqr([4n, 4n], T, 5n); // [16n, 32n, 16n]
FpXQX_red([16n, 32n, 16n], T, 5n); // [1n, 2n, 1n]
FpXQX_normalize([[1n, 0n, 1n], [1n]], T, 5n); // [[1n, 0n, 1n], 1n]
FlxqX_normalize([[1n, 0n, 1n], [1n]], T, 5n); // [[4n], [1n]]
FlxqX_mul([[1n, 1n], [1n]], [[1n, 1n], [1n]], T, 5n); // [[4n, 2n], [2n, 2n], [1n]]
F2xqX_sqr([3n, 1n], 7n); // [2n, 0n, 1n]
F2xqX_normalize([7n, 1n], 7n); // [0n, 1n]
```


PARI extension-polynomial division accepts the same coefficient representations:
generic P is ExtensionPolynomial, word P is bigint[][], binary P is bigint[].
An outer divisor may be P or a cached reduction object:
`ExtensionReduction<P> = { polynomial: P; inverse: P }`;
`ExtensionModulus<P> = P | ExtensionReduction<P>`. Obtain valid caches with
the matching `_get_red` function; it returns a plain polynomial below PARI's
threshold and retains an existing cache.

| Operation | Generic / word | Binary | Result |
|---|---|---|---|
| Quotient and remainder | `FpXQX_divrem(x,S,T,p)` / `FlxqX_divrem(x,S,T,p)` | `F2xqX_divrem(x,S,T)` | [P,P] |
| Remainder | `FpXQX_rem(x,S,T,p)` / `FlxqX_rem(x,S,T,p)` | `F2xqX_rem(x,S,T)` | P |
| Quotient | `FpXQX_div(x,S,T,p)` / `FlxqX_div(x,S,T,p)` | `F2xqX_div(x,S,T)` | P |
| Truncated reciprocal | `FpXQX_invBarrett(S,T,p)` / `FlxqX_invBarrett(S,T,p)` | `F2xqX_invBarrett(S,T)` | P |
| Prepare reduction | `FpXQX_get_red(S,T,p)` / `FlxqX_get_red(S,T,p)` | `F2xqX_get_red(S,T)` | ExtensionModulus<P> |

`_invBarrett` takes a plain polynomial. The reciprocal is truncated to
max(0, degree(S)-1) coefficients, with trailing zeros removed. The positive-degree inner modulus and coefficient guards
from extension arithmetic also apply here. Raw inner coefficients retain native
behavior, including the degree constraints in Kronecker/Newton packing paths.
Generic division delegates primes below 2^64 to word arithmetic, converting constant
polynomial coefficients to integers. Its early `_rem` return preserves their tags.

This bundled PARI version has a binary Barrett remainder branch that retains stale
high coefficients. The port preserves that behavior: for example, `F2xqX_rem` of
X^100+1 by X² over the field with T=7 returns the original input.
`F2xqX_divrem` retains its separate native behavior and gives remainder [1n] there.

```ts
import { FpXQX_divrem, FpXQX_rem, FlxqX_divrem, FlxqX_get_red,
  F2xqX_divrem } from '@sagemath-ts/parigp-ts';
const T = [3n, 0n, 1n], x = [[1n], [], [1n]], S = [[1n], [1n]];
FpXQX_divrem(x, S, T, 17n); // [[16n, 1n], [2n]]
FlxqX_divrem(x, S, T, 17n); // [[[16n], [1n]], [[2n]]]
FpXQX_rem([[1n]], S, T, 17n); // [[1n]]: native early-return tags
F2xqX_divrem([1n, 0n, 1n], [1n, 1n], 7n); // [[1n, 1n], []]
const divisor = Array.from({length: 16}, () => [1n]);
const cached = FlxqX_get_red(divisor, T, 17n);
FlxqX_divrem(divisor, cached, T, 17n); // [[[1n]], []]
```


PARI extension-polynomial GCD uses the same P, T and p representations:

| Operation | Generic / word | Binary | Result |
|---|---|---|---|
| Native unscaled GCD | FpXQX_gcd(x,y,T,p) / FlxqX_gcd(x,y,T,p) | F2xqX_gcd(x,y,T) | P |
| Bézout coefficients | FpXQX_extgcd(x,y,T,p) / FlxqX_extgcd(x,y,T,p) | F2xqX_extgcd(x,y,T) | [P,P,P] = [gcd,U,V] |
| Half-GCD matrix | FpXQX_halfgcd(x,y,T,p) / FlxqX_halfgcd(x,y,T,p) | F2xqX_halfgcd(x,y,T) | ExtensionMatrix<P> |

ExtensionMatrix<P> is [[P,P],[P,P]] in row-major order. GCD is not made monic.
The retained native binary remainder quirk can also make GCD fail to terminate:
F2xqX_gcd(X^100+1, X², 7n) cycles. Extgcd uses paired division and avoids that case.
Extgcd follows native output timing: a zero first polynomial raises a division
error when computing U. Half-GCD follows the native initial swap/remainder matrix
conventions and does not first reduce the inputs.

```ts
import { FpXQX_gcd, FpXQX_extgcd, FlxqX_extgcd, F2xqX_halfgcd,
  F2xqX_mul } from '@sagemath-ts/parigp-ts';
const T = [3n, 0n, 1n];
FpXQX_gcd([2n, 2n], [4n, 4n], T, 17n); // [4n, 4n]
FpXQX_extgcd([1n, 0n, 1n], [1n, 1n], T, 17n); // [[2n], [1n], [1n, 16n]]
FlxqX_extgcd([[1n], [], [1n]], [[1n], [1n]], T, 17n); // [[[2n]], [[1n]], [[1n], [16n]]]
F2xqX_halfgcd([1n], [8n, 1n], 7n); // [[[1n], []], [[1n, 1n], [1n]]]
F2xqX_mul([32n, 1n], [1n], 7n); // []: native XOR packing cancellation
```

PARI extension quotient operations take the same coefficient representation P and
an outer modulus S of type ExtensionModulus<P>:

| Operation | Generic / word | Binary | Result |
|---|---|---|---|
| Multiply | FpXQXQ_mul(x,y,S,T,p) / FlxqXQ_mul(x,y,S,T,p) | F2xqXQ_mul(x,y,S,T) | P |
| Square | FpXQXQ_sqr(x,S,T,p) / FlxqXQ_sqr(x,S,T,p) | F2xqXQ_sqr(x,S,T) | P |
| Safe inverse | FpXQXQ_invsafe(x,S,T,p) / FlxqXQ_invsafe(x,S,T,p) | F2xqXQ_invsafe(x,S,T) | P or null |
| Inverse | FpXQXQ_inv(x,S,T,p) / FlxqXQ_inv(x,S,T,p) | F2xqXQ_inv(x,S,T) | P |
| Divide | FpXQXQ_div(x,y,S,T,p) / FlxqXQ_div(x,y,S,T,p) | — | P |
| Signed power | FpXQXQ_pow(x,n,S,T,p) / FlxqXQ_pow(x,n,S,T,p) | F2xqXQ_pow(x,n,S,T) | P |
| Power table | FpXQXQ_powers(x,l,S,T,p) / FlxqXQ_powers(x,l,S,T,p) | F2xqXQ_powers(x,l,S,T) | P[] |

Signed exponents n are bigint. FlxqXQ_powu(x,n,S,T,p) additionally accepts bigint
exponents from 0 through 2^64-1 and preserves native unsigned-power shortcuts.
Power-table l is a number: a nonnegative integer less than 2^32-1. Tables include
the zeroth power and retain the original first power. Shared extension guards
apply; safe inverse may still propagate native arithmetic errors over composites.
Primitive operations accept an outer reciprocal cache; inner reciprocals are
prepared and reused internally at native powering cutoffs.

Generic powering at word primes retains a bundled conversion quirk: a cached
scalar-coefficient S is interpreted as one inner coefficient. Cached S containing
polynomial coefficients raises RangeError on that conversion to exclude invalid
native GEN reads. This conversion occurs for exponents other than 0 and ±1.
Use a plain S for generic powering, or the word entry point with word coefficients
when supplying a cache.

Inverse error payloads follow native GEN formatting and its 1,600-character cutoff.
Long payloads use the native omitted-object message, including its newline.

```ts
import { FpXQXQ_inv, FpXQXQ_invsafe, FpXQXQ_powers, FlxqXQ_powu,
  F2xqXQ_inv, F2xqXQ_powers } from '@sagemath-ts/parigp-ts';
const T = [3n, 0n, 1n], S = [1n, 0n, 1n];
FpXQXQ_inv([1n, 1n], S, T, 17n); // [9n, 8n]
FpXQXQ_invsafe(S, S, T, 17n); // null
FpXQXQ_powers([1n, 1n], 3, S, T, 17n); // [[1n], [1n, 1n], [0n, 2n], [15n, 2n]]
FlxqXQ_powu([[1n], [1n]], 3n, [[1n], [], [1n]], T, 17n); // [[15n], [2n]]
F2xqXQ_inv([2n, 1n], [3n, 1n, 1n], 7n); // [2n, 3n]
F2xqXQ_powers([2n, 1n], 3, [3n, 1n, 1n], 7n); // [[1n], [2n, 1n], [0n, 1n], [3n, 3n]]
```


### PARI extension composition

Import these names from `@sagemath-ts/parigp-ts`. `P` denotes an extension
polynomial: `(bigint | bigint[])[]` for generic, `bigint[][]` for word, and packed
`bigint[]` for binary coefficients. Outer modulus `S` accepts the same plain or
cached shape as the quotient operations above. Inner `T` is `bigint[]`, or packed
`bigint` for binary. Generic/word calls end in `p: bigint`.

| Function | Arguments before S, T, p | Result |
|---|---|---|
| FpXQX_FpXQXQ_eval / FlxqX_FlxqXQ_eval | Q: P, x: P | P |
| F2xqX_F2xqXQ_eval | Q: P, x: P (no p argument) | P |
| FpXQX_FpXQXQV_eval / FlxqX_FlxqXQV_eval | Q: P, powers: P[] | P |
| F2xqX_F2xqXQV_eval | Q: P, powers: P[] (no p argument) | P |

A supplied table may contain raw entries. Generic/binary blocks reduce inner
coefficients; word conversion truncates outer coefficients to deg(S) before its
matrix product. In the word APIs a nonzero Q requires a sufficient table and a
nonzero S when matrix conversion is reached. Existing coefficient/modulus guards
apply; see DEVIATIONS.md for source-order error handling.

`FlxqM_mul(A: bigint[][][], B: bigint[][][], T: bigint[], p: bigint)` returns the
same matrix type. Column and row slot zero are unused arrays; each real cell is
an ascending word polynomial. Matrices must be rectangular and compatible, with
reduced word coefficients and 2 <= p < 2^64.

```typescript
import { FpXQX_FpXQXQ_eval, FlxqX_FlxqXQ_eval,
  F2xqX_F2xqXQ_eval, FlxqM_mul } from '@sagemath-ts/parigp-ts';
FpXQX_FpXQXQ_eval([1n,1n,1n], [1n,1n], [1n,0n,1n], [3n,0n,1n], 17n); // [2n,3n]
FlxqX_FlxqXQ_eval([[1n],[1n],[1n]], [[1n],[1n]], [[1n],[],[1n]], [3n,0n,1n], 17n); // [[2n],[3n]]
F2xqX_F2xqXQ_eval([1n,2n,1n], [2n,1n], [3n,1n,1n], 7n); // [2n,3n]
const A = [[], [[],[1n],[2n]], [[],[3n],[4n]]];
FlxqM_mul(A, A, [3n,0n,1n], 17n); // [[],[[],[7n],[10n]],[[],[15n],[5n]]]
```


### PARI power-table cache boundaries

`FpXQ_powers` and `Flxq_powers` remove numeric trailing zeros while preserving the
first power's unreduced coefficient values. They prepare native reciprocal caches
before returning even a count-zero table. `FpX_FpXQV_eval` and `Flx_FlxqV_eval`
also prepare the reciprocal for a constant nonzero Q; canonical zero Q returns
zero before preparation. A nonunit leading coefficient can therefore raise a
native inverse error once the relevant modulus-degree threshold is reached.

```typescript
import { FpXQ_powers, FpX_FpXQV_eval } from '@sagemath-ts/parigp-ts';
FpXQ_powers([18n,0n], 1, [3n,0n,1n], 17n); // [[1n],[18n]]
const T = [1n, ...Array(35).fill(0n), 2n];
FpX_FpXQV_eval([0n,0n], [], T, 4n); // []
// FpXQ_powers([1n], 0, T, 4n) throws PariError:
// impossible inverse in Fp_inv: Mod(2, 4).
```


Automorphism storage and native inverse errors:

```typescript
import { FpXQ_autpowers, FpXQ_auttrace } from '@sagemath-ts/parigp-ts';
FpXQ_autpowers([18n, 0n], 1, [3n, 0n, 1n], 17n); // [[], [0n, 1n], [18n]]
const T = [1n, ...Array(35).fill(0n), 2n];
FpXQ_auttrace([[1n], [1n]], 1n, T, 4n); // throws PariError: impossible inverse in Fp_inv: Mod(2, 4).
```


Coefficient substitution in a polynomial over an extension field is exported as:

- `FpXY_FpXQ_evalx(P, x: bigint[], T: bigint[], p: bigint)` and
  `FpXY_FpXQV_evalx(P, V: bigint[][], T, p)` return `ExtensionPolynomial`.
  `P` is an `ExtensionPolynomial`: integer coefficients are copied unchanged,
  while polynomial coefficients are evaluated modulo `(T, p)`.
- `FlxY_Flxq_evalx(P: bigint[][], x: bigint[], T: bigint[], p: bigint)` and
  `FlxY_FlxqV_evalx(P, V: bigint[][], T, p)` return `bigint[][]` and use reduced
  word coefficients.
- `F2xY_F2xq_evalx(P: bigint[], x: bigint, T: bigint)` and
  `F2xY_F2xqV_evalx(P, V: bigint[], T)` return `bigint[]`, whose entries are
  packed binary coefficient polynomials.
- `F2x_F2xq_eval(Q: bigint, x: bigint, T: bigint)` and
  `F2x_F2xqV_eval(Q, V: bigint[], T)` return one packed binary polynomial.

`V` is a zero-indexed power table. Direct substitution constructs a table even
for an empty outer polynomial; a supplied-table generic call with only integer
coefficients copies them without using `T` or `p` for arithmetic. Existing polynomial/table
input contracts and native error ordering apply.

```typescript
import { FpXY_FpXQ_evalx, FlxY_Flxq_evalx, F2xY_F2xq_evalx, F2x_F2xq_eval } from '@sagemath-ts/parigp-ts';
FpXY_FpXQ_evalx([18n, [18n, 0n]], [1n], [3n, 0n, 1n], 17n); // [18n, [1n]]
FlxY_Flxq_evalx([[1n], [1n, 1n]], [2n], [3n, 0n, 1n], 17n); // [[1n], [3n]]
F2xY_F2xq_evalx([3n, 5n], 2n, 7n); // [3n, 2n]
F2x_F2xq_eval((1n << 4096n) - 1n, 2n, 7n); // 1n
```


### PARI extension automorphisms

All eight APIs are exported from `@sagemath-ts/parigp-ts`. Generic `P` is
`ExtensionPolynomial`, word `P` is `bigint[][]`, and binary `P` is `bigint[]`
with packed binary coefficients. `S: ExtensionModulus<P>` is the outer modulus
or its existing cache object. Generic/word `T` is `bigint[]`; binary `T` is a
packed `bigint`. Generic/word calls end in `p: bigint`; binary calls omit p.

| Function | Arguments before S, T, p | Result |
|---|---|---|
| FpXQXQ_autpow / FlxqXQ_autpow | aut: [bigint[], P], n: bigint | [bigint[], P] |
| FpXQXQ_auttrace / FlxqXQ_auttrace | aut: [P, P], n: bigint | [P, P] |
| FpXQXQ_autsum / FlxqXQ_autsum | aut: [bigint[], P, P], n: bigint | [bigint[], P, P] |
| F2xqXQ_autpow | aut: [bigint, P], n: bigint | [bigint, P] |
| F2xqXQ_auttrace | aut: [bigint, P, P], n: bigint | [bigint, P, P] |

Power tuples contain the inner generator image and outer generator image.
Generic/word trace tuples contain the outer image and additive aggregate, with
no inner image. Their autsum triples also contain the inner image and multiply
the aggregate. Binary trace triples include the inner image and add the aggregate.

The seven signed native entries accept nonzero signed 64-bit n. Negative n is
cast to its unsigned 64-bit bit pattern, as in C; it does not request an inverse.
`FlxqXQ_auttrace` accepts 1 <= n < 2^64. Zero and out-of-range n raise RangeError.
Existing extension-field and table-storage guards apply. Even n=1 prepares inner
and outer reciprocal caches before returning canonical storage with raw values.

```typescript
import { FpXQXQ_auttrace, FpXQXQ_autsum, F2xqXQ_auttrace } from '@sagemath-ts/parigp-ts';
const T = [1n, 0n, 1n], S = [1n, 0n, 1n], B = [0n, 1n];
FpXQXQ_auttrace([B, [2n, 1n]], 3n, S, T, 17n); // [[0n,1n],[6n,3n]]
FpXQXQ_autsum([[0n,1n], B, [2n,1n]], 3n, S, T, 17n); // [[0n,1n],[0n,1n],[2n,11n]]
FpXQXQ_auttrace([B, [1n]], -1n, S, T, 17n); // [[0n,1n],[]]
F2xqXQ_auttrace([2n, [0n,1n], [1n,1n]], 3n, [1n,0n,1n], 7n); // [2n,[0n,1n],[1n,1n]]
```


### PARI extension projections and truncated products

Import these eight APIs from `@sagemath-ts/parigp-ts`. Generic `P` below is
`ExtensionPolynomial`; word `P` is `bigint[][]`. Coefficients are ascending.

| Function | Arguments | Result |
|---|---|---|
| random_FpXQX / random_FlxqX | length: number, T: bigint[], p: bigint | P |
| FpXQX_dotproduct | x: P, y: P, T: bigint[], p: bigint | bigint or bigint[] |
| FlxqX_dotproduct | x: P, y: P, T: bigint[], p: bigint | bigint[] |
| FpXQXn_mul / FlxqXn_mul | x: P, y: P, n: number, T: bigint[], p: bigint | P |
| FpXQXn_sqr / FlxqXn_sqr | x: P, n: number, T: bigint[], p: bigint | P |

Random sampling consumes the PARI random state in ascending outer/inner order.
Only the canonical degree of T is used. Constant T gives zero coefficients;
length zero consumes no random values, even with zero T or p. Word T entries
must fit unsigned 64-bit storage, but need not be reduced modulo p for sampling.
A requested coefficient of positive degree requires a positive sampling bound.

Dot products use the shorter canonical outer length and reduce after summing.
An empty generic dot product is scalar `0n`; a polynomial-valued cancellation
retains the polynomial zero `[]`. Truncation keeps outer powers below n. Generic
packed products truncate before coefficient reduction; word products reduce
first, which can raise an inverse error even at n=0. Generic all-integer squares
retain exact integer coefficients, as the native routine does. Existing extension
arithmetic input guards apply to dot products and truncated arithmetic. Counts
must be nonnegative safe integers; native variable metadata is omitted.

```typescript
import { setrand, random_FpXQX, random_FlxqX, FpXQX_dotproduct,
  FlxqX_dotproduct, FpXQXn_mul, FpXQXn_sqr, FlxqXn_mul, FlxqXn_sqr } from '@sagemath-ts/parigp-ts';
const T = [1n, 0n, 1n];
setrand(1n);
random_FpXQX(2, T, 17n); // [[13n,11n],[11n,3n]]
setrand(1n);
random_FlxqX(2, T, 17n); // [[13n,11n],[11n,3n]]
FpXQX_dotproduct([1n,1n], [1n,-1n], T, 17n); // 0n
FlxqX_dotproduct([[1n],[2n]], [[3n],[4n]], T, 17n); // [11n]
FpXQXn_sqr([18n,19n], 2, T, 17n); // [324n,684n]
FpXQXn_mul([18n,19n], [18n,19n], 2, T, 17n); // [1n,4n]
FlxqXn_mul([[1n],[2n]], [[1n],[2n]], 2, T, 17n); // [[1n],[4n]]
FlxqXn_sqr([[1n],[2n]], 2, T, 17n); // [[1n],[4n]]
```


### PARI extension minimal polynomials

`FpXQXQ_minpoly(x: ExtensionPolynomial, S: ExtensionModulus<ExtensionPolynomial>,
T: bigint[], p: bigint): ExtensionPolynomial` and
`FlxqXQ_minpoly(x: bigint[][], S: ExtensionModulus<bigint[][]>, T: bigint[],
p: bigint): bigint[][]` are exported from `@sagemath-ts/parigp-ts`.
They compute the monic minimal polynomial over the coefficient field F_p[t]/T
of x in its outer quotient by S. S may have repeated factors.

Supply a coefficient field (prime p and irreducible T), positive-degree S, and
an outer representative of degree below deg(S). The adapter rejects zero or
constant S and oversized outer representatives with RangeError; existing field
and word-coefficient guards apply. Primality and irreducibility are assumed,
not checked. Generic results preserve native scalar/polynomial coefficient tags.

The native Shoup algorithm consumes PARI's random state. It prepares the outer
reciprocal with the original T before constructing the power table's inner
reciprocal, so nonunit inputs can expose different errors at those stages.

```typescript
import { setrand, FpXQXQ_minpoly, FlxqXQ_minpoly } from '@sagemath-ts/parigp-ts';
setrand(1n);
FpXQXQ_minpoly([0n,1n], [1n,0n,1n], [3n,0n,1n], 17n); // [[1n],[],1n]
setrand(1n);
FlxqXQ_minpoly([[],[1n]], [[1n],[],[1n]], [3n,0n,1n], 17n); // [[1n],[],[1n]]
setrand(1n);
FpXQXQ_minpoly([[1n,1n]], [1n,0n,1n], [3n,0n,1n], 17n); // [[16n,16n],1n]
setrand(1n);
FpXQXQ_minpoly([0n,1n], [1n,0n,0n,0n,1n], [1n,1n,1n], 2n); // [[1n],[],[],[],1n]
```


Prime quotient powers and inverses preserve complete native exceptions. Powers
with exponent zero or +/-1 keep their earlier shortcuts; other exponents prepare
one reciprocal after any negative-power inversion. A nonunit leading coefficient
at a cache threshold can therefore raise an error even for a constant base.
`FpXQ_inv(a,T,p)` uses the generic inverse, while the explicit four-argument
Hensel form uses its native unsigned-word initialization. An undefined fourth
argument selects the three-argument behavior.

```typescript
import { FpXQ_pow, FpXQ_inv } from '@sagemath-ts/parigp-ts';
const T = [1n, ...Array(90).fill(0n), 2n];
FpXQ_pow([1n], 0n, T, 4n); // [1n]
FpXQ_pow([1n], 2n, T, 4n); // throws PariError: impossible inverse in Fl_inv: Mod(2, 4).
FpXQ_inv([], [1n,0n,1n], 17n); // throws PariError: impossible inverse in FpXQ_inv: 0.
FpXQ_pow([], -2n, [1n,0n,1n], 17n); // throws PariError: impossible inverse in Flxq_inv: 0.
```


### PARI extension Frobenius

These exports are available from `@sagemath-ts/parigp-ts`. Outer coefficients are
ascending; generic `ExtensionPolynomial` coefficients retain integer/polynomial
storage tags. `ExtensionModulus<P>` accepts P or `{polynomial: P, inverse: P}`.

| Export | Signature |
|---|---|
| `Flx_Frobenius` | `(T: bigint[], p: bigint): bigint[]` |
| `FpXQX_Frobenius` | `(S: ExtensionModulus<ExtensionPolynomial>, T: bigint[], p: bigint): ExtensionPolynomial` |
| `FlxqX_Frobenius` | `(S: ExtensionModulus<bigint[][]>, T: bigint[], p: bigint): bigint[][]` |
| `F2xqX_Frobenius` | `(S: ExtensionModulus<bigint[]>, T: bigint): bigint[]` |
| `FpXQXQ_halfFrobenius` | `(a: ExtensionPolynomial, S: ExtensionModulus<ExtensionPolynomial>, T: bigint[], p: bigint): ExtensionPolynomial` |
| `FlxqXQ_halfFrobenius` | `(a: bigint[][], S: ExtensionModulus<bigint[][]>, T: bigint[], p: bigint): bigint[][]` |

`Flx_Frobenius` computes t^p modulo T. The outer Frobenius functions compute
X^(p^deg(T)) modulo S over the coefficient field. For odd prime p, half Frobenius
computes a^((p^deg(T)-1)/2). Characteristic two retains PARI's aggregate exponent
2^deg(T)-1. The native field assumptions are not checked by primality or
irreducibility tests. The existing extension/word guards apply; generic inputs
that become a constant coefficient modulus on word conversion raise RangeError.
Generic full Frobenius also inherits the documented cached-polynomial conversion
guard from FpXQXQ_pow. Half Frobenius converts both parts of the cache correctly.

```typescript
import {
  Flx_Frobenius, FpXQX_Frobenius, FlxqX_Frobenius, F2xqX_Frobenius,
  FpXQXQ_halfFrobenius, FlxqXQ_halfFrobenius,
} from '@sagemath-ts/parigp-ts';
const T = [3n,0n,1n], S = [1n,0n,1n], W = [[1n],[],[1n]];
Flx_Frobenius(T, 17n);                            // [0n,16n]
FpXQX_Frobenius(S, T, 17n);                       // [0n,1n]
FlxqX_Frobenius(W, T, 17n);                       // [[],[1n]]
F2xqX_Frobenius([1n,1n,1n], 7n);                  // [0n,1n]
FpXQXQ_halfFrobenius([1n,1n], S, T, 17n);          // [1n]
FlxqXQ_halfFrobenius([[1n],[1n]], W, T, 17n);       // [[1n]]
```


### PARI prime-polynomial linear boundaries

The existing `FpX_red`, `FpX_add`, `FpX_sub`, `FpX_neg` and `FpX_Fp_mul` exports
preserve native low-level signed-modulus reduction and zero-modulus shortcuts.
This does not make a negative or zero modulus a finite-field characteristic.
If an operation reaches division by zero, it raises
`PariError('impossible inverse in dvmdii: 0.')`.

```typescript
import { FpX_red, FpX_add, FpX_Fp_mul } from '@sagemath-ts/parigp-ts';
FpX_red([-1n,18n,0n], -17n);          // [16n,1n]
FpX_add([0n,1n], [0n,-1n], 0n);       // []: cancellation returns before division
FpX_Fp_mul([1n,1n,1n], 0n, 0n);       // []: raw scalar zero returns first
FpX_Fp_mul([0n,0n], 1n, 0n);          // []: canonical empty polynomial
```


### PARI extension root counts, split parts and derivatives

The following exports are available from `@sagemath-ts/parigp-ts`. Coefficients
are ascending. Generic `ExtensionPolynomial` coefficients are bigint or bigint[];
word coefficients are bigint[]; binary coefficients and T are packed bigint bits.

| Export | Signature |
|---|---|
| `FpXQX_split_part` | `(f: ExtensionPolynomial, T: bigint[], p: bigint): ExtensionPolynomial` |
| `FpXQX_nbroots` | `(f: ExtensionPolynomial, T: bigint[], p: bigint): number` |
| `FlxqX_nbroots` | `(f: bigint[][], T: bigint[], p: bigint): number` |
| `F2xqX_nbroots` | `(f: bigint[], T: bigint): number` |
| `FqX_nbroots` | `(f: ExtensionPolynomial, T: bigint[] \| null, p: bigint): number` |
| `FpXQX_is_squarefree` | `(f: ExtensionPolynomial, T: bigint[], p: bigint): boolean` |
| `FlxqX_is_squarefree` | `(f: bigint[][], T: bigint[], p: bigint): boolean` |
| `FpXX_deriv` | `(f: ExtensionPolynomial, p: bigint): ExtensionPolynomial` |
| `FlxX_deriv` | `(f: bigint[][], p: bigint): bigint[][]` |

Counts are distinct roots, with -1 for the native zero-polynomial sentinel.
`FpXQX_split_part` retains PARI's unscaled GCD result; it need not be monic.
`FqX_nbroots(f, null, p)` counts over the prime field and requires integer
coefficients; a non-null T selects the extension. Native degree-zero/one shortcuts
precede use of T. Word conversion can change the raw degree before that shortcut.
The remaining extension/GCD guards and native field assumptions still apply.

Derivatives act on the outer variable and reduce coefficients modulo p without
reducing by an inner modulus. Generic scalar/polynomial tags are preserved.
At p=0, zero integer coefficients return before reduction; a reached scalar
operation raises dvmdii, while a polynomial-coefficient operation raises umodui.

```typescript
import {
  FpXQX_split_part, FpXQX_nbroots, FlxqX_nbroots, F2xqX_nbroots, FqX_nbroots,
  FpXQX_is_squarefree, FlxqX_is_squarefree, FpXX_deriv, FlxX_deriv,
} from '@sagemath-ts/parigp-ts';
const T = [1n,0n,1n], S = [1n,0n,1n], W = [[1n],[],[1n]];
FpXQX_split_part(S, T, 3n);                       // [1n,0n,1n]
FpXQX_nbroots(S, T, 3n);                          // 2
FlxqX_nbroots(W, T, 3n);                          // 2
F2xqX_nbroots([1n,1n,1n], 7n);                    // 2
FqX_nbroots(S, null, 3n);                         // 0
FqX_nbroots(S, T, 3n);                            // 2
FpXQX_is_squarefree([1n,-2n,1n], T, 3n);           // false
FlxqX_is_squarefree(W, T, 3n);                    // true
FpXX_deriv([1n,[2n,3n],4n], 17n);                 // [[2n,3n],8n]
FlxX_deriv([[1n],[2n,3n],[4n]], 17n);              // [[2n,3n],[8n]]
```


PARI low-level matrix products and composition preserve signed modulus magnitudes.
These boundary behaviors do not make nonpositive moduli finite fields.

```typescript
import { FpXQ_powers, FpX_FpXQ_eval, FpM_mul } from '@sagemath-ts/parigp-ts';
FpXQ_powers([-1n], 2, [1n,1n], -2n); // [[1n],[-1n],[1n]]
FpXQ_powers([-1n], 3, [1n,1n], -2n); // [[1n],[1n],[1n],[1n]]
FpX_FpXQ_eval([1n,2n,3n], [1n], [1n,1n], -17n); // [6n]
FpM_mul([[]], [[]], 0n); // [[]]
FpM_mul([[],[0n,1n],[0n,1n]], [[],[0n,1n,-1n]], 0n); // [[],[0n,0n]]
```

At modulus zero, a nonzero entry in the integer matrix product raises
`PariError('impossible inverse in dvmdii: 0.')` during reduction.


PARI polynomial observation helpers accept ascending bigint coefficients:
`FpX_deriv(f,p)`, `FpX_eval(f,x,p)`, `FpX_center(f,p,half)` and
`FpX_div_by_X_x(f,x,p)` are exported from `@sagemath-ts/parigp-ts`.
The last returns only the quotient by X-x. Centering assumes reduced coefficients;
it compares magnitudes with half and conditionally subtracts p, without a modular
reduction or a final storage normalization.

```typescript
import { FpX_deriv, FpX_eval, FpX_center, FpX_div_by_X_x } from '@sagemath-ts/parigp-ts';
FpX_deriv([1n,0n,0n], 0n); // []
FpX_eval([-17n,1n], 17n, 0n); // 0n (cancellation before reduction)
FpX_eval([1n,2n,3n], -1n, -17n); // 2n
FpX_center([1n,16n], 17n, 8n); // [1n,-1n]
FpX_div_by_X_x([1n,1n], 17n, 0n); // [1n] (no remainder requested)
```

These nonpositive-modulus examples describe native low-level boundary behavior,
not finite fields. Reached division by zero raises the native PariError.


PARI scalar exports preserve native signed operands and low-level shortcuts:

```typescript
import { Fp_inv, Fp_mul, Fp_red, Fp_div } from '@sagemath-ts/parigp-ts';
Fp_inv(-1n, 17n); // 16n
Fp_mul(-1n, 1n, 17n); // 16n
Fp_red(0n, 0n); // 0n
Fp_div(0n, 17n, 17n); // 0n (native word-divisor shortcut)
Fp_div(1n, 2n, (1n<<64n)+13n); // 9223372036854775815n
```

`Fp_sqr(0n,0n)` raises `PariError('impossible inverse in dvmdii: 0.')`;
its native remainder primitive differs from Fp_red's zero shortcut. Failed
inversion reports PariError with PARI's exact payload. These boundary cases
should not be interpreted as general division rules for finite fields.


PARI's low-level Fp_pow preserves its zero-exponent divisibility shortcut:

```typescript
import { Fp_pow } from '@sagemath-ts/parigp-ts';
Fp_pow(17n, 0n, 17n); // 0n
Fp_pow(0n, 0n, 0n); // 0n
Fp_pow(1n, 0n, 0n); // 1n
Fp_pow(-1n, -1n, 17n); // 16n
```

For noninvertible bases, a negative exponent can report Fl_inv or Fp_inv in the
PariError payload, depending on the exponent and modulus word sizes.


PARI scalar order functions preserve native word-bound defaults:

```typescript
import { Fp_order, znorder, xgcd } from '@sagemath-ts/parigp-ts';
Fp_order(2n, 0n, 17n); // 8n: zero defaults to p-1 for a word modulus
Fp_order(2n, 1n << 64n, 17n); // 8n: a larger-than-word bound also defaults
znorder(2n, -17n, 0n); // 8n
xgcd(0n, 0n); // [0n, 0n, 0n]
```

Generic Fp_order moduli require a positive order bound. Failed scalar znorder
coprimality checks raise PariError with both normalized operands in its payload.


PARI rational polynomial helpers canonicalize the common denominator and reduce
fractions before finite-field conversion:

```typescript
import { QPoly_normalize, QPoly_to_fractions, QPoly_to_FpX } from '@sagemath-ts/parigp-ts';
QPoly_normalize({num:[2n,4n,0n],den:-6n}); // {num:[-1n,-2n],den:3n}
QPoly_to_fractions({num:[2n,4n,0n],den:-6n}); // [[-1n,3n],[-2n,3n]]
QPoly_to_FpX({num:[2n,4n,6n],den:2n},2n); // [1n,0n,1n]
```

A zero QPoly denominator raises PariError, including for the zero polynomial.


PARI Galois integer helpers retain native word arithmetic and checked logarithms:

```typescript
import { ulcm, factoru_small, logint, Forprime } from '@sagemath-ts/parigp-ts';
ulcm(3, 2**63); // 2**63 (native unsigned 64-bit overflow)
factoru_small(0); // [[0,1]]
logint(3n**8192n, 3n); // 8192
new Forprime(499979).next(); // 500009 (native initialized-table transition)
```

Forprime keeps exact internal progress, returns number values, and returns zero
when unsigned-word primes are exhausted. Integer-number arguments use native
magnitude conversion and reject magnitudes at least 2^64. Logarithm bases below
two and nonpositive operands raise PariError; the base is validated first.


PARI integer-polynomial observations normalize trailing zeros. Nonzero constants
are squarefree, and indexpartial includes native p-adic denominator refinement:

```typescript
import { ZX_neg, ZX_is_squarefree, indexpartial } from '@sagemath-ts/parigp-ts';
ZX_neg([1n, 0n]); // [-1n]
ZX_is_squarefree([2n]); // true
indexpartial([-2n, 0n, 0n, 0n, 1n]); // 8n
indexpartial([1n, -2n, 1n]); // 0n (zero discriminant)
```


PARI Hensel lifting retains one-indexed root/factor vectors and native monic
normalization. A single already monic factor retains its original coefficients:

```typescript
import { ZpX_liftfact, bezout_lift_fact, ZpX_ZpXQ_liftroot } from '@sagemath-ts/parigp-ts';
ZpX_liftfact([4n, 1n], [[], [0n, 1n]], 2n, 1); // [[], [4n, 1n]]
const f = [-9n, -2n, -1n], factors = [[], [0n, 1n], [2n, 1n]];
ZpX_liftfact(f, factors, 3n, 3); // [[], [18n, 1n], [11n, 1n]]
bezout_lift_fact(f, factors, 3n, 3); // [[], [10n, 23n], [18n, 4n]]
ZpX_ZpXQ_liftroot([-1n, 0n, 1n], [4n], [1n, 0n, 1n], 3n, 1); // [4n]
```

Root/factor inputs must satisfy the lifting preconditions; an inexact polynomial
division raises RangeError. For multiple factors, bezout_lift_fact requires
precision at least two. Scalar inverse failures use native PariError messages.


PARI inverse Vandermonde matrices retain one-indexed rows and columns:

```typescript
import { FpV_invVandermonde } from '@sagemath-ts/parigp-ts';
FpV_invVandermonde([0n, 0n, 1n], 1n, 5n);
// [[0n,0n,0n], [0n,1n,0n], [0n,4n,1n]]
```

A singular point set raises PariError even when the denominator multiplier is zero;
native batch inversion occurs before multiplication by that denominator.


PARI permutation helpers retain one-indexed vectors and checked word exponents:

```typescript
import { perm_powu, perm_cycles, cyc_pow, vecsmall_lexcmp } from '@sagemath-ts/parigp-ts';
const p = [0, 2, 3, 1];
perm_powu(p, -1); // [0,2,3,1] (native unsigned magnitude conversion)
cyc_pow(perm_cycles(p), -1); // [[],[0,1,3,2]] (signed exponent)
vecsmall_lexcmp([0], [0,1,2,3]); // -1
```

Unsigned exponents reject magnitudes at least 2^64; signed cycle exponents reject
magnitudes at least 2^63. Permutation orders round to Number after unsigned-word
arithmetic. group_order rounds once after exact signed products and raises
RangeError if an intermediate product leaves PARI's signed-word domain.


PARI symmetric-polynomial evaluation returns `bigint[] | bigint`. Vectors are
one-indexed. Each Newton sum is reduced modulo `p`; the weighted combination is
an integer column, with no final reduction. An empty or all-zero coefficient
list returns scalar `0n`. Cancellation of nonzero coefficients returns a column.

```typescript
import { sympol_eval } from '@sagemath-ts/parigp-ts';
const orbits = [[], [0n, 2n], [0n, 3n]];
sympol_eval({ v: [0, 4], w: [0, 1] }, orbits, 5n); // [0n,8n,12n]
sympol_eval({ v: [0, -1], w: [0, 1] }, orbits, 5n); // [0n,-2n,-3n]
sympol_eval({ v: [0, 0], w: [0, 1] }, orbits, 5n); // 0n
```

Callers that accept zero weights must narrow this result before using array
methods. Signed exponent entries follow PARI's unsigned-word cast. Internally,
`arith1.ts`'s `Fp_powu(a, exponent, modulus)` accepts bigint arguments, wraps the
exponent to 64 bits and delegates to the shared native powering schedule. It is
not added to the package root exports.


Unit subgroups preserve PARI's enumeration order, including ties in subgroup size:

```typescript
import { listznstarelts } from '@sagemath-ts/parigp-ts';
listznstarelts(8, 2); // [[0,1],[0,1,3],[0,1,7],[0,1,5]]
listznstarelts(9, -1); // [[0,1],[0,1,4,7]] (unsigned-word order argument)
```

The internal `char.ts` znstar, `subgroup.ts` subgrouplist and `hnf_snf.ts`
ZM_hnfmodid use zero-indexed bigint vectors/column matrices. These dependencies
are available by source subpath and are not added to package root exports.

```typescript
import { subgrouplist } from '@sagemath-ts/parigp-ts/src/subgroup.ts';
subgrouplist([4n]); // [[[4n]],[[2n]],[[1n]]] (all subgroups, native order)
subgrouplist([4n], 2n); // [[[2n]],[[1n]]] (index at most 2)
subgrouplist([4n], [2n]); // [[[2n]]] (index exactly 2)
```

`ifactor1.ts` provides `absZ_factor_limit_strict(n, limit)`, returning
`[factors, unresolvedOrNull]`, where each factor and unresolved value is a
`[bigint base, bigint exponent]` pair. Limit zero uses the configured default
500000. The unresolved base may itself be prime when the native branch has not
certified it. `ifactor.ts`'s default-limit helper delegates to this engine.


`galoisconj4(T, den?)` returns `QPoly[]`, with PARI's flag-four identity fallback
when its Galois search returns no structure. Generic linear input also returns
`x`; the `x - 1` and `x + 1` cyclotomic shortcuts return reduced constants.
`galoisinit` retains its nullable result. Validation raises `PariError`, checks
squarefreeness before monicity, and includes constant-polynomial values in errors.

```typescript
import { galoisconj4 } from '@sagemath-ts/parigp-ts';
galoisconj4([-2n, 1n]); // [{ num: [0n,1n], den: 1n }]
galoisconj4([-2n, 0n, 0n, 1n]); // [{ num: [0n,1n], den: 1n }]
galoisconj4([-1n, 1n]); // [{ num: [1n], den: 1n }]
```

The search routines require irreducible defining polynomials; squarefreeness
alone does not establish that precondition. Some reducible inputs return no
structure, while others cause the native and ported searches to continue.


The internal source-subpath helper `fixedfieldsurmer(l, NS, W)` searches for a
separating symmetric polynomial. `NS` is a one-indexed list of one-indexed
Newton-sum columns and `W` is a one-indexed list of 1–31 weights, matching its
native caller's range. It returns `{v,w}` or `null` after exhausting its native
base-four search. It uses an exact search counter, including beyond 32-bit limits.

```typescript
import { fixedfieldsurmer } from '@sagemath-ts/parigp-ts/src/galconj.ts';
fixedfieldsurmer(17n, [[], [0n,0n,1n]], [0,1]); // {v:[0,1],w:[0,1]}
```


Additional Galois helpers are available from `@sagemath-ts/parigp-ts/src/galconj.ts`:

- `radicalu(n: number)` and `eulerphiu(n: number)` return numbers after native
  unsigned-word magnitude conversion. The totient returns 2 at zero. The exact
  `eulerphiu(n: bigint): bigint` backend in `src/arith2.ts` takes `0 <= n < 2^64`.
- `findpsi(D, pstart, P, S, o, out)` selects a Frobenius prime and writes
  `{Tmod, psi, p}` to `out`. `D` is a nonzero bad-prime multiple; `P` defines a
  Galois field and `S: QPoly` generates a normal cyclic subgroup of order `o > 1`.
  `pstart` and `o` are numbers. `Tmod` and `psi` are one-indexed output vectors.
- `inittest(L, M, borne, ladic)` constructs the permutation-filter state from
  one-indexed roots and a one-indexed row-major matrix. `galois_test_perm(td, pf)`
  returns a boolean and updates the state's row order and cached test matrices.
- `galoisfindgroups(lo, sg, f)` filters subgroup lists by their sorted distinct
  residues modulo `f`. The outer `lo` list is zero-indexed; its rows and `sg` are
  one-indexed number vectors.

```typescript
import { eulerphiu, radicalu, findpsi, inittest, galois_test_perm, galoisfindgroups }
  from '@sagemath-ts/parigp-ts/src/galconj.ts';
eulerphiu(0); // 2
eulerphiu(-12); // 4
radicalu(12); // 6
const out = {Tmod: [] as bigint[][], psi: [] as number[], p: 0};
findpsi(-4n, 2, [1n,0n,1n], {num:[0n,-1n],den:1n}, 2, out); // 3
out.psi; // [0,1]
const td = inittest([0n,1n,2n], [[],[0n,1n,0n],[0n,0n,1n]], 10n, 1n << 96n);
galois_test_perm(td, [0,2,1]); // true
galoisfindgroups([[0,1,3],[0,1,2]], [0,1], 2); // [[0,1,3]]
```


The Galois source subpath also exposes `FpM_ker(M, nrows, ncols, p)` for zero-indexed
row matrices; it returns a zero-indexed list of kernel columns. `lcmBig(a,b)`
returns the nonnegative bigint LCM. The `alglin1.ts` backend's `FpM_ker(M,p)` and
`Flv.ts`'s word-prime `Flm_ker(M,p)` use one-indexed column matrices with an unused
entry in each column. Prime moduli are positive. Native generic-field kernels
can retain signed, unreduced coefficients once both dimensions reach five.

`F3v.ts: F3m_ker(columns, rows)` takes packed bigint columns (two bits per trit,
least-significant trit first) and returns packed kernel columns. The supplied row
count preserves zero padding. These helpers are available from their source
subpaths; package-root exports are unchanged.

```typescript
import { FpM_ker, lcmBig } from '@sagemath-ts/parigp-ts/src/galconj.ts';
import { FpM_ker as kernel } from '@sagemath-ts/parigp-ts/src/alglin1.ts';
import { Flm_ker } from '@sagemath-ts/parigp-ts/src/Flv.ts';
import { F3m_ker } from '@sagemath-ts/parigp-ts/src/F3v.ts';
FpM_ker([[1n,1n],[0n,0n]], 2, 2, 17n); // [[16n,1n]]
lcmBig(-6n,8n); // 24n
kernel([[],[0n,1n],[0n,1n]], 17n); // [[],[0n,16n,1n]]
Flm_ker([[],[0n,1n],[0n,1n]], 17n); // [[],[0n,16n,1n]]
F3m_ker([1n,1n], 1); // [6n], representing the kernel vector [2,1]
```


`src/ifactor.ts` exposes `squfof(n)`, `pollardbrent(n)` and
`Z_pollardbrent(n, rounds, seed)`, returning a factor list or `null`. The first two
now retain native stage-selection gates: SQUFOF declines at 2^46; Pollard–Brent
declines one-limb inputs and two-limb inputs whose low limb is below 2^32. Later
factorization stages handle declined inputs. Explicit round budgets and seeds
are integer numbers, computed exactly internally. Retry exhaustion raises
`PariError`; ordinary budget exhaustion returns `null`.

```typescript
import { Z_pollardbrent } from '@sagemath-ts/parigp-ts/src/ifactor.ts';
Z_pollardbrent(1022117n, 2**26, 0); // [1009n,1013n]
```

The internal `src/galconj.ts` helper
`galoisgenfixedfield0(O, L, sigma, T, bad, gf, gb)` is also available for direct
comparisons. `O` and `L` are one-indexed orbit/root vectors, `sigma` is a `QPoly`,
`bad` is a bigint or null, and `gf`/`gb` are the native-style Frobenius/bound states.
It updates the prime/factor/psi fields of `gf` and returns `{PG,Pg,V}` or null;
`V` contains the symmetric polynomial, its root values and the defining polynomial.
The prime is replaced when the fixed polynomial is not squarefree modulo it.


The PARI source subpath `src/ifactor.ts` also exports power-stage helpers:
`is_357_power(x: bigint, mask: number): [number, bigint, number]` takes a positive
integer and a mask from 0 through 7 (cube/fifth/seventh bits). It returns the
exponent, root, and updated mask; failure returns `[0, 0n, 0]`. Word-sized inputs
prioritize smaller exponents; larger inputs prioritize larger exponents, as in
native PARI. For example, `is_357_power(32n, 7)` is `[5, 2n, 2]`, and
`is_357_power(3n ** 21n, 7)` is `[3, 2187n, 5]`.
`is_pth_power(x: bigint, iterator: { next(): number | null }, cutoffbits: number)`
returns `[number, bigint]`. Supply an increasing prime iterator starting at 11
or later, a positive input and positive cutoff. Failure returns `[0, x]` and
consumes the first prime that fails the cutoff (or exhausts the iterator).
The iterator can wrap `forprime(start: number, end: number): Generator<number>`.


The internal PARI source module `src/_ecm.ts` exposes
`ECM_loop(N: bigint, nbc: number, seed: number, B1: number): bigint | null` for a
single native ECM round. Inputs follow the native stage preconditions: positive
odd N, 4 through 64 curves in multiples of four, positive integer seed and a
supported factor-stage B1 bound. For example,
`ECM_loop(4295229443n, 8, 1, 142)` returns `65537n`.
`new ECM(N, nbc, seed).round(B1)` retains state and advances the BigInt seed across
rounds. Normal consumers should use `Z_factor`; `ifactor.ts:ellfacteur` retains
its existing `maxRounds` limit in insisting mode.


`ifactor.ts:is_kth_power(x: bigint, n: number): bigint | null` takes a nonzero
integer and a positive integer exponent. Negative inputs may return null at a
modular filter or raise `PariError` with
`sorry, sqrtnr for x < 0 is not yet implemented.` if they pass those filters.
For example, `is_kth_power(-16n, 3)` is null, while `is_kth_power(-8n, 3)` raises
that error. `is_kth_power(2n, 17886697)` is null without a large allocation.

The test-only `src/mpqs.ts:mpqsInternals` object additionally exposes `newHandle`,
`mpqs_FB_ctor` and `combine_large_primes` for native relation comparisons. Its
`Fl_sqrt(a: number, p: number)` accepts reduced quadratic residues for word
primes in the factor-base range and selects the smaller native root:
`mpqsInternals.Fl_sqrt(4, 5)` returns `2`.


`src/F2v.ts` additionally exports `F2Ms_colelim(M: readonly (readonly number[])[],
rows: number): number[]` and `F2Ms_ker(M: readonly (readonly number[])[], rows:
number): bigint[]`. Each sparse column contains distinct one-based row indices;
the outer array is zero-based. `F2Ms_colelim` returns surviving one-based column
indices. `F2Ms_ker` returns bitsets with bit i for outer column i. Above 640 rows,
it uses native randomized block Lanczos and can return a proper subspace of the
kernel. Set the native RNG seed when the exact returned basis matters:

```ts
const { F2Ms_ker } = require('./src/F2v.ts');
const { setrand } = require('./src/random.ts');
setrand(1n);
F2Ms_ker([[], [], []], 641); // [4n, 6n, 5n]
```

The MPQS test-only `mpqsInternals.F2Ms_ker` keeps its existing one-based
Uint32Array result and delegates to this backend.


The test-only `src/mpqs.ts:mpqsInternals` object also exposes
`mpqs_factorback(handle, packedExponents): bigint` and
`mpqs_check_rel(handle, { Y, relp }, q, mode): void` for native debug comparisons.
The former supports signed packed exponents and signed nonzero moduli; the latter
raises native `PariError` diagnostics when a checked relation is wrong. Mode 0 is
factorization; mode 1 preserves the native class-group early returns.

Its word symbol helpers retain their Number word arguments (exact integers,
positive denominator): `krouu(1, 2 ** 32)` is `1` and `kroiu(-100n, 3)` is `-1`.


`mpqsInternals.Fl_inv(a: number, p: number): number` requires exact integer word
arguments with `0 <= a < p` and positive p. It returns the native modular inverse
or raises native `PariError` for a nonunit. For example, `Fl_inv(3, 97)` returns
`65`, while `Fl_inv(0, 2)` raises
`impossible inverse in Fl_inv: Mod(0, 2).` The valid case `Fl_inv(0, 1)` is `0`.

For source comparisons, `mpqsInternals` also exposes the existing
`mpqs_set_parameters`, `mpqs_create_FB`, `mpqs_sieve_array_ctor`, `mpqs_poly_ctor`,
`mpqs_set_sieve_threshold`, `mpqs_locate_A_range`, `mpqs_self_init`, `mpqs_sieve`
and `mpqs_eval_sieve` helpers. They use the typed handle returned by `newHandle`
and require the corresponding initialized native state; class-group entry
points remain unported.


The internal `src/_mpqs_hash.ts` exports `RelationData` (`{ Y: bigint; relp:
number[] }`), `relationHash(row): bigint`, and `RelationTable<T extends
RelationData = RelationData>`. Construct a table with a nonnegative integer
minimum size, insert immutable relation values with `add(row)`, read `size`,
and iterate `values()` in native bucket order. Duplicate relation values are
ignored. This helper implements the MPQS relation-key subset of PARI hashing.

For example, `relationHash({ Y: 1n, relp: [1048578] })` is
`6171988346495030103n`. With native seed 1, `mpqs(1000509270188293n)` returns
`[[10005089n, 1n], [100000037n, 1n]]` in that order.


When MPQS debugging is enabled, a failed consistency check after Gaussian
elimination emits `console.warn('MPQS: wrong relation found after Gauss')` and
continues, as PARI does. The earlier per-relation checker still raises its native
PariError when it fails.

The test-only `mpqsInternals` exposes `mpqs_eval_cand` and
`mpqs_solve_linear_system`; the latter takes an initialized handle and a
`RelationTable` of full relations. For source diagnostics, a handle with N=143,
the standard factor base [2,3,5,7,11,13], and relation {Y:10n, relp:[]} emits that
warning and returns `[[13n, 1n], [11n, 1n]]`. That deliberately inconsistent
relation is a diagnostic fixture, not a valid factorization relation.

### Seeded Pollard lambda with a custom hash

`groups.discrete_log_lambda(target, base, bounds, operation?, identity?, inverse?, op?, hashFunction?)`
accepts a hash returning `bigint`, including negative values. Step sizes use Sage's
CPython random stream. With a fixed seed and the same custom hash, the selected
logarithm and random walk match Sage, including intervals containing multiple logs.

```ts
import { groups, Mod } from 'sagemath-ts';
import { set_random_seed } from 'sagemath-ts/misc/randstate';
set_random_seed(0);
groups.discrete_log_lambda(
  Mod(300n, 17n), Mod(1n, 17n), [0n, 256n], '+',
  undefined, undefined, undefined, v => -(v.value ** 2n + 1n)
); // 232n
```

The default string hash remains an adapter to JavaScript objects; use an explicit
shared hash for identical native walks. Single-value bounds preserve Sage's
`ZeroDivisionError` after hashing, without random step draws. Both Pollard routines
call an element's `set_immutable()` hook at Sage's original walk positions; lambda
can make the input target immutable. Rho's stored-point lookup verifies element
equality when string keys collide.

`groups.discrete_log_rho(target, base, primeOrder, operation?, identity?, inverse?, op?)`
requires an explicit prime order and uses the modular ring's CPython stream for
random exponents. Its string hash remains the adapter documented in `DEVIATIONS.md`.

```ts
import { groups, Mod } from 'sagemath-ts';
import { set_random_seed } from 'sagemath-ts/misc/randstate';
set_random_seed(0);
const base = Mod(4n, 1019n); // multiplicative order 509
groups.discrete_log_rho(base.pow(250n), base, 509n, '*'); // 250n
```

`groups.order_from_multiple(element, multiple, factorization?, operation?, identity?, inverse?, op?, options?)`
accepts `options: { plist?: IntegerLike[], check?: boolean }`; `check` defaults to
true. Empty factorization and prime lists request automatic factorization.

```ts
import { groups, Mod } from 'sagemath-ts';
groups.order_from_multiple(Mod(1n, 12n), 12n, [], '+'); // 12n
```

`groups.multiples(step, n, start?, indexed?, operation?, op?)` validates and copies
its inputs immediately. Its callback runs before each value is returned; after a
callback error, a later `next()` can continue as in Sage. Copy hooks `__copy__()`
or `copy()` support custom objects with private/internal state; ordinary objects
use a shallow copy. The returned object retains the JavaScript generator protocol.

```ts
import { groups, Mod } from 'sagemath-ts';
let calls = 0;
const iterator = groups.multiples(Mod(2n, 11n), 2n, Mod(5n, 11n), false, 'other',
  (x, y) => { calls++; return x.add(y); });
const first = iterator.next();
if (!first.done) first.value.value; // 5n
calls; // 1
```

`NumberFieldElement.mul(other: NumberFieldElement | RationalLike)` accepts field
products and exact scalar actions; `RationalLike` is `bigint | Integer | Rational`.
Generic group routines recognize both Sage `is_zero()`/`is_one()` predicates and
the port's existing `isZero()`/`isOne()` names.

```ts
import { groups } from 'sagemath-ts';
import { QuadraticField, Rational } from 'sagemath-ts/rings';
const a = QuadraticField.create(2n, 'a').gen();
a.mul(new Rational(7n, 3n)).list().map(String); // ['0', '7/3']
groups.bsgs(a, a.mul(31n), [0n, 64n], '+'); // 31n
```

Number-field `add`, `sub`, `div`, and `eq` also accept
`NumberFieldElement | RationalLike`. The constructor
`K.__call__(x: NumberFieldElement | RationalLike | number)` accepts `Integer`
wrappers; a JavaScript `number` follows Sage's binary64-real-to-rational conversion.
Use `bigint` or `Rational` when the input must specify an exact value.

```ts
import { QuadraticField, Integer, Rational } from 'sagemath-ts/rings';
const K = QuadraticField.create(2n, 'a');
K.gen().add(new Integer(3n)).list().map(String); // ['3', '1']
K.gen().div(3n).list().map(String); // ['0', '1/3']
K.__call__(1.5).eq(new Rational(3n, 2n)); // true
```

`NumberFieldElement.__getitem__(index: number | IntegerLike)` follows Sage's
coefficient-indexing rules. Square-root/gaussian representations permit `-2` and
`-1`; general representations require `0 <= index < degree`. Out-of-range indices
raise `IndexError`. Integral host-number indices use the ordinary array-index
convention; exact bigint and Integer indices are also accepted.

```ts
import { QuadraticField } from 'sagemath-ts/rings';
const a = QuadraticField.create(2n, 'a').gen().add(2n);
a.__getitem__(-2).toString(); // '2'
a.__getitem__(-1n).toString(); // '1'
```

`Integer.add`, `sub`, and `mul` return `Integer` for `IntegerLike` operands and
`Rational` for rational operands, including rationals with denominator one.
`Integer.eq`, `lt`, `le`, `gt`, and `ge` also accept rational operands and compare
exactly. `Integer.div` retains its documented floor-division meaning.

Generic groups require `eq(other: this): boolean`. Standard operation dispatch
recognizes `__invert__()` and `inv()`, and derives missing parent identities from
an element's neutral power/action. Integer multiplicative `multiple` returns
`Integer | Rational`; standard `parseGroupOps` with an Integer sample returns
`GroupOps<Integer | Rational>` because reciprocals belong to QQ.

```ts
import { Integer, Rational, groups } from 'sagemath-ts';
new Integer(2n).add(new Rational(1n, 3n)).toString(); // '7/3'
new Integer(2n).eq(new Rational(2n)); // true
groups.bsgs(new Integer(2n), new Integer(8n), [0n, 10n], '*'); // 3n
const ops = groups.parseGroupOps('*', undefined, undefined, undefined, new Integer(2n));
ops.inverse(new Integer(2n)).toString(); // '1/2'
```

`Rational.eq(number)` compares through Sage's Real Double Field coercion, including
binary64 rounding, underflow, overflow and NaN. IntegerLike and Rational operands
use exact comparison.

```ts
import { Rational } from 'sagemath-ts';
new Rational(1n, 10n).eq(0.1); // true
const q = new Rational(2n ** 100n + 1n);
q.eq(2 ** 100); // true: both round to the same binary64 value
q.eq(2n ** 100n); // false: exact integer comparison
```


`RealDoubleElement.add`, `sub`, `mul`, `div`, and `eq` accept RDF elements,
IntegerLike, Rational, number, and boolean operands. Arithmetic coerces scalars
through RDF. Addition of bigint/bool zero and multiplication by bigint/bool one
reuse the original object, as Sage does for Python integers and booleans.
Generic `multiple` uses Sage's binary group algorithm, including its identity and
idempotence shortcuts; floating-point rounding follows that operation order.

```ts
import { RDF, Rational } from 'sagemath-ts/rings';
import { groups } from 'sagemath-ts';
RDF.__call__(1.5).add(2n).value; // 3.5
RDF.__call__(1.5).mul(new Rational(1n, 3n)).value; // 0.5
Object.is(groups.multiple(RDF.__call__(-0), 2n, '*').value, -0); // true
```


`NumberFieldElement.pow(n: IntegerLike)` accepts bigint and Integer exponents;
exponent one returns the same object. `new NumberFieldElement(K, coefficients)`
reduces a rational coefficient vector modulo K's defining polynomial and pads it
to the degree. `K.polynomial()` and `K.defining_polynomial()` retain the original
leading coefficient; degree-one `K.gen()` is the rational root.

```ts
import { NumberField, NumberFieldElement, RationalPolynomial, Rational, Integer } from 'sagemath-ts/rings';
const f = RationalPolynomial.fromBigInts([-1n, 2n]);
const K = new NumberField(f, 'a');
K.gen().list().map(String); // ['1/2']
K.polynomial() === f; // true
const a = new NumberFieldElement(K, [1n, 2n, 3n].map(c => new Rational(c)));
a.list().map(String); // ['11/4']
a.pow(new Integer(1n)) === a; // true
```

Number-field `minimal_polynomial()`, `absolute_trace()` and `relative_trace()`
are port convenience aliases for `minpoly()` and absolute-field `trace()`.
`numerator()` extends the quadratic Sage method to all supported absolute fields
as `a * a.denominator()`; these extra method names are documented in DEVIATIONS.md.


`NumberField.different()` returns and caches the maximal-order different ideal.
`NumberFieldIdeal.inverse()` supports nonzero fractional ideals in general absolute
fields; `pow(exponent: IntegerLike)` accepts bigint and Integer and returns the
original ideal for exponent one. `AbsoluteOrder.different()` and `codifferent()`
are port conveniences for the ambient maximal order, not a supplied nonmaximal basis.

```ts
import { QuadraticField, Integer } from 'sagemath-ts/rings';
const K = QuadraticField.create(-23n, 'a');
K.different().norm().toString(); // '23'
K.different() === K.different(); // true
const I = K.ideal(2n, K.gen().add(1n));
I.pow(new Integer(1n)) === I; // true
I.mul(I.inverse()).eq(K.ideal(1n)); // true
```

The internal `_pari_ideal_data(): NfIdealData` bridge caches PARI trace/ideal metadata;
`NfIdealData` lives at `@sagemath-ts/parigp-ts/src/base1.js`. Internal
`arith/power.ts` also exports `generic_power_pos<T extends { mul(other: T): T }>(a: T,
n: bigint): T`, requiring a positive exponent. These are dependency interfaces,
not package-root exports.


`NumberFieldIdeal.gens_two(): [NumberFieldElement, NumberFieldElement]` returns a
cached pair, both with the ambient field as parent. The first generates the ideal's
intersection with QQ, so it can be fractional. The second is zero exactly for a
rational principal ideal. **22.0.0 type migration:** the first entry was formerly
`bigint`; use its field-element operations or `.list()[0]` for its rational value.

`integral_basis()` and the port's `zk_basis()` return an HNF basis of the ideal in
the maximal order, with an empty list for the zero ideal. `free_module()` retains
its port record shape `{ basis: NumberFieldElement[], rank: number }` and caches
that record; it does not expose Sage's full FreeModule API.

```ts
import { QuadraticField, Rational } from 'sagemath-ts/rings';
const K = QuadraticField.create(-5n, 'a');
const I = K.ideal(new Rational(3n, 2n));
const [a, b] = I.gens_two();
a.toString(); // '3/2'
b.is_zero(); // true
a.parent() === K; // true
I.free_module() === I.free_module(); // true
K.ideal(0n).integral_basis().length; // 0
```


`NumberFieldIdeal.denominator(): NumberFieldIdeal` and `numerator()` return cached,
coprime integral ideals D and N with I = N/D. **23.0.0 type migration:** denominator
formerly returned bigint; use `I.denominator().smallest_integer()` when you need the
least positive integer clearing the ideal's denominators. Fractional `is_coprime()`
checks both numerator and denominator prime supports. Zero ideals have no Sage
fractional-ideal numerator, denominator, coprimality or divisibility methods; calling
those port methods raises AttributeError.

`K.ideal()` and all-zero generator lists produce a fresh zero ideal with one zero
generator. `intersection()` retains fractional coefficients. Multiplication and
division preserve Sage's generator-count dispatch, including its zero-ideal errors.

```ts
import { QuadraticField, Rational } from 'sagemath-ts/rings';
const K = QuadraticField.create(-1n, 'i');
const I = K.ideal(K.gen().mul(4n).add(3n).div(5n));
I.numerator().norm().toString(); // '5'
I.denominator().norm().toString(); // '5'
I.numerator().div(I.denominator()).eq(I); // true
const half = K.ideal(new Rational(1n, 2n));
half.intersection(half).eq(half); // true
K.ideal().is_zero(); // true
```

`NumberField.ideal(...gens: NumberFieldIdealInput[]): NumberFieldIdeal` and
`fractional_ideal(...gens: NumberFieldIdealInput[]): NumberFieldFractionalIdeal`
accept scalars, generator lists, and existing ideals. The aliases are exported from
`rings/number_field/number_field_ideal.ts` and the number-field barrel:
`NumberFieldIdealGenerator = RationalLike | number | NumberFieldElement`;
`NumberFieldIdealInput = NumberFieldIdealGenerator | NumberFieldIdealGenerator[] |
NumberFieldIdeal`. Host numbers follow the existing floating-point coercion rules.

Nonzero factory results have the fractional class. The fractional factory rejects
zero and preserves the identity of a same-field fractional ideal. Direct ideal
constructors take `(K, gens: (NumberFieldIdealGenerator |
NumberFieldIdealGenerator[])[])`. `intersection(other: NumberFieldIdealInput)`
coerces scalar/list/ideal inputs and delegates to PARI's LLL-kernel intersection.
String prefixes distinguish base and fractional ideals; full Sage generator display
reduction is still incomplete.

```ts
import { QuadraticField, Rational } from 'sagemath-ts/rings';
const K = QuadraticField.create(-1n, 'i');
const I = K.ideal([new Rational(1n, 2n)]);
K.fractional_ideal(I) === I; // true
I.intersection([new Rational(3n, 2n), 3n]).norm().toString(); // '9/4'
I.intersection([]).is_zero(); // true
K.ideal().gens().length; // 1
```

`NumberField.__call__(x: NumberFieldInput)` also accepts an array of exactly
`K.degree()` coefficients. `NumberFieldInput` is exported from the number-field
module/barrel and `sagemath-ts/rings`; it consists of RationalLike, number,
NumberFieldElement, or an array of those scalar types. Vector coefficients must
coerce to rationals. A wrong length raises ValueError.

`NumberFieldIdeal.contains(x: NumberFieldInput | NumberFieldIdeal): boolean` first
coerces x to the field; a TypeError returns false. `add`, `mul`, `div`, `divides`
and `is_coprime` accept `NumberFieldIdealInput`. Addition/multiplication accept
scalar and generator-list operands. Division rejects generator lists and nonzero
host floats; use Rational or Integer for exact scalar division. Foreign ideals
follow each operation's native coercion rules.

```ts
import { NumberField, RationalPolynomial, Rational } from 'sagemath-ts/rings';
const K = new NumberField(RationalPolynomial.fromBigInts([3n, 0n, 1n]), 'a');
K.__call__([1n, new Rational(1n, 2n)]).toString(); // '1/2*a + 1'
K._pari_integral_basis().map(String); // ['1', '1/2*a - 1/2']
const I = K.ideal(2n);
I.contains(2n); // true
I.contains(new Rational(1n, 2n)); // false
I.divides(6n); // true
I.is_coprime(3n); // true
I.add([3n]).norm().toString(); // '1'
I.div(new Rational(1n, 2n)).norm().toString(); // '16'
```

The internal PARI dependency exports `diviiround(x: bigint, y: bigint): bigint`
from `src/gen3.js` (nearest quotient; half ties toward positive infinity),
`ZM_hnfcenter(input: bigint[][]): bigint[][]` from `src/hnf_snf.js` (copies square
HNF columns and centers residues), and `idealmul(nf: NfIdealData, I: bigint[][],
iDen: bigint, J: bigint[][], jDen: bigint): [bigint[][], bigint]` from `src/base4.js`.
The product adapter uses the rational column representation below and checks native
matrix dimensions before removing primitive content. These are dependency interfaces.

The internal PARI `src/base4.js` adapter also exports `idealintersect(nf: NfIdealData,
I: bigint[][], iDen: bigint, J: bigint[][], jDen: bigint): [bigint[][], bigint]`,
using the same column/denominator representation as idealdiv below. It preserves
input matrices and normalizes the output denominator. Malformed shapes or nonpositive
denominators raise RangeError.

The internal PARI `src/base4.js` adapter exports `idealdiv(nf: NfIdealData,
I: bigint[][], iDen: bigint, J: bigint[][], jDen: bigint): [bigint[][], bigint]`.
It uses integral HNF columns plus positive denominators; `[]` represents the zero
ideal. The helper is a dependency interface, not a root package export.

The internal PARI `src/lll.js` dependency exports
`fplll_fast(B: bigint[][], delta = 0.99, eta = 0.51, keepfirst = false,
trackTransform = true): [number, bigint[][], bigint[][] | null]`. Matrices contain
columns. It copies B, starts U at the identity or null, and returns status plus B/U,
including partial changes on status -1. Inputs must be nonempty and rectangular.
This is the uncertified fast stage, with the AArch64 Clang floating profile documented
in DESIGN.md; it is not a complete LLL reducer and can loop on dependent keep-first
bases. It is not a package-root export.

The internal `src/_binary64.js` dependency exports `fma(a: number, b: number,
c: number): number`, `ldexp(x: number, exponent: number): number`, and
`frexp(x: number): [number, number]`. They preserve single-rounding binary64 fused
arithmetic, power-of-two scaling and mantissa/exponent decomposition. `ldexp` takes
an integer exponent; these helpers do not expose native floating exception flags.

```ts
import { fplll_fast } from '@sagemath-ts/parigp-ts/src/lll.js';
const [status, B, U] = fplll_fast([[1n, 1n], [0n, 2n]]);
status; // 0
B; // [[1n, 1n], [-1n, 1n]]
U; // [[1n, 0n], [-1n, 1n]]
```

The internal PARI `src/lll.js` dependency also exports
`fplll_dpe(B: bigint[][] | null, G: bigint[][] | null = null, delta = 0.99,
eta = 0.51, keepfirst = false, trackTransform = true, wantNorms = false):
[number, bigint[][] | null, bigint[][] | null, bigint[][] | null, MpReal[] | null]`.
It copies inputs and returns `[status, G, B, U, norms]`. A null G requests incremental
Gram construction; B may be null with a supplied square Gram. Only the updated
triangle `G[j][i]` for `i <= j` is authoritative. Incremental failure returns G=null;
a supplied Gram retains partial state. Requested norms use PARI `MpReal` records.

This is the DPE stage, not the complete LLL wrapper. Dependent keep-first bases can
loop or fail; unrepresentable Gram-Schmidt shifts raise the RangeError documented in
DEVIATIONS.md. The adapter requires nonempty compatible matrix shapes.

```ts
import { fplll_dpe } from '@sagemath-ts/parigp-ts/src/lll.js';
const [status, G, B, U] = fplll_dpe([[1n, 1n], [0n, 2n]]);
status; // 0
G; // [[2n, 0n], [0n, 2n]]
B; // [[1n, 1n], [-1n, 1n]]
U; // [[1n, 0n], [-1n, 1n]]
```

The internal PARI `src/lll.js` dependency exports two real-arithmetic stages:

- `fplll_heuristic(B: bigint[][], delta = 0.99, eta = 0.51, keepfirst = false,
  trackTransform = true, precision = 64, gramPrecision = precision):
  [number, bigint[][], bigint[][] | null]` returns `[status, B, U]`.
- `fplll(B: bigint[][] | null, G: bigint[][] | null = null, delta = 0.99,
  eta = 0.51, keepfirst = false, trackTransform = true, wantNorms = false,
  precision = 64): [number, bigint[][] | null, bigint[][] | null,
  bigint[][] | null, MpReal[] | null]` returns `[status, G, B, U, norms]`.

Precision counts bits and must be a positive multiple of 64. Inputs are copied and
matrices contain columns. Status -1 retains partial state; the proved stage's Gram
and optional norms follow the DPE adapter's conventions, with norms at the requested
real precision. These stages are not the complete LLL wrapper and can fail or loop
on dependent keep-first bases.

```ts
import { fplll, fplll_heuristic } from '@sagemath-ts/parigp-ts/src/lll.js';
const B = [[1n, 1n], [0n, 2n]];
fplll_heuristic(B, 0.99, 0.51, false, true, 128, 256)[1]; // [[1n, 1n], [-1n, 1n]]
fplll(B, null, 0.99, 0.51, false, true, true, 128)[4]?.map(r => r.p); // [128, 128]
```

Internal `src/gen3.js` exports `roundr_safe(x: MpReal): bigint`, rounding half-integer
ties toward positive infinity and allowing integer results beyond the available
mantissa precision. `src/kernel/none/cmp.js` exports `abscmprr(x: MpReal,
y: MpReal): number`, returning -1/0/1 for absolute-value comparison; a zero is zero
regardless of its accuracy exponent. Neither helper is a package-root export.


Internal `src/bibli1.js` exports `QrScalar = bigint | MpReal` and
`Householder = [MpReal, QrScalar[]]`, plus these column-matrix dependencies:

- `QR_init(x: QrScalar[][], precision = 64): [number, QrScalar[] | null,
  Householder[] | null, QrScalar[][] | null]` returns `[success, B, Q, L]`.
  B contains squared norms, Q has n-1 reflectors, and L is lower triangular.
  On precision failure success is zero and all three outputs are null.
- `R_from_QR(x: QrScalar[][], precision = 64): QrScalar[][] | null`
  returns the transpose of L.
- `gaussred_from_QR(x: QrScalar[][], precision = 64): QrScalar[][] | null`
  returns the native Gaussian reduction of x^T*x; input must be square and full rank.

Inputs must be nonempty rectangular columns with rows >= columns. Working precision
counts bits and must be a positive multiple of 64. Real input cells retain their
own precision except for native explicit conversions. Inputs are preserved. These
are internal dependencies, not root exports or the generic matqr API.

```ts
import { QR_init, R_from_QR, gaussred_from_QR } from '@sagemath-ts/parigp-ts/src/bibli1.js';
const I = [[1n, 0n], [0n, 1n]];
QR_init(I)[1]; // [1n, 1n]
R_from_QR(I)?.[1]?.[1]; // 1n
gaussred_from_QR(I)?.[0]; // [1n, 0n]
```


Internal `src/_real_matrix.js` defines `MatrixReal = MpReal | 0n`, where 0n is an
exact integer zero. The real FLATTER dependencies use zero-indexed columns:

- `src/alglin2.js`: `qfgaussred_positive(a: MatrixReal[][]): MatrixReal[][] | null`
  reads the upper triangle and returns null on a nonpositive pivot.
- `src/alglin2.js`: `RgM_Cholesky(M: MatrixReal[][], precision = 64):
  MatrixReal[][] | null` returns the native Cholesky factor, or null. Working precision
  counts bits and must be a positive multiple of 64; real cells determine actual
  square-root precision.
- `src/alglin1.js`: `RgM_inv_upper(A: MatrixReal[][]): MatrixReal[][]`
  performs backward substitution on an invertible upper-triangular matrix.

These preserve input state and accept empty matrices. They implement the real/zero
coefficient domain, not generic rational or complex matrices, and are not root exports.

```ts
import { RgM_Cholesky } from '@sagemath-ts/parigp-ts/src/alglin2.js';
import { RgM_inv_upper } from '@sagemath-ts/parigp-ts/src/alglin1.js';
import { itor } from '@sagemath-ts/parigp-ts/src/qfb.js';
const D = [[4n, 0n], [0n, 9n]].map(c => c.map(x => itor(x, 128)));
RgM_Cholesky(D)?.[0]?.[0]; // itor(2n, 128)
(RgM_inv_upper(D)[0]?.[0] as ReturnType<typeof itor>).e; // -2
```

Internal `src/_real_matrix.js` also exports `matrixRealMul(x, y)`,
`matrixRealSub(x, y)`, `matrixRealNeg(x)`, `matrixRealDiv(x, y)` and
`matrixRealInv(x)`. Arguments/results are MatrixReal, except matrixRealInv returns
MpReal. They implement only those PARI numeric branches, including native exact-zero
and finite-accuracy-zero error dispatch. Large allocated-zero inverses raise the
existing deterministic error documented in DEVIATIONS.md.


Internal `src/RgV.js` defines `RgScalar = bigint | MpReal` and exports:

- `RgV_dotsquare(x: RgScalar[]): RgScalar`
- `RgV_dotproduct(x: RgScalar[], y: RgScalar[]): RgScalar`
- `gram_matrix(x: RgScalar[][]): RgScalar[][]`
- `RgM_mul(x: RgScalar[][], y: RgScalar[][]): RgScalar[][]`

Matrices use zero-indexed columns. These preserve inputs, native integer/real
arithmetic order and integer multiplication delegation. Empty vectors have dot
product zero; dimensions must be compatible. The adapters support integer/real
coefficients, not all generic PARI types, and are not root exports.

Internal `src/polarit2.js` also exports `RgM_rescale_to_int(x: (bigint | MpReal)[][]):
bigint[][]`. It preserves all-integer input values. In mixed input, native scaling
can round integer entries too. An all-zero real matrix with negative accuracy raises
`PariError('overflow in expo()')`, matching PARI's sentinel shift.

```ts
import { RgV_dotsquare, RgV_dotproduct, gram_matrix, RgM_mul } from '@sagemath-ts/parigp-ts/src/RgV.js';
import { RgM_rescale_to_int } from '@sagemath-ts/parigp-ts/src/polarit2.js';
import { itor } from '@sagemath-ts/parigp-ts/src/qfb.js';
RgV_dotsquare([1n, 2n]); // 5n
RgV_dotproduct([1n, 2n], [3n, 4n]); // 11n
gram_matrix([[1n, 2n], [3n, 4n]]); // [[5n, 11n], [11n, 25n]]
RgM_mul([[1n, 2n], [3n, 4n]], [[1n, 0n], [0n, 1n]]); // [[1n, 2n], [3n, 4n]]
RgM_rescale_to_int([[3n, itor(16n, 64)]]); // [[2n, 8n]]
```


Internal `src/lll.js` re-exports these precision-selection dependencies from
`src/_lll_gso.js` (all matrices contain zero-indexed columns):

- `drop(R: QrScalar[][]): bigint`, `potential(R: QrScalar[][]): bigint`, and
  `spread(R: QrScalar[][]): bigint` return native diagonal-exponent statistics.
  drop requires matching integer/real types on the diagonal.
- `condition_bound(U: QrScalar[][], lower = false): bigint` and
  `GS_extraprec(L: QrScalar[][], lower = false): bigint` compute native condition
  and precision bounds for an invertible triangular matrix.
- `gramschmidt_upper(M: bigint[][]): MpReal[][]` converts an invertible upper
  triangular integer matrix at its native selected precision.
- `gramschmidt_dynprec(M: bigint[][]): QrScalar[][]` returns adaptive QR's upper R
  for an integer basis with full column rank.
- `RgM_Cholesky_dynprec(M: bigint[][]): MatrixReal[][]` returns adaptive Cholesky
  for a positive-definite integer Gram matrix.

Inputs are preserved and must be nonempty with rows >= columns; Cholesky requires
a square matrix. Full rank/positive definiteness are caller preconditions. These
routines preserve native precision escalation, with a RangeError for allocation
precision outside the exact JavaScript integer range. They are not root exports or
a complete LLL reducer.

```ts
import { drop, potential, spread, condition_bound, GS_extraprec,
  gramschmidt_upper, gramschmidt_dynprec, RgM_Cholesky_dynprec } from '@sagemath-ts/parigp-ts/src/lll.js';
const I = [[1n, 0n], [0n, 1n]];
[drop(I), potential(I), spread(I), condition_bound(I), GS_extraprec(I)]; // [0n, 0n, 0n, 0n, 4n]
gramschmidt_upper(I)[0]?.[0]?.p; // 64
(gramschmidt_dynprec(I)[0]?.[0] as { p: number }).p; // 64
(RgM_Cholesky_dynprec(I)[0]?.[0] as { p: number }).p; // 64
```


Internal modular linear algebra dependencies use an unused slot zero on **both**
matrix axes, as in `[[], [0n, a11, a21], [0n, a12, a22]]`:

- `src/Flv.js`: `Flm_pivots(x: bigint[][], p: bigint): [number[] | null, number]`
  returns `[pivotRows, nullity]`; `Flm_gauss(a: bigint[][], b: bigint[][],
  p: bigint): bigint[][] | null` returns a solve or null on rank failure.
  p must be a prime in the unsigned 64-bit range.
- `src/alglin1.js`: `ZM_pivots(x: bigint[][]): [number[] | null, number]` and
  `ZM_rank(x: bigint[][]): number` compute certified integer rank profiles/rank.
- `src/alglin1.js`: `ZM_gauss(a: bigint[][], b: bigint[][]):
  [bigint[][], bigint] | null` returns `[numerator, positiveCommonDenominator]`
  or null on rank failure. It accepts matrix RHS values and requires rows >= columns.

Nonnull public pivot vectors also have unused slot zero. Their other entries are
one-based pivot row indices, or zero for dependent columns. Inputs are preserved.
Overdetermined solves use selected independent rows; they do not verify every RHS
row. The native rank-certificate error for certain common trial-prime factors is
preserved. These helpers are internal subpath exports, not package-root exports.

```ts
import { Flm_pivots, Flm_gauss } from '@sagemath-ts/parigp-ts/src/Flv.js';
import { ZM_pivots, ZM_rank, ZM_gauss } from '@sagemath-ts/parigp-ts/src/alglin1.js';
const A = [[], [0n, 1n, 0n], [0n, 0n, 2n]], B = [[], [0n, 3n, 3n]];
Flm_pivots(A, 101n); // [[0, 1, 2], 0]
Flm_gauss(A, B, 101n); // [[], [0n, 3n, 52n]]
ZM_pivots(A); // [[0, 1, 2], 0]
ZM_rank(A); // 2
ZM_gauss(A, B); // [[[], [0n, 6n, 3n]], 2n]
```

`src/_matrix_inverse.js` additionally exports row-oriented internal kernels (no dummy
matrix slots): `wordMatrixPivots(A, m, n, p)` and `integerMatrixPivots(A, m, n)` return
`[number[] | null, number]`. Their pivot vectors contain n entries with one-based row
values and zero for dependent columns, with no dummy vector slot.
`wordMatrixSolve(A, B, m, n, k, p): bigint[][] | null` and
`integerMatrixSolve(A, B, m, n, k): [bigint[][], bigint] | null` take A with m rows/n
columns and B with m rows/k columns; integer results are numerator rows/denominator.
All size arguments are numbers and p is bigint. These kernels serve the public
column adapters above and retain their numeric domain and solve preconditions.


Internal `src/lll.js` also re-exports the adaptive wrapper from `src/_lll_wrapper.js`.
Here `Matrix = bigint[][]` and all matrices are zero-indexed **columns**:

- `ZM_lll(input: Matrix, delta = 0.99, flag = LLL_IM): LllResult` returns the native
  transformation/image, kernel, [kernel,image] pair, or reduced basis for INPLACE.
  `LllResult = Matrix | [Matrix, Matrix]`; choose the result shape using flags.
- `ZM_lll_norms(input: Matrix, delta = 0.99, flag = LLL_IM):
  [LllResult, MpReal[] | null]` also returns native squared Gram–Schmidt norms.
  Trivial and successful NOCERTIFY shortcuts leave norms null.
- `lllfp(input: (bigint | MpReal)[][], delta = 0.99, flag = LLL_IM)` additionally
  rescales integer/real input. Its result is LllResult, an inexact matrix/pair for
  a one-column shortcut, or null for failed inexact Gram Cholesky. Inexact Gram
  coefficients must be MpReal or exact zero; general mixed Gram coefficients are
  not yet supported.
- `flat(M: Matrix, flag = 0): [Matrix, Matrix | null, bigint, bigint]` returns one
  compression step's [basis, transformation, drop, potential]. INPLACE omits U.
- `ZM_flatter(M: Matrix, flag = 0): Matrix | null` returns accumulated U (null if
  no accepted step), or the basis for INPLACE.
- `ZM_flatter_rank(M: Matrix, rank: number, flag = 0): Matrix | null` includes
  rank augmentation. Supply the actual integer rank.
- `ZM_flattergram(M: Matrix, flag = 0): Matrix | null` compresses a positive-definite
  integer Gram matrix; the native implementation can return an identity matrix.

Flags are `LLL_KER=1`, `LLL_IM=2`, `LLL_ALL=4`, `LLL_GRAM=0x100`,
`LLL_KEEP_FIRST=0x200`, `LLL_INPLACE=0x400`, `LLL_COMPATIBLE=0x800`,
`LLL_UPPER=0x1000`, `LLL_NOCERTIFY=0x2000`, and `LLL_NOFLATTER=0x4000`.
Inputs are preserved. Direct flat/FLATTER routines require at least two columns;
flat/ZM_flatter require full column rank. The rank-deficient Gram route has a
64-augmentation resource limit; see DEVIATIONS.md. Requesting norms on the native
positive-definite 2×2 image shortcut raises its caught PariError; use INPLACE
without the IM bit to follow the general basis route. Native low-accuracy norm
records (including zero precision) are preserved, not mathematically corrected.

```ts
import { ZM_lll_norms, LLL_INPLACE } from '@sagemath-ts/parigp-ts/src/lll.js';
import { mpreal_to_frac } from '@sagemath-ts/parigp-ts/src/qfb.js';
const [NB, NN] = ZM_lll_norms([[1n, 0n], [0n, 1n]], 0.99, LLL_INPLACE);
NB; // [[1n, 0n], [0n, 1n]]
NN!.map(mpreal_to_frac); // [[9223372036854775808n, 9223372036854775808n], [9223372036854775808n, 9223372036854775808n]]
ZM_lll_norms([[1n]])[1]; // null
```

Internal factorization dependencies use zero-indexed polynomial coefficients and
zero-indexed matrix **columns**:

- `src/hnf_snf.js`: `ZM_hnfperm(A, withU = false, withPerm = false)` returns
  `[H, U | null, permutation | null]`. U is a full unimodular transformation;
  its initial columns span the kernel. The one-based row permutation maps the
  remaining columns of A*U to H. `hnfperm(A)` returns all three outputs.
  `ZM_hnf_knapsack(A): bigint[][] | null` tests the permuted HNF and restores row
  order. Its source accepts ±1 entries and empty rows; columns retain native order.
  Zero-rank inputs preserve the bundled source's caught PariError.
- `src/QX_factor.js`: `LLL_check_progress(Bnorm: MpReal, n0: number, A, final: boolean)`
  returns retained basis columns truncated to n0 rows, or null for an irreducibility
  decision. n0 must be within the input row count. Native invalid progress/no-factor
  errors are preserved. The debug timing output is omitted.
- That module also exports `Mignotte_bound(P): bigint | MpReal`,
  `Beauzamy_bound(P): MpReal`, `factor_bound(P): bigint`, and `root_bound(P): bigint`.
  P must be normalized and have positive degree. Bounds retain the native real
  operation order, type-sensitive comparisons and conservative integer rounding.
- `src/bibli2.js`: `vecbinomial(n: number): bigint[]`, for nonnegative integer n.
- `src/ZX.js`: `ZX_Z_eval(P: readonly bigint[], x: bigint): bigint`, using sparse
  coefficient gaps to avoid repeated polynomial multiplications.
- `src/gen3.js`: `ceil_safe(x: bigint | [bigint, bigint] | MpReal): bigint`.
  Rational pairs are normalized with native integer-division construction errors.
  A real value's accuracy affects this upper bound: exact real 2 returns 3, while
  integer 2 returns 2. This is PARI's conservative ceiling, not ordinary rounding.
- `src/trans1.js`: `powruhalf(x: MpReal, n: bigint): MpReal | MpComplex` computes
  x^(n/2) with native word-exponent powering. n is an unsigned 64-bit integer.
  Nonnegative-sign inputs have an MpReal overload. Negative x and odd n return complex.
- `src/kernel/none/level1.js`: `cmpir(x: bigint, y: MpReal): number` and
  `cmpri(x: MpReal, y: bigint): number` return -1, 0 or 1. Native comparison rounds
  the integer to the real's precision and respects finite-accuracy real zeros.

```ts
import { hnfperm, ZM_hnf_knapsack } from '@sagemath-ts/parigp-ts/src/hnf_snf.js';
import { LLL_check_progress, factor_bound, root_bound } from '@sagemath-ts/parigp-ts/src/QX_factor.js';
import { itor } from '@sagemath-ts/parigp-ts/src/qfb.js';
const HA = [[1n,0n],[2n,1n]];
hnfperm(HA); // [[[1n,0n],[0n,1n]],[[-2n,1n],[1n,0n]],[2,1]]
ZM_hnf_knapsack(HA); // [[0n,1n],[1n,0n]]
LLL_check_progress(itor(3n,64), 2, [[1n,0n,0n],[0n,1n,0n],[0n,0n,2n]], true);
// [[1n,0n],[0n,1n]]
factor_bound([-2n,0n,1n]); // 4n
root_bound([-2n,0n,1n]);   // 2n
```

```ts
import { vecbinomial } from '@sagemath-ts/parigp-ts/src/bibli2.js';
import { ZX_Z_eval } from '@sagemath-ts/parigp-ts/src/ZX.js';
import { ceil_safe } from '@sagemath-ts/parigp-ts/src/gen3.js';
import { powruhalf } from '@sagemath-ts/parigp-ts/src/trans1.js';
import { itor } from '@sagemath-ts/parigp-ts/src/qfb.js';
vecbinomial(5); // [1n,5n,10n,10n,5n,1n]
ZX_Z_eval([-2n,0n,0n,1n], 3n); // 25n
ceil_safe(2n);          // 2n
ceil_safe(itor(2n,64)); // 3n
ceil_safe([3n,2n]);    // 2n
powruhalf(itor(-4n,64), 1n); // {re:0n, im:{s:1,e:1,m:9223372036854775808n,p:64}}
```

Internal `src/qfb.js` exports `redimagsl2(q: Qfb): { Q: Qfb; U: bigint[][] }` for
native definite-form reduction. U is **row-oriented**, like qfbredsl2's U. Unlike
the intentionally corrected qfbredsl2 entry, this compatibility entry preserves
PARI's large-input negative-b sign behavior and may return determinant -1.

```ts
import { ZM_lll, lllfp, flat, ZM_flatter, ZM_flatter_rank,
  ZM_flattergram, LLL_ALL } from '@sagemath-ts/parigp-ts/src/lll.js';
import { itor, mkqfb, redimagsl2 } from '@sagemath-ts/parigp-ts/src/qfb.js';
const I = [[1n, 0n], [0n, 1n]], B = [[1n, 0n], [3n, 1n]];
ZM_lll(B); // [[1n, 0n], [-3n, 1n]]
ZM_lll(B, 0.99, LLL_ALL); // [[], [[1n, 0n], [-3n, 1n]]]
lllfp(B.map(c => c.map(x => itor(x, 64)))); // [[1n, 0n], [-3n, 1n]]
flat(I); // [I, I, 0n, 0n]
ZM_flatter(I); // null
ZM_flatter_rank(I, 2); // null
ZM_flattergram(I); // I
redimagsl2(mkqfb(1n, 6n, 10n, -4n)).U; // [[1n, -3n], [0n, 1n]]
```


### PARI Gram transforms and Hermite kernel columns

`lllgramint(G, delta?: {n: bigint, d: bigint}): bigint[][]` is exported by
`@sagemath-ts/parigp-ts`. Matrices are zero-indexed columns. The result is the native
image transformation, including fewer columns on rank-deficient input; callers
must check its width when they require a square transformation. The optional
legacy rational delta maps to native double precision; its default is 99/100.

`ZM_hnflll(A, wantB: boolean, remove: boolean)` returns `{H, B}` and uses dummy
slot zero on both matrix axes. `remove=true` removes zero HNF columns while keeping
the full transformation, including kernel columns. `B` is null if wantB=false.

```typescript
import {lllgramint, ZM_hnflll} from '@sagemath-ts/parigp-ts';
lllgramint([[1n,1n],[1n,1n]]); // [[1n,0n]]
lllgramint([[0n]]); // []
const {H,B}=ZM_hnflll([[],[0n,0n,0n],[0n,0n,0n]],true,true);
// H = [[]]; B = [[],[0n,1n,0n],[0n,0n,1n]]
```

Quaternion ideal `reduced_basis()` applies Sage's final-column sign normalization
so the basis transformation has determinant +1. `quadratic_form()` preserves the
native odd-diagonal error when the original normalized matrix remains over ZZ;
`minimal_element()` and theta-series callers preserve that error too.


### Direct base ideals and fractional methods

Use `K.ideal(gens)` for normal ideal construction. Nonzero results, including
`K.decomposition(p)` entries, are `NumberFieldFractionalIdeal` instances. Direct
`new NumberFieldIdeal(K, gens)` constructs Sage's base class even if gens are
nonzero: inversion and fractional-only methods raise native errors on that class.
The shared TypeScript interface does not imply that those operations are supported
on every instance.

```typescript
import {NumberField, RationalPolynomial, NumberFieldIdeal,
        NumberFieldFractionalIdeal} from 'sagemath-ts/rings';
const K=new NumberField(RationalPolynomial.fromBigInts([-2n,0n,1n]),'a');
const base=new NumberFieldIdeal(K,[K.gen()]);
// base.inverse() raises TypeError: bad operand type for unary ~: 'NumberFieldIdeal'
const P=K.decomposition(2n)[0][0];
P instanceof NumberFieldFractionalIdeal; // true
P.ramification_index(); // 2n
P.residue_class_degree(); // 1n
P.inverse().mul(P).norm().toString(); // '1'
```


### Number-field ideal valuations

Use `I.valuation(P)`, where I is an integral or fractional ideal and P is a prime
ideal (or an input convertible to one). This corrects the old reversed argument
meaning in version 24.0.0. The return type is `bigint | 'Infinity'`; zero ideals
have infinite valuation. `prime_below()` returns the rational prime below P.

```typescript
import { NumberField, RationalPolynomial } from 'sagemath-ts/rings';
const VK = new NumberField(RationalPolynomial.fromBigInts([-2n, 0n, 1n]), 'a');
const VP = VK.decomposition(2n)[0]![0];
const VI = VK.ideal(VK.gen().add(1n).div(VK.__call__(2n)));
VI.valuation(VP);                  // -2n
VK.ideal(0n).valuation(VP);        // 'Infinity'
VK.ideal(8n).valuation(VP);        // 6n
VP.prime_below();                 // 2n
```

The internal PARI imports `src/base2.js`, `src/base3.js` and `src/base4.js` expose
`primedec_end(nf, primeHNFs, p)`, `ZC_nfval(nonzeroVector, prime)` and
`idealval(nf, integerHNF, positiveDenominator, prime)`. Matrices are zero-indexed
columns. Prime data is `{ p, generator, e, f, tau }`, with bigint scalars, a bigint
vector generator and a multiplication matrix tau (or `1n` for inert primes).
`nfmaxord_ideal_data` supplies all NfIdealData fields, including `polynomial`,
`basis`, `basisDenominator` and `discriminant`. Direct low-level callers must supply
certified prime HNFs and a valid maximal-order multiplication table; ZC_nfval rejects
zero vectors. These are internal representation adapters, not Sage public APIs.


### Number-field ideal factorization

`K.factor(input: NumberFieldIdealInput)` accepts an ideal, scalar, field element,
or generator list. `I.factor()` returns a cached array of `[primeIdeal, bigintExponent]`.
Fractional denominator primes have negative exponents. Factors follow Sage's HNF
matrix order in the chosen integral basis; PARI's internal prime order differs.

```typescript
import { NumberField, RationalPolynomial, Rational } from 'sagemath-ts/rings';
const FK = new NumberField(RationalPolynomial.fromBigInts([1n, 0n, 1n]), 'i');
const FI = FK.ideal(FK.gen().mul(4n).add(3n).div(FK.__call__(5n)));
FI.norm().toString();                         // '1'
FI.factor().map(([P, e]) => [P.norm().toString(), e]); // [['5', 1n], ['5', -1n]]
FK.factor(new Rational(1n, 2n)).map(([, e]) => e);     // [-2n]
FI.factor() === FK.factor(FI);                // true
```

Internal PARI base2 exports `idealprimedec(nf,p)`,
`idealprimedec_limit_f(nf,p,limit)` (zero means unlimited),
`idealprimedec_kummer(nf,polynomial,e,p)` and `cmp_prime_ideal(P,Q)`.
Base4 exports `pr_hnf(nf,P)`, `idealfactor(nf,H,d)`, `idealnumden(nf,H,d)`,
`idealismaximal(nf,H,d)`, `idealhnf_principal(nf,coordinates,d)` and
`idealadd(nf,H,d,J,e)`. All use zero-indexed columns and positive denominators;
idealismaximal returns prime data or null. The principal-HNF and addition routines
return `[H,d]`; idealnumden returns `[numeratorHNF,denominatorHNF|null]`.
Base3 exports `zkmultable_inv(M): [bigint[],bigint]` and
`zkmultable_capZ(M): bigint` for nonsingular integral multiplication tables.
`NfIdealData.index?: bigint` is supplied by the factory; older direct data records
can omit it and use the exact basis-determinant fallback.

The internal alglin1 FpM image/suppl/inverse APIs use columns with unused slot zero
on both axes. `hnf_snf.ZM_hnfmodprime(columns,p)` omits those slots.
F2v's `F2m_gauss_pivot(columns,rows)` returns `[oneBasedPivots|null,nullity]`;
`F2m_gauss(A,B,rows)` returns packed BigInt solution columns or null.


### Absolute number-field trace and norm

No-argument `trace()` and `norm()` return Rational values. Quadratic elements use
direct formulas; other absolute fields delegate to PARI's dedicated trace/norm
routines. `absolute_norm()` and `relative_norm()` have the same value over these
absolute fields. `norm(K?: RationalField | NumberField | null)` also accepts QQ,
the element's own field, and degree-one fields. The own-field result is the
original NumberFieldElement, except in degree one where Sage returns Rational.
Generic `trace(K?: RationalField | NumberField | null)` returns an element of K
when K is a number field of degree one or the parent itself. No-argument/QQ/null
results are Rational. Quadratic `trace` rejects every explicit argument (including
null), matching Sage's quadratic override. Proper intermediate fields and
embeddings still need the relative-matrix dependencies.

```typescript
import { QQ, NumberField, RationalPolynomial } from 'sagemath-ts/rings';
const TB = new NumberField(RationalPolynomial.fromBigInts([-2n,0n,0n,1n]), 'a');
const tb = TB.gen().add(1n);
tb.trace(QQ).toString(); // '3'
tb.norm(QQ).toString();  // '3'
tb.trace(TB) === tb;     // true
tb.norm(TB) === tb;      // true
const T1 = new NumberField(RationalPolynomial.fromBigInts([-3n,1n]), 'b');
tb.trace(T1).parent() === T1; // true
tb.norm(T1).toString();      // '3'
```

```typescript
import { NumberField, RationalPolynomial } from 'sagemath-ts/rings';
const TQ = new NumberField(RationalPolynomial.fromBigInts([1n,-1n,2n]), 'a');
TQ.gen().trace().toString();          // '1/2'
TQ.gen().norm().toString();           // '1/2'
TQ.gen().add(1n).trace().toString();   // '5/2'
TQ.gen().add(1n).norm().toString();    // '2'
```

Internal PARI QQ polynomials use `RationalPolynomialData = [bigint[], bigint]`
(coefficients, common denominator); scalar RationalPair results are reduced
`[bigint,bigint]`. Input denominators are nonzero; output denominators are positive.

```typescript
import { RgXQ_trace, RgXQ_norm } from '@sagemath-ts/parigp-ts/src/RgX.js';
RgXQ_trace([[1n,1n],2n], [[-2n,0n,0n,1n],1n]); // [3n,2n]
RgXQ_norm([[1n,1n],2n], [[-2n,0n,0n,1n],1n]);  // [3n,8n]
```

RgX also exports `RgX_rem(x,y)` and `RgXQ_mul(x,y,T)` on those QQ tuples.
ZX exports `ZX_rem(integerCoefficients,monicDivisor)`, `QX_mul(x,y)` and
`QX_ZX_rem(x,monicDivisor)`. Polarit3 exports `QXQ_mul(x,y,monicModulus)` and
`QXQ_norm(x,integerModulus)`. Alglin2's `gnorm` and `gtrace` accept
`RationalPolMod = {value:RationalPolynomialData,modulus:RationalPolynomialData}`;
the value must already be reduced modulo a nonconstant modulus. These entries
cover the rational polynomial-modulus branch, not general PARI GEN values.


### PARI bounded integer-factor recombination

Internal imports from `@sagemath-ts/parigp-ts/src/QX_factor.js` provide:

- `ZX_divides_i(x, y, B = null): bigint[] | null` and `ZX_divides(x, y)`:
  exact integer quotient, or null. The normalized divisor must have positive
  degree. The optional coefficient bound skips the leading quotient coefficient.
- `cmbf_precs(q, A, B): [flag, a, b, qToA, qToB]`: native lifting precision;
  positive bounds, with `2n <= q < 2n ** 31n`.
- `cmbf(pol, famod, bound, p, a, b, klim)`:
  `[bigint[][], bigint[][][], number, boolean]`, containing factors, their modular
  groups, maximum subset size and completion. Inputs must satisfy the native
  squarefree/congruence/trace-bound contract; `false` leaves a group for lattice
  recombination. All polynomial arrays use ascending coefficients and are copied.

`cmbf_maxK(nb)` from `@sagemath-ts/parigp-ts/src/nffactor.js` gives the native subset cap.
`FpXV_prod(values, p)` from `@sagemath-ts/parigp-ts/src/FpX.js` returns a polynomial array,
except that an empty product is `1n`; singleton coefficients retain their raw values.
`centermodii(x, p, half)` from `@sagemath-ts/parigp-ts/src/polarit2.js` retains negative
half ties. A null half selects the source's uncentered remainder adjustment.

```typescript
import { ZX_divides_i, ZX_divides, cmbf_precs, cmbf } from '@sagemath-ts/parigp-ts/src/QX_factor.js';
ZX_divides_i([2n, 2n], [1n, 1n], 0n); // [2n]: leading quotient is not bound-checked
ZX_divides([-1n, 0n, 1n], [1n, 1n]); // [-1n, 1n]
cmbf_precs(5n, 100n, 100n); // [1, 16, 3, 152587890625n, 125n]
cmbf([-1n, 0n, 1n], [[1n, 1n], [124n, 1n]], 3n, 5n, 3, 1, 1);
// [[[1n, 1n], [-1n, 1n]], [[[1n, 1n]], [[124n, 1n]]], 0, true]
```

```typescript
import { cmbf_maxK } from '@sagemath-ts/parigp-ts/src/nffactor.js';
import { FpXV_prod } from '@sagemath-ts/parigp-ts/src/FpX.js';
import { centermodii } from '@sagemath-ts/parigp-ts/src/polarit2.js';
cmbf_maxK(11); // 3
FpXV_prod([[1n, 1n], [4n, 1n]], 5n); // [4n, 0n, 1n]
FpXV_prod([], 5n); // 1n
centermodii(-5n, 10n, 5n); // -5n
```


### PARI general integer Hermite forms

`@sagemath-ts/parigp-ts/src/hnf_snf.js` also exports these column-matrix entries:

- `ZM_hnf(A): bigint[][]`: native dispatch at seven/eight input columns.
- `ZM_hnfall_i(A, withU = false, remove = 1)` and `ZM_hnfall(A, withU = false, remove = 1)`:
  `[bigint[][], bigint[][] | null]`. `remove` is `0 | 1 | 2`: retain zero columns,
  remove them only from H, or remove the corresponding columns from H and U.
- `hnfall(A): [bigint[][], bigint[][]]`: remove zero columns only from H and return
  the full transformation U. The initial columns of A*U are zero; the rest equal H.

All inputs must have equally sized columns; inputs are preserved and outputs own
independent storage. Removing columns from U keeps only its image basis, which need
not be square. No transformation is computed unless requested.

```ts
import { ZM_hnf, ZM_hnfall, ZM_hnfall_i, hnfall } from '@sagemath-ts/parigp-ts/src/hnf_snf.js';
const hermiteInput = [[1n, 2n], [2n, 4n], [1n, 0n]];
ZM_hnf(hermiteInput); // [[1n, 0n], [0n, 2n]]
hnfall(hermiteInput);
// [[[1n, 0n], [0n, 2n]], [[-2n, 1n, 0n], [0n, 0n, 1n], [1n, 0n, -1n]]]
ZM_hnfall_i(hermiteInput, true, 2);
// [[[1n, 0n], [0n, 2n]], [[0n, 0n, 1n], [1n, 0n, -1n]]]
ZM_hnfall(hermiteInput, false, 0)[0]; // [[0n, 0n], [1n, 0n], [0n, 2n]]
```


### PARI integer polynomial factorization and Newton sums

`@sagemath-ts/parigp-ts/src/QX_factor.js` provides these typed entries:

- `ZX_factor(P): [bigint[], number][]` and `QX_factor([P, denominator])` return
  ordered primitive factors and multiplicities. Constant content is omitted;
  a nonzero constant gives `[]`, while the zero polynomial gives `[[[], 1]]`.
- `ZX_is_irred(P): boolean`; `ZX_squff(P): [bigint[][], number[]]` is the separate
  squarefree decomposition of a normalized nonzero integer polynomial.
- `ZX_gcd_all(A, B): [bigint[], bigint[]]` returns GCD and native scaled primitive
  quotient. Its second result is generally different from the raw quotient A/GCD.
- `ZX_DDF(P)` deflates before factorization and inflates factors in native order;
  `ZX_DDF_max(P, dmax)` uses the native degree-bound prime-selection route.
  DDF inputs must be squarefree, nonconstant, with nonzero constant and positive
  leading coefficient. These intermediate factors retain their unsorted order.
- `pick_prime(P, flag)` returns the chosen prime or `0n` for its modular stopping
  decision. Flag 0 counts factors, 1 counts roots, and other flags use the native
  degree test (including its negative-flag behavior).
- `combine_factors(P, famod, p, klim)` lifts and recombines monic modular factors.
  They must multiply to P/lc(P) modulo the prime p.
- `LLL_cmbf(P, famod, p, pa, bound, a, rec)` takes the unresolved lifted factor
  group produced by native bounded recombination. Supply its actual exponent a,
  modulus pa=p^a and maximum-subset rec; it is not a standalone arbitrary-lattice API.
- `chk_factors(P, M, bound, famod, pa)` checks the projected knapsack columns,
  returning factor arrays or null. M has one row per modular factor and at least
  two retained factor columns.
- `chk_factors_get(lt, famod, c, null, N)` selects nonzero entries of c and returns
  the modular product with optional leading scale lt. The empty product is `1n`.

All polynomial coefficients are ascending. Inputs are preserved; factorization
consumes the shared PARI random stream. Factor/knapsack helpers retain the native
caller preconditions. The extension-coefficient T argument is null in these adapters.

```ts
import { ZX_factor, QX_factor, ZX_is_irred, ZX_squff, ZX_gcd_all } from '@sagemath-ts/parigp-ts/src/QX_factor.js';
ZX_factor([-1n, 0n, 0n, 0n, 1n]); // [[[-1n, 1n], 1], [[1n, 1n], 1], [[1n, 0n, 1n], 1]]
QX_factor([[-7n, 0n, 7n], 3n]); // [[[-1n, 1n], 1], [[1n, 1n], 1]]
ZX_is_irred([-2n, 0n, 1n]); // true
ZX_squff([-1n, 3n, -3n, 1n]); // [[[-1n, 1n]], [3]]
ZX_gcd_all([1n, 3n, 2n], [-1n, -1n, 2n]); // [[1n, 2n], [2n, 2n]]
```

```ts
import { pick_prime, ZX_DDF, ZX_DDF_max, combine_factors, chk_factors, chk_factors_get } from '@sagemath-ts/parigp-ts/src/QX_factor.js';
const fourth = [-1n, 0n, 0n, 0n, 1n];
pick_prime(fourth, 0); // 3n
ZX_DDF_max(fourth, 0); // [[1n, 1n], [-1n, 1n], [1n, 0n, 1n]]
ZX_DDF(fourth); // [[1n, 0n, 1n], [1n, 1n], [-1n, 1n]]
combine_factors([-1n, 0n, 1n], [[1n, 1n], [2n, 1n]], 3n, 1); // [[1n, 1n], [-1n, 1n]]
chk_factors([-1n, 0n, 1n], [[1n, 0n], [0n, 1n]], 3n, [[1n, 1n], [124n, 1n]], 125n);
// [[-1n, 1n], [1n, 1n]]
chk_factors_get(null, [[1n, 1n], [4n, 1n]], [1n, 0n], null, 5n); // [1n, 1n]
```

`@sagemath-ts/parigp-ts/src/polarit2.js` exports `polsym(P, n)` and
`polsym_gen(P, prefix, n, null, N)`. P has integer coefficients. N=null returns
normalized rational pairs; a positive modulus returns bigints with the original
nonpositive residue storage. Results include the zeroth sum (degree) through n.
A null prefix starts from the degree; an existing prefix is copied and continued,
without recomputing or verifying its entries. It must fit within the requested result.

```ts
import { polsym, polsym_gen } from '@sagemath-ts/parigp-ts/src/polarit2.js';
polsym([1n, 0n, 2n], 4); // [[2n, 1n], [0n, 1n], [-1n, 1n], [0n, 1n], [1n, 2n]]
polsym_gen([-2n, 0n, 1n], [2n, 0n], 4, null, 5n); // [2n, 0n, -1n, 0n, -2n]
```

### Number-field constructor factorization window

`RationalPolynomial.isIrreducible()` and checked `NumberField` construction use
Sage's PARI factorization route for degrees 30 through 300, inclusively. Boolean
irreducibility results are cached on each immutable polynomial, including false.
The NTL ranges outside this window still use the legacy local backend and retain
its documented bounds and performance gaps.

```ts
import { NumberField, RationalPolynomial } from 'sagemath-ts/rings';
const largeEquation = RationalPolynomial.fromBigInts([-2n, ...Array<bigint>(199).fill(0n), 1n]);
largeEquation.isIrreducible(); // true
new NumberField(largeEquation, 'a').degree(); // 200
```


### NTL integer-factor trace dependencies

Import directly from `@sagemath-ts/ntl-ts/src/ZZXFactoring.js`, with
`@sagemath-ts/ntl-ts` installed in the importing package. In this workspace it is
a dependency of `packages/sagemath-ts`. These are native
factorization stages used by the complete integer driver described below.
Vectors use ordinary zero-indexed arrays and matrices use rows. All inputs are
copied for returned updates. Native long parameters are safe integer `number`s;
the prime `p` in `Compute_pb`, `Compute_pdelta` and `Compute_pb_eff` also accepts
`bigint` to preserve the full native word range;
coefficients, traces, moduli and bounds are `bigint`s. Supply valid dimensions,
positive prime powers, a prime `p >= 2`, and reduced trace residues for chopping.

| Entry | Result and input contract |
|---|---|
| `ComputeTrace(Tr, f, d, P)` | Copied `bigint[]` with entry `d-1` replaced; `f` monic of positive degree, `d > 0`, `Tr.length >= d`, `P > 1`. Earlier traces are supplied by the caller. |
| `ChopTraces(C, Tr, d, pb, pdelta, P, lc)` | Copied `bigint[]`; first `d` entries replaced, unused tail retained. `C`, `Tr`, `pb` have at least `d` entries. |
| `DenseChopTraces(C, Tr, d, d1, pb_eff, pdelta, P, lc, A)` | Copied `bigint[]`; first `d1` entries replaced using an `A` with `d1` rows and `d` columns. |
| `Compute_pb(b, pb, p, d, root_bound, n)` | `[number[], bigint[]]`; preserve the preceding `d-1` entries and resize both returned vectors to `d`. |
| `Compute_pdelta(delta, pdelta, p, bit_delta)` | `[number, bigint]`; advance until the modulus has more than `bit_delta` bits. |
| `BuildReductionMatrix(r, d, pdelta, chop_vec, B_L)` | `[bigint[][], number]`, the matrix and native scale `C`. `chop_vec` has `r` rows of `d` traces; `B_L` has `r` columns. |
| `Compute_pb_eff(p, d, root_bound, n, ran_bits)` | `[number, bigint]`, the compressed-trace exponent and prime power. |
| `d1_val(bit_delta, r, s)` | `number`, native compressed-trace count; `bit_delta > 0`. |

```ts
import { ComputeTrace, ChopTraces, DenseChopTraces }
  from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
ComputeTrace([0n, 0n, 91n], [-2n, 0n, 1n], 2, 125n); // [0n, 4n, 91n]
ChopTraces([0n, 0n, 99n], [4n, 8n], 2, [5n, 5n], 5n, 125n, 2n);
// [2n, 1n, 99n]
DenseChopTraces([0n, 0n, 99n], [4n, 8n], 2, 2, 5n, 5n, 125n, 2n,
  [[1n, 2n], [-1n, 1n]]); // [-1n, 0n, 99n]
```

```ts
import { Compute_pb, Compute_pdelta, Compute_pb_eff, d1_val, BuildReductionMatrix }
  from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
Compute_pb([], [], 3, 1, 2n, 4); // [[3], [27n]]
Compute_pdelta(1, 3n, 3, 4); // [3, 27n]
Compute_pb_eff(3, 2, 1n, 4, 1); // [4, 81n]
d1_val(8, 16, 8); // 5
BuildReductionMatrix(2, 2, 9n, [[1n, 2n], [3n, 4n]], [[1n, 0n], [1n, -1n]]);
// [[[2n,0n,1n,2n], [2n,-2n,-2n,-2n], [0n,0n,9n,0n], [0n,0n,0n,9n]], 2]
```


### NTL integer products and multifactor lifting

These direct module imports require `@sagemath-ts/ntl-ts` in the importing
package's dependencies (the workspace `packages/sagemath-ts` supplies it).
Coefficients are constant-first `bigint` arrays. Inputs remain independent of
returned arrays. The complete integer factorization driver uses these helpers.
`ZZX1.mul(a, b, state?: PolynomialProductState)` and
`ZZX1.sqr(a, state?: PolynomialProductState)` accept the shared FFT cache/stream
pair. `ChooseSS(da: number, maxbitsa: number, db: number, maxbitsb: number): boolean`
exposes NTL's transform-choice predicate; the counts are exactly represented native
integers, with the same numeric guards as SSRatio plus a coefficient-size overflow
guard. Omitting state retains the previous coefficient-only product behavior.


```ts
import { mul, sqr } from '@sagemath-ts/ntl-ts/src/ZZX1.js';
mul([-2n, 0n, 3n], [5n, -7n, 1n]); // [-10n, 14n, 13n, -21n, 3n]
sqr([-2n, 0n, 3n]); // [4n, 0n, -12n, 0n, 9n]
```

`MultiLift(a: readonly (readonly bigint[])[], f: readonly bigint[], e: number,
p: bigint, options?: zz_pXOptions): bigint[][]` lifts at least two monic factors modulo the word prime
`p` to modulus `p ** BigInt(e)`. The target `f` must be monic. For `e > 1`, the
factors must be pairwise coprime and their product congruent to `f` modulo `p`.
The exponent must be an exactly represented integer; the native overflow check
rejects `e >= 2 ** 60`. The exponent-one path copies the normalized factors
before checking coprimality or product congruence. During lifting, constant-factor
quotients are rejected by native modulus construction. Input factor order is retained.
Optional `zz_pXOptions` supplies the word `maxroot` and shared product `state`.

```ts
import { MultiLift } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
MultiLift([[1n, 1n], [2n, 1n]], [7n, 3n, 1n], 1, 5n);
// [[1n, 1n], [2n, 1n]]
MultiLift([[1n, 1n], [2n, 1n]], [7n, 3n, 1n], 4, 5n);
// [[281n, 1n], [347n, 1n]]
```


### NTL exact integer lattices

Direct imports require `@sagemath-ts/ntl-ts` in the importing package. Matrices
are rectangular row arrays of `bigint`; returned arrays are independent of inputs.
Native integer parameters use integer `number` values.

| Entry from `@sagemath-ts/ntl-ts/src/LLL.js` | Returned tuple |
|---|---|
| `LLL(B, a=3, b=4, transform=false)` | `[rank, det, basis, U]` |
| `LLL_plus(B, a=3, b=4, transform=false)` | `[rank, D, basis, U]` |
| `image(B, transform=false)` | `[rank, det, basis, U]` |
| `LatticeSolve(A, y, reduce=0, initial=[], columns=A[0]?.length ?? 0)` | `[success, x]`; success is 0 or 1. Failure returns a copy of `initial`. |

`D` is a `bigint[]` of ordered Gram determinants beginning with 1; `det` is its
final value (1 for rank zero). Basis matrices retain all rows, with dependent
zero rows first. A requested transformation satisfies `U * B = basis`; otherwise
`U` is `null`. LLL requires `a > 0`, `b > 0`, `a <= b` and `a > floor(b/4)`.
LatticeSolve uses `reduce=0`, `1` or `2`; its result satisfies `x * A = y` on
success. Supply `columns` to represent a matrix with zero rows and nonzero width.

```ts
import { LLL_plus, image, LatticeSolve } from '@sagemath-ts/ntl-ts/src/LLL.js';
LLL_plus([[1n,1n,1n], [-1n,0n,2n], [3n,5n,6n]], 3, 4, true);
// [3, [1n,1n,2n,9n], [[0n,1n,0n],[1n,0n,1n],[-1n,0n,2n]],
//  [[-4n,-1n,1n],[5n,1n,-1n],[0n,1n,0n]]]
image([[2n,0n], [0n,3n], [1n,1n]])[2]; // [[0n,0n],[1n,0n],[0n,1n]]
LatticeSolve([[2n,0n], [0n,3n], [1n,1n]], [1n,2n]); // [1, [1n,1n,-1n]]
LatticeSolve([[2n,0n], [0n,2n]], [1n,0n], 0, [19n]); // [0, [19n]]
```

```ts
import { _ntl_gexteucl } from '@sagemath-ts/ntl-ts/src/lip.js';
_ntl_gexteucl(0n, 0n); // [1n, 0n, 0n]
_ntl_gexteucl(3n, 2n); // [1n, -1n, 1n], returned as [s,t,d]
```


### NTL factor-recovery stages

Import these from `@sagemath-ts/ntl-ts/src/ZZXFactoring.js` in a package that
declares the NTL dependency. Inputs remain independent of returned arrays.

| Entry | Result and contract |
|---|---|
| `PolyEval(f: readonly bigint[], a: bigint)` | `bigint`; native Horner evaluation, rejecting the zero polynomial. |
| `RootBound(f: readonly bigint[])` | `bigint`; the native scaled root bound, requiring a nonzero constant coefficient. |
| `CutAway(D, M, C, r, d)` | `bigint[][]`; remove rows by Gram-determinant bounds and return the integer image of the first `r` columns after division by `C`. `M` uses rows; `D` has at least `M.length + 1` entries. Integer `C/r/d` use `number`. |
| `AdditionalLifting(P1, e1, w1, p, new_bound, f, doubling, state?: PolynomialProductState)` | `[bigint, number, bigint[][]]`, the new modulus, exponent and factors. `p/P1` and coefficients are `bigint`; `e1/new_bound` are integer numbers. `P1` is ignored as an input. |

Additional lifting uses `max(2*e1, new_bound)` when `doubling` is true and
`new_bound` otherwise. It normalizes the target modulo the new prime power and
reduces the old factor list modulo the word prime before relifting.

```ts
import { PolyEval, RootBound, CutAway }
  from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
PolyEval([-2n, 0n, 1n], 3n); // 7n
RootBound([-2n, 0n, 1n]); // 2n
RootBound([-2n, 3n]); // 3n
RootBound([-7n]); // 7n
CutAway([1n,1n,1n,100n], [[1n,0n,0n],[0n,1n,0n],[0n,0n,10n]], 1, 2, 1);
// [[1n,0n],[0n,1n]]
```

```ts
import { AdditionalLifting } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
AdditionalLifting(5n, 1, [[1n,1n],[2n,1n]], 5n, 4, [7n,3n,1n], false);
// [625n, 4, [[281n,1n],[347n,1n]]]
```


### NTL factor selection helpers

These entries are exported from `@sagemath-ts/ntl-ts/src/ZZXFactoring.js`.
Coefficient vectors contain `bigint`, degree/count/index vectors contain integer
`number`, and modular helpers take an explicit `bigint` modulus `p > 1`.
Returned vectors are independent of inputs. Factor indices must be valid,
sorted and distinct for complement products and removal.

| Entry | Result |
|---|---|
| `inplace_rev(f)` | Normalized reversed integer coefficient vector. |
| `RecordPattern(pat, fac, p)` | New count vector of `pat.length`; `fac` consists of `[polynomial, positive_degree]` pairs, with valid degree indices. Repeated degrees overwrite. Requires a word-prime `p`. |
| `NumFactors(pat)` | Sum of every count, including index zero. |
| `CalcPossibleDegrees(pat)` | `bigint` encoding all possible subset degrees as set bits. |
| `CalcPossibleDegrees(fac, k, p)` | `bigint[]` of native suffix bitsets for cardinality `k`, retaining trailing lower-cardinality slots. Empty factors return `[]` before validating `k`. |
| `ConstTermTest(W, I, ct, lc, prod, ProdLen, p)` | `[number, bigint[], number]`: divisibility result (0/1), copied prefix-product cache and updated prefix length. `I` must be nonempty and `prod` sufficiently long. |
| `BalCopy(G, p)` | Centered coefficient vector, retaining a positive coefficient at the even-modulus midpoint. |
| `mul(W, p, I?)` | Modular product of all factors or the selected indices; empty product is `[1n]`. |
| `InvMul(W, I, p, state?: PolynomialProductState)` | Modular product of the complementary factors. |
| `RemoveFactors(W, I, p)` | Remaining modular factors in their original order. |
| `unpack(a, n)` | Bits 0 through `n` of `abs(a)`, as numbers; `n = -1` returns `[]`. |
| `SubPattern(p1, p2)` | Elementwise difference; errors on unequal lengths or a negative result. |

```ts
import { CalcPossibleDegrees, NumFactors, unpack, SubPattern }
  from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
CalcPossibleDegrees([0, 2, 1]); // 31n (degrees 0 through 4)
NumFactors([0, 2, 1]); // 3
unpack(-13n, 4); // [1, 0, 1, 1, 0]
SubPattern([0, 2, 1], [0, 1, 0]); // [0, 1, 1]
```

```ts
import { mul, InvMul, ConstTermTest, BalCopy }
  from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
const W = [[1n, 1n], [2n, 1n], [3n, 1n]];
mul(W, 101n); // [6n, 11n, 6n, 1n]
InvMul(W, [1], 101n); // [3n, 4n, 1n]
ConstTermTest(W, [1, 2], 12n, 1n, [7n, 8n, 9n], 0, 101n);
// [1, [2n, 6n, 9n], 1]
BalCopy([8n, 9n, 16n], 16n); // [8n, -7n]
```


### NTL word polynomial quotient arithmetic

Import from `@sagemath-ts/ntl-ts/src/lzz_pX.js` in a package declaring the NTL
dependency. Polynomial arrays contain `bigint`; `p` is an explicit word modulus
with `1 < p < 2^60`. Arrays are normalized and copied. Use `build` when replacing
a modulus; its arithmetic, reciprocal and crossover fields are internal.

| Entry | Result and contract |
|---|---|
| `new zz_pXModulus(f: readonly bigint[] \| null, p: bigint, options?: zz_pXOptions)` | Context `F` with polynomial `f`, degree `n` and modulus `p`. `null` leaves it uninitialized. A polynomial must have positive degree. |
| `build(F, f)` | Rebuilds `F` and its caches; returns `void`. |
| `rem(a, F)` | `bigint[]` remainder; `a` may have arbitrary degree. |
| `MulMod(a, b, F)` | Reduced product; both operands must have degree less than `F.n`. |
| `SqrMod(a, F)` | Reduced square; `a` must have degree less than `F.n`. |
| `mul(a, b, p, options?: zz_pXOptions)` | Ordinary polynomial product as `bigint[]`; follows the native word FFT crossover and size checks. |
| `sqr(a, p, options?: zz_pXOptions)` | Ordinary polynomial square as `bigint[]`. |
| `MulByXMod(a, f, p)` | Reduced X times `a`; `a` must have degree less than `f`. |
| `InvMod(a, f, p, options?: zz_pXOptions)` | Inverse, or an error if noninvertible; requires reduced `a`. |
| `InvModStatus(a, f, p, options?: zz_pXOptions)` | `[number, bigint[]]`: status 0 and an inverse, or status 1 and the monic gcd. |
| `PowerXMod(e: bigint, F)` | X to the exponent `e`. |
| `PowerXPlusAMod(a: bigint, e: bigint, F)` | `(X+a)^e`. |
| `PowerMod(g, e: bigint, F)` | `g^e`; `g` must be reduced, including for exponent zero. |

Negative powers invert their final result with the options captured by `F`.
Inverse functions use the native word extended GCD schedule. Composite moduli are permitted;
required scalar or polynomial inverses may not exist. Construction preserves
NTL's default 60-bit context profile, including when those errors occur.

```ts
import { zz_pXModulus, PowerXMod, PowerXPlusAMod, PowerMod, InvModStatus }
  from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
const F = new zz_pXModulus([1n, 0n, 1n], 5n);
PowerXMod(3n, F); // [0n, 4n]
PowerXMod(-1n, F); // [0n, 4n]
PowerXPlusAMod(1n, 3n, F); // [3n, 2n]
PowerMod([1n, 1n], 3n, F); // [3n, 2n]
InvModStatus([2n, 1n], F.f, 5n); // [1, [2n, 1n]]
```

```ts
import { zz_pXModulus, rem, MulMod, SqrMod, build }
  from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
const F = new zz_pXModulus([1n, 0n, 1n], 5n);
rem([1n, 2n, 3n, 4n], F); // [3n, 3n]
MulMod([1n, 1n], [2n, 1n], F); // [1n, 3n]
SqrMod([1n, 1n], F); // [0n, 2n]
build(F, [-1n, 1n]);
rem([1n, 2n, 3n], F); // [1n]
```


### NTL word matrix multiplication

Import `mul` from `@sagemath-ts/ntl-ts/src/mat_lzz_p.js` in a package declaring
the NTL dependency:

```ts
mul(A: readonly (readonly bigint[])[], B: readonly (readonly bigint[])[],
    p: bigint, columns?: number, inner?: number): bigint[][]
```

Matrices use rows. `p` must satisfy `1 < p < 2^60`; composite moduli are
supported. Inputs must be rectangular and compatible. `columns` defaults to
`B[0]?.length ?? 0`; `inner` defaults to `A[0]?.length ?? B.length`. Supply them
when empty row lists need explicit dimensions. Results do not alias inputs.

```ts
import { mul } from '@sagemath-ts/ntl-ts/src/mat_lzz_p.js';
const A = [[1n, 2n], [3n, 4n]];
mul(A, [[5n, 6n], [7n, 8n]], 101n); // [[19n, 22n], [43n, 50n]]
mul(A, A, 5n); // [[2n, 0n], [0n, 2n]]
```

```ts
import { mul } from '@sagemath-ts/ntl-ts/src/mat_lzz_p.js';
mul([[], [], []], [], 17n, 4, 0);
// [[0n,0n,0n,0n], [0n,0n,0n,0n], [0n,0n,0n,0n]]
```


### NTL word polynomial composition and GCD

Import from `@sagemath-ts/ntl-ts/src/lzz_pX1.js` in a package declaring the NTL
dependency. `F` is a `zz_pXModulus` from `lzz_pX.js`. Coefficients use bigint;
integer size/iteration parameters use exactly represented numbers.

| Entry | Result and contract |
|---|---|
| `new zz_pXNewArgument()` | Empty composition cache with `mat: bigint[][]` and `poly: bigint[]`. |
| `build(H, h, F, m)` | Fill cache `H` for base polynomial `h`, with positive integer block size `m`; `h` must be reduced modulo `F`. Returns void. |
| `CompMod(g, h, F)` | `g(h)` reduced modulo `F`, using either a coefficient vector `h` or a prebuilt `zz_pXNewArgument`. |
| `Comp2Mod(g1, g2, h, F)` | `[bigint[], bigint[]]`, two compositions sharing a cache. |
| `Comp3Mod(g1, g2, g3, h, F)` | `[bigint[], bigint[], bigint[]]`, three compositions sharing a cache. |
| `reduce(H, F)` | Rebuild `H` for the reduced base modulo the new `F`; retains its block size. Returns void. |
| `GCD(u, v, p: bigint, options?: zz_pXOptions)` | Monic polynomial gcd as bigint coefficients, or `[]` when both are zero. |
| `XGCD(a, b, p: bigint, options?: zz_pXOptions)` | `[d, s, t]` with monic `d = a*s + b*t`; returns `[[], [1n], []]` for two zero polynomials. |

Caches and returned arrays are independent of the supplied polynomial arrays.
Use build/reduce to update cache state consistently with its coefficient context.
A zero outer polynomial in the direct CompMod form returns zero before building
the cache; the prebuilt form also returns constants before checking cache state.

```ts
import { zz_pXModulus } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
import { zz_pXNewArgument, build, CompMod, Comp2Mod, Comp3Mod }
  from '@sagemath-ts/ntl-ts/src/lzz_pX1.js';
const F = new zz_pXModulus([1n,0n,1n], 5n);
const H = new zz_pXNewArgument();
build(H, [1n,1n], F, 2);
H.mat; // [[1n,0n], [1n,1n]]
H.poly; // [0n,2n]
CompMod([1n,2n,3n], H, F); // [3n,3n]
Comp2Mod([1n,2n,3n], [0n,1n], [1n,1n], F); // [[3n,3n], [1n,1n]]
Comp3Mod([1n,2n,3n], [0n,1n], [1n], [1n,1n], F);
// [[3n,3n], [1n,1n], [1n]]
```

```ts
import { zz_pXModulus } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
import { zz_pXNewArgument, build, reduce, CompMod, GCD }
  from '@sagemath-ts/ntl-ts/src/lzz_pX1.js';
const F = new zz_pXModulus([1n,0n,1n], 5n), H = new zz_pXNewArgument();
build(H, [1n,1n], F, 2);
const G = new zz_pXModulus([-1n,1n], 5n);
reduce(H, G);
H.mat; // [[1n], [2n]]
H.poly; // [4n]
CompMod([1n,2n,3n], H, G); // [2n]
GCD([-2n,0n,2n], [-3n,3n], 5n); // [4n,1n]
```

`TraceMap(a, d, F, b)` and `PowerCompose(h, q, F)` are exported from
`@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js`. They implement NTL's binary
composition trace and iterated composition, returning bigint coefficient arrays.
Negative integer d/q values error. PowerCompose with q=0 returns the polynomial
X for an initialized modulus; an uninitialized modulus raises NTL's scratch
vector-length error before this shortcut.

```ts
import { zz_pXModulus } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
import { TraceMap, PowerCompose }
  from '@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js';
const F = new zz_pXModulus([1n,0n,1n], 3n), frobenius = [0n,2n];
TraceMap([1n,1n], 2, F, frobenius); // [2n]
PowerCompose(frobenius, 2, F); // [0n,1n]
```


### NTL distinct-degree factorization

Import `NewDDF` and `SFCanZass1` from
`@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js` in a package declaring the NTL
dependency. Coefficients and the explicit prime modulus `p` use bigint.

- `NewDDF(f, h, p, options?: zz_pXOptions): [bigint[], number][]` groups the monic squarefree input by
  irreducible-factor degree, in native order. Supply `h = X^p mod f`.
- `SFCanZass1(f, p, options?: zz_pXOptions): [[bigint[], number][], bigint[]]` computes that Frobenius
  polynomial and returns `[groups, h]`. It rejects constant inputs.

Returned arrays are independent of the inputs. Each group's polynomial may
still contain several irreducible factors of the reported degree.

```ts
import { SFCanZass1, NewDDF }
  from '@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js';
const f = [-2n,2n,-1n,1n]; // (X-1)(X^2+2), over F_5
const [groups, h] = SFCanZass1(f, 5n);
groups; // [[[4n,1n],1], [[2n,0n,1n],2]]
h; // [3n,4n,4n]
NewDDF(f, h, 5n); // [[[4n,1n],1], [[2n,0n,1n],2]]
```

```ts
import { NewDDF } from '@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js';
NewDDF([1n], [], 5n); // []
NewDDF([-1n,1n], [1n], 5n); // [[[4n,1n],1]]
```


### NTL sequence minimum polynomials

`MinPolySeq(a: readonly bigint[], m: number, p: bigint, options?: zz_pXOptions): bigint[]` is exported
from `@sagemath-ts/ntl-ts/src/lzz_pX1.js`. It reconstructs the minimum polynomial
of a linearly generated sequence over the prime field, with caller-supplied degree
bound `m`. Supply at least `2*m` entries; later entries are ignored. The native
Berlekamp–Massey/half-GCD crossover is preserved. The result is an independent
coefficient array. Invalid degree bounds or short sequences error.

```ts
import { MinPolySeq } from '@sagemath-ts/ntl-ts/src/lzz_pX1.js';
const fibonacci = [0n,1n,1n,2n,3n,5n];
MinPolySeq(fibonacci, 3, 101n); // [100n,100n,1n], X^2-X-1
MinPolySeq([...fibonacci,999n,-1n], 3, 101n); // [100n,100n,1n]
```

```ts
import { MinPolySeq } from '@sagemath-ts/ntl-ts/src/lzz_pX1.js';
MinPolySeq([0n,0n,0n,0n], 2, 101n); // [1n]
MinPolySeq([1n,0n,0n,0n], 2, 101n); // [0n,1n]
MinPolySeq([], 0, 101n); // [1n]
```


### NTL small-prime iterator

Import `PrimeSeq` from `@sagemath-ts/ntl-ts/src/ZZ.js` in a package declaring the
NTL dependency. `new PrimeSeq()` starts at two. `next(): bigint` returns the next
prime, or zero after the native bound 1073676289. `reset(b: bigint): void` moves
to the first prime at least b; values at most two restart the sequence.
Each instance has independent iterator state.

```ts
import { PrimeSeq } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const primes = new PrimeSeq();
Array.from({length:6}, () => primes.next()); // [2n,3n,5n,7n,11n,13n]
primes.reset(100n);
primes.next(); // 101n
primes.next(); // 103n
```

```ts
import { PrimeSeq } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const primes = new PrimeSeq();
primes.reset(1073676289n);
primes.next(); // 0n
primes.next(); // 0n
primes.reset(-5n);
primes.next(); // 2n
```


### NTL cached multipliers and transposed projections

Import `zz_pXMultiplier`, `build` and `MulMod` from
`@sagemath-ts/ntl-ts/src/lzz_pX.js`. `new zz_pXMultiplier()` starts with the zero
polynomial; `new zz_pXMultiplier(b, F)` builds a cache for reduced b.
`build(B, b, F)` rebuilds it, and `MulMod(a, B, F)` uses it. Existing modulus-build
and polynomial-product signatures are still available. `B.val()` returns an
independent coefficient array. Use build to keep cache state consistent with F.

`UpdateMap(a, B, F): bigint[]` and `ProjectPowers(a, count, h, F): bigint[]` are
exported from `lzz_pX1.js`. UpdateMap returns the inner products of a with
`b*X^i mod F` for i=0..deg(F)-1, stripping trailing zeros. ProjectPowers returns
exactly count inner products with `1, h, h^2, ... mod F`; h can be a polynomial
or a `zz_pXNewArgument`. ProjectPowers requires a.length <= deg(F); UpdateMap
strips trailing zeros before checking that limit.

```ts
import { zz_pXModulus, zz_pXMultiplier, MulMod }
  from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
import { UpdateMap } from '@sagemath-ts/ntl-ts/src/lzz_pX1.js';
const F = new zz_pXModulus([1n,0n,1n], 5n);
const B = new zz_pXMultiplier([1n,1n], F);
B.val(); // [1n,1n]
B.UseFFT; // 0, the native plain-product path for this small degree
MulMod([1n,2n], B, F); // [4n,3n]
UpdateMap([1n,2n], B, F); // [3n,1n]
```

```ts
import { zz_pXModulus } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
import { zz_pXNewArgument, build, ProjectPowers }
  from '@sagemath-ts/ntl-ts/src/lzz_pX1.js';
const F = new zz_pXModulus([1n,0n,1n], 5n), a = [1n,2n];
ProjectPowers(a, 5, [0n,1n], F); // [1n,2n,4n,3n,1n]
const H = new zz_pXNewArgument();
build(H, [0n,1n], F, 2);
ProjectPowers(a, 5, H, F); // [1n,2n,4n,3n,1n]
ProjectPowers(a, 0, [0n,1n], F); // []
```

The prepared form requires a positive count. A valid prepared cache with count=0
throws `ProjectPowers: prepared argument requires a positive count`; bundled NTL
crashes on this undefined input. Invalid vector lengths and uninitialized caches
retain their earlier errors. See DEVIATIONS.md for the recorded native probe.


### NTL deterministic byte streams

Import from `@sagemath-ts/ntl-ts/src/ZZ.js` in a package declaring that dependency:

- `sha256(data: Uint8Array, length = 32): Uint8Array` and
  `hmac_sha256(key: Uint8Array, data: Uint8Array, length = 32): Uint8Array` return
  the digest prefix, with length clamped to 0..32.
- `DeriveKey(data: Uint8Array, length: number): Uint8Array` uses NTL's HMAC counter
  construction. Negative lengths throw `DeriveKey: bad args`.
- `new RandomStream(key: Uint8Array | RandomStream)` uses the first 32 key bytes
  or copies the source stream's current state. `get(n: number): Uint8Array`
  returns independent bytes; negative counts throw `RandomStream::get: bad args`.
  `set_nonce(nonce: bigint): void` wraps modulo 2^64 and restarts the counter.
  `assign(other: RandomStream): this` copies state into the existing receiver;
  aliases see the change. Counts must be exactly represented integer numbers.
- The historical names `salsa20_init(key: Uint8Array): Uint32Array`,
  `salsa20_core(state: Uint32Array): void` and
  `salsa20_apply(state: Uint32Array): Uint32Array` expose NTL's **ChaCha20** kernel.
  Core mutates the first 16 words. Apply returns a new block and increments the
  input's 64-bit counter. Keys shorter than 32 bytes and states shorter than 16
  words throw adapter errors.

This is the generic NTL byte-stream profile. HMAC retains native 32-bit data-length
conversion. Global seeding and FFT-cache random consumption remain separate work.

```ts
import { sha256, RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const hex = (a: Uint8Array) => Array.from(a, x => x.toString(16).padStart(2, '0')).join('');
hex(sha256(new Uint8Array()));
// 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855'
hex(new RandomStream(new Uint8Array(32)).get(16));
// '76b8e0ada0f13d90405d6ae55386bd28'
```

```ts
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const stream = new RandomStream(new Uint8Array(32));
stream.get(17);
const copy = new RandomStream(stream);
Array.from(stream.get(5)); // [210, 25, 184, 160, 141]
Array.from(copy.get(5));   // [210, 25, 184, 160, 141]
stream.set_nonce(0n);
Array.from(stream.get(4)); // [118, 184, 224, 173]
```


### NTL integer sampling with an explicit stream

Import from `@sagemath-ts/ntl-ts/src/ZZ.js`. All of these helpers take a
`RandomStream`; they do not use a process-global RNG.

- `RandomWord(stream): bigint` draws an unsigned 64-bit word.
  `VectorRandomWord(k: number, stream): bigint[]` draws k such words (empty for k<=0).
- `RandomBits_long(l: number, stream): bigint` permits 0..63 bits;
  `RandomBits_ulong(l: number, stream): bigint` permits 0..64.
  Nonpositive lengths return zero; larger lengths throw `RandomBits: length too big`.
- `RandomLen_long(l: number, stream): bigint` produces an exact l-bit signed word;
  l>=64 throws `RandomLen: length too big`.
- `RandomBits(l: number, stream): bigint` and its alias `RandomBits_ZZ` produce
  arbitrary-size values below 2^l. `RandomLen(l: number, stream): bigint` and its
  alias `RandomLen_ZZ` produce exact l-bit values. Nonpositive lengths return zero;
  exact length one returns one without a draw. Lengths >=2^60 throw native overflow
  errors. Requests must fit available memory; counts are exactly represented integers.
- `RandomBnd(bound: bigint, stream, options?: {word?: boolean}): bigint` returns
  a value below bound, or zero for bound<=1. Its default selects the native ZZ
  overload. `word:true` selects the native long overload and requires a signed
  64-bit bound; an out-of-range bound throws a RangeError before sampling.

Native overloads consume different bytes even when they return the same integer:

```ts
import { RandomStream, RandomBnd } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const a = new RandomStream(new Uint8Array(32));
const b = new RandomStream(new Uint8Array(32));
RandomBnd(256n, a);               // 118n
RandomBnd(256n, b, {word:true});  // 118n
Array.from(a.get(4)); // [224, 173, 160, 241]
Array.from(b.get(4)); // [184, 224, 173, 160]
```

```ts
import { RandomStream, RandomLen_long, RandomLen_ZZ } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const a = new RandomStream(new Uint8Array(32));
const b = new RandomStream(new Uint8Array(32));
RandomLen_long(9, a); // 374n
RandomLen_ZZ(9, b);   // 374n
Array.from(a.get(4)); // [184, 224, 173, 160]
Array.from(b.get(4)); // [224, 173, 160, 241]
```


### NTL minimum polynomials of quotient elements

Import from `@sagemath-ts/ntl-ts/src/lzz_pX1.js`, using `zz_pXModulus` from
`@sagemath-ts/ntl-ts/src/lzz_pX.js` and `RandomStream` from its `ZZ.js` module:

- `DoMinPolyMod(g: readonly bigint[], F, m: number, R: readonly bigint[]): bigint[]`
  reconstructs from the sequence of R applied to successive powers of g.
- `IrredPolyMod(g: readonly bigint[], F, m=F.n): bigint[]` uses R=[1]; it requires
  that the element's minimum polynomial is irreducible for its native guarantee.
- `ProbMinPolyMod(g, F, stream): bigint[]` or `(g, F, m: number, stream): bigint[]`
  returns a probable divisor using random projection.
- `MinPolyMod(g, F, stream): bigint[]` or `(g, F, m: number, stream): bigint[]`
  verifies and retries to return the full minimum polynomial.

Use a prime coefficient modulus, reduced g, and a valid degree bound m (default
F.n). The bound must satisfy 1<=m<=F.n and bound the true degree for the native
minimum-polynomial guarantees. Arrays list coefficients from constant upward.
The randomized helpers consume the supplied stream and preserve native error
ordering, including errors after random draws.

```ts
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { zz_pXModulus } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
import { ProbMinPolyMod, MinPolyMod } from '@sagemath-ts/ntl-ts/src/lzz_pX1.js';
const F = new zz_pXModulus([1n,1n,1n], 2n);
const a = new RandomStream(new Uint8Array(32));
const b = new RandomStream(new Uint8Array(32));
ProbMinPolyMod([0n,1n], F, a); // [1n]
MinPolyMod([0n,1n], F, b);     // [1n,1n,1n]
Array.from(a.get(4)); // [224,173,160,241]
Array.from(b.get(4)); // [160,241,61,144]
```

```ts
import { zz_pXModulus } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
import { DoMinPolyMod, IrredPolyMod } from '@sagemath-ts/ntl-ts/src/lzz_pX1.js';
const F = new zz_pXModulus([1n,1n,1n], 2n);
DoMinPolyMod([0n,1n], F, 2, [1n]); // [1n,1n,1n]
IrredPolyMod([0n,1n], F);          // [1n,1n,1n]
```


### NTL cached samplers, random polynomials and root products

From `@sagemath-ts/ntl-ts/src/ZZ.js`:

- `new RandomBndGenerator(bound: bigint|null, stream: RandomStream)` caches a
  signed-word bound greater than one. Null constructs an uninitialized generator.
  `build(bound: bigint): void` updates the cache; an invalid bound leaves it intact.
  `next(): bigint` samples below that bound. An uninitialized next call throws
  `RandomBndGenerator::next: uninitialized generator`.
- `new RandomBndGenerator(other)` copies the cache and shares its stream reference.
  `assign(other): this` copies the bound and stream into the existing receiver.
  Public `p: bigint`, `nb: number` and `mask: bigint` expose the initialized cache.
- `VectorRandomBnd(k: number, bound: bigint, stream): bigint[]` returns k samples,
  an empty array for k<=0, or k zeros for bound<=1. Bounds must fit signed 64 bits.
  Counts must be exactly represented integers and fit available memory.

From `@sagemath-ts/ntl-ts/src/lzz_pX.js`,
`random(n: number, p: bigint, stream: RandomStream): bigint[]` draws n coefficients
modulo p and strips trailing zero coefficients. The result can have degree below
n-1, including the zero polynomial `[]`. Negative n throws a native length error.

From `@sagemath-ts/ntl-ts/src/lzz_pX1.js`,
`BuildFromRoots(roots: readonly bigint[], p: bigint, options?: zz_pXOptions): bigint[]` returns the monic
product of X-root, preserving multiplicities and zero roots. It returns `[1n]`
for an empty root list and uses a balanced product tree above the native crossover.

```ts
import { BuildFromRoots } from '@sagemath-ts/ntl-ts/src/lzz_pX1.js';
BuildFromRoots([1n,2n,0n], 5n); // [0n,2n,2n,1n]
```

```ts
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { random } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
const stream = new RandomStream(new Uint8Array(32));
random(3, 5n, stream);      // []
Array.from(stream.get(4)); // [241,61,144,64]
```


### NTL prime-field roots and square-free factor recovery

Import from `@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js`. Polynomials are dense
bigint coefficient arrays from constant upward. Supply a prime bigint modulus p
and an explicit `RandomStream` from the dependency's `ZZ.js` module.

- `FindRoots(f, p, stream, options?: zz_pXOptions): bigint[]` returns roots in native split order.
  `FindRoot(f, p, stream, options?: zz_pXOptions): bigint` finds one root using the smaller successful split.
  These require a monic polynomial that splits into distinct linear factors.
- `RootEDF(f, p, stream, options?: zz_pXOptions): bigint[][]` returns the corresponding linear factors.
- `FindFactors(f, g, roots, p, options?: zz_pXOptions): bigint[][]` recovers ordered groups for an element g
  whose distinct field values are supplied in roots. It is deterministic.
- `EDFSplit(f, b, d: number, p, stream, options?: zz_pXOptions): bigint[][]` performs one equal-degree split.
  It may return groups that require further splitting.
- `EDF(f, b, d: number, p, stream, options?: zz_pXOptions): bigint[][]` completes equal-degree factorization.
  Here f is monic and square-free with all irreducible factors of degree d, and
  b is X^p modulo f.
- `SFCanZass2(groups: readonly (readonly [readonly bigint[], number])[], h, p, stream, options?: zz_pXOptions)`
  completes distinct-degree groups using their common Frobenius h, returning
  `bigint[][]`. This accepts the groups/Frobenius returned by `SFCanZass1`.
- `SFCanZass(f, p, stream, options?: zz_pXOptions): bigint[][]` factors a monic square-free polynomial,
  preserving native factor order. It handles constants and linear polynomials
  before the DDF/EDF stages.

The native preconditions determine termination and factorization guarantees.
Zero-degree guards preserve earlier native validation or work in earlier groups.
Errors can occur after consuming random bytes; the supplied stream retains that
native state. The word-polynomial division error is `zz_pX: division by zero`.

```ts
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { FindRoots, FindRoot } from '@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js';
FindRoots([0n,1n,1n], 2n, new RandomStream(new Uint8Array(32))); // [1n,0n]
FindRoot([0n,1n,1n], 2n, new RandomStream(new Uint8Array(32)));  // 1n
```

```ts
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { EDF, SFCanZass } from '@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js';
const f = [2n,1n,0n,1n,1n];
EDF(f, [0n,0n,0n,1n], 2, 3n, new RandomStream(new Uint8Array(32)));
// [[1n,0n,1n],[2n,1n,1n]]
SFCanZass(f, 3n, new RandomStream(new Uint8Array(32)));
// [[1n,0n,1n],[2n,1n,1n]]
```


### NTL degree-dependent word-polynomial settings

`zz_pXOptions`, exported as a type from `@sagemath-ts/ntl-ts/src/lzz_pX.js`,
is `{maxroot?: number; fftPrime?: number; state?: PolynomialProductState}`.
The default ordinary-context maxroot is 25.
Supply a nonnegative safe integer;
it determines the native context's algorithm thresholds and FFT size limit.
The constructor captures the settings. Operations taking F inherit them; the raw
GCD, sequence, root-product and factorization entries accept options last.
Supply `state: {context: FFTPrimeContext, stream: RandomStream}` to model FFT-table
initialization and its random draws. Nested operations inherit the same state;
when an operation also takes a stream, supply the same stream in both positions.
Omitting state retains the coefficient-only warmed-context behavior.

```ts
import { zz_pXModulus } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
const f = [1n, ...Array<bigint>(63).fill(0n), 3n];
const F = new zz_pXModulus(f, 9n, {maxroot: 60});
F.n;     // 64
F.f[64]; // 3n; this context uses the plain construction branch
```

```ts
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { SFCanZass } from '@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js';
const stream = new RandomStream(new Uint8Array(32));
SFCanZass([2n,1n,0n,1n,1n], 3n, stream, {maxroot: 60});
// [[1n,0n,1n],[2n,1n,1n]]
Array.from(stream.get(4)); // [93,106,229,83]
```


### NTL stateful small-prime selection for integer polynomials

Import `LocalInfoT`, `SmallPrimeFactorization`, and the type
`SmallPrimeFactorizationOptions` from `@sagemath-ts/ntl-ts/src/ZZXFactoring.js`.

```ts
SmallPrimeFactorization(info: LocalInfoT, f: readonly bigint[], stream: RandomStream,
  options?: {InitNumPrimes?: number; MaxNumPrimes?: number; context?: FFTPrimeContext}): bigint[][] | null
```

Use a fresh `new LocalInfoT()` for a new primitive squarefree integer polynomial.
The defaults are seven initial good primes and a maximum of fifty; unsuitable
primes are skipped without a fifty-attempt limit. Tuning values must be exactly
represented integers. A null result takes the native irreducibility exit; otherwise
results are monic modular factors in native order over `info.context!.p`.

The state exposes number fields `n`, `NumPrimes`, `NumFactors`, bigint
`PossibleDegrees`, a readonly `s: PrimeSeq`, and
`context: Readonly<{p: bigint; maxroot: number}> | null`. Its prime-vector getter
returns `(bigint | undefined)[]`, and pattern getter returns `number[][]`; both are
snapshots. Assigning either property copies its visible prefix and retains older
initialized storage beyond that length. Read only the first NumPrimes entries as
active retained prime information. Unwritten native scalar slots have no portable
value and appear as undefined in the port. Repeated calls preserve native state;
changing the polynomial without resetting information need not preserve the
mathematical guarantees of fresh factorization. These are explicit contexts and
streams. Supply `options.context` to include native FFT initialization using the
same stream; omit it for the earlier warmed-context coefficient adapter.

```ts
import { LocalInfoT, SmallPrimeFactorization }
  from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const info = new LocalInfoT(), stream = new RandomStream(new Uint8Array(32));
SmallPrimeFactorization(info, [-1n,0n,1n], stream); // [[1n,1n],[2n,1n]]
info.context;   // {p: 3n, maxroot: 2}
info.NumPrimes; // 6
Array.from(stream.get(4)); // [184,224,173,160]
```

```ts
import { LocalInfoT, SmallPrimeFactorization }
  from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const info = new LocalInfoT();
info.p = [101n,103n,107n]; info.p = [2n];
info.pattern = [[9],[8],[7]]; info.pattern = [[]];
info.s.reset(3n);
SmallPrimeFactorization(info, [-2n,0n,1n], new RandomStream(new Uint8Array(32)),
  {InitNumPrimes: 3}); // null
info.p;       // [3n,103n,107n]
info.pattern; // [[0,0,1],[8],[7]]
```


### NTL retained factor-information updates

Import `UpdateLocalInfo` and its type `UpdateLocalInfoOptions` from
`@sagemath-ts/ntl-ts/src/ZZXFactoring.js`.

```ts
UpdateLocalInfo(info: LocalInfoT, W: readonly (readonly bigint[])[],
  factors: readonly (readonly bigint[])[], f: readonly bigint[], count: number,
  modulus: bigint, options?: {van_Hoeij?: number; MaxNumPrimes?: number;
    state?: PolynomialProductState}): bigint[] | null
```

W contains the remaining modular factors over modulus (at least two); factors
contains recovered integer factors, and f is the remaining integer polynomial.
The routine updates info and returns recomputed subset-degree suffixes, or null
when that cache is unchanged. Use `updated ?? previous` to retain a caller cache.
Options replace native integer tuning globals: van_Hoeij defaults to 1 (any
nonzero value enables it), and MaxNumPrimes defaults to 50. Inputs are copied;
partial info changes survive errors. The saved info.context is preserved.
The optional shared state includes native FFT-cache creation and random draws.

```ts
import { LocalInfoT, SmallPrimeFactorization, UpdateLocalInfo }
  from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const info = new LocalInfoT();
SmallPrimeFactorization(info, [-1n,0n,1n], new RandomStream(new Uint8Array(32)));
UpdateLocalInfo(info, [[1n,1n]], [[-1n,1n]], [1n,1n], 1, 125n); // [2n]
info.NumFactors;     // 1
info.PossibleDegrees; // 3n
info.context;        // {p: 3n, maxroot: 2}
```

```ts
import { LocalInfoT, UpdateLocalInfo }
  from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
const info = new LocalInfoT();
info.n = 2; info.PossibleDegrees = 7n;
UpdateLocalInfo(info, [], [], [-1n,0n,1n], 1, 125n,
  {van_Hoeij: 0, MaxNumPrimes: 3}); // null
info.p;         // [3n]
info.NumPrimes; // 1
info.context;  // null; the previous context was restored
```


### NTL integer-factor recombination by cardinality

Import `CardinalitySearch`, `CardinalitySearch1`, `FindTrueFactors`, and the type
`FactorRecombinationOptions` from `@sagemath-ts/ntl-ts/src/ZZXFactoring.js`.

```ts
CardinalitySearch(recovered: readonly (readonly bigint[])[], f: readonly bigint[],
  W: readonly (readonly bigint[])[], info: LocalInfoT, count: number, bound: number,
  modulus: bigint, options?: UpdateLocalInfoOptions): [bigint[][], bigint[], bigint[][]]
CardinalitySearch1(recovered: readonly (readonly bigint[])[], f: readonly bigint[],
  W: readonly (readonly bigint[])[], info: LocalInfoT, count: number, bound: number,
  modulus: bigint, options?: FactorRecombinationOptions): [bigint[][], bigint[], bigint[][]]
FindTrueFactors(f: readonly bigint[], W: readonly (readonly bigint[])[], modulus: bigint,
  info: LocalInfoT, bound: number, options?: FactorRecombinationOptions): bigint[][]
```

These are recombination helpers: W contains positive-degree lifted monic modular
factors over modulus (at least two), and info describes the original polynomial.
The search tuple gives recovered factors, remaining f, and remaining W. Inputs
are copied; info retains native partial updates. `bound` limits candidate coefficient
bit lengths. The optimized search requires count greater than one. Options use
native integer defaults `{van_Hoeij: 1, MaxNumPrimes: 50, MaxPrune: 10}`; MaxPrune
applies to the optimized search and FindTrueFactors. Both options types also accept
`state?: PolynomialProductState`, which is forwarded to word-context initialization,
modular products, exact division and local updates. Results retain native order.

```ts
import { LocalInfoT, CardinalitySearch } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
const info = new LocalInfoT(); info.n = 4; info.PossibleDegrees = 31n;
const W = [[-1n,1n],[1n,1n],[2n,1n],[3n,1n]];
CardinalitySearch([], [-6n,-5n,5n,5n,1n], W, info, 1, 100, 2n**127n-1n);
// [[[-1n,1n],[1n,1n],[2n,1n]], [3n,1n], [[3n,1n]]]
```

```ts
import { LocalInfoT, CardinalitySearch1 } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
const info = new LocalInfoT(); info.n = 4; info.PossibleDegrees = 31n;
const W = [[-1n,1n],[1n,1n],[2n,1n],[3n,1n]];
CardinalitySearch1([], [-6n,-5n,5n,5n,1n], W, info, 2, 100, 2n**127n-1n);
// [[[-1n,0n,1n]], [6n,5n,1n], [[2n,1n],[3n,1n]]]
```

```ts
import { LocalInfoT, FindTrueFactors } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
const info = new LocalInfoT(); info.n = 4; info.PossibleDegrees = 31n;
const W = [[-1n,1n],[1n,1n],[2n,1n],[3n,1n]];
FindTrueFactors([-6n,-5n,5n,5n,1n], W, 2n**127n-1n, info, 100);
// [[-1n,1n],[1n,1n],[2n,1n],[3n,1n]]
```


### NTL probable-prime checks and generators

Import these free functions from `@sagemath-ts/ntl-ts/src/ZZ.js`:

```ts
ComputePrimeBound(bits: number): bigint
ErrBoundTest(bits: number, trials: number, errorBits: number): number
ProbPrime(n: bigint, stream: RandomStream,
  options?: {NumTrials?: number; word?: boolean}): number
RandomPrime(bits: number, stream: RandomStream, options?: {NumTrials?: number}): bigint
OldRandomPrime(bits: number, stream: RandomStream, options?: {NumTrials?: number}): bigint
RandomPrime_long(bits: number, stream: RandomStream, options?: {NumTrials?: number}): bigint
GenPrime(bits: number, stream: RandomStream, options?: {err?: number}): bigint
GenPrime_long(bits: number, stream: RandomStream, options?: {err?: number}): bigint
```

ProbPrime returns 0 or 1 and consumes the explicit stream when it needs random
witnesses. NumTrials defaults to 10; negative values mean zero random trials,
with the initial base-two check retained. GenPrime uses err=80 by default and
clamps it to 1–512. ErrBoundTest returns native bound codes 0–5 (zero means the
bound is insufficient). These are probable-prime routines, not primality proofs.
All generators require at least two bits; word generators allow at most 63 bits.
GenPrime allows at most 2^20 bits. Number parameters must be exactly represented
native signed integers. The historical generator stays sequential; RandomPrime
uses seeded candidate blocks at 256 bits and above.

```ts
import { ComputePrimeBound, ErrBoundTest, ProbPrime, RandomStream }
  from '@sagemath-ts/ntl-ts/src/ZZ.js';
ComputePrimeBound(2305); // 23050n
ErrBoundTest(128, 8, 80); // 0
const stream = new RandomStream(new Uint8Array(32));
ProbPrime(17n, stream, {NumTrials: 2}); // 1
Array.from(stream.get(4)); // [64, 93, 106, 229]
```

```ts
import { RandomPrime, OldRandomPrime, RandomPrime_long, RandomStream }
  from '@sagemath-ts/ntl-ts/src/ZZ.js';
const fresh = () => new RandomStream(new Uint8Array(32));
RandomPrime(32, fresh());     // 2419978657n
OldRandomPrime(32, fresh());  // 2419978657n
RandomPrime_long(32, fresh()); // 2915387179n
```

```ts
import { GenPrime, GenPrime_long, RandomPrime, RandomStream }
  from '@sagemath-ts/ntl-ts/src/ZZ.js';
const fresh = () => new RandomStream(new Uint8Array(32));
GenPrime(32, fresh());      // 2419978657n
GenPrime_long(60, fresh()); // 1007410285329247537n
RandomPrime(256, fresh(), {NumTrials: 2});
// 99719095675712960086022384902837345215225263017301052027152607768417317746133n
```

The package root also re-exports ProbPrime and RandomPrime. Their older
`ProbPrime(n: ZZ, NumTrials?: number)` and
`RandomPrime(bits: number, NumTrials?: number): ZZ` overloads retain their existing
`NTL_NOT_IMPLEMENTED` errors; use the bigint/explicit-stream overloads above.


### NTL word matrix elimination, inverses and solving

Import from `@sagemath-ts/ntl-ts/src/mat_lzz_p.js`. Matrices are rectangular bigint
rows and p satisfies 1 < p < 2^60. Dimensions use exactly represented native
integers. Set columns for a zero-row matrix; returned arrays do not store empty
column metadata.

```ts
type WordMatrixOptions = {columns?: number};
type WordMatrixEliminationOptions = WordMatrixOptions & {w?: number};
type WordMatrixInverseOptions = WordMatrixOptions & {
  relax?: boolean; previous?: readonly (readonly bigint[])[]
};
type WordMatrixSolveOptions = WordMatrixOptions & {
  relax?: boolean; left?: boolean; previous?: readonly bigint[]
};
gauss(A: readonly (readonly bigint[])[], p: bigint,
  options?: WordMatrixEliminationOptions): [number, bigint[][]]
image(A: readonly (readonly bigint[])[], p: bigint, options?: WordMatrixOptions): bigint[][]
kernel(A: readonly (readonly bigint[])[], p: bigint, options?: WordMatrixOptions): bigint[][]
relaxed_inv(A: readonly (readonly bigint[])[], p: bigint,
  options?: WordMatrixInverseOptions): [bigint, bigint[][]]
inv(A: readonly (readonly bigint[])[], p: bigint, options?: WordMatrixOptions): bigint[][]
relaxed_determinant(A: readonly (readonly bigint[])[], p: bigint,
  options?: WordMatrixOptions & {relax?: boolean}): bigint
determinant(A: readonly (readonly bigint[])[], p: bigint, options?: WordMatrixOptions): bigint
relaxed_solve(A: readonly (readonly bigint[])[], b: readonly bigint[], p: bigint,
  options?: WordMatrixSolveOptions): [bigint, bigint[]]
solve(A: readonly (readonly bigint[])[], b: readonly bigint[], p: bigint,
  options?: Omit<WordMatrixSolveOptions, 'relax'>): [bigint, bigint[]]
```

The four option types are exported. gauss returns [rank, fullMatrix], using w
pivot columns (all columns by default). image returns pivot rows without pivot
normalization; kernel returns a **left** kernel basis. Inputs and outputs are
independent arrays.

```ts
import {gauss, image, kernel} from '@sagemath-ts/ntl-ts/src/mat_lzz_p.js';
const A = [[1n,2n],[2n,4n],[0n,1n]];
gauss(A, 5n);  // [2, [[1n,2n],[0n,1n],[0n,0n]]]
image(A, 5n);  // [[1n,2n],[0n,1n]]
kernel(A, 5n); // [[3n,1n,0n]]; each row v satisfies v*A=0
```

relaxed_inv returns [determinantStatus, inverseOrPrevious]. Relaxed routines skip
nonunit pivots by default; strict inv/determinant/solve throw at the first such
pivot. inv also throws for a singular matrix. Zero status in relaxed routines
preserves the previous output, defaulting to an empty array.

```ts
import {relaxed_inv, inv, relaxed_determinant, determinant}
  from '@sagemath-ts/ntl-ts/src/mat_lzz_p.js';
const A = [[1n,2n],[3n,4n]];
relaxed_inv(A, 5n); // [3n, [[3n,1n],[4n,2n]]]
inv(A, 5n); // [[3n,1n],[4n,2n]]
relaxed_determinant(A, 5n); // 3n
determinant(A, 5n); // 3n
```

Solve entries return [determinantStatus, solutionOrPrevious]. By default they solve
A*x=b; left=true solves x*A=b.

```ts
import {relaxed_solve, solve} from '@sagemath-ts/ntl-ts/src/mat_lzz_p.js';
const A = [[1n,2n],[3n,4n]], b = [1n,0n];
relaxed_solve(A, b, 5n); // [3n, [3n,4n]]
solve(A, b, 5n); // [3n, [3n,4n]]
relaxed_solve(A, b, 5n, {left: true}); // [3n, [3n,1n]]
```

Composite moduli preserve NTL's pivot behavior, including zero status when it
cannot find a unit pivot; that status need not equal the mathematical determinant
over a general composite ring.

```ts
import {relaxed_inv, inv, relaxed_solve} from '@sagemath-ts/ntl-ts/src/mat_lzz_p.js';
relaxed_inv([[3n,1n],[1n,1n]], 9n); // [2n, [[5n,4n],[4n,6n]]]
// inv([[3n,1n],[1n,1n]], 9n) throws Error('InvMod: inverse undefined').
relaxed_inv([[3n]], 9n, {previous: [[7n]]}); // [0n, [[7n]]]
relaxed_solve([[3n]], [1n], 9n, {previous: [11n,12n]}); // [0n, [2n,3n]]
```


### NTL FFT-prime caches

Import from `@sagemath-ts/ntl-ts/src/FFT.js` and import RandomStream from the
corresponding `ZZ.js` module. One FFTPrimeContext represents one native cache;
share it across calls to preserve cache reuse and stream consumption.

```ts
interface FFTPrimeInfo {
  readonly q: bigint;
  readonly qrecip: number;
  readonly RootTable: readonly [readonly bigint[], readonly bigint[]];
  readonly TwoInvTable: readonly bigint[];
}
class FFTPrimeContext {
  constructor();
  length(): number;
  get(index: number): FFTPrimeInfo;
}
IsFFTPrime(n: bigint, stream: RandomStream, previousRoot?: bigint): [number, bigint]
CalcMaxRoot(p: bigint): number
NextFFTPrime(index: number, context: FFTPrimeContext, stream: RandomStream): [bigint, bigint]
InitFFTPrimeInfo(q: bigint, root: bigint): FFTPrimeInfo
UseFFTPrime(index: number, context: FFTPrimeContext, stream: RandomStream): void
GetFFTPrime(index: number, context: FFTPrimeContext): bigint
GetFFTPrimeRecip(index: number, context: FFTPrimeContext): number
```

IsFFTPrime returns [1, root] on success and [0, previousRoot] on failure, with
previousRoot defaulting to 0n. Inputs and previous roots fit signed 64-bit words.
The cache uses primes below 2^60, at most 25 root levels, and indices 0–19,999.
Number indices must be exactly represented native integers. Reading an entry
before initialization raises RangeError. CalcMaxRoot rejects p=1, where the native
loop does not terminate, and signed-overflow inputs.

```ts
import {IsFFTPrime, CalcMaxRoot} from '@sagemath-ts/ntl-ts/src/FFT.js';
import {RandomStream} from '@sagemath-ts/ntl-ts/src/ZZ.js';
const stream = new RandomStream(new Uint8Array(32));
IsFFTPrime(97n, stream); // [1, 46n]
IsFFTPrime(21n, stream, 123n); // [0, 123n]
CalcMaxRoot(97n); // 5
```

InitFFTPrimeInfo accepts 1 < q < 2^60 and a canonical root 0 <= root < q. Native
nonunit inverses still raise errors. Information objects and nested arrays are
immutable; RootTable contains successive two-power roots and their inverses.

```ts
import {InitFFTPrimeInfo} from '@sagemath-ts/ntl-ts/src/FFT.js';
const info = InitFFTPrimeInfo(17n, 3n);
info.RootTable; // [[1n,16n,13n,9n,3n], [1n,16n,4n,2n,6n]]
info.TwoInvTable; // [1n,9n,13n,15n,16n]
```

```ts
import {FFTPrimeContext, UseFFTPrime, GetFFTPrime, GetFFTPrimeRecip}
  from '@sagemath-ts/ntl-ts/src/FFT.js';
import {RandomStream} from '@sagemath-ts/ntl-ts/src/ZZ.js';
const context = new FFTPrimeContext(), stream = new RandomStream(new Uint8Array(32));
context.length(); // 0
UseFFTPrime(0, context, stream);
context.length(); // 1
GetFFTPrime(0, context); // 882705526964617217n
GetFFTPrimeRecip(0, context); // 1.1328806373726086e-18
context.get(0).RootTable[0][0]; // 1n
UseFFTPrime(0, context, stream); // cached: consumes no random bytes
Array.from(stream.get(4)); // [152,186,151,124]
```

NextFFTPrime is the underlying candidate helper. It advances search state without
adding a cache entry; repeating its index rolls back the candidate search while
continuing to consume the caller's stream.

```ts
import {FFTPrimeContext, NextFFTPrime, UseFFTPrime, GetFFTPrime}
  from '@sagemath-ts/ntl-ts/src/FFT.js';
import {RandomStream} from '@sagemath-ts/ntl-ts/src/ZZ.js';
const context = new FFTPrimeContext(), stream = new RandomStream(new Uint8Array(32));
NextFFTPrime(0, context, stream)[0]; // 882705526964617217n
context.length(); // 0
NextFFTPrime(0, context, stream)[0]; // 882705526964617217n
UseFFTPrime(0, context, stream);
GetFFTPrime(0, context); // 882705526964617217n
```


### NTL incremental CRT and modular determinants

These lower-level NTL exports use BigInt values and matrices as arrays of rows.
`CRTInRange(g, a)` returns 0 or 1 for `-a/2 < g <= a/2` with `a > 0`.
Scalar `CRT(g, a, G, p, {word?: boolean})` returns
`[modified: number, residue: bigint, modulus: bigint]`. It requires `a > 0`,
`p > 1`, `0 <= G < p`; noncoprime moduli throw the native inverse error.
The word option uses signed 64-bit inputs for G and p.

```ts
import { CRT, CRTInRange } from '@sagemath-ts/ntl-ts/src/ZZ.js';
CRTInRange(-3n, 6n); // 0
CRTInRange(3n, 6n);  // 1
CRT(1n, 3n, 0n, 2n); // [1, -2n, 6n]
CRT(1n, 3n, 1n, 2n); // [0, 1n, 6n]
```

Matrix `CRT(g, a, G, p, {columns?, residueColumns?})` returns the same tuple with
a matrix residue. The word modulus satisfies `1 < p < 2^60`, entries of G are
normalized modulo p, and the two matrices must have equal dimensions. Optional
column counts preserve zero-row shapes; otherwise they come from the first row.
Input arrays stay unchanged, and output rows are independently owned.

```ts
import { CRT } from '@sagemath-ts/ntl-ts/src/mat_ZZ.js';
CRT([[1n, -1n]], 3n, [[0n, 1n]], 2n); // [1, [[-2n, -1n]], 6n]
```

`DetBound(A, {columns?})` returns the bit length of the native product of row-norm
bounds. It accepts rectangular integer matrices. An empty matrix has bound one;
a zero row makes the bound zero.

```ts
import { DetBound } from '@sagemath-ts/ntl-ts/src/mat_ZZ.js';
DetBound([[3n, 4n], [0n, 1n]]); // 3
```

`determinant(A, p, {columns?})` in `mat_ZZ_p.js` requires a square matrix and any
BigInt modulus `p > 1`. It returns a canonical residue, with empty determinant one.
It preserves the first nonzero pivot and throws if that pivot is a nonunit.

```ts
import { determinant } from '@sagemath-ts/ntl-ts/src/mat_ZZ_p.js';
const p = (1n << 127n) - 1n;
determinant([[p + 2n, -1n], [1n, p + 3n]], p); // 7n
```


### NTL integer determinant and inverse reconstruction

From `@sagemath-ts/ntl-ts/src/mat_ZZ.js`:

- `determinant(A, context, stream, {columns?, deterministic?})` returns a BigInt.
- `inv(A, context, stream, {columns?})` returns the integral inverse. It requires
  determinant 1 or -1 and otherwise throws `inv: non-invertible matrix`.
- `inv(A, context, stream, {status: true, columns?, previous?, deterministic?})`
  returns `[determinant, adjugateOrPrevious]`. Singular input retains a copy of
  previous (default `[]`); empty input returns `[1n, []]`.

A is an array of BigInt rows. The explicit `FFTPrimeContext` and `RandomStream`
are shared across calls; cached primes consume no additional initialization bytes.
By default these routines use NTL's probabilistic stabilization. Set
`deterministic: true` on determinant or status inverse to use the native
full-bound reconstruction path. Strict inverse follows the native probabilistic
default. Matrices and previous outputs remain unchanged; returned rows are owned.

```ts
import { determinant, inv } from '@sagemath-ts/ntl-ts/src/mat_ZZ.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const context = new FFTPrimeContext();
const stream = new RandomStream(new Uint8Array(32));
const A = [[1n, 2n], [3n, 4n]];
determinant(A, context, stream); // -2n
inv(A, context, stream, {status: true}); // [-2n, [[4n, -2n], [-3n, 1n]]]
```

```ts
import { inv } from '@sagemath-ts/ntl-ts/src/mat_ZZ.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const context = new FFTPrimeContext();
const stream = new RandomStream(new Uint8Array(32));
inv([[1n, 2n], [3n, 5n]], context, stream); // [[-5n, 2n], [3n, -1n]]
```

```ts
import { inv } from '@sagemath-ts/ntl-ts/src/mat_ZZ.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const context = new FFTPrimeContext();
const stream = new RandomStream(new Uint8Array(32));
inv([[1n, 2n], [2n, 4n]], context, stream, {status: true, previous: [[7n]]});
// [0n, [[7n]]]
```

```ts
import { determinant } from '@sagemath-ts/ntl-ts/src/mat_ZZ.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const context = new FFTPrimeContext();
const stream = new RandomStream(new Uint8Array(32));
determinant([[1n, 1n << 1024n], [0n, 1n]], context, stream, {deterministic: true}); // 1n
context.length(); // 18
```


### NTL certified integer row reduction

`gauss(M, context, stream, {columns?})` from
`@sagemath-ts/ntl-ts/src/ZZXFactoring.js` returns `[d: bigint, R: bigint[][]]`, where
`R/d` is reduced row echelon form. M must be nonempty with linearly independent
rows. Dependent nonempty input can retry indefinitely, matching the original
helper's precondition. This is the factorization helper; its output scale d is the
signed determinant of selected original pivot columns.

It shares an explicit FFTPrimeContext and RandomStream, preserves the native
prime sampling and retry schedule, and leaves M unchanged. Each call samples a
new word prime even when the FFT cache is already initialized.

```ts
import { gauss } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const context = new FFTPrimeContext();
const stream = new RandomStream(new Uint8Array(32));
gauss([[1n, 2n], [3n, 4n]], context, stream);
// [-2n, [[-2n, 0n], [0n, -2n]]]
```

```ts
import { gauss } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const context = new FFTPrimeContext();
const stream = new RandomStream(new Uint8Array(32));
gauss([[1n, 2n, 3n], [0n, 1n, 4n]], context, stream);
// [1n, [[1n, 0n, -5n], [0n, 1n, 4n]]]
```


### NTL modular products with explicit random state

`ZZ_pX.mul(a, b, p, state?)` accepts coefficient arrays and optional
`PolynomialProductState = { context: FFTPrimeContext, stream: RandomStream }`.
Sharing state preserves native FFT-prime initialization and subsequent random bytes.
Passing the same array twice selects native squaring thresholds. Three-argument
calls retain the existing coefficient-only behavior.

```ts
import { mul } from '@sagemath-ts/ntl-ts/src/ZZ_pX.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const context = new FFTPrimeContext();
const stream = new RandomStream(new Uint8Array(32));
const p = (1n << 127n) - 1n;
const a = Array<bigint>(80).fill(1n);
mul(a, a.slice(), p, { context, stream });
context.length(); // 0: distinct operands use the multiplication threshold
const square = mul(a, a, p, { context, stream });
[square.length, square[0], square[79], square[158]]; // [159, 1n, 80n, 1n]
context.length(); // 5
```

`ZZXFactoring.mul(W, p, indices?, state?)` shares the same state across its
minimum-degree product tree. It copies selected factors and preserves caller arrays.

```ts
import { mul } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const context = new FFTPrimeContext();
const stream = new RandomStream(new Uint8Array(32));
const W = Array.from({ length: 4 }, () => Array<bigint>(101).fill(1n));
const product = mul(W, (1n << 127n) - 1n, [3, 1, 3, 0], { context, stream });
[product.length, product[0], product[400]]; // [401, 1n, 1n]
context.length(); // 5
```

`ZZX1.SSRatio(na, maxa, nb, maxb)` returns the original floating padding ratio
used to select Schoenhage–Strassen multiplication. Arguments are exactly represented
signed native integers; undefined signed-overflow inputs have explicit adapter guards.

```ts
import { SSRatio } from '@sagemath-ts/ntl-ts/src/ZZX1.js';
SSRatio(1, 1, 1, 1); // 1.4
SSRatio(0, 127, 199, 127); // 0
```


### NTL integer polynomial exact division

`ZZX1.divide(a, b, options?)`, `PlainDivide(a, b, options?)` and
`HomDivide(a, b, options?)` return `[status, quotientOrPrevious]`, with status 0 or 1.
Coefficient arrays are constant term first. `divide` also accepts a bigint scalar
as `b`. `IntegerPolynomialDivisionOptions` has optional `previous` and
`state: { context: FFTPrimeContext, stream: RandomStream }`; outputs own their arrays.

The ordinary dispatcher uses PlainDivide when the divisor degree or quotient-degree
bound is at most eight. Its zero-numerator behavior follows the original, including
the difference from explicitly calling HomDivide:

```ts
import { divide, PlainDivide, HomDivide } from '@sagemath-ts/ntl-ts/src/ZZX1.js';
const options = { previous: [17n, -2n] };
divide([], [1n, 1n], options); // [0, [17n, -2n]]
PlainDivide([], [1n, 1n], options); // [0, [17n, -2n]]
HomDivide([], [1n, 1n], options); // [1, []]
divide([], [], options); // [1, []]
```

```ts
import { divide } from '@sagemath-ts/ntl-ts/src/ZZX1.js';
divide([2n, -4n], -2n); // [1, [-1n, 2n]]
divide([1n, 2n], 2n, { previous: [7n] }); // [0, [7n]]
```

HomDivide reconstructs the quotient with FFT primes and certifies the result using
the native coefficient bound. Supplying state preserves prime-search random draws;
coefficient-only callers may omit it.

```ts
import { HomDivide } from '@sagemath-ts/ntl-ts/src/ZZX1.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const context = new FFTPrimeContext();
const stream = new RandomStream(new Uint8Array(32));
HomDivide([1n, 2n, 1n], [1n, 1n], { state: { context, stream } }); // [1, [1n, 1n]]
context.length(); // 2
```


The internal `_ZZX_kernels.divide(a, b, state?)` adapter accepts polynomial arrays
and returns the quotient or `null`, delegating to the same ordinary dispatcher.
For example, a zero numerator with divisor `[1n, 1n]` returns `null`. Its optional
third argument is the same `PolynomialProductState` used above.


### NTL van Hoeij factor acceptance

`ZZXFactoring.GotThem(previous, B_L, W, f, bound, modulus, context, stream, options?)`
returns `[status, factors]`. It appends accepted factors to a copy of `previous`, or
returns an unchanged copy when status is 0. `B_L` uses rows and must have independent
rows; `options.columns` preserves an explicit column count for empty matrices. `W`
must contain at least one lifted polynomial per column. Bounds are exactly represented
signed native integers. Both polynomial products and integer division share the supplied
FFTPrimeContext and RandomStream.

```ts
import { GotThem } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const B = [[1n,1n,1n,1n,0n,0n,0n,0n], [0n,0n,0n,0n,1n,1n,1n,1n]];
const W = [...Array.from({ length: 4 }, () => [1n,1n]),
           ...Array.from({ length: 4 }, () => [-1n,1n])];
const f = [1n,0n,-4n,0n,6n,0n,-4n,0n,1n];
const context = new FFTPrimeContext();
const stream = new RandomStream(new Uint8Array(32));
GotThem([], B, W, f, 16, 65537n, context, stream);
// [1, [[1n,4n,6n,4n,1n], [1n,-4n,6n,-4n,1n]]]
context.length(); // 3
```

Groups with at most three lifted factors are rejected by this acceptance helper:

```ts
import { GotThem } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
GotThem([[7n]], [[1n,1n,1n]], [[1n,1n],[1n,1n],[1n,1n]],
        [1n,3n,3n,1n], 16, 65537n,
        new FFTPrimeContext(), new RandomStream(new Uint8Array(32)));
// [0, [[7n]]]
```

The original's group-degree keys are equal, so accepted factors retain certified row
order. Trial-division rejection emits `X` and a newline; the port uses `console.error`.
Verbose timing output is not exposed by this adapter.


### NTL cold word-context initialization

A shared `FFTPrimeContext` retains native primes, roots and reciprocal tables.
Initializing a word modulus consumes the supplied stream only for missing primes.
The native four-prime limit is checked after generating every required prime;
a failed constructor can therefore extend the cache and advance the stream.

```ts
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { LocalInfoT, SmallPrimeFactorization }
  from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
const context = new FFTPrimeContext();
const stream = new RandomStream(Uint8Array.from({length:32}, (_, i) => i));
SmallPrimeFactorization(new LocalInfoT(), [-1n,0n,1n], stream, {context});
// [[2n,1n],[1n,1n]]
context.length(); // 1
Array.from(stream.get(4)); // [246,94,222,93]
```

```ts
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { zz_pXModulus } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
const state = {context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32))};
try { new zz_pXModulus(null, 3n, {maxroot:240, state}); }
catch (e) { (e as Error).message; } // 'zz_pInit: too many primes'
state.context.length(); // 5
const F = new zz_pXModulus(null, 3n, {state});
F.PrimeCnt; // 1; reuses the initialized cache
```


### NTL factor recombination with shared cache and stream

`UpdateLocalInfo`, `CardinalitySearch`, `CardinalitySearch1` and `FindTrueFactors`
accept `options.state: PolynomialProductState`. `InvMul` accepts the same state
as its fourth argument. Reuse the state from small-prime selection to preserve
native random draws throughout factor recovery.

```ts
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { LocalInfoT, UpdateLocalInfo } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
const state = {context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32))};
const info = new LocalInfoT(); info.n = 2; info.PossibleDegrees = 7n;
UpdateLocalInfo(info, [], [], [-1n,0n,1n], 1, 125n,
  {van_Hoeij:0, MaxNumPrimes:3, state}); // null
info.p; // [3n]
state.context.length(); // 1
Array.from(state.stream.get(4)); // [152,186,151,124]
```

```ts
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { LocalInfoT, CardinalitySearch } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
const state = {context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32))};
const info = new LocalInfoT(); info.n = 20; info.PossibleDegrees = (1n<<21n)-1n;
const g = [1n,...Array<bigint>(9).fill(0n),1n];
const h = [2n,...Array<bigint>(9).fill(0n),1n];
const f = [2n,...Array<bigint>(9).fill(0n),3n,...Array<bigint>(9).fill(0n),1n];
CardinalitySearch([], f, [g,h], info, 1, 100, 65537n, {state}); // [[g],h,[h]]
state.context.length(); // 2; exact division initialized the shared primes
```

```ts
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { InvMul } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
const state = {context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32))};
const W = Array.from({length:4}, () => Array<bigint>(101).fill(1n));
const f = InvMul(W, [], 65537n, state);
[f.length, f[0], f[400]]; // [401,1n,1n]
state.context.length(); // 2
```


### NTL integer polynomial product state and transform choice

The explicit state follows NTL's plain/Karatsuba/transform dispatch. Homomorphic
products initialize enough FFT primes to certify the coefficient bound. Integer
Schönhage–Strassen products do not draw from the FFT-prime stream. These state
effects accompany the port's existing exact coefficient kernel.

```ts
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { mul } from '@sagemath-ts/ntl-ts/src/ZZX1.js';
const state = {context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32))};
const a = Array<bigint>(150).fill(1n);
const c = mul(a, a, state);
[c.length, c[0], c[149], c[298]]; // [299,1n,150n,1n]
state.context.length(); // 1
```

```ts
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { sqr } from '@sagemath-ts/ntl-ts/src/ZZX1.js';
const state = {context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32))};
sqr(Array<bigint>(80).fill(1n<<1919n), state).length; // 159
state.context.length(); // 0
```

```ts
import { ChooseSS } from '@sagemath-ts/ntl-ts/src/ZZX1.js';
ChooseSS(149, 1, 149, 1);          // false
ChooseSS(79, 1920, 79, 1920);      // true
ChooseSS(79, 3400, 8192, 3400);    // true
ChooseSS(79, 6000, 16384, 6000);   // true
ChooseSS(79, 10000, 65536, 10000); // true
ChooseSS(79, 17500, 65536, 17500); // true
```


### NTL truncated transforms and retained word-modulus caches

From `@sagemath-ts/ntl-ts/src/FFT.js`:

```ts
FFTFwd_trunc(a: readonly bigint[], k: number, info: FFTPrimeInfo, yn: number, xn: number): bigint[];
FFTRev1_trunc(a: readonly bigint[], k: number, info: FFTPrimeInfo, yn: number): bigint[];
```

From `@sagemath-ts/ntl-ts/src/FFT_impl.js`:

```ts
FFTRoundUp(xn: number, k: number): number;
```

Transforms use bit-reversed frequency order and return independent arrays of length
yn. Prefix lengths must equal FFTRoundUp(length,k), lie in [1,2^k], and fit the input
buffer and prime root table. Forward inputs have xn coefficients; inverse inputs
have yn frequency values. Coefficients are reduced modulo info.q. FFTRoundUp accepts
exactly represented signed native integer lengths and exponents from 0 through 62.

<!-- llm-test: ntl-truncated-transform-roundtrip -->
```ts
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { FFTPrimeContext, UseFFTPrime, FFTFwd_trunc, FFTRev1_trunc }
  from '@sagemath-ts/ntl-ts/src/FFT.js';
const context = new FFTPrimeContext();
UseFFTPrime(0, context, new RandomStream(new Uint8Array(32)));
const a = [1n, 2n, 3n, 4n], info = context.get(0);
const frequencies = FFTFwd_trunc(a, 2, info, 4, 4);
console.log(FFTRev1_trunc(frequencies, 2, info, 4)); // [1n, 2n, 3n, 4n]
```

<!-- llm-test: ntl-transform-prefix-rounding -->
```ts
import { FFTRoundUp } from '@sagemath-ts/ntl-ts/src/FFT_impl.js';
console.log(FFTRoundUp(401, 9)); // 416
```

Word modulus and multiplier build failures retain NTL's partially updated caches.
Later reductions may succeed using an old reciprocal, or raise `FFT rep mismatch`
when no compatible cache exists. Supply `{state:{context,stream}}` to zz_pXModulus
to match ambient FFT roots exactly. Without state, exceptional truncated interpolation
after failed rebuilds uses a private context seeded with a 32-byte zero key.

<!-- llm-test: ntl-word-failed-rebuild-remainder -->
```ts
import { zz_pXModulus, build, rem } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
const good = [1n, ...Array<bigint>(47).fill(0n), 1n];
const F = new zz_pXModulus(good, 6n);
try { build(F, [...good.slice(0, 48), 2n]); } catch {}
const a = [...Array<bigint>(94).fill(0n), 1n];
console.log(rem(a, F).flatMap((x, i) => x ? [[i, x]] : [])); // [[30, 5n], [46, 5n]]
```


### NTL word FFT residue reconstruction

From `@sagemath-ts/ntl-ts/src/lzz_pX.js`:

```ts
FromModularRep(residues: readonly (readonly bigint[])[], F: zz_pXModulus): bigint[];
```

Use one residue row per F.NumPrimes and one column per output coefficient. Signed
entries are reduced modulo the corresponding FFT prime. Equal row lengths are
required; zero columns are preserved. The one-prime case centers exactly; multiple
primes use the fused binary64 estimate of the bundled AArch64 native profile.
The return value owns its coefficient array.

<!-- llm-test: ntl-word-exact-crt-centering -->
```ts
import { zz_pXModulus, FromModularRep } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
const F = new zz_pXModulus(null, 6n);
console.log(FromModularRep([[441352763482308608n]], F)); // [2n]
```

### NTL transposed transforms and rebuilt projection caches

From `@sagemath-ts/ntl-ts/src/FFT.js`:

```typescript
FFTFwd_trans(a: readonly bigint[], k: number, info: FFTPrimeInfo): bigint[];
FFTRev1_trans(a: readonly bigint[], k: number, info: FFTPrimeInfo): bigint[];
```

Both require at least `2 ** k` input coefficients and return independent arrays of
that length. They are the transposes of the forward and normalized inverse FFT;
inputs are reduced modulo `info.q`. Invalid root exponents and short buffers throw
labeled `RangeError`s. Cached root powers belong to the supplied prime-info object.

```typescript
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
import { FFTPrimeContext, UseFFTPrime, FFTFwd_trans, FFTRev1_trans }
  from '@sagemath-ts/ntl-ts/src/FFT.js';
const context = new FFTPrimeContext();
UseFFTPrime(0, context, new RandomStream(new Uint8Array(32)));
const info = context.get(0);
console.log(FFTFwd_trans([1n, 0n, 0n, 0n], 2, info)); // [1n, 1n, 1n, 1n]
console.log(FFTRev1_trans([1n, 1n, 1n, 1n], 2, info)); // [1n, 0n, 0n, 0n]
```

`lzz_pX1.UpdateMap(a,B,F)` preserves native cached transforms across successful and
failed modulus/multiplier rebuilds, including native `FFT rep mismatch` errors.
For example, a plain rebuild still uses a previously cached FFT multiplier:

```typescript
import { zz_pXModulus, zz_pXMultiplier, build }
  from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
import { UpdateMap } from '@sagemath-ts/ntl-ts/src/lzz_pX1.js';
const F = new zz_pXModulus([1n, ...Array<bigint>(47).fill(0n), 1n], 6n);
const B = new zz_pXMultiplier([2n, ...Array<bigint>(46).fill(0n), 3n], F);
build(F, [1n, 1n]);
console.log(UpdateMap([1n], B, F)); // [2n]
```

The internal `lzz_pX._zz_pX_transform_kernels(F,B)` returns retained sizes `k,l` and
`RevTofftRep(a,k,offset)`, `mul(rep,'FRep'|'B1'|'B2')`,
`RevFromfftRep(rep,lo,hi)` and `AddExpand(target,small)`. A representation has shape
`{ k: number, rows: bigint[][] }`, with one full frequency row per FFT prime.
This adapter supports the mirrored projection module; callers should use UpdateMap.
Its returned arrays own their storage. Never mutate the internal cache fields.

### NTL arbitrary-modulus polynomial quotients

From `@sagemath-ts/ntl-ts/src/ZZ_pX.js`:

```typescript
new ZZ_pXModulus(f: readonly bigint[] | null, p: bigint, state?: PolynomialProductState);
new ZZ_pXMultiplier();
new ZZ_pXMultiplier(b: readonly bigint[], F: ZZ_pXModulus);
build(F: ZZ_pXModulus, f: readonly bigint[]): void;
build(B: ZZ_pXMultiplier, b: readonly bigint[], F: ZZ_pXModulus): void;
rem(a: readonly bigint[], F: ZZ_pXModulus): bigint[];
MulMod(a: readonly bigint[], b: readonly bigint[] | ZZ_pXMultiplier, F: ZZ_pXModulus): bigint[];
SqrMod(a: readonly bigint[], F: ZZ_pXModulus): bigint[];
FromModularRep(residues: readonly (readonly bigint[])[], F: ZZ_pXModulus): bigint[];
```

Coefficient arrays run from constant term upward. The explicit coefficient modulus
`p` must exceed one; composite moduli retain NTL's noninvertible-element errors.
Passing `null` creates an uninitialized quotient. `F.val()` and `B.val()` return
independent coefficient arrays. Builds preserve native partial changes after errors,
including the public polynomial and cached FFT representations. Degree 22 starts
quotient FFT initialization. Optional `state` shares `FFTPrimeContext` and
`RandomStream`; omitted state uses a private context with a 32-byte zero key.
The existing object-valued `MulMod(ZZ_pX, ZZ_pX, ZZ_pX)` overload remains separate.

<!-- llm-test: ntl-big-quotient-arithmetic -->
```typescript
import { ZZ_pXModulus, ZZ_pXMultiplier, rem, MulMod, SqrMod }
  from '@sagemath-ts/ntl-ts/src/ZZ_pX.js';
const F = new ZZ_pXModulus([1n, 0n, 1n], 7n);
const B = new ZZ_pXMultiplier([3n, 1n], F);
console.log(rem([0n, 0n, 0n, 1n], F)); // [0n, 6n]
console.log(MulMod([1n, 2n], B, F)); // [1n]
console.log(SqrMod([1n, 2n], F)); // [4n, 4n]
```

CRT input has one row per modulus-specific FFT prime and equal column counts.
Signed residues are normalized; empty columns are retained. Up to 800 primes,
reconstruction follows the bundled AArch64 fused floating estimate, including
its rounding at the midpoint even with one prime. Larger prime counts use
NTL's exact centered CRT tree. Calling CRT initializes the modulus's FFT state
before checking the row count.

<!-- llm-test: ntl-big-quotient-crt -->
```typescript
import { ZZ_pXModulus, FromModularRep } from '@sagemath-ts/ntl-ts/src/ZZ_pX.js';
const F = new ZZ_pXModulus(null, 6n);
console.log(FromModularRep([[441352763482308608n]], F)); // [3n]
```


Cached word and arbitrary-modulus multipliers retain their allocated FFT-prime count
across plain and failed rebuilds. A changed count raises `fftRep: inconsistent use`
(word) or `FFTRep: inconsistent use` (arbitrary modulus). Use a new multiplier object
when the new modulus requires a different count. Retained FFT data keeps the
coefficient modulus under which it was built during truncated inverse operations.

### NTL ordinary word products and extended GCD

The final options argument accepts `maxroot` and shared `state` as in
`zz_pXModulus`. Products preserve the native pre-product FFT size check, even
when leading coefficients cancel modulo a composite modulus. Extended GCD
returns copied coefficient arrays and retains the word half-GCD crossovers.

```ts
import { mul, sqr, InvMod, InvModStatus } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
mul([1n,2n], [3n,4n], 5n); // [3n,0n,3n]
sqr([1n,2n], 5n); // [1n,4n,4n]
InvMod([1n,1n], [1n,0n,1n], 5n, {maxroot: 0}); // [3n,2n]
InvModStatus([1n,1n], [1n,0n,1n], 5n, {maxroot: 0}); // [0,[3n,2n]]
```

```ts
import { XGCD } from '@sagemath-ts/ntl-ts/src/lzz_pX1.js';
XGCD([1n,0n,1n], [1n,1n], 5n); // [[1n],[3n],[3n,2n]]
XGCD([], [], 5n); // [[],[1n],[]]
```

### NTL Hensel lifting with shared FFT state

`MultiLift` passes the supplied state through its word factor tree, integer
products and arbitrary-modulus quotient operations. `AdditionalLifting` accepts
that state as its last argument and chooses its own degree-dependent word
`maxroot`, after normalizing the target's leading coefficient. Errors preserve
the FFT tables and random-stream effects of operations already completed.
Input arrays remain independent of the returned factors.

```ts
import { MultiLift } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const state = {context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32))};
const a = [1n,...Array<bigint>(21).fill(0n),1n];
const b = [2n,...Array<bigint>(21).fill(0n),1n];
const f = [2n,...Array<bigint>(21).fill(0n),3n,...Array<bigint>(21).fill(0n),1n];
MultiLift([a,b], f, 9, 17n, {state}); // [a,b], copied
state.context.length(); // 2
```

```ts
import { AdditionalLifting } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const state = {context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32))};
AdditionalLifting(5n, 1, [[1n,1n],[2n,1n]], 5n, 4, [7n,3n,1n], false, state);
// [625n,4,[[281n,1n],[347n,1n]]]
state.context.length(); // 1
```


### NTL FFT-prime coefficient contexts and integer GCD state

Set `fftPrime` to a cached FFT-prime index and pass its prime as p. This selects
native FFTInit behavior: `PrimeCnt=0`, `NumPrimes=1`, and a prime-derived
`MaxRoot`; `maxroot` is ignored. Operands and retained caches must use the same
coefficient context. A mismatched p or non-safe-integer index raises RangeError.
Without shared state this mode creates a private zero-key stream/cache.

```ts
import { zz_pXModulus, rem } from '@sagemath-ts/ntl-ts/src/lzz_pX.js';
const p = 882705526964617217n;
const F = new zz_pXModulus([1n,0n,1n], p, {fftPrime:0});
[F.PrimeCnt,F.NumPrimes,F.MaxRoot]; // [0,1,25]
rem([0n,0n,1n], F); // [882705526964617216n]
```

Integer `ZZX_GCD` and `ZZX_SquareFreeDecomp` accept shared state as their last
argument. Prime initialization, modular GCD and exact-division certification all
advance that state in native order. Without it, each call creates private state.

```ts
import { ZZX_GCD, ZZX_SquareFreeDecomp } from '@sagemath-ts/ntl-ts';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const state = {context:new FFTPrimeContext(), stream:new RandomStream(new Uint8Array(32))};
ZZX_GCD([-2n,0n,2n], [2n,-4n,2n], state); // [-2n,2n]
ZZX_SquareFreeDecomp([1n,4n,4n], state); // [[[1n,2n],2]]
state.context.length(); // 2
```


### NTL complete integer factorization

`ZZX_factor(f, bound = 0, options?)` and `ZZX_SFFactor(f, bound = 0, options?)`
are exported from `@sagemath-ts/ntl-ts`. Their native names `factor`, `SFFactor`,
`ll_SFFactor` and `FindTrueFactors_vH` are available from
`@sagemath-ts/ntl-ts/src/ZZXFactoring.js`.

| Entry | Result and contract |
|---|---|
| `factor(f: readonly bigint[], bound?: number, options?: IntegerFactorizationOptions)` | `[bigint, Array<[bigint[], number]>]`: signed content and ordered primitive factors with multiplicities. Zero/constants return their constant coefficient and no factors. |
| `SFFactor(f, bound?, options?)` | `bigint[][]`; f must be primitive and squarefree with positive leading coefficient. Zero raises `SFFactor: bad args`. Includes native deflation/inflation. |
| `ll_SFFactor(f, bound?, options?)` | `bigint[][]`; the primitive squarefree driver before deflation. Same positive-leading-coefficient precondition. |
| `FindTrueFactors_vH(f, lifted, P: bigint, p: bigint, exponent: number, info: LocalInfoT, bound: number, options?)` | `bigint[][]`; adaptive van Hoeij recombination. Requires the retained local information and lifted factors produced by native prime selection, with P = p**exponent satisfying the driver's coefficient and pruning bounds. It can update info. |

A zero bound selects the native bound. `IntegerFactorizationOptions` accepts
`state?: PolynomialProductState`, `InitNumPrimes`, `MaxNumPrimes`, `MaxPrune`,
`van_Hoeij` and `PowerHack`, with native defaults 7, 50, 10, 1 and 1. All control
numbers must be safe integers. The internal `ok_to_abandon` flag allows an
intermediate deflation stage to return unresolved factors; SFFactor manages it.
Native factor order is retained. Arrays are copied; shared state retains FFT
cache and random-stream effects. Omitting state creates a private zero-key stream
and cache for the call.

```ts
import { ZZX_factor, ZZX_SFFactor } from '@sagemath-ts/ntl-ts';
import { FFTPrimeContext } from '@sagemath-ts/ntl-ts/src/FFT.js';
import { RandomStream } from '@sagemath-ts/ntl-ts/src/ZZ.js';
const state = {context:new FFTPrimeContext(), stream:new RandomStream(new Uint8Array(32))};
ZZX_factor([-2n,-4n,-2n], 0, {state}); // [-2n, [[[1n,1n],2]]]
state.context.length(); // 2
ZZX_SFFactor([-1n,0n,0n,0n,1n]); // [[1n,0n,1n],[1n,1n],[-1n,1n]]
```

```ts
import { ll_SFFactor, Compute_pb, Compute_pdelta, Compute_pb_eff }
  from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
ll_SFFactor([1n,3n,2n]); // [[1n,1n],[1n,2n]]
const p = 882705526964617217n;
Compute_pb([], [], p, 1, p*p, 4); // [[3], [p**3n]]
Compute_pdelta(0, 1n, p, 120); // [3, p**3n]
Compute_pb_eff(p, 1, p*p, 4, 120); // [5, p**5n]
```

`RationalPolynomial.isIrreducible()` delegates primitive integer factorization
to PARI for degrees 30 through 300 and NTL otherwise. Both true and false results
are cached. NumberField construction uses the same check.

```ts
import { RationalPolynomial, NumberField } from 'sagemath-ts/rings';
const f = RationalPolynomial.fromBigInts([-2n,...Array<bigint>(300).fill(0n),1n]);
f.isIrreducible(); // true
new NumberField(f, 'a').degree(); // 301
```

### PARI conjugate-count bound and word derivatives

`numberofconjugates(T, pinit = 2n): number` computes PARI's modular upper bound
on the automorphism count of the monic irreducible integer polynomial `T`.
The starting prime is a signed 64-bit argument represented by bigint; PARI's
prime iterator receives its unsigned word value. This is an upper bound, not
an assertion that the field is Galois.

`Flx_deriv(T, p): bigint[]` and `Flx_is_squarefree(T, p): boolean` require
coefficients reduced modulo the word prime `p`. Zero has squarefree predicate
false; nonzero constants have predicate true. Inputs are not mutated.

```ts
import { numberofconjugates, Flx_deriv, Flx_is_squarefree }
  from '@sagemath-ts/parigp-ts';
numberofconjugates([-2n,0n,0n,1n]); // 1
numberofconjugates([-2n,0n,0n,1n], -1n); // 3: unsigned prime range is exhausted
Flx_deriv([1n,0n,1n], 3n); // [0n,2n]
Flx_is_squarefree([1n,0n,1n], 2n); // false
Flx_is_squarefree([1n], 2n); // true
```


### Extension-field scalar square roots

`FiniteFieldElement.sqrt({extend?: boolean, all?: boolean})` defaults to
extend=false and all=false. It returns an element, or an array for all=true.
A runtime boolean all returns the union. Nonsquares return [] for all=true;
otherwise they raise ValueError("element is not a square"). extend=true raises
NotImplementedError. Roots use PARI's state-dependent choice and native order.

```ts
import { GFpn } from 'sagemath-ts/rings/finite_rings';
import { setrand } from '@sagemath-ts/parigp-ts';
const F = GFpn(3n, 2, [1, 0], 'a');
setrand(17n);
F.__call__(2n).sqrt({ all: true }).map(String); // ['a', '2*a']
```

Dependency export: `FF_issquareall(x: PariFfelt): PariFfelt | null` from
`@sagemath-ts/parigp-ts` is the root-output-pointer form, on reduced finite-field
records. It preserves PARI random state and returns polynomial coefficient
arrays. It does not expose the omitted-pointer predicate.


### Extension-field square predicates and norms in the dependency layer

`FiniteFieldElement.is_square(): boolean` delegates to PARI's deterministic
square predicate without advancing its random stream.

```ts
import { GFpn } from 'sagemath-ts/rings/finite_rings';
const F = GFpn(3n, 2, [1, 0], 'a');
F.gen().is_square(); // true
F.gen().add(F.one()).is_square(); // false
```

Exports from `@sagemath-ts/parigp-ts`:

- `FpX_resultant(x: bigint[], y: bigint[], p: bigint): bigint`
- `Flx_resultant(x: bigint[], y: bigint[], p: bigint): bigint`
- `FpXQ_norm(x: bigint[], T: bigint[], p: bigint): bigint`
- `Flxq_norm(x: bigint[], T: bigint[], p: bigint): bigint`
- `FpXQ_issquare(x: bigint[], T: bigint[], p: bigint): boolean`
- `Flxq_issquare(x: bigint[], T: bigint[], p: bigint): boolean`
- `FF_norm(x: PariFfelt): bigint`
- `FF_issquare(x: PariFfelt): boolean`

Arrays use ascending coefficients. Word inputs are reduced with a positive
word modulus. Field predicates require irreducible defining polynomials.
Resultants follow the native Euclidean/half-GCD thresholds. The internal source
export `polynomialResultant(x, y, arithmetic, p)` shares that recursion.
These dependency functions do not change the separate Sage element norm API.


### Curve coordinate lifting over extension fields

For `EllipticCurveGeneric`, `E.lift_x(x, true)` returns all points, sorted by
Sage's y-coordinate order; `E.lift_x(x)` returns the first or raises ValueError.
The overloads accept Integer wrappers, bigint, existing integral-number inputs,
and field elements. A coordinate already in the base field retains point type
`EllipticCurvePoint<F>`; an element of another field `G` gives
`EllipticCurvePoint<F | G>` because the curve can be promoted. Unknown inputs have
a checked overload returning `EllipticCurvePoint<FieldElement>` (or an array with
`all: true`). The third positional flag is `extend` (default false):
`E.lift_x(x, false, true)` returns one point after adjoining y if necessary;
`E.lift_x(x, true, true)` returns all such points. These overloads use the broader
`EllipticCurvePoint<FieldElement>` result because the coefficient field can change.
Existing base-field roots return immediately. Otherwise QQ produces a quadratic
number field, a prime field produces a quadratic finite field, and an explicit
finite extension produces a polynomial quotient over that field, with generator y.

`lift_x` checks canonical parent maps: prime-field coordinates embed in an
extension curve, and extension coordinates promote a prime-field curve even when
the coordinate is constant. Rational coordinates are rejected on finite-field
curves. Differently named explicit extension parents have no implicit canonical
embedding. Strings and null raise AttributeError because they have no parent.
`E.base_extend(R)` and its alias `change_ring(R)` convert whole coefficient
elements through R; an unchanged field returns E itself.

`E.is_x_coord(x: unknown): boolean` converts x through its base-field constructor,
accepting scalar strings, null and convertible field elements. Conversion
TypeError becomes `x must be coercible into the base ring of the curve`; other
exceptions propagate. It tests existence without extracting odd-characteristic
square roots. Real coordinates, other starting parent families, and global
field/curve caching remain under audit. These notes concern the generic curve class; the
separate optimized public `EllipticCurve` factory uses these same coordinate
callers for its prime-field short models; other factory behavior remains under audit.

```ts
import { GFpn } from 'sagemath-ts/rings/finite_rings';
import { EllipticCurveGeneric } from 'sagemath-ts/schemes/elliptic_curves';
const K = GFpn(3n, 2, [1, 0], 'a');
const E = new EllipticCurveGeneric(K, [K.zero(), K.zero(), K.zero(), K.one(), K.zero()]);
E.is_x_coord(K.one()); // true
E.is_x_coord('1'); // true
E.is_x_coord(null); // true
E.lift_x(K.one(), true).map(P => String(P.y())); // ['a', '2*a']
```


```ts
import { GF, GFpn } from 'sagemath-ts/rings/finite_rings';
import { QQ } from 'sagemath-ts/rings/rational_field';
import { EllipticCurveGeneric } from 'sagemath-ts/schemes/elliptic_curves';
const F = GF(3n);
const K = GFpn(3n, 2, [1, 0], 'a');
const E = new EllipticCurveGeneric(F, [F.zero(), F.zero(), F.zero(), F.one(), F.zero()]);
const P = E.lift_x(K.one());
String(P.x()); // '1'
String(P.y()); // 'a'
P.curve.base_ring === K; // true
P.curve === E; // false
E.lift_x(F.zero()).curve === E; // true
E.change_ring(F) === E; // true
E.lift_x(QQ.__call__(1n)); // throws TypeError
E.lift_x('1'); // throws AttributeError
```


```ts
import { GF } from 'sagemath-ts/rings/finite_rings';
import { EllipticCurveGeneric } from 'sagemath-ts/schemes/elliptic_curves';
const F = GF(3n);
const E = new EllipticCurveGeneric(F, [F.zero(), F.zero(), F.zero(), F.one(), F.one()]);
E.lift_x(F.__call__(2n), true).length; // 0
const points = E.lift_x(F.__call__(2n), true, true);
points.map(P => String(P.y())); // ['y', '2*y']
String(points[0].curve.base_ring); // 'Finite Field in y of size 3^2'
points[0].curve.is_on_curve(points[0].x(), points[0].y()); // true
```

`RationalPolynomial(coefficients: Rational[], variableName = 'x')` retains the
variable name in arithmetic, quotients/remainders, derivatives and defining-field
representations. The readonly `variableName` records it. This coefficient helper
expects arithmetic operands in the same variable. `NumberFieldElement.isZero()`
is the generic ring-element spelling of its existing Sage method `is_zero()`.

```ts
import { RationalPolynomial, NumberField, Rational } from 'sagemath-ts/rings';
const f = new RationalPolynomial([new Rational(-2n), Rational.zero(), Rational.one()], 'y');
String(f); // 'y^2 - 2'
String(f.derivative()); // '2*y'
f.scale(Rational.zero()).variableName; // 'y'
const K = new NumberField(f, 'a');
String(K); // 'Number Field in a with defining polynomial y^2 - 2'
K.zero().isZero(); // true
K.one().isZero(); // false
```

The optimized `EllipticCurveFiniteField.is_x_coord(x: unknown)` uses the same base
field conversion. Its `lift_x(x, all = false, extend = false)` uses the generic
coercion, promotion and root-selection logic. Points over the original field
retain the optimized point class and their `.x`/`.y` properties. A promoted or
newly extended point has the generic point class and `.x()`/`.y()` methods.
Typed base-field/integer inputs with extend false retain the optimized result
type; unknown inputs or extend true return a union of the two point types
(an array of that union for all true). General model conversion and the public
factory's extension-field input support remain under audit.

```ts
import { EllipticCurve, GF, QQ } from 'sagemath-ts';
const F = GF(3n);
const E = EllipticCurve(F, [1n, 0n]);
E.is_x_coord('0'); // true
E.lift_x(0n).curve === E; // true
E.lift_x(1n, true).length; // 0
E.lift_x(1n, true, true).map(String); // ['(1 : y : 1)', '(1 : 2*y : 1)']
E.lift_x(QQ.zero()); // throws TypeError
E.lift_x('0'); // throws AttributeError
```


Finite extension elements' `trace(): PrimeFieldElement` delegates to PARI's
`FF_trace` and returns an element of the existing prime-subfield parent. It does
not compute a sequence of Frobenius powers. The PARI package also exports:

- `FF_trace(x: PariFfelt): bigint`
- `FpXQ_trace(x: bigint[], T: bigint[], p: bigint): bigint`
- `Flxq_trace(x: bigint[], T: bigint[], p: bigint): bigint`
- `F2xq_trace(x: bigint, T: bigint): bigint`

Polynomial arrays use ascending, reduced coefficients. Binary polynomials use
bits. These are native-kernel interfaces: use a nonconstant separable modulus,
reduced elements and a positive prime; word kernels require p < 2^64. `FF_trace`
selects the native binary/word/large-prime kernel. The low-level FpXQ kernel keeps
PARI's derivative-degree indexing; it is not interchangeable with the word kernel
when the characteristic divides the modulus degree. Public finite-field trace
selects the appropriate kernel.

```ts
import { GFpn } from 'sagemath-ts/rings/finite_rings';
import { FF_trace, FpXQ_trace, Flxq_trace, F2xq_trace, PariType } from '@sagemath-ts/parigp-ts';
const F = GFpn(3n, 2, [1, 0], 'a');
F.gen().add(F.one()).trace().toString(); // '2'
F.gen().trace().parent === F.baseField; // true
FpXQ_trace([1n, 1n], [1n, 0n, 1n], 3n); // 2n
Flxq_trace([1n, 1n], [1n, 0n, 1n], 3n); // 2n
F2xq_trace(3n, 11n); // 1n
FF_trace({type: PariType.t_FFELT, p: 3n, degree: 2, value: [1n, 1n], definingPoly: [1n, 0n, 1n]}); // 2n
```


Finite extension elements expose `charpoly(varName = 'x'): Polynomial<PrimeFieldElement>`
through PARI `FF_charpoly`. `norm(): PrimeFieldElement` returns the signed constant
coefficient of that polynomial, retaining the existing prime-subfield parent.
The characteristic polynomial has the full extension degree even for subfield elements.

The PARI package exports `FF_charpoly(x: PariFfelt): bigint[]`,
`FpXQ_charpoly(x: bigint[], T: bigint[], p: bigint): bigint[]` and
`Flxq_charpoly(x: bigint[], T: bigint[], p: bigint): bigint[]`.
The quotient kernels accept reduced elements and a monic positive-degree modulus.
`FpX_FpXY_resultant(T: bigint[], Q: bigint[][], p: bigint): bigint[]` and
`Flx_FlxY_resultant(T: bigint[], Q: bigint[][], p: bigint): bigint[]` also accept
nonmonic T. Q is indexed by eliminated-variable degree, then retained-variable
degree; arrays are ascending. Use nonzero T and Q for the arbitrary-prime route.
Word kernels require positive prime p < 2^64 and reduced coefficients.
`FpV_polint(x: bigint[], y: bigint[], p: bigint): bigint[]` and
`Flv_polint(x: bigint[], y: bigint[], p: bigint): bigint[]` interpolate equally sized,
nonempty vectors at distinct reduced points. Repeated points and malformed native
inputs are outside these low-level contracts. Degree-one binary native charpoly
preserves PARI's Fl_inv error caused by repeated interpolation points.

```ts
import { GFpn } from 'sagemath-ts/rings/finite_rings';
import { FpX_FpXY_resultant, FpV_polint } from '@sagemath-ts/parigp-ts';
const F = GFpn(3n, 2, [1, 0], 'a');
F.gen().charpoly('y').toString(); // 'y^2 + 1'
F.gen().norm().toString(); // '1'
F.gen().norm().parent === F.baseField; // true
F.one().charpoly().toString(); // 'x^2 + x + 1'
FpX_FpXY_resultant([1n, 0n, 2n], [[1n, 1n], [0n, 1n]], 5n); // [2n, 4n, 3n]
FpV_polint([0n, 1n, 2n], [1n, 3n, 7n], 17n); // [1n, 1n, 1n]
```


Generic group `groups.linear_relation(P, Q, operation = '+', identity?, inverse?, op?,
options?: {ord_p?: IntegerLike; ord_q?: IntegerLike}): [bigint, bigint]` follows
Sage's divisor/BSGS search. Standard groups supply `additive_order()` or
`multiplicative_order()`; custom operations require explicit orders. It preserves
bundled Sage's `[order(Q), 0n]` early return when the two orders are coprime.
`IntegerMod.additive_order(): bigint` returns the modulus divided by its gcd with
the residue. Generic elliptic points' `additive_order(): bigint` aliases their
existing `order()` method, with the same supported domains.

`EllipticCurveGeneric._p_primary_torsion_basis` uses that group-relation dependency
and returns Sage's selected basis, including extension-field point ordering.
`division_points(P, m, true)` applies Sage's repeated-factor reduction for nonzero
2-torsion targets. Point-list results propagate an already-known target order;
finite-field scalar multiplication also propagates known orders. Generic finite
point multiplication delegates to PARI over supported prime fields of characteristic
> 3; characteristics 2/3 and extension-field backend routing remain open.
`division_polynomial` accepts Sage's special indices 0, -1 and -2.

```ts
import { Mod, groups } from 'sagemath-ts';
Mod(2n, 20n).additive_order(); // 10n
Mod(2n, 90308402384902n).additive_order(); // 45154201192451n
groups.linear_relation(Mod(2n, 20n), Mod(4n, 20n)); // [2n, 1n]
groups.linear_relation(Mod(0n, 20n), Mod(2n, 20n)); // [10n, 0n]
```


Prime short-model curves expose `E.abelian_group(): AbelianGroupStructure`, the
same cached `{invariants, generators, order}` record as `abelian_group(E)`.
Computing the group updates `E.gens()` to its corrected direct-product basis.
`E._points_via_group_structure()` enumerates multiples of that basis;
`E.points()` sorts and caches the result. The returned array is immutable:
assignment raises Sage's ValueError. Use `.slice()` for a mutable copy.

```ts
import { EllipticCurve, GF } from 'sagemath-ts';
const E = EllipticCurve(GF(5n), [1n, 3n]);
E.points().map(String); // ['(0 : 1 : 0)', '(1 : 0 : 1)', '(4 : 1 : 1)', '(4 : 4 : 1)']
E.points() === E.points(); // true
E.abelian_group().invariants; // [4n]
E.abelian_group() === E.abelian_group(); // true
```


The default `EllipticCurve` has a five-coefficient overload:
`EllipticCurve<F>(field: FieldParent, coeffs: [unknown, unknown, unknown, unknown, unknown]): EllipticCurveGeneric<F>`.
It preserves the supplied model in characteristics 2, 3 and larger primes, over
explicit finite extensions and QQ. These results use the existing generic API:
`base_ring`, `a_invariants()`, `point([x,y])`, and point `x()`, `y()`, `is_zero()`.
The two-coefficient prime-field overload retains `EllipticCurveFiniteField` and
its coordinate properties. Its coefficient types now include IntegerLike and
field elements. Both constructors coerce every coefficient into the base ring.

Generic finite curves still lack the optimized class's full PARI cardinality and
group backend. Use the two-coefficient prime path when those operations are needed;
converting models without transforming points is not a valid workaround. General
finite-field backend integration remains an open implementation task.

```ts
import { EllipticCurve, GF } from 'sagemath-ts';
const K = GF(2n);
const E = EllipticCurve(K, [1n, 0n, 0n, 0n, 1n]);
E.a_invariants().map(String); // ['1', '0', '0', '0', '1']
E.point([K.zero(), K.one()]).mul(2n).toString(); // '(0 : 1 : 0)'
```


`EllipticCurveGeneric.pari_curve()` and `.__pari__()` return the same cached
PARI general-model record over PrimeField, FiniteFieldPrime and prime IntegerModRing
parents with characteristic > 3. Characteristic-two/three prime parents,
GF2Field and explicit FiniteFieldExtension parents of any characteristic are also
supported for scalar multiplication. The return type is now
`EllipticCurve | FFEllipticCurve`; `'field' in model` distinguishes the latter.
Other parents currently raise NotImplementedError.
Generic point `mul(n: IntegerLike | number)` and `rmul(n: IntegerLike | number)`
use `ZZ` coercion and delegate to PARI on those supported parents. Scalar results
preserve the original equation and propagate any known order. On the prime
parents of characteristic > 3, point `order()` and `order({algorithm: 'pari'})` now delegate to PARI,
caching the curve cardinality and point order. `order({algorithm: 'hybrid'})`
follows Sage's bounded-search/partial-factorization schedule; `generic_small`
uses generic bounds. `additive_order(options?)` accepts the same options as
`order(options?)`. `_compute_order(algorithm?)` is the corresponding uncached
nonzero-point calculation. Generic prime curves also expose `cardinality()`,
`order()` and `cardinality_pari()`; the last does not set the curve's `_order` cache.
Additional cardinality algorithms, extension-degree options, group APIs and other
parent backends remain open.

PARI exports `ell_to_a4a6_bc(E, p): [bigint, bigint, [bigint,bigint,bigint,bigint]]`
for general models over p > 3, and `FpE_changepoint(P, ch, p)` /
`FpE_changepointinv(P, ch, p)` for reduced prime-field coordinates with invertible
`ch[0]`. The latter use the existing EllipticPoint union and readonly four-term
`ch = [u,r,s,t]`. `ellmul` now also accepts general `ellinit` records over p > 3;
its existing short-model record API remains supported.
`ellcard(E)` and `ellorder(E, P, multiple?)` also accept general ellinit records
over p > 3; `ellorder` accepts either existing point representation. Omitting
`multiple` uses the cached factored group exponent for general records.

```ts
import { EllipticCurve, GF, Integer } from 'sagemath-ts';
import { ellmul, ellcard, ellorder, ell_to_a4a6_bc } from '@sagemath-ts/parigp-ts';
const K = GF(11n);
const E = EllipticCurve(K, [1n, 0n, 1n, 0n, 0n]);
E.pari_curve() === E.__pari__(); // true
const model = E.pari_curve();
if ('field' in model) throw new Error('expected a prime-field model');
ell_to_a4a6_bc(model, 11n); // [5n, 6n, [6n, 3n, 3n, 9n]]
ellmul(model, {isInfinity: false, x: 0n, y: 0n}, 2n); // {isInfinity: false, x: 0n, y: 10n}
E.point([K.zero(), K.zero()]).mul(new Integer(2n)).toString(); // '(0 : 10 : 1)'
ellcard(model); // 6n
ellorder(model, {isInfinity: false, x: 0n, y: 0n}); // 3n
E.point([K.zero(), K.zero()]).order({algorithm: 'pari'}); // 3n
```


`Integer.factor(options?: {limit?: IntegerLike})` accepts a trial-division bound.
With a bound, the returned `[base, exponent][]` may retain a composite cofactor;
without one it retains the existing complete-factorization behavior.
`factor_trial_division(m: IntegerLike, limit?: IntegerLike)` is available from
`sagemath-ts/rings`. Its default limit is the signed 64-bit maximum. Negative
inputs retain the established `[-1n, 1n]` unit entry. The free `factor(n)` API
still takes no options.

```ts
import { Integer, EllipticCurve, GF } from 'sagemath-ts';
import { factor_trial_division } from 'sagemath-ts/rings';
new Integer(143n).factor({limit: 9n}); // [[143n, 1n]]
new Integer(143n).factor({limit: new Integer(12n)}); // [[11n, 1n], [13n, 1n]]
factor_trial_division(-143n, 12n); // [[-1n, 1n], [11n, 1n], [13n, 1n]]
const K = GF(11n);
const E = EllipticCurve(K, [0n, 0n, 0n, 1n, 1n]);
E.cardinality_pari(); // 14n
E.cardinality(); // 14n
E.order(); // 14n
E.point([K.zero(), K.one()]).order({algorithm: 'hybrid'}); // 7n
```


```ts
import { EllipticCurve, GF } from 'sagemath-ts';
const K = GF(11n);
const E = EllipticCurve(K, [1n, 1n]);
E.point(0n, 1n).order({algorithm: 'generic_small'}); // 7n
E.point(0n, 1n).additive_order({algorithm: 'hybrid'}); // 7n
E.point(0n, 1n).order({algorithm: 'pari'}); // 7n
E._order; // 14n
```


```ts
import { EllipticCurve, GF, Integer } from 'sagemath-ts';
const E = EllipticCurve(GF(11n), [1n, 1n]);
const P = E.point(0n, 1n);
P.order(); // 7n
P.mul(new Integer(2n)).toString(); // '(3 : 3 : 1)'
P.mul(2).toString(); // '(3 : 3 : 1)'
E.pari_curve() === E.__pari__(); // true
E.pari_curve() === E.toPari(); // true
```


PARI binary elliptic kernels use the existing `EllipticPoint` union with packed
nonnegative bigint coefficient bits. Over GF(2)[y]/T, `T` is also a bit polynomial
(e.g. `7n` is y²+y+1). Curve operations require reduced coordinates and a valid
field modulus. `F2xqECoefficient` is either ordinary `a2: bigint` for
Y²+XY=X³+a2X²+a6, or supersingular `[a3,a4,inverse(a3)]` for
Y²+a3Y=X³+a4X+a6. The scalar kernels recover a6 from the point and do not take it.

Exports from `@sagemath-ts/parigp-ts`:

- `F2xqE_add(P,Q,a,T)`, `F2xqE_sub(P,Q,a,T)`, `F2xqE_dbl(P,a,T)`, `F2xqE_neg(P,a,T)`.
- `F2xqE_mul(P,n,a,T)` for signed bigint scalars, using PARI's powering windows.
- `F2xqE_changepoint(P,ch,T)`, `F2xqE_changepointinv(P,ch,T)`, with
  `F2xqEChange = readonly [u,r,s,t]` and nonzero u.
- `F2xq_invsafe(x,T): bigint | null`, `F2xq_inv(x,T): bigint`, and
  `F2xq_div(x,y,T): bigint`. The latter two throw native PariError on nonunits.

General Sage binary-curve scalar multiplication now delegates through the model
adapters below. Binary point-order/cardinality/group backends remain open.

```ts
import { F2xq_inv, F2xq_invsafe, F2xqE_dbl, F2xqE_mul } from '@sagemath-ts/parigp-ts';
const P = {isInfinity: false as const, x: 1n, y: 2n};
F2xq_inv(2n, 7n); // 3n
F2xq_invsafe(3n, 5n); // null
F2xqE_dbl(P, 1n, 7n); // {isInfinity: false, x: 0n, y: 1n}
F2xqE_mul(P, 4n, 1n, 7n); // {isInfinity: true}
```


Finite-field model adapters are exported from `@sagemath-ts/parigp-ts`:

- `ellinit_Fq(coefficients, field): FFEllipticCurve | null` accepts `[j]`, `[a4, a6]`
  or `[a1, a2, a3, a4, a6]` (`FFEllipticScalar` entries), plus a `PariFfelt` field
  descriptor. A j-invariant is reduced into that field before selecting the model.
  It returns null for a singular model, matching the native internal initializer.
- `FF_ellinit(invariants, field): FFEllipticCurve` accepts the twelve precomputed
  `FFEllipticInvariants` entries (a1 through disc). It allows singular records and
  sets j to zero for them. It does not recompute the supplied invariants.
- `FF_ellmul(E, P, n)` and the `ellmul(E, P, n)` overload accept a finite-field curve
  record, `FFEllipticInputPoint`, and signed bigint scalar. They return
  `FFEllipticPoint`, with PariFfelt coordinates, in the original model.

`FFEllipticScalar` is `bigint | PariFfelt`. Integers are prime-field constants;
use a PariFfelt coefficient array for an extension generator. `FFEllipticCurve`
is now `BinaryFFEllipticCurve | OddFFEllipticCurve`. Both store field-valued
invariants and `field`; narrow with `'binaryModel' in E` to access
`binaryModel: readonly [F2xqECoefficient, bigint, F2xqEChange]`. Otherwise use
`oddModel: readonly [FlxqECoefficient, bigint[], FqEllipticChange]`.
Only binaryModel uses packed bits. Inputs must belong to the supplied valid
field. `ellinit(coefficients, 2n/3n)` returns the same finite-field record type,
or null for a singular curve. With a general bigint domain, its return type is
`EllipticCurve | FFEllipticCurve | null`; narrow before calling prime-only APIs.
Without a domain it still returns the integer record. For other domains, singular
initialization retains the existing exception behavior (an open fidelity gap).
`ellj`, `elldisc` and `ellcoeffs` accept both record families and preserve their
scalar types. `ellisnonsingular` checks the field-valued discriminant correctly.

```ts
import { PariType, ellinit_Fq } from '@sagemath-ts/parigp-ts';
const field = {type: PariType.t_FFELT as const, p: 3n, degree: 2,
  definingPoly: [1n, 0n, 1n], value: [0n, 1n]};
const E = ellinit_Fq([field], field)!; // j is the extension generator
E.a2.value; // [0n, 1n]
E.a6.value; // [1n]
E.j.value;  // [0n, 1n]
ellinit_Fq([1n, 1n], field)!.a4.value; // [1n]
ellinit_Fq([0n, 0n], field); // null (singular short model)
```

```ts
import { PariType, ellinit_Fq, ellmul, FF_ellmul } from '@sagemath-ts/parigp-ts';
const field = {type: PariType.t_FFELT as const, p: 2n, degree: 2,
  definingPoly: [1n, 1n, 1n], value: [0n, 1n]};
const E = ellinit_Fq([1n, 1n, 0n, 0n, 1n], field)!;
const P = {isInfinity: false as const, x: 1n, y: field};
const Q = ellmul(E, P, 2n);
if (!Q.isInfinity) [Q.x.value, Q.y.value]; // [[], [1n]]
FF_ellmul(E, P, 4n); // {isInfinity: true}
ellinit_Fq([0n, 0n, 0n, 1n, 1n], field); // null
```

```ts
import { EllipticCurve } from 'sagemath-ts';
import { FiniteFieldExtension } from 'sagemath-ts/rings/finite_rings';
// This constructor's modulus array omits the monic leading coefficient.
const K = new FiniteFieldExtension(2n, 2, [1, 1], 'a');
const E = EllipticCurve(K, [1n, 1n, 0n, 0n, 1n]);
const P = E.point([K.one(), K.gen()]);
P.mul(2n).toString(); // '(0 : 1 : 1)'
P.mul(4n).is_zero(); // true
E.pari_curve() === E.__pari__(); // true
```


Odd-extension elliptic kernels are exported from `@sagemath-ts/parigp-ts`.
`FqEllipticPoint` is `{isInfinity: true}` or
`{isInfinity: false, x: bigint[], y: bigint[]}` with reduced ascending polynomial
coefficients; zero is `[]`. `FqEllipticChange` is the readonly tuple `[u,r,s,t]`
of four such polynomials, with nonzero u. T is an irreducible polynomial over p.

The following names are available with both `FlxqE_` and `FpXQE_` prefixes:

- `add(P,Q,a,T,p)`, `sub(P,Q,a,T,p)`, `dbl(P,a,T,p)`, `neg(P,T,p)`.
- `mul(P,n,a,T,p)` for signed bigint n, using native gen_pow_i windows.
- `changepoint(P,ch,T,p)` and `changepointinv(P,ch,T,p)`.

FlxqE uses the word polynomial backend (odd prime p < 2^64); FpXQE uses the
arbitrary-prime polynomial backend. The coefficient a is polynomial a4 for
Y²=X³+a4X+a6. FlxqE also accepts `FlxqECoefficient = bigint[] | readonly [bigint[]]`:
in characteristic three, `[a2]` denotes the ordinary model Y²=X³+a2X²+a6.
These kernels recover a6 from a point. Native FF initialization selects the
characteristic-three model and FF scalar multiplication restores original
coordinates. Generic Sage point orders/cardinality/groups over these additional
parents still need their native backend integration.

```ts
import { FlxqE_dbl, FlxqE_mul, FpXQE_mul } from '@sagemath-ts/parigp-ts';
const T = [1n, 0n, 1n];
const P = {isInfinity: false as const, x: [1n], y: [1n]};
FlxqE_dbl(P, [[1n]], T, 3n); // {isInfinity: false, x: [1n], y: [2n]}
FlxqE_mul(P, 3n, [[1n]], T, 3n); // {isInfinity: true}
FpXQE_mul({isInfinity: false, x: [1n], y: [2n]}, 2n, [2n], T, 7n);
// {isInfinity: false, x: [], y: [1n]}
```

```ts
import { EllipticCurve, GF } from 'sagemath-ts';
import { FiniteFieldExtension } from 'sagemath-ts/rings/finite_rings';
const K = GF(3n);
const E = EllipticCurve(K, [1n, 0n, 0n, 1n, 1n]);
E.point([K.zero(), K.one()]).mul(2n).toString(); // '(0 : 2 : 1)'
const L = new FiniteFieldExtension(3n, 2, [1, 0], 'a');
const C = EllipticCurve(L, [0n, 0n, 0n, 1n, 1n]);
const P = C.point([L.gen(), L.one()]);
P.mul(2n).toString(); // '(a + 1 : 0 : 1)'
P.mul(4n).is_zero(); // true
C.pari_curve() === C.__pari__(); // true
```


```ts
import { ellinit, ellj, elldisc, ellcoeffs, ellisnonsingular, ellmul } from '@sagemath-ts/parigp-ts';
const E = ellinit([1n, 0n, 0n, 1n, 1n], 3n)!;
ellj(E).value; // [2n]
elldisc(E).value; // [2n]
ellcoeffs(E)[0].value; // [1n]
ellisnonsingular(E); // true
const Q = ellmul(E, {isInfinity: false, x: 0n, y: 1n}, 2n);
if (!Q.isInfinity) Q.y.value; // [2n]
ellinit([0n, 0n], 3n); // null
```


Finite-field point orders with a known annihilating multiple:

- `FF_ellorder(E, P, order): bigint` accepts an `FFEllipticCurve` and
  `FFEllipticInputPoint`, converts coordinates to the cached native model and
  delegates to the relevant extension kernel.
- `F2xqE_order(P, order, a, T)`, `FlxqE_order(P, order, a, T, p)` and
  `FpXQE_order(P, order, a4, T, p)` use the same point/coefficient encodings as
  their scalar-multiplication counterparts.
- `GroupOrder` is a positive bigint, factor rows `[[prime, exponent], ...]`, or
  `[positiveBigint, factorRows]`. The bound must annihilate the point; supplied
  factors must describe that bound. The routines do not verify those assumptions.
- `gen_order<T>(a, order, power, isIdentity): bigint` exposes the shared native
  recursion. It also accepts null to report a missing-order error. Identity and
  power callbacks use the represented group's operations.

These supplied-bound APIs do not yet provide automatic curve cardinality or
finite-field group structure, nor default general-field `ellorder` dispatch.

```ts
import { ellinit, FF_ellorder, gen_order } from '@sagemath-ts/parigp-ts';
const E = ellinit([1n, 0n, 0n, 1n, 1n], 3n)!;
const P = {isInfinity: false as const, x: 0n, y: 1n};
FF_ellorder(E, P, 12n); // 3n
FF_ellorder(E, P, [[2n, 2n], [3n, 1n]]); // 3n
FF_ellorder(E, P, [12n, [[2n, 2n], [3n, 1n]]]); // 3n
gen_order(4n, 12n, (x, n) => x * n % 12n, x => x === 0n); // 3n
```


`elltrace_extension(t: bigint, n: number, q: bigint): bigint` computes the
extension trace by powering modulo `X² - tX + q`. The nonnegative degree `n`
uses the existing native unsigned-word powering adapter; degree zero returns 2.
`Fp_ffellcard(a4, a6, q, n, p): bigint` counts a curve defined over the prime
field over an extension, normally with `q = p ** BigInt(n)`. It counts over the
base field and extends the trace. This helper does not yet provide the general
extension-field cardinality dispatcher, and inherits the existing prime-counter
algorithm limitations documented in DEVIATIONS.md.

```ts
import { elltrace_extension, Fp_ffellcard } from '@sagemath-ts/parigp-ts';
elltrace_extension(3n, 2, 7n); // -5n
elltrace_extension(3n, 0, 7n); // 2n
Fp_ffellcard(1n, 1n, 49n, 2, 7n); // 55n
```


PARI word-field point kernels are available from `@sagemath-ts/parigp-ts`:
`Flj_dbl_pre(P,a4,p,pi)`, `Flj_add_pre(P,Q,a4,p,pi)`, `Flj_neg(P,p)`,
`Flj_mulu_pre(P,n,a4,p,pi)`, `Fle_to_Flj(P)`, `Flj_to_Fle_pre(P,p,pi)`,
`Fle_dbl(P,a4,p)`, `Fle_add(P,Q,a4,p)`, `Fle_mulu(P,n,a4,p)` and
`Fle_order(P,order,a4,p)`. Jacobian points use `{X,Y,Z}`; affine points use
`{isInfinity:false,x,y}` or `{isInfinity:true}`. Residues/scalars are bigint,
`n` is unsigned 64-bit, and `order` is a supplied `GroupOrder`.
`pi` preserves the native signature but is unused by exact BigInt reduction.
`Fl_ellcard_Shanks(a4,a6,p): bigint` is a module-level export of
`parigp-ts/src/elliptic/group.ts` for nonsingular curves with
`99 < p < 2^63 - 2^32`; normal callers use `ellcard`.

```typescript
import { Fle_mulu, Fle_order, Flj_dbl_pre } from '@sagemath-ts/parigp-ts';
const P = {isInfinity: false as const, x: 0n, y: 1n};
Fle_order(P, 5n, 1n, 7n); // => 5n
Fle_mulu(P, 5n, 1n, 7n).isInfinity; // => true
Flj_dbl_pre({X: 3n, Y: 4n, Z: 0n}, 1n, 7n, 0n).X; // => 3n
```
