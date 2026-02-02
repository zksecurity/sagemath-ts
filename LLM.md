# sagemath-ts LLM Reference

TypeScript port of SageMath number theory modules. Pure TypeScript, no native dependencies.

**Warning**: Educational only. Not constant-time, not audited.

## Installation

```bash
bun add sagemath-ts
# or
npm install sagemath-ts
```

## Quick Reference

### Arithmetic

```typescript
import { gcd, lcm, xgcd, factor, is_prime, power_mod, inverse_mod, euler_phi, crt } from 'sagemath-ts';

gcd(12n, 8n)                    // 4n
lcm(12n, 8n)                    // 24n
xgcd(15n, 6n)                   // { gcd: 3n, x: 1n, y: -2n }  (ax + by = gcd)
factor(60n)                     // [[2n, 2], [3n, 1], [5n, 1]]
is_prime(97n)                   // true
power_mod(2n, 100n, 1000000007n) // 976371285n
inverse_mod(3n, 7n)             // 5n (3*5 ≡ 1 mod 7)
euler_phi(60n)                  // 16n
crt([2n, 3n], [5n, 7n])         // 17n (x ≡ 2 mod 5, x ≡ 3 mod 7)
```

### Primes

```typescript
import { next_prime, previous_prime, prime_range, prime_factors, is_prime_power } from 'sagemath-ts';

next_prime(10n)                 // 11n
previous_prime(10n)             // 7n
prime_range(10n, 30n)           // [11n, 13n, 17n, 19n, 23n, 29n]
prime_factors(60n)              // [2n, 3n, 5n]
is_prime_power(8n)              // { base: 2n, exp: 3 }
```

### Divisibility

```typescript
import { divisors, sigma, valuation, is_square, isqrt, squarefree_part } from 'sagemath-ts';

divisors(12n)                   // [1n, 2n, 3n, 4n, 6n, 12n]
sigma(12n, 0)                   // 6n (number of divisors)
sigma(12n, 1)                   // 28n (sum of divisors)
valuation(12n, 2n)              // 2 (12 = 2^2 * 3)
is_square(16n)                  // true
isqrt(17n)                      // 4n
squarefree_part(12n)            // 3n
```

### Rings and Fields

```typescript
import { ZZ, QQ, Rational, Zmod, GF, Mod } from 'sagemath-ts';

// Integers
const n = ZZ(42n);

// Rationals
const q = new Rational(3n, 4n);
q.add(new Rational(1n, 2n))     // 5/4

// Modular integers (Z/nZ)
const Z7 = Zmod(7n);
const a = Z7(3n);
a.mul(Z7(5n))                   // 1 (mod 7)
Mod(10n, 7n)                    // 3

// Finite fields (GF(p))
const F13 = GF(13n);
const x = F13(5n);
x.pow(12n)                      // 1 (Fermat's little theorem)
```

### Matrices

```typescript
import { Matrix, ZZ, identity_matrix, zero_matrix } from 'sagemath-ts';

const A = Matrix(ZZ, [[1n, 2n], [3n, 4n]]);
const B = Matrix(ZZ, [[5n, 6n], [7n, 8n]]);

A.add(B)                        // element-wise addition
A.mul(B)                        // matrix multiplication
A.det()                         // determinant: -2n
A.transpose()
A.inverse()                     // over rationals
identity_matrix(ZZ, 3)          // 3x3 identity
zero_matrix(ZZ, 2, 3)           // 2x3 zeros
```

### Elliptic Curves

```typescript
import { EllipticCurve, GF } from 'sagemath-ts';

// Over finite field GF(p)
const F = GF(101n);
const E = EllipticCurve(F, [1n, 2n]);  // y² = x³ + x + 2

const P = E.point([3n, 6n]);
const Q = E.point([10n, 15n]);

P.add(Q)                        // point addition
P.mul(5n)                       // scalar multiplication
E.order()                       // curve order
E.is_on_curve([3n, 6n])         // true
E.random_point()
```

### Pairings

```typescript
import { weil_pairing, tate_pairing, ate_pairing, embedding_degree } from 'sagemath-ts';

weil_pairing(P, Q, n)           // Weil pairing e(P, Q)
tate_pairing(P, Q, n, k)        // Tate pairing
embedding_degree(E, n)          // smallest k where n | p^k - 1
```

### Lattices and LLL

```typescript
import { IntegerLattice, lllReduce, gramSchmidt, vector } from 'sagemath-ts';

const basis = [
  vector([1n, 0n, 0n]),
  vector([0n, 1n, 0n]),
  vector([0n, 0n, 1n])
];

const L = IntegerLattice(basis);
const reduced = lllReduce(basis);  // LLL-reduced basis
```

### LWE Cryptography

```typescript
import { LWE, Regev, DiscreteGaussianInteger } from 'sagemath-ts';

// Standard LWE
const lwe = new LWE({
  n: 256,
  q: 40961n,
  D: new DiscreteGaussianInteger(3.2)
});

const sample = lwe.sample();    // { a: [...], b: bigint }

// Regev encryption
const regev = new Regev({ n: 256, q: 40961n });
const { pk, sk } = regev.keygen();
const ct = regev.encrypt(pk, 1n);
const pt = regev.decrypt(sk, ct);
```

### Discrete Gaussian

```typescript
import { DiscreteGaussianInteger } from 'sagemath-ts';

const D = new DiscreteGaussianInteger(3.2);  // sigma = 3.2
D.sample()                      // sample from distribution
```

### Coding Theory

```typescript
import { ReedSolomonCode, createClassicalReedSolomonCode, GF } from 'sagemath-ts';

const F = GF(929n);
const RS = createClassicalReedSolomonCode(F, 10, 4);  // n=10, k=4

RS.encode([1n, 2n, 3n, 4n])     // encode message
RS.decode(codeword)             // decode/correct errors
RS.minimum_distance()           // d = n - k + 1 = 7
```

### Groups

```typescript
import { bsgs, pohlig_hellman, order_from_multiple } from 'sagemath-ts';

// Baby-step giant-step discrete log
bsgs(g, h, bounds, { op: mul, inverse: inv })  // find x: g^x = h

// Group order
order_from_multiple(g, multiple, { op: mul, identity: one })
```

## Type Coercion

Functions accept `IntegerLike = bigint | number | Integer`:

```typescript
gcd(12, 8)      // number
gcd(12n, 8n)    // bigint
gcd(ZZ(12), ZZ(8))  // Integer
```

Convert explicitly:

```typescript
import { toBigInt, toRational } from 'sagemath-ts';

toBigInt(42)              // 42n
toBigInt('123')           // 123n
toRational('3/4')         // Rational(3n, 4n)
```

## Namespaced Imports

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

## Common Patterns

### Modular arithmetic chain
```typescript
const p = 1000000007n;
const result = power_mod(
  inverse_mod(3n, p),
  euler_phi(p) - 1n,
  p
);
```

### Curve point operations
```typescript
const E = EllipticCurve(GF(101n), [1n, 2n]);
const G = E.random_point();
const n = E.order();
const P = G.mul(n);  // should be identity
```

### CRT for large moduli
```typescript
const x = crt([a1, a2, a3], [m1, m2, m3]);
// x ≡ a1 (mod m1), x ≡ a2 (mod m2), x ≡ a3 (mod m3)
```

## Related Packages

- `@sagemath-ts/parigp-ts` - PARI/GP algorithms
- `@sagemath-ts/flint-ts` - FLINT library port
- `@sagemath-ts/ntl-ts` - NTL library port
- `@zksecurity/cheatsheets` - Crypto curve parameters
