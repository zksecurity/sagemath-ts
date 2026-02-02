/**
 * Property-based tests for arithmetic functions using fast-check.
 *
 * These tests verify mathematical invariants that must hold for all inputs,
 * not just specific test vectors. When a counterexample is found, it's
 * automatically saved as a regression test.
 *
 * Run with: bun test tests/property/typescript/arith.property.test.ts
 */

import { beforeAll, describe, expect, it } from 'bun:test';
import * as fc from 'fast-check';
import {
  type Factorization,
  crt,
  divisors,
  euler_phi,
  factor,
  gcd,
  inverse_mod,
  is_prime,
  is_prime_power,
  is_square,
  is_squarefree,
  isqrt,
  kronecker_symbol,
  lcm,
  moebius,
  next_prime,
  power_mod,
  previous_prime,
  radical,
  sigma,
  xgcd,
} from '../../../packages/sagemath-ts/src/index.js';
import {
  anyBigInt,
  assertProperty,
  coprimePair,
  defaultFcParams,
  mediumPrime,
  nonNegativeBigInt,
  positiveBigInt,
  smallPrime,
} from './utils/property-test-utils.js';

// ============================================
// GCD Properties
// ============================================

describe('gcd properties', () => {
  it('gcd(a, b) === gcd(b, a) [commutativity]', () => {
    assertProperty(
      'arith',
      'gcd commutativity',
      fc.tuple(nonNegativeBigInt(), nonNegativeBigInt()),
      ([a, b]) => gcd(a, b) === gcd(b, a),
      defaultFcParams
    );
  });

  it('gcd(a, gcd(b, c)) === gcd(gcd(a, b), c) [associativity]', () => {
    assertProperty(
      'arith',
      'gcd associativity',
      fc.tuple(positiveBigInt(10n ** 10n), positiveBigInt(10n ** 10n), positiveBigInt(10n ** 10n)),
      ([a, b, c]) => gcd(a, gcd(b, c)) === gcd(gcd(a, b), c),
      defaultFcParams
    );
  });

  it('gcd(a, b) divides both a and b', () => {
    assertProperty(
      'arith',
      'gcd divides inputs',
      fc.tuple(positiveBigInt(), positiveBigInt()),
      ([a, b]) => {
        const g = gcd(a, b);
        return a % g === 0n && b % g === 0n;
      },
      defaultFcParams
    );
  });

  it('gcd(a, 0) === |a| [identity]', () => {
    assertProperty(
      'arith',
      'gcd identity with zero',
      anyBigInt(),
      (a) => gcd(a, 0n) === (a < 0n ? -a : a),
      defaultFcParams
    );
  });

  it('gcd(a*k, b*k) === |k| * gcd(a, b) [homogeneity]', () => {
    assertProperty(
      'arith',
      'gcd homogeneity',
      fc.tuple(positiveBigInt(10n ** 8n), positiveBigInt(10n ** 8n), positiveBigInt(10n ** 4n)),
      ([a, b, k]) => gcd(a * k, b * k) === k * gcd(a, b),
      defaultFcParams
    );
  });
});

// ============================================
// LCM Properties
// ============================================

describe('lcm properties', () => {
  it('lcm(a, b) === lcm(b, a) [commutativity]', () => {
    assertProperty(
      'arith',
      'lcm commutativity',
      fc.tuple(positiveBigInt(10n ** 8n), positiveBigInt(10n ** 8n)),
      ([a, b]) => lcm(a, b) === lcm(b, a),
      defaultFcParams
    );
  });

  it('a and b both divide lcm(a, b)', () => {
    assertProperty(
      'arith',
      'lcm is multiple',
      fc.tuple(positiveBigInt(10n ** 8n), positiveBigInt(10n ** 8n)),
      ([a, b]) => {
        const l = lcm(a, b);
        return l % a === 0n && l % b === 0n;
      },
      defaultFcParams
    );
  });

  it('gcd(a, b) * lcm(a, b) === a * b', () => {
    assertProperty(
      'arith',
      'gcd-lcm product identity',
      fc.tuple(positiveBigInt(10n ** 8n), positiveBigInt(10n ** 8n)),
      ([a, b]) => gcd(a, b) * lcm(a, b) === a * b,
      defaultFcParams
    );
  });
});

// ============================================
// Extended GCD Properties
// ============================================

describe('xgcd properties', () => {
  it('xgcd returns Bézout coefficients: g = s*a + t*b', () => {
    assertProperty(
      'arith',
      'xgcd Bézout identity',
      fc.tuple(positiveBigInt(10n ** 12n), positiveBigInt(10n ** 12n)),
      ([a, b]) => {
        const [g, s, t] = xgcd(a, b);
        return g === s * a + t * b;
      },
      defaultFcParams
    );
  });

  it('xgcd(a, b)[0] === gcd(a, b)', () => {
    assertProperty(
      'arith',
      'xgcd consistent with gcd',
      fc.tuple(nonNegativeBigInt(10n ** 12n), nonNegativeBigInt(10n ** 12n)),
      ([a, b]) => xgcd(a, b)[0] === gcd(a, b),
      defaultFcParams
    );
  });
});

// ============================================
// Factorization Properties
// ============================================

describe('factor properties', () => {
  it('product of factors equals original', () => {
    assertProperty(
      'arith',
      'factor product identity',
      fc.bigInt({ min: 2n, max: 10n ** 10n }),
      (n) => {
        const f = factor(n);
        const product = f.reduce((acc, [p, e]) => acc * p ** e, 1n);
        return product === n;
      },
      defaultFcParams
    );
  });

  it('all factors are prime', () => {
    assertProperty(
      'arith',
      'factor yields primes',
      fc.bigInt({ min: 2n, max: 10n ** 8n }),
      (n) => {
        const f = factor(n);
        return f.every(([p, _e]) => is_prime(p));
      },
      defaultFcParams
    );
  });

  it('exponents are positive', () => {
    assertProperty(
      'arith',
      'factor positive exponents',
      fc.bigInt({ min: 2n, max: 10n ** 8n }),
      (n) => {
        const f = factor(n);
        return f.every(([_p, e]) => e > 0n);
      },
      defaultFcParams
    );
  });

  it('factors are in ascending order', () => {
    assertProperty(
      'arith',
      'factor ascending order',
      fc.bigInt({ min: 2n, max: 10n ** 8n }),
      (n) => {
        const f = factor(n);
        for (let i = 1; i < f.length; i++) {
          if (f[i][0] <= f[i - 1][0]) return false;
        }
        return true;
      },
      defaultFcParams
    );
  });
});

// ============================================
// Primality Properties
// ============================================

describe('is_prime properties', () => {
  it('known small primes return true', () => {
    const primes = [
      2n,
      3n,
      5n,
      7n,
      11n,
      13n,
      17n,
      19n,
      23n,
      29n,
      31n,
      37n,
      41n,
      43n,
      47n,
      53n,
      59n,
      61n,
      67n,
      71n,
      73n,
      79n,
      83n,
      89n,
      97n,
    ];
    for (const p of primes) {
      expect(is_prime(p)).toBe(true);
    }
  });

  it('known small composites return false', () => {
    const composites = [
      4n,
      6n,
      8n,
      9n,
      10n,
      12n,
      14n,
      15n,
      16n,
      18n,
      20n,
      21n,
      22n,
      24n,
      25n,
      26n,
      27n,
      28n,
    ];
    for (const n of composites) {
      expect(is_prime(n)).toBe(false);
    }
  });

  it('0 and 1 are not prime', () => {
    expect(is_prime(0n)).toBe(false);
    expect(is_prime(1n)).toBe(false);
  });

  it('product of two distinct primes is composite', () => {
    const primes = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n];
    assertProperty(
      'arith',
      'prime product is composite',
      fc.tuple(
        fc.integer({ min: 0, max: primes.length - 1 }),
        fc.integer({ min: 0, max: primes.length - 1 })
      ),
      ([i, j]) => {
        if (i === j) return true; // Skip p*p
        return !is_prime(primes[i] * primes[j]);
      },
      defaultFcParams
    );
  });

  it('is_prime agrees with trial division for small n', () => {
    assertProperty(
      'arith',
      'is_prime correctness',
      fc.bigInt({ min: 2n, max: 10000n }),
      (n) => {
        // Trial division primality check
        if (n < 2n) return !is_prime(n);
        for (let d = 2n; d * d <= n; d++) {
          if (n % d === 0n) return !is_prime(n);
        }
        return is_prime(n);
      },
      { ...defaultFcParams, numRuns: 200 }
    );
  });
});

// ============================================
// Square Root Properties
// ============================================

describe('isqrt properties', () => {
  it('isqrt(n)² ≤ n < (isqrt(n)+1)²', () => {
    assertProperty(
      'arith',
      'isqrt bounds',
      nonNegativeBigInt(10n ** 30n),
      (n) => {
        const r = isqrt(n);
        return r * r <= n && n < (r + 1n) * (r + 1n);
      },
      defaultFcParams
    );
  });

  it('isqrt(n²) === n for perfect squares', () => {
    assertProperty(
      'arith',
      'isqrt perfect squares',
      nonNegativeBigInt(10n ** 15n),
      (n) => isqrt(n * n) === n,
      defaultFcParams
    );
  });

  it('isqrt is monotonic: a <= b implies isqrt(a) <= isqrt(b)', () => {
    assertProperty(
      'arith',
      'isqrt monotonic',
      fc.tuple(nonNegativeBigInt(10n ** 20n), nonNegativeBigInt(10n ** 20n)),
      ([a, b]) => {
        // Only test the implication: if a <= b, then isqrt(a) <= isqrt(b)
        if (a <= b) {
          return isqrt(a) <= isqrt(b);
        }
        // If a > b, test the other direction
        return isqrt(b) <= isqrt(a);
      },
      defaultFcParams
    );
  });
});

// ============================================
// Euler Phi Properties
// ============================================

describe('euler_phi properties', () => {
  it('euler_phi(p) === p - 1 for prime p', () => {
    assertProperty(
      'arith',
      'euler_phi of prime',
      mediumPrime(),
      (p) => euler_phi(p) === p - 1n,
      defaultFcParams
    );
  });

  it('euler_phi(p^k) === p^(k-1) * (p-1) for prime power', () => {
    assertProperty(
      'arith',
      'euler_phi of prime power',
      fc.tuple(smallPrime(), fc.integer({ min: 1, max: 5 })),
      ([p, k]) => {
        const pk = p ** BigInt(k);
        const expected = p ** BigInt(k - 1) * (p - 1n);
        return euler_phi(pk) === expected;
      },
      defaultFcParams
    );
  });

  it('euler_phi is multiplicative for coprime inputs', () => {
    assertProperty(
      'arith',
      'euler_phi multiplicative',
      coprimePair(10n ** 6n),
      ([a, b]) => euler_phi(a * b) === euler_phi(a) * euler_phi(b),
      defaultFcParams
    );
  });
});

// ============================================
// Modular Arithmetic Properties
// ============================================

describe('power_mod properties', () => {
  it('power_mod(a, 1, m) === a mod m', () => {
    assertProperty(
      'arith',
      'power_mod base case',
      fc.tuple(fc.bigInt({ min: 0n, max: 10n ** 10n }), fc.bigInt({ min: 2n, max: 10n ** 6n })),
      ([a, m]) => power_mod(a, 1n, m) === ((a % m) + m) % m,
      defaultFcParams
    );
  });

  it('power_mod(a, 0, m) === 1 (for m > 1)', () => {
    assertProperty(
      'arith',
      'power_mod zero exponent',
      fc.tuple(fc.bigInt({ min: 1n, max: 10n ** 10n }), fc.bigInt({ min: 2n, max: 10n ** 6n })),
      ([a, m]) => power_mod(a, 0n, m) === 1n,
      defaultFcParams
    );
  });

  it('power_mod(a, n+m, mod) === power_mod(a, n, mod) * power_mod(a, m, mod) mod mod', () => {
    assertProperty(
      'arith',
      'power_mod additive exponents',
      fc.tuple(
        fc.bigInt({ min: 2n, max: 100n }),
        fc.bigInt({ min: 1n, max: 50n }),
        fc.bigInt({ min: 1n, max: 50n }),
        fc.bigInt({ min: 2n, max: 10n ** 6n })
      ),
      ([a, n, m, mod]) => {
        const left = power_mod(a, n + m, mod);
        const right = (power_mod(a, n, mod) * power_mod(a, m, mod)) % mod;
        return left === right;
      },
      defaultFcParams
    );
  });
});

describe('inverse_mod properties', () => {
  it('a * inverse_mod(a, m) ≡ 1 (mod m) when gcd(a, m) = 1', () => {
    assertProperty(
      'arith',
      'inverse_mod correctness',
      fc.tuple(fc.bigInt({ min: 1n, max: 10n ** 8n }), mediumPrime()),
      ([a, p]) => {
        const aModP = ((a % p) + p) % p;
        if (aModP === 0n) return true; // Skip zero
        const inv = inverse_mod(aModP, p);
        return (aModP * inv) % p === 1n;
      },
      defaultFcParams
    );
  });
});

// ============================================
// CRT Properties
// ============================================

describe('crt properties', () => {
  it('crt(a, b, m, n) ≡ a (mod m) and ≡ b (mod n)', () => {
    assertProperty(
      'arith',
      'crt solution correctness',
      fc
        .tuple(
          fc.bigInt({ min: 0n, max: 1000n }),
          fc.bigInt({ min: 0n, max: 1000n }),
          smallPrime(),
          smallPrime()
        )
        .filter(([_a, _b, m, n]) => m !== n), // Require coprime moduli
      ([a, b, m, n]) => {
        const x = crt(a % m, b % n, m, n);
        return x % m === a % m && x % n === b % n;
      },
      defaultFcParams
    );
  });
});

// ============================================
// Divisor Function Properties
// ============================================

describe('divisors properties', () => {
  it('divisors(n) includes 1 and n', () => {
    assertProperty(
      'arith',
      'divisors contains 1 and n',
      positiveBigInt(10n ** 8n),
      (n) => {
        const divs = divisors(n);
        return divs.includes(1n) && divs.includes(n);
      },
      defaultFcParams
    );
  });

  it('all divisors actually divide n', () => {
    assertProperty(
      'arith',
      'divisors correctness',
      positiveBigInt(10n ** 6n),
      (n) => {
        const divs = divisors(n);
        return divs.every((d) => n % d === 0n);
      },
      defaultFcParams
    );
  });

  it('sigma(n, 0) === number of divisors', () => {
    assertProperty(
      'arith',
      'sigma(n, 0) counts divisors',
      positiveBigInt(10n ** 6n),
      (n) => sigma(n, 0n) === BigInt(divisors(n).length),
      defaultFcParams
    );
  });

  it('sigma(n, 1) === sum of divisors', () => {
    assertProperty(
      'arith',
      'sigma(n, 1) sums divisors',
      positiveBigInt(10n ** 6n),
      (n) => {
        const divs = divisors(n);
        const sum = divs.reduce((a, b) => a + b, 0n);
        return sigma(n, 1n) === sum;
      },
      defaultFcParams
    );
  });
});

// ============================================
// Prime Navigation Properties
// ============================================

describe('next_prime / previous_prime properties', () => {
  it('next_prime(n) > n', () => {
    assertProperty(
      'arith',
      'next_prime is greater',
      nonNegativeBigInt(10n ** 10n),
      (n) => next_prime(n) > n,
      defaultFcParams
    );
  });

  it('next_prime(n) is prime', () => {
    assertProperty(
      'arith',
      'next_prime yields prime',
      nonNegativeBigInt(10n ** 10n),
      (n) => is_prime(next_prime(n)),
      defaultFcParams
    );
  });

  it('no prime between n and next_prime(n)', () => {
    assertProperty(
      'arith',
      'next_prime is next',
      nonNegativeBigInt(10n ** 8n),
      (n) => {
        const np = next_prime(n);
        for (let k = n + 1n; k < np; k++) {
          if (is_prime(k)) return false;
        }
        return true;
      },
      { ...defaultFcParams, numRuns: 50 } // Slower test
    );
  });

  it('previous_prime(next_prime(p)) === p for prime p > 2', () => {
    assertProperty(
      'arith',
      'previous_prime inverse of next_prime',
      mediumPrime().filter((p) => p > 2n),
      (p) => previous_prime(next_prime(p)) === p,
      defaultFcParams
    );
  });
});

// ============================================
// Radical Properties
// ============================================

describe('radical properties', () => {
  it('radical(n) divides n', () => {
    assertProperty(
      'arith',
      'radical divides n',
      positiveBigInt(10n ** 8n),
      (n) => n % radical(n) === 0n,
      defaultFcParams
    );
  });

  it('radical(n) is squarefree', () => {
    assertProperty(
      'arith',
      'radical is squarefree',
      positiveBigInt(10n ** 8n),
      (n) => is_squarefree(radical(n)),
      defaultFcParams
    );
  });

  it('radical(p) === p for prime p', () => {
    assertProperty(
      'arith',
      'radical of prime',
      mediumPrime(),
      (p) => radical(p) === p,
      defaultFcParams
    );
  });
});
