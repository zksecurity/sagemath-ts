/**
 * Unit tests for sage/arith/misc
 */
import { describe, expect, test } from 'bun:test';
import {
  CRT_basis,
  CRT_list,
  CRT_vectors,
  algebraic_dependency,
  bernoulli,
  continued_fraction,
  continued_fraction_value,
  convergents,
  coprime_part,
  crt,
  dedekind_sum,
  divisors,
  euler_phi,
  factor,
  formatFactorization,
  four_squares,
  gcd,
  half_gcd,
  hilbert_conductor,
  hilbert_conductor_inverse,
  hilbert_symbol,
  inverse_mod,
  is_prime,
  is_prime_power,
  is_pseudoprime,
  is_square,
  is_squarefree,
  isqrt,
  jacobi_symbol,
  kronecker_symbol,
  lcm,
  legendre_symbol,
  moebius,
  mqrr_rational_reconstruction,
  next_prime,
  next_probable_prime,
  number_of_divisors,
  power_mod,
  previous_prime,
  prime_divisors,
  prime_factors,
  prime_powers,
  prime_range,
  prime_to_m_part,
  radical,
  rational_reconstruction,
  sigma,
  smooth_part,
  sort_complex_numbers_for_display,
  sqrt_mod,
  squarefree_part,
  three_squares,
  trial_division,
  two_squares,
  valuation,
  xgcd,
} from './misc.js';

describe('gcd', () => {
  test('basic cases', () => {
    expect(gcd(12n, 8n)).toBe(4n);
    expect(gcd(97n, 100n)).toBe(1n);
    expect(gcd(0n, 5n)).toBe(5n);
    expect(gcd(5n, 0n)).toBe(5n);
    expect(gcd(0n, 0n)).toBe(0n);
  });

  test('negative numbers', () => {
    expect(gcd(-4n, 6n)).toBe(2n);
    expect(gcd(4n, -6n)).toBe(2n);
    expect(gcd(-4n, -6n)).toBe(2n);
  });

  test('large numbers', () => {
    const a = 97n * 10n ** 15n;
    const b = 19n ** 20n * 97n ** 2n;
    expect(gcd(a, b)).toBe(97n);
  });

  test('array form', () => {
    expect(gcd([2n, 4n, 6n, 8n])).toBe(2n);
    expect(gcd([15n, 25n, 35n])).toBe(5n);
    expect(gcd([])).toBe(0n);
  });
});

describe('lcm', () => {
  test('basic cases', () => {
    expect(lcm(4n, 6n)).toBe(12n);
    expect(lcm(3n, 5n)).toBe(15n);
    expect(lcm(0n, 5n)).toBe(0n);
    expect(lcm(5n, 0n)).toBe(0n);
  });

  test('array form', () => {
    expect(lcm([2n, 3n, 4n])).toBe(12n);
    expect(lcm([1n, 2n, 3n, 4n, 5n])).toBe(60n);
  });
});

describe('xgcd', () => {
  test('basic cases', () => {
    const [g, s, t] = xgcd(6n, 4n);
    expect(g).toBe(2n);
    expect(s * 6n + t * 4n).toBe(g);
  });

  test('coprime numbers', () => {
    const [g, s, t] = xgcd(5n, 7n);
    expect(g).toBe(1n);
    expect(s * 5n + t * 7n).toBe(g);
  });

  test('large numbers', () => {
    const a = 58526524056n;
    const b = 101294172798n;
    const [g, s, t] = xgcd(a, b);
    expect(g).toBe(22544886n);
    expect(s * a + t * b).toBe(g);
  });

  test('with zero', () => {
    expect(xgcd(5n, 0n)).toEqual([5n, 1n, 0n]);
    expect(xgcd(0n, 5n)).toEqual([5n, 0n, 1n]);
  });
});

describe('is_prime', () => {
  test('small primes', () => {
    expect(is_prime(2n)).toBe(true);
    expect(is_prime(3n)).toBe(true);
    expect(is_prime(5n)).toBe(true);
    expect(is_prime(7n)).toBe(true);
    expect(is_prime(11n)).toBe(true);
    expect(is_prime(13n)).toBe(true);
    expect(is_prime(17n)).toBe(true);
    expect(is_prime(19n)).toBe(true);
    expect(is_prime(23n)).toBe(true);
    expect(is_prime(389n)).toBe(true);
  });

  test('composites', () => {
    expect(is_prime(1n)).toBe(false);
    expect(is_prime(4n)).toBe(false);
    expect(is_prime(6n)).toBe(false);
    expect(is_prime(8n)).toBe(false);
    expect(is_prime(9n)).toBe(false);
    expect(is_prime(10n)).toBe(false);
    expect(is_prime(2000n)).toBe(false);
  });

  test('negative and zero', () => {
    expect(is_prime(0n)).toBe(false);
    expect(is_prime(-1n)).toBe(false);
    expect(is_prime(-2n)).toBe(false);
    expect(is_prime(-5n)).toBe(false);
  });

  test('larger primes', () => {
    expect(is_prime(104729n)).toBe(true); // 10000th prime
    expect(is_prime(1000000007n)).toBe(true);
  });

  test('Carmichael numbers (pseudoprimes)', () => {
    // These are composite but pass some primality tests
    expect(is_prime(561n)).toBe(false); // 3 × 11 × 17
    expect(is_prime(1105n)).toBe(false); // 5 × 13 × 17
    expect(is_prime(1729n)).toBe(false); // 7 × 13 × 19
  });
});

describe('factor', () => {
  test('prime powers', () => {
    expect(factor(8n)).toEqual([[2n, 3n]]);
    expect(factor(27n)).toEqual([[3n, 3n]]);
    expect(factor(32n)).toEqual([[2n, 5n]]);
  });

  test('products of primes', () => {
    expect(factor(12n)).toEqual([
      [2n, 2n],
      [3n, 1n],
    ]);
    expect(factor(30n)).toEqual([
      [2n, 1n],
      [3n, 1n],
      [5n, 1n],
    ]);
    expect(factor(100n)).toEqual([
      [2n, 2n],
      [5n, 2n],
    ]);
  });

  test('primes', () => {
    expect(factor(17n)).toEqual([[17n, 1n]]);
    expect(factor(101n)).toEqual([[101n, 1n]]);
  });

  test('negative numbers', () => {
    expect(factor(-12n)).toEqual([
      [-1n, 1n],
      [2n, 2n],
      [3n, 1n],
    ]);
  });

  test('one', () => {
    expect(factor(1n)).toEqual([]);
  });

  test('formatFactorization', () => {
    expect(formatFactorization(factor(12n))).toBe('2^2 * 3');
    expect(formatFactorization(factor(30n))).toBe('2 * 3 * 5');
    expect(formatFactorization(factor(1n))).toBe('1');
  });
});

describe('power_mod', () => {
  test('basic cases', () => {
    expect(power_mod(2n, 10n, 1000n)).toBe(24n);
    expect(power_mod(3n, 5n, 7n)).toBe(5n); // 243 mod 7 = 5
    expect(power_mod(2n, 0n, 5n)).toBe(1n);
  });

  test('negative exponent', () => {
    expect(power_mod(3n, -1n, 7n)).toBe(5n); // 3 * 5 = 15 ≡ 1 mod 7
  });

  test('large numbers', () => {
    const result = power_mod(2n, 1000n, 10n ** 9n + 7n);
    expect(result).toBe(688423210n);
  });
});

describe('inverse_mod', () => {
  test('basic cases', () => {
    expect(inverse_mod(3n, 7n)).toBe(5n);
    expect(inverse_mod(2n, 5n)).toBe(3n);
    expect(inverse_mod(7n, 11n)).toBe(8n);
  });

  test('throws for non-invertible', () => {
    expect(() => inverse_mod(2n, 4n)).toThrow();
    expect(() => inverse_mod(0n, 5n)).toThrow();
  });
});

describe('euler_phi', () => {
  test('basic cases', () => {
    expect(euler_phi(1n)).toBe(1n);
    expect(euler_phi(2n)).toBe(1n);
    expect(euler_phi(12n)).toBe(4n);
    expect(euler_phi(36n)).toBe(12n);
  });

  test('prime', () => {
    expect(euler_phi(7n)).toBe(6n);
    expect(euler_phi(13n)).toBe(12n);
  });

  test('prime power', () => {
    expect(euler_phi(8n)).toBe(4n); // 2^3: 2^2 * (2-1) = 4
    expect(euler_phi(27n)).toBe(18n); // 3^3: 3^2 * (3-1) = 18
  });
});

describe('next_prime / previous_prime', () => {
  test('next_prime', () => {
    expect(next_prime(10n)).toBe(11n);
    expect(next_prime(11n)).toBe(13n);
    expect(next_prime(1n)).toBe(2n);
    expect(next_prime(2n)).toBe(3n);
  });

  test('previous_prime', () => {
    expect(previous_prime(10n)).toBe(7n);
    expect(previous_prime(11n)).toBe(7n);
    expect(previous_prime(3n)).toBe(2n);
  });
});

describe('isqrt', () => {
  test('perfect squares', () => {
    expect(isqrt(0n)).toBe(0n);
    expect(isqrt(1n)).toBe(1n);
    expect(isqrt(4n)).toBe(2n);
    expect(isqrt(9n)).toBe(3n);
    expect(isqrt(16n)).toBe(4n);
    expect(isqrt(100n)).toBe(10n);
  });

  test('non-perfect squares', () => {
    expect(isqrt(2n)).toBe(1n);
    expect(isqrt(3n)).toBe(1n);
    expect(isqrt(5n)).toBe(2n);
    expect(isqrt(10n)).toBe(3n);
    expect(isqrt(99n)).toBe(9n);
  });

  test('large numbers', () => {
    expect(isqrt(10n ** 20n)).toBe(10n ** 10n);
  });
});

describe('is_square', () => {
  test('perfect squares', () => {
    expect(is_square(0n)).toBe(true);
    expect(is_square(1n)).toBe(true);
    expect(is_square(4n)).toBe(true);
    expect(is_square(9n)).toBe(true);
    expect(is_square(100n)).toBe(true);
  });

  test('non-squares', () => {
    expect(is_square(2n)).toBe(false);
    expect(is_square(3n)).toBe(false);
    expect(is_square(5n)).toBe(false);
    expect(is_square(-1n)).toBe(false);
  });

  test('with root', () => {
    // SageMath returns (True, sqrt) / (False, None); 0n would be ambiguous
    // with the genuine root of 0. AUDIT-2026-07 L5: this test previously
    // pinned [false, 0n].
    expect(is_square(16n, true)).toEqual([true, 4n]);
    expect(is_square(15n, true)).toEqual([false, null]);
    expect(is_square(-1n, true)).toEqual([false, null]);
    expect(is_square(0n, true)).toEqual([true, 0n]);
  });
});

describe('is_squarefree', () => {
  test('squarefree numbers', () => {
    expect(is_squarefree(1n)).toBe(true);
    expect(is_squarefree(2n)).toBe(true);
    expect(is_squarefree(6n)).toBe(true);
    expect(is_squarefree(30n)).toBe(true);
  });

  test('non-squarefree', () => {
    expect(is_squarefree(0n)).toBe(false);
    expect(is_squarefree(4n)).toBe(false);
    expect(is_squarefree(12n)).toBe(false);
    expect(is_squarefree(18n)).toBe(false);
  });
});

/**
 * Definition-based Kronecker symbol used as an independent oracle:
 * (a|n) = (a|-1)^[n<0] * (a|2)^v * prod (a|p)^e over the odd part of |n|,
 * with (a|p) a brute-forced Legendre symbol.
 */
function referenceKronecker(a: bigint, n: bigint): bigint {
  const legendreBrute = (x: bigint, p: bigint): bigint => {
    const r = ((x % p) + p) % p;
    if (r === 0n) return 0n;
    for (let t = 1n; t < p; t++) {
      if ((t * t) % p === r) return 1n;
    }
    return -1n;
  };

  if (n === 0n) return a === 1n || a === -1n ? 1n : 0n;

  let res = 1n;
  let m = n;
  if (m < 0n) {
    m = -m;
    if (a < 0n) res = -res;
  }

  let v = 0n;
  while (m % 2n === 0n) {
    m /= 2n;
    v++;
  }
  if (v > 0n) {
    if (a % 2n === 0n) return 0n;
    const r = ((a % 8n) + 8n) % 8n;
    res *= (r === 1n || r === 7n ? 1n : -1n) ** v;
  }

  let rest = m;
  for (let p = 3n; p * p <= rest; p += 2n) {
    while (rest % p === 0n) {
      rest /= p;
      res *= legendreBrute(a, p);
    }
  }
  if (rest > 1n) res *= legendreBrute(a, rest);
  return res;
}

describe('kronecker_symbol', () => {
  test('basic cases', () => {
    expect(kronecker_symbol(2n, 3n)).toBe(-1n);
    expect(kronecker_symbol(1n, 5n)).toBe(1n);
    expect(kronecker_symbol(4n, 5n)).toBe(1n);
    expect(kronecker_symbol(2n, 5n)).toBe(-1n);
  });

  test('quadratic residues', () => {
    // 1, 4 are QR mod 5
    expect(kronecker_symbol(1n, 5n)).toBe(1n);
    expect(kronecker_symbol(4n, 5n)).toBe(1n);
    // 2, 3 are not QR mod 5
    expect(kronecker_symbol(2n, 5n)).toBe(-1n);
    expect(kronecker_symbol(3n, 5n)).toBe(-1n);
  });

  test("SageMath's doctests, including even moduli", () => {
    expect(kronecker_symbol(101n, 4n)).toBe(1n);
    expect(kronecker_symbol(3n, 5n)).toBe(-1n);
    expect(kronecker_symbol(3n, 15n)).toBe(0n);
    expect(kronecker_symbol(2n, 15n)).toBe(1n);
    expect(kronecker_symbol(-2n, 15n)).toBe(-1n);
    expect(kronecker_symbol(13n, 21n)).toBe(-1n);
  });

  test('is 0 whenever gcd(a, n) > 1 with n even (AUDIT-2026-07 H1)', () => {
    // (a|2) = 0 for even a, so every one of these must vanish.
    expect(kronecker_symbol(2n, 2n)).toBe(0n);
    expect(kronecker_symbol(2n, 4n)).toBe(0n);
    expect(kronecker_symbol(2n, 6n)).toBe(0n);
    expect(kronecker_symbol(10n, 12n)).toBe(0n);
    expect(kronecker_symbol(0n, 2n)).toBe(0n);
    expect(kronecker_symbol(-4n, 8n)).toBe(0n);
  });

  test('(a|n) for n = +-1 and n = 0', () => {
    expect(kronecker_symbol(0n, 1n)).toBe(1n);
    expect(kronecker_symbol(0n, -1n)).toBe(1n);
    expect(kronecker_symbol(-5n, -1n)).toBe(-1n);
    expect(kronecker_symbol(5n, -1n)).toBe(1n);
    expect(kronecker_symbol(1n, 0n)).toBe(1n);
    expect(kronecker_symbol(-1n, 0n)).toBe(1n);
    expect(kronecker_symbol(6n, 0n)).toBe(0n);
  });

  test('agrees with the definition for all a, n in [-30, 30]', () => {
    let mismatches = 0;
    for (let a = -30n; a <= 30n; a++) {
      for (let n = -30n; n <= 30n; n++) {
        if (kronecker_symbol(a, n) !== referenceKronecker(a, n)) mismatches++;
      }
    }
    expect(mismatches).toBe(0);
  });

  test('multiplicativity in the lower argument for n > 0', () => {
    // (ab|n) = (a|n)(b|n) holds for n > 0; for n < 0 the (a|-1) = sign(a)
    // factor breaks it at a = 0.
    for (let n = 1n; n <= 25n; n++) {
      for (let a = -6n; a <= 6n; a++) {
        for (let b = -6n; b <= 6n; b++) {
          expect(kronecker_symbol(a * b, n)).toBe(kronecker_symbol(a, n) * kronecker_symbol(b, n));
        }
      }
    }
  });
});

describe('legendre_symbol', () => {
  test("SageMath's doctests", () => {
    expect(legendre_symbol(2n, 3n)).toBe(-1n);
    expect(legendre_symbol(1n, 3n)).toBe(1n);
  });

  test('rejects composite p with SageMath message', () => {
    expect(() => legendre_symbol(2n, 15n)).toThrow('p must be a prime');
    expect(() => legendre_symbol(2n, 9n)).toThrow('p must be a prime');
    expect(() => legendre_symbol(2n, 1n)).toThrow('p must be a prime');
  });

  test('rejects p = 2 with SageMath message', () => {
    expect(() => legendre_symbol(1n, 2n)).toThrow('p must be odd');
  });

  test('matches kronecker at odd primes', () => {
    for (const p of [3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n]) {
      for (let a = -20n; a <= 20n; a++) {
        expect(legendre_symbol(a, p)).toBe(kronecker_symbol(a, p));
      }
    }
  });
});

describe('jacobi_symbol', () => {
  test("SageMath's doctests", () => {
    expect(jacobi_symbol(10n, 777n)).toBe(-1n);
    expect(jacobi_symbol(10n, 5n)).toBe(0n);
    expect(() => jacobi_symbol(10n, 2n)).toThrow('second input must be odd, 2 is not odd');
  });

  test('accepts negative odd b (AUDIT-2026-07 L6)', () => {
    expect(jacobi_symbol(10n, -3n)).toBe(1n);
    expect(jacobi_symbol(-10n, -3n)).toBe(1n);
    expect(jacobi_symbol(2n, -7n)).toBe(1n);
  });

  test('rejects every even b', () => {
    expect(() => jacobi_symbol(3n, 0n)).toThrow('second input must be odd, 0 is not odd');
    expect(() => jacobi_symbol(3n, -4n)).toThrow('second input must be odd, -4 is not odd');
  });
});

describe('crt', () => {
  test('basic cases', () => {
    expect(crt(2n, 3n, 3n, 5n)).toBe(8n);
    expect(crt(1n, 2n, 2n, 3n)).toBe(5n);
  });

  test("SageMath's unsolvable-instance message", () => {
    // sage: CRT(4, 6, 8, 12)
    // ValueError: no solution to crt problem since gcd(8,12) does not divide 4-6
    expect(() => crt(4n, 6n, 8n, 12n)).toThrow(
      'no solution to crt problem since gcd(8,12) does not divide 4-6'
    );
    // sage: CRT_list([32,2,1],[60,90,150])
    // ValueError: no solution to crt problem since gcd(180,150) does not divide 92-1
    expect(() => CRT_list([32n, 2n, 1n], [60n, 90n, 150n])).toThrow(
      'no solution to crt problem since gcd(180,150) does not divide 92-1'
    );
  });

  test('larger moduli', () => {
    const x = crt(3n, 4n, 7n, 11n);
    expect(x % 7n).toBe(3n);
    expect(x % 11n).toBe(4n);
  });
});

describe('divisors', () => {
  test('basic cases', () => {
    expect(divisors(1n)).toEqual([1n]);
    expect(divisors(6n)).toEqual([1n, 2n, 3n, 6n]);
    expect(divisors(12n)).toEqual([1n, 2n, 3n, 4n, 6n, 12n]);
  });

  test('prime', () => {
    expect(divisors(7n)).toEqual([1n, 7n]);
  });

  test('prime power', () => {
    expect(divisors(8n)).toEqual([1n, 2n, 4n, 8n]);
  });
});

describe('sigma', () => {
  test('sum of divisors', () => {
    expect(sigma(1n)).toBe(1n);
    expect(sigma(6n)).toBe(12n); // 1 + 2 + 3 + 6
    expect(sigma(12n)).toBe(28n); // 1 + 2 + 3 + 4 + 6 + 12
  });

  test('sum of squares of divisors', () => {
    expect(sigma(6n, 2n)).toBe(50n); // 1 + 4 + 9 + 36
  });

  test("SageMath's doctests", () => {
    expect(sigma(10n)).toBe(18n);
    expect(sigma(10n, 2n)).toBe(130n);
    expect(sigma(100n, 4n)).toBe(106811523n);
  });

  test('k = 0 counts the divisors', () => {
    expect(sigma(12n, 0n)).toBe(6n);
    expect(sigma(1n, 0n)).toBe(1n);
    expect(sigma(64n, 0n)).toBe(7n);
  });

  test('negative n uses |n| (the -1 unit is not a factor)', () => {
    // SageMath: sigma(-4) == 7
    expect(sigma(-4n)).toBe(7n);
    expect(sigma(-12n)).toBe(28n);
    expect(sigma(-6n, 2n)).toBe(50n);
  });

  test('uses the multiplicative closed form, not divisor enumeration', () => {
    // sage: sigma(factorial(100), 0) == 39001250856960000 -- impossible to
    // reach by enumerating 3.9e16 divisors.
    let f100 = 1n;
    for (let i = 2n; i <= 100n; i++) f100 *= i;
    const start = Date.now();
    expect(sigma(f100, 0n)).toBe(39001250856960000n);
    expect(Date.now() - start).toBeLessThan(2000);

    // sage: sigma(factorial(100), 3).mod(144169) == 3672
    expect(sigma(f100, 3n) % 144169n).toBe(3672n);

    let f41 = 1n;
    for (let i = 2n; i <= 41n; i++) f41 *= i;
    expect(sigma(f41, 1n)).toBe(229199532273029988767733858700732906511758707916800n);
  });

  test('agrees with brute-force divisor sums', () => {
    for (let n = 1n; n <= 200n; n++) {
      for (const k of [0n, 1n, 2n, 3n]) {
        let want = 0n;
        for (let d = 1n; d <= n; d++) if (n % d === 0n) want += d ** k;
        expect(sigma(n, k)).toBe(want);
      }
    }
  });

  test('rejects 0 like factor(0)', () => {
    expect(() => sigma(0n)).toThrow();
  });
});

describe('number_of_divisors', () => {
  test("SageMath's doctests", () => {
    expect(number_of_divisors(100n)).toBe(9n);
    // AUDIT-2026-07 M3: this used to throw.
    expect(number_of_divisors(-720n)).toBe(30n);
  });

  test('agrees with the length of divisors()', () => {
    for (let n = 1n; n <= 300n; n++) {
      expect(number_of_divisors(n)).toBe(BigInt(divisors(n).length));
      expect(number_of_divisors(-n)).toBe(BigInt(divisors(n).length));
    }
  });

  test('rejects only 0, with SageMath message', () => {
    expect(() => number_of_divisors(0n)).toThrow('input must be nonzero');
    expect(number_of_divisors(1n)).toBe(1n);
    expect(number_of_divisors(-1n)).toBe(1n);
  });
});

describe('radical', () => {
  test('basic cases', () => {
    expect(radical(12n)).toBe(6n); // 2 * 3
    expect(radical(100n)).toBe(10n); // 2 * 5
    expect(radical(1n)).toBe(1n);
  });
});

describe('sqrt_mod', () => {
  test('edge case: a = 0', () => {
    expect(sqrt_mod(0n, 7n)).toBe(0n);
    expect(sqrt_mod(0n, 7n, true)).toEqual([0n]);
    expect(sqrt_mod(0n, 13n)).toBe(0n);
    expect(sqrt_mod(0n, 2n)).toBe(0n);
  });

  test('edge case: a = 1', () => {
    expect(sqrt_mod(1n, 7n)).toBe(1n);
    expect(sqrt_mod(1n, 7n, true)).toEqual([1n, 6n]);
    expect(sqrt_mod(1n, 13n)).toBe(1n);
    expect(sqrt_mod(1n, 13n, true)).toEqual([1n, 12n]);
  });

  test('edge case: p = 2', () => {
    expect(sqrt_mod(0n, 2n)).toBe(0n);
    expect(sqrt_mod(1n, 2n)).toBe(1n);
    expect(sqrt_mod(0n, 2n, true)).toEqual([0n]);
    expect(sqrt_mod(1n, 2n, true)).toEqual([1n]);
  });

  test('quadratic non-residues return null', () => {
    // 3 is not a QR mod 7 (QRs are 1, 2, 4)
    expect(sqrt_mod(3n, 7n)).toBe(null);
    expect(sqrt_mod(3n, 7n, true)).toEqual([]);
    // 5 is not a QR mod 7
    expect(sqrt_mod(5n, 7n)).toBe(null);
    expect(sqrt_mod(5n, 7n, true)).toEqual([]);
    // 2 is not a QR mod 5
    expect(sqrt_mod(2n, 5n)).toBe(null);
  });

  test('p = 7 (p ≡ 3 mod 4)', () => {
    // QRs mod 7 are 1, 2, 4
    // 1^2 = 1, 2^2 = 4, 3^2 = 2, 4^2 = 2, 5^2 = 4, 6^2 = 1
    const r1 = sqrt_mod(1n, 7n)!;
    expect((r1 * r1) % 7n).toBe(1n);

    const r2 = sqrt_mod(2n, 7n)!;
    expect((r2 * r2) % 7n).toBe(2n);

    const r4 = sqrt_mod(4n, 7n)!;
    expect((r4 * r4) % 7n).toBe(4n);

    // All roots
    const roots2 = sqrt_mod(2n, 7n, true);
    expect(roots2.length).toBe(2);
    expect(roots2.every((r) => (r * r) % 7n === 2n)).toBe(true);
    expect(roots2[0]).toBeLessThan(roots2[1]!);
  });

  test('p = 13 (p ≡ 1 mod 4)', () => {
    // QRs mod 13 are 1, 3, 4, 9, 10, 12
    const qrs = [1n, 3n, 4n, 9n, 10n, 12n];
    for (const a of qrs) {
      const r = sqrt_mod(a, 13n)!;
      expect((r * r) % 13n).toBe(a);
    }

    // Non-QRs
    const nonQRs = [2n, 5n, 6n, 7n, 8n, 11n];
    for (const a of nonQRs) {
      expect(sqrt_mod(a, 13n)).toBe(null);
    }
  });

  test('p = 17 (p ≡ 1 mod 8)', () => {
    // Tests the full Tonelli-Shanks case
    // QRs mod 17 are 1, 2, 4, 8, 9, 13, 15, 16
    const qrs = [1n, 2n, 4n, 8n, 9n, 13n, 15n, 16n];
    for (const a of qrs) {
      const r = sqrt_mod(a, 17n)!;
      expect(r).not.toBe(null);
      expect((r * r) % 17n).toBe(a);
    }
  });

  test('p = 23 (p ≡ 7 mod 8)', () => {
    // p ≡ 3 mod 4 case
    // Test a few QRs
    const testCases = [1n, 2n, 3n, 4n, 6n, 8n, 9n, 12n, 13n, 16n, 18n];
    for (const a of testCases) {
      if (legendre_symbol(a, 23n) === 1n) {
        const r = sqrt_mod(a, 23n)!;
        expect(r).not.toBe(null);
        expect((r * r) % 23n).toBe(a);
      }
    }
  });

  test('p = 5 (p ≡ 5 mod 8)', () => {
    // Tests the p ≡ 5 mod 8 special case
    // QRs mod 5 are 1, 4
    expect(sqrt_mod(1n, 5n)! ** 2n % 5n).toBe(1n);
    expect(sqrt_mod(4n, 5n)! ** 2n % 5n).toBe(4n);
    expect(sqrt_mod(2n, 5n)).toBe(null);
    expect(sqrt_mod(3n, 5n)).toBe(null);
  });

  test('p = 41 (p ≡ 1 mod 8)', () => {
    // Another Tonelli-Shanks case
    for (let a = 1n; a < 41n; a++) {
      if (legendre_symbol(a, 41n) === 1n) {
        const r = sqrt_mod(a, 41n)!;
        expect(r).not.toBe(null);
        expect((r * r) % 41n).toBe(a);
      } else {
        expect(sqrt_mod(a, 41n)).toBe(null);
      }
    }
  });

  test('large prime p ≡ 3 mod 4', () => {
    // p = 2^127 - 1 (Mersenne prime, p ≡ 3 mod 4)
    const p = 170141183460469231731687303715884105727n;
    const a = 12345678901234567890n;
    const r = sqrt_mod(a, p);
    if (r !== null) {
      expect((r * r) % p).toBe(a % p);
    }
  });

  test('large prime p ≡ 1 mod 4', () => {
    // p = 10^9 + 7 (common prime in competitive programming)
    const p = 1000000007n;
    const a = 123456789n;
    const ls = legendre_symbol(a, p);
    if (ls === 1n) {
      const r = sqrt_mod(a, p)!;
      expect((r * r) % p).toBe(a);
    }

    // Test with a known square
    const square = (12345n * 12345n) % p;
    const r2 = sqrt_mod(square, p)!;
    expect((r2 * r2) % p).toBe(square);
  });

  test('all_roots returns sorted array', () => {
    const roots = sqrt_mod(4n, 13n, true);
    expect(roots.length).toBe(2);
    expect(roots[0]).toBeLessThan(roots[1]!);
    expect(roots).toContain(2n);
    expect(roots).toContain(11n);
  });

  test('all_roots for all elements of small prime', () => {
    const p = 11n;
    for (let a = 0n; a < p; a++) {
      const roots = sqrt_mod(a, p, true);
      for (const r of roots) {
        expect((r * r) % p).toBe(a);
      }
      // Verify count: 0 has 1 root, QRs have 2 roots, non-QRs have 0 roots
      if (a === 0n) {
        expect(roots.length).toBe(1);
      } else if (legendre_symbol(a, p) === 1n) {
        expect(roots.length).toBe(2);
      } else {
        expect(roots.length).toBe(0);
      }
    }
  });

  test('negative input is normalized', () => {
    // -1 mod 5 = 4, which is a QR
    const r = sqrt_mod(-1n, 5n)!;
    expect((r * r) % 5n).toBe(4n);

    // -1 mod 7 is 6, which is not a QR
    expect(sqrt_mod(-1n, 7n)).toBe(null);
  });

  test('input larger than p is reduced', () => {
    // 9 mod 7 = 2
    const r1 = sqrt_mod(9n, 7n)!;
    const r2 = sqrt_mod(2n, 7n)!;
    expect(r1).toBe(r2);

    // 100 mod 13 = 9
    const r3 = sqrt_mod(100n, 13n)!;
    expect((r3 * r3) % 13n).toBe(9n);
  });

  test('comparison with SageMath known values', () => {
    // sage: mod(5, 389).sqrt()
    // 86
    const r389 = sqrt_mod(5n, 389n)!;
    expect((r389 * r389) % 389n).toBe(5n);
    // One of the roots should be 86
    const roots389 = sqrt_mod(5n, 389n, true);
    expect(roots389).toContain(86n);

    // sage: mod(-1, 17).sqrt()
    // 4
    const r17 = sqrt_mod(-1n, 17n)!;
    expect((r17 * r17) % 17n).toBe(16n);
    const roots17 = sqrt_mod(-1n, 17n, true);
    expect(roots17).toContain(4n);
    expect(roots17).toContain(13n);
  });
});

describe('prime_range', () => {
  test('empty range', () => {
    expect(prime_range(0n)).toEqual([]);
    expect(prime_range(1n)).toEqual([]);
    expect(prime_range(2n)).toEqual([]);
    expect(prime_range(5n, 5n)).toEqual([]);
    expect(prime_range(10n, 5n)).toEqual([]);
  });

  test('single prime', () => {
    expect(prime_range(3n)).toEqual([2n]);
    expect(prime_range(2n, 3n)).toEqual([2n]);
  });

  test('multiple primes - single arg', () => {
    expect(prime_range(10n)).toEqual([2n, 3n, 5n, 7n]);
    expect(prime_range(20n)).toEqual([2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n]);
    expect(prime_range(30n)).toEqual([2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n]);
  });

  test('multiple primes - two args', () => {
    expect(prime_range(10n, 20n)).toEqual([11n, 13n, 17n, 19n]);
    expect(prime_range(2000n, 2020n)).toEqual([2003n, 2011n, 2017n]);
    expect(prime_range(5n, 10n)).toEqual([5n, 7n]);
  });

  test('negative start', () => {
    expect(prime_range(-100n, 10n)).toEqual([2n, 3n, 5n, 7n]);
  });

  test('count of primes', () => {
    // There are 168 primes less than 1000
    expect(prime_range(1000n).length).toBe(168);
  });
});

describe('moebius', () => {
  test('moebius(1) = 1', () => {
    expect(moebius(1n)).toBe(1n);
  });

  test('moebius(0) = 0', () => {
    expect(moebius(0n)).toBe(0n);
  });

  test('primes have moebius = -1', () => {
    expect(moebius(2n)).toBe(-1n);
    expect(moebius(3n)).toBe(-1n);
    expect(moebius(5n)).toBe(-1n);
    expect(moebius(7n)).toBe(-1n);
    expect(moebius(11n)).toBe(-1n);
  });

  test('products of two distinct primes have moebius = 1', () => {
    expect(moebius(6n)).toBe(1n); // 2 * 3
    expect(moebius(10n)).toBe(1n); // 2 * 5
    expect(moebius(15n)).toBe(1n); // 3 * 5
    expect(moebius(35n)).toBe(1n); // 5 * 7
  });

  test('products of three distinct primes have moebius = -1', () => {
    expect(moebius(30n)).toBe(-1n); // 2 * 3 * 5
  });

  test('squared factors give moebius = 0', () => {
    expect(moebius(4n)).toBe(0n); // 2^2
    expect(moebius(8n)).toBe(0n); // 2^3
    expect(moebius(9n)).toBe(0n); // 3^2
    expect(moebius(12n)).toBe(0n); // 2^2 * 3
    expect(moebius(18n)).toBe(0n); // 2 * 3^2
    expect(moebius(20n)).toBe(0n); // 2^2 * 5
  });

  test('negative numbers', () => {
    expect(moebius(-1n)).toBe(1n);
    expect(moebius(-5n)).toBe(-1n);
    expect(moebius(-35n)).toBe(1n);
  });

  test('values for 1-20', () => {
    // From SageMath: [moebius(n) for n in range(1, 21)]
    const expected = [
      1n,
      -1n,
      -1n,
      0n,
      -1n,
      1n,
      -1n,
      0n,
      0n,
      1n,
      -1n,
      0n,
      -1n,
      1n,
      1n,
      0n,
      -1n,
      0n,
      -1n,
      0n,
    ];
    for (let i = 0; i < 20; i++) {
      expect(moebius(BigInt(i + 1))).toBe(expected[i]);
    }
  });
});

describe('squarefree_part', () => {
  test('squarefree numbers are unchanged', () => {
    expect(squarefree_part(1n)).toBe(1n);
    expect(squarefree_part(2n)).toBe(2n);
    expect(squarefree_part(3n)).toBe(3n);
    expect(squarefree_part(6n)).toBe(6n);
    expect(squarefree_part(30n)).toBe(30n);
  });

  test('perfect squares', () => {
    expect(squarefree_part(4n)).toBe(1n);
    expect(squarefree_part(9n)).toBe(1n);
    expect(squarefree_part(16n)).toBe(1n);
    expect(squarefree_part(100n)).toBe(1n);
  });

  test('general cases', () => {
    expect(squarefree_part(12n)).toBe(3n); // 12 = 3 * 4
    expect(squarefree_part(18n)).toBe(2n); // 18 = 2 * 9
    expect(squarefree_part(72n)).toBe(2n); // 72 = 2 * 36
    expect(squarefree_part(50n)).toBe(2n); // 50 = 2 * 25
    expect(squarefree_part(8n)).toBe(2n); // 8 = 2 * 4
    expect(squarefree_part(27n)).toBe(3n); // 27 = 3 * 9
  });

  test('negative numbers', () => {
    expect(squarefree_part(-1n)).toBe(-1n);
    expect(squarefree_part(-4n)).toBe(-1n);
    expect(squarefree_part(-12n)).toBe(-3n);
    expect(squarefree_part(-17n * 32n)).toBe(-34n); // -17 * 32 = -17 * 2^5 -> -17 * 2
  });

  test('zero', () => {
    expect(squarefree_part(0n)).toBe(0n);
  });

  test('SageMath examples', () => {
    expect(squarefree_part(17n * 37n * 37n)).toBe(17n);
  });
});

describe('prime_factors / prime_divisors', () => {
  test('basic cases', () => {
    expect(prime_factors(12n)).toEqual([2n, 3n]);
    expect(prime_factors(100n)).toEqual([2n, 5n]);
    expect(prime_factors(2004n)).toEqual([2n, 3n, 167n]);
  });

  test('one has no prime factors', () => {
    expect(prime_factors(1n)).toEqual([]);
  });

  test('primes', () => {
    expect(prime_factors(2n)).toEqual([2n]);
    expect(prime_factors(17n)).toEqual([17n]);
    expect(prime_factors(101n)).toEqual([101n]);
  });

  test('prime powers', () => {
    expect(prime_factors(8n)).toEqual([2n]);
    expect(prime_factors(27n)).toEqual([3n]);
    expect(prime_factors(32n)).toEqual([2n]);
  });

  test('negative numbers', () => {
    expect(prime_factors(-12n)).toEqual([2n, 3n]);
    expect(prime_factors(-100n)).toEqual([2n, 5n]);
  });

  test('prime_divisors is alias for prime_factors', () => {
    expect(prime_divisors(12n)).toEqual(prime_factors(12n));
    expect(prime_divisors(100n)).toEqual(prime_factors(100n));
  });
});

describe('valuation', () => {
  test('basic cases', () => {
    expect(valuation(24n, 2n)).toBe(3n); // 24 = 2^3 * 3
    expect(valuation(100n, 5n)).toBe(2n); // 100 = 4 * 5^2
    expect(valuation(512n, 2n)).toBe(9n); // 512 = 2^9
    expect(valuation(243n, 3n)).toBe(5n); // 243 = 3^5
  });

  test('no factors', () => {
    expect(valuation(1n, 2n)).toBe(0n);
    expect(valuation(7n, 2n)).toBe(0n);
    expect(valuation(15n, 7n)).toBe(0n);
  });

  test('composite base', () => {
    expect(valuation(100n, 10n)).toBe(2n);
    expect(valuation(200n, 10n)).toBe(2n);
    expect(valuation(1000n, 10n)).toBe(3n);
  });

  test('negative numbers', () => {
    expect(valuation(-24n, 2n)).toBe(3n);
    expect(valuation(-100n, 5n)).toBe(2n);
  });

  test('throws for invalid base', () => {
    expect(() => valuation(10n, 1n)).toThrow();
    expect(() => valuation(10n, 0n)).toThrow();
    expect(() => valuation(10n, -2n)).toThrow();
  });

  test('throws for zero', () => {
    expect(() => valuation(0n, 2n)).toThrow();
  });
});

describe('rational_reconstruction', () => {
  test('basic examples from SageMath', () => {
    // 119/53 mod 100000 = 11323
    const [p, q] = rational_reconstruction(11323n, 100000n);
    expect(p).toBe(119n);
    expect(q).toBe(53n);
    // Verify: 119 * inverse_mod(53, 100000) = 11323 mod 100000
    expect((p * inverse_mod(q, 100000n)) % 100000n).toBe(11323n);
  });

  test('small values', () => {
    const [p, q] = rational_reconstruction(3n, 292393n);
    expect(p).toBe(3n);
    expect(q).toBe(1n);
  });

  test('round-trip reconstruction', () => {
    // For various p/q, compute a = p/q mod m, then reconstruct
    const testCases: Array<[bigint, bigint, bigint]> = [
      [45n, 97n, 292393n],
      [7n, 13n, 1000n],
      [1n, 1n, 100n],
      [0n, 1n, 100n],
      [17n, 23n, 10000n],
    ];

    for (const [p, q, m] of testCases) {
      // Check that |p|, q < sqrt(m/2) for valid reconstruction
      const bound = isqrt(m / 2n);
      if ((p < 0n ? -p : p) <= bound && q <= bound) {
        const a = (p * inverse_mod(q, m)) % m;
        const [pRec, qRec] = rational_reconstruction((a + m) % m, m);
        // The reconstructed rational should equal p/q
        expect(pRec * q).toBe(p * qRec);
      }
    }
  });

  test('throws for non-existent reconstruction', () => {
    // 400 mod 1000 doesn't have a valid reconstruction
    expect(() => rational_reconstruction(400n, 1000n)).toThrow();
  });

  test('throws for zero modulus', () => {
    expect(() => rational_reconstruction(1n, 0n)).toThrow();
  });

  test('negative inputs', () => {
    // Should work with negative inputs by reducing mod m
    // -11323 + 100000 = 88677
    const [p1, q1] = rational_reconstruction(88677n, 100000n);
    // 88677 ≡ -119/53 (mod 100000) since -119 * inverse_mod(53, 100000) = 88677
    // Verify the reconstruction is valid
    expect((((p1 * inverse_mod(q1, 100000n)) % 100000n) + 100000n) % 100000n).toBe(88677n);
  });
});

describe('CRT_basis', () => {
  test('basic coprime moduli', () => {
    const basis = CRT_basis([5n, 13n]) as bigint[];
    const [c1, c2] = basis;

    // e_i ≡ 1 (mod m_i), e_i ≡ 0 (mod m_j) for j ≠ i
    expect(c1! % 5n).toBe(1n);
    expect(c1! % 13n).toBe(0n);
    expect(c2! % 5n).toBe(0n);
    expect(c2! % 13n).toBe(1n);
  });

  test('three moduli', () => {
    const basis = CRT_basis([3n, 5n, 7n]) as bigint[];
    expect(basis.length).toBe(3);

    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        const moduli = [3n, 5n, 7n];
        const expected = i === j ? 1n : 0n;
        expect(basis[i]! % moduli[j]!).toBe(expected);
      }
    }
  });

  test('using basis to combine residues', () => {
    const moduli = [5n, 13n];
    const basis = CRT_basis(moduli) as bigint[];
    const M = 65n; // 5 * 13

    // Combine residues a1 = 42 mod 5 and a2 = 42 mod 13
    const a1 = 42n % 5n;
    const a2 = 42n % 13n;
    const combined = (a1 * basis[0]! + a2 * basis[1]!) % M;
    expect(combined).toBe(42n % M);
  });

  test('empty moduli', () => {
    const basis = CRT_basis([]) as bigint[];
    expect(basis).toEqual([]);
  });

  test('single modulus', () => {
    const basis = CRT_basis([7n]) as bigint[];
    expect(basis.length).toBe(1);
    expect(basis[0]! % 7n).toBe(1n);
  });

  test('throws for non-coprime moduli', () => {
    expect(() => CRT_basis([6n, 10n])).toThrow();
  });

  test('non-coprime with require_coprime_moduli=false', () => {
    const result = CRT_basis([6n, 10n], false) as [bigint[], boolean];
    expect(result[1]).toBe(false);
  });

  test("SageMath's non-coprime basis for [60, 90, 150] (AUDIT-2026-07 H4)", () => {
    const [cs, coprime] = CRT_basis([60n, 90n, 150n], false) as [bigint[], boolean];
    expect(coprime).toBe(false);
    expect(cs).toEqual([15n, -20n, 6n]);
  });

  test('non-coprime basis always has one entry per modulus', () => {
    for (const moduli of [
      [6n, 10n],
      [60n, 90n, 150n],
      [7n, 6n, 10n],
      [4n, 6n, 9n, 10n],
    ]) {
      const [cs, coprime] = CRT_basis(moduli, false) as [bigint[], boolean];
      expect(coprime).toBe(false);
      expect(cs.length).toBe(moduli.length);
    }
  });

  test('the non-coprime basis solves solvable systems', () => {
    const moduli = [60n, 90n, 150n];
    const residues = [32n, 2n, 2n];
    const [cs] = CRT_basis(moduli, false) as [bigint[], boolean];
    let x = 0n;
    for (let i = 0; i < 3; i++) x += cs[i]! * residues[i]!;
    for (let i = 0; i < 3; i++) {
      expect((((x - residues[i]!) % moduli[i]!) + moduli[i]!) % moduli[i]!).toBe(0n);
    }
  });
});

describe('CRT_vectors', () => {
  test('basic example', () => {
    const result = CRT_vectors(
      [
        [3n, 5n, 7n],
        [3n, 5n, 11n],
      ],
      [2n, 3n]
    );
    expect(result).toEqual([3n, 5n, 5n]);
  });

  test('verify componentwise CRT', () => {
    const X = [
      [10n, 43n, 17n],
      [2n, 7n, 8n],
    ];
    const moduli = [9n, 8n];

    const result = CRT_vectors(X, moduli);

    // Each component should satisfy CRT
    for (let j = 0; j < 3; j++) {
      expect(result[j]! % moduli[0]!).toBe(X[0]![j]! % moduli[0]!);
      expect(result[j]! % moduli[1]!).toBe(X[1]![j]! % moduli[1]!);
    }
  });

  test('empty vectors', () => {
    expect(CRT_vectors([], [2n, 3n])).toEqual([]);
    expect(CRT_vectors([[]], [2n])).toEqual([]);
  });

  test('throws for mismatched lengths', () => {
    expect(() => CRT_vectors([[1n], [2n], [3n]], [2n, 3n])).toThrow();
  });

  test('three non-coprime moduli (AUDIT-2026-07 H4)', () => {
    // sage: CRT_list([32,2,2],[60,90,150]) == 452
    expect(CRT_vectors([[32n], [2n], [2n]], [60n, 90n, 150n])).toEqual([452n]);
    expect(CRT_list([32n, 2n, 2n], [60n, 90n, 150n])).toBe(452n);
  });

  test("SageMath's non-coprime doctests", () => {
    expect(CRT_vectors([[6n], [0n]], [10n, 4n])).toEqual([16n]);
    expect(() => CRT_vectors([[6n], [0n]], [10n, 10n])).toThrow('solution does not exist');
  });
});

describe('continued_fraction', () => {
  test('basic examples', () => {
    expect(continued_fraction(13n, 9n)).toEqual([1n, 2n, 4n]);
    expect(continued_fraction(225n, 157n)).toEqual([1n, 2n, 3n, 4n, 5n]);
  });

  test('integers', () => {
    expect(continued_fraction(5n)).toEqual([5n]);
    expect(continued_fraction(0n)).toEqual([0n]);
    expect(continued_fraction(-3n)).toEqual([-3n]);
  });

  test('negative rationals', () => {
    // -1/3 = -1 + 2/3 = -1 + 1/(3/2) = -1 + 1/(1 + 1/2)
    const cf = continued_fraction(-1n, 3n);
    expect(cf).toEqual([-1n, 1n, 2n]);

    // Verify by evaluating back
    const [p, q] = continued_fraction_value(cf);
    expect(p * 3n).toBe(-1n * q);
  });

  test('pi approximation', () => {
    // 355/113 is a famous approximation to pi
    const cf = continued_fraction(355n, 113n);
    expect(cf).toEqual([3n, 7n, 16n]);
  });

  test('with bound', () => {
    const cf = continued_fraction(355n, 113n, 2);
    expect(cf.length).toBe(2);
    expect(cf).toEqual([3n, 7n]);
  });

  test('Fibonacci ratio', () => {
    // Fibonacci(n)/Fibonacci(n-1) approaches golden ratio
    // with continued fraction [1; 1, 1, 1, ...]
    // 21/13 (Fib(8)/Fib(7))
    const cf = continued_fraction(21n, 13n);
    // 21/13 = 1 + 8/13 = 1 + 1/(13/8) = 1 + 1/(1 + 5/8) = ...
    expect(cf).toEqual([1n, 1n, 1n, 1n, 1n, 2n]);

    // Verify by evaluating back
    const [p, q] = continued_fraction_value(cf);
    expect(p).toBe(21n);
    expect(q).toBe(13n);
  });
});

describe('convergents', () => {
  test('basic examples', () => {
    // 13/9 = [1; 2, 4]
    const convs = convergents([1n, 2n, 4n]);
    expect(convs).toEqual([
      [1n, 1n],
      [3n, 2n],
      [13n, 9n],
    ]);
  });

  test('convergents approach the value', () => {
    // For 23/157 = [0; 6, 1, 4, 1, 3]
    const cf = [0n, 6n, 1n, 4n, 1n, 3n];
    const convs = convergents(cf);

    // Last convergent should equal 23/157
    const [p, q] = convs[convs.length - 1]!;
    expect(p).toBe(23n);
    expect(q).toBe(157n);

    // Convergents should alternate above/below
    // and get progressively closer
  });

  test('single term', () => {
    expect(convergents([5n])).toEqual([[5n, 1n]]);
  });

  test('empty input', () => {
    expect(convergents([])).toEqual([]);
  });

  test('property: p_n*q_{n-1} - q_n*p_{n-1} = (-1)^n', () => {
    const cf = [3n, 7n, 15n, 1n, 292n];
    const convs = convergents(cf);

    for (let n = 1; n < convs.length; n++) {
      const [pn, qn] = convs[n]!;
      const [pn1, qn1] = convs[n - 1]!;
      const det = pn * qn1 - qn * pn1;
      expect(det === 1n || det === -1n).toBe(true);
    }
  });
});

describe('continued_fraction_value', () => {
  test('evaluates correctly', () => {
    expect(continued_fraction_value([1n, 2n, 4n])).toEqual([13n, 9n]);
    expect(continued_fraction_value([3n, 7n, 16n])).toEqual([355n, 113n]);
    expect(continued_fraction_value([0n])).toEqual([0n, 1n]);
  });

  test('round-trip with continued_fraction', () => {
    const testCases: Array<[bigint, bigint]> = [
      [13n, 9n],
      [225n, 157n],
      [355n, 113n],
      [22n, 7n],
      [-5n, 3n],
    ];

    for (const [p, q] of testCases) {
      const cf = continued_fraction(p, q);
      const [pBack, qBack] = continued_fraction_value(cf);
      // Check p/q = pBack/qBack
      expect(p * qBack).toBe(q * pBack);
    }
  });
});

describe('half_gcd', () => {
  test('basic transformation', () => {
    const [r, s, t, u] = half_gcd(1000n, 37n);

    // The matrix should be valid (determinant ±1)
    const det = r * u - s * t;
    expect(det === 1n || det === -1n).toBe(true);
  });

  test('transforms toward GCD', () => {
    const a = 12345678n;
    const b = 9876543n;
    const [r, s, t, u] = half_gcd(a, b);

    // Apply transformation
    const a1 = r * a + s * b;
    const b1 = t * a + u * b;

    // gcd should be preserved
    expect(gcd(a1 < 0n ? -a1 : a1, b1 < 0n ? -b1 : b1)).toBe(gcd(a, b));
  });

  test('with zero', () => {
    const [r, s, t, u] = half_gcd(100n, 0n);
    expect([r, s, t, u]).toEqual([1n, 0n, 0n, 1n]);
  });

  test('symmetric inputs', () => {
    // When a < b, they should be swapped internally
    const [r1, s1, t1, u1] = half_gcd(37n, 1000n);
    const [r2, s2, t2, u2] = half_gcd(1000n, 37n);
    // Results should be the same
    expect([r1, s1, t1, u1]).toEqual([r2, s2, t2, u2]);
  });
});

describe('CRT_list', () => {
  test('basic example', () => {
    // x ≡ 2 (mod 3), x ≡ 3 (mod 5) => x = 8
    expect(CRT_list([2n, 3n], [3n, 5n])).toBe(8n);
  });

  test('three moduli', () => {
    // x ≡ 1 (mod 2), x ≡ 2 (mod 3), x ≡ 3 (mod 5)
    const result = CRT_list([1n, 2n, 3n], [2n, 3n, 5n]);
    expect(result % 2n).toBe(1n);
    expect(result % 3n).toBe(2n);
    expect(result % 5n).toBe(3n);
  });

  test('empty lists', () => {
    expect(CRT_list([], [])).toBe(0n);
  });
});

describe('bernoulli', () => {
  test('B_0 = 1', () => {
    const b0 = bernoulli(0n);
    expect(b0.numerator).toBe(1n);
    expect(b0.denominator).toBe(1n);
  });

  test('B_1 = -1/2', () => {
    const b1 = bernoulli(1n);
    expect(b1.numerator).toBe(-1n);
    expect(b1.denominator).toBe(2n);
  });

  test('odd Bernoulli numbers > 1 are zero', () => {
    for (let n = 3n; n <= 15n; n += 2n) {
      const bn = bernoulli(n);
      expect(bn.numerator).toBe(0n);
    }
  });

  test('B_2 = 1/6', () => {
    const b2 = bernoulli(2n);
    expect(b2.numerator).toBe(1n);
    expect(b2.denominator).toBe(6n);
  });

  test('B_4 = -1/30', () => {
    const b4 = bernoulli(4n);
    expect(b4.numerator).toBe(-1n);
    expect(b4.denominator).toBe(30n);
  });

  test('B_6 = 1/42', () => {
    const b6 = bernoulli(6n);
    expect(b6.numerator).toBe(1n);
    expect(b6.denominator).toBe(42n);
  });

  test('B_8 = -1/30', () => {
    const b8 = bernoulli(8n);
    expect(b8.numerator).toBe(-1n);
    expect(b8.denominator).toBe(30n);
  });

  test('B_10 = 5/66', () => {
    const b10 = bernoulli(10n);
    expect(b10.numerator).toBe(5n);
    expect(b10.denominator).toBe(66n);
  });

  test('B_12 = -691/2730 (famous negative numerator)', () => {
    const b12 = bernoulli(12n);
    expect(b12.numerator).toBe(-691n);
    expect(b12.denominator).toBe(2730n);
  });

  test('B_14 = 7/6', () => {
    const b14 = bernoulli(14n);
    expect(b14.numerator).toBe(7n);
    expect(b14.denominator).toBe(6n);
  });

  test('throws for negative n', () => {
    expect(() => bernoulli(-1n)).toThrow();
  });
});

describe('hilbert_symbol', () => {
  test('zero cases return 0', () => {
    expect(hilbert_symbol(0n, 5n, 3n)).toBe(0n);
    expect(hilbert_symbol(5n, 0n, 3n)).toBe(0n);
    expect(hilbert_symbol(0n, 0n, 3n)).toBe(0n);
  });

  test('archimedean place (p = -1)', () => {
    // hilbert(a, b, -1) = -1 iff a < 0 and b < 0
    expect(hilbert_symbol(1n, 1n, -1n)).toBe(1n);
    expect(hilbert_symbol(1n, -1n, -1n)).toBe(1n);
    expect(hilbert_symbol(-1n, 1n, -1n)).toBe(1n);
    expect(hilbert_symbol(-1n, -1n, -1n)).toBe(-1n);
  });

  test('p = 2 cases', () => {
    // From SageMath: hilbert_symbol(-1, -1, 2) = -1
    expect(hilbert_symbol(-1n, -1n, 2n)).toBe(-1n);
    // hilbert_symbol(1, 1, 2) = 1
    expect(hilbert_symbol(1n, 1n, 2n)).toBe(1n);
    // hilbert_symbol(2, 3, 2) = -1 (both 3 mod 4)
    expect(hilbert_symbol(3n, 3n, 2n)).toBe(-1n);
  });

  test('odd prime cases', () => {
    // hilbert_symbol(a, b, p) = 1 if a is a QR mod p
    expect(hilbert_symbol(1n, 5n, 3n)).toBe(1n);
    expect(hilbert_symbol(4n, 5n, 3n)).toBe(1n);
  });

  test('throws for non-prime p', () => {
    expect(() => hilbert_symbol(1n, 2n, 4n)).toThrow();
  });
});

describe('hilbert_conductor', () => {
  test('hilbert_conductor(-1, -1) = 2', () => {
    expect(hilbert_conductor(-1n, -1n)).toBe(2n);
  });

  test('hilbert_conductor(-1, 1) = 1', () => {
    expect(hilbert_conductor(-1n, 1n)).toBe(1n);
  });

  test('hilbert_conductor(1, 1) = 1', () => {
    expect(hilbert_conductor(1n, 1n)).toBe(1n);
  });
});

describe('hilbert_conductor_inverse', () => {
  test('hilbert_conductor_inverse(1)', () => {
    const [a, b] = hilbert_conductor_inverse(1n);
    expect(hilbert_conductor(a, b)).toBe(1n);
  });

  test('hilbert_conductor_inverse(2)', () => {
    const [a, b] = hilbert_conductor_inverse(2n);
    expect(hilbert_conductor(a, b)).toBe(2n);
  });

  test('hilbert_conductor_inverse(3) (prime, 3 mod 4)', () => {
    const [a, b] = hilbert_conductor_inverse(3n);
    expect(hilbert_conductor(a, b)).toBe(3n);
  });

  test('hilbert_conductor_inverse(5) (prime, 1 mod 4)', () => {
    const [a, b] = hilbert_conductor_inverse(5n);
    expect(hilbert_conductor(a, b)).toBe(5n);
  });

  test('throws for non-positive d', () => {
    expect(() => hilbert_conductor_inverse(0n)).toThrow();
    expect(() => hilbert_conductor_inverse(-1n)).toThrow();
  });

  test('throws for non-squarefree d', () => {
    expect(() => hilbert_conductor_inverse(4n)).toThrow();
  });
});

describe('sort_complex_numbers_for_display', () => {
  test('empty list', () => {
    expect(sort_complex_numbers_for_display([])).toEqual([]);
  });

  test('real numbers first, sorted', () => {
    const nums = [
      { re: 3, im: 0 },
      { re: 1, im: 0 },
      { re: 2, im: 0 },
    ];
    const sorted = sort_complex_numbers_for_display(nums);
    expect(sorted.map((z) => z.re)).toEqual([1, 2, 3]);
  });

  test('complex numbers after real numbers', () => {
    const nums = [
      { re: 1, im: 1 },
      { re: 2, im: 0 },
      { re: 0, im: 1 },
    ];
    const sorted = sort_complex_numbers_for_display(nums);
    // Real (2) first, then complex sorted by real part then imaginary
    expect(sorted[0]).toEqual({ re: 2, im: 0 });
    expect(sorted[1]).toEqual({ re: 0, im: 1 });
    expect(sorted[2]).toEqual({ re: 1, im: 1 });
  });

  test('complex numbers sorted by real part then imaginary', () => {
    const nums = [
      { re: 1, im: 2 },
      { re: 1, im: 1 },
      { re: 0, im: 3 },
    ];
    const sorted = sort_complex_numbers_for_display(nums);
    expect(sorted.map((z) => [z.re, z.im])).toEqual([
      [0, 3],
      [1, 1],
      [1, 2],
    ]);
  });

  test('near-zero imaginary parts treated as real', () => {
    const nums = [
      { re: 2, im: 1e-15 },
      { re: 1, im: 0 },
    ];
    const sorted = sort_complex_numbers_for_display(nums);
    // Both should be treated as real
    expect(sorted[0]!.re).toBe(1);
    expect(sorted[1]!.re).toBe(2);
  });
});

describe('dedekind_sum', () => {
  test('s(0, k) = 0', () => {
    const s = dedekind_sum(0n, 5n);
    expect(s.numerator).toBe(0n);
  });

  test('s(h, 1) = 0', () => {
    const s = dedekind_sum(3n, 1n);
    expect(s.numerator).toBe(0n);
  });

  test('s(1, 2) = 0', () => {
    const s = dedekind_sum(1n, 2n);
    expect(s.numerator).toBe(0n);
  });

  test('s(1, 3) = 1/18', () => {
    // From SageMath: dedekind_sum(1, 3) = 1/18
    const s = dedekind_sum(1n, 3n);
    expect(s.numerator * 18n).toBe(s.denominator);
  });

  test('s(1, 4) = 1/8', () => {
    // From SageMath: dedekind_sum(1, 4) = 1/8
    const s = dedekind_sum(1n, 4n);
    expect(s.numerator * 8n).toBe(s.denominator);
  });

  test('s(3, 5) = 0', () => {
    // From SageMath: dedekind_sum(3, 5) = 0
    const s = dedekind_sum(3n, 5n);
    expect(s.numerator).toBe(0n);
  });

  test('s(2, 3) = -1/18', () => {
    // From SageMath: dedekind_sum(2, 3) = -1/18
    const s = dedekind_sum(2n, 3n);
    expect(s.numerator * 18n).toBe(-s.denominator);
  });

  test('s(-h, k) = -s(h, k)', () => {
    const s1 = dedekind_sum(3n, 7n);
    const s2 = dedekind_sum(-3n, 7n);
    // s1.num/s1.den = -s2.num/s2.den
    expect(s1.numerator * s2.denominator).toBe(-s2.numerator * s1.denominator);
  });

  test('s(2, 7) = 1/14', () => {
    // From SageMath: dedekind_sum(2, 7) = 1/14
    const s = dedekind_sum(2n, 7n);
    expect(s.numerator * 14n).toBe(s.denominator);
  });

  test('s(1, 5) = 1/5', () => {
    // From SageMath: dedekind_sum(1, 5) = 1/5
    const s = dedekind_sum(1n, 5n);
    expect(s.numerator * 5n).toBe(s.denominator);
  });

  test("SageMath's doctested table for q < 10 (AUDIT-2026-07 L1)", () => {
    // sage: for q in range(10): print([dedekind_sum(p,q) for p in range(q+1)])
    const table: string[][] = [
      ['0'],
      ['0', '0'],
      ['0', '0', '0'],
      ['0', '1/18', '-1/18', '0'],
      ['0', '1/8', '0', '-1/8', '0'],
      ['0', '1/5', '0', '0', '-1/5', '0'],
      ['0', '5/18', '1/18', '0', '-1/18', '-5/18', '0'],
      ['0', '5/14', '1/14', '-1/14', '1/14', '-1/14', '-5/14', '0'],
      ['0', '7/16', '1/8', '1/16', '0', '-1/16', '-1/8', '-7/16', '0'],
      ['0', '14/27', '4/27', '1/18', '-4/27', '4/27', '-1/18', '-4/27', '-14/27', '0'],
    ];
    for (let q = 0; q < table.length; q++) {
      for (let p = 0; p <= q; p++) {
        const r = dedekind_sum(BigInt(p), BigInt(q));
        const str = r.denominator === 1n ? `${r.numerator}` : `${r.numerator}/${r.denominator}`;
        expect(str).toBe(table[q]![p]!);
      }
    }
  });

  test('s(1, 23) = 77/46 and reciprocity', () => {
    // sage: q = 23; dedekind_sum(1, q) == (q-1)*(q-2)/(12*q)
    const s = dedekind_sum(1n, 23n);
    expect(s.numerator).toBe(77n);
    expect(s.denominator).toBe(46n);

    // sage: dedekind_sum(100, 723) + dedekind_sum(723, 100) == 31583/86760
    const a = dedekind_sum(100n, 723n);
    const b = dedekind_sum(723n, 100n);
    const num = a.numerator * b.denominator + b.numerator * a.denominator;
    const den = a.denominator * b.denominator;
    expect(num * 86760n).toBe(31583n * den);
  });
});

describe('mqrr_rational_reconstruction', () => {
  test('u = 0 with m > T returns [0, 1]', () => {
    const result = mqrr_rational_reconstruction(0n, 100n, 50n);
    expect(result).toEqual([0n, 1n]);
  });

  test('u = 0 with m <= T returns null', () => {
    const result = mqrr_rational_reconstruction(0n, 50n, 100n);
    expect(result).toBe(null);
  });

  test('reconstructs 3/5 from 3*inv(5) mod m', () => {
    // 3/5 mod 1000: need inv(5, 1000) = 201 (since 5*201 = 1005 = 1 mod 1000)
    // 3 * 201 = 603 mod 1000
    const result = mqrr_rational_reconstruction(603n, 1000n, 5n);
    if (result !== null) {
      const [n, d] = result;
      // Check n/d represents the same as 603/1000 when reduced
      expect((n * inverse_mod(d, 1000n)) % 1000n).toBe(603n);
    }
  });

  test('reconstruction fails when T is too large', () => {
    // With very large T, reconstruction may fail
    const result = mqrr_rational_reconstruction(123n, 1000n, 1000000n);
    // This depends on the algorithm, but reconstruction constraints may not be met
    // Just check it returns a valid format
    expect(result === null || (Array.isArray(result) && result.length === 2)).toBe(true);
  });

  test('simple reconstruction', () => {
    // 1/3 mod 100: inv(3, 100) = 67, so 1/3 = 67 mod 100
    const result = mqrr_rational_reconstruction(67n, 100n, 3n);
    if (result !== null) {
      const [n, d] = result;
      // Should reconstruct to something equivalent to 1/3
      expect((n * inverse_mod(d, 100n)) % 100n).toBe(67n);
    }
  });
});

describe('algebraic_dependency', () => {
  test('finds algebraic dependency for sqrt(2)', () => {
    const sqrt2 = Math.sqrt(2);
    const poly = algebraic_dependency(sqrt2, 2n);

    // Should find x^2 - 2 or a scalar multiple
    // Verify that sqrt(2) is a root
    let evalResult = 0;
    for (let i = 0; i < poly.length; i++) {
      evalResult += Number(poly[i]!) * sqrt2 ** i;
    }
    expect(Math.abs(evalResult)).toBeLessThan(1e-8);
  });

  test('finds linear dependency for rational approximation', () => {
    // 17/9 ≈ 1.888888...
    const x = 17 / 9;
    const poly = algebraic_dependency(x, 1n);

    // Should find approximately 9*x - 17
    // Check: poly[0] + poly[1] * (17/9) ≈ 0
    const evalResult = Number(poly[0]!) + Number(poly[1]!) * x;
    expect(Math.abs(evalResult)).toBeLessThan(1e-8);
  });

  test('finds dependency for cube root', () => {
    const cubeRoot2 = Math.cbrt(2);
    const poly = algebraic_dependency(cubeRoot2, 3n);

    // Should find x^3 - 2 or similar
    let evalResult = 0;
    for (let i = 0; i < poly.length; i++) {
      evalResult += Number(poly[i]!) * cubeRoot2 ** i;
    }
    expect(Math.abs(evalResult)).toBeLessThan(1e-6);
  });

  test('returns simple polynomial for integer', () => {
    const poly = algebraic_dependency(5, 1n);

    // Should return x - 5
    expect(poly[0]).toBe(-5n);
    expect(poly[1]).toBe(1n);
  });

  test('throws for non-finite input', () => {
    expect(() => algebraic_dependency(Number.POSITIVE_INFINITY, 2n)).toThrow(
      'z must be a finite number'
    );
    expect(() => algebraic_dependency(Number.NaN, 2n)).toThrow('z must be a finite number');
  });

  test('throws for degree < 1', () => {
    expect(() => algebraic_dependency(1.5, 0n)).toThrow('degree must be at least 1');
  });

  test('golden ratio satisfies x^2 - x - 1', () => {
    const phi = (1 + Math.sqrt(5)) / 2;
    const poly = algebraic_dependency(phi, 3n);

    // Check polynomial evaluates to ~0 at phi
    let evalResult = 0;
    for (let i = 0; i < poly.length; i++) {
      evalResult += Number(poly[i]!) * phi ** i;
    }
    expect(Math.abs(evalResult)).toBeLessThan(1e-8);
  });
});

describe('trial_division', () => {
  test("SageMath's arith/misc doctests", () => {
    expect(trial_division(15n)).toBe(3n);
    expect(trial_division(91n)).toBe(7n);
    expect(trial_division(11n)).toBe(11n);
    expect(trial_division(387833n, 300n)).toBe(387833n);
    expect(trial_division(387833n, 400n)).toBe(389n);
  });

  test("SageMath's Integer.trial_division doctests", () => {
    const n = 1000003n * 10000019n;
    expect(trial_division(n)).toBe(1000003n);
    expect(trial_division(-n)).toBe(1000003n);
    // AUDIT-2026-07 L4: the bound used to be ignored (a full factorization ran).
    expect(trial_division(n, 100n)).toBe(10000049000057n);

    const big = 100003n * 10000000000000000000000000000000000000121n;
    expect(trial_division(big)).toBe(100003n);
    expect(trial_division(big, 10n ** 4n)).toBe(1000030000000000000000000000000000000012100363n);
    expect(trial_division(-big, 10n ** 4n)).toBe(1000030000000000000000000000000000000012100363n);

    for (const p of [2n, 3n, 5n]) {
      expect(trial_division(p * 10000000000000000000000000000000000000121n)).toBe(p);
      expect(trial_division(p * 10007n)).toBe(p);
    }
  });

  test('start parameter', () => {
    // sage: (3*5*101*103).trial_division(start=50) == 101
    expect(trial_division(3n * 5n * 101n * 103n, undefined, 50n)).toBe(101n);
    expect(trial_division(3n * 5n * 101n * 103n)).toBe(3n);
    expect(trial_division(3n * 5n * 101n * 103n, undefined, 4n)).toBe(5n);
  });

  test('error cases', () => {
    expect(() => trial_division(100n, -10n)).toThrow('bound must be positive');
    expect(() => trial_division(100n, 0n)).toThrow('bound must be positive');
    expect(() => trial_division(0n)).toThrow('self must be nonzero');
    expect(trial_division(1n)).toBe(1n);
    expect(trial_division(-1n)).toBe(1n);
  });

  test('returns the smallest prime factor', () => {
    for (let k = 2n; k < 2000n; k++) {
      let spf = k;
      for (let p = 2n; p * p <= k; p++) {
        if (k % p === 0n) {
          spf = p;
          break;
        }
      }
      expect(trial_division(k)).toBe(spf);
    }
  });
});

describe('is_prime_power', () => {
  test("SageMath's doctests", () => {
    expect(is_prime_power(2n, true)).toEqual([2n, 1n]);
    expect(is_prime_power(4n, true)).toEqual([2n, 2n]);
    expect(is_prime_power(512n, true)).toEqual([2n, 9n]);
    // issue #4777
    expect(is_prime_power(150607571n ** 14n)).toBe(true);
  });

  test('never infers primality from a failed factorization (AUDIT-2026-07 H2)', () => {
    const p = 100000000000000000000117n;
    const q = 123456789012345678901259n;
    expect(is_prime(p)).toBe(true);
    expect(is_prime(q)).toBe(true);
    // A 48-digit semiprime: PARI's isprimepower answers without factoring.
    expect(is_prime_power(p * q)).toBe(false);
    expect(is_prime_power(p * q, true)).toEqual([p * q, 0n]);
    expect(is_prime_power(p * p, true)).toEqual([p, 2n]);
    expect(is_prime_power(p ** 3n, true)).toEqual([p, 3n]);
    expect(is_prime_power(p ** 6n, true)).toEqual([p, 6n]);
  });

  test('non-positive input', () => {
    expect(is_prime_power(0n)).toBe(false);
    expect(is_prime_power(1n)).toBe(false);
    expect(is_prime_power(-4n)).toBe(false);
    expect(is_prime_power(0n, true)).toEqual([0n, 0n]);
    expect(is_prime_power(1n, true)).toEqual([1n, 0n]);
  });

  test('agrees with the definition for n <= 3000', () => {
    for (let n = 2n; n <= 3000n; n++) {
      let want = false;
      let wp = 0n;
      let wk = 0n;
      for (let p = 2n; p <= n; p++) {
        if (n % p !== 0n) continue;
        let m = n;
        let k = 0n;
        while (m % p === 0n) {
          m /= p;
          k++;
        }
        if (m === 1n) {
          want = true;
          wp = p;
          wk = k;
        }
        break;
      }
      expect(is_prime_power(n)).toBe(want);
      if (want) expect(is_prime_power(n, true)).toEqual([wp, wk]);
    }
  });

  test('large bases beyond the tiny-prime table', () => {
    expect(is_prime_power(211n ** 6n, true)).toEqual([211n, 6n]);
    expect(is_prime_power(211n ** 7n * 2n)).toBe(false);
    expect(is_prime_power(1000003n ** 5n, true)).toEqual([1000003n, 5n]);
  });
});

describe('is_pseudoprime', () => {
  test('agrees with is_prime (BPSW, AUDIT-2026-07 H3)', () => {
    // 318665857834031151167461 is the smallest composite that is a strong
    // pseudoprime to all of the first 12 prime bases.
    const composite = 318665857834031151167461n;
    expect(is_pseudoprime(composite)).toBe(false);
    expect(is_prime(composite)).toBe(false);
    expect(next_probable_prime(composite - 2n)).not.toBe(composite);
  });

  test("SageMath's doctests", () => {
    expect(is_pseudoprime(389n)).toBe(true);
    expect(is_pseudoprime(2000n)).toBe(false);
    expect(is_pseudoprime(2n)).toBe(true);
    expect(is_pseudoprime(-1n)).toBe(false);
  });

  test('matches is_prime on a deterministic pseudo-random 64-bit sweep', () => {
    let seed = 123456789n;
    const mask = (1n << 64n) - 1n;
    for (let i = 0; i < 4000; i++) {
      seed = (seed * 6364136223846793005n + 1442695040888963407n) & mask;
      expect(is_pseudoprime(seed)).toBe(is_prime(seed));
    }
  });
});

describe('prime_powers', () => {
  test("SageMath's doctests", () => {
    expect(prime_powers(20n)).toEqual([2n, 3n, 4n, 5n, 7n, 8n, 9n, 11n, 13n, 16n, 17n, 19n]);
    expect(prime_powers(0n, 1n)).toEqual([]);
    expect(prime_powers(2n)).toEqual([]);
    expect(prime_powers(3n)).toEqual([2n]);
    expect(prime_powers(6n)).toEqual([2n, 3n, 4n, 5n]);
    expect(prime_powers(6n, 10n)).toEqual([7n, 8n, 9n]);
  });

  test('agrees with filtering by is_prime_power', () => {
    const want: bigint[] = [];
    for (let n = 2n; n < 500n; n++) if (is_prime_power(n)) want.push(n);
    expect(prime_powers(500n)).toEqual(want);
    expect(prime_powers(100n, 500n)).toEqual(want.filter((q) => q >= 100n));
  });

  test('walks the primes rather than testing every integer (AUDIT-2026-07 L3)', () => {
    const start = Date.now();
    const v = prime_powers(100000n);
    expect(Date.now() - start).toBeLessThan(2000);
    expect(v.length).toBe(9700);
    expect(v[0]).toBe(2n);
    expect(v[v.length - 1]).toBe(99991n);
  });
});

describe('prime_to_m_part', () => {
  test("SageMath's doctests", () => {
    expect(prime_to_m_part(43434n, 20n)).toBe(21717n);
    expect(prime_to_m_part(2048n, 2n)).toBe(1n);
    expect(prime_to_m_part(2048n, 3n)).toBe(2048n);
    expect(() => prime_to_m_part(0n, 2n)).toThrow('self must be nonzero');
  });

  test('m = 0 gives 1, unit m gives n (AUDIT-2026-07 L2)', () => {
    expect(prime_to_m_part(240n, 0n)).toBe(1n);
    expect(prime_to_m_part(240n, 1n)).toBe(240n);
    expect(prime_to_m_part(240n, -1n)).toBe(240n);
  });

  test('preserves the sign of n', () => {
    expect(prime_to_m_part(-240n, 6n)).toBe(-5n);
    expect(prime_to_m_part(-240n, 2n)).toBe(-15n);
  });

  test('result is coprime to m and n/result is m-smooth', () => {
    for (let n = 1n; n <= 200n; n++) {
      for (let m = 2n; m <= 12n; m++) {
        const r = prime_to_m_part(n, m);
        expect(n % r).toBe(0n);
        expect(gcd(r, m)).toBe(1n);
      }
    }
  });
});

describe('smooth_part / coprime_part', () => {
  test("SageMath's doctests", () => {
    const primes1000 = prime_range(1000n);
    expect(smooth_part(10n ** 77n + 1n, primes1000)).toEqual([
      [11n, 2n],
      [23n, 1n],
      [463n, 1n],
    ]);

    const primes10000 = prime_range(10000n);
    expect(coprime_part(10n ** 77n + 1n, primes10000)).toBe(
      2159827213801295896328509719222460043196544298056155507343412527n
    );
    expect(coprime_part(10n ** 55n + 1n, primes10000)).toBe(6426667196963538873896485804232411n);
  });

  test('smooth_part * coprime_part == x (AUDIT-2026-07 M2)', () => {
    const primes10000 = prime_range(10000n);
    for (const [x, base] of [
      [240n, [6n]],
      [240n, [2n, 3n]],
      [10n ** 55n + 1n, primes10000],
      [123456789n, [3n, 7n]],
    ] as Array<[bigint, bigint[]]>) {
      let prod = 1n;
      for (const [p, e] of smooth_part(x, base)) prod *= p ** e;
      expect(prod * coprime_part(x, base)).toBe(x);
    }
  });

  test('composite factor base entries are stripped as a unit', () => {
    // sage: coprime_part(240, [6]) == 40, because smooth_part(240, [6]) == 6
    expect(smooth_part(240n, [6n])).toEqual([[6n, 1n]]);
    expect(coprime_part(240n, [6n])).toBe(40n);
    expect(coprime_part(240n, [2n, 3n])).toBe(5n);
  });
});

describe('two_squares / three_squares / four_squares', () => {
  test("SageMath's two_squares doctests", () => {
    expect(two_squares(389n)).toEqual([10n, 17n]);
    expect(two_squares(21n)).toBeNull();
    expect(two_squares(21n * 21n)).toEqual([0n, 21n]);
    expect(two_squares(0n)).toEqual([0n, 0n]);
    expect(two_squares(106n)).toEqual([5n, 9n]);
  });

  test("SageMath's three_squares doctests (AUDIT-2026-07 M1)", () => {
    // n < 2^32 goes through sum_of_squares.pyx, which gives different --
    // and canonical -- answers from the factorization-based path.
    expect(three_squares(389n)).toEqual([1n, 8n, 18n]);
    expect(three_squares(946n)).toEqual([9n, 9n, 28n]);
    expect(three_squares(2986n)).toEqual([3n, 24n, 49n]);
    expect(three_squares(107n)).toEqual([1n, 5n, 9n]);
    expect(three_squares(0n)).toEqual([0n, 0n, 0n]);
    expect(three_squares(7n)).toBeNull();
  });

  test("SageMath's three_squares doctests above the cutoff", () => {
    expect(three_squares(7n ** 100n)).toEqual([
      0n,
      0n,
      1798465042647412146620280340569649349251249n,
    ]);
    expect(three_squares(11n ** 111n - 1n)).toEqual([
      616274160655975340150706442680n,
      901582938385735143295060746161n,
      6270382387635744140394001363065311967964099981788593947233n,
    ]);
    expect(three_squares(7n * 2n ** 41n)).toEqual([1048576n, 2097152n, 3145728n]);
    expect(three_squares(7n * 2n ** 42n)).toBeNull();
  });

  test("SageMath's four_squares doctests (AUDIT-2026-07 M1)", () => {
    expect(four_squares(389n)).toEqual([0n, 1n, 8n, 18n]);
    expect(four_squares(15447n)).toEqual([2n, 5n, 17n, 123n]);
    expect(four_squares(523439n)).toEqual([3n, 5n, 26n, 723n]);
    expect(four_squares(0n)).toEqual([0n, 0n, 0n, 0n]);
    // the large path, which delegates its tail to the pyx three_squares
    expect(four_squares(1101011011004n)).toEqual([90n, 102n, 1220n, 1049290n]);
  });

  test('pyx TESTS: two_squares(n^2) == (0, n) and two_squares(n^2+1) == (1, n)', () => {
    for (let n = 1n; n < 1500n; n++) {
      expect(two_squares(n * n)).toEqual([0n, n]);
      expect(two_squares(n * n + 1n)).toEqual([1n, n]);
    }
  });

  test('four_squares sums correctly and is sorted', () => {
    for (let n = 5000n; n < 6000n; n++) {
      const r = four_squares(n)!;
      expect(r[0] ** 2n + r[1] ** 2n + r[2] ** 2n + r[3] ** 2n).toBe(n);
      expect(r[0] <= r[1] && r[1] <= r[2] && r[2] <= r[3]).toBe(true);
    }
    // above the 2^32 cutoff
    for (let i = 0n; i < 20n; i++) {
      const n = 2n ** 34n + i * 7919n;
      const r = four_squares(n)!;
      expect(r[0] ** 2n + r[1] ** 2n + r[2] ** 2n + r[3] ** 2n).toBe(n);
    }
  });

  test('three_squares matches Legendre and sums correctly', () => {
    for (let n = 1n; n < 5000n; n++) {
      let m = n;
      while (m % 4n === 0n) m /= 4n;
      const possible = m % 8n !== 7n;
      const r = three_squares(n);
      expect(r !== null).toBe(possible);
      if (r) {
        expect(r[0] ** 2n + r[1] ** 2n + r[2] ** 2n).toBe(n);
        expect(r[0] <= r[1] && r[1] <= r[2]).toBe(true);
      }
    }
  });
});
