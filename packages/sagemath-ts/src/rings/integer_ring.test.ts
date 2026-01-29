/**
 * Tests for integer_ring.ts (Integer class methods)
 */

import { describe, test, expect } from 'bun:test';
import { Integer, ZZ } from './integer_ring.js';

describe('Integer', () => {
  describe('nth_root', () => {
    test('perfect cubes', () => {
      expect(new Integer(125n).nth_root(3n).value).toBe(5n);
      expect(new Integer(8n).nth_root(3n).value).toBe(2n);
      expect(new Integer(1n).nth_root(3n).value).toBe(1n);
      expect(new Integer(0n).nth_root(3n).value).toBe(0n);
    });

    test('perfect squares', () => {
      expect(new Integer(100n).nth_root(2n).value).toBe(10n);
      expect(new Integer(144n).nth_root(2n).value).toBe(12n);
    });

    test('negative cubes', () => {
      expect(new Integer(-125n).nth_root(3n).value).toBe(-5n);
      expect(new Integer(-8n).nth_root(3n).value).toBe(-2n);
    });

    test('throws for non-perfect powers', () => {
      expect(() => new Integer(124n).nth_root(3n)).toThrow('is not a 3rd power');
      expect(() => new Integer(10n).nth_root(2n)).toThrow('is not a 2nd power');
    });

    test('throws for even root of negative', () => {
      expect(() => new Integer(-4n).nth_root(2n)).toThrow('cannot take even root');
    });

    test('truncate mode', () => {
      const [root, exact] = new Integer(124n).nth_root(3n, true);
      expect(root.value).toBe(4n);
      expect(exact).toBe(false);

      const [root2, exact2] = new Integer(125n).nth_root(3n, true);
      expect(root2.value).toBe(5n);
      expect(exact2).toBe(true);
    });
  });

  describe('exact_log', () => {
    test('powers of 2', () => {
      expect(new Integer(1024n).exact_log(2n)).toBe(10n);
      expect(new Integer(256n).exact_log(2n)).toBe(8n);
      expect(new Integer(1n).exact_log(2n)).toBe(0n);
    });

    test('powers of 10', () => {
      expect(new Integer(1000n).exact_log(10n)).toBe(3n);
      expect(new Integer(10000n).exact_log(10n)).toBe(4n);
    });

    test('non-exact logs', () => {
      expect(new Integer(1023n).exact_log(2n)).toBe(9n); // floor(log2(1023)) = 9
      expect(new Integer(1025n).exact_log(2n)).toBe(10n); // floor(log2(1025)) = 10
    });

    test('throws for invalid base', () => {
      expect(() => new Integer(10n).exact_log(1n)).toThrow('base must be >= 2');
    });

    test('throws for non-positive self', () => {
      expect(() => new Integer(0n).exact_log(2n)).toThrow('self must be positive');
      expect(() => new Integer(-5n).exact_log(2n)).toThrow('self must be positive');
    });
  });

  describe('is_perfect_power', () => {
    test('perfect squares', () => {
      expect(new Integer(4n).is_perfect_power()).toBe(true);
      expect(new Integer(9n).is_perfect_power()).toBe(true);
      expect(new Integer(144n).is_perfect_power()).toBe(true);
    });

    test('perfect cubes', () => {
      expect(new Integer(8n).is_perfect_power()).toBe(true);
      expect(new Integer(27n).is_perfect_power()).toBe(true);
      expect(new Integer(125n).is_perfect_power()).toBe(true);
    });

    test('negative perfect powers', () => {
      expect(new Integer(-8n).is_perfect_power()).toBe(true);
      expect(new Integer(-27n).is_perfect_power()).toBe(true);
    });

    test('non-perfect powers', () => {
      expect(new Integer(10n).is_perfect_power()).toBe(false);
      expect(new Integer(12n).is_perfect_power()).toBe(false);
      expect(new Integer(-4n).is_perfect_power()).toBe(false);
    });

    test('special cases', () => {
      expect(new Integer(0n).is_perfect_power()).toBe(true);
      expect(new Integer(1n).is_perfect_power()).toBe(true);
      expect(new Integer(-1n).is_perfect_power()).toBe(true);
    });
  });

  describe('is_prime_power', () => {
    test('primes', () => {
      expect(new Integer(2n).is_prime_power()).toBe(true);
      expect(new Integer(3n).is_prime_power()).toBe(true);
      expect(new Integer(17n).is_prime_power()).toBe(true);
    });

    test('prime powers', () => {
      expect(new Integer(4n).is_prime_power()).toBe(true); // 2^2
      expect(new Integer(8n).is_prime_power()).toBe(true); // 2^3
      expect(new Integer(27n).is_prime_power()).toBe(true); // 3^3
      expect(new Integer(125n).is_prime_power()).toBe(true); // 5^3
    });

    test('non-prime powers', () => {
      expect(new Integer(1n).is_prime_power()).toBe(false);
      expect(new Integer(6n).is_prime_power()).toBe(false); // 2*3
      expect(new Integer(12n).is_prime_power()).toBe(false); // 2^2*3
      expect(new Integer(-7n).is_prime_power()).toBe(false); // negative
    });
  });

  describe('binomial', () => {
    test('basic binomial coefficients', () => {
      expect(new Integer(5n).binomial(2n).value).toBe(10n);
      expect(new Integer(10n).binomial(3n).value).toBe(120n);
      expect(new Integer(20n).binomial(10n).value).toBe(184756n);
    });

    test('edge cases', () => {
      expect(new Integer(5n).binomial(0n).value).toBe(1n);
      expect(new Integer(5n).binomial(5n).value).toBe(1n);
      expect(new Integer(5n).binomial(6n).value).toBe(0n);
      expect(new Integer(5n).binomial(-1n).value).toBe(0n);
    });

    test('negative n', () => {
      expect(new Integer(-2n).binomial(3n).value).toBe(-4n);
    });
  });

  describe('factorial', () => {
    test('basic factorials', () => {
      expect(new Integer(0n).factorial().value).toBe(1n);
      expect(new Integer(1n).factorial().value).toBe(1n);
      expect(new Integer(5n).factorial().value).toBe(120n);
      expect(new Integer(10n).factorial().value).toBe(3628800n);
    });

    test('throws for negative', () => {
      expect(() => new Integer(-1n).factorial()).toThrow();
    });
  });

  describe('euler_phi', () => {
    test('basic euler phi', () => {
      expect(new Integer(1n).euler_phi().value).toBe(1n);
      expect(new Integer(12n).euler_phi().value).toBe(4n);
      expect(new Integer(10n).euler_phi().value).toBe(4n);
      expect(new Integer(7n).euler_phi().value).toBe(6n); // prime
    });
  });

  describe('sigma', () => {
    test('sum of divisors', () => {
      expect(new Integer(12n).sigma(1n).value).toBe(28n); // 1+2+3+4+6+12 = 28
      expect(new Integer(6n).sigma(1n).value).toBe(12n); // 1+2+3+6 = 12
    });

    test('number of divisors', () => {
      expect(new Integer(12n).sigma(0n).value).toBe(6n); // divisors: 1,2,3,4,6,12
    });
  });

  describe('moebius', () => {
    test('moebius function', () => {
      expect(new Integer(1n).moebius()).toBe(1n);
      expect(new Integer(2n).moebius()).toBe(-1n);
      expect(new Integer(6n).moebius()).toBe(1n); // 2*3
      expect(new Integer(4n).moebius()).toBe(0n); // 2^2
      expect(new Integer(30n).moebius()).toBe(-1n); // 2*3*5
    });
  });

  describe('radical', () => {
    test('basic radical', () => {
      expect(new Integer(12n).radical().value).toBe(6n); // 2*3
      expect(new Integer(100n).radical().value).toBe(10n); // 2*5
      expect(new Integer(72n).radical().value).toBe(6n); // 2*3
    });
  });

  describe('prime_divisors', () => {
    test('basic prime divisors', () => {
      const divs = new Integer(12n).prime_divisors().map(d => d.value);
      expect(divs).toEqual([2n, 3n]);
    });

    test('prime has itself as only divisor', () => {
      const divs = new Integer(17n).prime_divisors().map(d => d.value);
      expect(divs).toEqual([17n]);
    });
  });

  describe('divisors', () => {
    test('basic divisors', () => {
      const divs = new Integer(12n).divisors().map(d => d.value);
      expect(divs.sort((a, b) => Number(a - b))).toEqual([1n, 2n, 3n, 4n, 6n, 12n]);
    });
  });

  describe('squarefree_part', () => {
    test('basic squarefree part', () => {
      expect(new Integer(12n).squarefree_part().value).toBe(3n); // 12 = 3*4
      expect(new Integer(72n).squarefree_part().value).toBe(2n); // 72 = 2*36
      expect(new Integer(100n).squarefree_part().value).toBe(1n); // 100 = 1*100
    });
  });

  describe('jacobi and kronecker', () => {
    test('jacobi symbol', () => {
      expect(new Integer(2n).jacobi(15n)).toBe(1n);
      expect(new Integer(7n).jacobi(15n)).toBe(-1n);
    });

    test('kronecker symbol', () => {
      expect(new Integer(5n).kronecker(3n)).toBe(-1n);
      expect(new Integer(4n).kronecker(3n)).toBe(1n);
    });
  });

  describe('is_squarefree', () => {
    test('squarefree numbers', () => {
      expect(new Integer(6n).is_squarefree()).toBe(true);
      expect(new Integer(10n).is_squarefree()).toBe(true);
      expect(new Integer(1n).is_squarefree()).toBe(true);
    });

    test('non-squarefree numbers', () => {
      expect(new Integer(4n).is_squarefree()).toBe(false);
      expect(new Integer(12n).is_squarefree()).toBe(false);
      expect(new Integer(0n).is_squarefree()).toBe(false);
    });
  });

  describe('is_irreducible', () => {
    test('primes are irreducible', () => {
      expect(new Integer(7n).is_irreducible()).toBe(true);
      expect(new Integer(-7n).is_irreducible()).toBe(true);
    });

    test('composites are not irreducible', () => {
      expect(new Integer(6n).is_irreducible()).toBe(false);
      expect(new Integer(1n).is_irreducible()).toBe(false);
    });
  });

  describe('is_discriminant', () => {
    test('valid discriminants', () => {
      expect(new Integer(-3n).is_discriminant()).toBe(true); // -3 mod 4 = 1
      expect(new Integer(-4n).is_discriminant()).toBe(true); // -4 mod 4 = 0
      expect(new Integer(5n).is_discriminant()).toBe(true); // 5 mod 4 = 1
    });

    test('invalid discriminants', () => {
      expect(new Integer(2n).is_discriminant()).toBe(false); // 2 mod 4 = 2
      expect(new Integer(3n).is_discriminant()).toBe(false); // 3 mod 4 = 3
      expect(new Integer(4n).is_discriminant()).toBe(false); // perfect square
    });
  });

  describe('inverse_mod and powermod', () => {
    test('inverse_mod', () => {
      expect(new Integer(3n).inverse_mod(7n).value).toBe(5n); // 3*5 = 15 = 1 mod 7
    });

    test('powermod', () => {
      expect(new Integer(2n).powermod(10n, 1000n).value).toBe(24n);
      expect(new Integer(3n).powermod(5n, 7n).value).toBe(5n); // 3^5 = 243 = 5 mod 7
    });
  });

  describe('multiplicative_order', () => {
    test('basic multiplicative order', () => {
      expect(new Integer(2n).multiplicative_order(7n)).toBe(3n); // 2^3 = 8 = 1 mod 7
      expect(new Integer(3n).multiplicative_order(7n)).toBe(6n); // 3 is primitive root mod 7
    });
  });

  describe('fibonacci and lucas', () => {
    test('fibonacci', () => {
      expect(new Integer(0n).fibonacci().value).toBe(0n);
      expect(new Integer(1n).fibonacci().value).toBe(1n);
      expect(new Integer(10n).fibonacci().value).toBe(55n);
    });

    test('lucas', () => {
      expect(new Integer(0n).lucas_number().value).toBe(2n);
      expect(new Integer(1n).lucas_number().value).toBe(1n);
      expect(new Integer(5n).lucas_number().value).toBe(11n);
    });
  });

  describe('bell_number and catalan_number', () => {
    test('bell numbers', () => {
      expect(new Integer(0n).bell_number().value).toBe(1n);
      expect(new Integer(1n).bell_number().value).toBe(1n);
      expect(new Integer(4n).bell_number().value).toBe(15n);
      expect(new Integer(5n).bell_number().value).toBe(52n);
    });

    test('catalan numbers', () => {
      expect(new Integer(0n).catalan_number().value).toBe(1n);
      expect(new Integer(1n).catalan_number().value).toBe(1n);
      expect(new Integer(4n).catalan_number().value).toBe(14n);
      expect(new Integer(5n).catalan_number().value).toBe(42n);
    });
  });

  describe('carmichael_lambda', () => {
    test('basic carmichael lambda', () => {
      expect(new Integer(1n).carmichael_lambda().value).toBe(1n);
      expect(new Integer(8n).carmichael_lambda().value).toBe(2n);
      expect(new Integer(12n).carmichael_lambda().value).toBe(2n);
    });
  });

  describe('sqrt_mod', () => {
    test('basic sqrt_mod', () => {
      const r = new Integer(2n).sqrt_mod(7n);
      expect(r).not.toBeNull();
      expect((r!.value * r!.value) % 7n).toBe(2n);
    });

    test('no square root', () => {
      const r = new Integer(3n).sqrt_mod(7n);
      expect(r).toBeNull();
    });
  });

  describe('sqrtrem', () => {
    test('basic sqrtrem', () => {
      const [s, r] = new Integer(10n).sqrtrem();
      expect(s.value).toBe(3n);
      expect(r.value).toBe(1n); // 10 = 9 + 1
    });

    test('perfect square', () => {
      const [s, r] = new Integer(16n).sqrtrem();
      expect(s.value).toBe(4n);
      expect(r.value).toBe(0n);
    });
  });

  describe('core', () => {
    test('squarefree core (t=2)', () => {
      expect(new Integer(72n).core(2n).value).toBe(2n); // 72 = 2 * 36
      expect(new Integer(12n).core(2n).value).toBe(3n); // 12 = 3 * 4
    });

    test('cubefree core (t=3)', () => {
      expect(new Integer(54n).core(3n).value).toBe(2n); // 54 = 2 * 27
    });
  });

  describe('number_of_partitions', () => {
    test('basic partition numbers', () => {
      expect(new Integer(0n).number_of_partitions().value).toBe(1n);
      expect(new Integer(1n).number_of_partitions().value).toBe(1n);
      expect(new Integer(5n).number_of_partitions().value).toBe(7n);
      expect(new Integer(10n).number_of_partitions().value).toBe(42n);
    });
  });
});

describe('Integer class_number', () => {
  test('known class numbers for imaginary quadratic discriminants', () => {
    expect(new Integer(-3n).class_number()).toBe(1n);
    expect(new Integer(-4n).class_number()).toBe(1n);
    expect(new Integer(-7n).class_number()).toBe(1n);
    expect(new Integer(-8n).class_number()).toBe(1n);
    expect(new Integer(-11n).class_number()).toBe(1n);
    expect(new Integer(-15n).class_number()).toBe(2n);
    expect(new Integer(-23n).class_number()).toBe(3n);
    expect(new Integer(-163n).class_number()).toBe(1n); // Heegner number
  });

  test('known class numbers for real quadratic discriminants', () => {
    expect(new Integer(5n).class_number()).toBe(1n);
    expect(new Integer(8n).class_number()).toBe(1n);
    expect(new Integer(12n).class_number()).toBe(1n);
    expect(new Integer(40n).class_number()).toBe(2n);
    expect(new Integer(65n).class_number()).toBe(2n);
  });

  test('throws for square integers', () => {
    expect(() => new Integer(4n).class_number()).toThrow('not defined for square integers');
    expect(() => new Integer(9n).class_number()).toThrow('not defined for square integers');
    expect(() => new Integer(0n).class_number()).toThrow('not defined for square integers');
  });

  test('throws for non-discriminants', () => {
    expect(() => new Integer(2n).class_number()).toThrow('congruent to 0 or 1 modulo 4');
    expect(() => new Integer(3n).class_number()).toThrow('congruent to 0 or 1 modulo 4');
    expect(() => new Integer(-2n).class_number()).toThrow('congruent to 0 or 1 modulo 4');
  });
});

describe('Integer nth_root_mod', () => {
  test('square roots (n=2)', () => {
    // 3^2 = 9 ≡ 2 (mod 7)
    const r1 = new Integer(2n).nth_root_mod(2n, 7n);
    expect((r1.value * r1.value) % 7n).toBe(2n);

    // 4^2 = 16 ≡ 2 (mod 7) as well
    const r2 = new Integer(4n).nth_root_mod(2n, 7n);
    expect((r2.value * r2.value) % 7n).toBe(4n);
  });

  test('cube roots (n=3)', () => {
    // Find x such that x^3 ≡ 8 (mod 13)
    const r = new Integer(8n).nth_root_mod(3n, 13n);
    expect(_power_mod(r.value, 3n, 13n)).toBe(8n % 13n);
  });

  test('fourth roots (n=4)', () => {
    // Find x such that x^4 ≡ 16 (mod 17)
    const r = new Integer(16n).nth_root_mod(4n, 17n);
    expect(_power_mod(r.value, 4n, 17n)).toBe(16n % 17n);
  });

  test('trivial cases', () => {
    expect(new Integer(0n).nth_root_mod(3n, 7n).value).toBe(0n);
    expect(new Integer(1n).nth_root_mod(5n, 11n).value).toBe(1n);
  });

  test('throws when no root exists', () => {
    // 3 is not a quadratic residue mod 7
    expect(() => new Integer(3n).nth_root_mod(2n, 7n)).toThrow('no n-th root');
  });

  test('n=1 returns the value', () => {
    expect(new Integer(5n).nth_root_mod(1n, 7n).value).toBe(5n);
    expect(new Integer(10n).nth_root_mod(1n, 11n).value).toBe(10n);
  });
});

// Import power_mod for testing
import { power_mod as _power_mod } from '../arith/misc.js';

describe('Integer __invert__', () => {
  test('units are invertible', () => {
    expect(new Integer(1n).__invert__().value).toBe(1n);
    expect(new Integer(-1n).__invert__().value).toBe(-1n);
  });

  test('non-units throw', () => {
    expect(() => new Integer(2n).__invert__()).toThrow('not invertible');
    expect(() => new Integer(0n).__invert__()).toThrow('not invertible');
    expect(() => new Integer(-5n).__invert__()).toThrow('not invertible');
  });
});

describe('IntegerRing ZZ', () => {
  test('coercion', () => {
    expect(ZZ.__call__(5)).toBe(5n);
    expect(ZZ.__call__(5n)).toBe(5n);
    expect(ZZ.__call__('123')).toBe(123n);
  });

  test('zero and one', () => {
    expect(ZZ.zero()).toBe(0n);
    expect(ZZ.one()).toBe(1n);
  });

  test('properties', () => {
    expect(ZZ.is_field()).toBe(false);
    expect(ZZ.is_ring()).toBe(true);
    expect(ZZ.is_integral_domain()).toBe(true);
    expect(ZZ.characteristic()).toBe(0n);
  });
});
