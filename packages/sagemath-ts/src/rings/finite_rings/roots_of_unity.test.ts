/**
 * Tests for roots of unity, multiplicative subgroups, and FFT domains
 */

import { describe, expect, it } from 'bun:test';
import { euler_phi } from '../../arith/misc.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import { GF } from './finite_field_extension.js';
import {
  CosetDomain,
  FFTDomain,
  createFFTDomain,
  cyclotomic_polynomial,
  elementOfOrder,
  findMultiplicativeGenerator,
  has_primitive_root,
  maxFFTSize,
  multiplicative_order,
  primitive_nth_root,
  roots_of_unity,
  twoAdicity,
  validFFTSizes,
} from './roots_of_unity.js';

describe('has_primitive_root', () => {
  it('should return true for divisors of p-1 in prime fields', () => {
    const F17 = GF(17);

    // 17 - 1 = 16 = 2^4
    expect(has_primitive_root(F17, 1n)).toBe(true);
    expect(has_primitive_root(F17, 2n)).toBe(true);
    expect(has_primitive_root(F17, 4n)).toBe(true);
    expect(has_primitive_root(F17, 8n)).toBe(true);
    expect(has_primitive_root(F17, 16n)).toBe(true);
  });

  it('should return false for non-divisors of p-1', () => {
    const F17 = GF(17);

    // 3, 5, 6, 7 do not divide 16
    expect(has_primitive_root(F17, 3n)).toBe(false);
    expect(has_primitive_root(F17, 5n)).toBe(false);
    expect(has_primitive_root(F17, 6n)).toBe(false);
    expect(has_primitive_root(F17, 7n)).toBe(false);
  });

  it('should work for various prime fields', () => {
    // F7: 7 - 1 = 6 = 2 * 3
    const F7 = GF(7);
    expect(has_primitive_root(F7, 1n)).toBe(true);
    expect(has_primitive_root(F7, 2n)).toBe(true);
    expect(has_primitive_root(F7, 3n)).toBe(true);
    expect(has_primitive_root(F7, 6n)).toBe(true);
    expect(has_primitive_root(F7, 4n)).toBe(false);

    // F13: 13 - 1 = 12 = 2^2 * 3
    const F13 = GF(13);
    expect(has_primitive_root(F13, 1n)).toBe(true);
    expect(has_primitive_root(F13, 2n)).toBe(true);
    expect(has_primitive_root(F13, 3n)).toBe(true);
    expect(has_primitive_root(F13, 4n)).toBe(true);
    expect(has_primitive_root(F13, 6n)).toBe(true);
    expect(has_primitive_root(F13, 12n)).toBe(true);
    expect(has_primitive_root(F13, 5n)).toBe(false);
  });

  it('should handle edge cases', () => {
    const F7 = GF(7);
    expect(has_primitive_root(F7, 0n)).toBe(false);
    expect(has_primitive_root(F7, -1n)).toBe(false);
  });
});

describe('primitive_nth_root', () => {
  it('should find primitive roots in GF(7)', () => {
    const F7 = GF(7);

    // 1st root is always 1
    const omega1 = primitive_nth_root(F7, 1n);
    expect(omega1.isOne()).toBe(true);

    // 2nd root: omega^2 = 1, omega != 1
    const omega2 = primitive_nth_root(F7, 2n);
    expect(omega2.pow(2n).isOne()).toBe(true);
    expect(omega2.isOne()).toBe(false);

    // 3rd root: omega^3 = 1, omega != 1, omega^2 != 1
    const omega3 = primitive_nth_root(F7, 3n);
    expect(omega3.pow(3n).isOne()).toBe(true);
    expect(omega3.isOne()).toBe(false);
    expect(omega3.pow(2n).isOne()).toBe(false);

    // 6th root (generator)
    const omega6 = primitive_nth_root(F7, 6n);
    expect(omega6.pow(6n).isOne()).toBe(true);
    for (let k = 1n; k < 6n; k++) {
      expect(omega6.pow(k).isOne()).toBe(false);
    }
  });

  it('should find primitive roots in GF(17)', () => {
    const F17 = GF(17);

    // 8th root
    const omega8 = primitive_nth_root(F17, 8n);
    expect(omega8.pow(8n).isOne()).toBe(true);
    expect(omega8.pow(4n).isOne()).toBe(false);
    expect(omega8.pow(2n).isOne()).toBe(false);
    expect(omega8.isOne()).toBe(false);

    // 16th root (generator)
    const omega16 = primitive_nth_root(F17, 16n);
    expect(omega16.pow(16n).isOne()).toBe(true);
    expect(omega16.pow(8n).isOne()).toBe(false);
  });

  it('should throw for non-existent roots', () => {
    const F7 = GF(7);

    expect(() => primitive_nth_root(F7, 4n)).toThrow();
    expect(() => primitive_nth_root(F7, 5n)).toThrow();
    expect(() => primitive_nth_root(F7, 0n)).toThrow();
  });

  it('should satisfy omega = g^((p-1)/n)', () => {
    const F17 = GF(17);
    const g = findMultiplicativeGenerator(F17);

    for (const n of [2n, 4n, 8n, 16n]) {
      const omega = primitive_nth_root(F17, n);
      const expected = g.pow((17n - 1n) / n);
      expect(omega.eq(expected)).toBe(true);
    }
  });
});

describe('roots_of_unity', () => {
  it('should return all n-th roots in GF(7)', () => {
    const F7 = GF(7);

    // 3rd roots of unity
    const roots3 = roots_of_unity(F7, 3n);
    expect(roots3.length).toBe(3);

    // All should satisfy x^3 = 1
    for (const r of roots3) {
      expect(r.pow(3n).isOne()).toBe(true);
    }

    // Should include 1
    expect(roots3[0]!.isOne()).toBe(true);

    // Should be distinct
    const values = roots3.map((r) => {
      if ('value' in r && typeof r.value === 'bigint') {
        return r.value;
      }
      return r.toString();
    });
    const uniqueValues = new Set(values);
    expect(uniqueValues.size).toBe(3);
  });

  it('should return all n-th roots in GF(13)', () => {
    const F13 = GF(13);

    // 4th roots of unity
    const roots4 = roots_of_unity(F13, 4n);
    expect(roots4.length).toBe(4);

    // All should satisfy x^4 = 1
    for (const r of roots4) {
      expect(r.pow(4n).isOne()).toBe(true);
    }

    // The roots should be 1, omega, omega^2, omega^3
    const omega = primitive_nth_root(F13, 4n);
    expect(roots4[0]!.isOne()).toBe(true);
    expect(roots4[1]!.eq(omega)).toBe(true);
    expect(roots4[2]!.eq(omega.pow(2n))).toBe(true);
    expect(roots4[3]!.eq(omega.pow(3n))).toBe(true);
  });

  it('should return [1] for n=1', () => {
    const F7 = GF(7);
    const roots1 = roots_of_unity(F7, 1n);
    expect(roots1.length).toBe(1);
    expect(roots1[0]!.isOne()).toBe(true);
  });
});

describe('multiplicative_order', () => {
  it('should compute order correctly in GF(7)', () => {
    const F7 = GF(7);

    // 1 has order 1
    expect(multiplicative_order(F7.__call__(1n))).toBe(1n);

    // -1 = 6 has order 2
    expect(multiplicative_order(F7.__call__(6n))).toBe(2n);

    // 2 has order 3 (2^3 = 8 = 1 mod 7)
    expect(multiplicative_order(F7.__call__(2n))).toBe(3n);

    // 4 has order 3 (4 = 2^2, so 4^3 = 2^6 = 1)
    expect(multiplicative_order(F7.__call__(4n))).toBe(3n);

    // 3 is a generator (order 6)
    expect(multiplicative_order(F7.__call__(3n))).toBe(6n);

    // 5 is also a generator
    expect(multiplicative_order(F7.__call__(5n))).toBe(6n);
  });

  it('should compute order correctly in GF(13)', () => {
    const F13 = GF(13);

    expect(multiplicative_order(F13.__call__(1n))).toBe(1n);
    expect(multiplicative_order(F13.__call__(12n))).toBe(2n); // -1

    // All orders divide 12
    for (let i = 1n; i < 13n; i++) {
      const order = multiplicative_order(F13.__call__(i));
      expect(12n % order).toBe(0n);
    }
  });

  it('should throw for zero', () => {
    const F7 = GF(7);
    // sage: k(0).multiplicative_order()
    // ArithmeticError: Multiplicative order of 0 not defined.
    expect(() => multiplicative_order(F7.__call__(0n))).toThrow(
      'Multiplicative order of 0 not defined.'
    );
  });

  it('should use one exponentiation per prime power, not per divisor (L34)', () => {
    // Goldilocks: p - 1 = 2^32 * 3 * 5 * 17 * 257 * 65537 has 1056 divisors,
    // so the old divisor enumeration needed up to 1056 exponentiations where
    // element_base.pyx:709-714 needs 6 plus a few squarings.
    const p = 18446744069414584321n;
    const F = GF(p);
    const start = Date.now();
    // 7 is a multiplicative generator of the Goldilocks field.
    expect(multiplicative_order(F.__call__(7n))).toBe(p - 1n);
    // the order of a proper power divides it properly
    expect(multiplicative_order(F.__call__(7n).pow((p - 1n) / 3n))).toBe(3n);
    expect(multiplicative_order(F.__call__(7n).pow((p - 1n) / 65537n))).toBe(65537n);
    expect(multiplicative_order(F.__call__(7n).pow(2n ** 32n))).toBe((p - 1n) / 2n ** 32n);
    expect(multiplicative_order(F.__call__(p - 1n))).toBe(2n);
    expect(Date.now() - start).toBeLessThan(20000);
  });
});

describe('findMultiplicativeGenerator', () => {
  it('should find generator in GF(7)', () => {
    const F7 = GF(7);
    const g = findMultiplicativeGenerator(F7);

    // g should have order 6
    expect(multiplicative_order(g)).toBe(6n);
    expect(g.pow(6n).isOne()).toBe(true);
  });

  it('should find generator in GF(17)', () => {
    const F17 = GF(17);
    const g = findMultiplicativeGenerator(F17);

    // g should have order 16
    expect(multiplicative_order(g)).toBe(16n);
  });

  it('should find generator in GF(101)', () => {
    const F101 = GF(101);
    const g = findMultiplicativeGenerator(F101);

    // g should have order 100
    expect(multiplicative_order(g)).toBe(100n);
  });
});

describe('cyclotomic_polynomial', () => {
  it('should compute Phi_1(x) = x - 1', () => {
    const coeffs = cyclotomic_polynomial(1) as number[];
    expect(coeffs).toEqual([-1, 1]);
  });

  it('should compute Phi_2(x) = x + 1', () => {
    const coeffs = cyclotomic_polynomial(2) as number[];
    expect(coeffs).toEqual([1, 1]);
  });

  it('should compute Phi_3(x) = x^2 + x + 1', () => {
    const coeffs = cyclotomic_polynomial(3) as number[];
    expect(coeffs).toEqual([1, 1, 1]);
  });

  it('should compute Phi_4(x) = x^2 + 1', () => {
    const coeffs = cyclotomic_polynomial(4) as number[];
    expect(coeffs).toEqual([1, 0, 1]);
  });

  it('should compute Phi_5(x) = x^4 + x^3 + x^2 + x + 1', () => {
    const coeffs = cyclotomic_polynomial(5) as number[];
    expect(coeffs).toEqual([1, 1, 1, 1, 1]);
  });

  it('should compute Phi_6(x) = x^2 - x + 1', () => {
    const coeffs = cyclotomic_polynomial(6) as number[];
    expect(coeffs).toEqual([1, -1, 1]);
  });

  it('should compute Phi_8(x) = x^4 + 1', () => {
    const coeffs = cyclotomic_polynomial(8) as number[];
    expect(coeffs).toEqual([1, 0, 0, 0, 1]);
  });

  it('should compute Phi_12(x) = x^4 - x^2 + 1', () => {
    const coeffs = cyclotomic_polynomial(12) as number[];
    expect(coeffs).toEqual([1, 0, -1, 0, 1]);
  });

  it('should have degree phi(n) for cyclotomic polynomial', () => {
    for (const n of [1n, 2n, 3n, 4n, 5n, 6n, 7n, 8n, 9n, 10n, 12n, 15n, 20n]) {
      const coeffs = cyclotomic_polynomial(n) as number[];
      const degree = coeffs.length - 1;
      expect(BigInt(degree)).toBe(euler_phi(n));
    }
  });

  it('should have primitive n-th roots as roots of cyclotomic polynomial', () => {
    // Test in GF(13) which has primitive 12th roots
    const F13 = GF(13);
    const polyRing = new PolynomialRing(F13, 'x');

    // Get Phi_3 and verify primitive 3rd roots are roots
    const phi3Coeffs = cyclotomic_polynomial(3) as number[];
    const phi3 = polyRing.__call__(phi3Coeffs.map((c) => F13.__call__(c)));

    const omega3 = primitive_nth_root(F13, 3n);
    // Evaluate phi3 at omega3
    const value = phi3.evaluate(omega3);
    expect(value.isZero()).toBe(true);

    // omega3^2 should also be a root
    const omega3_2 = omega3.pow(2n);
    const value2 = phi3.evaluate(omega3_2);
    expect(value2.isZero()).toBe(true);
  });
});

describe('FFTDomain', () => {
  describe('construction', () => {
    it('should create domain with valid power-of-2 size', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 8n);

      expect(domain.size).toBe(8n);
      expect(domain.field).toBe(F17);
    });

    it('should reject non-power-of-2 sizes', () => {
      const F17 = GF(17);
      expect(() => new FFTDomain(F17, 3n)).toThrow();
      expect(() => new FFTDomain(F17, 6n)).toThrow();
      expect(() => new FFTDomain(F17, 7n)).toThrow();
    });

    it('should reject sizes that do not divide |F*|', () => {
      const F7 = GF(7); // |F*| = 6
      expect(() => new FFTDomain(F7, 8n)).toThrow(); // 8 does not divide 6
    });

    it('should work with Fermat prime', () => {
      // 65537 = 2^16 + 1, so |F*| = 2^16
      const F = GF(65537);
      const domain = new FFTDomain(F, 65536n);

      expect(domain.size).toBe(65536n);
    });
  });

  describe('generator', () => {
    it('should have generator that is primitive n-th root', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 8n);

      // generator^8 = 1
      expect(domain.generator.pow(8n).isOne()).toBe(true);

      // generator^4 != 1
      expect(domain.generator.pow(4n).isOne()).toBe(false);
    });
  });

  describe('elements', () => {
    it('should return correct number of elements', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 8n);

      const elems = domain.elements();
      expect(elems.length).toBe(8);
    });

    it('should return [1, omega, omega^2, ..., omega^(n-1)]', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 4n);

      const elems = domain.elements();
      const omega = domain.generator;

      expect(elems[0]!.isOne()).toBe(true);
      expect(elems[1]!.eq(omega)).toBe(true);
      expect(elems[2]!.eq(omega.pow(2n))).toBe(true);
      expect(elems[3]!.eq(omega.pow(3n))).toBe(true);
    });

    it('should all satisfy x^n = 1', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 8n);

      for (const elem of domain.elements()) {
        expect(elem.pow(8n).isOne()).toBe(true);
      }
    });
  });

  describe('contains', () => {
    it('should return true for domain elements', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 4n);

      for (const elem of domain.elements()) {
        expect(domain.contains(elem)).toBe(true);
      }
    });

    it('should return false for non-domain elements', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 4n);

      // 2^4 = 16 != 1, so 2 is not a 4th root of unity
      const two = F17.__call__(2n);
      expect(domain.contains(two)).toBe(false);
    });
  });

  describe('coset', () => {
    it('should create coset with offset', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 4n);
      const offset = F17.__call__(3n);

      const coset = domain.coset(offset);
      expect(coset.size).toBe(4n);
      expect(coset.offset).toBe(offset);
    });

    it('should have coset elements be offset * domain', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 4n);
      const offset = F17.__call__(3n);

      const coset = domain.coset(offset);
      const cosetElems = coset.elements();
      const domainElems = domain.elements();

      for (let i = 0; i < 4; i++) {
        const expected = offset.mul(domainElems[i]!);
        expect(cosetElems[i]!.eq(expected)).toBe(true);
      }
    });
  });

  describe('iteration', () => {
    it('should be iterable', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 4n);

      const elems = [...domain];
      expect(elems.length).toBe(4);
    });
  });

  describe('log2Size', () => {
    it('should return log base 2 of size', () => {
      const F17 = GF(17);

      expect(new FFTDomain(F17, 1n).log2Size).toBe(0);
      expect(new FFTDomain(F17, 2n).log2Size).toBe(1);
      expect(new FFTDomain(F17, 4n).log2Size).toBe(2);
      expect(new FFTDomain(F17, 8n).log2Size).toBe(3);
      expect(new FFTDomain(F17, 16n).log2Size).toBe(4);
    });
  });
});

describe('CosetDomain', () => {
  describe('elements', () => {
    it('should return offset * [1, omega, omega^2, ...]', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 4n);
      const offset = F17.__call__(5n);

      const coset = new CosetDomain(F17, 4n, domain.generator, offset);
      const elems = coset.elements();

      expect(elems[0]!.eq(offset)).toBe(true);
      expect(elems[1]!.eq(offset.mul(domain.generator))).toBe(true);
    });
  });

  describe('fold', () => {
    it('should create domain of half size', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 8n);
      const offset = F17.__call__(3n);

      const coset = domain.coset(offset);
      const challenge = F17.__call__(7n);
      const folded = coset.fold(challenge);

      expect(folded.size).toBe(4n);
    });

    it('should have squared generator', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 8n);
      const offset = F17.__call__(3n);

      const coset = domain.coset(offset);
      const challenge = F17.__call__(7n);
      const folded = coset.fold(challenge);

      // Folded generator should be omega^2
      expect(folded.generator.eq(domain.generator.pow(2n))).toBe(true);
    });

    it('should throw for size 1 domain', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 1n);
      const offset = F17.__call__(3n);

      const coset = domain.coset(offset);
      expect(() => coset.fold(F17.__call__(7n))).toThrow();
    });
  });

  describe('contains', () => {
    it('should detect coset membership', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 4n);
      const offset = F17.__call__(3n);

      const coset = domain.coset(offset);

      for (const elem of coset.elements()) {
        expect(coset.contains(elem)).toBe(true);
      }
    });
  });
});

describe('maxFFTSize', () => {
  it('should return largest power of 2 dividing |F*|', () => {
    // F7: |F*| = 6 = 2 * 3, so max FFT size is 2
    expect(maxFFTSize(GF(7))).toBe(2n);

    // F13: |F*| = 12 = 4 * 3, so max FFT size is 4
    expect(maxFFTSize(GF(13))).toBe(4n);

    // F17: |F*| = 16 = 2^4, so max FFT size is 16
    expect(maxFFTSize(GF(17))).toBe(16n);

    // F97: |F*| = 96 = 32 * 3, so max FFT size is 32
    expect(maxFFTSize(GF(97))).toBe(32n);

    // F65537: |F*| = 65536 = 2^16
    expect(maxFFTSize(GF(65537))).toBe(65536n);
  });
});

describe('validFFTSizes', () => {
  it('should return all powers of 2 up to max', () => {
    const F17 = GF(17);
    const sizes = validFFTSizes(F17);

    expect(sizes).toEqual([1n, 2n, 4n, 8n, 16n]);
  });

  it('should work for fields with limited FFT support', () => {
    const F7 = GF(7);
    const sizes = validFFTSizes(F7);

    expect(sizes).toEqual([1n, 2n]);
  });
});

describe('twoAdicity', () => {
  it('should return log2 of max FFT size', () => {
    expect(twoAdicity(GF(7))).toBe(1); // 6 = 2^1 * 3
    expect(twoAdicity(GF(13))).toBe(2); // 12 = 2^2 * 3
    expect(twoAdicity(GF(17))).toBe(4); // 16 = 2^4
    expect(twoAdicity(GF(97))).toBe(5); // 96 = 2^5 * 3
    expect(twoAdicity(GF(65537))).toBe(16); // 65536 = 2^16
  });
});

describe('createFFTDomain', () => {
  it('should create domain with exact power-of-2 size', () => {
    const F17 = GF(17);
    const domain = createFFTDomain(F17, 8n);

    expect(domain.size).toBe(8n);
  });

  it('should round down to nearest power of 2', () => {
    const F17 = GF(17);

    const domain = createFFTDomain(F17, 10n);
    expect(domain.size).toBe(8n);

    const domain2 = createFFTDomain(F17, 15n);
    expect(domain2.size).toBe(8n);
  });

  it('should cap at max supported size', () => {
    const F7 = GF(7); // max FFT size is 2
    const domain = createFFTDomain(F7, 16n);

    expect(domain.size).toBe(2n);
  });

  it('should throw in exact mode for non-power-of-2', () => {
    const F17 = GF(17);

    expect(() => createFFTDomain(F17, 10n, true)).toThrow();
  });

  it('should throw in exact mode for unsupported size', () => {
    const F7 = GF(7);

    expect(() => createFFTDomain(F7, 8n, true)).toThrow();
  });
});

describe('elementOfOrder', () => {
  it('should find elements of specified order', () => {
    const F13 = GF(13);

    // Orders dividing 12: 1, 2, 3, 4, 6, 12
    for (const order of [1n, 2n, 3n, 4n, 6n, 12n]) {
      const elem = elementOfOrder(F13, order);
      expect(multiplicative_order(elem)).toBe(order);
    }
  });

  it('should throw for orders that do not divide |F*|', () => {
    const F13 = GF(13);

    expect(() => elementOfOrder(F13, 5n)).toThrow();
    expect(() => elementOfOrder(F13, 7n)).toThrow();
  });
});

describe('ZK-specific tests', () => {
  describe('FRI-style domain folding', () => {
    it('should support iterative folding', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 16n);
      const offset = F17.__call__(3n);
      let coset: CosetDomain = domain.coset(offset);

      // Fold 4 times: 16 -> 8 -> 4 -> 2 -> 1
      for (let i = 0; i < 4; i++) {
        const challenge = F17.__call__(BigInt(i + 2));
        coset = coset.fold(challenge);
        expect(coset.size).toBe(BigInt(16 >> (i + 1)));
      }

      expect(coset.size).toBe(1n);
    });
  });

  describe('STARK-style domain properties', () => {
    it('should have disjoint cosets', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 4n);

      // Find an offset not in domain
      const offset = F17.__call__(2n); // 2 is not a 4th root of unity

      const coset = domain.coset(offset);

      // Check no overlap
      const domainElems = new Set(domain.elements().map((e) => e.toString()));
      for (const elem of coset.elements()) {
        expect(domainElems.has(elem.toString())).toBe(false);
      }
    });

    it('should compute vanishing polynomial correctly', () => {
      const F17 = GF(17);
      const domain = new FFTDomain(F17, 4n);
      const polyRing = new PolynomialRing(F17, 'x');

      const Z = domain.vanishingPolynomial(polyRing);

      // Z(x) = x^4 - 1
      // All domain elements should be roots
      for (const elem of domain.elements()) {
        const value = Z.evaluate(elem);
        expect(value.isZero()).toBe(true);
      }

      // A non-domain element should not be a root
      const nonMember = F17.__call__(2n);
      const value = Z.evaluate(nonMember);
      expect(value.isZero()).toBe(false);
    });
  });

  describe('BLS12-381 style field', () => {
    it('should work with realistic ZK field parameters', () => {
      // BLS12-381 scalar field has r-1 = 2^32 * t for some odd t
      // We'll use a smaller example: p such that p-1 = 2^20 * odd
      // p = 2^20 * 7 + 1 = 7340033 (this is prime)
      const p = 7340033n;
      const F = GF(p);

      // Should support FFT up to size 2^20
      expect(twoAdicity(F)).toBe(20);
      expect(maxFFTSize(F)).toBe(1048576n);

      // Create a moderately sized domain
      const domain = new FFTDomain(F, 1024n);
      expect(domain.size).toBe(1024n);

      // Verify generator properties
      expect(domain.generator.pow(1024n).isOne()).toBe(true);
      expect(domain.generator.pow(512n).isOne()).toBe(false);
    });
  });
});
