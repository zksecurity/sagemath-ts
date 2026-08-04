/**
 * Tests for prime finite fields GF(p)
 */

import { describe, expect, it } from 'bun:test';
import {
  FiniteFieldElement,
  FiniteFieldPrime,
  GF,
  analyzeFiniteFieldOrder,
} from './finite_field_constructor.js';
import { IntegerMod, Mod } from './integer_mod.js';
import { IntegerModRing, Zmod } from './integer_mod_ring.js';

describe('IntegerMod', () => {
  describe('basic operations', () => {
    it('should create elements correctly', () => {
      const a = Mod(7n, 10n);
      expect(a.value).toBe(7n);
      expect(a.modulus).toBe(10n);
    });

    it('should reduce values modulo n', () => {
      expect(Mod(15n, 10n).value).toBe(5n);
      expect(Mod(-3n, 10n).value).toBe(7n);
      expect(Mod(-13n, 10n).value).toBe(7n);
    });

    it('should add elements', () => {
      const a = Mod(7n, 10n);
      const b = Mod(8n, 10n);
      expect(a.add(b).value).toBe(5n); // 7 + 8 = 15 ≡ 5
    });

    it('should subtract elements', () => {
      const a = Mod(3n, 10n);
      const b = Mod(8n, 10n);
      expect(a.sub(b).value).toBe(5n); // 3 - 8 = -5 ≡ 5
    });

    it('should multiply elements', () => {
      const a = Mod(7n, 10n);
      const b = Mod(3n, 10n);
      expect(a.mul(b).value).toBe(1n); // 7 * 3 = 21 ≡ 1
    });

    it('should negate elements', () => {
      expect(Mod(3n, 10n).neg().value).toBe(7n);
      expect(Mod(0n, 10n).neg().value).toBe(0n);
    });

    it('should compute inverses when they exist', () => {
      const a = Mod(3n, 7n);
      const inv = a.inv();
      expect(a.mul(inv).value).toBe(1n);
    });

    it('should throw when computing inverse that does not exist', () => {
      const a = Mod(2n, 10n); // gcd(2, 10) = 2 ≠ 1
      expect(() => a.inv()).toThrow();
    });

    it('should compute powers', () => {
      const a = Mod(2n, 7n);
      expect(a.pow(3).value).toBe(1n); // 2^3 = 8 ≡ 1
      expect(a.pow(6).value).toBe(1n); // Fermat's little theorem
    });

    it('should check equality', () => {
      const a = Mod(3n, 7n);
      expect(a.eq(3n)).toBe(true);
      expect(a.eq(10n)).toBe(true); // 10 ≡ 3 mod 7
      expect(a.eq(4n)).toBe(false);
    });

    it('should identify zero and one', () => {
      expect(Mod(0n, 7n).isZero()).toBe(true);
      expect(Mod(7n, 7n).isZero()).toBe(true);
      expect(Mod(1n, 7n).isOne()).toBe(true);
      expect(Mod(8n, 7n).isOne()).toBe(true);
    });

    it('should identify units', () => {
      // In Z/10Z, units are 1, 3, 7, 9
      expect(Mod(1n, 10n).isUnit()).toBe(true);
      expect(Mod(3n, 10n).isUnit()).toBe(true);
      expect(Mod(7n, 10n).isUnit()).toBe(true);
      expect(Mod(9n, 10n).isUnit()).toBe(true);
      expect(Mod(2n, 10n).isUnit()).toBe(false);
      expect(Mod(5n, 10n).isUnit()).toBe(false);
    });
  });
});

describe('IntegerModRing (Zmod)', () => {
  describe('ring construction', () => {
    it('should create Z/nZ', () => {
      const Z10 = Zmod(10n);
      expect(Z10.modulus).toBe(10n);
      expect(Z10.characteristic).toBe(10n);
      expect(Z10.order).toBe(10n);
    });

    it('should reject non-positive modulus', () => {
      expect(() => Zmod(0n)).toThrow();
      expect(() => Zmod(-5n)).toThrow();
    });

    it('should create elements via __call__', () => {
      const Z7 = Zmod(7n);
      const a = Z7.__call__(10n);
      expect(a.value).toBe(3n);
    });

    it('should return correct zero and one', () => {
      const Z7 = Zmod(7n);
      expect(Z7.zero().value).toBe(0n);
      expect(Z7.one().value).toBe(1n);
    });

    it('should detect if ring is a field', () => {
      expect(Zmod(7n).is_field()).toBe(true); // 7 is prime
      expect(Zmod(10n).is_field()).toBe(false); // 10 is not prime
      expect(Zmod(4n).is_field()).toBe(false); // 4 is not prime
    });
  });

  describe('iteration', () => {
    it('should iterate over all elements', () => {
      const Z5 = Zmod(5n);
      const elements = [...Z5];
      expect(elements.length).toBe(5);
      expect(elements.map((e) => e.value)).toEqual([0n, 1n, 2n, 3n, 4n]);
    });

    it('should list all elements', () => {
      const Z3 = Zmod(3n);
      const elements = Z3.list();
      expect(elements.map((e) => e.value)).toEqual([0n, 1n, 2n]);
    });

    it('should list all units', () => {
      const Z10 = Zmod(10n);
      const units = Z10.units();
      expect(units.map((e) => e.value).sort((a, b) => Number(a - b))).toEqual([1n, 3n, 7n, 9n]);
    });
  });
});

describe('GF (FiniteFieldPrime)', () => {
  describe('field construction', () => {
    it('should create GF(p) for prime p', () => {
      const F7 = GF(7n);
      expect(F7.characteristic).toBe(7n);
      expect(F7.order).toBe(7n);
      expect(F7.degree).toBe(1);
    });

    it('should reject non-prime order', () => {
      expect(() => GF(10n)).toThrow();
      expect(() => GF(1n)).toThrow();
    });

    it('should reject non-prime powers correctly', () => {
      expect(() => GF(6n)).toThrow(); // 6 = 2 * 3
      expect(() => GF(12n)).toThrow(); // 12 = 2^2 * 3
    });

    it('should accept large primes', () => {
      // Mersenne prime 2^31 - 1
      const p = 2147483647n;
      const F = GF(p);
      expect(F.characteristic).toBe(p);
    });

    it('should allow skipping primality check', () => {
      // This would fail the check but we skip it
      const F = GF(7n, { check: false });
      expect(F.characteristic).toBe(7n);
    });
  });

  describe('field arithmetic in GF(7)', () => {
    const F7 = GF(7n);

    it('should add elements', () => {
      const a = F7.__call__(3n);
      const b = F7.__call__(5n);
      expect(a.add(b).value).toBe(1n); // 3 + 5 = 8 ≡ 1
    });

    it('should subtract elements', () => {
      const a = F7.__call__(2n);
      const b = F7.__call__(5n);
      expect(a.sub(b).value).toBe(4n); // 2 - 5 = -3 ≡ 4
    });

    it('should multiply elements', () => {
      const a = F7.__call__(3n);
      const b = F7.__call__(5n);
      expect(a.mul(b).value).toBe(1n); // 3 * 5 = 15 ≡ 1
    });

    it('should divide elements', () => {
      const a = F7.__call__(6n);
      const b = F7.__call__(2n);
      expect(a.div(b).value).toBe(3n); // 6 / 2 = 3
    });

    it('should negate elements', () => {
      const a = F7.__call__(3n);
      expect(a.neg().value).toBe(4n); // -3 ≡ 4
    });

    it('should compute inverses', () => {
      // Check all non-zero elements have inverses
      for (let i = 1n; i < 7n; i++) {
        const elem = F7.__call__(i);
        const inv = elem.inv();
        expect(elem.mul(inv).isOne()).toBe(true);
      }
    });

    it('should throw on division by zero', () => {
      const a = F7.__call__(3n);
      expect(() => a.div(0n)).toThrow();
      expect(() => F7.zero().inv()).toThrow();
    });
  });

  describe('field arithmetic in GF(13)', () => {
    const F13 = GF(13n);

    it('should compute all inverses', () => {
      for (let i = 1n; i < 13n; i++) {
        const elem = F13.__call__(i);
        const inv = elem.inv();
        expect(elem.mul(inv).isOne()).toBe(true);
      }
    });

    it('should satisfy Fermat little theorem', () => {
      for (let i = 1n; i < 13n; i++) {
        const elem = F13.__call__(i);
        // a^(p-1) = 1 for all a != 0
        expect(elem.pow(12n).isOne()).toBe(true);
      }
    });

    it('should compute powers correctly', () => {
      const a = F13.__call__(2n);
      expect(a.pow(0).value).toBe(1n);
      expect(a.pow(1).value).toBe(2n);
      expect(a.pow(2).value).toBe(4n);
      expect(a.pow(3).value).toBe(8n);
      expect(a.pow(4).value).toBe(3n); // 16 ≡ 3 mod 13
    });

    it('should compute negative powers via inverse', () => {
      const a = F13.__call__(2n);
      const invA = a.inv();
      expect(a.pow(-1n).value).toBe(invA.value);
      expect(a.pow(-2n).value).toBe(invA.pow(2n).value);
    });
  });

  describe('field arithmetic in GF(101)', () => {
    const F101 = GF(101n);

    it('should handle arithmetic near modulus', () => {
      const a = F101.__call__(100n);
      const b = F101.__call__(2n);
      expect(a.add(b).value).toBe(1n); // 100 + 2 = 102 ≡ 1
    });

    it('should compute all inverses', () => {
      // Test a sampling of inverses
      const testValues = [1n, 2n, 10n, 50n, 99n, 100n];
      for (const i of testValues) {
        const elem = F101.__call__(i);
        const inv = elem.inv();
        expect(elem.mul(inv).isOne()).toBe(true);
      }
    });

    it('should find multiplicative generator', () => {
      const g = F101.multiplicative_generator();
      // g should generate the entire multiplicative group
      // The multiplicative group has order 100
      expect(g.multiplicative_order()).toBe(100n);
    });
  });

  describe('field axioms', () => {
    const F7 = GF(7n);

    it('should satisfy additive commutativity', () => {
      for (let i = 0n; i < 7n; i++) {
        for (let j = 0n; j < 7n; j++) {
          const a = F7.__call__(i);
          const b = F7.__call__(j);
          expect(a.add(b).value).toBe(b.add(a).value);
        }
      }
    });

    it('should satisfy multiplicative commutativity', () => {
      for (let i = 0n; i < 7n; i++) {
        for (let j = 0n; j < 7n; j++) {
          const a = F7.__call__(i);
          const b = F7.__call__(j);
          expect(a.mul(b).value).toBe(b.mul(a).value);
        }
      }
    });

    it('should satisfy additive associativity', () => {
      const a = F7.__call__(2n);
      const b = F7.__call__(3n);
      const c = F7.__call__(5n);
      expect(a.add(b).add(c).value).toBe(a.add(b.add(c)).value);
    });

    it('should satisfy multiplicative associativity', () => {
      const a = F7.__call__(2n);
      const b = F7.__call__(3n);
      const c = F7.__call__(5n);
      expect(a.mul(b).mul(c).value).toBe(a.mul(b.mul(c)).value);
    });

    it('should satisfy distributivity', () => {
      const a = F7.__call__(2n);
      const b = F7.__call__(3n);
      const c = F7.__call__(5n);
      // a * (b + c) = a*b + a*c
      expect(a.mul(b.add(c)).value).toBe(a.mul(b).add(a.mul(c)).value);
    });

    it('should have additive identity', () => {
      for (let i = 0n; i < 7n; i++) {
        const a = F7.__call__(i);
        expect(a.add(F7.zero()).value).toBe(a.value);
      }
    });

    it('should have multiplicative identity', () => {
      for (let i = 0n; i < 7n; i++) {
        const a = F7.__call__(i);
        expect(a.mul(F7.one()).value).toBe(a.value);
      }
    });

    it('should have additive inverses', () => {
      for (let i = 0n; i < 7n; i++) {
        const a = F7.__call__(i);
        expect(a.add(a.neg()).isZero()).toBe(true);
      }
    });

    it('should have multiplicative inverses for non-zero elements', () => {
      for (let i = 1n; i < 7n; i++) {
        const a = F7.__call__(i);
        expect(a.mul(a.inv()).isOne()).toBe(true);
      }
    });
  });

  describe('quadratic residues', () => {
    const F7 = GF(7n);

    it('should identify quadratic residues', () => {
      // In GF(7), squares are: 0^2=0, 1^2=1, 2^2=4, 3^2=2, 4^2=2, 5^2=4, 6^2=1
      // So residues are 0, 1, 2, 4
      expect(F7.__call__(0n).is_square()).toBe(true);
      expect(F7.__call__(1n).is_square()).toBe(true);
      expect(F7.__call__(2n).is_square()).toBe(true);
      expect(F7.__call__(3n).is_square()).toBe(false);
      expect(F7.__call__(4n).is_square()).toBe(true);
      expect(F7.__call__(5n).is_square()).toBe(false);
      expect(F7.__call__(6n).is_square()).toBe(false);
    });

    it('should compute square roots', () => {
      const a = F7.__call__(2n);
      const sqrtA = a.sqrt();
      expect(sqrtA.mul(sqrtA).value).toBe(2n);
    });

    it('should throw for non-residues', () => {
      const a = F7.__call__(3n);
      expect(() => a.sqrt()).toThrow();
    });

    it('should find quadratic non-residue', () => {
      const nr = F7.quadratic_non_residue();
      expect(nr.is_square()).toBe(false);
    });
  });

  describe('iteration', () => {
    it('should iterate over all elements', () => {
      const F5 = GF(5n);
      const elements = [...F5];
      expect(elements.length).toBe(5);
      expect(elements.map((e) => e.value)).toEqual([0n, 1n, 2n, 3n, 4n]);
    });

    it('should list all elements', () => {
      const F3 = GF(3n);
      const elements = F3.list();
      expect(elements.map((e) => e.value)).toEqual([0n, 1n, 2n]);
    });
  });

  describe('large prime fields', () => {
    // Curve25519 prime: 2^255 - 19
    const p25519 = 2n ** 255n - 19n;

    it('should create field with Curve25519 prime', () => {
      const Fp = GF(p25519);
      expect(Fp.characteristic).toBe(p25519);
    });

    it('should perform arithmetic with large values', () => {
      const Fp = GF(p25519);
      const a = Fp.__call__(12345678901234567890123456789n);
      const b = Fp.__call__(98765432109876543210987654321n);

      // Test basic operations
      const sum = a.add(b);
      const diff = a.sub(b);
      const prod = a.mul(b);
      const quot = a.div(b);

      // Verify inverse property
      expect(quot.mul(b).value).toBe(a.value);
    });

    it('should compute inverses for large values', () => {
      const Fp = GF(p25519);
      const a = Fp.__call__(2n ** 200n + 1n);
      const inv = a.inv();
      expect(a.mul(inv).isOne()).toBe(true);
    });

    it('should compute large powers efficiently', () => {
      const Fp = GF(p25519);
      const a = Fp.__call__(7n);

      // Fermat's little theorem: a^(p-1) = 1
      expect(a.pow(p25519 - 1n).isOne()).toBe(true);
    });

    it('should handle values larger than prime', () => {
      const Fp = GF(p25519);
      const large = p25519 + 100n;
      const elem = Fp.__call__(large);
      expect(elem.value).toBe(100n);
    });
  });

  describe('special primes', () => {
    it('should work with GF(2)', () => {
      const F2 = GF(2n);
      expect(F2.cardinality()).toBe(2n);

      const zero = F2.zero();
      const one = F2.one();

      expect(zero.add(zero).value).toBe(0n);
      expect(zero.add(one).value).toBe(1n);
      expect(one.add(one).value).toBe(0n);
      expect(one.mul(one).value).toBe(1n);
    });

    it('should work with GF(3)', () => {
      const F3 = GF(3n);
      const one = F3.one();
      const two = F3.__call__(2n);

      expect(two.add(two).value).toBe(1n); // 2+2=4≡1
      expect(two.mul(two).value).toBe(1n); // 2*2=4≡1
      expect(two.inv().value).toBe(2n); // 2*2=4≡1, so 2^(-1)=2
    });

    it('should work with Mersenne prime 2^31-1', () => {
      const p = 2n ** 31n - 1n;
      const Fp = GF(p);

      const a = Fp.__call__(1234567890n);
      const b = Fp.__call__(987654321n);

      const prod = a.mul(b);
      const quot = a.div(b);

      expect(quot.mul(b).value).toBe(a.value);
    });

    it('should work with 64-bit prime', () => {
      // A large 64-bit prime
      const p = 18446744073709551557n; // 2^64 - 59
      const Fp = GF(p);

      const a = Fp.__call__(2n ** 63n);
      const b = a.inv();
      expect(a.mul(b).isOne()).toBe(true);
    });
  });
});

describe('analyzeFiniteFieldOrder', () => {
  it('should analyze prime orders', () => {
    const info = analyzeFiniteFieldOrder(7n);
    expect(info).not.toBeNull();
    expect(info!.order).toBe(7n);
    expect(info!.characteristic).toBe(7n);
    expect(info!.degree).toBe(1n);
    expect(info!.isPrimeField).toBe(true);
  });

  it('should analyze prime power orders', () => {
    const info = analyzeFiniteFieldOrder(8n); // 2^3
    expect(info).not.toBeNull();
    expect(info!.order).toBe(8n);
    expect(info!.characteristic).toBe(2n);
    expect(info!.degree).toBe(3n);
    expect(info!.isPrimeField).toBe(false);
  });

  it('should reject non-prime powers', () => {
    expect(analyzeFiniteFieldOrder(6n)).toBeNull(); // 2 * 3
    expect(analyzeFiniteFieldOrder(12n)).toBeNull(); // 2^2 * 3
    expect(analyzeFiniteFieldOrder(1n)).toBeNull(); // too small
  });
});

describe('multiplicative order', () => {
  it('should compute multiplicative order in GF(7)', () => {
    const F7 = GF(7n);

    // 1 has order 1
    expect(F7.__call__(1n).multiplicative_order()).toBe(1n);

    // 6 = -1 has order 2 (since 6^2 = 36 ≡ 1)
    expect(F7.__call__(6n).multiplicative_order()).toBe(2n);

    // 2 has order 3 (2^3 = 8 ≡ 1)
    expect(F7.__call__(2n).multiplicative_order()).toBe(3n);

    // 3 is a generator (order 6 = p-1)
    expect(F7.__call__(3n).multiplicative_order()).toBe(6n);
  });

  it('should compute multiplicative order in GF(13)', () => {
    const F13 = GF(13n);

    // All non-zero elements should have orders dividing 12
    for (let i = 1n; i < 13n; i++) {
      const elem = F13.__call__(i);
      const order = elem.multiplicative_order();
      expect(12n % order).toBe(0n);
      expect(elem.pow(order).isOne()).toBe(true);
    }
  });
});
