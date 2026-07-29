/**
 * Tests for Finite Field Extensions GF(p^n)
 *
 * Verifies the implementation against field axioms and SageMath behavior.
 */
import { describe, expect, test } from 'bun:test';
import {
  available_characteristics,
  available_degrees,
  conway_polynomial,
  has_conway_polynomial,
} from './conway_polynomials.js';
import {
  FiniteFieldElement,
  type FiniteFieldExtension,
  GF,
  GFpn,
  PrimeField,
  PrimeFieldElement,
} from './finite_field_extension.js';

describe('Conway Polynomials', () => {
  test('has Conway polynomials for common cases', () => {
    expect(has_conway_polynomial(2, 2)).toBe(true);
    expect(has_conway_polynomial(2, 8)).toBe(true);
    expect(has_conway_polynomial(3, 2)).toBe(true);
    expect(has_conway_polynomial(5, 2)).toBe(true);
    expect(has_conway_polynomial(7, 2)).toBe(true);
  });

  test('Conway polynomial for GF(4)', () => {
    // x^2 + x + 1
    const coeffs = conway_polynomial(2, 2);
    expect(coeffs).toEqual([1, 1]); // 1 + x (monic, so x^2 implicit)
  });

  test('Conway polynomial for GF(8)', () => {
    // x^3 + x + 1
    const coeffs = conway_polynomial(2, 3);
    expect(coeffs).toEqual([1, 1, 0]); // 1 + x
  });

  test('Conway polynomial for GF(9)', () => {
    // x^2 + 2x + 2
    const coeffs = conway_polynomial(3, 2);
    expect(coeffs).toEqual([2, 2]);
  });

  test('available characteristics', () => {
    const chars = available_characteristics();
    expect(chars).toContain(2);
    expect(chars).toContain(3);
    expect(chars).toContain(5);
    expect(chars).toContain(7);
  });

  test('available degrees for characteristic 2', () => {
    const degrees = available_degrees(2);
    expect(degrees).toContain(2);
    expect(degrees).toContain(3);
    expect(degrees).toContain(4);
    expect(degrees).toContain(8);
    expect(degrees).toContain(16);
  });
});

describe('Prime Fields GF(p)', () => {
  test('GF(2) construction', () => {
    const F2 = GF(2) as PrimeField;
    expect(F2.characteristic).toBe(2n);
    expect(F2.order).toBe(2n);
    expect(F2.cardinality()).toBe(2n);
  });

  test('GF(3) construction', () => {
    const F3 = GF(3) as PrimeField;
    expect(F3.characteristic).toBe(3n);
    expect(F3.order).toBe(3n);
  });

  test('GF(5) arithmetic', () => {
    const F5 = GF(5) as PrimeField;

    const a = F5.__call__(3);
    const b = F5.__call__(4);

    expect(a.add(b).value).toBe(2n); // 3 + 4 = 7 = 2 mod 5
    expect(a.mul(b).value).toBe(2n); // 3 * 4 = 12 = 2 mod 5
    expect(a.sub(b).value).toBe(4n); // 3 - 4 = -1 = 4 mod 5
    expect(a.inv().value).toBe(2n); // 3 * 2 = 6 = 1 mod 5
  });

  test('GF(7) iteration', () => {
    const F7 = GF(7) as PrimeField;
    const elements = [...F7];
    expect(elements.length).toBe(7);
  });
});

describe('GF(4) = GF(2^2)', () => {
  const F4 = GF(4) as FiniteFieldExtension;

  test('field construction', () => {
    expect(F4.characteristic).toBe(2n);
    expect(F4.degree).toBe(2);
    expect(F4.order).toBe(4n);
    expect(F4.cardinality()).toBe(4n);
  });

  test('has 4 elements', () => {
    const elements = [...F4];
    expect(elements.length).toBe(4);
  });

  test('generator satisfies modulus', () => {
    const a = F4.gen();

    // a^2 + a + 1 = 0 (Conway polynomial for GF(4))
    const result = a.pow(2).add(a).add(F4.one());
    expect(result.isZero()).toBe(true);
  });

  test('arithmetic', () => {
    const zero = F4.zero();
    const one = F4.one();
    const a = F4.gen();

    // Addition
    expect(one.add(one).isZero()).toBe(true); // 1 + 1 = 0 in char 2
    expect(a.add(a).isZero()).toBe(true); // a + a = 0 in char 2

    // Multiplication
    expect(one.mul(a).eq(a)).toBe(true);
    expect(a.mul(a).eq(a.add(one))).toBe(true); // a^2 = a + 1

    // Subtraction (same as addition in char 2)
    expect(a.sub(one).eq(a.add(one))).toBe(true);
  });

  test('multiplicative inverses', () => {
    let nonZeroCount = 0;
    for (const elem of F4) {
      if (!elem.isZero()) {
        const inv = elem.inv();
        expect(elem.mul(inv).isOne()).toBe(true);
        nonZeroCount++;
      }
    }
    expect(nonZeroCount).toBe(3); // 4 - 1 non-zero elements
  });

  test('multiplicative group is cyclic of order 3', () => {
    const a = F4.gen();

    // a^3 = 1
    expect(a.pow(3).isOne()).toBe(true);

    // Elements: 1, a, a^2 = a+1
    const one = F4.one();
    const aSq = a.pow(2);

    // Check these are all distinct
    expect(one.eq(a)).toBe(false);
    expect(one.eq(aSq)).toBe(false);
    expect(a.eq(aSq)).toBe(false);
  });

  test('Frobenius automorphism', () => {
    const a = F4.gen();

    // Frobenius: x -> x^2
    const aFrob = a.frobenius();
    expect(aFrob.eq(a.pow(2))).toBe(true);

    // a^2 = a + 1 in GF(4)
    expect(aFrob.eq(a.add(F4.one()))).toBe(true);

    // Frobenius squared is identity on GF(4)
    expect(a.frobenius(2).eq(a)).toBe(true);
  });
});

describe('GF(8) = GF(2^3)', () => {
  const F8 = GF(8) as FiniteFieldExtension;

  test('field construction', () => {
    expect(F8.characteristic).toBe(2n);
    expect(F8.degree).toBe(3);
    expect(F8.order).toBe(8n);
  });

  test('has 8 elements', () => {
    const elements = [...F8];
    expect(elements.length).toBe(8);
  });

  test('generator satisfies modulus', () => {
    const a = F8.gen();

    // a^3 + a + 1 = 0 (Conway polynomial for GF(8))
    const result = a.pow(3).add(a).add(F8.one());
    expect(result.isZero()).toBe(true);
  });

  test('multiplicative inverses', () => {
    let nonZeroCount = 0;
    for (const elem of F8) {
      if (!elem.isZero()) {
        const inv = elem.inv();
        expect(elem.mul(inv).isOne()).toBe(true);
        nonZeroCount++;
      }
    }
    expect(nonZeroCount).toBe(7);
  });

  test('multiplicative group has order 7 (prime)', () => {
    const a = F8.gen();

    // a^7 = 1
    expect(a.pow(7).isOne()).toBe(true);

    // a is a primitive element (generates all non-zero elements)
    const seen = new Set<string>();
    let elem = F8.one();
    for (let i = 0; i < 7; i++) {
      seen.add(elem.toString());
      elem = elem.mul(a);
    }
    expect(seen.size).toBe(7);
  });

  test('Frobenius automorphism', () => {
    const a = F8.gen();

    // Frobenius: x -> x^2
    const aFrob = a.frobenius();
    expect(aFrob.eq(a.pow(2))).toBe(true);

    // Frobenius^3 = identity on GF(8)
    expect(a.frobenius(3).eq(a)).toBe(true);
  });
});

describe('GF(16) = GF(2^4)', () => {
  const F16 = GF(16) as FiniteFieldExtension;

  test('field construction', () => {
    expect(F16.characteristic).toBe(2n);
    expect(F16.degree).toBe(4);
    expect(F16.order).toBe(16n);
  });

  test('has 16 elements', () => {
    const elements = [...F16];
    expect(elements.length).toBe(16);
  });

  test('generator satisfies modulus', () => {
    const a = F16.gen();

    // a^4 + a + 1 = 0 (Conway polynomial for GF(16))
    const result = a.pow(4).add(a).add(F16.one());
    expect(result.isZero()).toBe(true);
  });

  test('multiplicative inverses', () => {
    for (const elem of F16) {
      if (!elem.isZero()) {
        const inv = elem.inv();
        expect(elem.mul(inv).isOne()).toBe(true);
      }
    }
  });

  test('multiplicative group has order 15', () => {
    const a = F16.gen();
    expect(a.pow(15).isOne()).toBe(true);
  });
});

describe('GF(9) = GF(3^2)', () => {
  const F9 = GF(9) as FiniteFieldExtension;

  test('field construction', () => {
    expect(F9.characteristic).toBe(3n);
    expect(F9.degree).toBe(2);
    expect(F9.order).toBe(9n);
  });

  test('has 9 elements', () => {
    const elements = [...F9];
    expect(elements.length).toBe(9);
  });

  test('generator satisfies modulus', () => {
    const a = F9.gen();

    // a^2 + 2a + 2 = 0 (Conway polynomial for GF(9))
    const result = a.pow(2).add(F9.__call__(2).mul(a)).add(F9.__call__(2));
    expect(result.isZero()).toBe(true);
  });

  test('characteristic 3 property', () => {
    // In GF(9), 3 * x = 0 for all x
    for (const elem of F9) {
      expect(elem.add(elem).add(elem).isZero()).toBe(true);
    }
  });

  test('multiplicative inverses', () => {
    for (const elem of F9) {
      if (!elem.isZero()) {
        const inv = elem.inv();
        expect(elem.mul(inv).isOne()).toBe(true);
      }
    }
  });

  test('multiplicative group has order 8', () => {
    const a = F9.gen();
    expect(a.pow(8).isOne()).toBe(true);
  });

  test('Frobenius automorphism', () => {
    const a = F9.gen();

    // Frobenius: x -> x^3
    const aFrob = a.frobenius();
    expect(aFrob.eq(a.pow(3))).toBe(true);

    // Frobenius^2 = identity on GF(9)
    expect(a.frobenius(2).eq(a)).toBe(true);
  });
});

describe('GF(27) = GF(3^3)', () => {
  const F27 = GF(27) as FiniteFieldExtension;

  test('field construction', () => {
    expect(F27.characteristic).toBe(3n);
    expect(F27.degree).toBe(3);
    expect(F27.order).toBe(27n);
  });

  test('has 27 elements', () => {
    const elements = [...F27];
    expect(elements.length).toBe(27);
  });

  test('generator satisfies modulus', () => {
    const a = F27.gen();

    // a^3 + 2a + 1 = 0 (Conway polynomial for GF(27))
    const result = a.pow(3).add(F27.__call__(2).mul(a)).add(F27.one());
    expect(result.isZero()).toBe(true);
  });

  test('multiplicative inverses', () => {
    for (const elem of F27) {
      if (!elem.isZero()) {
        const inv = elem.inv();
        expect(elem.mul(inv).isOne()).toBe(true);
      }
    }
  });

  test('multiplicative group has order 26', () => {
    const a = F27.gen();
    expect(a.pow(26).isOne()).toBe(true);
  });

  test('Frobenius automorphism', () => {
    const a = F27.gen();

    // Frobenius^3 = identity on GF(27)
    expect(a.frobenius(3).eq(a)).toBe(true);
  });
});

describe('GF(25) = GF(5^2)', () => {
  const F25 = GF(25) as FiniteFieldExtension;

  test('field construction', () => {
    expect(F25.characteristic).toBe(5n);
    expect(F25.degree).toBe(2);
    expect(F25.order).toBe(25n);
  });

  test('has 25 elements', () => {
    const elements = [...F25];
    expect(elements.length).toBe(25);
  });

  test('generator satisfies modulus', () => {
    const a = F25.gen();

    // a^2 + 4a + 2 = 0 (Conway polynomial for GF(25))
    const result = a.pow(2).add(F25.__call__(4).mul(a)).add(F25.__call__(2));
    expect(result.isZero()).toBe(true);
  });

  test('multiplicative inverses', () => {
    for (const elem of F25) {
      if (!elem.isZero()) {
        const inv = elem.inv();
        expect(elem.mul(inv).isOne()).toBe(true);
      }
    }
  });

  test('multiplicative group has order 24', () => {
    const a = F25.gen();
    expect(a.pow(24).isOne()).toBe(true);
  });
});

describe('Field Axioms', () => {
  const fields = [
    { name: 'GF(4)', F: GF(4) as FiniteFieldExtension },
    { name: 'GF(8)', F: GF(8) as FiniteFieldExtension },
    { name: 'GF(9)', F: GF(9) as FiniteFieldExtension },
  ];

  for (const { name, F } of fields) {
    describe(name, () => {
      test('additive identity', () => {
        const zero = F.zero();
        for (const a of F) {
          expect(a.add(zero).eq(a)).toBe(true);
        }
      });

      test('additive inverse', () => {
        for (const a of F) {
          const negA = a.neg();
          expect(a.add(negA).isZero()).toBe(true);
        }
      });

      test('multiplicative identity', () => {
        const one = F.one();
        for (const a of F) {
          expect(a.mul(one).eq(a)).toBe(true);
        }
      });

      test('multiplicative inverse (non-zero)', () => {
        for (const a of F) {
          if (!a.isZero()) {
            const invA = a.inv();
            expect(a.mul(invA).isOne()).toBe(true);
          }
        }
      });

      test('commutativity of addition', () => {
        const elements = [...F];
        for (let i = 0; i < Math.min(elements.length, 10); i++) {
          for (let j = 0; j < Math.min(elements.length, 10); j++) {
            const a = elements[i]!;
            const b = elements[j]!;
            expect(a.add(b).eq(b.add(a))).toBe(true);
          }
        }
      });

      test('commutativity of multiplication', () => {
        const elements = [...F];
        for (let i = 0; i < Math.min(elements.length, 10); i++) {
          for (let j = 0; j < Math.min(elements.length, 10); j++) {
            const a = elements[i]!;
            const b = elements[j]!;
            expect(a.mul(b).eq(b.mul(a))).toBe(true);
          }
        }
      });

      test('associativity of addition', () => {
        const elements = [...F].slice(0, 5);
        for (const a of elements) {
          for (const b of elements) {
            for (const c of elements) {
              expect(a.add(b.add(c)).eq(a.add(b).add(c))).toBe(true);
            }
          }
        }
      });

      test('associativity of multiplication', () => {
        const elements = [...F].slice(0, 5);
        for (const a of elements) {
          for (const b of elements) {
            for (const c of elements) {
              expect(a.mul(b.mul(c)).eq(a.mul(b).mul(c))).toBe(true);
            }
          }
        }
      });

      test('distributivity', () => {
        const elements = [...F].slice(0, 5);
        for (const a of elements) {
          for (const b of elements) {
            for (const c of elements) {
              // a * (b + c) = a*b + a*c
              expect(a.mul(b.add(c)).eq(a.mul(b).add(a.mul(c)))).toBe(true);
            }
          }
        }
      });
    });
  }
});

describe('Frobenius Automorphism', () => {
  test('GF(4): Frobenius is non-trivial automorphism', () => {
    const F4 = GF(4) as FiniteFieldExtension;
    const a = F4.gen();

    // a and a^2 are conjugates
    const aFrob = a.frobenius();
    expect(aFrob.eq(a)).toBe(false);
    expect(aFrob.eq(a.pow(2))).toBe(true);
  });

  test('GF(8): Frobenius has order 3', () => {
    const F8 = GF(8) as FiniteFieldExtension;
    const a = F8.gen();

    expect(a.frobenius(1).eq(a)).toBe(false);
    expect(a.frobenius(2).eq(a)).toBe(false);
    expect(a.frobenius(3).eq(a)).toBe(true);
  });

  test('GF(9): Frobenius fixes GF(3)', () => {
    const F9 = GF(9) as FiniteFieldExtension;

    // Elements of GF(3) embedded in GF(9) are fixed by Frobenius
    for (let i = 0; i < 3; i++) {
      const elem = F9.__call__(i);
      expect(elem.frobenius().eq(elem)).toBe(true);
    }
  });

  test('Frobenius is a field homomorphism', () => {
    const F8 = GF(8) as FiniteFieldExtension;
    const elements = [...F8].slice(0, 5);

    for (const a of elements) {
      for (const b of elements) {
        // Frobenius(a + b) = Frobenius(a) + Frobenius(b)
        expect(a.add(b).frobenius().eq(a.frobenius().add(b.frobenius()))).toBe(true);

        // Frobenius(a * b) = Frobenius(a) * Frobenius(b)
        expect(a.mul(b).frobenius().eq(a.frobenius().mul(b.frobenius()))).toBe(true);
      }
    }
  });
});

describe('Trace and Norm', () => {
  test('GF(4) trace', () => {
    const F4 = GF(4) as FiniteFieldExtension;

    // Trace: x + x^2
    // For elements of GF(2), trace is 0
    expect(F4.zero().trace().value).toBe(0n);
    expect(F4.one().trace().value).toBe(0n); // 1 + 1 = 0 in char 2

    // For generator a: Tr(a) = a + a^2 = a + (a+1) = 1 in char 2
    const a = F4.gen();
    expect(a.trace().value).toBe(1n);
  });

  test('GF(9) trace', () => {
    const F9 = GF(9) as FiniteFieldExtension;

    // Trace: x + x^3
    // For elements of GF(3), trace is 2x (since x + x = 2x in char 3)
    expect(F9.zero().trace().value).toBe(0n);
    expect(F9.one().trace().value).toBe(2n); // 1 + 1 = 2
    expect(F9.__call__(2).trace().value).toBe(1n); // 2 + 2 = 4 = 1 mod 3
  });
});

describe('Element Representation', () => {
  test('GF(4) element coefficients', () => {
    const F4 = GF(4) as FiniteFieldExtension;
    const a = F4.gen();

    // a = 0 + 1*a
    const coeffs = a.coefficients();
    expect(coeffs[0]!.value).toBe(0n);
    expect(coeffs[1]!.value).toBe(1n);
  });

  test('GF(8) element integer representation', () => {
    const F8 = GF(8) as FiniteFieldExtension;

    // Elements can be represented as integers 0-7
    const seen = new Set<bigint>();
    for (const elem of F8) {
      const intRep = elem.integer_representation();
      expect(intRep >= 0n && intRep < 8n).toBe(true);
      seen.add(intRep);
    }
    expect(seen.size).toBe(8);
  });

  test('fromInteger round-trip', () => {
    const F8 = GF(8) as FiniteFieldExtension;

    for (let i = 0n; i < 8n; i++) {
      const elem = F8.fromInteger(i);
      expect(elem.integer_representation()).toBe(i);
    }
  });
});

describe('GFpn explicit construction', () => {
  test('GFpn(2, 4) equals GF(16)', () => {
    const F16a = GFpn(2, 4);
    const F16b = GF(16) as FiniteFieldExtension;

    expect(F16a.order).toBe(F16b.order);
    expect(F16a.characteristic).toBe(F16b.characteristic);
    expect(F16a.degree).toBe(F16b.degree);
  });

  test('GFpn with custom modulus', () => {
    // x^2 + x + 1 over GF(2)
    const F4 = GFpn(2, 2, [1, 1]);
    expect(F4.order).toBe(4n);

    const a = F4.gen();
    // a^2 + a + 1 = 0
    expect(a.pow(2).add(a).add(F4.one()).isZero()).toBe(true);
  });

  test('custom variable name', () => {
    const F4 = GFpn(2, 2, undefined, 'alpha');
    expect(F4.toString()).toContain('alpha');
    expect(F4.gen().toString()).toBe('alpha');
  });
});

describe('Error handling', () => {
  test('GF rejects non-prime-powers', () => {
    expect(() => GF(6)).toThrow();
    expect(() => GF(10)).toThrow();
    expect(() => GF(15)).toThrow();
  });

  test('GF rejects values < 2', () => {
    expect(() => GF(1)).toThrow();
    expect(() => GF(0)).toThrow();
  });

  test('division by zero throws', () => {
    const F4 = GF(4) as FiniteFieldExtension;
    expect(() => F4.zero().inv()).toThrow();
    expect(() => F4.one().div(F4.zero())).toThrow();
  });

  test('PrimeField rejects non-primes', () => {
    expect(() => new PrimeField(4)).toThrow();
    expect(() => new PrimeField(1)).toThrow();
  });
});

describe('Large fields', () => {
  test('GF(256) = GF(2^8)', () => {
    const F256 = GF(256) as FiniteFieldExtension;
    expect(F256.order).toBe(256n);
    expect(F256.characteristic).toBe(2n);
    expect(F256.degree).toBe(8);

    // Don't enumerate all elements, just check structure
    const a = F256.gen();
    expect(a.pow(255).isOne()).toBe(true);
  });

  test('GF(49) = GF(7^2)', () => {
    const F49 = GF(49) as FiniteFieldExtension;
    expect(F49.order).toBe(49n);
    expect(F49.characteristic).toBe(7n);
    expect(F49.degree).toBe(2);

    const a = F49.gen();
    expect(a.pow(48).isOne()).toBe(true);
  });
});

describe('Default modulus is irreducible', () => {
  // The Conway table only covers p <= 31, so every prime beyond it exercises
  // findIrreducible(). No test above p = 7 existed, which is why a reducible
  // default modulus went unnoticed.
  //
  // A degree-2 or degree-3 polynomial over GF(p) is irreducible iff it has no
  // root in GF(p), so root-finding is a complete check at these degrees and is
  // independent of the implementation under test.
  const hasRootInBaseField = (F: FiniteFieldExtension): bigint | null => {
    for (let r = 0n; r < F.characteristic; r++) {
      if (F.modulus.evaluate(F.baseField.__call__(r)).value === 0n) {
        return r;
      }
    }
    return null;
  };

  // p = 787 = 1 (mod 3), so x^2 + x + 1 = (x - 379)(x - 407) over GF(787).
  // A polynomial that splits into distinct linear factors satisfies
  // x^p = x mod f, the exact case the irreducibility test must reject.
  test('GF(787^2) does not use the reducible x^2 + x + 1', () => {
    const F = GFpn(787, 2);
    expect(hasRootInBaseField(F)).toBe(null);
  });

  test('every default modulus past the Conway table is irreducible', () => {
    const primes = [37, 41, 43, 61, 67, 73, 79, 97, 101, 127, 151, 181, 199, 787, 1009];
    for (const p of primes) {
      for (const n of [2, 3]) {
        const F = GFpn(p, n);
        const root = hasRootInBaseField(F);
        expect(`GF(${p}^${n}) modulus=${F.modulus.toString()} root=${root}`).toBe(
          `GF(${p}^${n}) modulus=${F.modulus.toString()} root=null`
        );
      }
    }
  });

  // Degree-independent check: a reducible modulus makes the quotient a ring
  // with zero divisors, so some nonzero element fails to be invertible.
  test('nonzero elements are invertible', () => {
    for (const [p, n] of [
      [787, 2],
      [37, 2],
      [67, 3],
      [1009, 2],
    ] as Array<[number, number]>) {
      const F = GFpn(p, n);
      const a = F.gen();
      // Walk a handful of elements spread across the field.
      for (let k = 1n; k < 40n; k++) {
        const elt = a.mul(F.__call__(k)).add(F.__call__(k * k + 1n));
        if (elt.isZero()) continue;
        expect(elt.mul(elt.inv()).isOne()).toBe(true);
      }
    }
  });
});

describe('Power operations', () => {
  test('negative exponents', () => {
    const F8 = GF(8) as FiniteFieldExtension;
    const a = F8.gen();

    // a^(-1) * a = 1
    const aInv = a.pow(-1);
    expect(a.mul(aInv).isOne()).toBe(true);

    // a^(-2) = (a^2)^(-1)
    expect(a.pow(-2).eq(a.pow(2).inv())).toBe(true);
  });

  test('large exponents', () => {
    const F16 = GF(16) as FiniteFieldExtension;
    const a = F16.gen();

    // a^15 = 1 (multiplicative group order)
    expect(a.pow(15).isOne()).toBe(true);

    // a^16 = a
    expect(a.pow(16).eq(a)).toBe(true);

    // a^1000 = a^(1000 mod 15) = a^10
    expect(a.pow(1000).eq(a.pow(1000n % 15n))).toBe(true);
  });
});

describe('PrimeField generators (M23)', () => {
  // sage: GF(13).gen()  ->  1     (finite_field_prime_modn.py:277)
  // The generator of GF(p) over its prime field is a root of the modulus x - 1,
  // *not* a generator of the multiplicative group.
  test('gen() is one, not a primitive root', () => {
    for (const p of [2n, 3n, 13n, 1009n, 65537n]) {
      expect(new PrimeField(p).gen().value).toBe(1n);
    }
  });

  test('multiplicative_generator() is a primitive root', () => {
    // sage: GF(997).multiplicative_generator()  ->  7
    expect(new PrimeField(997n).multiplicative_generator().value).toBe(7n);
    // sage: GF(1009, modulus='primitive').gen()  ->  11
    expect(new PrimeField(1009n).multiplicative_generator().value).toBe(11n);
    expect(new PrimeField(13n).multiplicative_generator().value).toBe(2n);
    expect(new PrimeField(2n).multiplicative_generator().value).toBe(1n);
  });

  test('primitive_element() is an alias', () => {
    expect(new PrimeField(997n).primitive_element().value).toBe(7n);
  });
});

describe('Conway moduli give genuine fields (C2/M21)', () => {
  // Every one of these entries was reducible or non-Conway before; a reducible
  // modulus makes the quotient a ring with zero divisors.
  const cases: Array<[number, number]> = [
    [29, 2],
    [31, 2],
    [17, 2],
    [13, 4],
    [17, 3],
    [23, 3],
    [31, 3],
    [3, 8],
    [5, 5],
    [5, 6],
    [7, 5],
  ];

  test('every nonzero element is invertible', () => {
    for (const [p, n] of cases) {
      const F = GFpn(p, n);
      let checked = 0;
      for (const elt of F) {
        if (elt.isZero()) continue;
        expect(elt.mul(elt.inv()).isOne()).toBe(true);
        checked++;
      }
      expect(`${p}^${n}:${checked}`).toBe(`${p}^${n}:${Number(F.order) - 1}`);
    }
  });

  test('the generator is a multiplicative generator', () => {
    // Conway polynomials are primitive, so the class of x generates F*.
    for (const [p, n] of [
      [17, 2],
      [31, 3],
      [29, 2],
      [7, 5],
    ] as Array<[number, number]>) {
      const F = GFpn(p, n);
      const g = F.gen();
      const order = F.order - 1n;
      expect(g.pow(order).isOne()).toBe(true);
      // g^(order/q) != 1 for every prime q | order
      const primes = new Set<bigint>();
      let m = order;
      for (let d = 2n; d * d <= m; d++) {
        while (m % d === 0n) {
          primes.add(d);
          m /= d;
        }
      }
      if (m > 1n) primes.add(m);
      for (const q of primes) {
        expect(g.pow(order / q).isOne()).toBe(false);
      }
      // ... and primitiveElement() returns it immediately
      expect(F.primitiveElement().eq(g)).toBe(true);
    }
  });
});

describe('default modulus search (H120/M24)', () => {
  test('GF(2^n) is constructible for every n in [2, 64]', () => {
    const failures: number[] = [];
    for (let n = 2; n <= 64; n++) {
      try {
        GFpn(2, n);
      } catch {
        failures.push(n);
      }
    }
    expect(failures).toEqual([]);
  });

  test("algorithm='first_lexicographic' reproduces the SageMath doctest", () => {
    // sage: GF(19)['x'].irreducible_element(21, algorithm='first_lexicographic')
    // x^21 + x + 5
    //
    // This is no longer the *default*: since the NTL/PARI delegation landed,
    // `irreducible_element(21)` over GF(19) goes to PARI's ffinit, exactly as
    // in SageMath (see 'default modulus delegates ...' below).
    const F = GFpn(19, 21, 'first_lexicographic', 'x');
    expect(F.modulus.toString()).toBe('x^21 + x + 5');
  });

  test('GF(p^2) builds for 64-bit p', () => {
    // Used to trial divide up to sqrt(p^2) = p, i.e. never finished.
    const p = 18446744073709551557n;
    const start = Date.now();
    const F = GF(p * p) as FiniteFieldExtension;
    expect(F.degree).toBe(2);
    expect(F.characteristic).toBe(p);
    expect(Date.now() - start).toBeLessThan(5000);
  });

  test('GF() rejects orders that are not prime powers', () => {
    // 1000003 * 1000033: PARI's trial division gives up on it, so the old
    // is_prime_power-based test called it prime.
    expect(() => GF(1000003n * 1000033n)).toThrow('is not a prime power');
    expect((GF(1000003n * 1000003n) as FiniteFieldExtension).degree).toBe(2);
  });
});

describe('irreducible_element delegates to NTL/PARI (H120)', () => {
  // Port of PolynomialRing_dense_mod_p.irreducible_element
  // (reference/sage/src/sage/rings/polynomial/polynomial_ring.py:3560-3626).
  // Every expected value below was produced by running the real SageMath 10.3.
  const list = (F: FiniteFieldExtension): number[] =>
    F.modulus.coeffs.map((c) => Number(c.value));

  test("default: Conway when available (SageMath's GF(2^8).modulus())", () => {
    // sage: k.<a> = GF(2**8); k.modulus()
    // x^8 + x^4 + x^3 + x^2 + 1            (finite_field_givaro.py:69)
    expect(GFpn(2, 8, undefined, 'x').modulus.toString()).toBe('x^8 + x^4 + x^3 + x^2 + 1');
    // sage: GF(4,'a').modulus()  ->  x^2 + x + 1   (finite_field_constructor.py:399)
    expect(GFpn(2, 2, undefined, 'x').modulus.toString()).toBe('x^2 + x + 1');
    // sage: GF(5)['x'].irreducible_element(4, algorithm='conway').list()
    expect(list(GFpn(5, 4))).toEqual([2, 4, 4, 0, 1]);
  });

  test("default modulus equals SageMath's for every (p, n) our Conway table covers", () => {
    // sage: [[int(c) for c in GF(p**n,'a').modulus().list()] for ...]
    const expected: Array<[number, number, number[]]> = [
      [2, 2, [1, 1, 1]],
      [2, 3, [1, 1, 0, 1]],
      [2, 4, [1, 1, 0, 0, 1]],
      [2, 5, [1, 0, 1, 0, 0, 1]],
      [2, 8, [1, 0, 1, 1, 1, 0, 0, 0, 1]],
      [3, 2, [2, 2, 1]],
      [3, 3, [1, 2, 0, 1]],
      [3, 4, [2, 0, 0, 2, 1]],
      [3, 5, [1, 2, 0, 0, 0, 1]],
      [3, 8, [2, 2, 2, 0, 1, 2, 0, 0, 1]],
      [5, 2, [2, 4, 1]],
      [5, 3, [3, 3, 0, 1]],
      [5, 4, [2, 4, 4, 0, 1]],
      [5, 5, [3, 4, 0, 0, 0, 1]],
      [5, 8, [2, 4, 3, 0, 1, 0, 0, 0, 1]],
      [7, 2, [3, 6, 1]],
      [7, 3, [4, 0, 6, 1]],
      [7, 4, [3, 4, 5, 0, 1]],
      [7, 5, [4, 1, 0, 0, 0, 1]],
      [7, 8, [3, 2, 6, 4, 0, 0, 0, 0, 1]],
      [11, 2, [2, 7, 1]],
      [11, 3, [9, 2, 0, 1]],
      [11, 4, [2, 10, 8, 0, 1]],
      [13, 2, [2, 12, 1]],
      [13, 4, [2, 12, 3, 0, 1]],
      [17, 2, [3, 16, 1]],
      [19, 3, [17, 4, 0, 1]],
      [23, 5, [18, 3, 0, 0, 0, 1]],
      [29, 2, [2, 24, 1]],
      [31, 3, [28, 1, 0, 1]],
    ];
    for (const [p, n, want] of expected) {
      expect(`${p}^${n}:${list(GFpn(p, n))}`).toBe(`${p}^${n}:${want}`);
    }
  });

  test("default: GF(1009^8) — SageMath has no Conway entry either, so both use PARI's ffinit", () => {
    // sage: exists_conway_polynomial(1009, 8)  ->  False
    // sage: [int(c) for c in GF(1009**8,'a').modulus().list()]
    expect(list(GFpn(1009, 8))).toEqual([1, 1005, 999, 10, 15, 1003, 1002, 1, 1]);
  });

  test("algorithm='adleman-lenstra' matches PARI's ffinit exactly", () => {
    // sage: [int(c) for c in GF(p)['x'].irreducible_element(n, algorithm='adleman-lenstra').list()]
    const expected: Array<[number, number, number[]]> = [
      [3, 2, [2, 1, 1]],
      [3, 3, [2, 2, 0, 1]],
      [3, 4, [1, 1, 1, 1, 1]],
      [3, 7, [1, 0, 2, 1, 2, 0, 1, 1]],
      [3, 12, [2, 0, 0, 1, 1, 0, 2, 0, 0, 2, 2, 0, 1]],
      [5, 2, [1, 1, 1]],
      [5, 3, [4, 3, 1, 1]],
      [5, 4, [3, 1, 2, 1, 1]],
      [5, 7, [1, 1, 4, 3, 3, 3, 1, 1]],
      [5, 12, [2, 0, 4, 2, 1, 2, 0, 1, 2, 2, 3, 2, 1]],
      [97, 2, [96, 1, 1]],
      [97, 3, [1, 93, 1, 1]],
      [97, 4, [1, 1, 1, 1, 1]],
      [97, 7, [1, 88, 14, 28, 90, 85, 1, 1]],
      [97, 12, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]],
      [101, 3, [100, 99, 1, 1]],
      [101, 4, [23, 20, 4, 1, 1]],
      [101, 12, [75, 80, 34, 39, 73, 40, 38, 86, 33, 5, 24, 7, 1]],
      [1009, 2, [3, 1, 1]],
      [1009, 7, [1, 1000, 14, 28, 1002, 997, 1, 1]],
      [1009, 12, [969, 284, 694, 118, 911, 41, 405, 402, 772, 899, 2, 7, 1]],
      [2, 21, [1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]],
      [
        19,
        21,
        [3, 17, 16, 12, 8, 12, 1, 13, 15, 14, 13, 11, 11, 5, 9, 6, 14, 6, 4, 13, 10, 1],
      ],
    ];
    for (const [p, n, want] of expected) {
      expect(`${p}^${n}:${list(GFpn(p, n, 'adleman-lenstra'))}`).toBe(`${p}^${n}:${want}`);
    }
  });

  test("GF(19^21)'s default modulus is exactly PARI's ffinit(19, 21)", () => {
    // Our Conway table stops at 19^10, so this takes SageMath's
    // 'adleman-lenstra' branch.  (SageMath itself has a Conway polynomial
    // here -- see the "Finite Fields — Conway Table and Minimal Polynomials" deviation.)
    expect(list(GFpn(19, 21))).toEqual(list(GFpn(19, 21, 'adleman-lenstra')));
  });

  test("algorithm='minimal_weight' matches NTL's GF2X_BuildSparseIrred", () => {
    // sage: [int(c) for c in GF(2)['x'].irreducible_element(n, algorithm='minimal_weight').list()]
    expect(GFpn(2, 8, 'minimal_weight', 'x').modulus.toString()).toBe('x^8 + x^4 + x^3 + x + 1');
    expect(GFpn(2, 33, 'minimal_weight', 'x').modulus.toString()).toBe('x^33 + x^10 + 1');
    expect(GFpn(2, 16, 'minimal_weight', 'x').modulus.toString()).toBe(
      'x^16 + x^5 + x^3 + x + 1'
    );
    expect(GFpn(2, 32, 'minimal_weight', 'x').modulus.toString()).toBe(
      'x^32 + x^7 + x^3 + x^2 + 1'
    );
    expect(GFpn(2, 64, 'minimal_weight', 'x').modulus.toString()).toBe(
      'x^64 + x^4 + x^3 + x + 1'
    );
    expect(GFpn(2, 100, 'minimal_weight', 'x').modulus.toString()).toBe('x^100 + x^15 + 1');
  });

  test("algorithm='first_lexicographic' over GF(2) matches NTL's GF2X_BuildIrred", () => {
    // sage: GF(2)['x'].irreducible_element(33, algorithm='first_lexicographic')
    // x^33 + x^6 + x^3 + x + 1                       (polynomial_ring.py:3541)
    expect(GFpn(2, 33, 'first_lexicographic', 'x').modulus.toString()).toBe(
      'x^33 + x^6 + x^3 + x + 1'
    );
    expect(GFpn(2, 100, 'first_lexicographic', 'x').modulus.toString()).toBe(
      'x^100 + x^6 + x^5 + x^2 + 1'
    );
  });

  test('degree 1 gives x - 1', () => {
    // sage: GF(5)['x'].irreducible_element(1)  ->  x + 4
    expect(GFpn(5, 1, undefined, 'x').modulus.toString()).toBe('x + 4');
    expect(GFpn(97, 1, undefined, 'x').modulus.toString()).toBe('x + 96');
    expect(GFpn(2, 1, undefined, 'x').modulus.toString()).toBe('x + 1');
  });

  test("SageMath's error messages are preserved", () => {
    // sage: GF(5)['x'].irreducible_element(3, algorithm='minimal_weight')
    // NotImplementedError: 'minimal_weight' option only implemented for p = 2
    expect(() => GFpn(5, 3, 'minimal_weight')).toThrow(
      "'minimal_weight' option only implemented for p = 2"
    );
    // sage: GF(5)['x'].irreducible_element(3, algorithm='nosuch')
    // ValueError: no such algorithm for finding an irreducible polynomial: nosuch
    expect(() => GFpn(5, 3, 'nosuch' as never)).toThrow(
      'no such algorithm for finding an irreducible polynomial: nosuch'
    );
    // 'ffprimroot' needs PARI's ffgen/ffprimroot/charpoly, which parigp-ts
    // does not have; it must say so rather than silently returning something else.
    expect(() => GFpn(1009, 4, 'primitive')).toThrow('ffprimroot');
  });

  test("algorithm='primitive' uses the Conway polynomial when there is one", () => {
    // Sage: `if exists_conway_polynomial(p, n): algorithm = "conway"`.
    expect(list(GFpn(5, 4, 'primitive'))).toEqual(list(GFpn(5, 4, 'conway')));
    expect(list(GFpn(5, 4, 'primitive'))).toEqual([2, 4, 4, 0, 1]);
  });

  test("algorithm='random' returns an irreducible polynomial of the right degree", () => {
    for (const [p, n] of [
      [2, 6],
      [3, 5],
      [101, 3],
    ] as Array<[number, number]>) {
      const F = GFpn(p, n, 'random');
      expect(F.modulus.degree()).toBe(n);
      // A reducible modulus makes the quotient a ring with zero divisors.
      for (let k = 1n; k < 30n; k++) {
        const elt = F.gen().mul(F.__call__(k)).add(F.__call__(k * k + 1n));
        if (elt.isZero()) continue;
        expect(elt.mul(elt.inv()).isOne()).toBe(true);
      }
    }
  });

  test('every GF(2^n) default modulus is irreducible for n in [2, 64]', () => {
    // Independent oracle: a monic f of degree n over F_2 is irreducible iff
    // x^(2^n) == x mod f and gcd(x^(2^(n/q)) - x, f) == 1 for each prime q | n.
    // Implemented here directly on bit-packed GF(2) polynomials, so it shares
    // no code with the construction under test.
    const clmul = (a: bigint, b: bigint): bigint => {
      let r = 0n;
      let x = a;
      let y = b;
      while (y !== 0n) {
        if (y & 1n) r ^= x;
        x <<= 1n;
        y >>= 1n;
      }
      return r;
    };
    const bitlen = (x: bigint): number => (x === 0n ? 0 : x.toString(2).length);
    const rem = (a: bigint, b: bigint): bigint => {
      let r = a;
      const db = bitlen(b) - 1;
      while (bitlen(r) - 1 >= db && r !== 0n) r ^= b << BigInt(bitlen(r) - 1 - db);
      return r;
    };
    const gcd = (a: bigint, b: bigint): bigint => {
      while (b !== 0n) [a, b] = [b, rem(a, b)];
      return a;
    };
    const powmod = (e: bigint, f: bigint): bigint => {
      let result = 1n;
      let base = rem(2n, f);
      let k = e;
      while (k > 0n) {
        if (k & 1n) result = rem(clmul(result, base), f);
        base = rem(clmul(base, base), f);
        k >>= 1n;
      }
      return result;
    };
    const primeFactors = (n: number): number[] => {
      const s = new Set<number>();
      let m = n;
      for (let d = 2; d * d <= m; d++)
        while (m % d === 0) {
          s.add(d);
          m /= d;
        }
      if (m > 1) s.add(m);
      return [...s];
    };

    const bad: string[] = [];
    for (let n = 2; n <= 64; n++) {
      const F = GFpn(2, n);
      let f = 0n;
      F.modulus.coeffs.forEach((c, i) => {
        if (c.value === 1n) f |= 1n << BigInt(i);
      });
      if (bitlen(f) - 1 !== n) {
        bad.push(`n=${n}: degree ${bitlen(f) - 1}`);
        continue;
      }
      let ok = powmod(1n << BigInt(n), f) === 2n;
      for (const q of primeFactors(n)) {
        if (gcd(powmod(1n << BigInt(n / q), f) ^ 2n, f) !== 1n) ok = false;
      }
      if (!ok) bad.push(`n=${n}: reducible`);
    }
    expect(bad).toEqual([]);
  });
});
