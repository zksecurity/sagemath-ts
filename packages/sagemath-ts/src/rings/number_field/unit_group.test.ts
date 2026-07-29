/**
 * Tests for unit groups of number fields.
 *
 * The audit of 2026-07 found that this module had no test file at all, which
 * let M28 (`regulator`/`log_embedding` embedding a + b*alpha as a + b*sqrt|disc|)
 * and M29 (`quadraticUnitGroup` branching on disc == -3/-4 and assuming
 * alpha = sqrt(d)) ship.
 */

import { describe, expect, it } from 'bun:test';
import { NotImplementedError } from '../../errors.js';
import { Rational } from '../rational.js';
import { CyclotomicField, NumberFieldConstructor, QuadraticField } from './number_field.js';
import { RationalPolynomial } from './number_field.js';
import {
  S_UnitGroup,
  UnitGroup,
  quadraticUnitGroup,
  realQuadraticFundamentalUnit,
} from './unit_group.js';

describe('quadraticUnitGroup torsion (audit M29)', () => {
  it('finds the 6 roots of unity of Q(sqrt(-3)) whatever the model', () => {
    // x^2 + 3: disc(Z[alpha]) = -12 but disc(K) = -3, so the field is Q(zeta_3)
    // and the torsion subgroup has order 6 (Sage: K.unit_group().zeta_order() == 6).
    for (const poly of [
      [3n, 0n, 1n], // x^2 + 3,   alpha = sqrt(-3)
      [1n, 1n, 1n], // x^2 + x + 1, alpha = zeta_3
      [12n, 0n, 1n], // x^2 + 12,  alpha = 2*sqrt(-3)
    ]) {
      const K = NumberFieldConstructor(poly, 'a');
      const U = quadraticUnitGroup(K);
      expect(U.torsion_order()).toBe(6n);
      const zeta = U.torsion_generator();
      expect(zeta.pow(6n).is_one()).toBe(true);
      expect(zeta.pow(3n).is_one()).toBe(false);
      expect(zeta.pow(2n).is_one()).toBe(false);
      // a root of unity is a unit of O_K: norm +/- 1
      expect(zeta.norm().eq(Rational.one())).toBe(true);
      expect(zeta.is_integral_unit()).toBe(true);
    }
  });

  it('finds the 4 roots of unity of Q(i) whatever the model', () => {
    for (const poly of [
      [1n, 0n, 1n], // x^2 + 1,  alpha = i
      [4n, 0n, 1n], // x^2 + 4,  alpha = 2i
      [9n, 0n, 1n], // x^2 + 9,  alpha = 3i
    ]) {
      const K = NumberFieldConstructor(poly, 'a');
      const U = quadraticUnitGroup(K);
      expect(U.torsion_order()).toBe(4n);
      const zeta = U.torsion_generator();
      expect(zeta.pow(4n).is_one()).toBe(true);
      expect(zeta.pow(2n).is_one()).toBe(false);
      expect(zeta.is_integral_unit()).toBe(true);
    }
  });

  it('has only +/-1 in the other imaginary quadratic fields', () => {
    for (const D of [-5n, -7n, -11n, -15n, -23n]) {
      const K = QuadraticField.create(D);
      const U = quadraticUnitGroup(K);
      expect(U.torsion_order()).toBe(2n);
      expect(U.torsion_generator().eq(K.__call__(-1n))).toBe(true);
      expect(U.rank()).toBe(0);
    }
  });

  it('has rank 1 and torsion 2 for real quadratic fields', () => {
    const K = QuadraticField.create(2n);
    const U = quadraticUnitGroup(K);
    expect(U.rank()).toBe(1);
    expect(U.torsion_order()).toBe(2n);
  });

  it('enumerates the roots of unity', () => {
    const K = QuadraticField.create(-1n);
    const roots = quadraticUnitGroup(K).roots_of_unity();
    expect(roots.length).toBe(4);
    for (const r of roots) {
      expect(r.pow(4n).is_one()).toBe(true);
    }
  });
});

describe('regulator and log embedding (audit M28)', () => {
  it('embeds a + b*alpha at the root of the defining polynomial', () => {
    // Q(sqrt(2)) has fundamental unit 1 + sqrt(2) and regulator
    // log(1 + sqrt(2)) = 0.8813735870...  The port used to evaluate at
    // sqrt(|disc|) = sqrt(8) and report log(1 + sqrt(8)) = 1.3425.
    const K = QuadraticField.create(2n);
    const eps = K.one().add(K.gen());
    const U = new UnitGroup(K, 2n, K.__call__(-1n), [eps]);
    expect(U.regulator()).toBeCloseTo(Math.log(1 + Math.SQRT2), 10);
    expect(U.log_embedding()[0]![0]!).toBeCloseTo(0.8813735870195429, 10);
  });

  it('is independent of the model of the field', () => {
    // Q(sqrt(2)) = Q(sqrt(8)); the fundamental unit 1 + sqrt(2) is
    // 1 + alpha/2 when alpha = sqrt(8).
    const K = QuadraticField.create(8n);
    const eps = K.one().add(K.gen().scalarMul(new Rational(1n, 2n)));
    const U = new UnitGroup(K, 2n, K.__call__(-1n), [eps]);
    expect(U.regulator()).toBeCloseTo(Math.log(1 + Math.SQRT2), 10);
  });

  it('matches Sage for Q(sqrt(5))', () => {
    // Sage: QuadraticField(5).regulator() == 0.481211825059603 = log((1+sqrt5)/2)
    const K = QuadraticField.create(5n);
    const eps = K.one().add(K.gen()).scalarMul(new Rational(1n, 2n));
    const U = new UnitGroup(K, 2n, K.__call__(-1n), [eps]);
    expect(U.regulator()).toBeCloseTo(0.4812118250596035, 10);
  });

  it('returns 1 for a rank-0 unit group', () => {
    const K = QuadraticField.create(-1n);
    expect(quadraticUnitGroup(K).regulator()).toBe(1);
    expect(quadraticUnitGroup(K).log_embedding()).toEqual([]);
  });
});

describe('UnitGroup structure', () => {
  it('reports the Dirichlet rank', () => {
    expect(quadraticUnitGroup(QuadraticField.create(-5n)).rank()).toBe(0);
    expect(quadraticUnitGroup(QuadraticField.create(5n)).rank()).toBe(1);
    // Q(zeta_5): r1 = 0, r2 = 2, rank = 1
    expect(CyclotomicField.create(5n).unit_group().rank()).toBe(1);
  });

  it('exposes zeta(n) for the available roots of unity', () => {
    const K = QuadraticField.create(-1n);
    const U = quadraticUnitGroup(K);
    expect(U.zeta(2n).eq(K.__call__(-1n))).toBe(true);
    expect(U.zeta(4n).pow(4n).is_one()).toBe(true);
    expect(() => U.zeta(3n)).toThrow('no 3-th root of unity in this field');
  });

  it('S_UnitGroup adds |S| to the rank', () => {
    const K = QuadraticField.create(-1n);
    const S = [K.prime_above(5n)];
    const SU = new S_UnitGroup(K, S, 4n, quadraticUnitGroup(K).torsion_generator(), []);
    expect(SU.rank()).toBe(1);
  });
});

describe('fundamental unit of a real quadratic field (audit item 15/M28)', () => {
  // epsilon = u + v*w_D with w_D = quadgen(disc), from PARI's quadunit.
  it('reproduces the classical fundamental units', () => {
    const expected: Array<[bigint, string]> = [
      [2n, 'a + 1'], // 1 + sqrt2
      [3n, 'a + 2'], // 2 + sqrt3
      [5n, '1/2*a + 1/2'], // (1 + sqrt5)/2
      [6n, '2*a + 5'],
      [7n, '3*a + 8'],
      [13n, '1/2*a + 3/2'],
      [17n, 'a + 4'], // PARI: quadunit(17) = 3 + 2*w = 4 + sqrt17
      [94n, '221064*a + 2143295'],
    ];
    for (const [d, str] of expected) {
      const K = QuadraticField.create(d);
      const eps = realQuadraticFundamentalUnit(K);
      expect(eps.toString()).toBe(str);
      // a unit of O_K: an algebraic integer of norm +/-1
      const n = eps.norm();
      expect(n.denominator).toBe(1n);
      expect(n.numerator === 1n || n.numerator === -1n).toBe(true);
      expect(eps.is_integral()).toBe(true);
      expect(eps.is_integral_unit()).toBe(true);
    }
  });

  it('is > 1 at the embedding alpha |-> the larger real root', () => {
    for (const d of [2n, 3n, 5n, 6n, 7n, 11n, 13n, 94n, 1000003n]) {
      const K = QuadraticField.create(d);
      const U = quadraticUnitGroup(K);
      // regulator = log|epsilon| > 0 exactly when |epsilon| > 1
      expect(U.regulator()).toBeGreaterThan(0);
    }
  });

  it('is the *fundamental* unit: no smaller unit > 1 exists', () => {
    // brute force over u + v*w_D for v below the computed v
    // brute force is only feasible while v stays small, so 94 (v = 221064)
    // is covered by the exact-value test above instead
    for (const d of [2n, 3n, 5n, 6n, 7n, 10n, 11n, 13n, 14n, 15n, 21n, 22n, 23n, 46n]) {
      const K = QuadraticField.create(d);
      const D = K.discriminant();
      const eps = realQuadraticFundamentalUnit(K);
      // epsilon = c0 + c1*sqrt(d);  write it over the basis 1, sqrt(d)
      const c0 = eps.list()[0]!;
      const c1 = eps.list()[1]!;
      const target =
        Number(c0.numerator) / Number(c0.denominator) +
        (Number(c1.numerator) / Number(c1.denominator)) * Math.sqrt(Number(d));
      const half = (D & 1n) === 1n;
      const w = half ? (1 + Math.sqrt(Number(D))) / 2 : Math.sqrt(Number(D)) / 2;
      let best = Number.POSITIVE_INFINITY;
      const vmax = Math.ceil(target / w) + 2;
      expect(vmax).toBeLessThan(20000);
      for (let v = 1; v <= vmax; v++) {
        const centre = Math.round(
          v * (half ? (Math.sqrt(Number(D)) - 1) / 2 : Math.sqrt(Number(D)) / 2)
        );
        for (let u = centre - 3; u <= centre + 3; u++) {
          const U = BigInt(u);
          const V = BigInt(v);
          const nrm = half ? U * U + U * V - V * V * ((D - 1n) / 4n) : U * U - V * V * (D / 4n);
          if (nrm !== 1n && nrm !== -1n) continue;
          const val = u + v * w;
          if (val > 1 + 1e-9 && val < best) best = val;
        }
      }
      expect(Number.isFinite(best)).toBe(true);
      expect(Math.abs(best - target)).toBeLessThan(1e-6 * Math.max(1, target));
    }
  });

  it('matches the Sage doctest NumberField(1/2*x^2 - 1/6).units() == (3*a + 2,)', () => {
    // sage/rings/number_field/number_field.py:7207
    const K = NumberFieldConstructor(
      new RationalPolynomial([new Rational(-1n, 6n), Rational.zero(), new Rational(1n, 2n)]),
      'a'
    );
    expect(realQuadraticFundamentalUnit(K).toString()).toBe('3*a + 2');
  });

  it('feeds the unit group, so fundamental_units() no longer throws', () => {
    const K = QuadraticField.create(2n);
    const U = quadraticUnitGroup(K);
    expect(U.rank()).toBe(1);
    expect(U.fundamental_units().length).toBe(1);
    expect(U.gens().length).toBe(2);
  });
});

describe('regulator of a real quadratic field', () => {
  it('matches the Sage doctest NumberField(x^2 - 2).regulator() == 0.881373587019543', () => {
    // sage/rings/number_field/number_field.py:6948
    expect(NumberFieldConstructor([-2n, 0n, 1n], 'a').regulator()).toBeCloseTo(
      0.881373587019543,
      12
    );
  });

  it('matches Sage for several real quadratic fields', () => {
    const expected: Array<[bigint, number]> = [
      [3n, 1.31695789692482],
      [5n, 0.481211825059603],
      [6n, 2.29243166956118],
      [7n, 2.76865938331357],
      [94n, 15.2710021030312],
    ];
    for (const [d, r] of expected) {
      expect(QuadraticField.create(d).regulator()).toBeCloseTo(r, 11);
    }
  });

  it('does not overflow for a fundamental unit with hundreds of digits', () => {
    // eps for Q(sqrt(1000003)) has ~250 digits; evaluating it as a double
    // would give Infinity.
    const R = QuadraticField.create(1000003n).regulator();
    expect(Number.isFinite(R)).toBe(true);
    // R = log(trace(eps)) - log(1 + N*exp(-2R)) with N = +-1
    const eps = quadraticUnitGroup(QuadraticField.create(1000003n)).fundamental_units()[0]!;
    const tr = eps.trace();
    const bits = tr.numerator.toString(2).length;
    const logTr = Math.log(Number(tr.numerator >> BigInt(bits - 53))) + (bits - 53) * Math.LN2;
    const N = Number(eps.norm().numerator);
    expect(R + Math.log(1 + N * Math.exp(-2 * R))).toBeCloseTo(logTr, 9);
  });

  it('is 1 for imaginary quadratic fields, as Sage reports', () => {
    expect(QuadraticField.create(-1n).regulator()).toBe(1);
    expect(QuadraticField.create(-5n).regulator()).toBe(1);
  });
});

describe('UnitGroup.log / exp for rank 1 (real quadratic)', () => {
  it('inverts exp exactly, including large negative exponents', () => {
    // Sage: all(UK.log(u^k) == (0,k) for k in range(10)) is True
    for (const d of [2n, 3n, 5n, 6n, 7n, 13n, 94n]) {
      const K = QuadraticField.create(d);
      const U = quadraticUnitGroup(K);
      const eps = U.fundamental_units()[0]!;
      for (let k = -8; k <= 8; k++) {
        for (const negate of [false, true]) {
          const base = eps.pow(BigInt(k));
          const u = negate ? base.neg() : base;
          const v = U.log(u);
          expect(v).toEqual([negate ? 1n : 0n, BigInt(k)]);
          expect(U.exp(v).eq(u)).toBe(true);
        }
      }
    }
  });

  it('rejects non-units', () => {
    const K = QuadraticField.create(2n);
    const U = quadraticUnitGroup(K);
    expect(() => U.log(K.__call__(2n))).toThrow('element is not a unit');
  });
});

describe('roots of unity of a field of degree > 2 (NumberField.nfrootsof1)', () => {
  // sage: x = polygen(QQ); K.<a> = NumberField(x^4 - 8*x^2 + 36)
  // sage: UK = UnitGroup(K); UK
  // Unit group with structure C4 x Z of
  //  Number Field in a with defining polynomial x^4 - 8*x^2 + 36
  // sage: UK.zeta_order()
  // 4
  // sage: UK.roots_of_unity()
  // [1/12*a^3 - 1/6*a, -1, -1/12*a^3 + 1/6*a, 1]
  it('x^4 - 8x^2 + 36 has C4 torsion', () => {
    const K = NumberFieldConstructor([36n, 0n, -8n, 0n, 1n], 'a');
    const U = K.unit_group();
    expect(U.rank()).toBe(1);
    expect(U.zeta_order()).toBe(4n);
    expect(U.toString()).toBe(
      'Unit group with structure C4 x Z of Number Field in a with defining polynomial x^4 - 8*x^2 + 36'
    );
    expect(U.roots_of_unity().map((z) => z.toString()).sort()).toEqual(
      ['1/12*a^3 - 1/6*a', '-1', '-1/12*a^3 + 1/6*a', '1'].sort()
    );
  }, 60000);

  // sage: K.<a> = NumberField(x^4 - x^2 + 4); U = UnitGroup(K)
  // sage: U.zeta_order()
  // 6
  // sage: U.torsion_generator().value() # random
  // -1/4*a^3 - 1/4*a + 1/2
  it('x^4 - x^2 + 4 has C6 torsion generated by -1/4*a^3 - 1/4*a + 1/2', () => {
    const K = NumberFieldConstructor([4n, 0n, -1n, 0n, 1n], 'a');
    const U = K.unit_group();
    expect(U.zeta_order()).toBe(6n);
    expect(U.torsion_generator().toString()).toBe('-1/4*a^3 - 1/4*a + 1/2');
  }, 60000);

  it('cyclotomic fields: w = n for even n and 2n for odd n', () => {
    const cases: Array<[bigint[], bigint]> = [
      [[1n, 1n, 1n, 1n, 1n], 10n], // Phi_5
      [[1n, 1n, 1n, 1n, 1n, 1n, 1n], 14n], // Phi_7
      [[1n, 0n, 0n, 0n, 1n], 8n], // Phi_8
      [[1n, 0n, 0n, 1n, 0n, 0n, 1n], 18n], // Phi_9
      [[1n, 0n, -1n, 0n, 1n], 12n], // Phi_12
    ];
    for (const [poly, w] of cases) {
      const K = NumberFieldConstructor(poly, 'a');
      const r = K.nfrootsof1()!;
      expect(r.order).toBe(w);
      // the generator really has order exactly w
      expect(r.generator.pow(w).eq(K.one())).toBe(true);
      for (const q of [2n, 3n, 5n, 7n]) {
        if (w % q === 0n) expect(r.generator.pow(w / q).eq(K.one())).toBe(false);
      }
    }
  }, 120000);

  it('a field with a real embedding has exactly {1, -1}', () => {
    for (const poly of [
      [-2n, 0n, 0n, 1n],
      [2n, 0n, 0n, 1n],
      [-1n, -1n, 0n, 0n, 0n, 1n],
    ] as bigint[][]) {
      const K = NumberFieldConstructor(poly, 'a');
      const r = K.nfrootsof1()!;
      expect(r.order).toBe(2n);
      expect(r.generator.eq(K.__call__(-1n))).toBe(true);
    }
  });

  it('never claims a torsion order it has not proved', () => {
    // fundamental units still need bnfinit
    const K = NumberFieldConstructor([-2n, 0n, 0n, 1n], 'a');
    expect(() => K.unit_group().fundamental_units()).toThrow(NotImplementedError);
  });
});
