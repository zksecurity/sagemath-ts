/**
 * Tests for unit groups of number fields.
 *
 * The audit of 2026-07 found that this module had no test file at all, which
 * let M28 (`regulator`/`log_embedding` embedding a + b*alpha as a + b*sqrt|disc|)
 * and M29 (`quadraticUnitGroup` branching on disc == -3/-4 and assuming
 * alpha = sqrt(d)) ship.
 */

import { describe, expect, it } from 'bun:test';
import { Rational } from '../rational.js';
import { CyclotomicField, NumberFieldConstructor, QuadraticField } from './number_field.js';
import { S_UnitGroup, UnitGroup, quadraticUnitGroup } from './unit_group.js';

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
