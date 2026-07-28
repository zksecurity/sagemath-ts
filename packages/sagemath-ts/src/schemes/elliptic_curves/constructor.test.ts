/**
 * Tests for the elliptic curve constructors.
 *
 * Expected values come from SageMath's `sage/schemes/elliptic_curves/constructor.py`
 * (`coefficients_from_j`, `EllipticCurve_from_j`, `EllipticCurve_from_c4c6`).
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import type { FiniteFieldElement } from '../../rings/finite_rings/finite_field_prime.js';
import { EllipticCurve, EllipticCurve_from_c4c6, EllipticCurve_from_j } from './constructor.js';

describe('EllipticCurve', () => {
  it('accepts the short and the long form', () => {
    const K = GF(23n);
    const short = EllipticCurve<FiniteFieldElement>(K, [1n, 1n]);
    expect(short.a_invariants().map(String)).toEqual(['0', '0', '0', '1', '1']);

    const long = EllipticCurve<FiniteFieldElement>(K, [1n, 2n, 3n, 4n, 5n]);
    expect(long.a_invariants().map(String)).toEqual(['1', '2', '3', '4', '5']);
  });

  it('rejects a wrong number of coefficients', () => {
    const K = GF(23n);
    expect(() => EllipticCurve<FiniteFieldElement>(K, [1n, 2n, 3n] as never)).toThrow(
      'Invalid number of coefficients'
    );
  });

  it('rejects singular curves', () => {
    const K = GF(7n);
    expect(() => EllipticCurve<FiniteFieldElement>(K, [0n, 0n])).toThrow(
      'defines a singular curve'
    );
  });
});

describe('EllipticCurve_from_j', () => {
  it('produces a curve with the requested j-invariant over GF(p), p > 3', () => {
    for (const p of [5n, 7n, 11n, 13n, 101n, 1009n]) {
      const K = GF(p);
      for (let j = 0n; j < p; j++) {
        const E = EllipticCurve_from_j<FiniteFieldElement>(K, j);
        expect(E.j_invariant().eq(K.__call__(j))).toBe(true);
      }
    }
  });

  it("matches Sage's coefficients_from_j for the special j-invariants", () => {
    // sage: (general field, char != 2,3)
    //   j == 0    -> [0, 0, 0, 0, 1]
    //   j == 1728 -> [0, 0, 0, 1, 0]
    //   else      -> [0, 0, 0, -3*j*k, -2*j*k^2] with k = j - 1728
    const K = GF(101n);
    expect(EllipticCurve_from_j<FiniteFieldElement>(K, 0n).a_invariants().map(String)).toEqual([
      '0',
      '0',
      '0',
      '0',
      '1',
    ]);
    expect(EllipticCurve_from_j<FiniteFieldElement>(K, 1728n).a_invariants().map(String)).toEqual([
      '0',
      '0',
      '0',
      '1',
      '0',
    ]);

    const j = K.__call__(5n);
    const k = j.sub(K.__call__(1728n));
    const E = EllipticCurve_from_j<FiniteFieldElement>(K, j);
    expect(E.a4().eq(j.mul(k).mul(K.__call__(-3n)))).toBe(true);
    expect(E.a6().eq(j.mul(k).mul(k).mul(K.__call__(-2n)))).toBe(true);
  });

  it("matches Sage's characteristic-2 branch (constructor.py:712-716)", () => {
    // sage: char == 2: j == 0 -> [0, 0, 1, 0, 0]; else -> [1, 0, 0, 0, 1/j]
    const K = GF(2n);
    expect(EllipticCurve_from_j<FiniteFieldElement>(K, 0n).a_invariants().map(String)).toEqual([
      '0',
      '0',
      '1',
      '0',
      '0',
    ]);
    const E1 = EllipticCurve_from_j<FiniteFieldElement>(K, 1n);
    expect(E1.a_invariants().map(String)).toEqual(['1', '0', '0', '0', '1']);
    expect(E1.j_invariant().eq(K.one())).toBe(true);
  });

  it("matches Sage's characteristic-3 branch (constructor.py:717-721)", () => {
    // sage: char == 3: j == 0 -> [0, 0, 0, 1, 0]; else -> [0, j, 0, 0, -j^2]
    const K = GF(3n);
    expect(EllipticCurve_from_j<FiniteFieldElement>(K, 0n).a_invariants().map(String)).toEqual([
      '0',
      '0',
      '0',
      '1',
      '0',
    ]);
    for (const j of [1n, 2n]) {
      const E = EllipticCurve_from_j<FiniteFieldElement>(K, j);
      expect(E.a_invariants().map(String)).toEqual([
        '0',
        String(j),
        '0',
        '0',
        K.__call__(-(j * j)).toString(),
      ]);
      expect(E.j_invariant().eq(K.__call__(j))).toBe(true);
    }
  });
});

describe('EllipticCurve_from_c4c6', () => {
  it('round-trips the c-invariants', () => {
    const K = GF(101n);
    const E = EllipticCurve<FiniteFieldElement>(K, [1n, 2n, 3n, 4n, 5n]);
    const [c4, c6] = E.c_invariants();
    const F = EllipticCurve_from_c4c6<FiniteFieldElement>(K, c4, c6);
    expect(F.c4().eq(c4)).toBe(true);
    expect(F.c6().eq(c6)).toBe(true);
    expect(F.j_invariant().eq(E.j_invariant())).toBe(true);
  });
});
