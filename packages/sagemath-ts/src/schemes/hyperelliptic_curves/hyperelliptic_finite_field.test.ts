/**
 * Tests for `HyperellipticCurve_finite_field`.
 *
 * Every expected value below was obtained by running the installed SageMath on
 * the same input; the doctest values of
 * `sage/schemes/hyperelliptic_curves/hyperelliptic_finite_field.py` are
 * reproduced verbatim where they exist.
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import { GFpn } from '../../rings/finite_rings/finite_field_extension.js';
import type { FiniteFieldElement } from '../../rings/finite_rings/finite_field_prime.js';
import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import { HyperellipticCurve } from './constructor.js';
import { zz_poly_repr } from './hyperelliptic_finite_field.js';

type FF = FiniteFieldElement;

// biome-ignore lint/suspicious/noExplicitAny: the constructor returns a dynamic subclass
type AnyCurve = any;

function fp(p: bigint, name = 't') {
  const K = GF(p);
  const R = new PolynomialRing<FF>(K, name);
  return { K, R, t: R.gen(), c: (n: bigint) => R.__call__(K.__call__(n)) };
}

describe('_frobenius_coefficient_bound_charpoly (py:73-125)', () => {
  it('matches the doctests over GF(37)', () => {
    const { t, c } = fp(37n);
    const bounds = [3, 5, 7].map((d) =>
      (
        HyperellipticCurve(t.pow(d).add(t).add(c(1n))) as AnyCurve
      )._frobenius_coefficient_bound_charpoly()
    );
    expect(bounds).toEqual([1, 2, 3]);
  });

  it('matches the doctests over GF(next_prime(10^9))', () => {
    const { t, c } = fp(1000000007n);
    const bounds = [3, 5, 7, 9, 11, 13].map((d) =>
      (
        HyperellipticCurve(t.pow(d).add(t).add(c(1n))) as AnyCurve
      )._frobenius_coefficient_bound_charpoly()
    );
    expect(bounds).toEqual([1, 2, 2, 3, 3, 4]);
  });
});

describe('_frobenius_coefficient_bound_traces (py:127-180)', () => {
  it('matches the doctests over GF(37)', () => {
    const { t, c } = fp(37n);
    const bounds = [3, 5, 7].map((d) =>
      (
        HyperellipticCurve(t.pow(d).add(t).add(c(1n))) as AnyCurve
      )._frobenius_coefficient_bound_traces()
    );
    expect(bounds).toEqual([1, 2, 2]);
  });

  it('matches the doctests over GF(next_prime(10^9))', () => {
    const { t, c } = fp(1000000007n);
    const H = (d: number) => HyperellipticCurve(t.pow(d).add(t).add(c(1n))) as AnyCurve;
    expect([3, 5, 7].map((d) => H(d)._frobenius_coefficient_bound_traces())).toEqual([1, 1, 1]);
    expect(H(9)._frobenius_coefficient_bound_traces(3)).toBe(2);
    expect(H(11)._frobenius_coefficient_bound_traces(3)).toBe(2);
    expect(H(13)._frobenius_coefficient_bound_traces(5)).toBe(3);
  });

  it('matches the issue 18831 doctest over GF(11)', () => {
    const { t, c } = fp(11n);
    const H = HyperellipticCurve(t.pow(5).sub(t).add(c(1n))) as AnyCurve;
    expect(H._frobenius_coefficient_bound_traces()).toBe(2);
  });
});

describe('frobenius_polynomial and zeta_function', () => {
  it('y^2 = t^5 + t + 2 over GF(37) (py:316-325, 1456-1466)', () => {
    const { R, K, t, c } = fp(37n);
    const H = HyperellipticCurve(t.pow(5).add(t).add(c(2n))) as AnyCurve;
    expect(zz_poly_repr(H.frobenius_polynomial())).toBe('x^4 + x^3 - 52*x^2 + 37*x + 1369');
    expect(String(H.zeta_function())).toBe(
      '(1369*x^4 + 37*x^3 - 52*x^2 + x + 1)/(37*x^2 - 38*x + 1)'
    );
    const two = K.__call__(2n);
    const twist = HyperellipticCurve(
      t.pow(5).scalar_mul(two).add(t.scalar_mul(two)).add(c(4n))
    ) as AnyCurve;
    void R;
    expect(zz_poly_repr(twist.frobenius_polynomial())).toBe('x^4 - x^3 - 52*x^2 - 37*x + 1369');
    expect(String(twist.zeta_function())).toBe(
      '(1369*x^4 - 37*x^3 - 52*x^2 - x + 1)/(37*x^2 - 38*x + 1)'
    );
  });

  it('the Jacobian cardinality is frobenius_polynomial(1)', () => {
    const { t, c } = fp(37n);
    const H = HyperellipticCurve(t.pow(5).add(t).add(c(2n))) as AnyCurve;
    // 1 + 1 - 52 + 37 + 1369
    expect(H.jacobian().cardinality()).toBe(1356n);
    expect(H.jacobian().dimension()).toBe(2n);
  });

  it('y^2 = t^9 + t^3 + 1 over GF(5)', () => {
    const { t, c } = fp(5n);
    const H = HyperellipticCurve(t.pow(9).add(t.pow(3)).add(c(1n))) as AnyCurve;
    expect(zz_poly_repr(H.frobenius_polynomial())).toBe(
      'x^8 + 3*x^7 + 5*x^6 + 125*x^2 + 375*x + 625'
    );
  });

  it('char 2: y^2 + t^4 y = t^9 + t over GF(2) (py:1441-1444)', () => {
    const K = GF(2n);
    const R = new PolynomialRing<FF>(K, 't');
    const t = R.gen();
    const H = HyperellipticCurve(t.pow(9).add(t), t.pow(4)) as AnyCurve;
    expect(String(H.zeta_function())).toBe(
      '(16*x^8 + 8*x^7 + 8*x^6 + 4*x^5 + 6*x^4 + 2*x^3 + 2*x^2 + x + 1)/(2*x^2 - 3*x + 1)'
    );
  });

  it('char 2 over GF(4) (py:1446-1449)', () => {
    const F4 = GFpn(2n, 2, undefined, 'a');
    const R = new PolynomialRing(F4, 't');
    const t = R.gen();
    const c = (n: bigint) => R.__call__(F4.__call__(n));
    const H = HyperellipticCurve(
      t.pow(5).add(t.pow(3)).add(t.pow(2)).add(t).add(c(1n)),
      t.pow(2).add(t).add(c(1n))
    ) as AnyCurve;
    expect(String(H.zeta_function())).toBe('(16*x^4 + 8*x^3 + x^2 + 2*x + 1)/(4*x^2 - 5*x + 1)');
  });

  it('rejects an unknown algorithm and reports the unported ones', () => {
    const { t, c } = fp(37n);
    const H = HyperellipticCurve(t.pow(5).add(t).add(c(2n))) as AnyCurve;
    expect(() => H.frobenius_polynomial({ algorithm: 'bogus' })).toThrow('unknown algorithm bogus');
    expect(() => H.frobenius_polynomial({ algorithm: 'matrix' })).toThrow(
      'SAGE_NOT_IMPLEMENTED: frobenius_polynomial(algorithm="matrix")'
    );
    expect(() => H.frobenius_polynomial({ algorithm: 'pari' })).toThrow(
      'SAGE_NOT_IMPLEMENTED: frobenius_polynomial(algorithm="pari")'
    );
    expect(() => H.frobenius_matrix()).toThrow('SAGE_NOT_IMPLEMENTED: frobenius_matrix');
    expect(() => H.frobenius_matrix_hypellfrob()).toThrow(
      'SAGE_NOT_IMPLEMENTED: frobenius_matrix_hypellfrob'
    );
    expect(() => H.count_points_hypellfrob()).toThrow(
      'SAGE_NOT_IMPLEMENTED: count_points_hypellfrob'
    );
    expect(() => H.cardinality_hypellfrob()).toThrow(
      'SAGE_NOT_IMPLEMENTED: cardinality_hypellfrob'
    );
    expect(() => H.count_points_matrix_traces()).toThrow(
      'SAGE_NOT_IMPLEMENTED: count_points_matrix_traces'
    );
  });
});

describe('count_points and cardinality', () => {
  it('y^2 = t^5 + t + 2 over GF(37)', () => {
    const { t, c } = fp(37n);
    const H = HyperellipticCurve(t.pow(5).add(t).add(c(2n))) as AnyCurve;
    expect(H.count_points(4)).toEqual([39n, 1265n, 50922n, 1873873n]);
    expect(H.cardinality()).toBe(39n);
    expect(H.cardinality(2)).toBe(1265n);
    expect(H.cardinality(3)).toBe(50922n);
  });

  it('count_points_exhaustive over GF(5) (py:1021-1054)', () => {
    const { t, c } = fp(5n);
    const H = HyperellipticCurve(t.pow(9).add(t.pow(3)).add(c(1n))) as AnyCurve;
    expect(H.count_points_exhaustive(5)).toEqual([9n, 27n, 108n, 675n, 3069n]);
    expect(H.count_points_exhaustive(15)).toEqual([
      9n,
      27n,
      108n,
      675n,
      3069n,
      16302n,
      78633n,
      389475n,
      1954044n,
      9768627n,
      48814533n,
      244072650n,
      1220693769n,
      6103414827n,
      30517927308n,
    ]);
  });

  it('count_points over GF(3) and GF(9) (py:1178-1216)', () => {
    const { t, c } = fp(3n, 'x');
    const H = HyperellipticCurve(t.pow(3).add(t.pow(2)).add(c(1n))) as AnyCurve;
    expect(H.count_points(4)).toEqual([6n, 12n, 18n, 96n]);

    const F9 = GFpn(3n, 2, undefined, 'a');
    const R9 = new PolynomialRing(F9, 'x');
    const x9 = R9.gen();
    const c9 = (n: bigint) => R9.__call__(F9.__call__(n));
    const H9 = HyperellipticCurve(x9.pow(5).add(x9.pow(2)).add(c9(1n))) as AnyCurve;
    // Sage: H.count_points(5) == [18, 78, 738, 6366, 60018]
    expect(H9.count_points(3)).toEqual([18n, 78n, 738n]);
  });

  it('count_points over GF(4) with h nonzero (py:1213-1216)', () => {
    const F4 = GFpn(2n, 2, undefined, 'a');
    const a = F4.gen();
    const R = new PolynomialRing(F4, 'x');
    const x = R.gen();
    const c = (n: bigint) => R.__call__(F4.__call__(n));
    const H = HyperellipticCurve(
      x.pow(5).add(x.pow(2).scalar_mul(a)).add(c(1n)),
      x.add(R.__call__(a.add(F4.one())))
    ) as AnyCurve;
    expect(H.count_points(6)).toEqual([2n, 24n, 74n, 256n, 1082n, 4272n]);
  });

  it('char 2 over GF(2^5)', () => {
    const K = GFpn(2n, 5, undefined, 'e');
    const e = K.gen();
    const R = new PolynomialRing(K, 't');
    const t = R.gen();
    const H = HyperellipticCurve(
      t
        .pow(5)
        .add(t.scalar_mul(e))
        .add(R.__call__(e.pow(3))),
      t
    ) as AnyCurve;
    expect(H.cardinality()).toBe(32n);
    expect(H.count_points(2)).toEqual([32n, 1056n]);
    // py:591-595
    expect(zz_poly_repr(H.frobenius_polynomial())).toBe('x^4 - x^3 + 16*x^2 - 32*x + 1024');
  });

  it('issue 20391: even degree models', () => {
    const { t, c } = fp(23n, 'x');
    expect((HyperellipticCurve(t.pow(8).add(c(1n))) as AnyCurve).cardinality()).toBe(24n);
    const { t: u, c: cu } = fp(4099n, 'x');
    expect((HyperellipticCurve(u.pow(6).add(u).add(cu(1n))) as AnyCurve).count_points(1)).toEqual([
      4106n,
    ]);
  });

  it('issue 19122', () => {
    const { K, t, c } = fp(19n, 'x');
    const f = t
      .pow(4)
      .scalar_mul(K.__call__(15n))
      .add(t.pow(3).scalar_mul(K.__call__(7n)))
      .add(t.pow(2).scalar_mul(K.__call__(3n)))
      .add(t.scalar_mul(K.__call__(7n)))
      .add(c(18n));
    expect((HyperellipticCurve(f) as AnyCurve).cardinality_exhaustive(1)).toBe(19n);
  });

  it('issue 21195: points at infinity on genus 1 curves', () => {
    for (const p of [2n, 3n]) {
      const { t } = fp(p, 'z');
      const H = HyperellipticCurve(t.pow(2).neg().add(t), t.pow(2)) as AnyCurve;
      expect(H.count_points_exhaustive()).toEqual([5n]);
    }
  });
});

describe('points (py:632-909)', () => {
  it('y^2 = x^7 - x^2 - 1 over GF(7) (py:847-851)', () => {
    const { t, c } = fp(7n, 'x');
    const C = HyperellipticCurve(t.pow(7).sub(t.pow(2)).sub(c(1n))) as AnyCurve;
    expect(C.points().map(String)).toEqual([
      '(0 : 1 : 0)',
      '(2 : 5 : 1)',
      '(2 : 2 : 1)',
      '(3 : 0 : 1)',
      '(4 : 6 : 1)',
      '(4 : 1 : 1)',
      '(5 : 0 : 1)',
      '(6 : 5 : 1)',
      '(6 : 2 : 1)',
    ]);
  });

  it('_points_cache_sqrt over GF(7) (py:764-770)', () => {
    const { t, c } = fp(7n, 'x');
    const C = HyperellipticCurve(t.pow(3).add(t.pow(2)).sub(c(1n))) as AnyCurve;
    expect(C._points_cache_sqrt().map(String)).toEqual([
      '(0 : 1 : 0)',
      '(1 : 6 : 1)',
      '(1 : 1 : 1)',
      '(2 : 5 : 1)',
      '(2 : 2 : 1)',
      '(3 : 0 : 1)',
      '(4 : 4 : 1)',
      '(4 : 3 : 1)',
      '(5 : 4 : 1)',
      '(5 : 3 : 1)',
    ]);
    const bf = C._points_cache_sqrt({ brute_force: true }).map(String);
    expect(new Set(bf)).toEqual(new Set(C._points_cache_sqrt().map(String)));
  });

  it('conics are allowed (py:860-867)', () => {
    const { K, t, c } = fp(7n, 'x');
    const H = HyperellipticCurve(
      t
        .pow(2)
        .scalar_mul(K.__call__(3n))
        .add(t.scalar_mul(K.__call__(5n)))
        .add(c(1n))
    ) as AnyCurve;
    expect(H.points().map(String)).toEqual([
      '(0 : 6 : 1)',
      '(0 : 1 : 1)',
      '(1 : 4 : 1)',
      '(1 : 3 : 1)',
      '(2 : 4 : 1)',
      '(2 : 3 : 1)',
      '(3 : 6 : 1)',
      '(3 : 1 : 1)',
    ]);
  });

  it('a genus 2 curve over GF(11) (py:880-884)', () => {
    const { t, c } = fp(11n, 'x');
    let f = t;
    for (let i = 1n; i <= 5n; i++) {
      f = f.mul(t.add(c(i)));
    }
    const H = HyperellipticCurve(f) as AnyCurve;
    expect(H.points().map(String)).toEqual([
      '(0 : 1 : 0)',
      '(0 : 0 : 1)',
      '(1 : 7 : 1)',
      '(1 : 4 : 1)',
      '(5 : 7 : 1)',
      '(5 : 4 : 1)',
      '(6 : 0 : 1)',
      '(7 : 0 : 1)',
      '(8 : 0 : 1)',
      '(9 : 0 : 1)',
      '(10 : 0 : 1)',
    ]);
  });

  it('counts points on extension fields (py:645-649, 855-858)', () => {
    const F49 = GFpn(7n, 2, undefined, 'a');
    const R = new PolynomialRing(F49, 'x');
    const x = R.gen();
    const c = (n: bigint) => R.__call__(F49.__call__(n));
    const C = HyperellipticCurve(
      x.pow(5).sub(x.pow(2)).sub(c(1n)),
      x.pow(2).add(R.__call__(F49.gen()))
    ) as AnyCurve;
    expect(C._points_fast_sqrt().length).toBe(31);

    const F121 = GFpn(11n, 2, undefined, 'a');
    const R2 = new PolynomialRing(F121, 'x');
    const x2 = R2.gen();
    const c2 = (n: bigint) => R2.__call__(F121.__call__(n));
    const C2 = HyperellipticCurve(x2.pow(5).add(x2).sub(c2(1n)), x2.pow(2).add(c2(2n))) as AnyCurve;
    expect(C2.points().length).toBe(122);
  });

  it('cardinality_exhaustive over GF(9) (py:1250-1254)', () => {
    const F9 = GFpn(3n, 2, undefined, 'a');
    const R = new PolynomialRing(F9, 'x');
    const x = R.gen();
    const c = (n: bigint) => R.__call__(F9.__call__(n));
    const C = HyperellipticCurve(
      x.pow(7).sub(c(1n)),
      x.pow(2).add(R.__call__(F9.gen()))
    ) as AnyCurve;
    expect(C.cardinality_exhaustive()).toBe(7n);
  });

  it('cardinality_exhaustive over GF(1031) (py:1256-1260)', () => {
    const { K, t, c } = fp(1031n);
    const H = HyperellipticCurve(
      t
        .pow(7)
        .add(t.pow(5).scalar_mul(K.__call__(3n)))
        .add(c(5n))
    ) as AnyCurve;
    expect(H.cardinality_exhaustive()).toBe(1025n);
  });
});

describe('Cartier matrix, Hasse-Witt, a-number and p-rank (py:1475-1947)', () => {
  const mat = (M: unknown[][]) => M.map((row) => row.map(String).join(' '));

  it('y^2 = x^7 - 1 over GF(9) (py:1629-1634, 1827-1832)', () => {
    const F9 = GFpn(3n, 2, undefined, 'x');
    const R = new PolynomialRing(F9, 'y');
    const y = R.gen();
    const c = (n: bigint) => R.__call__(F9.__call__(n));
    const C = HyperellipticCurve(y.pow(7).sub(c(1n)), R.zero()) as AnyCurve;
    expect(mat(C.Cartier_matrix())).toEqual(['0 0 2', '0 0 0', '0 1 0']);
    expect(mat(C.Hasse_Witt())).toEqual(['0 0 0', '0 0 0', '0 0 0']);
    expect(C.a_number()).toBe(1);
    expect(C.p_rank()).toBe(0);
  });

  it('y^2 = x^5 + 1 over GF(49) (py:1636-1640, 1834-1838)', () => {
    const F49 = GFpn(7n, 2, undefined, 'x');
    const R = new PolynomialRing(F49, 'y');
    const y = R.gen();
    const c = (n: bigint) => R.__call__(F49.__call__(n));
    const C = HyperellipticCurve(y.pow(5).add(c(1n)), R.zero()) as AnyCurve;
    expect(mat(C.Cartier_matrix())).toEqual(['0 3', '0 0']);
    expect(mat(C.Hasse_Witt())).toEqual(['0 0', '0 0']);
    expect(C.a_number()).toBe(1);
    expect(C.p_rank()).toBe(0);
  });

  it('y^2 = x^7 - 1 over GF(5)', () => {
    const { R, t, c } = fp(5n, 'y');
    const C = HyperellipticCurve(t.pow(7).sub(c(1n)), R.zero()) as AnyCurve;
    expect(mat(C.Cartier_matrix())).toEqual(['0 0 0', '0 0 3', '1 0 0']);
    expect(mat(C.Hasse_Witt())).toEqual(['0 0 0', '0 0 0', '0 0 0']);
    expect(C.a_number()).toBe(1);
    expect(C.p_rank()).toBe(0);
  });

  it('y^2 = x^5 + x + 1 over GF(11) is ordinary', () => {
    const { t, c } = fp(11n, 'y');
    const C = HyperellipticCurve(t.pow(5).add(t).add(c(1n))) as AnyCurve;
    expect(mat(C.Cartier_matrix())).toEqual(['10 5', '5 5']);
    expect(mat(C.Hasse_Witt())).toEqual(['4 9', '9 6']);
    expect(C.a_number()).toBe(0);
    expect(C.p_rank()).toBe(2);
  });

  it('rejects the cases upstream rejects (py:1530-1558)', () => {
    const K2 = GF(2n);
    const R2 = new PolynomialRing<FF>(K2, 'x');
    const x2 = R2.gen();
    const C2 = HyperellipticCurve(x2.pow(7).sub(R2.one()), x2) as AnyCurve;
    expect(() => C2.Cartier_matrix()).toThrow('p must be odd');

    const { R, t, c } = fp(5n, 'x');
    const C3 = HyperellipticCurve(t.pow(7).sub(c(1n)), c(4n)) as AnyCurve;
    expect(() => C3.Cartier_matrix()).toThrow('E must be of the form y^2 = f(x)');

    const C4 = HyperellipticCurve(t.pow(8).sub(c(1n)), R.zero()) as AnyCurve;
    expect(() => C4.Cartier_matrix()).toThrow('In this implementation the degree of f must be odd');

    const C5 = HyperellipticCurve(t.pow(5).add(c(1n)), R.zero(), {
      check_squarefree: false,
    }) as AnyCurve;
    expect(() => C5.Cartier_matrix()).toThrow('curve is not smooth');
  });
});
