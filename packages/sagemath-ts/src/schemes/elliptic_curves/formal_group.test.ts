/**
 * Tests for the formal group of an elliptic curve.
 *
 * The expected values are SageMath's own doctests from
 * `sage/schemes/elliptic_curves/formal_group.py`. Doctests stated over QQ are
 * reproduced over a large prime field (GF(1000003)) by reducing the rational
 * coefficients; doctests already stated over a finite field are used verbatim.
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import { EllipticCurve } from './constructor.js';

/** A prime large enough that the small rational doctest coefficients reduce injectively. */
const BIGP = 1000003n;

function reduce(a: bigint, b: bigint = 1n): string {
  const K = GF(BIGP);
  return K.__call__(a).div(K.__call__(b)).toString();
}

/** Coefficients [c_0, ..., c_{prec-1}] of a power series, as strings. */
function coeffs(f: { list(): Array<{ toString(): string }> }, prec: number): string[] {
  const list = f.list();
  const out: string[] = [];
  for (let i = 0; i < prec; i++) {
    out.push((list[i] ?? { toString: () => '0' }).toString());
  }
  return out;
}

describe('EllipticCurveFormalGroup', () => {
  describe('w', () => {
    it("matches Sage's doctest for EllipticCurve([0,0,1,-1,0]).formal_group().w(10)", () => {
      // sage: EllipticCurve([0,0,1,-1,0]).formal_group().w(10)
      //  t^3 + t^6 - t^7 + 2*t^9 + O(t^10)
      const K = GF(BIGP);
      const E = EllipticCurve(K, [0n, 0n, 1n, -1n, 0n]);
      const w = E.formal_group().w(10);
      const c = coeffs(w, 10);
      const expected = ['0', '0', '0', '1', '0', '0', '1', reduce(-1n), '0', '2'];
      expect(c).toEqual(expected);
    });
  });

  describe('x and y', () => {
    it("matches Sage's x(10) for EllipticCurve([0,0,1,-1,0])", () => {
      // sage: EllipticCurve([0, 0, 1, -1, 0]).formal_group().x(10)
      //  t^-2 - t + t^2 - t^4 + 2*t^5 - t^6 - 2*t^7 + 6*t^8 - 6*t^9 + O(t^10)
      const K = GF(BIGP);
      const E = EllipticCurve(K, [0n, 0n, 1n, -1n, 0n]);
      const fg = E.formal_group();
      const x = fg.x(10);
      expect(x.valuation()).toBe(-2);
      expect(x.residue().toString()).toBe('0'); // no t^-1 term
      // x(t) = t^-2 * u(t) with u = (w/t^3)^-1; the port's LaurentSeriesElement
      // exposes no coefficient accessor for negative valuations, so compare the
      // unit part against Sage's coefficients of t^-2, t^-1, t^0, ...
      const u = fg.w(20)._shiftRight(3).add_bigoh(12).inv();
      const expected = [1n, 0n, 0n, -1n, 1n, 0n, -1n, 2n, -1n, -2n, 6n, -6n];
      expect(coeffs(u, expected.length)).toEqual(expected.map((v) => reduce(v)));
    });

    it("matches Sage's y(10) for EllipticCurve([0,0,1,-1,0])", () => {
      // sage: EllipticCurve([0, 0, 1, -1, 0]).formal_group().y(10)
      //  -t^-3 + 1 - t + t^3 - 2*t^4 + t^5 + 2*t^6 - 6*t^7 + 6*t^8 + 3*t^9 + O(t^10)
      const K = GF(BIGP);
      const E = EllipticCurve(K, [0n, 0n, 1n, -1n, 0n]);
      const fg = E.formal_group();
      const y = fg.y(10);
      expect(y.valuation()).toBe(-3);
      expect(y.residue().toString()).toBe('0'); // no t^-1 term
      // y(t) = -t^-3 * u(t); compare the unit part against Sage's coefficients
      // of t^-3, t^-2, t^-1, t^0, ...
      const u = fg.w(20)._shiftRight(3).add_bigoh(13).inv().neg();
      const expected = [-1n, 0n, 0n, 1n, -1n, 0n, 1n, -2n, 1n, 2n, -6n, 6n, 3n];
      expect(coeffs(u, expected.length)).toEqual(expected.map((v) => reduce(v)));
    });
  });

  describe('differential', () => {
    it("matches Sage's doctest over GF(53) for EllipticCurve([-1, 1/4])", () => {
      // sage: EllipticCurve(Integers(53), [-1, 1/4]).formal_group().differential(15)
      //  1 + 51*t^4 + 14*t^6 + 6*t^8 + 48*t^10 + 24*t^12 + 13*t^14 + O(t^15)
      const K = GF(53n);
      const E = EllipticCurve(K, [K.__call__(-1n), K.__call__(4n).inv()]);
      const d = E.formal_group().differential(15);
      const c = coeffs(d, 15);
      expect(c).toEqual([
        '1',
        '0',
        '0',
        '0',
        '51',
        '0',
        '14',
        '0',
        '6',
        '0',
        '48',
        '0',
        '24',
        '0',
        '13',
      ]);
    });

    it("matches Sage's rational doctest for EllipticCurve([-1, 1/4]) reduced mod p", () => {
      // sage: EllipticCurve([-1, 1/4]).formal_group().differential(15)
      //  1 - 2*t^4 + 3/4*t^6 + 6*t^8 - 5*t^10 - 305/16*t^12 + 105/4*t^14 + O(t^15)
      const K = GF(BIGP);
      const E = EllipticCurve(K, [K.__call__(-1n), K.__call__(4n).inv()]);
      const c = coeffs(E.formal_group().differential(15), 15);
      expect(c[0]).toBe('1');
      expect(c[4]).toBe(reduce(-2n));
      expect(c[6]).toBe(reduce(3n, 4n));
      expect(c[8]).toBe(reduce(6n));
      expect(c[10]).toBe(reduce(-5n));
      expect(c[12]).toBe(reduce(-305n, 16n));
      expect(c[14]).toBe(reduce(105n, 4n));
      for (const odd of [1, 2, 3, 5, 7, 9, 11, 13]) {
        expect(c[odd]).toBe('0');
      }
    });

    it('is not zero-padded beyond t^3 (audit H91)', () => {
      // The old implementation emitted only four hardcoded coefficients and
      // padded with zeros; over GF(7) the true series is 1 + 6*t^4 + 5*t^6 + ...
      const K = GF(7n);
      const E = EllipticCurve(K, [3n, 4n]);
      const c = coeffs(E.formal_group().differential(10), 10);
      expect(c[4]).toBe('6');
      expect(c[6]).toBe('5');
    });
  });

  describe('log', () => {
    it("matches Sage's doctest for EllipticCurve([-1, 1/4]).formal_group().log(15)", () => {
      // sage: EllipticCurve([-1, 1/4]).formal_group().log(15)
      //  t - 2/5*t^5 + 3/28*t^7 + 2/3*t^9 - 5/11*t^11 - 305/208*t^13 + O(t^15)
      const K = GF(BIGP);
      const E = EllipticCurve(K, [K.__call__(-1n), K.__call__(4n).inv()]);
      const c = coeffs(E.formal_group().log(15), 15);
      expect(c[1]).toBe('1');
      expect(c[5]).toBe(reduce(-2n, 5n));
      expect(c[7]).toBe(reduce(3n, 28n));
      expect(c[9]).toBe(reduce(2n, 3n));
      expect(c[11]).toBe(reduce(-5n, 11n));
      expect(c[13]).toBe(reduce(-305n, 208n));
    });
  });

  describe('inverse', () => {
    it("matches Sage's generic doctest i(6) for E = [a1,a2,a3,a4,a6]", () => {
      // sage: E.formal_group().inverse(6)
      //  -t - a1*t^2 - a1^2*t^3 + (-a1^3 - a3)*t^4 + (-a1^4 - 3*a1*a3)*t^5 + O(t^6)
      const K = GF(BIGP);
      const [a1, a3] = [1n, 3n];
      const E = EllipticCurve(K, [a1, 2n, a3, 4n, 6n]);
      const c = coeffs(E.formal_group().inverse(6), 6);
      expect(c).toEqual([
        '0',
        reduce(-1n),
        reduce(-a1),
        reduce(-(a1 * a1)),
        reduce(-(a1 ** 3n + a3)),
        reduce(-(a1 ** 4n + 3n * a1 * a3)),
      ]);
    });

    it('satisfies F(t, i(t)) = 0', () => {
      const K = GF(7n);
      const E = EllipticCurve(K, [3n, 4n]);
      const fg = E.formal_group();
      const prec = 7;
      const F = fg.group_law(prec);
      const i = fg.inverse(prec);
      // Substitute t1 -> t, t2 -> i(t) and check the result vanishes.
      const R = i.parent();
      const t = R.gen().add_bigoh(prec);
      let acc = R.zero().add_bigoh(prec);
      for (const [key, c] of F.terms) {
        const [p, q] = key.split(',').map(Number) as [number, number];
        acc = acc.add(t.pow(p).mul(i.pow(q))._scalarMul(c)).add_bigoh(prec);
      }
      expect(acc.is_zero()).toBe(true);
    });
  });

  describe('group_law', () => {
    it("matches Sage's doctests for EllipticCurve(GF(7), [3, 4])", () => {
      // sage: ehat.group_law(3) -> t1 + t2 + O(t1, t2)^3
      // sage: ehat.group_law(7) ->
      //   t1 + t2 + t1^4*t2 + 2*t1^3*t2^2 + 2*t1^2*t2^3 + t1*t2^4 + O(t1, t2)^7
      const K = GF(7n);
      const E = EllipticCurve(K, [3n, 4n]);
      const fg = E.formal_group();
      expect(fg.group_law(3).toString()).toBe('t1 + t2 + O(t1, t2)^3');
      expect(fg.group_law(4).toString()).toBe('t1 + t2 + O(t1, t2)^4');
      expect(fg.group_law(7).toString()).toBe(
        't1 + t2 + t1^4*t2 + 2*t1^3*t2^2 + 2*t1^2*t2^3 + t1*t2^4 + O(t1, t2)^7'
      );
    });

    it("matches Sage's doctest for EllipticCurve([1, 2]).group_law(6)", () => {
      // sage: e.formal_group().group_law(6)
      //  t1 + t2 - 2*t1^4*t2 - 4*t1^3*t2^2 - 4*t1^2*t2^3 - 2*t1*t2^4 + O(t1, t2)^6
      const K = GF(BIGP);
      const E = EllipticCurve(K, [1n, 2n]);
      const F = E.formal_group().group_law(6);
      expect(F.coefficient(1, 0).toString()).toBe('1');
      expect(F.coefficient(0, 1).toString()).toBe('1');
      expect(F.coefficient(4, 1).toString()).toBe(reduce(-2n));
      expect(F.coefficient(3, 2).toString()).toBe(reduce(-4n));
      expect(F.coefficient(2, 3).toString()).toBe(reduce(-4n));
      expect(F.coefficient(1, 4).toString()).toBe(reduce(-2n));
    });

    it("matches Sage's doctest for curve 14a1 group_law(5)", () => {
      // sage: ehat.group_law(5)
      //  t1 + t2 - t1*t2 - 2*t1^3*t2 - 3*t1^2*t2^2 - 2*t1*t2^3 + O(t1, t2)^5
      const K = GF(BIGP);
      const E = EllipticCurve(K, [1n, 0n, 1n, 4n, -6n]);
      const F = E.formal_group().group_law(5);
      expect(F.coefficient(1, 1).toString()).toBe(reduce(-1n));
      expect(F.coefficient(3, 1).toString()).toBe(reduce(-2n));
      expect(F.coefficient(2, 2).toString()).toBe(reduce(-3n));
      expect(F.coefficient(1, 3).toString()).toBe(reduce(-2n));
    });

    it('is symmetric and satisfies F(t1, 0) = t1, F(0, t2) = t2', () => {
      const K = GF(BIGP);
      const E = EllipticCurve(K, [1n, 2n, 3n, 4n, 6n]);
      const F = E.formal_group().group_law(7);
      for (const [key, c] of F.terms) {
        const [i, j] = key.split(',').map(Number) as [number, number];
        expect(F.coefficient(j, i).toString()).toBe(c.toString());
        // F(t1, 0) = t1 means the only pure-t1 term is t1 itself.
        if (j === 0) expect(i).toBe(1);
        if (i === 0) expect(j).toBe(1);
      }
      expect(F.coefficient(1, 0).toString()).toBe('1');
      expect(F.coefficient(0, 1).toString()).toBe('1');
    });

    it('rejects a nonpositive precision', () => {
      const K = GF(7n);
      const E = EllipticCurve(K, [3n, 4n]);
      expect(() => E.formal_group().group_law(0)).toThrow('The precision must be positive.');
    });
  });

  describe('mult_by_n', () => {
    it("matches Sage's doctest for EllipticCurve(GF(17), [1, 1])", () => {
      // sage: F.mult_by_n(10, 50)
      //  10*t + 5*t^5 + 7*t^7 + 13*t^9 + t^11 + 16*t^13 + 13*t^15 + 9*t^17 + 16*t^19 + ...
      const K = GF(17n);
      const E = EllipticCurve(K, [1n, 1n]);
      const c = coeffs(E.formal_group().mult_by_n(10, 20), 20);
      const expected = new Map<number, string>([
        [1, '10'],
        [5, '5'],
        [7, '7'],
        [9, '13'],
        [11, '1'],
        [13, '16'],
        [15, '13'],
        [17, '9'],
        [19, '16'],
      ]);
      for (let i = 0; i < 20; i++) {
        expect(c[i]).toBe(expected.get(i) ?? '0');
      }
    });

    it("matches Sage's doctest for EllipticCurve(GF(101), [1, 1]).mult_by_n(100, 20)", () => {
      // sage: F.mult_by_n(100, 20) -> 100*t + O(t^20)
      const K = GF(101n);
      const E = EllipticCurve(K, [1n, 1n]);
      const c = coeffs(E.formal_group().mult_by_n(100, 20), 20);
      expect(c[1]).toBe('100');
      for (let i = 0; i < 20; i++) {
        if (i !== 1) expect(c[i]).toBe('0');
      }
    });

    it("matches Sage's doctest for EllipticCurve(QQ, [1,2,3,4,6]).mult_by_n(2, 5)", () => {
      // sage: E.formal().mult_by_n(2, prec=5) -> 2*t - t^2 - 4*t^3 - 19*t^4 + O(t^5)
      const K = GF(BIGP);
      const E = EllipticCurve(K, [1n, 2n, 3n, 4n, 6n]);
      const c = coeffs(E.formal_group().mult_by_n(2, 5), 5);
      expect(c).toEqual(['0', '2', reduce(-1n), reduce(-4n), reduce(-19n)]);
    });

    it('handles n = 0, 1 and -1', () => {
      const K = GF(BIGP);
      const E = EllipticCurve(K, [1n, 2n, 3n, 4n, 6n]);
      const fg = E.formal_group();
      expect(fg.mult_by_n(0, 5).is_zero()).toBe(true);
      expect(coeffs(fg.mult_by_n(1, 5), 5)).toEqual(['0', '1', '0', '0', '0']);
      expect(coeffs(fg.mult_by_n(-1, 5), 5)).toEqual(coeffs(fg.inverse(5), 5));
    });

    it("verifies Sage's low-degree identity [-2] = [-1] o [2] = [2] o [-1]", () => {
      const K = GF(BIGP);
      const E = EllipticCurve(K, [1n, 2n, 3n, 4n, 6n]);
      const fg = E.formal_group();
      const none = fg.mult_by_n(-1, 5);
      const two = fg.mult_by_n(2, 5);
      const ntwo = fg.mult_by_n(-2, 5);
      expect(ntwo.sub(none.__call__(two)).is_zero()).toBe(true);
      expect(ntwo.sub(two.__call__(none)).is_zero()).toBe(true);
    });
  });

  describe('sigma', () => {
    it("matches Sage's doctest for curve 14a", () => {
      // sage: E = EllipticCurve('14a'); F = E.formal_group(); F.sigma(5)
      //  t + 1/2*t^2 + 1/3*t^3 + 3/4*t^4 + O(t^5)
      const K = GF(BIGP);
      const E = EllipticCurve(K, [1n, 0n, 1n, 4n, -6n]);
      const c = coeffs(E.formal_group().sigma(5), 5);
      expect(c).toEqual(['0', '1', reduce(1n, 2n), reduce(1n, 3n), reduce(3n, 4n)]);
    });
  });
});
