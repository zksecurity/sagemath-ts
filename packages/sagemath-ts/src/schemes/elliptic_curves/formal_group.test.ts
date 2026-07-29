/**
 * Tests for the formal group of an elliptic curve.
 *
 * The expected values are SageMath's own doctests from
 * `sage/schemes/elliptic_curves/formal_group.py`.
 *
 * The first block reproduces the doctests stated over QQ inside a large prime
 * field (GF(1000003)) by reducing the rational coefficients, and runs the
 * finite-field doctests verbatim. The second block ("over QQ") runs the same
 * doctests over QQ itself, coefficient for coefficient as Sage prints them.
 * The third block exercises the characteristic-zero branch of `mult_by_n`
 * (the formal point on E base-changed to the Laurent series ring). The fourth
 * block pins the formal group with algebraic identities that hold
 * independently of any doctest, including the three-variable associativity
 * `F(x, F(y, z)) == F(F(x, y), z)` of Sage's own TESTS block.
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import { MPowerSeriesRing } from '../../rings/power_series_ring.js';
import { QQ } from '../../rings/rational_field.js';
import { EllipticCurve } from './constructor.js';
import type { BivariatePowerSeries } from './formal_group.js';

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

// ---------------------------------------------------------------------------
// Helpers for the doctests stated over QQ, which the port can run verbatim.
// ---------------------------------------------------------------------------

// `QQ` and `Rational` do not structurally satisfy the elliptic-curve
// `FieldRing`/`FieldElement` interfaces (QQ exposes `characteristic()` as a
// method, the interface wants a property), so the QQ doctests below go through
// `any`. The runtime path is exercised end to end by the assertions.

/** A rational number a/b in QQ. */
function Q(a: bigint, b: bigint = 1n): any {
  return (QQ as any).__call__(a).div((QQ as any).__call__(b));
}

/** An elliptic curve over QQ from its a-invariants. */
function curveQQ(ainvs: any[]): any {
  return EllipticCurve(QQ as any, ainvs as any) as any;
}

/** Coefficients c_0 .. c_{prec-1} of a power series, as strings. */
function qcoeffs(f: any, prec: number): string[] {
  const list = f.list();
  const out: string[] = [];
  for (let i = 0; i < prec; i++) {
    out.push(list[i] === undefined ? '0' : String(list[i]));
  }
  return out;
}

/**
 * `f(a)` for a univariate power series `f` and a multivariate argument `a`,
 * i.e. SageMath's plain function application (the coercion framework sends it
 * to `PowerSeries_poly.__call__`, power_series_poly.pyx:176).
 */
function evalAt(f: any, a: any): any {
  return f.__call__(a);
}

/** `F(A, B)` -- SageMath's `MPowerSeries._subs_formal`. */
function subs(F: BivariatePowerSeries, A: any, B: any): any {
  return (F as any)._subs_formal(A, B);
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
      // Sage's coefficients of t^-2, t^-1, t^0, ..., t^9
      const expected = [1n, 0n, 0n, -1n, 1n, 0n, -1n, 2n, -1n, -2n, 6n, -6n];
      expect(fg.x_list(10).map(String)).toEqual(expected.map((v) => reduce(v)));
      expect(String(x.__getitem__(-2))).toBe('1');
      expect(String(x.__getitem__(9))).toBe(reduce(-6n));
      expect(x.prec()).toBe(10);
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
      // Sage's coefficients of t^-3, t^-2, t^-1, t^0, ..., t^9
      const expected = [-1n, 0n, 0n, 1n, -1n, 0n, 1n, -2n, 1n, 2n, -6n, 6n, 3n];
      expect(fg.y_list(10).map(String)).toEqual(expected.map((v) => reduce(v)));
      expect(String(y.__getitem__(-3))).toBe(reduce(-1n));
      expect(String(y.__getitem__(9))).toBe('3');
      expect(y.prec()).toBe(10);
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

// ---------------------------------------------------------------------------
// The same doctests run verbatim over QQ, exactly as they appear in
// reference/sage/src/sage/schemes/elliptic_curves/formal_group.py.
// ---------------------------------------------------------------------------

describe('EllipticCurveFormalGroup over QQ (Sage doctests verbatim)', () => {
  it('w(10) for EllipticCurve([0,0,1,-1,0])', () => {
    // sage: e = EllipticCurve([0, 0, 1, -1, 0])
    // sage: e.formal_group().w(10)
    //  t^3 + t^6 - t^7 + 2*t^9 + O(t^10)
    const F = curveQQ([0n, 0n, 1n, -1n, 0n]).formal_group();
    expect(qcoeffs(F.w(10), 10)).toEqual(['0', '0', '0', '1', '0', '0', '1', '-1', '0', '2']);
  });

  it('w(7), w(20) and w(35) for EllipticCurve([3,2,-4,-2,5])', () => {
    // sage: e = EllipticCurve([3, 2, -4, -2, 5]); e.formal_group().w(20) / w(7) / w(35)
    const F = curveQQ([3n, 2n, -4n, -2n, 5n]).formal_group();

    const w20 = [
      0n,
      0n,
      0n,
      1n,
      3n,
      11n,
      35n,
      101n,
      237n,
      312n,
      -949n,
      -10389n,
      -57087n,
      -244092n,
      -865333n,
      -2455206n,
      -4366196n,
      6136610n,
      109938783n,
      688672497n,
    ];
    expect(qcoeffs(F.w(20), 20)).toEqual(w20.map(String));

    // Cached lower precision must truncate, not recompute something else.
    expect(qcoeffs(F.w(7), 7)).toEqual(w20.slice(0, 7).map(String));

    const w35 = w20.concat([
      3219525807n,
      12337076504n,
      38106669615n,
      79452618700n,
      -33430470002n,
      -1522228110356n,
      -10561222329021n,
      -52449326572178n,
      -211701726058446n,
      -693522772940043n,
      -1613471639599050n,
      -421817906421378n,
      23651687753515182n,
      181817896829144595n,
      950887648021211163n,
    ]);
    expect(qcoeffs(F.w(35), 35)).toEqual(w35.map(String));
  });

  it('x(10) and y(10) for EllipticCurve([0,0,1,-1,0])', () => {
    // sage: EllipticCurve([0, 0, 1, -1, 0]).formal_group().x(10)
    //  t^-2 - t + t^2 - t^4 + 2*t^5 - t^6 - 2*t^7 + 6*t^8 - 6*t^9 + O(t^10)
    // sage: EllipticCurve([0, 0, 1, -1, 0]).formal_group().y(10)
    //  -t^-3 + 1 - t + t^3 - 2*t^4 + t^5 + 2*t^6 - 6*t^7 + 6*t^8 + 3*t^9 + O(t^10)
    const F = curveQQ([0n, 0n, 1n, -1n, 0n]).formal_group();

    // The Laurent series print exactly as Sage prints them.
    expect(F.x(10).toString()).toBe(
      't^-2 - t + t^2 - t^4 + 2*t^5 - t^6 - 2*t^7 + 6*t^8 - 6*t^9 + O(t^10)'
    );
    expect(F.y(10).toString()).toBe(
      '-t^-3 + 1 - t + t^3 - 2*t^4 + t^5 + 2*t^6 - 6*t^7 + 6*t^8 + 3*t^9 + O(t^10)'
    );

    expect(F.x(10).valuation()).toBe(-2);
    // coefficients of t^-2, t^-1, t^0, ..., t^9
    expect(F.x_list(10).map(String)).toEqual(
      [1n, 0n, 0n, -1n, 1n, 0n, -1n, 2n, -1n, -2n, 6n, -6n].map(String)
    );

    expect(F.y(10).valuation()).toBe(-3);
    // coefficients of t^-3, t^-2, t^-1, t^0, ..., t^9
    expect(F.y_list(10).map(String)).toEqual(
      [-1n, 0n, 0n, 1n, -1n, 0n, 1n, -2n, 1n, 2n, -6n, 6n, 3n].map(String)
    );
  });

  it('differential(15) and log(15) for EllipticCurve([-1, 1/4])', () => {
    // sage: EllipticCurve([-1, 1/4]).formal_group().differential(15)
    //  1 - 2*t^4 + 3/4*t^6 + 6*t^8 - 5*t^10 - 305/16*t^12 + 105/4*t^14 + O(t^15)
    // sage: EllipticCurve([-1, 1/4]).formal_group().log(15)
    //  t - 2/5*t^5 + 3/28*t^7 + 2/3*t^9 - 5/11*t^11 - 305/208*t^13 + O(t^15)
    const F = curveQQ([Q(-1n), Q(1n, 4n)]).formal_group();

    const d = qcoeffs(F.differential(15), 15);
    expect(d).toEqual([
      '1',
      '0',
      '0',
      '0',
      '-2',
      '0',
      '3/4',
      '0',
      '6',
      '0',
      '-5',
      '0',
      '-305/16',
      '0',
      '105/4',
    ]);

    const l = qcoeffs(F.log(15), 15);
    expect(l).toEqual([
      '0',
      '1',
      '0',
      '0',
      '0',
      '-2/5',
      '0',
      '3/28',
      '0',
      '2/3',
      '0',
      '-5/11',
      '0',
      '-305/208',
      '0',
    ]);
  });

  it('inverse(6) matches the generic formula -t - a1*t^2 - a1^2*t^3 + ...', () => {
    // sage: P.<a1, a2, a3, a4, a6> = ZZ[]; E = EllipticCurve(list(P.gens()))
    // sage: E.formal_group().inverse(6)
    //  -t - a1*t^2 - a1^2*t^3 + (-a1^3 - a3)*t^4 + (-a1^4 - 3*a1*a3)*t^5 + O(t^6)
    // (the port has no elliptic curve over a multivariate polynomial ring, so
    //  the generic formula is instantiated at a1..a6 = 1,2,3,4,6 over QQ)
    const [a1, a3] = [1n, 3n];
    const F = curveQQ([a1, 2n, a3, 4n, 6n]).formal_group();
    expect(qcoeffs(F.inverse(6), 6)).toEqual([
      '0',
      '-1',
      String(-a1),
      String(-(a1 * a1)),
      String(-(a1 ** 3n + a3)),
      String(-(a1 ** 4n + 3n * a1 * a3)),
    ]);
  });

  it('group_law(6) for EllipticCurve([1, 2])', () => {
    // sage: e = EllipticCurve([1, 2]); e.formal_group().group_law(6)
    //  t1 + t2 - 2*t1^4*t2 - 4*t1^3*t2^2 - 4*t1^2*t2^3 - 2*t1*t2^4 + O(t1, t2)^6
    const F = curveQQ([1n, 2n]).formal_group().group_law(6);
    const expected = new Map<string, string>([
      ['1,0', '1'],
      ['0,1', '1'],
      ['4,1', '-2'],
      ['3,2', '-4'],
      ['2,3', '-4'],
      ['1,4', '-2'],
    ]);
    expect(new Set(F.terms.keys())).toEqual(new Set(expected.keys()));
    for (const [key, v] of expected) {
      const [i, j] = key.split(',').map(Number) as [number, number];
      expect(String(F.coefficient(i, j))).toBe(v);
    }
  });

  it('group_law(3) and group_law(5) for curve 14a1', () => {
    // sage: e = EllipticCurve('14a1'); ehat = e.formal()
    // sage: ehat.group_law(3)  ->  t1 + t2 - t1*t2 + O(t1, t2)^3
    // sage: ehat.group_law(5)
    //  t1 + t2 - t1*t2 - 2*t1^3*t2 - 3*t1^2*t2^2 - 2*t1*t2^3 + O(t1, t2)^5
    const fg = curveQQ([1n, 0n, 1n, 4n, -6n]).formal_group();

    const F3 = fg.group_law(3);
    expect(new Set(F3.terms.keys())).toEqual(new Set(['1,0', '0,1', '1,1']));
    expect(String(F3.coefficient(1, 1))).toBe('-1');

    const F5 = fg.group_law(5);
    const expected = new Map<string, string>([
      ['1,0', '1'],
      ['0,1', '1'],
      ['1,1', '-1'],
      ['3,1', '-2'],
      ['2,2', '-3'],
      ['1,3', '-2'],
    ]);
    expect(new Set(F5.terms.keys())).toEqual(new Set(expected.keys()));
    for (const [key, v] of expected) {
      const [i, j] = key.split(',').map(Number) as [number, number];
      expect(String(F5.coefficient(i, j))).toBe(v);
    }
  });

  it('mult_by_n(2, 5) for EllipticCurve(QQ, [1,2,3,4,6])', () => {
    // sage: E = EllipticCurve(QQ, [1,2,3,4,6]); E.formal().mult_by_n(2, prec=5)
    //  2*t - t^2 - 4*t^3 - 19*t^4 + O(t^5)
    const F = curveQQ([1n, 2n, 3n, 4n, 6n]).formal_group();
    expect(qcoeffs(F.mult_by_n(2, 5), 5)).toEqual(['0', '2', '-1', '-4', '-19']);
  });

  it('mult_by_n(100, 20) for curve 37a', () => {
    // sage: E = EllipticCurve("37a"); F = E.formal_group(); F.mult_by_n(100, 20)
    const F = curveQQ([0n, 0n, 1n, -1n, 0n]).formal_group();
    const m = F.mult_by_n(100, 20);

    // Character for character, the line printed by Sage's doctest.
    expect(String(m)).toBe(
      '100*t - 49999950*t^4 + 3999999960*t^5 + 14285614285800*t^7 - 2999989920000150*t^8 + ' +
        '133333325333333400*t^9 - 3571378571674999800*t^10 + 1402585362624965454000*t^11 - ' +
        '146666057066712847999500*t^12 + 5336978000014213190385000*t^13 - ' +
        '519472790950932256570002000*t^14 + 93851927683683567270392002800*t^15 - ' +
        '6673787211563812368630730325175*t^16 + 320129060335050875009191524993000*t^17 - ' +
        '45670288869783478472872833214986000*t^18 + ' +
        '5302464956134111125466184947310391600*t^19 + O(t^20)'
    );

    const c = qcoeffs(m, 20);
    const expected = new Map<number, string>([
      [1, '100'],
      [4, '-49999950'],
      [5, '3999999960'],
      [7, '14285614285800'],
      [8, '-2999989920000150'],
      [9, '133333325333333400'],
      [10, '-3571378571674999800'],
      [11, '1402585362624965454000'],
      [12, '-146666057066712847999500'],
      [13, '5336978000014213190385000'],
      [14, '-519472790950932256570002000'],
      [15, '93851927683683567270392002800'],
      [16, '-6673787211563812368630730325175'],
      [17, '320129060335050875009191524993000'],
      [18, '-45670288869783478472872833214986000'],
      [19, '5302464956134111125466184947310391600'],
    ]);
    for (let i = 0; i < 20; i++) {
      expect(c[i]).toBe(expected.get(i) ?? '0');
    }
  });

  it('mult_by_n identity [-2] = [-1] o [2] = [2] o [-1] for e = [1,2,3,4,6]', () => {
    // sage: none = e.formal_group().mult_by_n(-1, 5); two = ...mult_by_n(2, 5)
    // sage: ntwo = ...mult_by_n(-2, 5); ntwo - none(two)  ->  O(t^5)
    const F = curveQQ([1n, 2n, 3n, 4n, 6n]).formal_group();
    const none = F.mult_by_n(-1, 5);
    const two = F.mult_by_n(2, 5);
    const ntwo = F.mult_by_n(-2, 5);
    expect(ntwo.sub(none.__call__(two)).is_zero()).toBe(true);
    expect(ntwo.sub(two.__call__(none)).is_zero()).toBe(true);
  });

  it('sigma(5) for curve 14a', () => {
    // sage: E = EllipticCurve('14a'); E.formal_group().sigma(5)
    //  t + 1/2*t^2 + 1/3*t^3 + 3/4*t^4 + O(t^5)
    const F = curveQQ([1n, 0n, 1n, 4n, -6n]).formal_group();
    expect(qcoeffs(F.sigma(5), 5)).toEqual(['0', '1', '1/2', '1/3', '3/4']);
  });

  it('mult_by_n(10, 50) for EllipticCurve(GF(17), [1, 1])', () => {
    // sage: F = EllipticCurve(GF(17), [1, 1]).formal_group(); F.mult_by_n(10, 50)
    const F = EllipticCurve(GF(17n), [1n, 1n]).formal_group();
    const c = coeffs(F.mult_by_n(10, 50), 50);
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
      [23, '15'],
      [25, '15'],
      [27, '2'],
      [29, '10'],
      [31, '8'],
      [33, '15'],
      [35, '6'],
      [37, '7'],
      [39, '9'],
      [41, '10'],
      [43, '5'],
      [45, '4'],
      [47, '6'],
      [49, '13'],
    ]);
    for (let i = 0; i < 50; i++) {
      expect(c[i]).toBe(expected.get(i) ?? '0');
    }
  }, 30000);

  it('mult_by_n(0, 5) and mult_by_n(1, 5) for e = EllipticCurve([1,2,3,4,6])', () => {
    // sage: e = EllipticCurve([1, 2, 3, 4, 6])
    // sage: e.formal_group().mult_by_n(0, 5)  ->  O(t^5)
    // sage: e.formal_group().mult_by_n(1, 5)  ->  t + O(t^5)
    const F = curveQQ([1n, 2n, 3n, 4n, 6n]).formal_group();
    expect(String(F.mult_by_n(0, 5))).toBe('O(t^5)');
    expect(String(F.mult_by_n(1, 5))).toBe('t + O(t^5)');
  });
});

// ---------------------------------------------------------------------------
// The characteristic-zero branch of ``mult_by_n``
// (formal_group.py:644-665): a formal point on E base-changed to the Laurent
// series ring, multiplied by n with the curve's own group law.
// ---------------------------------------------------------------------------

describe('mult_by_n over a field of characteristic zero', () => {
  it('takes the fast branch over QQ and the group-law branch over GF(17)', () => {
    // The two branches are observably different at low precision.  The general
    // branch answers at exactly the requested precision (it composes the group
    // law, which carries O(t1,t2)^prec).  The characteristic-zero branch
    // instead builds the formal point from ``x(prec-3)`` and ``y(prec-4)``;
    // with prec = 2 those are ``x(0)`` and ``y(0)`` (Sage's ``x``/``y`` clamp a
    // negative precision to 0), of relative precision 2 and 3, so -Q[0]/Q[1]
    // comes out with valuation 1 and relative precision 2, i.e. O(t^3).
    const q = curveQQ([1n, 2n, 3n, 4n, 6n]).formal_group().mult_by_n(2, 2);
    expect(q.prec()).toBe(3);

    const f = EllipticCurve(GF(17n), [1n, 1n]).formal_group().mult_by_n(2, 2);
    expect(f.prec()).toBe(2);

    // Whatever the branch, the coefficients are those of the doctest-verified
    // series 2*t - t^2 - 4*t^3 - 19*t^4 + O(t^5).
    const full = curveQQ([1n, 2n, 3n, 4n, 6n]).formal_group().mult_by_n(2, 5);
    expect(q.sub(full.add_bigoh(3)).is_zero()).toBe(true);
    expect(String(q)).toBe('2*t - t^2 + O(t^3)');
  });

  it('agrees with the group-law double-and-add it replaces (QQ)', () => {
    // Independent oracle: iterate the (doctest-verified) formal group law,
    // [n+1](t) = F([n](t), t), which is exactly what upstream's *general*
    // branch does -- and which the characteristic-zero branch never touches.
    const F = curveQQ([1n, 2n, 3n, 4n, 6n]).formal_group();
    const prec = 9;
    const GL = F.group_law(prec);
    const t = F.log(prec).parent().gen().add_bigoh(prec);

    let byGroupLaw = t;
    for (let n = 1; n <= 8; n++) {
      if (n > 1) byGroupLaw = subs(GL, byGroupLaw, t);
      const byFastPath = F.mult_by_n(n, prec);
      expect(byFastPath.sub(byGroupLaw).add_bigoh(prec).is_zero()).toBe(true);
      // and the negative multiple is the formal inverse of the positive one
      const neg = F.mult_by_n(-n, prec);
      expect(subs(GL, byFastPath, neg).add_bigoh(prec).is_zero()).toBe(true);
    }
    // Not vacuous: [8](t) = 8*t + ... is a genuine series, not O(t^9).
    expect(String(F.mult_by_n(8, prec).list()[1])).toBe('8');
  });

  it('[m] o [n] = [m*n] and [n](t) = n*t + O(t^2) over QQ', () => {
    const F = curveQQ([0n, 0n, 1n, -1n, 0n]).formal_group(); // 37a
    const prec = 10;
    for (const [m, n] of [
      [2, 3],
      [3, 4],
      [-2, 5],
      [7, -1],
      [6, 6],
    ] as [number, number][]) {
      const lhs = F.mult_by_n(m, prec).__call__(F.mult_by_n(n, prec)).add_bigoh(prec);
      const rhs = F.mult_by_n(m * n, prec).add_bigoh(prec);
      expect(lhs.sub(rhs).is_zero()).toBe(true);
      expect(String(rhs.list()[1])).toBe(String(m * n));
    }
  });

  it('the formal point it multiplies lies on the curve, and -x/y = t', () => {
    // The base change E -> E x_QQ QQ((t)) and the point (x(t), y(t)) that
    // upstream's branch forms: the Weierstrass equation holds as Laurent
    // series, and the formal parameter is recovered as -x/y.
    const E = curveQQ([1n, 2n, 3n, 4n, 6n]);
    const [a1, a2, a3, a4, a6] = E.ainvs();
    const F = E.formal_group();
    const prec = 15;
    const x = F.x(prec);
    const y = F.y(prec);
    const K = x.parent();
    const lhs = y.mul(y).add(x.mul(y).scalar_mul(a1)).add(y.scalar_mul(a3));
    const rhs = x
      .mul(x)
      .mul(x)
      .add(x.mul(x).scalar_mul(a2))
      .add(x.scalar_mul(a4))
      .add(K.__call__(a6));
    expect(lhs.sub(rhs).is_zero()).toBe(true);
    expect(String(x.neg().div(y).add_bigoh(prec))).toBe('t + O(t^15)');
    // and with the precisions upstream's branch actually uses, -x/y is [1](t)
    expect(
      String(
        F.x(prec - 3)
          .neg()
          .div(F.y(prec - 4))
      )
    ).toBe(String(F.mult_by_n(1, prec)));
  });
});

// ---------------------------------------------------------------------------
// Algebraic identities that pin the formal group independently of the doctests.
// ---------------------------------------------------------------------------

describe('EllipticCurveFormalGroup identities', () => {
  it('x(t) and y(t) satisfy the Weierstrass equation as Laurent series', () => {
    const E = curveQQ([1n, 2n, 3n, 4n, 6n]);
    const F = E.formal_group();
    const [a1, a2, a3, a4, a6] = E.ainvs();
    const prec = 12;

    // Multiply y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6 through by t^6 and
    // work with X = t^2*x(t) and Y = t^3*y(t), which are honest power series.
    const R = F.w(1).parent();
    const xs = F.x_list(prec); // c_{-2}, c_{-1}, ...
    const ys = F.y_list(prec); // c_{-3}, c_{-2}, ...
    const X = R.__call__(xs, xs.length); // = t^2 * x(t)
    const Y = R.__call__(ys, ys.length); // = t^3 * y(t)
    const n = Math.min(xs.length, ys.length);

    const lhs = Y.mul(Y)
      .add(X.mul(Y)._scalarMul(a1)._shiftLeft(1))
      .add(Y._scalarMul(a3)._shiftLeft(3));
    const rhs = X.mul(X)
      .mul(X)
      .add(X.mul(X)._scalarMul(a2)._shiftLeft(2))
      .add(X._scalarMul(a4)._shiftLeft(4))
      .add(R.__call__([a6], Number.POSITIVE_INFINITY)._shiftLeft(6));

    expect(lhs.sub(rhs).add_bigoh(n).is_zero()).toBe(true);
  });

  it('differential starts 1 + a1*t + (a1^2 + a2)*t^2 (Silverman p. 113)', () => {
    const E = curveQQ([1n, 2n, 3n, 4n, 6n]);
    const [a1, a2] = E.ainvs();
    const d = E.formal_group().differential(5).list();
    expect(String(d[0])).toBe('1');
    expect(String(d[1])).toBe(String(a1));
    expect(String(d[2])).toBe(String(a1.mul(a1).add(a2)));
  });

  it('log and its reversion (the formal exponential) are mutually inverse', () => {
    const F = curveQQ([1n, 2n, 3n, 4n, 6n]).formal_group();
    const prec = 12;
    const log = F.log(prec);
    const exp = log.reversion(prec);
    const t = log.parent().gen();
    expect(log.__call__(exp).sub(t).add_bigoh(prec).is_zero()).toBe(true);
    expect(exp.__call__(log).sub(t).add_bigoh(prec).is_zero()).toBe(true);
  });

  it('log is a homomorphism: log(F(t1,t2)) = log(t1) + log(t2)', () => {
    const E = curveQQ([1n, 2n, 3n, 4n, 6n]);
    const F = E.formal_group();
    const prec = 8;
    const GL = F.group_law(prec);
    const [T1, T2] = GL.parent().gens();
    const log = F.log(prec);
    const lhs = evalAt(log, GL);
    const rhs = evalAt(log, T1).add(evalAt(log, T2));
    expect(lhs.sub(rhs).add_bigoh(prec).is_zero()).toBe(true);
    // Not vacuous: log(F) really is t1 + t2 + ... to precision 8.
    expect(lhs.prec()).toBeGreaterThanOrEqual(prec);
    expect(String(lhs.__getitem__([1, 0]))).toBe('1');
    expect(String(lhs.__getitem__([0, 1]))).toBe('1');
  });

  it('log([n](t)) = n*log(t) for several n', () => {
    const F = curveQQ([1n, 2n, 3n, 4n, 6n]).formal_group();
    const prec = 8;
    const log = F.log(prec);
    for (const n of [2, 3, 5, -3]) {
      const lhs = log.__call__(F.mult_by_n(n, prec)).add_bigoh(prec);
      const rhs = log._scalarMul(Q(BigInt(n))).add_bigoh(prec);
      expect(lhs.sub(rhs).is_zero()).toBe(true);
    }
  });

  it("Sage's TESTS block: F(x, i(x)) = 0, F(x,y) = F(y,x), F(x,F(y,z)) = F(F(x,y),z)", () => {
    // sage: e = EllipticCurve(GF(7), [3, 4]); ehat = e.formal()
    // sage: F = ehat.group_law(7); F
    // t1 + t2 + t1^4*t2 + 2*t1^3*t2^2 + 2*t1^2*t2^3 + t1*t2^4 + O(t1, t2)^7
    // sage: R.<x,y,z> = GF(7)[[]]
    // sage: F(x, ehat.inverse()(x))
    // 0 + O(x, y, z)^7
    // sage: F(x, y) == F(y, x)
    // True
    // sage: F(x, F(y, z)) == F(F(x, y), z)
    // True
    const K = GF(7n);
    const ehat = EllipticCurve(K, [3n, 4n]).formal_group();
    const prec = 7;
    const F = ehat.group_law(prec);
    expect(F.toString()).toBe(
      't1 + t2 + t1^4*t2 + 2*t1^3*t2^2 + 2*t1^2*t2^3 + t1*t2^4 + O(t1, t2)^7'
    );

    // R.<x,y,z> = GF(7)[[]]
    const R3 = new MPowerSeriesRing(K as any, 'x,y,z');
    const [x, y, z] = R3.gens();

    // F(x, ehat.inverse()(x)) -> 0 + O(x, y, z)^7
    expect(subs(F, x, evalAt(ehat.inverse(prec), x)).toString()).toBe('0 + O(x, y, z)^7');

    // F(x, y) == F(y, x)
    expect(subs(F, x, y).eq(subs(F, y, x))).toBe(true);

    // F(x, F(y, z)) == F(F(x, y), z) -- the genuine three-variable identity
    const lhs = subs(F, x, subs(F, y, z));
    const rhs = subs(F, subs(F, x, y), z);
    expect(lhs.eq(rhs)).toBe(true);

    // The identity is not vacuous: both sides are x + y + z + ... to O(...)^7,
    // and they really do involve all three variables.
    expect(lhs.prec()).toBe(prec);
    expect(String(lhs.__getitem__([1, 0, 0]))).toBe('1');
    expect(String(lhs.__getitem__([0, 1, 0]))).toBe('1');
    expect(String(lhs.__getitem__([0, 0, 1]))).toBe('1');
    expect(lhs.is_zero()).toBe(false);
    expect(
      lhs.monomial_coefficients().some(([e]: [number[], unknown]) => e[0]! > 0 && e[2]! > 0)
    ).toBe(true);
  });

  it('the group law is commutative and associative in three variables (GF(17), QQ)', () => {
    // Same three-variable identities as Sage's GF(7) TESTS block, on two more
    // curves: one in characteristic 17, one over QQ with a1, a2, a3, a4, a6 all
    // nonzero (so every term of the group-law formula contributes).
    for (const [E, prec] of [
      [EllipticCurve(GF(17n), [1n, 1n]) as any, 10],
      [curveQQ([1n, 2n, 3n, 4n, 6n]), 8],
    ] as [any, number][]) {
      const k = E.base_ring;
      const F = E.formal_group().group_law(prec);

      // F(t1, t2) == F(t2, t1) coefficientwise
      for (const [e, c] of F.monomial_coefficients()) {
        const [i, j] = e as [number, number];
        expect(String(F.coefficient(j, i))).toBe(String(c));
      }

      const R3 = new MPowerSeriesRing(k, 'x,y,z');
      const [x, y, z] = R3.gens();
      expect(subs(F, x, y).eq(subs(F, y, x))).toBe(true);
      const lhs = subs(F, x, subs(F, y, z));
      const rhs = subs(F, subs(F, x, y), z);
      expect(lhs.eq(rhs)).toBe(true);
      expect(lhs.prec()).toBe(prec);
      expect(String(lhs.__getitem__([1, 0, 0]))).toBe('1');
      expect(String(lhs.__getitem__([0, 0, 1]))).toBe('1');
      // F(x, i(x)) = 0
      expect(subs(F, x, evalAt(E.formal_group().inverse(prec), x)).is_zero()).toBe(true);
    }
  });

  it('[m+n](t) = F([m](t), [n](t)) over GF(17)', () => {
    const E = EllipticCurve(GF(17n), [1n, 1n]);
    const F = E.formal_group();
    const prec = 12;
    const GL = F.group_law(prec);
    for (const [a, b] of [
      [2, 3],
      [4, 5],
      [7, 1],
      [3, -3],
    ] as [number, number][]) {
      const lhs = F.mult_by_n(a + b, prec);
      const rhs = subs(GL, F.mult_by_n(a, prec), F.mult_by_n(b, prec));
      expect(lhs.sub(rhs).add_bigoh(prec).is_zero()).toBe(true);
    }
  }, 30000);

  it('log over GF(7) is exact below t^7 and refuses to divide by 7', () => {
    // In characteristic p the formal logarithm only makes sense before t^p:
    // integrating the t^{p-1} coefficient divides by p. Sage raises there too.
    const F = EllipticCurve(GF(7n), [3n, 4n]).formal_group();
    // differential is 1 + 6*t^4 + 5*t^6 + ... -- not a zero-padded placeholder
    const d = coeffs(F.differential(10), 9);
    expect(d[0]).toBe('1');
    expect(d[4]).toBe('6');
    expect(d[6]).toBe('5');
    // log(7) integrates only t^0 .. t^5, so no division by 7 occurs
    expect(coeffs(F.log(7), 6)).toEqual(['0', '1', '0', '0', '0', '4']);
    // log(10) would need to divide the t^6 coefficient by 7
    // Verified: SageMath 10.3 raises
    // `ZeroDivisionError: inverse of Mod(0, 7) does not exist` here.
    // (This assertion previously pinned the port's own `division by zero in GF(p)`.)
    expect(() => F.log(10)).toThrow('inverse of Mod(0, 7) does not exist');
  });

  it("reproduces Sage's Bernardi sigma function for curve 14a", () => {
    // sage: E = EllipticCurve('14a'); L = E.padic_lseries(5)
    // sage: L.bernardi_sigma_function(prec=5)
    //  z + 1/24*z^3 + 29/384*z^5 - 8399/322560*z^7 - 291743/92897280*z^9 + O(z^10)
    //
    // padic_lseries.py:bernardi_sigma_function is formal_group.sigma with the
    // (a1^2 + 4*a2)/12 shift dropped and expressed in z = log(t); dropping that
    // shift multiplies sigma by exp(r*z^2/2) with r = (a1^2 + 4*a2)/12.
    const E = curveQQ([1n, 0n, 1n, 4n, -6n]);
    const [a1, a2] = E.ainvs();
    const F = E.formal_group();
    const N = 10;
    const sigmaOfT = F.sigma(N);
    const log = F.log(N);
    const sigmaOfZ = sigmaOfT.__call__(log.reversion(N)).add_bigoh(N);
    const r = a1
      .mul(a1)
      .add(a2.mul(Q(4n)))
      .div(Q(12n));
    const R = sigmaOfT.parent();
    const quad = R.gen()
      .pow(2)
      ._scalarMul(r.div(Q(2n)))
      .add_bigoh(N);
    const bernardi = sigmaOfZ.mul(quad.exp()).add_bigoh(N);
    expect(qcoeffs(bernardi, N)).toEqual([
      '0',
      '1',
      '0',
      '1/24',
      '0',
      '29/384',
      '0',
      '-8399/322560',
      '0',
      '-291743/92897280',
    ]);
  });
});
