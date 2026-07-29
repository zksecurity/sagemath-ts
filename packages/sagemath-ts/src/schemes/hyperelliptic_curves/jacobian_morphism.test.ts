/**
 * Tests for Cantor's algorithm and the Mumford representation.
 *
 * Every expected value below was produced by running the real SageMath
 * (`sage`, version 10.3) on the corresponding input; the docstring examples in
 * `sage/schemes/hyperelliptic_curves/jacobian_morphism.py` and
 * `jacobian_homset.py` are reproduced verbatim.
 *
 * In addition the group axioms (associativity, identity, inverse) and the
 * group order are checked against an exhaustive enumeration of all reduced
 * Mumford divisors over small finite fields.
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import type { FiniteFieldElement } from '../../rings/finite_rings/finite_field_prime.js';
import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import { Rational } from '../../rings/rational.js';
import { QQ } from '../../rings/rational_field.js';
import { HyperellipticCurve } from './constructor.js';
import type { HyperellipticCurve_finite_field } from './hyperelliptic_finite_field.js';
import type { HyperellipticJacobian_generic } from './jacobian_generic.js';
import type { JacobianHomset_divisor_classes } from './jacobian_homset.js';
import type { JacobianMorphism_divisor_class_field } from './jacobian_morphism.js';

type FF = FiniteFieldElement;

/** Build `GF(p)[x]` together with a few helpers. */
function fp(p: bigint, name = 'x') {
  const K = GF(p);
  const R = new PolynomialRing<FF>(K, name);
  return {
    K,
    R,
    x: R.gen(),
    /** constant polynomial */
    c: (n: bigint) => R.__call__(K.__call__(n)),
    /** field element */
    e: (n: bigint) => K.__call__(n),
  };
}

function qq(name = 'x') {
  const R = new PolynomialRing<Rational>(QQ as never, name);
  return {
    R,
    x: R.gen(),
    c: (n: bigint) => R.__call__(new Rational(n)),
    e: (n: bigint) => new Rational(n),
  };
}

// biome-ignore lint/suspicious/noExplicitAny: the constructor returns the dynamic subclass
type AnyCurve = any;
// The repo's `FiniteFieldElement` does not structurally satisfy `RingElement`
// in explicit type positions, so the generic classes of this module are
// referred to through these aliases in the tests.
// biome-ignore lint/suspicious/noExplicitAny: see above
type AnyJacobian = any;
// biome-ignore lint/suspicious/noExplicitAny: see above
type AnyHomset = any;
// biome-ignore lint/suspicious/noExplicitAny: see above
type AnyDivisor = any;

describe('jacobian_morphism module docstring (GF(37))', () => {
  const { R, x, c } = fp(37n);
  const f = x
    .pow(5)
    .add(x.pow(4).scalar_mul(R.base_ring.__call__(12n)))
    .add(x.pow(3).scalar_mul(R.base_ring.__call__(13n)))
    .add(x.pow(2).scalar_mul(R.base_ring.__call__(15n)))
    .add(x.scalar_mul(R.base_ring.__call__(33n)));
  const H: AnyCurve = HyperellipticCurve(f);
  const J = H.jacobian() as AnyJacobian;
  const JK = J.point_homset();
  void c;

  it('prints the curve and its Jacobian like Sage', () => {
    expect(String(H)).toBe(
      'Hyperelliptic Curve over Finite Field of size 37 defined by ' +
        'y^2 = x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x'
    );
    expect(String(J)).toBe(
      'Jacobian of Hyperelliptic Curve over Finite Field of size 37 defined by ' +
        'y^2 = x^5 + 12*x^4 + 13*x^3 + 15*x^2 + 33*x'
    );
  });

  it('reproduces the module docstring computation', () => {
    const P1 = H.lift_x(2n);
    const Q1 = H.lift_x(10n);
    expect(String(P1)).toBe('(2 : 11 : 1)');
    expect(String(Q1)).toBe('(10 : 18 : 1)');

    const P = JK.__call__(P1);
    const Q = JK.__call__(Q1);
    expect(String(P)).toBe('(x + 35, y + 26)');
    expect(String(Q)).toBe('(x + 27, y + 19)');
    expect(String(P.add(Q))).toBe('(x^2 + 25*x + 20, y + 13*x)');

    // 1904*P == 0, 34*P == 0, 35*P == P, 33*P == -P
    expect(String(P.mul(1904n))).toBe('(1)');
    expect(P.mul(34n).is_zero()).toBe(true);
    expect(P.mul(35n).eq(P)).toBe(true);
    expect(P.mul(33n).eq(P.neg())).toBe(true);

    expect(String(Q.mul(1904n))).toBe('(1)');
    expect(Q.mul(238n).is_zero()).toBe(true);
    expect(Q.mul(239n).eq(Q)).toBe(true);
    expect(Q.mul(237n).eq(Q.neg())).toBe(true);
  });

  it('reproduces the __neg__, _add_ and _sub_ docstrings', () => {
    const P1 = JK.__call__(H.lift_x(2n));
    expect(String(P1)).toBe('(x + 35, y + 26)');
    expect(String(P1.neg())).toBe('(x + 35, y + 11)');
    expect(String(P1.add(P1.neg()))).toBe('(1)');
    expect(String(P1.add(P1))).toBe('(x^2 + 33*x + 4, y + 13*x)');

    const P2 = JK.__call__(H.lift_x(4n));
    expect(String(P2)).toBe('(x + 33, y + 34)');
    expect(String(P1.sub(P2))).toBe('(x^2 + 31*x + 8, y + 7*x + 12)');
    expect(String(P1.add(P2))).toBe('(x^2 + 31*x + 8, y + 4*x + 18)');
    expect(String(P1.sub(P2).sub(P1.add(P2)).add(P2.mul(2n)))).toBe('(1)');
    expect(String(P1.sub(P1))).toBe('(1)');
  });

  it('has the second curve of the __neg__ docstring (h nonzero)', () => {
    const H2: AnyCurve = HyperellipticCurve(f, x);
    const J2 = (H2.jacobian() as AnyJacobian).point_homset();
    const P2 = J2.__call__(H2.lift_x(2n));
    expect(String(P2)).toBe('(x + 35, y + 24)');
    expect(String(P2.neg())).toBe('(x + 35, y + 15)');
    expect(String(P2.add(P2.neg()))).toBe('(1)');
  });
});

describe('cantor_reduction_simple / cantor_reduction docstrings (QQ)', () => {
  it('cantor_reduction_simple: 2-torsion on y^2 = x^5 - x', () => {
    const { x } = qq();
    const H: AnyCurve = HyperellipticCurve(x.pow(5).sub(x));
    expect(String(H)).toBe('Hyperelliptic Curve over Rational Field defined by y^2 = x^5 - x');
    const J = (H.jacobian() as AnyJacobian).point_homset();
    const P = J.__call__(H.lift_x(new Rational(-1n)));
    expect(String(P)).toBe('(x + 1, y)');
    expect(String(P.mul(2n))).toBe('(1)');
  });

  it('cantor_reduction: y^2 + x*y = x^5 - x', () => {
    const { x } = qq();
    const H: AnyCurve = HyperellipticCurve(x.pow(5).sub(x), x);
    expect(String(H)).toBe(
      'Hyperelliptic Curve over Rational Field defined by y^2 + x*y = x^5 - x'
    );
    const J = (H.jacobian() as AnyJacobian).point_homset();
    const Q = J.__call__(H.lift_x(new Rational(0n)));
    expect(String(Q)).toBe('(x, y)');
    expect(String(Q.mul(2n))).toBe('(1)');

    const P = J.__call__(H.lift_x(new Rational(-1n)));
    expect(String(P)).toBe('(x + 1, y)');
    expect(String(P.mul(2n))).toBe('(x^2 + 2*x + 1, y + 4*x + 4)');
    expect(String(P.mul(3n))).toBe('(x^2 - 487*x - 324, y + 10755*x + 7146)');
  });

  it('jacobian_generic docstring: y^2 + u*v = u^5 - u + 1 over QQ', () => {
    const { R, x } = qq('u');
    const H: AnyCurve = HyperellipticCurve(x.pow(5).sub(x).add(R.one()), x, { names: 'u,v' });
    expect(String(H)).toBe(
      'Hyperelliptic Curve over Rational Field defined by v^2 + u*v = u^5 - u + 1'
    );
    const J = H.jacobian() as AnyJacobian;
    const P = H.point([
      QQ.__call__(0n) as Rational,
      QQ.__call__(1n) as Rational,
      QQ.__call__(1n) as Rational,
    ]);
    expect(String(P)).toBe('(0 : 1 : 1)');
    const Q = J.point_homset().__call__(P);
    const expected = ['(1)', '(u, v - 1)', '(u^2, v + u - 1)', '(u^2, v + 1)', '(u, v + 1)', '(1)'];
    for (let i = 0; i < 6; i++) {
      expect(String(Q.mul(BigInt(i)))).toBe(expected[i]);
    }
  });

  it('jacobian_homset docstring: y^2 = x^5 + x + 1 over QQ', () => {
    const { x, c } = qq();
    const H: AnyCurve = HyperellipticCurve(x.pow(5).add(x).add(c(1n)));
    const J = (H.jacobian() as AnyJacobian).point_homset();
    const P = H.point([
      QQ.__call__(0n) as Rational,
      QQ.__call__(1n) as Rational,
      QQ.__call__(1n) as Rational,
    ]);
    const Q = J.__call__(P);
    expect(String(Q)).toBe('(x, y - 1)');
    expect(String(Q.add(Q))).toBe('(x^2, y - 1/2*x - 1)');
    expect(String(Q.mul(3n))).toBe('(x^2 - 1/64*x + 1/8, y + 255/512*x + 65/64)');
  });

  it('jacobian_homset docstring: y^2 = x^5 - 1 over GF(3)', () => {
    const { R, x, c, e } = fp(3n);
    const H: AnyCurve = HyperellipticCurve(x.pow(5).sub(c(1n)));
    const X = (H.jacobian() as AnyJacobian).point_homset();
    const a = x.pow(2).sub(x).add(c(1n));
    const b = x.neg().add(c(1n));
    const cc = x.sub(c(1n));
    const d = R.zero();
    void e;
    const D1 = X.__call__([a, b]);
    const D2 = X.__call__([cc, d]);
    expect(String(D1)).toBe('(x^2 + 2*x + 1, y + x + 2)');
    expect(String(D2)).toBe('(x + 2, y)');
    expect(String(D1.add(D2))).toBe('(x^2 + 2*x + 2, y + 2*x + 1)');
  });

  it('jacobian_generic docstring: J(FF) over GF(2003)', () => {
    const { R, x, c } = fp(2003n);
    const f = x
      .pow(5)
      .add(x.pow(3).scalar_mul(R.base_ring.__call__(1184n)))
      .add(x.pow(2).scalar_mul(R.base_ring.__call__(1846n)))
      .add(x.scalar_mul(R.base_ring.__call__(956n)))
      .add(c(560n));
    const H: AnyCurve = HyperellipticCurve(f);
    const J = H.jacobian() as AnyJacobian;
    const a = x
      .pow(2)
      .add(x.scalar_mul(R.base_ring.__call__(376n)))
      .add(c(245n));
    const b = x.scalar_mul(R.base_ring.__call__(1015n)).add(c(1368n));
    const X = J.point_homset();
    const D = X.__call__([a, b]);
    expect(String(D)).toBe('(x^2 + 376*x + 245, y + 988*x + 635)');
    expect(String(J.point(0n))).toBe('(1)');
    expect(D.eq(J.point([a, b]) as AnyDivisor)).toBe(true);
    expect(D.eq(D.add(J.point(0n) as AnyDivisor))).toBe(true);
  });
});

describe('an even degree model exercises the point-at-infinity branch', () => {
  // f, h found by instrumenting Sage's cantor_reduction; the sequences below
  // are exactly what Sage prints for i*D, i = 0..12.
  const cases: Array<{ f: bigint[]; h: bigint[]; seq: string[] }> = [
    {
      f: [3n, 6n, 9n, 4n, 7n, 9n, 6n],
      h: [10n, 2n, 9n, 10n],
      seq: [
        '(1)',
        '(x + 10, y)',
        '(x^2 + 9*x + 1, y + x + 10)',
        '(x^2 + 7*x + 7, y + 10)',
        '(x^2 + 3*x + 4, y + 7*x + 7)',
        '(x^2 + 8, y + 3*x + 8)',
        '(x^2 + 4*x + 7, y + 6*x + 3)',
        '(x^2 + 10*x + 8, y + 4*x + 1)',
        '(x^2 + 4*x + 2, y + 9*x + 7)',
        '(x^2 + 2*x + 7, y + 3*x + 8)',
        '(x + 7, y + 5)',
        '(x + 8, y + 9)',
        '(x^2 + 10*x + 2, y + 6*x + 4)',
      ],
    },
    {
      f: [5n, 9n, 10n, 3n, 6n, 8n, 3n],
      h: [7n, 9n, 10n, 2n],
      seq: [
        '(1)',
        '(x, y + 1)',
        '(x^2, y + 3*x + 1)',
        '(x^2 + 6*x + 4, y + 7*x + 4)',
        '(x^2 + 8*x + 1, y + 5*x + 5)',
        '(x^2 + 6*x + 4, y + 6*x)',
        '(x^2 + 7*x + 3, y + 3*x + 8)',
        '(x, y + 6)',
        '(1)',
        '(x, y + 1)',
        '(x^2, y + 3*x + 1)',
        '(x^2 + 6*x + 4, y + 7*x + 4)',
        '(x^2 + 8*x + 1, y + 5*x + 5)',
      ],
    },
  ];

  for (const [idx, tc] of cases.entries()) {
    it(`matches SageMath on the multiples of a point (case ${idx})`, () => {
      const { R, K } = fp(11n);
      const f = R.__call__(tc.f.map((v) => K.__call__(v)));
      const h = R.__call__(tc.h.map((v) => K.__call__(v)));
      const H: AnyCurve = HyperellipticCurve(f, h);
      expect(H.genus()).toBe(2);
      const pts = (H as AnyCurve).points();
      const JK = (H.jacobian() as AnyJacobian).point_homset();
      const D = JK.__call__(pts[1]!);
      for (let i = 0; i <= 12; i++) {
        expect(String(D.mul(BigInt(i)))).toBe(tc.seq[i]);
      }
    });
  }
});

/* -------------------------------------------------------------------------- */
/* Group-law verification against brute force                                 */
/* -------------------------------------------------------------------------- */

/**
 * Enumerate every reduced Mumford divisor `(a, b)` on the Jacobian of `C`:
 * `a` monic of degree at most `g`, `deg b < deg a`, and `a | b^2 + h b - f`.
 */
function enumerate_jacobian(H: AnyCurve, JK: AnyHomset): AnyDivisor[] {
  const K = H.base_ring();
  const R = JK.polynomial_ring();
  const [f, h] = H.hyperelliptic_polynomials();
  const g = H.genus();
  const q = Number(K.cardinality!());
  const elements: FF[] = [];
  for (let i = 0n; i < BigInt(q); i++) {
    elements.push(K.__call__(i));
  }

  const out: AnyDivisor[] = [];
  // deg a = 0: the identity
  out.push(JK.__call__(0n));

  const monic = (coeffs: FF[]) => R.__call__([...coeffs, K.one()]);

  for (let d = 1; d <= g; d++) {
    const aIdx = new Array(d).fill(0);
    while (true) {
      const a = monic(aIdx.map((i) => elements[i]!));
      const bIdx = new Array(d).fill(0);
      while (true) {
        const b = R.__call__(bIdx.map((i) => elements[i]!));
        if (b.mul(b).add(h.mul(b)).sub(f).mod(a).isZero()) {
          out.push(JK.__call__([a, b]));
        }
        let k = 0;
        while (k < d) {
          bIdx[k] += 1;
          if (bIdx[k] < q) break;
          bIdx[k] = 0;
          k += 1;
        }
        if (k === d) break;
      }
      let k = 0;
      while (k < d) {
        aIdx[k] += 1;
        if (aIdx[k] < q) break;
        aIdx[k] = 0;
        k += 1;
      }
      if (k === d) break;
    }
  }
  return out;
}

describe('the Jacobian group law over small finite fields', () => {
  const curves: Array<{ p: bigint; f: bigint[]; h: bigint[]; label: string }> = [
    // y^2 = x^5 + x + 1 over GF(5)
    { p: 5n, f: [1n, 1n, 0n, 0n, 0n, 1n], h: [], label: 'y^2 = x^5 + x + 1 / GF(5)' },
    // y^2 = x^5 - 1 over GF(7)
    { p: 7n, f: [6n, 0n, 0n, 0n, 0n, 1n], h: [], label: 'y^2 = x^5 - 1 / GF(7)' },
    // y^2 + x*y = x^5 + 2 over GF(5)
    { p: 5n, f: [2n, 0n, 0n, 0n, 0n, 1n], h: [0n, 1n], label: 'y^2 + xy = x^5 + 2 / GF(5)' },
    // y^2 + (x^2+1) y = x^5 + x over GF(4) is not a prime field; use GF(3)
    {
      p: 3n,
      f: [0n, 1n, 0n, 0n, 0n, 1n],
      h: [1n, 0n, 1n],
      label: 'y^2 + (x^2+1)y = x^5 + x / GF(3)',
    },
  ];

  for (const cv of curves) {
    describe(cv.label, () => {
      const { R, K } = fp(cv.p);
      const f = R.__call__(cv.f.map((v) => K.__call__(v)));
      const h = R.__call__(cv.h.map((v) => K.__call__(v)));
      const H = HyperellipticCurve(f, h.isZero() ? null : h) as AnyCurve;
      const J = H.jacobian() as AnyJacobian;
      const JK = J.point_homset();
      const all = enumerate_jacobian(H as AnyCurve, JK);

      it('has the cardinality predicted by the Frobenius polynomial', () => {
        expect(J.cardinality()).toBe(BigInt(all.length));
      });

      it('is closed, and the enumeration has no duplicates', () => {
        const keys = new Set(all.map(String));
        expect(keys.size).toBe(all.length);
        for (const A of all) {
          for (const B of all) {
            expect(keys.has(String(A.add(B)))).toBe(true);
          }
        }
      });

      it('satisfies identity, inverse and associativity', () => {
        const zero = JK.__call__(0n);
        for (const A of all) {
          expect(A.add(zero).eq(A)).toBe(true);
          expect(zero.add(A).eq(A)).toBe(true);
          expect(A.add(A.neg()).is_zero()).toBe(true);
          expect(A.neg().neg().eq(A)).toBe(true);
        }
        // Associativity is cubic in |J|; check it on every triple for small
        // groups and on a deterministic stride otherwise.
        const stride = all.length <= 30 ? 1 : Math.ceil(all.length / 12);
        const sample = all.filter((_, i) => i % stride === 0);
        for (const A of all) {
          for (const B of all) {
            expect(A.add(B).eq(B.add(A))).toBe(true);
          }
        }
        for (const A of sample) {
          for (const B of sample) {
            for (const Cc of sample) {
              expect(
                A.add(B)
                  .add(Cc)
                  .eq(A.add(B.add(Cc)))
              ).toBe(true);
            }
          }
        }
      });

      it('kills every element with the group order, and orders divide it', () => {
        const N = BigInt(all.length);
        for (const A of all) {
          expect(A.mul(N).is_zero()).toBe(true);
          // brute-force order
          let ord = 1n;
          let cur = A;
          while (!cur.is_zero()) {
            cur = cur.add(A);
            ord += 1n;
            expect(ord).toBeLessThanOrEqual(N);
          }
          expect(N % ord).toBe(0n);
          expect(A.mul(ord).is_zero()).toBe(true);
          expect(A.mul(ord + 1n).eq(A)).toBe(true);
          expect(A.mul(-1n).eq(A.neg())).toBe(true);
          expect(A.mul(ord - 1n).eq(A.neg())).toBe(true);
        }
      });
    });
  }
});

describe('divisor comparison follows SageMath exactly', () => {
  const { R, x, c } = fp(37n);
  void R;
  const f = x.pow(5).add(x).add(c(2n));
  const H: AnyCurve = HyperellipticCurve(f);

  it('two separately built Jacobians give incomparable divisors', () => {
    // sage: J1 = H.jacobian()(GF(37)); J2 = H.jacobian()(GF(37))
    // sage: J1(H.lift_x(...)) == J2(H.lift_x(...))
    // False
    const pt = H.points().find((p: { coords: FF[] }) => !p.coords[2]!.isZero());
    const J1 = (H.jacobian() as AnyJacobian).point_homset();
    const J2 = (H.jacobian() as AnyJacobian).point_homset();
    expect(J1.__call__(pt).eq(J2.__call__(pt))).toBe(false);
    // ... while the two Jacobians themselves are equal, as in Sage
    expect((H.jacobian() as AnyJacobian).eq(H.jacobian())).toBe(true);
  });

  it('the same "point" on different curves is not equal', () => {
    // sage: P1 == P2 -> False   (jacobian_morphism.py:705-719)
    const { x: y } = qq();
    const Ha: AnyCurve = HyperellipticCurve(y.pow(5).sub(y));
    const Hb: AnyCurve = HyperellipticCurve(y.pow(5).add(y));
    const Pa = (Ha.jacobian() as AnyJacobian).point_homset().__call__(Ha.lift_x(new Rational(0n)));
    const Pb = (Hb.jacobian() as AnyJacobian).point_homset().__call__(Hb.lift_x(new Rational(0n)));
    expect(String(Pa)).toBe('(x, y)');
    expect(String(Pb)).toBe('(x, y)');
    expect(Pa.eq(Pb)).toBe(false);
  });
});
