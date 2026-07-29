/**
 * Tests for `HyperellipticCurve_generic`.
 *
 * Expected values are the docstring examples of
 * `sage/schemes/hyperelliptic_curves/hyperelliptic_generic.py`, re-run against
 * the installed SageMath.
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import { GFpn } from '../../rings/finite_rings/finite_field_extension.js';
import type { FiniteFieldElement } from '../../rings/finite_rings/finite_field_prime.js';
import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import { Rational } from '../../rings/rational.js';
import { QQ } from '../../rings/rational_field.js';
import { HyperellipticCurve } from './constructor.js';
import { HyperellipticPoint, sage_poly_repr } from './hyperelliptic_generic.js';

type FF = FiniteFieldElement;

function qq(name = 'x') {
  const R = new PolynomialRing<Rational>(QQ as never, name);
  return { R, x: R.gen(), c: (n: bigint) => R.__call__(new Rational(n)) };
}
function fp(p: bigint, name = 'x') {
  const K = GF(p);
  const R = new PolynomialRing<FF>(K, name);
  return { K, R, x: R.gen(), c: (n: bigint) => R.__call__(K.__call__(n)) };
}

// biome-ignore lint/suspicious/noExplicitAny: the constructor returns a dynamic subclass
type AnyCurve = any;

describe('_repr_ (hyperelliptic_generic.py:145-167)', () => {
  it('renames the coordinates', () => {
    const { x, c } = qq();
    const f = x
      .pow(5)
      .scalar_mul(new Rational(4n))
      .sub(x.pow(3).scalar_mul(new Rational(30n)))
      .add(x.scalar_mul(new Rational(45n)))
      .sub(c(22n));
    expect(String(HyperellipticCurve(f))).toBe(
      'Hyperelliptic Curve over Rational Field defined by y^2 = 4*x^5 - 30*x^3 + 45*x - 22'
    );
    expect(String(HyperellipticCurve(f, null, { names: 'u,v' }))).toBe(
      'Hyperelliptic Curve over Rational Field defined by v^2 = 4*u^5 - 30*u^3 + 45*u - 22'
    );
    expect(String(HyperellipticCurve(x.pow(5).add(c(1n)), x.pow(3).add(c(2n))))).toBe(
      'Hyperelliptic Curve over Rational Field defined by y^2 + (x^3 + 2)*y = x^5 + 1'
    );
    expect(HyperellipticCurve(f).genus()).toBe(2);
  });
});

describe('sage_poly_repr reproduces Sage polynomial printing', () => {
  it('handles signs, unit coefficients and parentheses', () => {
    // (-t^3 + 1)*x^3 - t^2*x^2 - x + 1  (polynomial_element.pyx:3106-3110)
    expect(sage_poly_repr(['1', '-1', '-t^2', '-t^3 + 1'], 'x')).toBe(
      '(-t^3 + 1)*x^3 - t^2*x^2 - x + 1'
    );
    expect(sage_poly_repr([], 'x')).toBe('0');
    expect(sage_poly_repr([null, '1'], 'y')).toBe('y');
    expect(sage_poly_repr(['26', '1'], 'y')).toBe('y + 26');
    expect(sage_poly_repr(['-1', '1'], 'y')).toBe('y - 1');
  });
});

describe('is_singular / is_smooth (hyperelliptic_generic.py:212-267)', () => {
  it('always reports a smooth curve', () => {
    const { x, c } = qq();
    const H: AnyCurve = HyperellipticCurve(x.pow(5).add(c(1n)));
    expect(H.is_singular()).toBe(false);
    expect(H.is_smooth()).toBe(true);
  });
});

describe('is_x_coord and lift_x (hyperelliptic_generic.py:269-524)', () => {
  it('works over QQ for h = 0', () => {
    const { x, c } = qq();
    const H: AnyCurve = HyperellipticCurve(x.pow(5).add(x.pow(3)).add(c(1n)));
    expect(H.is_x_coord(0n)).toBe(true);
    expect(H.is_x_coord(3n)).toBe(false);
    expect(String(H.lift_x(0n))).toBe('(0 : -1 : 1)');
    expect(H.lift_x(4n, { all: true }).map(String)).toEqual(['(4 : -33 : 1)', '(4 : 33 : 1)']);
    expect(H.lift_x(3n, { all: true })).toEqual([]);
    expect(() => H.lift_x(3n)).toThrow(
      'No point with x-coordinate 3 on Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x^3 + 1'
    );
  });

  it('works over QQ for h nonzero', () => {
    const { x, c } = qq();
    const H: AnyCurve = HyperellipticCurve(x.pow(5).add(x.pow(3)).add(c(1n)), x.add(c(1n)));
    expect(H.is_x_coord(1n)).toBe(true);
    expect(String(H.lift_x(1n))).toBe('(1 : -3 : 1)');
  });

  it('works over GF(163)', () => {
    const { x, c } = fp(163n);
    const H: AnyCurve = HyperellipticCurve(x.pow(7).add(x).add(c(1n)));
    expect(H.is_x_coord(13n)).toBe(true);
    expect(String(H.lift_x(13n))).toBe('(13 : 41 : 1)');
  });

  it('works in characteristic two (issue 37097)', () => {
    const F = GFpn(2n, 4, undefined, 'z4');
    const R = new PolynomialRing(F, 'x');
    const x = R.gen();
    const c = (n: bigint) => R.__call__(F.__call__(n));
    const H: AnyCurve = HyperellipticCurve(x.pow(7).add(x.pow(3)).add(c(1n)), x.add(c(1n)));
    const z4 = F.gen();
    const x0 = z4.pow(3).add(z4.pow(2)).add(z4);
    expect(H.is_x_coord(x0)).toBe(true);
    expect(H.lift_x(x0, { all: true }).map(String)).toEqual([
      '(z4^3 + z4^2 + z4 : z4^2 + z4 + 1 : 1)',
      '(z4^3 + z4^2 + z4 : z4^3 : 1)',
    ]);
  });
});

describe('odd_degree_model (hyperelliptic_generic.py:533-635)', () => {
  it('returns the curve unchanged when it is already odd', () => {
    const { x, c } = qq();
    const H: AnyCurve = HyperellipticCurve(x.pow(5).add(c(1n)), null, { names: 'U, V' });
    expect(String(H.odd_degree_model())).toBe(
      'Hyperelliptic Curve over Rational Field defined by V^2 = U^5 + 1'
    );
  });

  it('moves a rational root to infinity', () => {
    const { x } = qq();
    const H: AnyCurve = HyperellipticCurve(x.pow(6).add(x));
    expect(String(H.odd_degree_model())).toBe(
      'Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + 1'
    );
  });

  it('raises when no rational root exists', () => {
    const { x, c } = qq();
    const sextic = x
      .pow(2)
      .add(c(2n))
      .mul(x.pow(2).add(c(3n)))
      .mul(x.pow(2).add(c(5n)));
    const H: AnyCurve = HyperellipticCurve(sextic);
    expect(() => H.odd_degree_model()).toThrow(
      'No odd degree model exists over field of definition'
    );
    expect(H.has_odd_degree_model()).toBe(false);
  });

  it('is not implemented for a nonzero h', () => {
    const { x, c } = qq();
    const H: AnyCurve = HyperellipticCurve(x.pow(5).add(c(1n)), c(1n));
    expect(() => H.odd_degree_model()).toThrow(
      'odd_degree_model only implemented for curves in Weierstrass form'
    );
  });

  it('has_odd_degree_model matches the docstring', () => {
    const { x, c } = qq();
    expect((HyperellipticCurve(x.pow(5).add(x)) as AnyCurve).has_odd_degree_model()).toBe(true);
    expect((HyperellipticCurve(x.pow(6).add(x)) as AnyCurve).has_odd_degree_model()).toBe(true);
    expect(
      (HyperellipticCurve(x.pow(6).add(x).add(c(1n))) as AnyCurve).has_odd_degree_model()
    ).toBe(false);
  });
});

describe('change_ring / base_extend (hyperelliptic_generic.py:113-143)', () => {
  it('base extends GF(7) to GF(49)', () => {
    const { x, c } = fp(7n);
    const H: AnyCurve = HyperellipticCurve(x.pow(8).add(x).add(c(5n)));
    const F49 = GFpn(7n, 2, undefined, 'a');
    expect(String(H.base_extend(F49))).toBe(
      'Hyperelliptic Curve over Finite Field in a of size 7^2 defined by y^2 = x^8 + x + 5'
    );
  });
});

describe('projective points', () => {
  it('normalises by the last nonzero coordinate', () => {
    const { K, x, c } = fp(11n);
    // y^2 + (2x^4+1) y = 3x^5 + x has the infinity point (8 : 1 : 0) in Sage
    const f = x.pow(5).scalar_mul(K.__call__(3n)).add(x);
    const h = x.pow(4).scalar_mul(K.__call__(2n)).add(c(1n));
    const H: AnyCurve = HyperellipticCurve(f, h);
    const pts = H.points().map(String);
    expect(pts).toEqual([
      '(0 : 1 : 0)',
      '(8 : 1 : 0)',
      '(0 : 0 : 1)',
      '(0 : 10 : 1)',
      '(1 : 1 : 1)',
      '(1 : 7 : 1)',
      '(4 : 2 : 1)',
      '(9 : 10 : 1)',
      '(9 : 1 : 1)',
      '(10 : 3 : 1)',
      '(10 : 5 : 1)',
    ]);
    expect(H.genus()).toBe(3);
  });

  it('rejects [0, 0, 0]', () => {
    const { K } = fp(7n);
    expect(() => HyperellipticPoint.normalize([K.zero(), K.zero(), K.zero()])).toThrow(
      '[0, 0, 0] does not define a valid projective point'
    );
  });
});

describe('unported generic methods throw NotImplementedError', () => {
  const { x, c } = qq();
  const H: AnyCurve = HyperellipticCurve(x.pow(5).add(c(1n)));
  for (const name of [
    'monsky_washnitzer_gens',
    'invariant_differential',
    'local_coordinates_at_nonweierstrass',
    'local_coordinates_at_weierstrass',
    'local_coordinates_at_infinity',
    'local_coord',
    'rational_points',
  ]) {
    it(`${name}`, () => {
      expect(() => H[name]()).toThrow('SAGE_NOT_IMPLEMENTED');
    });
  }
});
