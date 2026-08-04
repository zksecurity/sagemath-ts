/**
 * Tests for the hyperelliptic curve constructor.
 *
 * Expected values come from `sage/schemes/hyperelliptic_curves/constructor.py`
 * and were re-run against the installed SageMath.
 */

import { describe, expect, it } from 'bun:test';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import { GFpn } from '../../rings/finite_rings/finite_field_extension.js';
import type { FiniteFieldElement } from '../../rings/finite_rings/finite_field_prime.js';
import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import { Rational } from '../../rings/rational.js';
import { QQ } from '../../rings/rational_field.js';
import { HyperellipticCurve, _parse_multivariate_defining_equation } from './constructor.js';
import { HyperellipticCurve_finite_field } from './hyperelliptic_finite_field.js';
import {
  HyperellipticCurve_g2,
  HyperellipticCurve_g2_FiniteField,
  HyperellipticCurve_g2_RationalField,
} from './hyperelliptic_g2.js';
import { HyperellipticCurve_generic } from './hyperelliptic_generic.js';
import { HyperellipticCurve_rational_field } from './hyperelliptic_rational_field.js';

type FF = FiniteFieldElement;

function qq(name = 'x') {
  const R = new PolynomialRing<Rational>(QQ as never, name);
  return { R, x: R.gen(), c: (n: bigint) => R.__call__(new Rational(n)) };
}

describe('HyperellipticCurve basic examples (constructor.py:126-200)', () => {
  it('builds curves over QQ', () => {
    const { x, c } = qq();
    expect(String(HyperellipticCurve(x.pow(5).add(x).add(c(1n))))).toBe(
      'Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x + 1'
    );
    expect(String(HyperellipticCurve(x.pow(19).add(x).add(c(1n)), x.sub(c(2n))))).toBe(
      'Hyperelliptic Curve over Rational Field defined by y^2 + (x - 2)*y = x^19 + x + 1'
    );
  });

  it('builds curves over GF(9)', () => {
    const k = GFpn(3n, 2, undefined, 'a');
    const R = new PolynomialRing(k, 'x');
    const x = R.gen();
    const a = R.__call__(k.gen());
    const c = (n: bigint) => R.__call__(k.__call__(n));
    expect(String(HyperellipticCurve(x.pow(3).add(x).sub(c(1n)), x.add(a)))).toBe(
      'Hyperelliptic Curve over Finite Field in a of size 3^2 defined by y^2 + (x + a)*y = x^3 + x + 2'
    );
    expect(
      String(HyperellipticCurve(x.pow(3).add(x).sub(c(1n)), x.add(a), { names: ['X', 'Y'] }))
    ).toBe(
      'Hyperelliptic Curve over Finite Field in a of size 3^2 defined by Y^2 + (X + a)*Y = X^3 + X + 2'
    );
  });

  it('accepts characteristic two with nonzero h', () => {
    const k = GFpn(2n, 3, undefined, 'a');
    const R = new PolynomialRing(k, 'x');
    const x = R.gen();
    const c = (n: bigint) => R.__call__(k.__call__(n));
    expect(String(HyperellipticCurve(x.pow(7).add(c(1n)), x))).toBe(
      'Hyperelliptic Curve over Finite Field in a of size 2^3 defined by y^2 + x*y = x^7 + 1'
    );
    expect(String(HyperellipticCurve(x.pow(8).add(x.pow(7)).add(c(1n)), x.pow(4).add(c(1n))))).toBe(
      'Hyperelliptic Curve over Finite Field in a of size 2^3 defined by y^2 + (x^4 + 1)*y = x^8 + x^7 + 1'
    );
    expect(() => HyperellipticCurve(x.pow(8).add(c(1n)), x)).toThrow(
      'not a hyperelliptic curve: highly singular at infinity'
    );
    expect(() => HyperellipticCurve(x.pow(8).add(x.pow(7)).add(c(1n)), x.pow(4))).toThrow(
      'not a hyperelliptic curve: singularity in the provided affine patch'
    );
  });

  it('rejects h = 0 in characteristic two', () => {
    const K = GF(2n);
    const R = new PolynomialRing<FF>(K, 'x');
    const w = R.gen();
    expect(() => HyperellipticCurve(w.pow(7).add(R.one()))).toThrow(
      'for characteristic 2, argument h = 0 must be nonzero'
    );
  });

  it('allows genus 0 and 1 curves', () => {
    const { x, c } = qq();
    expect(String(HyperellipticCurve(x.pow(2).add(c(1n))))).toBe(
      'Hyperelliptic Curve over Rational Field defined by y^2 = x^2 + 1'
    );
    expect(String(HyperellipticCurve(x.pow(4).sub(c(1n))))).toBe(
      'Hyperelliptic Curve over Rational Field defined by y^2 = x^4 - 1'
    );
    expect(
      String(
        HyperellipticCurve(
          x
            .pow(3)
            .add(x.scalar_mul(new Rational(2n)))
            .add(c(2n))
        )
      )
    ).toBe('Hyperelliptic Curve over Rational Field defined by y^2 = x^3 + 2*x + 2');
  });

  it('rejects double roots but accepts check_squarefree = false', () => {
    const K = GF(7n);
    const R = new PolynomialRing<FF>(K, 'x');
    const x = R.gen();
    const c = (n: bigint) => R.__call__(K.__call__(n));
    const g = x.pow(3).sub(x).add(c(2n));
    const f = g.mul(g).mul(x.pow(6).sub(c(1n)));
    expect(() => HyperellipticCurve(f)).toThrow(
      'not a hyperelliptic curve: singularity in the provided affine patch'
    );
    expect(String(HyperellipticCurve(f, null, { check_squarefree: false }))).toBe(
      'Hyperelliptic Curve over Finite Field of size 7 defined by ' +
        'y^2 = x^12 + 5*x^10 + 4*x^9 + x^8 + 3*x^7 + 3*x^6 + 2*x^4 + 3*x^3 + 6*x^2 + 4*x + 3'
    );
  });

  it('rejects a model of far too high degree (constructor.py:202-218)', () => {
    const { x, c } = qq();
    const F = x.pow(6).add(c(1n));
    const h = x.pow(100);
    const f = F.sub(h.mul(h).scalar_mul(new Rational(1n, 4n)));
    expect(() => HyperellipticCurve(f, h)).toThrow(
      'not a hyperelliptic curve: highly singular at infinity'
    );
    expect(String(HyperellipticCurve(F))).toBe(
      'Hyperelliptic Curve over Rational Field defined by y^2 = x^6 + 1'
    );
  });

  it('accepts a constant f (constructor.py:241-245, issue 15516)', () => {
    const R = new PolynomialRing<Rational>(QQ as never, 'u');
    const u = R.gen();
    const c = (n: bigint) => R.__call__(new Rational(n));
    expect(String(HyperellipticCurve(c(-12n), u.pow(4).add(c(7n))))).toBe(
      'Hyperelliptic Curve over Rational Field defined by y^2 + (x^4 + 7)*y = -12'
    );
  });
});

describe('class specialisation (constructor.py:335-368)', () => {
  it('picks the genus 2 finite field class', () => {
    const K = GF(37n);
    const R = new PolynomialRing<FF>(K, 't');
    const t = R.gen();
    const H = HyperellipticCurve(t.pow(5).add(t).add(R.one()));
    expect(H).toBeInstanceOf(HyperellipticCurve_g2_FiniteField);
    expect(H).toBeInstanceOf(HyperellipticCurve_finite_field);
    expect(H).toBeInstanceOf(HyperellipticCurve_generic);
    expect(H.genus()).toBe(2);
  });

  it('picks the genus 2 rational field class', () => {
    const { x, c } = qq();
    const H = HyperellipticCurve(x.pow(5).add(x).add(c(1n)));
    expect(H).toBeInstanceOf(HyperellipticCurve_g2_RationalField);
    expect(H).toBeInstanceOf(HyperellipticCurve_rational_field);
  });

  it('picks the plain finite field class outside genus 2', () => {
    const K = GF(37n);
    const R = new PolynomialRing<FF>(K, 't');
    const t = R.gen();
    const H = HyperellipticCurve(t.pow(7).add(t).add(R.one()));
    expect(H.genus()).toBe(3);
    expect(H).toBeInstanceOf(HyperellipticCurve_finite_field);
    expect(H).not.toBeInstanceOf(HyperellipticCurve_g2_FiniteField);
    expect(H).not.toBeInstanceOf(HyperellipticCurve_g2);
  });
});

describe('unported entry points', () => {
  it('_parse_multivariate_defining_equation throws NotImplementedError', () => {
    expect(() => _parse_multivariate_defining_equation(null)).toThrow(
      'SAGE_NOT_IMPLEMENTED: _parse_multivariate_defining_equation'
    );
  });
});
