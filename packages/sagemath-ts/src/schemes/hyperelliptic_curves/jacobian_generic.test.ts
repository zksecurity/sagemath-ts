/**
 * Tests for `HyperellipticJacobian_generic`, `JacobianHomset_divisor_classes`
 * and the unported entry points of `hyperelliptic_rational_field`.
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import type { FiniteFieldElement } from '../../rings/finite_rings/finite_field_prime.js';
import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import { Rational } from '../../rings/rational.js';
import { QQ } from '../../rings/rational_field.js';
import { HyperellipticCurve } from './constructor.js';
import { zz_poly_repr } from './hyperelliptic_finite_field.js';
import { HyperellipticJacobian_g2 } from './jacobian_g2.js';
import { HyperellipticJacobian_generic } from './jacobian_generic.js';
import { JacobianHomset_divisor_classes } from './jacobian_homset.js';

type FF = FiniteFieldElement;

// biome-ignore lint/suspicious/noExplicitAny: the constructor returns a dynamic subclass
type AnyCurve = any;
// The repo's `FiniteFieldElement` does not structurally satisfy `RingElement`
// in explicit type positions, so the generic classes are referred to through
// these aliases in the tests.
// biome-ignore lint/suspicious/noExplicitAny: see above
type AnyJacobian = any;
// biome-ignore lint/suspicious/noExplicitAny: see AnyJacobian
type AnyHomset = any;
// biome-ignore lint/suspicious/noExplicitAny: see AnyJacobian
type AnyDivisor = any;

function fp(p: bigint, name = 'x') {
  const K = GF(p);
  const R = new PolynomialRing<FF>(K, name);
  return { K, R, x: R.gen(), c: (n: bigint) => R.__call__(K.__call__(n)) };
}

describe('HyperellipticJacobian_generic', () => {
  const { K, x, c } = fp(37n);
  const H: AnyCurve = HyperellipticCurve(x.pow(5).add(x).add(c(2n)));
  const J = H.jacobian();

  it('is the genus-2 Jacobian and prints like Sage', () => {
    expect(J).toBeInstanceOf(HyperellipticJacobian_g2);
    expect(J).toBeInstanceOf(HyperellipticJacobian_generic);
    expect(String(J)).toBe(`Jacobian of ${H}`);
    expect(J.dimension()).toBe(2n);
    expect(J.curve()).toBe(H);
    expect(J.base_ring()).toBe(K);
  });

  it('J(K) is a homset and is cached', () => {
    const X = J.__call__(K);
    expect(X).toBeInstanceOf(JacobianHomset_divisor_classes);
    expect(J.__call__(K)).toBe(X);
    expect((X as AnyHomset).value_ring()).toBe(K);
    expect((X as AnyHomset).curve()).toBe(H);
    expect((X as AnyHomset).codomain()).toBe(J);
  });

  it('J(0) and J([0]) are the identity', () => {
    expect(String(J.point(0n))).toBe('(1)');
    expect(String(J.__call__(0n))).toBe('(1)');
    expect(String(J.point_homset().__call__([0]))).toBe('(1)');
    expect(J.point_homset().zero().is_zero()).toBe(true);
  });

  it('accepts a pair of curve points as a difference', () => {
    const X = J.point_homset();
    const pts = H.points().filter((p: { coords: FF[] }) => !p.coords[2]!.isZero());
    const P = pts[0]!;
    const Q = pts[2]!;
    expect(String(X.__call__([P, Q]))).toBe(String(X.__call__(P).sub(X.__call__(Q))));
  });

  it('rejects invalid Mumford data with Sage error messages', () => {
    const X = J.point_homset();
    const a = x.pow(2).add(c(1n));
    const b = x; // b^2 + h b - f is not divisible by a
    expect(() => X.__call__([a, b])).toThrow('must be divisor on curve');
    expect(() => X.__call__(3n)).toThrow('does not determine a divisor class');
    expect(() => X.base_extend(null)).toThrow(
      'Jacobian point sets viewed as modules over rings other than ZZ not implemented'
    );
  });

  it('rejects the unported endomorphism helpers', () => {
    expect(() => J.geometric_endomorphism_algebra_is_field()).toThrow(
      'SAGE_NOT_IMPLEMENTED: geometric_endomorphism_algebra_is_field'
    );
    expect(() => J.geometric_endomorphism_ring_is_ZZ()).toThrow(
      'SAGE_NOT_IMPLEMENTED: geometric_endomorphism_ring_is_ZZ'
    );
  });
});

describe('cardinality over QQ is not implemented', () => {
  it('throws NotImplementedError', () => {
    const R = new PolynomialRing<Rational>(QQ as never, 'x');
    const x = R.gen();
    const H: AnyCurve = HyperellipticCurve(x.pow(5).add(x).add(R.one()));
    expect(() => H.jacobian().cardinality()).toThrow(
      'SAGE_NOT_IMPLEMENTED: cardinality of a Jacobian over a non-finite base ring'
    );
  });
});

describe('HyperellipticCurve_rational_field', () => {
  const R = new PolynomialRing<Rational>(QQ as never, 'x');
  const x = R.gen();
  const H: AnyCurve = HyperellipticCurve(
    x
      .pow(5)
      .sub(x.scalar_mul(new Rational(2n)))
      .add(R.__call__(new Rational(3n)))
  );

  it('matrix_of_frobenius is not ported', () => {
    expect(() => H.matrix_of_frobenius(5n)).toThrow('SAGE_NOT_IMPLEMENTED: matrix_of_frobenius');
  });

  it('lseries is not ported', () => {
    expect(() => H.lseries()).toThrow('SAGE_NOT_IMPLEMENTED: lseries');
  });
});

describe('Jacobian cardinality (jacobian_generic.py:420-436)', () => {
  it('matches SageMath over GF(101)', () => {
    // sage: R.<x> = GF(101)[]
    // sage: H = HyperellipticCurve(x^5 + 17*x^4 + 3*x^3 + 11)
    // sage: H.frobenius_polynomial()
    // x^4 - 11*x^3 + 68*x^2 - 1111*x + 10201
    // sage: H.frobenius_polynomial()(1)
    // 9148
    const { K, x, c } = fp(101n);
    const f = x
      .pow(5)
      .add(x.pow(4).scalar_mul(K.__call__(17n)))
      .add(x.pow(3).scalar_mul(K.__call__(3n)))
      .add(c(11n));
    const H: AnyCurve = HyperellipticCurve(f);
    expect(zz_poly_repr(H.frobenius_polynomial())).toBe('x^4 - 11*x^3 + 68*x^2 - 1111*x + 10201');
    expect(H.jacobian().cardinality()).toBe(9148n);
  });

  it('matches SageMath for a genus 3 curve over GF(7)', () => {
    // sage: R.<y> = GF(7)[]
    // sage: H = HyperellipticCurve(y^7 + y + 1)
    // sage: H.frobenius_polynomial()
    // x^6 + 21*x^4 + 147*x^2 + 343
    // sage: H.frobenius_polynomial()(1)
    // 512
    const { x, c } = fp(7n, 'y');
    const H: AnyCurve = HyperellipticCurve(x.pow(7).add(x).add(c(1n)));
    expect(H.genus()).toBe(3);
    expect(zz_poly_repr(H.frobenius_polynomial())).toBe('x^6 + 21*x^4 + 147*x^2 + 343');
    expect(H.jacobian().cardinality()).toBe(512n);
    expect(H.jacobian().dimension()).toBe(3n);
  });
});
