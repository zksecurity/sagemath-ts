/**
 * Unit tests for polynomial interpolation methods
 */
import { describe, expect, test } from 'bun:test';
import { ZeroDivisionError } from '../../errors.js';
import { FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import { GF2 } from '../finite_rings/gf2.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';
import { PolynomialRingConstructor } from './polynomial_ring.js';

// Create test fields
const F7 = new FiniteFieldPrime(7n);
const F13 = new FiniteFieldPrime(13n);
const F101 = new FiniteFieldPrime(101n);

describe('Lagrange interpolation over F7', () => {
  const [R, x] = PolynomialRingConstructor(F7, 'x');

  test('empty points returns zero', () => {
    const p = R.lagrange_polynomial([]);
    expect(p.isZero()).toBe(true);
  });

  test('single point returns constant', () => {
    const p = R.lagrange_polynomial([[F7.__call__(3), F7.__call__(5)]]);
    expect(p.degree()).toBe(0);
    expect(p.getCoeff(0).eq(5)).toBe(true);
  });

  test('two points gives line', () => {
    // Points (0, 1) and (1, 3)
    const p = R.lagrange_polynomial([
      [F7.__call__(0), F7.__call__(1)],
      [F7.__call__(1), F7.__call__(3)],
    ]);

    // Line through these points: y = 2x + 1
    expect(p.degree()).toBe(1);
    expect(p.evaluate(F7.__call__(0)).eq(1)).toBe(true);
    expect(p.evaluate(F7.__call__(1)).eq(3)).toBe(true);
  });

  test('three points gives parabola', () => {
    // Points (0, 1), (1, 2), (2, 5) -> y = x^2 + 1
    const p = R.lagrange_polynomial([
      [F7.__call__(0), F7.__call__(1)],
      [F7.__call__(1), F7.__call__(2)],
      [F7.__call__(2), F7.__call__(5)],
    ]);

    expect(p.degree()).toBe(2);
    expect(p.evaluate(F7.__call__(0)).eq(1)).toBe(true);
    expect(p.evaluate(F7.__call__(1)).eq(2)).toBe(true);
    expect(p.evaluate(F7.__call__(2)).eq(5)).toBe(true);

    // Verify it's x^2 + 1 (mod 7)
    expect(p.evaluate(F7.__call__(3)).eq(3 * 3 + 1)).toBe(true); // 10 mod 7 = 3
  });

  test('interpolation round-trip', () => {
    // Create a polynomial and verify we can recover it from point evaluations
    // p(x) = x^3 + 2x + 4
    const p = x
      .pow(3)
      .add(x.scalar_mul(F7.__call__(2)))
      .add(R.__call__(F7.__call__(4)));

    // Evaluate at 4 points (degree 3 polynomial needs 4 points)
    const points: Array<
      [
        typeof F7 extends CoefficientRing<infer T> ? T : never,
        typeof F7 extends CoefficientRing<infer T> ? T : never,
      ]
    > = [];
    for (let i = 0; i < 4; i++) {
      const xi = F7.__call__(i);
      points.push([xi, p.evaluate(xi)]);
    }

    // Interpolate and verify we get the same polynomial back
    const q = R.lagrange_polynomial(points);
    expect(q.eq(p)).toBe(true);
  });

  test('duplicate x values throws error', () => {
    // Sage has no precheck: a repeated x makes ``x_i - x_j`` zero and the
    // division in the divided-difference table raises ZeroDivisionError
    // (`polynomial_ring.py:2295` requires x_i - x_j to be invertible).
    expect(() => {
      R.lagrange_polynomial([
        [F7.__call__(1), F7.__call__(2)],
        [F7.__call__(1), F7.__call__(3)], // duplicate x = 1
      ]);
    }).toThrow(ZeroDivisionError);
  });

  test('divided_difference algorithm matches neville', () => {
    const points: Array<
      [
        typeof F7 extends CoefficientRing<infer T> ? T : never,
        typeof F7 extends CoefficientRing<infer T> ? T : never,
      ]
    > = [
      [F7.__call__(0), F7.__call__(1)],
      [F7.__call__(2), F7.__call__(2)],
      [F7.__call__(3), F7.__call__(6)],
      [F7.__call__(5), F7.__call__(1)],
    ];

    const p1 = R.lagrange_polynomial(points, 'divided_difference');
    // ``algorithm='neville'`` returns the whole last row of the Neville table
    // (a list), exactly as Sage does; the interpolating polynomial is its last
    // entry (`polynomial_ring.py:2373-2378`).
    const row = R.lagrange_polynomial(points, 'neville');
    expect(Array.isArray(row)).toBe(true);
    expect(row.length).toBe(points.length);

    expect(p1.eq(row[row.length - 1]!)).toBe(true);
  });
});

describe('Newton interpolation over F13', () => {
  const [R, x] = PolynomialRingConstructor(F13, 'x');

  test('newton_interpolation equals lagrange_polynomial', () => {
    const points: Array<
      [
        typeof F13 extends CoefficientRing<infer T> ? T : never,
        typeof F13 extends CoefficientRing<infer T> ? T : never,
      ]
    > = [
      [F13.__call__(0), F13.__call__(1)],
      [F13.__call__(1), F13.__call__(5)],
      [F13.__call__(2), F13.__call__(11)],
      [F13.__call__(3), F13.__call__(7)],
    ];

    const newton = R.newton_interpolation(points);
    const lagrange = R.lagrange_polynomial(points);

    expect(newton.eq(lagrange)).toBe(true);
  });

  test('divided_difference computation', () => {
    // For points (1, 2), (3, 4), the divided differences are:
    // F[0,0] = 2
    // F[1,0] = 4, F[1,1] = (4-2)/(3-1) = 2/2 = 1
    const points: Array<
      [
        typeof F13 extends CoefficientRing<infer T> ? T : never,
        typeof F13 extends CoefficientRing<infer T> ? T : never,
      ]
    > = [
      [F13.__call__(1), F13.__call__(2)],
      [F13.__call__(3), F13.__call__(4)],
    ];

    const coeffs = R.divided_difference(points);
    expect(coeffs.length).toBe(2);
    expect(coeffs[0]!.eq(2)).toBe(true);
    expect(coeffs[1]!.eq(1)).toBe(true);
  });
});

describe('Vanishing polynomial', () => {
  const [R, x] = PolynomialRingConstructor(F7, 'x');

  test('empty domain returns 1', () => {
    const Z = R.vanishing_polynomial([]);
    expect(Z.eq(R.one())).toBe(true);
  });

  test('single element domain', () => {
    // Z_H(x) = x - 3
    const Z = R.vanishing_polynomial([F7.__call__(3)]);
    expect(Z.degree()).toBe(1);
    expect(Z.evaluate(F7.__call__(3)).isZero()).toBe(true);
  });

  test('evaluates to zero on domain', () => {
    const domain = [F7.__call__(1), F7.__call__(2), F7.__call__(4)];
    const Z = R.vanishing_polynomial(domain);

    expect(Z.degree()).toBe(3);

    // Should evaluate to 0 on all domain elements
    for (const h of domain) {
      expect(Z.evaluate(h).isZero()).toBe(true);
    }

    // Should be non-zero on other elements
    expect(Z.evaluate(F7.__call__(0)).isZero()).toBe(false);
    expect(Z.evaluate(F7.__call__(3)).isZero()).toBe(false);
  });

  test('is monic', () => {
    const domain = [F7.__call__(1), F7.__call__(2), F7.__call__(3)];
    const Z = R.vanishing_polynomial(domain);

    expect(Z.is_monic()).toBe(true);
  });

  test('equals product of linear factors', () => {
    // (x - 1)(x - 2)(x - 3)
    const Z = R.vanishing_polynomial([F7.__call__(1), F7.__call__(2), F7.__call__(3)]);

    const expected = x
      .sub(R.__call__(F7.__call__(1)))
      .mul(x.sub(R.__call__(F7.__call__(2))))
      .mul(x.sub(R.__call__(F7.__call__(3))));

    expect(Z.eq(expected)).toBe(true);
  });
});

describe('Barycentric interpolation', () => {
  const [R, x] = PolynomialRingConstructor(F7, 'x');

  test('empty points returns zero', () => {
    const val = R.barycentric_interpolation([], F7.__call__(3));
    expect(val.isZero()).toBe(true);
  });

  test('single point returns y value', () => {
    const val = R.barycentric_interpolation([[F7.__call__(2), F7.__call__(5)]], F7.__call__(3));
    // For a single point, the interpolant is constant
    expect(val.eq(5)).toBe(true);
  });

  test('evaluation at data point returns y directly', () => {
    const points: Array<
      [
        typeof F7 extends CoefficientRing<infer T> ? T : never,
        typeof F7 extends CoefficientRing<infer T> ? T : never,
      ]
    > = [
      [F7.__call__(1), F7.__call__(3)],
      [F7.__call__(2), F7.__call__(5)],
      [F7.__call__(4), F7.__call__(1)],
    ];

    // Evaluating at x = 2 should return y = 5
    const val = R.barycentric_interpolation(points, F7.__call__(2));
    expect(val.eq(5)).toBe(true);
  });

  test('matches polynomial evaluation', () => {
    const points: Array<
      [
        typeof F7 extends CoefficientRing<infer T> ? T : never,
        typeof F7 extends CoefficientRing<infer T> ? T : never,
      ]
    > = [
      [F7.__call__(0), F7.__call__(1)],
      [F7.__call__(1), F7.__call__(2)],
      [F7.__call__(2), F7.__call__(5)],
    ];

    // Build full polynomial
    const p = R.lagrange_polynomial(points);

    // Evaluate using barycentric at multiple points
    for (let i = 0; i < 7; i++) {
      const evalPoint = F7.__call__(i);
      const baryVal = R.barycentric_interpolation(points, evalPoint);
      const polyVal = p.evaluate(evalPoint);

      expect(baryVal.eq(polyVal)).toBe(true);
    }
  });
});

describe('from_roots', () => {
  const [R, x] = PolynomialRingConstructor(F7, 'x');

  test('empty roots returns 1', () => {
    const p = R.from_roots([]);
    expect(p.eq(R.one())).toBe(true);
  });

  test('single root', () => {
    // (x - 3)
    const p = R.from_roots([F7.__call__(3)]);
    expect(p.degree()).toBe(1);
    expect(p.evaluate(F7.__call__(3)).isZero()).toBe(true);
  });

  test('multiple roots', () => {
    const roots = [F7.__call__(1), F7.__call__(2), F7.__call__(5)];
    const p = R.from_roots(roots);

    expect(p.degree()).toBe(3);
    expect(p.is_monic()).toBe(true);

    // Verify all roots
    for (const r of roots) {
      expect(p.evaluate(r).isZero()).toBe(true);
    }
  });

  test('equals vanishing_polynomial', () => {
    const roots = [F7.__call__(1), F7.__call__(3), F7.__call__(4)];
    const p1 = R.from_roots(roots);
    const p2 = R.vanishing_polynomial(roots);

    expect(p1.eq(p2)).toBe(true);
  });
});

describe('Cyclotomic polynomials', () => {
  const [R, x] = PolynomialRingConstructor(F101, 'x');

  test('Phi_1(x) = x - 1', () => {
    const phi1 = R.cyclotomic_polynomial(1);
    const expected = x.sub(R.one());

    expect(phi1.eq(expected)).toBe(true);
  });

  test('Phi_2(x) = x + 1', () => {
    const phi2 = R.cyclotomic_polynomial(2);
    const expected = x.add(R.one());

    expect(phi2.eq(expected)).toBe(true);
  });

  test('Phi_3(x) = x^2 + x + 1', () => {
    const phi3 = R.cyclotomic_polynomial(3);
    const expected = x.pow(2).add(x).add(R.one());

    expect(phi3.eq(expected)).toBe(true);
  });

  test('Phi_4(x) = x^2 + 1', () => {
    const phi4 = R.cyclotomic_polynomial(4);
    const expected = x.pow(2).add(R.one());

    expect(phi4.eq(expected)).toBe(true);
  });

  test('Phi_6(x) = x^2 - x + 1', () => {
    const phi6 = R.cyclotomic_polynomial(6);
    // x^2 - x + 1 = x^2 + (-1)*x + 1
    const expected = x.pow(2).sub(x).add(R.one());

    expect(phi6.eq(expected)).toBe(true);
  });

  test('Phi_8(x) = x^4 + 1', () => {
    const phi8 = R.cyclotomic_polynomial(8);
    const expected = x.pow(4).add(R.one());

    expect(phi8.eq(expected)).toBe(true);
  });

  test('degree of Phi_n equals euler_phi(n)', () => {
    const testCases: [number, number][] = [
      [1, 1],
      [2, 1],
      [3, 2],
      [4, 2],
      [5, 4],
      [6, 2],
      [7, 6],
      [8, 4],
      [9, 6],
      [10, 4],
      [12, 4],
    ];

    for (const [n, expectedDegree] of testCases) {
      const phi = R.cyclotomic_polynomial(n);
      expect(phi.degree()).toBe(expectedDegree);
    }
  });

  test('x^n - 1 = product of Phi_d for d | n', () => {
    // Test for n = 6: x^6 - 1 = Phi_1 * Phi_2 * Phi_3 * Phi_6
    const n = 6;
    const divisors = [1, 2, 3, 6];

    let product = R.one();
    for (const d of divisors) {
      product = product.mul(R.cyclotomic_polynomial(d));
    }

    const xNMinus1 = x.pow(n).sub(R.one());

    expect(product.eq(xNMinus1)).toBe(true);
  });

  test('Phi_n(x) is monic', () => {
    for (let n = 1; n <= 15; n++) {
      const phi = R.cyclotomic_polynomial(n);
      expect(phi.is_monic()).toBe(true);
    }
  });

  test('invalid n throws error', () => {
    // sage: R.cyclotomic_polynomial(0)
    // ArithmeticError: n=0 must be positive          (polynomial_ring.py:1188)
    expect(() => R.cyclotomic_polynomial(0)).toThrow('n=0 must be positive');
    expect(() => R.cyclotomic_polynomial(-1)).toThrow('n=-1 must be positive');
    expect(() => R.cyclotomic_polynomial(1.5)).toThrow('must be an integer');
  });

  test('Phi_30 (product of 3 distinct primes)', () => {
    // Phi_30(x) = x^8 + x^7 - x^5 - x^4 - x^3 + x + 1
    const phi30 = R.cyclotomic_polynomial(30);
    expect(phi30.degree()).toBe(8); // euler_phi(30) = 8

    // Verify it divides x^30 - 1
    const x30Minus1 = x.pow(30).sub(R.one());
    const [q, r] = x30Minus1.quo_rem(phi30);
    expect(r.isZero()).toBe(true);
  });
});

describe('Cyclotomic polynomials over GF(2)', () => {
  const [R, x] = PolynomialRingConstructor(GF2, 'x');

  test('Phi_1(x) = x + 1 (since -1 = 1 in GF(2))', () => {
    const phi1 = R.cyclotomic_polynomial(1);
    // x - 1 = x + 1 in GF(2)
    const expected = x.add(R.one());
    expect(phi1.eq(expected)).toBe(true);
  });

  test('Phi_3(x) = x^2 + x + 1', () => {
    const phi3 = R.cyclotomic_polynomial(3);
    const expected = x.pow(2).add(x).add(R.one());
    expect(phi3.eq(expected)).toBe(true);
  });

  test('Phi_7(x) = x^6 + x^5 + x^4 + x^3 + x^2 + x + 1', () => {
    const phi7 = R.cyclotomic_polynomial(7);
    expect(phi7.degree()).toBe(6);

    // In GF(2), Phi_7(x) = 1 + x + x^2 + x^3 + x^4 + x^5 + x^6
    let expected = R.one();
    for (let i = 1; i <= 6; i++) {
      expected = expected.add(x.pow(i));
    }
    expect(phi7.eq(expected)).toBe(true);
  });
});

describe('Interpolation edge cases over larger field', () => {
  const [R, x] = PolynomialRingConstructor(F101, 'x');

  test('high degree interpolation', () => {
    // Create a polynomial of degree 10
    let p = R.zero();
    for (let i = 0; i <= 10; i++) {
      p = p.add(x.pow(i).scalar_mul(F101.__call__(i + 1)));
    }

    // Evaluate at 11 points
    const points: Array<
      [
        typeof F101 extends CoefficientRing<infer T> ? T : never,
        typeof F101 extends CoefficientRing<infer T> ? T : never,
      ]
    > = [];
    for (let i = 0; i < 11; i++) {
      const xi = F101.__call__(i);
      points.push([xi, p.evaluate(xi)]);
    }

    // Interpolate
    const q = R.lagrange_polynomial(points);
    expect(q.eq(p)).toBe(true);
  });

  test('interpolation preserves zero polynomial', () => {
    // All y values are 0
    const points: Array<
      [
        typeof F101 extends CoefficientRing<infer T> ? T : never,
        typeof F101 extends CoefficientRing<infer T> ? T : never,
      ]
    > = [
      [F101.__call__(1), F101.__call__(0)],
      [F101.__call__(2), F101.__call__(0)],
      [F101.__call__(3), F101.__call__(0)],
    ];

    const p = R.lagrange_polynomial(points);
    expect(p.isZero()).toBe(true);
  });

  test('interpolation preserves constant polynomial', () => {
    // All y values are the same
    const c = F101.__call__(42);
    const points: Array<
      [
        typeof F101 extends CoefficientRing<infer T> ? T : never,
        typeof F101 extends CoefficientRing<infer T> ? T : never,
      ]
    > = [
      [F101.__call__(1), c],
      [F101.__call__(2), c],
      [F101.__call__(3), c],
    ];

    const p = R.lagrange_polynomial(points);
    expect(p.degree()).toBe(0);
    expect(p.getCoeff(0).eq(42)).toBe(true);
  });
});
