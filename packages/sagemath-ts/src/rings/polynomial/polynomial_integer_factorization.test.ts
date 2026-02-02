/**
 * Tests for polynomial factorization over ZZ (integers).
 *
 * These tests verify that polynomial factorization, root finding, and
 * irreducibility testing work correctly for integer polynomials.
 */

import { describe, expect, test } from 'bun:test';
import { Integer, ZZ } from '../integer_ring.js';
import type { Polynomial } from './polynomial_element.js';
import { PolynomialRing } from './polynomial_ring.js';

// Create a wrapper ring for Integer that works with the polynomial ring
const ZZRing = {
  zero: () => new Integer(0n),
  one: () => new Integer(1n),
  __call__: (x: unknown): Integer => {
    if (x instanceof Integer) return x;
    if (typeof x === 'bigint') return new Integer(x);
    if (typeof x === 'number') return new Integer(BigInt(x));
    if (typeof x === 'string') return new Integer(BigInt(x));
    throw new Error(`Cannot convert ${x} to Integer`);
  },
  is_field: () => false,
  is_integral_domain: () => true,
  characteristic: () => 0n,
  toString: () => 'Integer Ring',
};

// Create polynomial ring ZZ[x]
const ZZx = new PolynomialRing(ZZRing, 'x');

describe('Polynomial factorization over ZZ', () => {
  test('factor monic linear polynomial', () => {
    // f = x + 1
    const f = ZZx.__call__([new Integer(1n), new Integer(1n)]);
    const factors = f.factor();

    expect(factors.length).toBe(1);
    expect(factors[0]![0].degree()).toBe(1);
  });

  test('factor monic quadratic x^2 - 1', () => {
    // f = x^2 - 1 = (x - 1)(x + 1)
    const f = ZZx.__call__([new Integer(-1n), new Integer(0n), new Integer(1n)]);
    const factors = f.factor();

    // Should factor into two linear factors
    const linearFactors = factors.filter(([fac, _]) => fac.degree() === 1);
    expect(linearFactors.length).toBe(2);
  });

  test('factor quadratic that is irreducible over ZZ', () => {
    // f = x^2 + 1 (irreducible over ZZ)
    const f = ZZx.__call__([new Integer(1n), new Integer(0n), new Integer(1n)]);
    const factors = f.factor();

    // Should remain as single factor
    const nonConstantFactors = factors.filter(([fac, _]) => fac.degree() > 0);
    expect(nonConstantFactors.length).toBe(1);
    expect(nonConstantFactors[0]![0].degree()).toBe(2);
  });

  test('factor x^2 - 2x + 1 = (x - 1)^2', () => {
    // f = x^2 - 2x + 1 = (x - 1)^2
    const f = ZZx.__call__([new Integer(1n), new Integer(-2n), new Integer(1n)]);
    const factors = f.factor();

    // Should have (x - 1) with multiplicity 2
    const linearFactors = factors.filter(([fac, _]) => fac.degree() === 1);
    expect(linearFactors.length).toBeGreaterThanOrEqual(1);

    // Check total multiplicity of linear factors
    const totalLinearMult = linearFactors.reduce((sum, [_, mult]) => sum + mult, 0);
    expect(totalLinearMult).toBe(2);
  });
});

describe('Integer polynomial roots', () => {
  test('roots of x', () => {
    const f = ZZx.gen(); // x
    const roots = f.roots();
    expect(roots.length).toBe(1);
    expect((roots[0]![0] as Integer).value).toBe(0n);
    expect(roots[0]![1]).toBe(1);
  });

  test('roots of x^2', () => {
    const x = ZZx.gen();
    const f = x.mul(x); // x^2
    const roots = f.roots();
    expect(roots.length).toBe(1);
    expect((roots[0]![0] as Integer).value).toBe(0n);
    expect(roots[0]![1]).toBe(2); // multiplicity 2
  });

  test('roots of x^2 - 1', () => {
    // x^2 - 1 = (x - 1)(x + 1)
    const f = ZZx.__call__([new Integer(-1n), new Integer(0n), new Integer(1n)]);
    const roots = f.roots();

    expect(roots.length).toBe(2);

    const rootValues = roots.map(([r, _]) => Number((r as Integer).value)).sort((a, b) => a - b);
    expect(rootValues).toEqual([-1, 1]);
  });

  test('roots of polynomial with no integer roots', () => {
    // x^2 + 1 has no integer roots
    const f = ZZx.__call__([new Integer(1n), new Integer(0n), new Integer(1n)]);
    const roots = f.roots();
    expect(roots.length).toBe(0);
  });

  test('roots of x^3 - x', () => {
    // x^3 - x = x(x - 1)(x + 1), roots: 0, 1, -1
    const f = ZZx.__call__([new Integer(0n), new Integer(-1n), new Integer(0n), new Integer(1n)]);
    const roots = f.roots();

    expect(roots.length).toBe(3);

    const rootValues = roots.map(([r, _]) => Number((r as Integer).value)).sort((a, b) => a - b);
    expect(rootValues).toEqual([-1, 0, 1]);
  });
});

describe('Integer polynomial irreducibility', () => {
  test('linear polynomial is irreducible', () => {
    // x + 1 is irreducible
    const f = ZZx.__call__([new Integer(1n), new Integer(1n)]);
    expect(f.is_irreducible()).toBe(true);
  });

  test('non-primitive polynomial is not irreducible', () => {
    // 2x + 2 = 2(x + 1) is not irreducible
    const f = ZZx.__call__([new Integer(2n), new Integer(2n)]);
    expect(f.is_irreducible()).toBe(false);
  });

  test('x^2 + 1 is irreducible over ZZ', () => {
    const f = ZZx.__call__([new Integer(1n), new Integer(0n), new Integer(1n)]);
    expect(f.is_irreducible()).toBe(true);
  });

  test('x^2 - 1 is reducible over ZZ', () => {
    // x^2 - 1 = (x - 1)(x + 1)
    const f = ZZx.__call__([new Integer(-1n), new Integer(0n), new Integer(1n)]);
    expect(f.is_irreducible()).toBe(false);
  });

  test('x^2 - 2 is irreducible over ZZ', () => {
    // x^2 - 2 is irreducible (sqrt(2) is irrational)
    const f = ZZx.__call__([new Integer(-2n), new Integer(0n), new Integer(1n)]);
    expect(f.is_irreducible()).toBe(true);
  });

  test('cyclotomic polynomial Phi_3 = x^2 + x + 1 is irreducible', () => {
    const f = ZZx.__call__([new Integer(1n), new Integer(1n), new Integer(1n)]);
    expect(f.is_irreducible()).toBe(true);
  });
});

// Helper function to extract leading coefficient
function extractLC(poly: Polynomial<Integer>): bigint {
  const lc = poly.leading_coefficient();
  return (lc as Integer).value;
}
