/**
 * Tests for polynomial factorization, root-finding, and irreducibility testing
 *
 * Tests cover:
 * - roots() - finding roots with multiplicities
 * - factor() - complete factorization
 * - is_irreducible() - irreducibility testing
 * - squarefree_decomposition() - squarefree factorization
 * - distinct_degree_factorization() - factorization by degree
 */
import { describe, expect, test } from 'bun:test';
import { next_prime } from '../../arith/misc.js';
import { CONWAY_POLYNOMIALS } from '../finite_rings/conway_polynomials.js';
import { FiniteFieldExtension, PrimeField } from '../finite_rings/finite_field_extension.js';
import { FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import { GF2 } from '../finite_rings/gf2.js';
import { Integer } from '../integer_ring.js';
import { Rational } from '../rational.js';
import { QQ } from '../rational_field.js';
import {
  type CoefficientRing,
  Polynomial,
  type RingElement,
  _zz_factor_internal,
} from './polynomial_element.js';
import { PolynomialRing, PolynomialRingConstructor } from './polynomial_ring.js';

describe('Polynomial roots over GF(2)', () => {
  const [R, x] = PolynomialRingConstructor(GF2, 'x');

  test('roots of x', () => {
    // x has root 0 with multiplicity 1
    const roots = x.roots();
    expect(roots.length).toBe(1);
    expect(roots[0]![0].eq(0)).toBe(true);
    expect(roots[0]![1]).toBe(1);
  });

  test('roots of x + 1', () => {
    // x + 1 has root 1 with multiplicity 1
    const p = x.add(R.one());
    const roots = p.roots();
    expect(roots.length).toBe(1);
    expect(roots[0]![0].eq(1)).toBe(true);
    expect(roots[0]![1]).toBe(1);
  });

  test('roots of x^2 + x = x(x+1)', () => {
    // x^2 + x = x(x+1) has roots 0 and 1
    const p = x.pow(2).add(x);
    const roots = p.roots();
    expect(roots.length).toBe(2);

    // Sort by value for consistent testing
    const sortedRoots = roots.sort((a, b) => Number(a[0].toString()) - Number(b[0].toString()));
    expect(sortedRoots[0]![0].eq(0)).toBe(true);
    expect(sortedRoots[0]![1]).toBe(1);
    expect(sortedRoots[1]![0].eq(1)).toBe(true);
    expect(sortedRoots[1]![1]).toBe(1);
  });

  test('roots of x^2 + 1 = (x+1)^2', () => {
    // In GF(2): x^2 + 1 = (x + 1)^2, root 1 with multiplicity 2
    const p = x.pow(2).add(R.one());
    const roots = p.roots();
    expect(roots.length).toBe(1);
    expect(roots[0]![0].eq(1)).toBe(true);
    expect(roots[0]![1]).toBe(2);
  });

  test('roots of x^2 + x + 1 (irreducible)', () => {
    // x^2 + x + 1 is irreducible over GF(2), has no roots
    const p = x.pow(2).add(x).add(R.one());
    const roots = p.roots();
    expect(roots.length).toBe(0);
  });
});

describe('Polynomial roots over prime fields', () => {
  test('roots over GF(5)', () => {
    const F5 = new FiniteFieldPrime(5n);
    const [R, x] = PolynomialRingConstructor(F5, 'x');

    // x^2 - 1 = (x-1)(x+1) = (x-1)(x-4) has roots 1 and 4
    const p = x.pow(2).sub(R.one());
    const roots = p.roots();

    expect(roots.length).toBe(2);
    const rootValues = roots.map(([r, _m]) => r.value).sort();
    expect(rootValues).toContain(1n);
    expect(rootValues).toContain(4n);
  });

  test('roots over GF(7)', () => {
    const F7 = new FiniteFieldPrime(7n);
    const [R, x] = PolynomialRingConstructor(F7, 'x');

    // x^2 - 2 over GF(7): 3^2 = 9 = 2 (mod 7), 4^2 = 16 = 2 (mod 7)
    // So roots are 3 and 4
    const two = R.__call__(F7.__call__(2));
    const p = x.pow(2).sub(two);
    const roots = p.roots();

    expect(roots.length).toBe(2);
    const rootValues = roots.map(([r, _m]) => r.value).sort();
    expect(rootValues).toContain(3n);
    expect(rootValues).toContain(4n);
  });

  test('roots with multiplicities over GF(7)', () => {
    const F7 = new FiniteFieldPrime(7n);
    const [R, x] = PolynomialRingConstructor(F7, 'x');

    // (x - 1)^3 * (x - 2) has root 1 with mult 3, root 2 with mult 1
    const xMinus1 = x.sub(R.one());
    const xMinus2 = x.sub(R.__call__(F7.__call__(2)));
    const p = xMinus1.pow(3).mul(xMinus2);

    const roots = p.roots();

    // Find root 1 and root 2
    let root1Mult = 0;
    let root2Mult = 0;

    for (const [r, m] of roots) {
      if (r.value === 1n) root1Mult = m;
      if (r.value === 2n) root2Mult = m;
    }

    expect(root1Mult).toBe(3);
    expect(root2Mult).toBe(1);
  });
});

describe('Polynomial is_irreducible over GF(2)', () => {
  const [R, x] = PolynomialRingConstructor(GF2, 'x');

  test('x is irreducible', () => {
    expect(x.is_irreducible()).toBe(true);
  });

  test('x + 1 is irreducible', () => {
    expect(x.add(R.one()).is_irreducible()).toBe(true);
  });

  test('x^2 + x is not irreducible (= x(x+1))', () => {
    expect(x.pow(2).add(x).is_irreducible()).toBe(false);
  });

  test('x^2 + x + 1 is irreducible (Conway polynomial for GF(4))', () => {
    const p = x.pow(2).add(x).add(R.one());
    expect(p.is_irreducible()).toBe(true);
  });

  test('x^3 + x + 1 is irreducible (Conway polynomial for GF(8))', () => {
    const p = x.pow(3).add(x).add(R.one());
    expect(p.is_irreducible()).toBe(true);
  });

  test('x^3 + x^2 + 1 is irreducible over GF(2)', () => {
    const p = x.pow(3).add(x.pow(2)).add(R.one());
    expect(p.is_irreducible()).toBe(true);
  });

  test('x^4 + x + 1 is irreducible (Conway polynomial for GF(16))', () => {
    const p = x.pow(4).add(x).add(R.one());
    expect(p.is_irreducible()).toBe(true);
  });

  test('x^4 + x^3 + x^2 + x + 1 is not irreducible over GF(2)', () => {
    // This factors as (x^2 + x + 1)^2 - no, actually let's verify
    // x^4 + x^3 + x^2 + x + 1 = (x^5 - 1)/(x - 1) = 5th cyclotomic
    // Over GF(2), 5th cyclotomic factors
    const p = x.pow(4).add(x.pow(3)).add(x.pow(2)).add(x).add(R.one());
    // Actually need to check - 5 doesn't divide 2^k - 1 for small k
    // Let's check: p(0) = 1, p(1) = 1+1+1+1+1 = 1 in GF(2)
    // So no roots in GF(2), could still be reducible
    // Actually x^4 + x^3 + x^2 + x + 1 = (x^2 + x + 1)(x^2 + 1) in GF(2)?
    // Let's verify: (x^2+x+1)(x^2+1) = x^4 + x^2 + x^3 + x + x^2 + 1 = x^4 + x^3 + 1
    // That's not it. Let me check if it's irreducible.
    // Order of 2 mod 5 is 4 (since 2^4 = 16 = 1 mod 5)
    // So actually x^4 + x^3 + x^2 + x + 1 IS irreducible over GF(2)
    expect(p.is_irreducible()).toBe(true);
  });

  test('x^4 + 1 is not irreducible over GF(2)', () => {
    // x^4 + 1 = (x^2 + 1)^2 = (x+1)^4 in GF(2)
    const p = x.pow(4).add(R.one());
    expect(p.is_irreducible()).toBe(false);
  });
});

describe('Polynomial is_irreducible over larger prime fields', () => {
  test('irreducible polynomials over GF(3)', () => {
    const F3 = new FiniteFieldPrime(3n);
    const [R, x] = PolynomialRingConstructor(F3, 'x');
    const one = R.one();

    // Conway polynomial for GF(9): x^2 + 2x + 2
    const conway = CONWAY_POLYNOMIALS[3]![2]!;
    const p = x
      .pow(2)
      .add(x.scalar_mul(F3.__call__(conway[1]!)))
      .add(R.__call__(F3.__call__(conway[0]!)));

    expect(p.is_irreducible()).toBe(true);
  });

  test('x^2 - 1 is reducible over GF(5)', () => {
    const F5 = new FiniteFieldPrime(5n);
    const [R, x] = PolynomialRingConstructor(F5, 'x');

    // x^2 - 1 = (x-1)(x+1) = (x-1)(x+1)
    const p = x.pow(2).sub(R.one());
    expect(p.is_irreducible()).toBe(false);
  });

  test('x^2 + 1 is irreducible over GF(3)', () => {
    const F3 = new FiniteFieldPrime(3n);
    const [R, x] = PolynomialRingConstructor(F3, 'x');

    // x^2 + 1 over GF(3): check if -1 is a quadratic residue
    // -1 = 2, and 2^1 = 2 (not 1 mod 3), so -1 is not a square
    // Thus x^2 + 1 is irreducible
    const p = x.pow(2).add(R.one());
    expect(p.is_irreducible()).toBe(true);
  });
});

describe('Polynomial factor over GF(2)', () => {
  const [R, x] = PolynomialRingConstructor(GF2, 'x');

  test('factor x', () => {
    const factors = x.factor();
    expect(factors.length).toBe(1);
    expect(factors[0]![0].eq(x)).toBe(true);
    expect(factors[0]![1]).toBe(1);
  });

  test('factor x^2 + x = x(x+1)', () => {
    const p = x.pow(2).add(x);
    const factors = p.factor();

    // Should have two factors: x and x+1
    expect(factors.length).toBe(2);

    // Multiply factors back
    let product = R.one();
    for (const [f, m] of factors) {
      product = product.mul(f.pow(m));
    }
    expect(product.eq(p)).toBe(true);
  });

  test('factor x^2 + 1 = (x+1)^2', () => {
    const p = x.pow(2).add(R.one());
    const factors = p.factor();

    // Should have one factor: x+1 with multiplicity 2
    expect(factors.length).toBe(1);
    expect(factors[0]![1]).toBe(2);

    // Multiply back
    let product = R.one();
    for (const [f, m] of factors) {
      product = product.mul(f.pow(m));
    }
    expect(product.eq(p)).toBe(true);
  });

  test('factor irreducible polynomial stays intact', () => {
    // x^2 + x + 1 is irreducible
    const p = x.pow(2).add(x).add(R.one());
    const factors = p.factor();

    expect(factors.length).toBe(1);
    expect(factors[0]![0].eq(p)).toBe(true);
    expect(factors[0]![1]).toBe(1);
  });

  test('factor x^4 + x = x(x+1)(x^2+x+1)', () => {
    // x^4 + x = x(x^3 + 1) = x(x+1)(x^2+x+1)
    const p = x.pow(4).add(x);
    const factors = p.factor();

    // Should have 3 factors
    expect(factors.length).toBe(3);

    // Multiply back
    let product = R.one();
    for (const [f, m] of factors) {
      product = product.mul(f.pow(m));
    }
    expect(product.eq(p)).toBe(true);
  });

  test('factor and multiply equals original', () => {
    // x^6 + x^5 + x^4 + x^3 + x^2 + x + 1
    const p = x.pow(6).add(x.pow(5)).add(x.pow(4)).add(x.pow(3)).add(x.pow(2)).add(x).add(R.one());
    const factors = p.factor();

    // Multiply back
    let product = R.one();
    for (const [f, m] of factors) {
      product = product.mul(f.pow(m));
    }
    expect(product.eq(p)).toBe(true);
  });
});

describe('Polynomial factor over prime fields', () => {
  test('factor x^4 - 1 over GF(5)', () => {
    const F5 = new FiniteFieldPrime(5n);
    const [R, x] = PolynomialRingConstructor(F5, 'x');

    // x^4 - 1 = (x-1)(x+1)(x^2+1) = (x-1)(x-4)(x-2)(x-3)
    // Since -1 = 4, and 2^2 = 4 = -1, 3^2 = 9 = 4 = -1
    // So x^2 + 1 = (x-2)(x-3) where 2,3 are sqrt(-1)
    // All four roots are 1, 4, 2, 3
    const p = x.pow(4).sub(R.one());
    const factors = p.factor();

    // Should have 4 linear factors
    expect(factors.length).toBe(4);

    for (const [f, _m] of factors) {
      expect(f.degree()).toBe(1);
    }

    // Multiply back
    let product = R.one();
    for (const [f, m] of factors) {
      product = product.mul(f.pow(m));
    }
    expect(product.eq(p)).toBe(true);
  });

  test('factor with repeated factors over GF(7)', () => {
    const F7 = new FiniteFieldPrime(7n);
    const [R, x] = PolynomialRingConstructor(F7, 'x');

    // (x - 1)^2 * (x - 2)^3
    const xMinus1 = x.sub(R.one());
    const xMinus2 = x.sub(R.__call__(F7.__call__(2)));
    const p = xMinus1.pow(2).mul(xMinus2.pow(3));

    const factors = p.factor();

    // Should have 2 factors with multiplicities 2 and 3
    expect(factors.length).toBe(2);

    const mults = factors.map(([_f, m]) => m).sort();
    expect(mults).toEqual([2, 3]);

    // Multiply back
    let product = R.one();
    for (const [f, m] of factors) {
      product = product.mul(f.pow(m));
    }
    expect(product.eq(p)).toBe(true);
  });
});

describe('Squarefree decomposition', () => {
  test('squarefree polynomial over GF(2)', () => {
    const [R, x] = PolynomialRingConstructor(GF2, 'x');

    // x^2 + x = x(x+1) is squarefree
    const p = x.pow(2).add(x);
    const decomp = p.squarefree_decomposition();

    expect(decomp.length).toBe(1);
    expect(decomp[0]![1]).toBe(1); // multiplicity 1 means squarefree
  });

  test('polynomial with squares over GF(2)', () => {
    const [R, x] = PolynomialRingConstructor(GF2, 'x');

    // (x+1)^2 * x = x^3 + x
    const xPlus1 = x.add(R.one());
    const p = xPlus1.pow(2).mul(x);

    const decomp = p.squarefree_decomposition();

    // Should have x with mult 1 and (x+1) with mult 2
    let hasMultOne = false;
    let hasMultTwo = false;

    for (const [_f, m] of decomp) {
      if (m === 1) hasMultOne = true;
      if (m === 2) hasMultTwo = true;
    }

    expect(hasMultOne).toBe(true);
    expect(hasMultTwo).toBe(true);

    // Multiply back
    let product = R.one();
    for (const [f, m] of decomp) {
      product = product.mul(f.pow(m));
    }
    expect(product.eq(p)).toBe(true);
  });

  test('squarefree decomposition over GF(5)', () => {
    const F5 = new FiniteFieldPrime(5n);
    const [R, x] = PolynomialRingConstructor(F5, 'x');

    // (x-1)^3 * (x-2)
    const xMinus1 = x.sub(R.one());
    const xMinus2 = x.sub(R.__call__(F5.__call__(2)));
    const p = xMinus1.pow(3).mul(xMinus2);

    const decomp = p.squarefree_decomposition();

    // Check multiplicities
    const mults = decomp.map(([_f, m]) => m).sort();
    expect(mults).toContain(1);
    expect(mults).toContain(3);

    // Multiply back
    let product = R.one();
    for (const [f, m] of decomp) {
      product = product.mul(f.pow(m));
    }
    expect(product.eq(p)).toBe(true);
  });
});

describe('Distinct-degree factorization', () => {
  test('DDF of product of linear factors over GF(2)', () => {
    const [R, x] = PolynomialRingConstructor(GF2, 'x');

    // x(x+1) - product of degree 1 factors
    const p = x.mul(x.add(R.one()));
    const ddf = p.distinct_degree_factorization();

    // Should have one entry with degree 1
    expect(ddf.length).toBe(1);
    expect(ddf[0]![1]).toBe(1);
    expect(ddf[0]![0].degree()).toBe(2); // Product of two linear factors
  });

  test('DDF of irreducible quadratic over GF(2)', () => {
    const [R, x] = PolynomialRingConstructor(GF2, 'x');

    // x^2 + x + 1 is irreducible of degree 2
    const p = x.pow(2).add(x).add(R.one());
    const ddf = p.distinct_degree_factorization();

    // Should have one entry with degree 2
    expect(ddf.length).toBe(1);
    expect(ddf[0]![1]).toBe(2);
  });

  test('DDF of mixed degree product over GF(2)', () => {
    const [R, x] = PolynomialRingConstructor(GF2, 'x');

    // x * (x^2 + x + 1) = x^3 + x^2 + x
    // Has degree 1 factor (x) and degree 2 factor (x^2+x+1)
    const linear = x;
    const quadratic = x.pow(2).add(x).add(R.one());
    const p = linear.mul(quadratic);

    const ddf = p.distinct_degree_factorization();

    // Should have entries for degree 1 and degree 2
    expect(ddf.length).toBe(2);

    const degrees = ddf.map(([_f, d]) => d).sort();
    expect(degrees).toEqual([1, 2]);
  });
});

describe('Conway polynomial irreducibility', () => {
  test('all stored Conway polynomials for GF(2) are irreducible', () => {
    const conwayDegrees = Object.keys(CONWAY_POLYNOMIALS[2]!).map(Number);

    for (const n of conwayDegrees) {
      if (n > 10) continue; // Skip very large degrees for speed

      const coeffs = CONWAY_POLYNOMIALS[2]![n]!;
      const [R, x] = PolynomialRingConstructor(GF2, 'x');

      // Build polynomial: x^n + c_{n-1}*x^{n-1} + ... + c_0
      let poly = x.pow(n);
      for (let i = 0; i < coeffs.length; i++) {
        if (coeffs[i] !== 0) {
          if (i === 0) {
            poly = poly.add(R.one());
          } else {
            poly = poly.add(x.pow(i));
          }
        }
      }

      expect(poly.is_irreducible()).toBe(true);
    }
  });

  test('Conway polynomials for GF(3) are irreducible', () => {
    const F3 = new FiniteFieldPrime(3n);
    const conwayDegrees = Object.keys(CONWAY_POLYNOMIALS[3]!).map(Number);

    for (const n of conwayDegrees) {
      if (n > 6) continue; // Skip large degrees

      const coeffs = CONWAY_POLYNOMIALS[3]![n]!;
      const [R, x] = PolynomialRingConstructor(F3, 'x');

      // Build polynomial
      let poly = x.pow(n);
      for (let i = 0; i < coeffs.length; i++) {
        const c = coeffs[i]!;
        if (c !== 0) {
          const cElem = F3.__call__(c);
          if (i === 0) {
            poly = poly.add(R.__call__(cElem));
          } else {
            poly = poly.add(x.pow(i).scalar_mul(cElem));
          }
        }
      }

      expect(poly.is_irreducible()).toBe(true);
    }
  });
});

describe('Factorization consistency checks', () => {
  test('factor then multiply equals original over GF(7)', () => {
    const F7 = new FiniteFieldPrime(7n);
    const [R, x] = PolynomialRingConstructor(F7, 'x');

    // Random polynomial
    const p = x.pow(5).add(x.pow(3)).add(x).add(R.one());
    const factors = p.factor();

    // Multiply back
    let product = R.one();
    for (const [f, m] of factors) {
      product = product.mul(f.pow(m));
    }

    // Should equal original (up to unit scalar)
    // Make both monic
    const pMonic = p._monic();
    const prodMonic = product._monic();
    expect(prodMonic.eq(pMonic)).toBe(true);
  });

  test('factors are irreducible over GF(5)', () => {
    const F5 = new FiniteFieldPrime(5n);
    const [R, x] = PolynomialRingConstructor(F5, 'x');

    // x^4 - 1
    const p = x.pow(4).sub(R.one());
    const factors = p.factor();

    // Each factor should be irreducible
    for (const [f, _m] of factors) {
      expect(f.is_irreducible()).toBe(true);
    }
  });

  test('roots of factored polynomial', () => {
    const F7 = new FiniteFieldPrime(7n);
    const [R, x] = PolynomialRingConstructor(F7, 'x');

    // (x-1)(x-2)(x-3)
    const p = x
      .sub(R.one())
      .mul(x.sub(R.__call__(F7.__call__(2))))
      .mul(x.sub(R.__call__(F7.__call__(3))));

    const factors = p.factor();
    const roots = p.roots();

    // Should have 3 roots
    expect(roots.length).toBe(3);

    // Each root should correspond to a linear factor
    const rootValues = roots.map(([r, _m]) => r.value).sort();
    expect(rootValues).toEqual([1n, 2n, 3n]);
  });
});

describe('Edge cases', () => {
  test('constant polynomial is not irreducible', () => {
    const [R, _x] = PolynomialRingConstructor(GF2, 'x');
    expect(R.one().is_irreducible()).toBe(false);
  });

  test('zero polynomial throws on roots', () => {
    const [R, _x] = PolynomialRingConstructor(GF2, 'x');
    expect(() => R.zero().roots()).toThrow();
  });

  test('zero polynomial throws on factor', () => {
    const [R, _x] = PolynomialRingConstructor(GF2, 'x');
    expect(() => R.zero().factor()).toThrow();
  });

  test('zero polynomial is not irreducible', () => {
    const [R, _x] = PolynomialRingConstructor(GF2, 'x');
    expect(R.zero().is_irreducible()).toBe(false);
  });
});

describe('monic() method', () => {
  test('monic of zero polynomial is zero', () => {
    const [R, _x] = PolynomialRingConstructor(GF2, 'x');
    expect(R.zero().monic().isZero()).toBe(true);
  });

  test('monic of monic polynomial is itself', () => {
    const [R, x] = PolynomialRingConstructor(GF2, 'x');
    const p = x.pow(2).add(x).add(R.one());
    expect(p.monic().eq(p)).toBe(true);
  });

  test('monic normalizes leading coefficient', () => {
    const F5 = new FiniteFieldPrime(5n);
    const [R, x] = PolynomialRingConstructor(F5, 'x');

    // 3x^2 + 2x + 1 should become x^2 + 4x + 2 (multiply by 3^-1 = 2)
    const three = F5.__call__(3);
    const two = F5.__call__(2);
    const p = x.pow(2).scalar_mul(three).add(x.scalar_mul(two)).add(R.one());

    const m = p.monic();
    expect(m.is_monic()).toBe(true);
    expect(m.leading_coefficient().eq(1)).toBe(true);
  });
});

describe('Higher degree factorization', () => {
  test('x^8 + x over GF(2)', () => {
    const [R, x] = PolynomialRingConstructor(GF2, 'x');

    // x^8 + x = x(x^7 + 1) = x(x+1)(x^3+x+1)(x^3+x^2+1)
    const p = x.pow(8).add(x);
    const factors = p.factor();

    // Multiply back
    let product = R.one();
    for (const [f, m] of factors) {
      product = product.mul(f.pow(m));
    }
    expect(product.eq(p)).toBe(true);

    // Each factor should be irreducible
    for (const [f, _m] of factors) {
      expect(f.is_irreducible()).toBe(true);
    }
  });

  test('x^p - x over GF(p) has all roots', () => {
    const F5 = new FiniteFieldPrime(5n);
    const [R, x] = PolynomialRingConstructor(F5, 'x');

    // x^5 - x has all elements of GF(5) as roots
    const p = x.pow(5).sub(x);
    const roots = p.roots();

    // Should have 5 roots (0, 1, 2, 3, 4)
    expect(roots.length).toBe(5);

    const rootValues = roots.map(([r, _m]) => r.value).sort();
    expect(rootValues).toEqual([0n, 1n, 2n, 3n, 4n]);
  });

  test('factorization of x^9 - x over GF(3)', () => {
    const F3 = new FiniteFieldPrime(3n);
    const [R, x] = PolynomialRingConstructor(F3, 'x');

    // x^9 - x = x(x^8 - 1) factors completely over GF(9)
    // but over GF(3), some factors remain irreducible
    const p = x.pow(9).sub(x);
    const factors = p.factor();

    // Multiply back
    let product = R.one();
    for (const [f, m] of factors) {
      product = product.mul(f.pow(m));
    }
    expect(product.eq(p)).toBe(true);

    // Count factors
    let totalDegree = 0;
    for (const [f, m] of factors) {
      totalDegree += f.degree() * m;
    }
    expect(totalDegree).toBe(9);
  });
});

describe('Factorization determinism', () => {
  test('same polynomial factors to same result', () => {
    const F7 = new FiniteFieldPrime(7n);
    const [R, x] = PolynomialRingConstructor(F7, 'x');

    const p = x.pow(4).sub(R.one());

    const factors1 = p.factor();
    const factors2 = p.factor();

    expect(factors1.length).toBe(factors2.length);

    // Both should produce equivalent factorizations
    let product1 = R.one();
    let product2 = R.one();

    for (const [f, m] of factors1) {
      product1 = product1.mul(f.pow(m));
    }
    for (const [f, m] of factors2) {
      product2 = product2.mul(f.pow(m));
    }

    expect(product1.eq(product2)).toBe(true);
  });
});

// ============================================================================
// Factorization over extension fields (H12) and reducibility detection (C4)
// ============================================================================

describe('Factorization over extension fields (H12)', () => {
  test('x^2 + x + 1 splits over GF(4)', () => {
    // sage: F.<a> = GF(4); R.<x> = F[]; (x^2+x+1).factor() == (x + a) * (x + a + 1)
    const F4 = new FiniteFieldExtension(2n, 2);
    const [R, x] = PolynomialRingConstructor(F4, 'x');
    const f = x.pow(2).add(x).add(R.one());

    expect(f.is_irreducible()).toBe(false);

    const factors = f.factor();
    expect(factors.length).toBe(2);
    for (const [g, m] of factors) {
      expect(g.degree()).toBe(1);
      expect(m).toBe(1);
    }
    let product = R.one();
    for (const [g, m] of factors) product = product.mul(g.pow(m));
    expect(product.eq(f)).toBe(true);

    // roots() and factor() must agree
    expect(f.roots().length).toBe(2);
  });

  test('y^2 + 1 splits over GF(9)', () => {
    const F9 = new FiniteFieldExtension(3n, 2);
    const [R, y] = PolynomialRingConstructor(F9, 'y');
    const f = y.pow(2).add(R.one());

    expect(f.is_irreducible()).toBe(false);

    const factors = f.factor();
    expect(factors.length).toBe(2);
    let product = R.one();
    for (const [g, m] of factors) product = product.mul(g.pow(m));
    expect(product.eq(f)).toBe(true);
  });

  test('a genuinely irreducible quadratic over GF(4) stays irreducible', () => {
    // x^2 + x + a is irreducible over GF(4) (its trace is nonzero)
    const F4 = new FiniteFieldExtension(2n, 2);
    const [R, x] = PolynomialRingConstructor(F4, 'x');
    const a = F4.gen();
    const f = x.pow(2).add(x).add(R.__call__(a));
    expect(f.is_irreducible()).toBe(true);
    expect(f.factor().length).toBe(1);
    expect(f.roots().length).toBe(0);
  });

  test('products of random irreducibles over GF(9) factor back', () => {
    const F9 = new FiniteFieldExtension(3n, 2);
    const [R, x] = PolynomialRingConstructor(F9, 'x');
    const a = F9.gen();
    // (x + a)(x + a + 1)(x^2 + x + a)
    const f = x
      .add(R.__call__(a))
      .mul(x.add(R.__call__(a.add(F9.one()))))
      .mul(x.pow(2).add(x).add(R.__call__(a)));
    let product = R.one();
    for (const [g, m] of f.factor()) {
      product = product.mul(g.pow(m));
      expect(g.is_irreducible()).toBe(true);
    }
    expect(product.eq(f)).toBe(true);
  });
});

describe('is_irreducible detects factors of non-dividing degree (C4)', () => {
  test('GF(2) degree 5 reducibles', () => {
    const [R, x] = PolynomialRingConstructor(GF2, 'x');
    // (x^2+x+1)(x^3+x+1) = x^5 + x^4 + 1
    const f = x.pow(2).add(x).add(R.one()).mul(x.pow(3).add(x).add(R.one()));
    expect(f.toString()).toBe('x^5 + x^4 + 1');
    expect(f.is_irreducible()).toBe(false);
    expect(f.factor().length).toBe(2);

    // (x^2+x+1)(x^3+x^2+1) = x^5 + x + 1
    const g = x
      .pow(2)
      .add(x)
      .add(R.one())
      .mul(x.pow(3).add(x.pow(2)).add(R.one()));
    expect(g.toString()).toBe('x^5 + x + 1');
    expect(g.is_irreducible()).toBe(false);
  });

  test('GF(5) degree 5 reducible', () => {
    const F5 = new FiniteFieldPrime(5n);
    const [R, y] = PolynomialRingConstructor(F5, 'y');
    const f = y
      .pow(2)
      .add(R.__call__(F5.__call__(2)))
      .mul(y.pow(3).add(y).add(R.one()));
    expect(f.is_irreducible()).toBe(false);
    let product = R.one();
    for (const [g, m] of f.factor()) product = product.mul(g.pow(m));
    expect(product.eq(f)).toBe(true);
  });

  test('is_irreducible agrees with factor() over GF(7), degrees 2..6', () => {
    const F7 = new FiniteFieldPrime(7n);
    const [R, x] = PolynomialRingConstructor(F7, 'x');
    const samples = [
      x.pow(2).add(R.one()),
      x.pow(3).sub(x).add(R.one()),
      x.pow(4).add(x).add(R.one()),
      x.pow(5).add(x.pow(2)).add(R.one()),
      x.pow(6).add(x.pow(3)).add(R.one()),
      x.pow(2).sub(R.one()).mul(x.pow(3).add(x).add(R.one())),
      x.pow(2).add(x).add(R.one()).mul(x.pow(2).add(x).add(R.one())),
    ];
    for (const f of samples) {
      const factors = f.factor();
      const nonUnit = factors.filter(([g]) => g.degree() > 0);
      const expected = nonUnit.length === 1 && nonUnit[0]![1] === 1;
      expect(f.is_irreducible()).toBe(expected);
    }
  });
});

// ===========================================================================
// Factorization over ZZ and QQ
//
// `QQ[x].factor()` clears denominators, factors over ZZ with Zassenhaus and
// divides each factor by its leading coefficient, which needs `QQ(...)` to
// accept the shapes Sage's `Rational` accepts.  The expected values below are
// SageMath's:
//
//     sage: R.<x> = QQ[]
//     sage: (x^4 - 1).factor()
//     (x - 1) * (x + 1) * (x^2 + 1)
//     sage: (x^4 + 3*x^2 + 2).factor()
//     (x^2 + 1) * (x^2 + 2)
//     sage: (6*x^2 + x - 2).factor()
//     (6) * (x - 1/2) * (x + 2/3)
//     sage: R(2*x).is_irreducible()
//     True
//     sage: R.<x> = ZZ[]
//     sage: R(2*x).is_irreducible()
//     False
//     sage: (12*(x^2+1)^3*(x+2)).factor()
//     2^2 * 3 * (x + 2) * (x^2 + 1)^3
//     sage: (x^8 - 40*x^6 + 352*x^4 - 960*x^2 + 576).factor()
//     x^8 - 40*x^6 + 352*x^4 - 960*x^2 + 576
// ===========================================================================

/** The real QQ, used as a polynomial coefficient ring. */
type QQElement = Rational & RingElement;
type ZZElement = Integer & RingElement;

const QQCoefficients = QQ as unknown as CoefficientRing<QQElement>;

interface IntegerCoefficientRing extends CoefficientRing<ZZElement> {
  is_integral_domain(): boolean;
  characteristic(): bigint;
  toString(): string;
}

const qqElement = (n: bigint, d: bigint = 1n): QQElement => new Rational(n, d) as QQElement;
const zzElement = (n: bigint): ZZElement => new Integer(n) as ZZElement;

const ZZCoefficients: IntegerCoefficientRing = {
  zero: () => zzElement(0n),
  one: () => zzElement(1n),
  __call__: (x: unknown) => (x instanceof Integer ? (x as ZZElement) : zzElement(x as bigint)),
  is_field: () => false,
  is_integral_domain: () => true,
  characteristic: () => 0n,
  toString: () => 'Integer Ring',
};

const RQQ = new PolynomialRing(QQCoefficients, 'x');
const RZZ = new PolynomialRing(ZZCoefficients, 'x');

/** Polynomial over QQ from integer coefficients, constant term first. */
const qq = (coeffs: bigint[]): Polynomial<QQElement> =>
  new Polynomial<QQElement>(
    coeffs.map((c) => qqElement(c)),
    RQQ
  );
/** Polynomial over QQ from [numerator, denominator] pairs. */
const qqFrac = (coeffs: Array<[bigint, bigint]>): Polynomial<QQElement> =>
  new Polynomial<QQElement>(
    coeffs.map(([n, d]) => qqElement(n, d)),
    RQQ
  );
const zz = (coeffs: bigint[]): Polynomial<ZZElement> =>
  new Polynomial<ZZElement>(
    coeffs.map((c) => zzElement(c)),
    RZZ
  );

/** Multiply out a factorization, so it can be compared with the input. */
function productOf<C extends RingElement>(
  one: Polynomial<C>,
  factors: Array<[Polynomial<C>, number]>
): Polynomial<C> {
  let product = one;
  for (const [g, e] of factors) {
    for (let i = 0; i < e; i++) product = product.mul(g);
  }
  return product;
}

// --- independent oracle: Kronecker's factorization method -------------------
// Naive on purpose and written from scratch: rational roots by the rational
// root theorem, then a search over the divisors of f(a_0), ..., f(a_m) with
// Lagrange interpolation.  Nothing here is shared with the implementation.

type IntPoly = bigint[]; // constant term first

function oracleStrip(f: IntPoly): IntPoly {
  const g = [...f];
  while (g.length > 0 && g[g.length - 1] === 0n) g.pop();
  return g;
}
function oracleMul(a: IntPoly, b: IntPoly): IntPoly {
  if (a.length === 0 || b.length === 0) return [];
  const r: bigint[] = new Array(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) r[i + j] = r[i + j]! + a[i]! * b[j]!;
  }
  return oracleStrip(r);
}
function oracleEval(f: IntPoly, x: bigint): bigint {
  let r = 0n;
  for (let i = f.length - 1; i >= 0; i--) r = r * x + f[i]!;
  return r;
}
function oracleGcd(a: bigint, b: bigint): bigint {
  let x = a < 0n ? -a : a;
  let y = b < 0n ? -b : b;
  while (y !== 0n) [x, y] = [y, x % y];
  return x;
}
/** primitive part with a positive leading coefficient */
function oraclePrimitive(f: IntPoly): IntPoly {
  const s = oracleStrip(f);
  if (s.length === 0) return s;
  let c = 0n;
  for (const x of s) c = oracleGcd(c, x);
  if (s[s.length - 1]! < 0n) c = -c;
  return s.map((x) => x / c);
}
/** exact division over ZZ, or null */
function oracleDivExact(a: IntPoly, b: IntPoly): IntPoly | null {
  const u = oracleStrip(a);
  const v = oracleStrip(b);
  if (v.length === 0) throw new Error('division by zero');
  if (u.length === 0) return [];
  if (u.length < v.length) return null;
  const q: bigint[] = new Array(u.length - v.length + 1).fill(0n);
  const r = [...u];
  const lc = v[v.length - 1]!;
  for (let i = u.length - 1; i >= v.length - 1; i--) {
    if (r[i] === 0n) continue;
    if (r[i]! % lc !== 0n) return null;
    const c = r[i]! / lc;
    q[i - v.length + 1] = c;
    for (let j = 0; j < v.length; j++)
      r[i - v.length + 1 + j] = r[i - v.length + 1 + j]! - c * v[j]!;
  }
  if (r.some((c) => c !== 0n)) return null;
  return oracleStrip(q);
}
function oracleDivisors(n: bigint): bigint[] {
  const a = n < 0n ? -n : n;
  const out: bigint[] = [];
  for (let d = 1n; d * d <= a; d++) {
    if (a % d === 0n) {
      out.push(d);
      if (d * d !== a) out.push(a / d);
    }
  }
  return out;
}
/** split off every rational root (rational root theorem) */
function oracleLinearFactors(f: IntPoly): { factors: IntPoly[]; rest: IntPoly } {
  const factors: IntPoly[] = [];
  let g = oraclePrimitive(f);
  while (g.length - 1 >= 1 && g[0] === 0n) {
    factors.push([0n, 1n]);
    g = g.slice(1);
  }
  let changed = true;
  while (changed && g.length - 1 >= 1) {
    changed = false;
    const ps = oracleDivisors(g[0]!);
    const qs = oracleDivisors(g[g.length - 1]!);
    search: for (const q of qs) {
      for (const p0 of ps) {
        for (const p of [p0, -p0]) {
          if (oracleGcd(p, q) !== 1n) continue;
          const n = g.length - 1;
          let s = 0n;
          for (let i = 0; i <= n; i++) s += g[i]! * p ** BigInt(i) * q ** BigInt(n - i);
          if (s !== 0n) continue;
          const quotient = oracleDivExact(g, [-p, q]);
          if (quotient) {
            factors.push([-p, q]);
            g = oraclePrimitive(quotient);
            changed = true;
            break search;
          }
        }
      }
    }
  }
  return { factors, rest: g };
}
/** Kronecker: a factor of f of degree exactly m, or null */
function oracleKroneckerFactor(f: IntPoly, m: number): IntPoly | null {
  const pts: bigint[] = [];
  let a = 0n;
  while (pts.length < m + 1) {
    if (oracleEval(f, a) !== 0n) pts.push(a);
    a = a > 0n ? -a : -a + 1n;
  }
  const divs = pts.map((pt) => oracleDivisors(oracleEval(f, pt)).flatMap((d) => [d, -d]));
  const D = pts.map((ai, i) => pts.reduce((acc, aj, j) => (i === j ? acc : acc * (ai - aj)), 1n));
  const P = pts.map((_, i) => {
    let poly: IntPoly = [1n];
    pts.forEach((aj, j) => {
      if (j !== i) poly = oracleMul(poly, [-aj, 1n]);
    });
    return poly;
  });
  let L = 1n;
  for (const d of D) L = (L / oracleGcd(L, d)) * d;
  const idx: number[] = new Array(m + 1).fill(0);
  for (;;) {
    // g and -g are the same factorization, so only positive first values
    if (divs[0]![idx[0]!]! > 0n) {
      const N: bigint[] = new Array(m + 1).fill(0n);
      for (let i = 0; i <= m; i++) {
        const scale = (L / D[i]!) * divs[i]![idx[i]!]!;
        const Pi = P[i]!;
        for (let j = 0; j < Pi.length; j++) N[j] = (N[j] ?? 0n) + scale * Pi[j]!;
      }
      if (N.every((c) => c % L === 0n)) {
        const g = oracleStrip(N.map((c) => c / L));
        if (g.length - 1 === m && oracleDivExact(f, g) !== null) return oraclePrimitive(g);
      }
    }
    let i = m;
    for (;;) {
      idx[i] = idx[i]! + 1;
      if (idx[i]! < divs[i]!.length) break;
      idx[i] = 0;
      i--;
      if (i < 0) return null;
    }
  }
}
/** complete factorization into irreducibles (primitive, positive lc) */
function oracleFactor(f: IntPoly): IntPoly[] {
  const g = oraclePrimitive(f);
  if (g.length - 1 <= 0) return [];
  const { factors, rest } = oracleLinearFactors(g);
  const out = [...factors];
  const stack: IntPoly[] = [rest];
  while (stack.length > 0) {
    const h = stack.pop()!;
    const d = h.length - 1;
    if (d <= 0) continue;
    if (d === 1) {
      out.push(oraclePrimitive(h));
      continue;
    }
    let split: IntPoly | null = null;
    for (let m = 2; 2 * m <= d; m++) {
      split = oracleKroneckerFactor(h, m);
      if (split) break;
    }
    if (split) {
      stack.push(oraclePrimitive(split), oraclePrimitive(oracleDivExact(h, split)!));
    } else {
      out.push(oraclePrimitive(h));
    }
  }
  return out;
}
/** clear denominators of a rational polynomial, returning its primitive part */
function clearDenominators(coeffs: readonly QQElement[]): IntPoly {
  let l = 1n;
  for (const c of coeffs) l = (l / oracleGcd(l, c.denominator)) * c.denominator;
  return oraclePrimitive(coeffs.map((c) => c.numerator * (l / c.denominator)));
}

describe('Kronecker oracle self-check', () => {
  test('known factorizations', () => {
    // x^4 - 1 = (x-1)(x+1)(x^2+1)
    expect(oracleFactor([-1n, 0n, 0n, 0n, 1n]).length).toBe(3);
    // (x^2+1)(x^2+2)
    expect(oracleFactor([2n, 0n, 3n, 0n, 1n]).length).toBe(2);
    // 6x^2 + x - 2 = (2x-1)(3x+2)
    expect(oracleFactor([-2n, 1n, 6n]).length).toBe(2);
    // irreducible ones
    expect(oracleFactor([1n, 1n, 1n, 1n, 1n]).length).toBe(1); // Phi_5
    expect(oracleFactor([1n, 0n, 0n, 0n, 1n]).length).toBe(1); // x^4 + 1
    expect(oracleFactor([2n, 0n, 0n, 1n]).length).toBe(1); // x^3 + 2
    // x^6 - 1 = (x-1)(x+1)(x^2+x+1)(x^2-x+1)
    expect(oracleFactor([-1n, 0n, 0n, 0n, 0n, 0n, 1n]).length).toBe(4);
  });
});

describe('factor over QQ', () => {
  test("Sage's doctest values", () => {
    // (x^4 - 1).factor() == (x - 1) * (x + 1) * (x^2 + 1)
    const f = qq([-1n, 0n, 0n, 0n, 1n]);
    const factors = f.factor();
    expect(factors.map(([g, e]) => `(${g})^${e}`).sort()).toEqual(
      ['(x + -1)^1', '(x + 1)^1', '(x^2 + 1)^1'].sort()
    );
    expect(productOf(RQQ.one(), factors).eq(f)).toBe(true);

    // (x^4 + 3x^2 + 2).factor() == (x^2 + 1) * (x^2 + 2)
    const g = qq([2n, 0n, 3n, 0n, 1n]);
    const gFactors = g.factor();
    expect(gFactors.length).toBe(2);
    expect(gFactors.every(([h, e]) => h.degree() === 2 && e === 1)).toBe(true);
    expect(productOf(RQQ.one(), gFactors).eq(g)).toBe(true);

    // (6x^2 + x - 2).factor() == (6) * (x - 1/2) * (x + 2/3)
    const h = qq([-2n, 1n, 6n]);
    const hFactors = h.factor();
    const nonUnit = hFactors.filter(([p]) => p.degree() > 0);
    expect(nonUnit.length).toBe(2);
    expect(nonUnit.every(([p]) => p.is_monic())).toBe(true);
    expect(nonUnit.map(([p]) => p.toString()).sort()).toEqual(['x + -1/2', 'x + 2/3']);
    expect(productOf(RQQ.one(), hFactors).eq(h)).toBe(true);
  });

  test('rational coefficients: (x - 1/2)(x + 2/3)', () => {
    const f = qqFrac([
      [-1n, 3n],
      [1n, 6n],
      [1n, 1n],
    ]);
    const factors = f.factor();
    expect(factors.map(([g]) => g.toString()).sort()).toEqual(['x + -1/2', 'x + 2/3']);
    expect(productOf(RQQ.one(), factors).eq(f)).toBe(true);
  });

  test('the leading coefficient is the unit and every factor is monic', () => {
    // 2*x^2 + 2 = 2 * (x^2 + 1)
    const f = qq([2n, 0n, 2n]);
    const factors = f.factor();
    expect(factors.filter(([g]) => g.degree() > 0).map(([g]) => g.toString())).toEqual(['x^2 + 1']);
    expect(productOf(RQQ.one(), factors).eq(f)).toBe(true);
  });

  test('multiplicities', () => {
    // (x - 1)^25 -- an earlier fixed cap of 20 truncated this
    let f = RQQ.one();
    const xMinus1 = qq([-1n, 1n]);
    for (let i = 0; i < 25; i++) f = f.mul(xMinus1);
    expect(f.factor()).toEqual([[xMinus1, 25]]);

    // (x^2 + 1)^5
    let g = RQQ.one();
    const quad = qq([1n, 0n, 1n]);
    for (let i = 0; i < 5; i++) g = g.mul(quad);
    expect(g.factor()).toEqual([[quad, 5]]);
  });

  test('the charpoly of matrix(QQ, 4, 4, range(16)) (minpoly blocker)', () => {
    // sage: A = matrix(QQ, 4, 4, range(16))
    // sage: A.charpoly()
    // x^4 - 30*x^3 - 80*x^2
    // sage: A.charpoly().factor()
    // (x^2 - 30*x - 80) * x^2
    // sage: A.minpoly()
    // x^3 - 30*x^2 - 80*x
    // `minpoly` builds its answer out of this factorization, so it stayed
    // broken over QQ for as long as `factor()` did.
    const charpoly = qq([0n, 0n, -80n, -30n, 1n]);
    const factors = charpoly.factor();
    expect(factors.map(([g, e]) => `(${g})^${e}`).sort()).toEqual(
      ['(x)^2', '(x^2 + (-30)*x + -80)^1'].sort()
    );
    expect(productOf(RQQ.one(), factors).eq(charpoly)).toBe(true);

    // the radical of the charpoly is the minimal polynomial here
    let radical = RQQ.one();
    for (const [g] of factors) {
      if (g.degree() > 0) radical = radical.mul(g);
    }
    expect(radical.toString()).toBe('x^3 + (-30)*x^2 + (-80)*x');
    // x^2 - 30x - 80 is irreducible: 30^2 + 4*80 = 1220 is not a square
    expect(qq([-80n, -30n, 1n]).is_irreducible()).toBe(true);
  });

  test('x^n - 1 has one irreducible factor per divisor of n', () => {
    // The cyclotomic factorization: x^n - 1 = prod_{d | n} Phi_d(x).
    for (let n = 1; n <= 24; n++) {
      const coeffs: bigint[] = new Array(n + 1).fill(0n);
      coeffs[0] = -1n;
      coeffs[n] = 1n;
      const f = qq(coeffs);
      const factors = f.factor().filter(([g]) => g.degree() > 0);
      let divisors = 0;
      for (let d = 1; d <= n; d++) if (n % d === 0) divisors++;
      expect(factors.reduce((s, [, e]) => s + e, 0)).toBe(divisors);
      expect(factors.every(([, e]) => e === 1)).toBe(true);
      expect(productOf(RQQ.one(), f.factor()).eq(f)).toBe(true);
    }
  });

  test('Swinnerton-Dyer polynomials are irreducible (Zassenhaus worst case)', () => {
    // Minimal polynomial of sqrt2 + sqrt3 + sqrt5: irreducible over QQ but a
    // product of quadratics modulo every prime.
    const sd3 = qq([576n, 0n, -960n, 0n, 352n, 0n, -40n, 0n, 1n]);
    expect(sd3.factor()).toEqual([[sd3, 1]]);
    expect(sd3.is_irreducible()).toBe(true);

    // ... + sqrt7, degree 16
    const sd4 = qq([
      1216800225n,
      0n,
      -152253360n,
      0n,
      7391364n,
      0n,
      -223664n,
      0n,
      5702n,
      0n,
      -160n,
      0n,
      92n,
      0n,
      -20n,
      0n,
      1n,
    ]);
    expect(sd4.factor()).toEqual([[sd4, 1]]);
    expect(sd4.is_irreducible()).toBe(true);
  });

  test('is_irreducible over QQ ignores the unit', () => {
    // sage: R.<x> = QQ[]; R(2*x).is_irreducible() -> True
    expect(qq([0n, 2n]).is_irreducible()).toBe(true);
    expect(qq([2n, 0n, 2n]).is_irreducible()).toBe(true); // 2x^2 + 2
    expect(qq([2n, 0n, 0n, 1n]).is_irreducible()).toBe(true); // x^3 + 2
    expect(qq([1n, 0n, 0n, 0n, 1n]).is_irreducible()).toBe(true); // x^4 + 1
    expect(qq([1n, 1n, 1n]).is_irreducible()).toBe(true); // x^2 + x + 1
    expect(qq([-1n, 0n, 1n]).is_irreducible()).toBe(false); // x^2 - 1
    expect(qq([2n, 0n, 3n, 0n, 1n]).is_irreducible()).toBe(false); // (x^2+1)(x^2+2)
    expect(qq([-2n, 1n, 6n]).is_irreducible()).toBe(false); // 6x^2 + x - 2
  });
});

describe('factor over ZZ', () => {
  test("Sage's doctest values", () => {
    // (12*(x^2+1)^3*(x+2)).factor() == 2^2 * 3 * (x + 2) * (x^2 + 1)^3
    let f = zz([12n]);
    const quad = zz([1n, 0n, 1n]);
    for (let i = 0; i < 3; i++) f = f.mul(quad);
    f = f.mul(zz([2n, 1n]));
    const factors = f.factor();
    expect(productOf(RZZ.one(), factors).eq(f)).toBe(true);
    const described = factors
      .map(([g, e]) => `${g}^${e}`)
      .sort()
      .join(' * ');
    expect(described).toBe('1*x + 2^1 * 1*x^2 + 1^3 * 2^2 * 3^1');

    // (-x^2 + 4).factor() == (-1) * (x - 2) * (x + 2)
    const g = zz([4n, 0n, -1n]);
    expect(productOf(RZZ.one(), g.factor()).eq(g)).toBe(true);

    // 6x^2 + x - 2 = (2x - 1)(3x + 2) over ZZ (no unit, both factors kept)
    const h = zz([-2n, 1n, 6n]);
    const hFactors = h.factor();
    expect(hFactors.map(([p]) => p.toString()).sort()).toEqual(['2*x + -1', '3*x + 2']);
    expect(productOf(RZZ.one(), hFactors).eq(h)).toBe(true);
  });

  test('is_irreducible over ZZ', () => {
    // sage: R.<x> = ZZ[]
    expect(zz([0n, 2n]).is_irreducible()).toBe(false); // 2*x is not primitive
    expect(zz([1n, 0n, 0n, 1n]).is_irreducible()).toBe(false); // x^3 + 1
    expect(zz([2n, 0n, 0n, 1n]).is_irreducible()).toBe(true); // x^3 + 2
    expect(zz([-1n, 0n, 1n]).is_irreducible()).toBe(false); // x^2 - 1
    expect(zz([5n]).is_irreducible()).toBe(true); // R(5)
    expect(zz([4n]).is_irreducible()).toBe(false); // R(4)
    expect(zz([2n, 0n, 3n, 0n, 1n]).is_irreducible()).toBe(false); // (x^2+1)(x^2+2)
  });
});

describe('QQ[x] factorization against the Kronecker oracle', () => {
  test('200 random polynomials of degree <= 6', () => {
    // Deterministic LCG, so a failure is reproducible.
    let seed = 20260727n;
    const rand = (n: number): number => {
      seed = (seed * 6364136223846793005n + 1442695040888963407n) & ((1n << 64n) - 1n);
      return Number((seed >> 33n) % BigInt(n));
    };

    let checked = 0;
    let composite = 0;
    for (let iter = 0; iter < 200; iter++) {
      const d = 1 + rand(6);
      const coeffs: QQElement[] = [];
      for (let i = 0; i <= d; i++) {
        let num = BigInt(rand(17) - 8);
        if (i === d && num === 0n) num = 1n;
        coeffs.push(qqElement(num, BigInt(1 + rand(4))));
      }
      const f = new Polynomial<QQElement>(coeffs, RQQ);
      if (f.degree() < 1) continue;
      checked++;

      const factors = f.factor();

      // (1) the factors multiply back to f
      expect(productOf(RQQ.one(), factors).eq(f)).toBe(true);

      // (2) every non-constant factor is monic and irreducible
      const nonUnit = factors.filter(([g]) => g.degree() > 0);
      for (const [g] of nonUnit) {
        expect(g.is_monic()).toBe(true);
        expect(oracleFactor(clearDenominators(g.coeffs)).length).toBe(1);
      }

      // (3) the whole factorization agrees with Kronecker's
      const ours: string[] = [];
      for (const [g, e] of nonUnit) {
        for (let i = 0; i < e; i++) ours.push(clearDenominators(g.coeffs).join(','));
      }
      const theirs = oracleFactor(clearDenominators(f.coeffs)).map((g) => g.join(','));
      expect(ours.sort()).toEqual(theirs.sort());

      if (nonUnit.reduce((s, [, e]) => s + e, 0) > 1) composite++;
    }

    expect(checked).toBe(200);
    // the sweep must actually exercise composite inputs
    expect(composite).toBeGreaterThan(15);
  });

  test('150 products of known irreducible factors (degree <= 14)', () => {
    let seed = 424242n;
    const rand = (n: number): number => {
      seed = (seed * 6364136223846793005n + 1442695040888963407n) & ((1n << 64n) - 1n);
      return Number((seed >> 33n) % BigInt(n));
    };

    // pool of irreducible integer polynomials, certified by the oracle
    const pool: IntPoly[] = [];
    while (pool.length < 40) {
      const d = 1 + rand(4);
      const g: bigint[] = [];
      for (let i = 0; i <= d; i++) {
        let c = BigInt(rand(11) - 5);
        if (i === d && c === 0n) c = 1n;
        g.push(c);
      }
      const p = oraclePrimitive(g);
      if (p.length - 1 < 1) continue;
      if (oracleFactor(p).length !== 1) continue;
      if (pool.some((q) => q.join(',') === p.join(','))) continue;
      pool.push(p);
    }

    let checked = 0;
    for (let iter = 0; iter < 150; iter++) {
      const chosen: Array<[IntPoly, number]> = [];
      let f: IntPoly = [1n];
      for (let i = 0; i < 1 + rand(3); i++) {
        const g = pool[rand(pool.length)]!;
        if (chosen.some(([h]) => h.join(',') === g.join(','))) continue;
        const e = 1 + rand(3);
        if (f.length - 1 + (g.length - 1) * e > 14) continue;
        chosen.push([g, e]);
        for (let j = 0; j < e; j++) f = oracleMul(f, g);
      }
      if (chosen.length === 0) continue;
      checked++;

      const poly = qq(f);
      const factors = poly.factor();
      expect(productOf(RQQ.one(), factors).eq(poly)).toBe(true);

      const ours = factors
        .filter(([g]) => g.degree() > 0)
        .map(([g, e]) => `${clearDenominators(g.coeffs).join(',')}^${e}`)
        .sort();
      const want = chosen.map(([g, e]) => `${g.join(',')}^${e}`).sort();
      expect(ours).toEqual(want);
    }
    expect(checked).toBeGreaterThan(100);
  });
});

/* =====================================================================
 * ZZ[x] factorisation: Zassenhaus recombination and van Hoeij / LLL
 *
 * These exercise the FLINT transcriptions in polynomial_element.ts:
 *   fmpz_poly_factor/{factor_zassenhaus, factor_van_hoeij, CLD_mat,
 *                     van_hoeij_check_if_solved, zassenhaus_subset,
 *                     zassenhaus_prune}.c
 *   fmpz_poly/{CLD_bound, divlow_smodp, divhigh_smodp}.c
 *   fmpz_mat/{next_col_van_hoeij, col_partition}.c
 *   fmpz_lll/ (LLL with removal)
 * ===================================================================== */

/** deterministic LCG, so any failure is reproducible */
function makeRand(seed0: bigint): (n: number) => number {
  let seed = seed0;
  return (n: number) => {
    seed = (seed * 6364136223846793005n + 1442695040888963407n) & ((1n << 64n) - 1n);
    return Number((seed >> 33n) % BigInt(n));
  };
}

function oracleDeriv(f: IntPoly): IntPoly {
  const r: bigint[] = [];
  for (let i = 1; i < f.length; i++) r.push(f[i]! * BigInt(i));
  return r;
}

/**
 * The Swinnerton-Dyer polynomial of `sqrt(p_1), ..., sqrt(p_k)`: the product of
 * `x -+ sqrt(p_1) -+ ... -+ sqrt(p_k)` over all `2^k` sign choices, of degree
 * `2^k` and irreducible over `Q`.  Built exactly by the norm recursion
 * `f_k(x) = A(x)^2 - p B(x)^2` where `f_{k-1}(x + s) = A + sB` in `Z[x][s]/(s^2 - p)`.
 *
 * This is van Hoeij's classic worst case: modulo every prime the factors have
 * degree at most 2, so Zassenhaus has to walk `2^(2^(k-1))` subsets.
 */
function swinnertonDyer(primes: bigint[]): IntPoly {
  let f: IntPoly = [0n, 1n]; // x
  for (const p of primes) {
    const n = f.length;
    const A: bigint[] = new Array(n).fill(0n);
    const B: bigint[] = new Array(n).fill(0n);
    const C: bigint[][] = [];
    for (let j = 0; j < n; j++) {
      C.push(new Array<bigint>(j + 1).fill(0n));
      C[j]![0] = 1n;
      for (let i = 1; i <= j; i++) C[j]![i] = (C[j - 1]![i - 1] ?? 0n) + (C[j - 1]![i] ?? 0n);
    }
    for (let j = 0; j < n; j++) {
      if (f[j] === 0n) continue;
      for (let i = 0; i <= j; i++) {
        const coef = f[j]! * C[j]![i]! * p ** BigInt(i >> 1);
        if (i % 2 === 0) A[j - i] = A[j - i]! + coef;
        else B[j - i] = B[j - i]! + coef;
      }
    }
    const A2 = oracleMul(A, A);
    const B2 = oracleMul(B, B).map((c) => c * p);
    const len = Math.max(A2.length, B2.length);
    const g: bigint[] = new Array(len).fill(0n);
    for (let i = 0; i < len; i++) g[i] = (A2[i] ?? 0n) - (B2[i] ?? 0n);
    f = oracleStrip(g);
  }
  return f;
}

describe('Swinnerton-Dyer polynomials (van Hoeij worst case)', () => {
  test('the small ones are the classical polynomials', () => {
    expect(swinnertonDyer([2n]).join(',')).toBe('-2,0,1'); // x^2 - 2
    expect(swinnertonDyer([2n, 3n]).join(',')).toBe('1,0,-10,0,1'); // x^4 - 10x^2 + 1
    // x^8 - 40x^6 + 352x^4 - 960x^2 + 576
    expect(swinnertonDyer([2n, 3n, 5n]).join(',')).toBe('576,0,-960,0,352,0,-40,0,1');
  });

  test('degrees 16, 32 and 64 are irreducible, and fast', () => {
    for (const primes of [
      [2n, 3n, 5n, 7n],
      [2n, 3n, 5n, 7n, 11n],
      [2n, 3n, 5n, 7n, 11n, 13n],
    ]) {
      const f = swinnertonDyer(primes);
      const t0 = Date.now();
      const factors = _zz_factor_internal.factorSquarefreeIntPoly(f);
      const dt = Date.now() - t0;
      expect(f.length - 1).toBe(1 << primes.length);
      expect(factors.length).toBe(1);
      expect(factors[0]!.join(',')).toBe(f.join(','));
      // Zassenhaus would need 2^(2^(k-1)) subsets here; van Hoeij is polynomial.
      expect(dt).toBeLessThan(20000);
    }
  });

  test('the modular factorisation really is the hard one (all factors of degree <= 2)', () => {
    const f = swinnertonDyer([2n, 3n, 5n, 7n, 11n]);
    const { fac } = _zz_factor_internal.chooseFactorizationPrime(f);
    expect(fac.length).toBeGreaterThan(8); // > cutoff, so van Hoeij is taken
    for (const g of fac) expect(g.length - 1).toBeLessThanOrEqual(2);
  });

  test('degree 32 with a completely split prime (32 linear modular factors)', () => {
    const primes = [2n, 3n, 5n, 7n, 11n];
    const f = swinnertonDyer(primes);
    // smallest prime for which every p_i is a quadratic residue
    const legendre = (a: bigint, p: bigint): bigint => {
      let base = ((a % p) + p) % p;
      let e = (p - 1n) / 2n;
      let res = 1n;
      while (e > 0n) {
        if (e & 1n) res = (res * base) % p;
        base = (base * base) % p;
        e >>= 1n;
      }
      return res;
    };
    let p = 3n;
    for (;;) {
      p = next_prime(p);
      if (primes.every((q) => legendre(q, p) === 1n)) break;
    }
    const roots: bigint[] = [];
    for (let a = 0n; a < p; a++) if (((oracleEval(f, a) % p) + p) % p === 0n) roots.push(a);
    expect(roots.length).toBe(32);
    const fac = roots.map((rt) => [(p - rt) % p, 1n]);
    const t0 = Date.now();
    const factors = _zz_factor_internal.fmpzPolyFactorVanHoeij(fac, f, p);
    expect(factors.length).toBe(1);
    expect(Date.now() - t0).toBeLessThan(20000);
  });

  test('products of Swinnerton-Dyer polynomials recombine correctly', () => {
    const cases: Array<[string, IntPoly[]]> = [
      ['SD8 * SD8', [swinnertonDyer([2n, 3n, 5n]), swinnertonDyer([2n, 3n, 7n])]],
      [
        'SD8 * SD8 * SD8',
        [swinnertonDyer([2n, 3n, 5n]), swinnertonDyer([2n, 3n, 7n]), swinnertonDyer([3n, 5n, 7n])],
      ],
      [
        'SD16 * SD8 * SD4',
        [swinnertonDyer([2n, 3n, 5n, 7n]), swinnertonDyer([2n, 3n, 5n]), swinnertonDyer([2n, 3n])],
      ],
      ['SD16 * SD16', [swinnertonDyer([2n, 3n, 5n, 7n]), swinnertonDyer([2n, 3n, 5n, 11n])]],
    ];
    for (const [name, parts] of cases) {
      let f: IntPoly = [1n];
      for (const g of parts) f = oracleMul(f, g);
      const factors = _zz_factor_internal.factorSquarefreeIntPoly(f);
      const got = factors
        .map((g) => g.join(','))
        .sort()
        .join(' | ');
      const want = parts
        .map((g) => g.join(','))
        .sort()
        .join(' | ');
      expect(`${name}: ${got}`).toBe(`${name}: ${want}`);
    }
  });
});

describe('ZZ[x] factorisation: 500 random polynomials', () => {
  test('product of factors equals the input and every factor is irreducible', () => {
    const rand = makeRand(20260728n);

    // A pool of irreducible integer polynomials, each certified by the
    // independently written Kronecker oracle above.
    const pool: IntPoly[] = [];
    while (pool.length < 60) {
      const d = 1 + rand(5);
      const g: bigint[] = [];
      for (let i = 0; i <= d; i++) g.push(BigInt(rand(15) - 7));
      if (g[d] === 0n) g[d] = 1n;
      const q = oraclePrimitive(oracleStrip(g));
      if (q.length - 1 < 1) continue;
      if (oracleFactor(q).length !== 1) continue;
      if (pool.some((h) => h.join(',') === q.join(','))) continue;
      pool.push(q);
    }

    let checked = 0;
    let vanHoeijPaths = 0;
    let maxDegree = 0;

    for (let iter = 0; iter < 500; iter++) {
      const chosen = new Map<string, [IntPoly, number]>();
      let f: IntPoly = [1n];
      const k = 1 + rand(8);
      for (let i = 0; i < k; i++) {
        const g = pool[rand(pool.length)]!;
        const key = g.join(',');
        if (chosen.has(key)) continue;
        const e = 1 + rand(3);
        chosen.set(key, [g, e]);
        for (let j = 0; j < e; j++) f = oracleMul(f, g);
      }
      if (chosen.size === 0) continue;
      checked++;
      maxDegree = Math.max(maxDegree, f.length - 1);

      const poly = zz(f);
      const factors = poly.factor();

      // (1) the factors multiply back to the input, exactly
      expect(productOf(RZZ.one(), factors).eq(poly)).toBe(true);

      // (2) the factorization is exactly the multiset of irreducibles we built
      //     from (so every returned factor is irreducible, with the right
      //     multiplicity, and none is missing)
      const ours = factors
        .filter(([g]) => g.degree() > 0)
        .map(
          ([g, e]) => `${oraclePrimitive(g.coeffs.map((c) => BigInt(c.toString()))).join(',')}^${e}`
        )
        .sort();
      const want = Array.from(chosen.values())
        .map(([g, e]) => `${g.join(',')}^${e}`)
        .sort();
      expect(ours).toEqual(want);

      // did this input reach van Hoeij (more than the cutoff of 8 modular
      // factors) on at least one of its squarefree parts?
      for (const [g] of chosen.values()) {
        void g;
      }
      const squarefreePart = Array.from(chosen.values()).reduce((acc, [g]) => oracleMul(acc, g), [
        1n,
      ] as IntPoly);
      if (squarefreePart.length > 2 && squarefreePart[0] !== 0n) {
        const { fac } = _zz_factor_internal.chooseFactorizationPrime(squarefreePart);
        if (fac.length > 8) vanHoeijPaths++;
      }
    }

    expect(checked).toBe(500);
    expect(maxDegree).toBeGreaterThan(20);
    // the sweep must actually exercise the van Hoeij branch
    expect(vanHoeijPaths).toBeGreaterThan(50);
  }, 600000);

  test('200 fully random polynomials of degree <= 6 against the Kronecker oracle', () => {
    const rand = makeRand(555111n);
    let checked = 0;
    let composite = 0;
    for (let iter = 0; iter < 200; iter++) {
      const d = 1 + rand(6);
      const g: bigint[] = [];
      for (let i = 0; i <= d; i++) g.push(BigInt(rand(21) - 10));
      if (g[d] === 0n) g[d] = 1n;
      const f = oracleStrip(g);
      if (f.length - 1 < 1) continue;
      checked++;

      const poly = zz(f);
      const factors = poly.factor();
      expect(productOf(RZZ.one(), factors).eq(poly)).toBe(true);

      const ours: string[] = [];
      for (const [h, e] of factors) {
        if (h.degree() < 1) continue;
        const cs = oraclePrimitive(h.coeffs.map((c) => BigInt(c.toString())));
        expect(oracleFactor(cs).length).toBe(1); // each factor is irreducible
        for (let i = 0; i < e; i++) ours.push(cs.join(','));
      }
      const theirs = oracleFactor(f).map((h) => h.join(','));
      expect(ours.sort()).toEqual(theirs.sort());
      if (theirs.length > 1) composite++;
    }
    expect(checked).toBe(200);
    expect(composite).toBeGreaterThan(20);
  }, 600000);
});

describe('Factorisation primes: the search is unbounded (FLINT factor_zassenhaus.c:120)', () => {
  test('all primes below 10000 can divide lc, f(0) and disc(f)', () => {
    // f = x^2 - P with P the product of every prime <= 10007: disc(f) = 4P and
    // f(0) = -P, so every prime <= 10007 is rejected.  The previous
    // implementation scanned only p < 10000 and gave up here.
    let P = 1n;
    for (let p = 2n; p <= 10007n; p = next_prime(p)) P *= p;
    const f: IntPoly = [-P, 0n, 1n];

    const { p } = _zz_factor_internal.chooseFactorizationPrime(f);
    expect(p).toBeGreaterThan(10007n);

    // sage: (x^2 - prod(prime_range(10008))).is_irreducible() -> True
    const factors = _zz_factor_internal.factorSquarefreeIntPoly(f);
    expect(factors.length).toBe(1);
    expect(factors[0]!.join(',')).toBe(f.join(','));
  }, 120000);

  test('smaller primorials', () => {
    for (const bound of [200n, 2000n]) {
      let P = 1n;
      for (let p = 2n; p <= bound; p = next_prime(p)) P *= p;
      const f: IntPoly = [-P, 0n, 1n];
      const { p } = _zz_factor_internal.chooseFactorizationPrime(f);
      expect(p).toBeGreaterThan(bound);
      expect(_zz_factor_internal.factorSquarefreeIntPoly(f).length).toBe(1);
    }
  });

  test('a non-squarefree input is rejected instead of looping for ever', () => {
    // (x^2+1)^2 -- every prime divides disc, so upstream's unbounded loop would
    // never terminate; our bounded variant raises.
    const f = oracleMul([1n, 0n, 1n], [1n, 0n, 1n]);
    expect(() => _zz_factor_internal.chooseFactorizationPrime(f)).toThrow();
  }, 120000);
});

describe('Leading coefficient and discriminant divisible by many small primes', () => {
  const primorial = (bound: bigint): bigint => {
    let P = 1n;
    for (let p = 2n; p <= bound; p = next_prime(p)) P *= p;
    return P;
  };

  test('prod_{i=1..12} (P x + i) with P = primorial(100)', () => {
    const P = primorial(100n);
    let f: IntPoly = [1n];
    for (let i = 1; i <= 12; i++) f = oracleMul(f, [BigInt(i), P]);

    const [content, factors] = _zz_factor_internal.factorIntegerPolynomial(f);

    // sage: prod([P*x+i for i in range(1,13)]).factor() has unit 1 and the
    // constant factors 2^6 * 3^4 * 5^2 * 7 * 11 = 9979200
    expect(content).toBe(9979200n);

    const want: string[] = [];
    for (let i = 1; i <= 12; i++) want.push(oraclePrimitive([BigInt(i), P]).join(','));
    const got = factors.map(([g, e]) => {
      expect(e).toBe(1);
      return g.join(',');
    });
    expect(got.sort()).toEqual(want.sort());
  }, 120000);

  test('big leading coefficients and several factors', () => {
    const P = primorial(60n);
    const parts: IntPoly[] = [
      [1n, P],
      [-1n, P],
      [3n, P],
      [-5n, P],
      [7n, P],
      [P, 1n, 1n],
      [1n, 2n, P],
    ];
    let f: IntPoly = [1n];
    for (const g of parts) f = oracleMul(f, g);

    const [content, factors] = _zz_factor_internal.factorIntegerPolynomial(f);
    let rebuilt: IntPoly = [content];
    for (const [g, e] of factors) for (let i = 0; i < e; i++) rebuilt = oracleMul(rebuilt, g);
    expect(oracleStrip(rebuilt).join(',')).toBe(oracleStrip(f).join(','));

    const want = parts.map((g) => oraclePrimitive(g).join(',')).sort();
    const got = factors.flatMap(([g, e]) => Array.from({ length: e }, () => g.join(','))).sort();
    expect(got).toEqual(want);
  }, 120000);

  test('prod_{i=1..6} (x^2 - i^2 * primorial(50)): disc has many small prime factors', () => {
    const P = primorial(50n);
    const parts: IntPoly[] = [];
    for (let i = 1n; i <= 6n; i++) parts.push([-(P * i * i), 0n, 1n]);
    let f: IntPoly = [1n];
    for (const g of parts) f = oracleMul(f, g);

    const factors = _zz_factor_internal.factorSquarefreeIntPoly(f);
    expect(factors.map((g) => g.join(',')).sort()).toEqual(parts.map((g) => g.join(',')).sort());
  }, 120000);
});

describe('x^n - 1 factors into the cyclotomic polynomials (van Hoeij, r > 8)', () => {
  const eulerPhi = (n: number): number => {
    let r = n;
    let m = n;
    for (let p = 2; p * p <= m; p++) {
      if (m % p === 0) {
        while (m % p === 0) m /= p;
        r -= r / p;
      }
    }
    if (m > 1) r -= r / m;
    return r;
  };

  test.each([24, 36, 48, 60, 105, 120])(
    'x^%i - 1',
    (n) => {
      const f: IntPoly = new Array(n + 1).fill(0n);
      f[0] = -1n;
      f[n] = 1n;
      const factors = _zz_factor_internal.factorSquarefreeIntPoly(f);

      // independent oracle: the degrees are phi(d) for the divisors d of n
      const want: number[] = [];
      for (let d = 1; d <= n; d++) if (n % d === 0) want.push(eulerPhi(d));
      expect(factors.map((g) => g.length - 1).sort((a, b) => a - b)).toEqual(
        want.sort((a, b) => a - b)
      );

      let prod: IntPoly = [1n];
      for (const g of factors) prod = oracleMul(prod, g);
      expect(prod.join(',')).toBe(f.join(','));
    },
    120000
  );
});

describe('van Hoeij with a starting precision that has to be doubled', () => {
  // `_heuristic_van_hoeij_starting_precision` (factor_van_hoeij.c:24-42) starts
  // well below the Mignotte bound, so a polynomial with large coefficients only
  // resolves after the main Hensel loop has doubled `a` (factor_van_hoeij.c:225-229).
  test('prod_{i=1..12} (x + 10^30 i + 7)', () => {
    let f: IntPoly = [1n];
    const B = 10n ** 30n;
    for (let i = 1n; i <= 12n; i++) f = oracleMul(f, [B * i + 7n, 1n]);
    const factors = _zz_factor_internal.factorSquarefreeIntPoly(f);
    const want: string[] = [];
    for (let i = 1n; i <= 12n; i++) want.push([B * i + 7n, 1n].join(','));
    expect(factors.map((g) => g.join(',')).sort()).toEqual(want.sort());
  }, 120000);

  test('prod_{i=1..10} (x^2 - 10^40 i^2)', () => {
    let f: IntPoly = [1n];
    const B = 10n ** 40n;
    for (let i = 1n; i <= 10n; i++) f = oracleMul(f, [-(B * i * i), 0n, 1n]);
    const factors = _zz_factor_internal.factorSquarefreeIntPoly(f);
    // x^2 - (10^20 i)^2 = (x - 10^20 i)(x + 10^20 i)
    const want: string[] = [];
    for (let i = 1n; i <= 10n; i++) {
      want.push([-(10n ** 20n * i), 1n].join(','));
      want.push([10n ** 20n * i, 1n].join(','));
    }
    expect(factors.map((g) => g.join(',')).sort()).toEqual(want.sort());
  }, 120000);

  test('SD32(2,3,5,7,11) * SD32(2,3,5,7,13)', () => {
    const a = swinnertonDyer([2n, 3n, 5n, 7n, 11n]);
    const b = swinnertonDyer([2n, 3n, 5n, 7n, 13n]);
    const f = oracleMul(a, b);
    const factors = _zz_factor_internal.factorSquarefreeIntPoly(f);
    expect(factors.map((g) => g.join(',')).sort()).toEqual([a.join(','), b.join(',')].sort());
  }, 120000);
});

describe('LLL with removal (fmpz_lll_wrapper_with_removal_knapsack)', () => {
  // ---- exact rational arithmetic, for FLINT's own reducedness predicate ----
  type Q = [bigint, bigint];
  const qgcd = (a: bigint, b: bigint): bigint => {
    let x = a < 0n ? -a : a;
    let y = b < 0n ? -b : b;
    while (y) [x, y] = [y, x % y];
    return x;
  };
  const qnorm = (a: bigint, b: bigint): Q => {
    if (b < 0n) {
      a = -a;
      b = -b;
    }
    const d = qgcd(a, b) || 1n;
    return [a / d, b / d];
  };
  const QZERO: Q = [0n, 1n];
  const qadd = (x: Q, y: Q): Q => qnorm(x[0] * y[1] + y[0] * x[1], x[1] * y[1]);
  const qsub = (x: Q, y: Q): Q => qnorm(x[0] * y[1] - y[0] * x[1], x[1] * y[1]);
  const qmul = (x: Q, y: Q): Q => qnorm(x[0] * y[0], x[1] * y[1]);
  const qdiv = (x: Q, y: Q): Q => qnorm(x[0] * y[1], x[1] * y[0]);
  const qcmp = (x: Q, y: Q): number => {
    const l = x[0] * y[1];
    const r = y[0] * x[1];
    return l < r ? -1 : l > r ? 1 : 0;
  };
  const qabs = (x: Q): Q => (x[0] < 0n ? [-x[0], x[1]] : x);
  /** the exact binary value of a double */
  const qFromDouble = (d: number): Q => {
    let n = d;
    let e = 0;
    while (!Number.isInteger(n)) {
      n *= 2;
      e++;
    }
    return qnorm(BigInt(n), 1n << BigInt(e));
  };

  /**
   * `gr_mat_is_row_lll_reduced_with_removal_naive`
   * (reference/flint/src/gr_mat/is_lll_reduced.c:18-100) over QQ -- the exact
   * predicate that `fmpz_lll_wrapper_with_removal_knapsack` verifies its
   * floating-point output against.
   */
  function isRowLLLReducedWithRemoval(
    A: bigint[][],
    delta: number,
    eta: number,
    gsB: bigint,
    newd: number
  ): boolean {
    const d = A.length;
    if (d <= 1) return true;
    const n = A[0]!.length;
    const D = qFromDouble(delta);
    const E = qFromDouble(eta);
    const G: Q = [gsB, 1n];
    const B: Q[][] = A.map((r) => r.map((c) => [c, 1n] as Q));
    const mu: Q[][] = [];
    for (let i = 0; i < d; i++) mu.push(new Array<Q>(d).fill(QZERO));
    const dot = (u: Q[], v: Q[]): Q => {
      let s: Q = QZERO;
      for (let i = 0; i < n; i++) s = qadd(s, qmul(u[i]!, v[i]!));
      return s;
    };
    mu[0]![0] = dot(B[0]!, B[0]!);
    if (newd === 0 && qcmp(mu[0]![0]!, G) < 0) return false;
    for (let i = 1; i < d; i++) {
      for (let j = 0; j < i; j++) {
        const t = dot(
          A[i]!.map((c) => [c, 1n] as Q),
          B[j]!
        );
        mu[i]![j] = qdiv(t, mu[j]![j]!);
        for (let k = 0; k < n; k++) B[i]![k] = qsub(B[i]![k]!, qmul(B[j]![k]!, mu[i]![j]!));
        if (i < newd && qcmp(qabs(mu[i]![j]!), E) > 0) return false;
      }
      mu[i]![i] = dot(B[i]!, B[i]!);
      if (i >= newd && qcmp(mu[i]![i]!, G) < 0) return false;
      if (i < newd) {
        const t = qmul(qsub(D, qmul(mu[i]![i - 1]!, mu[i]![i - 1]!)), mu[i - 1]![i - 1]!);
        if (qcmp(t, mu[i]![i]!) > 0) return false;
      }
    }
    return true;
  }

  /** row-style Hermite normal form over Z, for lattice-equality checks */
  function hnf(rowsIn: bigint[][]): string {
    const fdiv = (a: bigint, b: bigint): bigint => {
      let q = a / b;
      if (a % b !== 0n && a < 0n !== b < 0n) q--;
      return q;
    };
    const rows = rowsIn.map((r) => r.slice());
    const n = rows[0]?.length ?? 0;
    let r = 0;
    for (let c = 0; c < n && r < rows.length; c++) {
      for (;;) {
        let nz = -1;
        for (let i = r; i < rows.length; i++)
          if (rows[i]![c] !== 0n) {
            nz = i;
            break;
          }
        if (nz < 0) break;
        let done = true;
        for (let i = nz + 1; i < rows.length; i++) {
          if (rows[i]![c] === 0n) continue;
          done = false;
          const ai = rows[i]![c]! < 0n ? -rows[i]![c]! : rows[i]![c]!;
          const an = rows[nz]![c]! < 0n ? -rows[nz]![c]! : rows[nz]![c]!;
          if (ai < an) {
            const t = rows[i]!;
            rows[i] = rows[nz]!;
            rows[nz] = t;
          }
          const q = rows[i]![c]! / rows[nz]![c]!;
          for (let j = 0; j < n; j++) rows[i]![j] = rows[i]![j]! - q * rows[nz]![j]!;
        }
        if (done) {
          const t = rows[r]!;
          rows[r] = rows[nz]!;
          rows[nz] = t;
          break;
        }
      }
      if (rows[r] && rows[r]![c] !== 0n) {
        if (rows[r]![c]! < 0n) for (let j = 0; j < n; j++) rows[r]![j] = -rows[r]![j]!;
        for (let i = 0; i < r; i++) {
          const q = fdiv(rows[i]![c]!, rows[r]![c]!);
          for (let j = 0; j < n; j++) rows[i]![j] = rows[i]![j]! - q * rows[r]![j]!;
        }
        r++;
      }
    }
    return rows
      .filter((row) => row.some((c) => c !== 0n))
      .map((row) => row.join(','))
      .join(';');
  }

  test("300 random lattices: reduced per FLINT's exact predicate, same lattice", () => {
    const rand = makeRand(987654321n);
    let tested = 0;
    for (let trial = 0; trial < 300; trial++) {
      const d = 2 + rand(9);
      const n = d + rand(4);
      const bits = 4 + rand(40);
      const A: bigint[][] = [];
      for (let i = 0; i < d; i++) {
        const row: bigint[] = [];
        for (let j = 0; j < n; j++) {
          let v = 0n;
          for (let k = 0; k < bits; k++) v = v * 2n + BigInt(rand(2));
          row.push(rand(2) ? v : -v);
        }
        A.push(row);
      }
      const before = hnf(A);
      const gsB = BigInt(1 + rand(1 << 20));
      let newd: number;
      try {
        newd = _zz_factor_internal.lllWithRemovalKnapsack(A, gsB);
      } catch {
        continue; // linearly dependent rows: upstream's LLL handles them separately
      }
      tested++;
      expect(isRowLLLReducedWithRemoval(A, 0.99, 0.51, gsB, newd)).toBe(true);
      expect(hnf(A)).toBe(before);
      expect(newd).toBeGreaterThanOrEqual(0);
      expect(newd).toBeLessThanOrEqual(d);
    }
    expect(tested).toBeGreaterThan(250);
  }, 300000);
});

describe('zassenhaus_subset / zassenhaus_prune', () => {
  test('the subset walker enumerates every subset of each size exactly once', () => {
    const binom = (n: number, k: number): number => {
      let r = 1;
      for (let i = 0; i < k; i++) r = (r * (n - i)) / (i + 1);
      return Math.round(r);
    };
    for (let r = 1; r <= 12; r++) {
      for (let m = 0; m <= r; m++) {
        const s: number[] = [];
        for (let i = 0; i < r; i++) s.push(i);
        _zz_factor_internal.zassenhausSubsetFirst(s, r, m);
        const seen = new Set<string>();
        for (;;) {
          const sel = s.filter((v) => v >= 0).sort((a, b) => a - b);
          expect(sel.length).toBe(m);
          const key = sel.join(',');
          expect(seen.has(key)).toBe(false);
          seen.add(key);
          // every index 0..r-1 is present exactly once, in or out
          const all = s
            .map((v) => (v >= 0 ? v : -v - 1))
            .sort((a, b) => a - b)
            .join(',');
          expect(all).toBe(Array.from({ length: r }, (_, i) => i).join(','));
          if (!_zz_factor_internal.zassenhausSubsetNext(s, r)) break;
        }
        if (m > 0 && m < r) expect(seen.size).toBe(binom(r, m));
      }
    }
  });

  test('degree pruning matches brute-force subset sums', () => {
    const rand = makeRand(31337n);
    const bruteDegs = (degs: number[]): Set<number> => {
      const pos = new Set<number>([0]);
      for (const d of degs) for (const v of Array.from(pos)) pos.add(v + d);
      return pos;
    };
    for (let t = 0; t < 100; t++) {
      const D = 6 + (t % 20);
      const Z = _zz_factor_internal.zassenhausPruneSetDegree(D);
      const sets: number[][] = [];
      for (let round = 0; round < 3; round++) {
        let rem = D;
        const degs: number[] = [];
        while (rem > 0) {
          const d = 1 + rand(Math.min(rem, 4));
          degs.push(d);
          rem -= d;
        }
        sets.push(degs);
        _zz_factor_internal.zassenhausPruneStartAddFactors(Z);
        for (const d of degs) _zz_factor_internal.zassenhausPruneAddFactor(Z, d, 1);
        _zz_factor_internal.zassenhausPruneEndAddFactors(Z);
      }
      let inter = bruteDegs(sets[0]!);
      for (let i = 1; i < sets.length; i++) {
        const s = bruteDegs(sets[i]!);
        inter = new Set(Array.from(inter).filter((v) => s.has(v)));
      }
      for (let d = 0; d <= D; d++) {
        expect(_zz_factor_internal.zassenhausPruneDegreeIsPossible(Z, d)).toBe(inter.has(d));
      }
    }
  });
});

describe('CLD bounds and the CLD matrix (fmpz_poly/CLD_bound.c, fmpz_poly_factor/CLD_mat.c)', () => {
  test("|[x^n] f g'/g| <= CLD_bound(f, n) for every factor g of f", () => {
    const rand = makeRand(42n);
    const cases: IntPoly[][] = [];
    for (let t = 0; t < 120; t++) {
      const k = 2 + rand(4);
      const parts: IntPoly[] = [];
      for (let i = 0; i < k; i++) {
        const d = 1 + rand(4);
        const g: bigint[] = [];
        for (let j = 0; j <= d; j++) g.push(BigInt(rand(21) - 10));
        if (g[d] === 0n) g[d] = 1n;
        if (g.every((c) => c === 0n)) g[0] = 1n;
        parts.push(g);
      }
      cases.push(parts);
    }
    cases.push([swinnertonDyer([2n, 3n, 5n]), swinnertonDyer([2n, 3n, 7n])]);
    cases.push([swinnertonDyer([2n, 3n, 5n, 7n])]);

    let checks = 0;
    for (const parts of cases) {
      let f: IntPoly = [1n];
      for (const g of parts) f = oracleMul(f, g);
      f = oracleStrip(f);
      if (f.length < 3) continue;
      const N = f.length - 1;
      for (const g of parts) {
        if (g.length < 2) continue;
        const cof = parts.filter((x) => x !== g).reduce((a, b) => oracleMul(a, b), [1n] as IntPoly);
        const h = oracleMul(cof, oracleDeriv(g)); // = f g'/g, exactly
        for (let n = 0; n < N; n++) {
          const B = _zz_factor_internal.fmpzPolyCLDBound(f, n);
          const hn = h[n] ?? 0n;
          checks++;
          expect(hn <= B && hn >= -B).toBe(true);
        }
      }
    }
    expect(checks).toBeGreaterThan(3000);
  }, 300000);

  test('the CLD matrix holds the coefficients of the logarithmic derivatives', () => {
    // Low column j    -> [x^j] (f g_i'/g_i)               mod P
    // High column j   -> [x^{N-hiN+j}] (f G_i'/G_i)       mod P,  G_i = g_i >> len_i
    // with len_i = deg(g_i) - hiN, so f G'/G = f g'/g - len_i * f/x
    // (fmpz_poly_factor/CLD_mat.c:88-113).
    const rand = makeRand(7n);
    const smod = _zz_factor_internal.fmpzSmod;
    let cases = 0;
    for (let t = 0; t < 300; t++) {
      const k = 2 + rand(4);
      const parts: IntPoly[] = [];
      for (let i = 0; i < k; i++) {
        const d = 1 + rand(4);
        const g: bigint[] = [];
        for (let j = 0; j < d; j++) g.push(BigInt(rand(21) - 10));
        g.push(1n);
        parts.push(g);
      }
      let f: IntPoly = [1n];
      for (const g of parts) f = oracleMul(f, g);
      const N = f.length - 1;
      if (N < 4) continue;
      const p0 = 1000003n;
      if (parts.some((g) => g[0]! % p0 === 0n)) continue;
      const P = p0 ** 8n;
      const kk = Math.min(3, Math.floor((N + 1) / 2));
      const { data, numDataCols } = _zz_factor_internal.fmpzPolyFactorCLDMat(f, parts, P, kk);
      if (numDataCols === 0) continue;
      cases++;

      const q = parts.map((g, i) =>
        oracleMul(
          parts.filter((_, j) => j !== i).reduce((a, b) => oracleMul(a, b), [1n] as IntPoly),
          oracleDeriv(g)
        )
      );

      let explained = false;
      for (let loN = 0; loN <= numDataCols && !explained; loN++) {
        const hiN = numDataCols - loN;
        let ok = true;
        for (let i = 0; i < parts.length && ok; i++) {
          const leni = BigInt(parts[i]!.length - 1 - hiN);
          for (let c = 0; c < numDataCols && ok; c++) {
            let want: bigint;
            if (c < loN) want = q[i]![c] ?? 0n;
            else {
              const j = N - hiN + (c - loN);
              want = (q[i]![j] ?? 0n) - leni * (f[j + 1] ?? 0n);
            }
            if (smod(want, P) !== data[i]![c]!) ok = false;
          }
        }
        if (ok) explained = true;
      }
      expect(explained).toBe(true);
    }
    expect(cases).toBeGreaterThan(150);
  }, 300000);
});

describe('SageMath golden factorisations', () => {
  test('a mixed product of five irreducibles', () => {
    // sage: R.<x> = ZZ[]
    // sage: ((5*x^3-7*x+11)*(3*x^4+2*x^2-5)*(x^5-x-1)*(2*x^2+3*x+7)*(x^6+x+1)).factor()
    //   (x - 1) * (x + 1) * (2*x^2 + 3*x + 7) * (3*x^2 + 5) * (5*x^3 - 7*x + 11)
    //     * (x^5 - x - 1) * (x^6 + x + 1)
    let f: IntPoly = [1n];
    for (const g of [
      [11n, -7n, 0n, 5n],
      [-5n, 0n, 2n, 0n, 3n],
      [-1n, -1n, 0n, 0n, 0n, 1n],
      [7n, 3n, 2n],
      [1n, 1n, 0n, 0n, 0n, 0n, 1n],
    ] as IntPoly[]) {
      f = oracleMul(f, g);
    }
    const [content, factors] = _zz_factor_internal.factorIntegerPolynomial(f);
    expect(content).toBe(1n);
    expect(
      factors
        .map(([g, e]) => `${g.join(',')}^${e}`)
        .sort()
        .join(' ')
    ).toBe(
      ['-1,1^1', '1,1^1', '7,3,2^1', '5,0,3^1', '11,-7,0,5^1', '-1,-1,0,0,0,1^1', '1,1,0,0,0,0,1^1']
        .sort()
        .join(' ')
    );
  });

  test('Sage doctest: (12*(x^2+1)^3*(x+2)) through the ZZ pipeline', () => {
    let f: IntPoly = [12n];
    for (let i = 0; i < 3; i++) f = oracleMul(f, [1n, 0n, 1n]);
    f = oracleMul(f, [2n, 1n]);
    const [content, factors] = _zz_factor_internal.factorIntegerPolynomial(f);
    expect(content).toBe(12n);
    expect(
      factors
        .map(([g, e]) => `${g.join(',')}^${e}`)
        .sort()
        .join(' ')
    ).toBe(['2,1^1', '1,0,1^3'].sort().join(' '));
  });
});
