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
import { CONWAY_POLYNOMIALS } from '../finite_rings/conway_polynomials.js';
import { FiniteFieldExtension, PrimeField } from '../finite_rings/finite_field_extension.js';
import { FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import { GF2 } from '../finite_rings/gf2.js';
import { Integer } from '../integer_ring.js';
import { Rational } from '../rational.js';
import { QQ } from '../rational_field.js';
import { type CoefficientRing, Polynomial, type RingElement } from './polynomial_element.js';
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
