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
import { Polynomial } from './polynomial_element.js';
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
