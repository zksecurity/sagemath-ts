/**
 * Unit tests for multivariate polynomial module
 *
 * Tests cover:
 * - Basic arithmetic (add, sub, mul, neg, pow)
 * - Degree and coefficient computations
 * - Leading term, monomial, and coefficient (lt, lm, lc)
 * - Evaluation and partial evaluation
 * - Term orderings (lex, deglex, degrevlex)
 * - Monomials, coefficients, and exponents extraction
 * - Derivative computation
 */
import { describe, expect, test } from 'bun:test';
import { FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import {
  MPolynomial,
  addExponents,
  exponentToKey,
  keyToExponent,
  totalDegree,
} from './multi_polynomial_element.js';
import { MPolynomialRing, MPolynomialRingConstructor } from './multi_polynomial_ring.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';

/**
 * A simple wrapper around bigint that implements RingElement.
 * Used for testing polynomials over integers.
 */
class IntegerElement implements RingElement {
  readonly value: bigint;

  constructor(value: bigint | number) {
    this.value = typeof value === 'number' ? BigInt(value) : value;
  }

  add(other: IntegerElement): IntegerElement {
    return new IntegerElement(this.value + other.value);
  }

  sub(other: IntegerElement): IntegerElement {
    return new IntegerElement(this.value - other.value);
  }

  mul(other: IntegerElement): IntegerElement {
    return new IntegerElement(this.value * other.value);
  }

  div(other: IntegerElement): IntegerElement {
    return new IntegerElement(this.value / other.value);
  }

  neg(): IntegerElement {
    return new IntegerElement(-this.value);
  }

  eq(other: IntegerElement | number): boolean {
    if (typeof other === 'number') {
      return this.value === BigInt(other);
    }
    return this.value === other.value;
  }

  isZero(): boolean {
    return this.value === 0n;
  }

  toString(): string {
    return this.value.toString();
  }
}

/**
 * Integer ring that produces IntegerElement objects.
 */
const IntegerRing: CoefficientRing<IntegerElement> = {
  zero(): IntegerElement {
    return new IntegerElement(0n);
  },
  one(): IntegerElement {
    return new IntegerElement(1n);
  },
  __call__(x: unknown): IntegerElement {
    if (x instanceof IntegerElement) {
      return x;
    }
    if (typeof x === 'bigint') {
      return new IntegerElement(x);
    }
    if (typeof x === 'number') {
      return new IntegerElement(BigInt(x));
    }
    throw new Error(`Cannot convert ${typeof x} to IntegerElement`);
  },
  is_field(): boolean {
    return false;
  },
};

describe('MPolynomialRing construction', () => {
  test('create ring with variable names', () => {
    const R = new MPolynomialRing(IntegerRing, ['x', 'y', 'z']);
    expect(R.ngens()).toBe(3);
    expect(R.variable_names()).toEqual(['x', 'y', 'z']);
  });

  test('reject duplicate variable names', () => {
    expect(() => new MPolynomialRing(IntegerRing, ['x', 'x', 'y'])).toThrow();
  });

  test('reject empty variable list', () => {
    expect(() => new MPolynomialRing(IntegerRing, [])).toThrow();
  });

  test('constructor function returns ring and generators', () => {
    const [R, x, y] = MPolynomialRingConstructor(IntegerRing, ['x', 'y']);
    expect(R.ngens()).toBe(2);
    expect(x.isTerm()).toBe(true);
    expect(y.isTerm()).toBe(true);
  });
});

describe('MPolynomial basic operations', () => {
  const [R, x, y, z] = MPolynomialRingConstructor(IntegerRing, ['x', 'y', 'z']);

  test('zero and one', () => {
    expect(R.zero().isZero()).toBe(true);
    expect(R.one().isConstant()).toBe(true);
    expect(R.one().constantCoefficient().eq(1)).toBe(true);
  });

  test('generators', () => {
    const gens = R.gens();
    expect(gens.length).toBe(3);

    expect(x.isTerm()).toBe(true);
    expect(x.isMonomial()).toBe(true);
    expect(x.degree()).toBe(1);
    expect(x.degreeIn(0)).toBe(1);
    expect(x.degreeIn(1)).toBe(0);
    expect(x.degreeIn(2)).toBe(0);
  });

  test('gens_dict', () => {
    const { x: xGen, y: yGen, z: zGen } = R.gens_dict();
    expect(xGen.eq(x)).toBe(true);
    expect(yGen.eq(y)).toBe(true);
    expect(zGen.eq(z)).toBe(true);
  });

  test('addition', () => {
    // x + y + z
    const sum = x.add(y).add(z);
    expect(sum.numberOfTerms()).toBe(3);
    expect(sum.degree()).toBe(1);
  });

  test('subtraction', () => {
    // x - x = 0
    const diff = x.sub(x);
    expect(diff.isZero()).toBe(true);
  });

  test('negation', () => {
    const negX = x.neg();
    expect(x.add(negX).isZero()).toBe(true);
  });

  test('multiplication of variables', () => {
    // x * y
    const prod = x.mul(y);
    expect(prod.isTerm()).toBe(true);
    expect(prod.degree()).toBe(2);
    expect(prod.coefficient([1, 1, 0]).eq(1)).toBe(true);
  });

  test('multiplication of polynomials', () => {
    // (x + 1) * (y + 1) = xy + x + y + 1
    const p1 = x.add(R.one());
    const p2 = y.add(R.one());
    const prod = p1.mul(p2);

    expect(prod.numberOfTerms()).toBe(4);
    expect(prod.coefficient([1, 1, 0]).eq(1)).toBe(true); // xy
    expect(prod.coefficient([1, 0, 0]).eq(1)).toBe(true); // x
    expect(prod.coefficient([0, 1, 0]).eq(1)).toBe(true); // y
    expect(prod.coefficient([0, 0, 0]).eq(1)).toBe(true); // 1
  });

  test('scalar multiplication', () => {
    // 3*x
    const scaled = x.scalarMul(new IntegerElement(3));
    expect(scaled.coefficient([1, 0, 0]).eq(3)).toBe(true);
  });

  test('power', () => {
    // (x + y)^2 = x^2 + 2xy + y^2
    const sum = x.add(y);
    const squared = sum.pow(2);

    expect(squared.coefficient([2, 0, 0]).eq(1)).toBe(true); // x^2
    expect(squared.coefficient([1, 1, 0]).eq(2)).toBe(true); // 2xy
    expect(squared.coefficient([0, 2, 0]).eq(1)).toBe(true); // y^2
  });
});

describe('MPolynomial degree operations', () => {
  const [R, x, y, z] = MPolynomialRingConstructor(IntegerRing, ['x', 'y', 'z']);

  test('total degree', () => {
    // x^2 + y*z + 1
    const f = x.pow(2).add(y.mul(z)).add(R.one());
    expect(f.degree()).toBe(2);
  });

  test('degree in specific variable', () => {
    // x^3*y^2 + z^5
    const f = x.pow(3).mul(y.pow(2)).add(z.pow(5));
    expect(f.degreeIn(0)).toBe(3); // x
    expect(f.degreeIn(1)).toBe(2); // y
    expect(f.degreeIn(2)).toBe(5); // z
  });

  test('zero polynomial degree', () => {
    expect(R.zero().degree()).toBe(-1);
    expect(R.zero().degreeIn(0)).toBe(-1);
  });

  test('constant polynomial degree', () => {
    const c = R.__call__(5);
    expect(c.degree()).toBe(0);
    expect(c.isConstant()).toBe(true);
  });
});

describe('MPolynomial coefficients and monomials', () => {
  const [R, x, y] = MPolynomialRingConstructor(IntegerRing, ['x', 'y']);

  test('coefficient extraction', () => {
    // 3x^2y + 2xy^2 + 5
    const f = x
      .pow(2)
      .mul(y)
      .scalarMul(new IntegerElement(3))
      .add(x.mul(y.pow(2)).scalarMul(new IntegerElement(2)))
      .add(R.__call__(5));

    expect(f.coefficient([2, 1]).eq(3)).toBe(true);
    expect(f.coefficient([1, 2]).eq(2)).toBe(true);
    expect(f.coefficient([0, 0]).eq(5)).toBe(true);
    expect(f.coefficient([1, 1]).eq(0)).toBe(true); // not present
  });

  test('monomials list', () => {
    // x^2 + xy + 1
    const f = x.pow(2).add(x.mul(y)).add(R.one());
    const monomials = f.monomials();

    expect(monomials.length).toBe(3);
    // All should be monomials (coefficient 1)
    for (const m of monomials) {
      expect(m.isMonomial()).toBe(true);
    }
  });

  test('coefficients list', () => {
    // 3x^2 + 2y + 5
    const f = x
      .pow(2)
      .scalarMul(new IntegerElement(3))
      .add(y.scalarMul(new IntegerElement(2)))
      .add(R.__call__(5));

    const coeffs = f.coefficients();
    expect(coeffs.length).toBe(3);
  });

  test('exponents list', () => {
    // x^2y + xy^2
    const f = x
      .pow(2)
      .mul(y)
      .add(x.mul(y.pow(2)));
    const exps = f.exponents();

    expect(exps.length).toBe(2);
    // Both have total degree 3
    for (const exp of exps) {
      expect(totalDegree(exp)).toBe(3);
    }
  });

  test('constant coefficient', () => {
    // x + y + 7
    const f = x.add(y).add(R.__call__(7));
    expect(f.constantCoefficient().eq(7)).toBe(true);

    // x + y (no constant term)
    const g = x.add(y);
    expect(g.constantCoefficient().eq(0)).toBe(true);
  });
});

describe('MPolynomial leading term operations', () => {
  const F = new FiniteFieldPrime(101n);

  test('degrevlex order (default)', () => {
    const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

    // x^2 + xy + y^2 + z^3
    const f = x.pow(2).add(x.mul(y)).add(y.pow(2)).add(z.pow(3));

    // In degrevlex: z^3 is leading (degree 3 > 2)
    const lt = f.lt();
    const lm = f.lm();
    const lc = f.lc();

    expect(f.degree()).toBe(3);
    expect(lm.coefficient([0, 0, 3]).eq(1)).toBe(true);
    expect(lc.eq(1)).toBe(true);
    expect(lt.eq(lm)).toBe(true);
  });

  test('lex order', () => {
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y'], 'lex');

    // y^3 + x
    const f = y.pow(3).add(x);

    // In lex: x is leading (x > y)
    const lm = f.lm();
    expect(lm.coefficient([1, 0]).eq(1)).toBe(true);
  });

  test('deglex order', () => {
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y'], 'deglex');

    // x + y^2
    const f = x.add(y.pow(2));

    // In deglex: y^2 is leading (degree 2 > 1)
    const lm = f.lm();
    expect(lm.coefficient([0, 2]).eq(1)).toBe(true);
  });
});

describe('MPolynomial evaluation', () => {
  const F = new FiniteFieldPrime(101n);
  const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

  test('evaluate constant polynomial', () => {
    const c = R.__call__(42);
    const result = c.evaluate({ x: F.__call__(1), y: F.__call__(2), z: F.__call__(3) });
    expect(result.eq(42)).toBe(true);
  });

  test('evaluate single variable', () => {
    // x^2 at x=5
    const f = x.pow(2);
    const result = f.evaluate({ x: F.__call__(5), y: F.__call__(0), z: F.__call__(0) });
    expect(result.eq(25)).toBe(true);
  });

  test('evaluate multivariate polynomial', () => {
    // x + 2y + 3z at (1, 2, 3) = 1 + 4 + 9 = 14
    const f = x.add(y.scalarMul(F.__call__(2))).add(z.scalarMul(F.__call__(3)));

    const result = f.evaluate({ x: F.__call__(1), y: F.__call__(2), z: F.__call__(3) });
    expect(result.eq(14)).toBe(true);
  });

  test('evaluate with array input', () => {
    // xy + z at (2, 3, 4) = 6 + 4 = 10
    const f = x.mul(y).add(z);
    const result = f.evaluate([F.__call__(2), F.__call__(3), F.__call__(4)]);
    expect(result.eq(10)).toBe(true);
  });

  test('partial evaluation', () => {
    // xy + z, substitute y=2
    const f = x.mul(y).add(z);
    const g = f.partialEvaluate({ y: F.__call__(2) });

    // g = 2x + z
    expect(g.coefficient([1, 0, 0]).eq(2)).toBe(true);
    expect(g.coefficient([0, 0, 1]).eq(1)).toBe(true);
  });

  test('partial evaluation with multiple variables', () => {
    // x + y + z, substitute x=1, z=3
    const f = x.add(y).add(z);
    const g = f.partialEvaluate({ x: F.__call__(1), z: F.__call__(3) });

    // g = y + 4
    expect(g.coefficient([0, 1, 0]).eq(1)).toBe(true);
    expect(g.constantCoefficient().eq(4)).toBe(true);
  });
});

describe('MPolynomial derivative', () => {
  const [R, x, y, z] = MPolynomialRingConstructor(IntegerRing, ['x', 'y', 'z']);

  test('derivative of constant', () => {
    const c = R.__call__(5);
    expect(c.derivative(0).isZero()).toBe(true);
    expect(c.derivative(1).isZero()).toBe(true);
  });

  test('derivative of variable', () => {
    // d/dx(x) = 1
    const dx = x.derivative(0);
    expect(dx.isConstant()).toBe(true);
    expect(dx.constantCoefficient().eq(1)).toBe(true);

    // d/dy(x) = 0
    expect(x.derivative(1).isZero()).toBe(true);
  });

  test('derivative of power', () => {
    // d/dx(x^3) = 3x^2
    const f = x.pow(3);
    const df = f.derivative(0);

    expect(df.degree()).toBe(2);
    expect(df.coefficient([2, 0, 0]).eq(3)).toBe(true);
  });

  test('derivative of multivariate polynomial', () => {
    // d/dx(x^2y + xy^2) = 2xy + y^2
    const f = x
      .pow(2)
      .mul(y)
      .add(x.mul(y.pow(2)));
    const df = f.derivative(0);

    expect(df.coefficient([1, 1, 0]).eq(2)).toBe(true); // 2xy
    expect(df.coefficient([0, 2, 0]).eq(1)).toBe(true); // y^2
  });
});

describe('MPolynomial homogeneity', () => {
  const [R, x, y, z] = MPolynomialRingConstructor(IntegerRing, ['x', 'y', 'z']);

  test('zero is homogeneous', () => {
    expect(R.zero().isHomogeneous()).toBe(true);
  });

  test('constant is homogeneous', () => {
    expect(R.one().isHomogeneous()).toBe(true);
  });

  test('single variable is homogeneous', () => {
    expect(x.isHomogeneous()).toBe(true);
  });

  test('sum of same degree is homogeneous', () => {
    // x^2 + xy + y^2 (all degree 2)
    const f = x.pow(2).add(x.mul(y)).add(y.pow(2));
    expect(f.isHomogeneous()).toBe(true);
  });

  test('mixed degrees is not homogeneous', () => {
    // x^2 + x + 1 (degrees 2, 1, 0)
    const f = x.pow(2).add(x).add(R.one());
    expect(f.isHomogeneous()).toBe(false);
  });
});

describe('MPolynomial equality', () => {
  const F = new FiniteFieldPrime(101n);
  const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y']);

  test('same polynomial is equal', () => {
    const f = x.pow(2).add(y);
    const g = x.pow(2).add(y);
    expect(f.eq(g)).toBe(true);
  });

  test('different polynomials are not equal', () => {
    const f = x.pow(2).add(y);
    const g = x.pow(2).add(y.pow(2));
    expect(f.eq(g)).toBe(false);
  });

  test('equality with constant', () => {
    const c = R.__call__(5);
    expect(c.eq(5)).toBe(true);
    expect(c.eq(6)).toBe(false);
  });

  test('zero equals zero', () => {
    expect(R.zero().eq(R.zero())).toBe(true);
    expect(R.zero().eq(0)).toBe(true);
  });
});

describe('MPolynomial string representation', () => {
  const [R, x, y] = MPolynomialRingConstructor(IntegerRing, ['x', 'y']);

  test('zero polynomial', () => {
    expect(R.zero().toString()).toBe('0');
  });

  test('constant polynomial', () => {
    expect(R.__call__(5).toString()).toBe('5');
  });

  test('single variable', () => {
    expect(x.toString()).toBe('x');
  });

  test('variable power', () => {
    expect(x.pow(3).toString()).toBe('x^3');
  });

  test('polynomial with multiple terms', () => {
    const f = x.pow(2).add(x.mul(y)).add(R.one());
    const str = f.toString();

    // Should contain all terms
    expect(str).toContain('x^2');
    expect(str).toContain('x*y');
    // Constant might be just "1"
  });
});

describe('MPolynomialRing element creation', () => {
  const F = new FiniteFieldPrime(101n);
  const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y']);

  test('create from number', () => {
    const c = R.__call__(42);
    expect(c.isConstant()).toBe(true);
    expect(c.constantCoefficient().eq(42)).toBe(true);
  });

  test('create from bigint', () => {
    const c = R.__call__(42n);
    expect(c.constantCoefficient().eq(42)).toBe(true);
  });

  test('create from base ring element', () => {
    const c = R.__call__(F.__call__(7));
    expect(c.constantCoefficient().eq(7)).toBe(true);
  });

  test('create from dictionary', () => {
    // { '2,0': 3, '0,1': 2 } = 3x^2 + 2y
    const f = R.__call__({ '2,0': 3, '0,1': 2 });
    expect(f.coefficient([2, 0]).eq(3)).toBe(true);
    expect(f.coefficient([0, 1]).eq(2)).toBe(true);
  });

  test('monomial creation', () => {
    const m = R.monomial([2, 3], F.__call__(5));
    expect(m.isTerm()).toBe(true);
    expect(m.coefficient([2, 3]).eq(5)).toBe(true);
  });

  test('fromTerms', () => {
    const f = R.fromTerms([
      [F.__call__(3), [2, 0]], // 3x^2
      [F.__call__(2), [0, 1]], // 2y
      [F.__call__(1), [0, 0]], // 1
    ]);

    expect(f.coefficient([2, 0]).eq(3)).toBe(true);
    expect(f.coefficient([0, 1]).eq(2)).toBe(true);
    expect(f.constantCoefficient().eq(1)).toBe(true);
  });
});

describe('Term order comparisons', () => {
  test('lex: x > y > z', () => {
    const F = new FiniteFieldPrime(101n);
    const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z'], 'lex');

    // x + y^2 + z^3: x should be leading
    const f = x.add(y.pow(2)).add(z.pow(3));
    const lm = f.lm();
    expect(lm.coefficient([1, 0, 0]).eq(1)).toBe(true);
  });

  test('deglex: higher degree first, then lex', () => {
    const F = new FiniteFieldPrime(101n);
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y'], 'deglex');

    // x^2 + y^2: same degree, x > y in lex
    const f = x.pow(2).add(y.pow(2));
    const lm = f.lm();
    expect(lm.coefficient([2, 0]).eq(1)).toBe(true);
  });

  test('degrevlex: higher degree first, then revlex', () => {
    const F = new FiniteFieldPrime(101n);
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y'], 'degrevlex');

    // x^2 + y^2: same degree, in degrevlex x^2 > y^2
    const f = x.pow(2).add(y.pow(2));
    const lm = f.lm();
    // In degrevlex, y^2 would be "smaller" so x^2 is leading
    expect(lm.coefficient([2, 0]).eq(1)).toBe(true);
  });
});

describe('Utility functions', () => {
  test('totalDegree', () => {
    expect(totalDegree([2, 3, 1])).toBe(6);
    expect(totalDegree([0, 0, 0])).toBe(0);
    expect(totalDegree([])).toBe(0);
  });

  test('addExponents', () => {
    expect(addExponents([1, 2], [3, 4])).toEqual([4, 6]);
    expect(addExponents([1, 2, 3], [0, 0])).toEqual([1, 2, 3]);
    expect(addExponents([], [1, 2])).toEqual([1, 2]);
  });

  test('exponentToKey and keyToExponent', () => {
    const exp = [2, 0, 3, 1];
    const key = exponentToKey(exp);
    const recovered = keyToExponent(key);
    expect(recovered).toEqual(exp);
  });
});

describe('MPolynomial sumcheck/GKR protocol methods', () => {
  const F = new FiniteFieldPrime(101n);

  describe('degrees()', () => {
    test('returns array of max degrees per variable', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      // f = x^2*y + y^3*z
      const f = x.pow(2).mul(y).add(y.pow(3).mul(z));

      expect(f.degrees()).toEqual([2, 3, 1]);
    });

    test('handles single variable polynomial', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      // f = x^5
      const f = x.pow(5);

      expect(f.degrees()).toEqual([5, 0, 0]);
    });

    test('handles constant polynomial', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      const f = R.__call__(42);

      expect(f.degrees()).toEqual([0, 0, 0]);
    });

    test('handles zero polynomial', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      expect(R.zero().degrees()).toEqual([-1, -1, -1]);
    });

    test('handles polynomial with all variables', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      // f = x^3 + y^2 + z^4
      const f = x.pow(3).add(y.pow(2)).add(z.pow(4));

      expect(f.degrees()).toEqual([3, 2, 4]);
    });

    test('matches degreeIn for each variable', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      const f = x.pow(2).mul(y.pow(3)).add(z.pow(4));
      const degs = f.degrees();

      expect(degs[0]).toBe(f.degreeIn(0));
      expect(degs[1]).toBe(f.degreeIn(1));
      expect(degs[2]).toBe(f.degreeIn(2));
    });
  });

  describe('variables()', () => {
    test('returns only appearing variables', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      // f = x^2 + z (y does not appear)
      const f = x.pow(2).add(z);

      const vars = f.variables();
      expect(vars.length).toBe(2);
      expect(vars[0]!.eq(x)).toBe(true);
      expect(vars[1]!.eq(z)).toBe(true);
    });

    test('returns all variables when all appear', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      const f = x.mul(y).mul(z);

      const vars = f.variables();
      expect(vars.length).toBe(3);
      expect(vars[0]!.eq(x)).toBe(true);
      expect(vars[1]!.eq(y)).toBe(true);
      expect(vars[2]!.eq(z)).toBe(true);
    });

    test('returns empty array for constant polynomial', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      const f = R.__call__(5);

      expect(f.variables()).toEqual([]);
    });

    test('returns empty array for zero polynomial', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      expect(R.zero().variables()).toEqual([]);
    });

    test('respects variable order', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      // f = z + x (y missing, order should still be x, z)
      const f = z.add(x);

      const vars = f.variables();
      expect(vars.length).toBe(2);
      expect(vars[0]!.eq(x)).toBe(true);
      expect(vars[1]!.eq(z)).toBe(true);
    });
  });

  describe('nvariables()', () => {
    test('counts appearing variables', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      // f = x^2 + z (y does not appear)
      const f = x.pow(2).add(z);

      expect(f.nvariables()).toBe(2);
    });

    test('equals length of variables()', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      const f = x.mul(y).add(z.pow(2));

      expect(f.nvariables()).toBe(f.variables().length);
    });

    test('returns 0 for constant polynomial', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      expect(R.__call__(42).nvariables()).toBe(0);
    });

    test('returns 0 for zero polynomial', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      expect(R.zero().nvariables()).toBe(0);
    });
  });

  describe('subs()', () => {
    test('works same as partialEvaluate', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      // f = xy + z
      const f = x.mul(y).add(z);

      const subsResult = f.subs({ y: F.__call__(2) });
      const partialResult = f.partialEvaluate({ y: F.__call__(2) });

      expect(subsResult.eq(partialResult)).toBe(true);
    });

    test('substitutes multiple variables', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      // f = x + y + z
      const f = x.add(y).add(z);

      // g = x + 3 + 5 = x + 8
      const g = f.subs({ y: F.__call__(3), z: F.__call__(5) });

      expect(g.coefficient([1, 0, 0]).eq(1)).toBe(true);
      expect(g.constantCoefficient().eq(8)).toBe(true);
    });

    test('substitutes single variable', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      // f = x^2 * y
      const f = x.pow(2).mul(y);

      // Substitute x = 3: result = 9 * y
      const g = f.subs({ x: F.__call__(3) });

      expect(g.coefficient([0, 1, 0]).eq(9)).toBe(true);
    });

    test('full substitution reduces to constant', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      // f = x + 2*y + z
      const f = x.add(y.scalarMul(F.__call__(2))).add(z);

      // Substitute all: 1 + 2*2 + 3 = 8
      const g = f.subs({ x: F.__call__(1), y: F.__call__(2), z: F.__call__(3) });

      expect(g.isConstant()).toBe(true);
      expect(g.constantCoefficient().eq(8)).toBe(true);
    });

    test('empty substitution returns same polynomial', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      const f = x.pow(2).add(y);
      const g = f.subs({});

      expect(g.eq(f)).toBe(true);
    });
  });

  describe('args()', () => {
    test('returns all ring generators', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      const f = x.pow(2).add(z);

      const args = f.args();
      expect(args.length).toBe(3);
      expect(args[0]!.eq(x)).toBe(true);
      expect(args[1]!.eq(y)).toBe(true);
      expect(args[2]!.eq(z)).toBe(true);
    });

    test('returns same as parent.gens()', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      const f = x.add(y);

      const args = f.args();
      const gens = R.gens();

      expect(args.length).toBe(gens.length);
      for (let i = 0; i < args.length; i++) {
        expect(args[i]!.eq(gens[i]!)).toBe(true);
      }
    });

    test('works for zero polynomial', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      const args = R.zero().args();
      expect(args.length).toBe(3);
    });

    test('works for constant polynomial', () => {
      const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

      const args = R.__call__(42).args();
      expect(args.length).toBe(3);
    });
  });

  describe('sumcheck protocol pattern', () => {
    test('SageMath sumcheck example pattern', () => {
      // Reference from sumcheck.sage:
      // (x0, x1, x2) = poly.args()
      // p1 = sum([poly.subs(x1=a, x2=b) for a in [0, 1] for b in [0, 1]])
      // assert p1.degree() == 1

      const [R, x0, x1, x2] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);

      // Create a simple multilinear polynomial
      // f = x0 + x1 + x2
      const poly = x0.add(x1).add(x2);

      // Destructure args
      const [v0, v1, v2] = poly.args();
      expect(v0!.eq(x0)).toBe(true);
      expect(v1!.eq(x1)).toBe(true);
      expect(v2!.eq(x2)).toBe(true);

      // Sum over boolean hypercube for x1 and x2
      // p1 = sum([poly.subs(x1=a, x2=b) for a in [0, 1] for b in [0, 1]])
      let p1 = R.zero();
      for (const a of [0, 1]) {
        for (const b of [0, 1]) {
          p1 = p1.add(poly.subs({ x1: F.__call__(a), x2: F.__call__(b) }));
        }
      }

      // p1 should be 4*x0 + (0+0+1+1+0+1+0+1) = 4*x0 + 4
      // (x0+0+0) + (x0+0+1) + (x0+1+0) + (x0+1+1) = 4*x0 + 4
      expect(p1.degree()).toBe(1);
      expect(p1.variables().length).toBe(1);
      expect(p1.variables()[0]!.eq(x0)).toBe(true);
    });

    test('degrees check for sumcheck', () => {
      const [R, x0, x1, x2] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);

      // Multilinear polynomial: x0*x1 + x1*x2 + x0
      const poly = x0.mul(x1).add(x1.mul(x2)).add(x0);

      // Each variable degree should be at most 1 for multilinear
      const degs = poly.degrees();
      expect(degs[0]).toBe(1); // x0
      expect(degs[1]).toBe(1); // x1
      expect(degs[2]).toBe(1); // x2

      // Sum of degrees array
      expect(degs.reduce((a, b) => a + b, 0)).toBeLessThanOrEqual(3);
    });
  });
});

describe('ZK constraint system example', () => {
  test('R1CS-like constraint encoding', () => {
    // Example: x * y = z constraint
    const F = new FiniteFieldPrime(101n);
    const [R, x, y, z] = MPolynomialRingConstructor(F, ['x', 'y', 'z']);

    // Constraint: x*y - z = 0
    const constraint = x.mul(y).sub(z);

    // Valid witness: x=3, y=7, z=21
    const witness1 = constraint.evaluate({
      x: F.__call__(3),
      y: F.__call__(7),
      z: F.__call__(21),
    });
    expect(witness1.isZero()).toBe(true);

    // Invalid witness: x=3, y=7, z=20
    const witness2 = constraint.evaluate({
      x: F.__call__(3),
      y: F.__call__(7),
      z: F.__call__(20),
    });
    expect(witness2.isZero()).toBe(false);
  });

  test('Boolean constraint encoding', () => {
    // Constraint: x * (1 - x) = 0 (x must be 0 or 1)
    const F = new FiniteFieldPrime(101n);
    const [R, x] = MPolynomialRingConstructor(F, ['x']);

    const constraint = x.mul(R.one().sub(x));

    // x = 0 satisfies
    expect(constraint.evaluate([F.__call__(0)]).isZero()).toBe(true);

    // x = 1 satisfies
    expect(constraint.evaluate([F.__call__(1)]).isZero()).toBe(true);

    // x = 2 does not satisfy
    expect(constraint.evaluate([F.__call__(2)]).isZero()).toBe(false);
  });

  test('Range check composition', () => {
    // Constraint: x^2 - y = 0 (prove x^2 = y)
    const F = new FiniteFieldPrime(101n);
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y']);

    const constraint = x.pow(2).sub(y);

    // Valid: x=5, y=25
    expect(constraint.evaluate({ x: F.__call__(5), y: F.__call__(25) }).isZero()).toBe(true);

    // Invalid: x=5, y=24
    expect(constraint.evaluate({ x: F.__call__(5), y: F.__call__(24) }).isZero()).toBe(false);
  });
});
