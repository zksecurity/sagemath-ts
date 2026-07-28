/**
 * Unit tests for multilinear extension module
 *
 * Tests cover:
 * - Binary conversion utilities (closestPowerOfTwo, intToBinary, binaryToInt)
 * - Boolean hypercube generation
 * - Equality polynomial (selector polynomial)
 * - MLE interpolation (dense and sparse)
 * - Multilinearity verification
 */
import { describe, expect, test } from 'bun:test';
import { FiniteFieldPrime } from '../rings/finite_rings/finite_field_prime.js';
import { MPolynomialRingConstructor } from '../rings/polynomial/multi_polynomial_ring.js';
import {
  closestPowerOfTwo,
  intToBinary,
  binaryToInt,
  booleanHypercube,
  MAX_HYPERCUBE_DIM,
  eqPolynomial,
  multilinearExtension,
  sparseMultilinearExtension,
  isMultilinear,
  verifyMLEEvaluation,
} from './multilinear.js';

// Test field
const F = new FiniteFieldPrime(101n);

describe('closestPowerOfTwo', () => {
  test('returns 0 for n <= 1', () => {
    expect(closestPowerOfTwo(0)).toBe(0);
    expect(closestPowerOfTwo(1)).toBe(0);
  });

  test('returns correct values for powers of 2', () => {
    expect(closestPowerOfTwo(2)).toBe(1);
    expect(closestPowerOfTwo(4)).toBe(2);
    expect(closestPowerOfTwo(8)).toBe(3);
    expect(closestPowerOfTwo(16)).toBe(4);
    expect(closestPowerOfTwo(32)).toBe(5);
  });

  test('returns ceiling for non-powers of 2', () => {
    expect(closestPowerOfTwo(3)).toBe(2); // ceil(log2(3)) = 2
    expect(closestPowerOfTwo(5)).toBe(3); // ceil(log2(5)) = 3
    expect(closestPowerOfTwo(6)).toBe(3);
    expect(closestPowerOfTwo(7)).toBe(3);
    expect(closestPowerOfTwo(9)).toBe(4);
    expect(closestPowerOfTwo(100)).toBe(7);
  });
});

describe('intToBinary', () => {
  test('converts 0 correctly', () => {
    expect(intToBinary(0, 1)).toEqual([0]);
    expect(intToBinary(0, 4)).toEqual([0, 0, 0, 0]);
  });

  test('converts small integers correctly', () => {
    expect(intToBinary(1, 1)).toEqual([1]);
    expect(intToBinary(1, 4)).toEqual([0, 0, 0, 1]);
    expect(intToBinary(5, 3)).toEqual([1, 0, 1]);
    expect(intToBinary(5, 4)).toEqual([0, 1, 0, 1]);
    expect(intToBinary(5, 8)).toEqual([0, 0, 0, 0, 0, 1, 0, 1]);
  });

  test('converts powers of 2 correctly', () => {
    expect(intToBinary(2, 3)).toEqual([0, 1, 0]);
    expect(intToBinary(4, 3)).toEqual([1, 0, 0]);
    expect(intToBinary(8, 4)).toEqual([1, 0, 0, 0]);
  });

  test('converts all-ones correctly', () => {
    expect(intToBinary(7, 3)).toEqual([1, 1, 1]);
    expect(intToBinary(15, 4)).toEqual([1, 1, 1, 1]);
  });

  test('throws for negative values', () => {
    expect(() => intToBinary(-1, 4)).toThrow();
  });

  test('throws for negative numBits', () => {
    expect(() => intToBinary(5, -1)).toThrow();
  });
});

describe('binaryToInt', () => {
  test('converts empty array to 0', () => {
    expect(binaryToInt([])).toBe(0n);
  });

  test('converts single bits correctly', () => {
    expect(binaryToInt([0])).toBe(0n);
    expect(binaryToInt([1])).toBe(1n);
  });

  test('converts multi-bit arrays correctly', () => {
    expect(binaryToInt([1, 0, 1])).toBe(5n);
    expect(binaryToInt([0, 1, 0, 1])).toBe(5n);
    expect(binaryToInt([1, 1, 1])).toBe(7n);
    expect(binaryToInt([1, 0, 0, 0])).toBe(8n);
  });

  test('roundtrips with intToBinary', () => {
    for (let i = 0; i < 20; i++) {
      const bits = intToBinary(i, 5);
      expect(binaryToInt(bits)).toBe(BigInt(i));
    }
  });

  test('is exact above 2^53', () => {
    // 64 bits: a `number` result would round here (parseInt on the joined
    // string loses the low bits).
    const bits = Array.from({ length: 64 }, (_, i) => (i % 3 === 0 ? 1 : 0));
    const value = binaryToInt(bits);
    expect(value.toString(2).padStart(64, '0')).toBe(bits.join(''));
    expect(value).toBe(0b1001001001001001001001001001001001001001001001001001001001001001n);
  });

  test('rejects malformed bit arrays', () => {
    // parseInt('121', 2) silently returns 1; we must not accept that.
    expect(() => binaryToInt([1, 2, 1])).toThrow('bits must contain only 0 and 1 (got 2)');
    expect(() => binaryToInt([0, -1])).toThrow();
    expect(() => binaryToInt([0.5])).toThrow();
  });
});

describe('booleanHypercube', () => {
  test('n=0 returns single empty point', () => {
    expect(booleanHypercube(0)).toEqual([[]]);
  });

  test('n=1 returns two points', () => {
    expect(booleanHypercube(1)).toEqual([[0], [1]]);
  });

  test('n=2 returns four points in order', () => {
    expect(booleanHypercube(2)).toEqual([
      [0, 0],
      [0, 1],
      [1, 0],
      [1, 1],
    ]);
  });

  test('n=3 returns eight points in order', () => {
    expect(booleanHypercube(3)).toEqual([
      [0, 0, 0],
      [0, 0, 1],
      [0, 1, 0],
      [0, 1, 1],
      [1, 0, 0],
      [1, 0, 1],
      [1, 1, 0],
      [1, 1, 1],
    ]);
  });

  test('returns 2^n points', () => {
    for (let n = 0; n <= 5; n++) {
      expect(booleanHypercube(n).length).toBe(2 ** n);
    }
  });

  test('throws for negative n', () => {
    expect(() => booleanHypercube(-1)).toThrow();
  });

  test('rejects dimensions it cannot materialize instead of wrapping around', () => {
    // With `1 << n` the 32-bit shift made n = 31 produce 0 points and n = 32
    // produce 1 point, silently.
    expect(() => booleanHypercube(31)).toThrow(`n must be at most ${MAX_HYPERCUBE_DIM}`);
    expect(() => booleanHypercube(32)).toThrow(`n must be at most ${MAX_HYPERCUBE_DIM}`);
    expect(() => booleanHypercube(MAX_HYPERCUBE_DIM + 1)).toThrow();
  });

  test('every point is distinct and has n coordinates', () => {
    const n = 4;
    const points = booleanHypercube(n);
    expect(points.length).toBe(16);
    const seen = new Set(points.map((p) => p.join('')));
    expect(seen.size).toBe(16);
    for (const p of points) {
      expect(p.length).toBe(n);
    }
  });
});

describe('eqPolynomial', () => {
  test('evaluates to 1 at matching point', () => {
    const [R, x0, x1] = MPolynomialRingConstructor(F, ['x0', 'x1']);
    const vars = [x0, x1];

    // eq([0, 0], x) at (0, 0) = 1
    const eq00 = eqPolynomial([0, 0], vars, R);
    expect(eq00.evaluate({ x0: F.__call__(0), x1: F.__call__(0) }).eq(1)).toBe(true);

    // eq([1, 0], x) at (1, 0) = 1
    const eq10 = eqPolynomial([1, 0], vars, R);
    expect(eq10.evaluate({ x0: F.__call__(1), x1: F.__call__(0) }).eq(1)).toBe(true);

    // eq([1, 1], x) at (1, 1) = 1
    const eq11 = eqPolynomial([1, 1], vars, R);
    expect(eq11.evaluate({ x0: F.__call__(1), x1: F.__call__(1) }).eq(1)).toBe(true);
  });

  test('evaluates to 0 at non-matching points', () => {
    const [R, x0, x1] = MPolynomialRingConstructor(F, ['x0', 'x1']);
    const vars = [x0, x1];

    const eq00 = eqPolynomial([0, 0], vars, R);

    // Should be 0 at all other points
    expect(eq00.evaluate({ x0: F.__call__(0), x1: F.__call__(1) }).eq(0)).toBe(true);
    expect(eq00.evaluate({ x0: F.__call__(1), x1: F.__call__(0) }).eq(0)).toBe(true);
    expect(eq00.evaluate({ x0: F.__call__(1), x1: F.__call__(1) }).eq(0)).toBe(true);
  });

  test('accepts integer base for bit decomposition', () => {
    const [R, x0, x1, x2] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);
    const vars = [x0, x1, x2];

    // 5 in binary is [1, 0, 1]
    const eq5 = eqPolynomial(5, vars, R);
    expect(eq5.evaluate({ x0: F.__call__(1), x1: F.__call__(0), x2: F.__call__(1) }).eq(1)).toBe(
      true
    );
    expect(eq5.evaluate({ x0: F.__call__(0), x1: F.__call__(0), x2: F.__call__(0) }).eq(0)).toBe(
      true
    );
  });

  test('is multilinear', () => {
    const [R, x0, x1, x2] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);
    const vars = [x0, x1, x2];

    for (let i = 0; i < 8; i++) {
      const eq = eqPolynomial(i, vars, R);
      expect(isMultilinear(eq)).toBe(true);
    }
  });

  test('throws for mismatched lengths', () => {
    const [R, x0, x1] = MPolynomialRingConstructor(F, ['x0', 'x1']);
    const vars = [x0, x1];

    expect(() => eqPolynomial([0, 0, 0], vars, R)).toThrow();
  });
});

describe('multilinearExtension', () => {
  test('interpolates single value as constant', () => {
    const values = [F.__call__(42n)];
    const poly = multilinearExtension(values, F);

    expect(poly.isConstant()).toBe(true);
    expect(poly.constantCoefficient().eq(42)).toBe(true);
  });

  test('interpolates two values correctly', () => {
    const values = [F.__call__(3n), F.__call__(7n)];
    const poly = multilinearExtension(values, F);

    // At x=0, should return 3
    expect(poly.evaluate({ x0: F.__call__(0) }).eq(3)).toBe(true);
    // At x=1, should return 7
    expect(poly.evaluate({ x0: F.__call__(1) }).eq(7)).toBe(true);
  });

  test('interpolates four values correctly', () => {
    const values = [F.__call__(1n), F.__call__(2n), F.__call__(3n), F.__call__(4n)];
    const poly = multilinearExtension(values, F);

    // Verify at all boolean hypercube points
    expect(poly.evaluate({ x0: F.__call__(0), x1: F.__call__(0) }).eq(1)).toBe(true);
    expect(poly.evaluate({ x0: F.__call__(0), x1: F.__call__(1) }).eq(2)).toBe(true);
    expect(poly.evaluate({ x0: F.__call__(1), x1: F.__call__(0) }).eq(3)).toBe(true);
    expect(poly.evaluate({ x0: F.__call__(1), x1: F.__call__(1) }).eq(4)).toBe(true);
  });

  test('handles non-power-of-2 length (missing values treated as 0)', () => {
    // 6 values requires 3 variables (8 points), last 2 treated as 0
    const values = [
      F.__call__(9n),
      F.__call__(2n),
      F.__call__(5n),
      F.__call__(4n),
      F.__call__(5n),
      F.__call__(7n),
    ];
    const poly = multilinearExtension(values, F);

    // Verify at provided values
    for (let i = 0; i < values.length; i++) {
      const bits = intToBinary(i, 3);
      const point = {
        x0: F.__call__(bits[0]!),
        x1: F.__call__(bits[1]!),
        x2: F.__call__(bits[2]!),
      };
      expect(poly.evaluate(point).eq(values[i]!)).toBe(true);
    }

    // Missing values (indices 6 and 7) evaluate to 0
    expect(
      poly.evaluate({ x0: F.__call__(1), x1: F.__call__(1), x2: F.__call__(0) }).eq(0)
    ).toBe(true);
    expect(
      poly.evaluate({ x0: F.__call__(1), x1: F.__call__(1), x2: F.__call__(1) }).eq(0)
    ).toBe(true);
  });

  test('result is multilinear', () => {
    const values = [
      F.__call__(1n),
      F.__call__(2n),
      F.__call__(3n),
      F.__call__(4n),
      F.__call__(5n),
      F.__call__(6n),
      F.__call__(7n),
      F.__call__(8n),
    ];
    const poly = multilinearExtension(values, F);

    expect(isMultilinear(poly)).toBe(true);
  });

  test('throws for empty values', () => {
    expect(() => multilinearExtension([], F)).toThrow();
  });

  test('accepts custom ring', () => {
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y']);
    const values = [F.__call__(1n), F.__call__(2n), F.__call__(3n), F.__call__(4n)];
    const poly = multilinearExtension(values, F, R);

    expect(poly.parent).toBe(R);
    expect(poly.evaluate({ x: F.__call__(0), y: F.__call__(0) }).eq(1)).toBe(true);
    expect(poly.evaluate({ x: F.__call__(1), y: F.__call__(1) }).eq(4)).toBe(true);
  });

  test('rejects a degenerate variable list (blueprint sanity check)', () => {
    // interpolate(values, [x0, x0]) in mle.sage:59-62. Verified in SageMath:
    //   interpolate([1,0,0,0], [x0, x0]) -> ValueError("Polynomial is not linear")
    //   interpolate([1,2,3,4], [x0, x0]) -> 3*x0 + 1 (accidentally multilinear)
    const [R, x0] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);
    const one = [F.__call__(1n), F.__call__(0n), F.__call__(0n), F.__call__(0n)];
    expect(() => multilinearExtension(one, F, R, [x0, x0])).toThrow('Polynomial is not linear');

    const ramp = [F.__call__(1n), F.__call__(2n), F.__call__(3n), F.__call__(4n)];
    const collapsed = multilinearExtension(ramp, F, R, [x0, x0]);
    expect(collapsed.toString()).toBe('3*x0 + 1');
  });
});

describe('sparseMultilinearExtension', () => {
  describe('selector polynomial (array of indices)', () => {
    test('empty array returns zero polynomial', () => {
      const [R, x0, x1] = MPolynomialRingConstructor(F, ['x0', 'x1']);
      const poly = sparseMultilinearExtension([], [x0, x1], R);

      expect(poly.isZero()).toBe(true);
    });

    test('single index creates selector for that point', () => {
      const [R, x0, x1] = MPolynomialRingConstructor(F, ['x0', 'x1']);
      const poly = sparseMultilinearExtension([2], [x0, x1], R);

      // Index 2 = (1, 0) in binary
      expect(poly.evaluate({ x0: F.__call__(1), x1: F.__call__(0) }).eq(1)).toBe(true);
      expect(poly.evaluate({ x0: F.__call__(0), x1: F.__call__(0) }).eq(0)).toBe(true);
      expect(poly.evaluate({ x0: F.__call__(0), x1: F.__call__(1) }).eq(0)).toBe(true);
      expect(poly.evaluate({ x0: F.__call__(1), x1: F.__call__(1) }).eq(0)).toBe(true);
    });

    test('multiple indices create sum of selectors', () => {
      const [R, x0, x1, x2] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);
      const vars = [x0, x1, x2];
      const poly = sparseMultilinearExtension([0, 2, 5], vars, R);

      // Should be 1 at indices 0, 2, 5 and 0 elsewhere
      for (let i = 0; i < 8; i++) {
        const bits = intToBinary(i, 3);
        const point = {
          x0: F.__call__(bits[0]!),
          x1: F.__call__(bits[1]!),
          x2: F.__call__(bits[2]!),
        };
        const expected = [0, 2, 5].includes(i) ? 1 : 0;
        expect(poly.evaluate(point).eq(expected)).toBe(true);
      }
    });
  });

  describe('sparse values (Map)', () => {
    test('empty map returns zero polynomial', () => {
      const [R, x0, x1] = MPolynomialRingConstructor(F, ['x0', 'x1']);
      const poly = sparseMultilinearExtension(new Map(), [x0, x1], R);

      expect(poly.isZero()).toBe(true);
    });

    test('single entry creates scaled selector', () => {
      const [R, x0, x1] = MPolynomialRingConstructor(F, ['x0', 'x1']);
      const entries = new Map([[1, F.__call__(5n)]]);
      const poly = sparseMultilinearExtension(entries, [x0, x1], R);

      // Index 1 = (0, 1) in binary, value 5
      expect(poly.evaluate({ x0: F.__call__(0), x1: F.__call__(1) }).eq(5)).toBe(true);
      expect(poly.evaluate({ x0: F.__call__(0), x1: F.__call__(0) }).eq(0)).toBe(true);
      expect(poly.evaluate({ x0: F.__call__(1), x1: F.__call__(0) }).eq(0)).toBe(true);
      expect(poly.evaluate({ x0: F.__call__(1), x1: F.__call__(1) }).eq(0)).toBe(true);
    });

    test('multiple entries create weighted sum', () => {
      const [R, x0, x1, x2] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);
      const vars = [x0, x1, x2];
      const entries = new Map([
        [0, F.__call__(9n)],
        [2, F.__call__(5n)],
        [5, F.__call__(2n)],
      ]);
      const poly = sparseMultilinearExtension(entries, vars, R);

      // Verify at specified entries
      expect(
        poly.evaluate({ x0: F.__call__(0), x1: F.__call__(0), x2: F.__call__(0) }).eq(9)
      ).toBe(true);
      expect(
        poly.evaluate({ x0: F.__call__(0), x1: F.__call__(1), x2: F.__call__(0) }).eq(5)
      ).toBe(true);
      expect(
        poly.evaluate({ x0: F.__call__(1), x1: F.__call__(0), x2: F.__call__(1) }).eq(2)
      ).toBe(true);

      // Verify zeros at unspecified entries
      expect(
        poly.evaluate({ x0: F.__call__(0), x1: F.__call__(0), x2: F.__call__(1) }).eq(0)
      ).toBe(true);
    });
  });

  test('result is multilinear', () => {
    const [R, x0, x1, x2] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);
    const vars = [x0, x1, x2];

    const selectorPoly = sparseMultilinearExtension([0, 3, 5, 7], vars, R);
    expect(isMultilinear(selectorPoly)).toBe(true);

    const valuePoly = sparseMultilinearExtension(
      new Map([
        [1, F.__call__(3n)],
        [4, F.__call__(7n)],
      ]),
      vars,
      R
    );
    expect(isMultilinear(valuePoly)).toBe(true);
  });

  test('throws for empty variables array', () => {
    const [R] = MPolynomialRingConstructor(F, ['x']);
    expect(() => sparseMultilinearExtension([0], [], R)).toThrow();
  });

  test('rejects a degenerate variable list (blueprint sanity check)', () => {
    // interpolate_sparse in mle.sage:127-129. Verified in SageMath:
    //   interpolate_sparse([0, 3], [x0, x0]) -> ValueError("Polynomial is not linear")
    const [R, x0] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);
    expect(() => sparseMultilinearExtension([0, 3], [x0, x0], R)).toThrow(
      'Polynomial is not linear'
    );
    expect(() =>
      sparseMultilinearExtension(new Map([[0, F.__call__(1n)]]), [x0, x0], R)
    ).toThrow('Polynomial is not linear');
  });

  test('single index is the selector, not the constant (documented deviation)', () => {
    // The blueprint short-circuits `interpolate_sparse([5], vars)` to the
    // constant 5 (mle.sage:108-112, a leftover debug branch that still prints).
    // We return eq(5, x) instead; verified against SageMath that the blueprint
    // really does return the constant 5 there.
    const [R, x0, x1, x2] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);
    const vars = [x0, x1, x2];
    const poly = sparseMultilinearExtension([5], vars, R);

    for (let i = 0; i < 8; i++) {
      const bits = intToBinary(i, 3);
      const point = {
        x0: F.__call__(bits[0]!),
        x1: F.__call__(bits[1]!),
        x2: F.__call__(bits[2]!),
      };
      expect(poly.evaluate(point).eq(i === 5 ? 1 : 0)).toBe(true);
    }
  });
});

describe('isMultilinear', () => {
  test('zero polynomial is multilinear', () => {
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y']);
    expect(isMultilinear(R.zero())).toBe(true);
  });

  test('constant polynomial is multilinear', () => {
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y']);
    expect(isMultilinear(R.__call__(42))).toBe(true);
  });

  test('single variable is multilinear', () => {
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y']);
    expect(isMultilinear(x)).toBe(true);
    expect(isMultilinear(y)).toBe(true);
  });

  test('x + y + xy is multilinear', () => {
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y']);
    const poly = x.add(y).add(x.mul(y));
    expect(isMultilinear(poly)).toBe(true);
  });

  test('x^2 is not multilinear', () => {
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y']);
    expect(isMultilinear(x.pow(2))).toBe(false);
  });

  test('xy^2 is not multilinear', () => {
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y']);
    expect(isMultilinear(x.mul(y.pow(2)))).toBe(false);
  });

  test('x^2 + y is not multilinear', () => {
    const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y']);
    expect(isMultilinear(x.pow(2).add(y))).toBe(false);
  });
});

describe('verifyMLEEvaluation', () => {
  test('returns true for correct evaluation', () => {
    const values = [F.__call__(1n), F.__call__(2n), F.__call__(3n), F.__call__(4n)];
    const poly = multilinearExtension(values, F);

    expect(verifyMLEEvaluation(poly, { x0: F.__call__(0), x1: F.__call__(0) }, F.__call__(1n))).toBe(
      true
    );
    expect(verifyMLEEvaluation(poly, { x0: F.__call__(1), x1: F.__call__(1) }, F.__call__(4n))).toBe(
      true
    );
  });

  test('returns false for incorrect evaluation', () => {
    const values = [F.__call__(1n), F.__call__(2n), F.__call__(3n), F.__call__(4n)];
    const poly = multilinearExtension(values, F);

    expect(verifyMLEEvaluation(poly, { x0: F.__call__(0), x1: F.__call__(0) }, F.__call__(2n))).toBe(
      false
    );
    expect(verifyMLEEvaluation(poly, { x0: F.__call__(1), x1: F.__call__(1) }, F.__call__(3n))).toBe(
      false
    );
  });
});

describe('MLE integration tests', () => {
  test('MLE matches original values at all hypercube points', () => {
    // Test with various value arrays
    const testCases = [
      [1n, 2n], // 2 values
      [1n, 2n, 3n, 4n], // 4 values
      [9n, 2n, 5n, 4n, 5n, 7n, 3n, 1n], // 8 values
    ];

    for (const vals of testCases) {
      const values = vals.map((v) => F.__call__(v));
      const numVars = closestPowerOfTwo(values.length);
      const poly = multilinearExtension(values, F);

      for (let i = 0; i < values.length; i++) {
        const bits = intToBinary(i, numVars);
        const point: Record<string, typeof values[0]> = {};
        for (let j = 0; j < numVars; j++) {
          point[`x${j}`] = F.__call__(bits[j]!);
        }
        expect(poly.evaluate(point).eq(values[i]!)).toBe(true);
      }
    }
  });

  test('sparse MLE equals dense MLE for same values', () => {
    const [R, x0, x1, x2] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);
    const vars = [x0, x1, x2];

    // Create sparse version
    const sparseEntries = new Map([
      [0, F.__call__(9n)],
      [2, F.__call__(5n)],
      [5, F.__call__(2n)],
    ]);
    const sparsePoly = sparseMultilinearExtension(sparseEntries, vars, R);

    // Create dense version (with zeros)
    const denseValues = [
      F.__call__(9n),
      F.__call__(0n),
      F.__call__(5n),
      F.__call__(0n),
      F.__call__(0n),
      F.__call__(2n),
      F.__call__(0n),
      F.__call__(0n),
    ];
    const densePoly = multilinearExtension(denseValues, F, R, vars);

    // Should evaluate to same values at all points
    for (let i = 0; i < 8; i++) {
      const bits = intToBinary(i, 3);
      const point = {
        x0: F.__call__(bits[0]!),
        x1: F.__call__(bits[1]!),
        x2: F.__call__(bits[2]!),
      };
      const sparseVal = sparsePoly.evaluate(point);
      const denseVal = densePoly.evaluate(point);
      expect(sparseVal.eq(denseVal)).toBe(true);
    }
  });

  test('evaluates correctly at non-boolean points', () => {
    // MLE extends to all of F^n, not just {0,1}^n
    const values = [F.__call__(1n), F.__call__(2n), F.__call__(3n), F.__call__(4n)];
    const poly = multilinearExtension(values, F);

    // Test at some non-boolean points
    // The polynomial is linear in each variable, so we can verify manually
    // For 2 variables: f(x,y) = f(0,0)(1-x)(1-y) + f(0,1)(1-x)y + f(1,0)x(1-y) + f(1,1)xy
    //                = 1*(1-x)(1-y) + 2*(1-x)*y + 3*x*(1-y) + 4*x*y
    //                = 1 - x - y + xy + 2y - 2xy + 3x - 3xy + 4xy
    //                = 1 + 2x + y

    // At x=2, y=3: 1 + 2*2 + 3 = 8
    const result = poly.evaluate({ x0: F.__call__(2), x1: F.__call__(3) });
    expect(result.eq(8)).toBe(true);
  });

  test('sum of all eq polynomials equals 1 on hypercube', () => {
    const [R, x0, x1] = MPolynomialRingConstructor(F, ['x0', 'x1']);
    const vars = [x0, x1];

    // Sum of all eq polynomials should be 1 everywhere on the hypercube
    let sumPoly = R.zero();
    for (let i = 0; i < 4; i++) {
      sumPoly = sumPoly.add(eqPolynomial(i, vars, R));
    }

    // Verify at all hypercube points
    for (let i = 0; i < 4; i++) {
      const bits = intToBinary(i, 2);
      const point = { x0: F.__call__(bits[0]!), x1: F.__call__(bits[1]!) };
      expect(sumPoly.evaluate(point).eq(1)).toBe(true);
    }
  });
});

describe('edge cases', () => {
  test('large field characteristic', () => {
    const bigF = new FiniteFieldPrime(
      8444461749428370424248824938781546531375899335154063827935233455917409239041n
    );
    const values = [bigF.__call__(1n), bigF.__call__(2n), bigF.__call__(3n), bigF.__call__(4n)];
    const poly = multilinearExtension(values, bigF);

    expect(poly.evaluate({ x0: bigF.__call__(0), x1: bigF.__call__(0) }).eq(1)).toBe(true);
    expect(poly.evaluate({ x0: bigF.__call__(1), x1: bigF.__call__(1) }).eq(4)).toBe(true);
  });

  test('values with zeros', () => {
    const values = [F.__call__(0n), F.__call__(0n), F.__call__(5n), F.__call__(0n)];
    const poly = multilinearExtension(values, F);

    expect(poly.evaluate({ x0: F.__call__(0), x1: F.__call__(0) }).eq(0)).toBe(true);
    expect(poly.evaluate({ x0: F.__call__(1), x1: F.__call__(0) }).eq(5)).toBe(true);
  });

  test('all-zero values', () => {
    const values = [F.__call__(0n), F.__call__(0n), F.__call__(0n), F.__call__(0n)];
    const poly = multilinearExtension(values, F);

    expect(poly.isZero()).toBe(true);
  });

  test('all-same values gives constant polynomial', () => {
    const values = [F.__call__(7n), F.__call__(7n), F.__call__(7n), F.__call__(7n)];
    const poly = multilinearExtension(values, F);

    expect(poly.isConstant()).toBe(true);
    expect(poly.constantCoefficient().eq(7)).toBe(true);
  });
});
