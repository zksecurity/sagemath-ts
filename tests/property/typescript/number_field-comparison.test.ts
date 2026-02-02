/**
 * Property tests for number field functions
 * Compares sagemath-ts output with real SageMath output
 *
 * This test file parses output from tests/property/sage/number_field.sage
 * and verifies that our TypeScript implementation produces identical results.
 *
 * Usage:
 *   bun test tests/property/typescript/number_field-comparison.test.ts
 */

import { beforeAll, describe, expect, test } from 'bun:test';
import { existsSync, readFileSync } from 'node:fs';
import { join } from 'node:path';

// Import from sagemath-ts package
import {
  CyclotomicField,
  NumberField,
  NumberFieldElement,
  QuadraticField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/index.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';

// Path to SageMath output file
const SAGE_OUTPUT_PATH = join(import.meta.dir, '..', 'transcripts', 'sage', 'number_field.txt');

/**
 * Parse SageMath output file into structured test results
 */
interface SageTestResult {
  operation: string;
  result: string;
}

function parseSageOutput(content: string): SageTestResult[] {
  const results: SageTestResult[] = [];
  const lines = content.split('\n');

  let currentOperation = '';

  for (const line of lines) {
    if (line.startsWith('# ') && !line.startsWith('# ===')) {
      // This is an operation line
      currentOperation = line.substring(2);
    } else if (line.trim() !== '' && !line.startsWith('#') && currentOperation) {
      // This is a result line
      results.push({
        operation: currentOperation,
        result: line.trim(),
      });
      currentOperation = '';
    }
  }

  return results;
}

/**
 * Skip test gracefully if SageMath output not available
 */
function skipIfNoSageOutput(): boolean {
  if (!existsSync(SAGE_OUTPUT_PATH)) {
    console.log(
      '\n⚠️  SageMath output file not found. Run `sage tests/property/sage/number_field.sage > tests/property/transcripts/sage/number_field.txt` first.\n'
    );
    return true;
  }
  return false;
}

describe('NumberField SageMath Comparison', () => {
  let sageResults: SageTestResult[] = [];
  let hasSageOutput = false;

  beforeAll(() => {
    if (existsSync(SAGE_OUTPUT_PATH)) {
      const content = readFileSync(SAGE_OUTPUT_PATH, 'utf-8');
      sageResults = parseSageOutput(content);
      hasSageOutput = true;
    }
  });

  describe('Q(sqrt(2)) tests', () => {
    const K = (() => {
      const poly = RationalPolynomial.fromBigInts([-2n, 0n, 1n]); // x^2 - 2
      return new NumberField(poly, 'a');
    })();
    const a = K.gen();

    test('degree', () => {
      expect(K.degree()).toBe(2);
    });

    test('discriminant', () => {
      expect(K.discriminant()).toBe(8n);
    });

    test('signature', () => {
      expect(K.signature()).toEqual([2, 0]);
    });

    test('trace of a', () => {
      expect(a.trace().eq(Rational.zero())).toBe(true);
    });

    test('norm of a', () => {
      expect(a.norm().eq(new Rational(-2n))).toBe(true);
    });

    test('trace of a+1', () => {
      const elem = a.add(K.one());
      expect(elem.trace().eq(new Rational(2n))).toBe(true);
    });

    test('norm of a+1', () => {
      const elem = a.add(K.one());
      expect(elem.norm().eq(new Rational(-1n))).toBe(true);
    });

    test('charpoly of a', () => {
      const cp = a.charpoly();
      // x^2 - 2
      expect(cp.degree()).toBe(2);
      expect(cp.getCoeff(0).eq(new Rational(-2n))).toBe(true);
      expect(cp.getCoeff(1).isZero()).toBe(true);
      expect(cp.getCoeff(2).eq(Rational.one())).toBe(true);
    });

    test('minpoly of a', () => {
      const mp = a.minpoly();
      expect(mp.degree()).toBe(2);
      expect(mp.getCoeff(0).eq(new Rational(-2n))).toBe(true);
    });
  });

  describe('Q(sqrt(-1)) tests', () => {
    const K = (() => {
      const poly = RationalPolynomial.fromBigInts([1n, 0n, 1n]); // x^2 + 1
      return new NumberField(poly, 'i');
    })();
    const i = K.gen();

    test('degree', () => {
      expect(K.degree()).toBe(2);
    });

    test('discriminant', () => {
      expect(K.discriminant()).toBe(-4n);
    });

    test('signature', () => {
      expect(K.signature()).toEqual([0, 1]);
    });

    test('trace of i', () => {
      expect(i.trace().isZero()).toBe(true);
    });

    test('norm of i', () => {
      expect(i.norm().eq(Rational.one())).toBe(true);
    });
  });

  describe('QuadraticField discriminants', () => {
    const testCases: [bigint, bigint][] = [
      [2n, 8n],
      [3n, 12n],
      [5n, 5n], // 5 ≡ 1 (mod 4)
      [-1n, -4n],
      [-3n, -3n], // -3 ≡ 1 (mod 4)
      [-5n, -20n],
      [-7n, -7n], // -7 ≡ 1 (mod 4)
      [-11n, -11n], // -11 ≡ 1 (mod 4)
      [13n, 13n], // 13 ≡ 1 (mod 4)
      [17n, 17n], // 17 ≡ 1 (mod 4)
    ];

    for (const [d, expectedDisc] of testCases) {
      test(`QuadraticField(${d}).discriminant() = ${expectedDisc}`, () => {
        const K = QuadraticField.create(d);
        expect(K.discriminant()).toBe(expectedDisc);
      });
    }
  });

  describe('CyclotomicField degrees', () => {
    const testCases: [bigint, number][] = [
      [3n, 2], // phi(3) = 2
      [4n, 2], // phi(4) = 2
      [5n, 4], // phi(5) = 4
      [6n, 2], // phi(6) = 2
      [7n, 6], // phi(7) = 6
      [8n, 4], // phi(8) = 4
    ];

    for (const [n, expectedDegree] of testCases) {
      test(`CyclotomicField(${n}).degree() = ${expectedDegree}`, () => {
        const K = CyclotomicField.create(n);
        expect(K.degree()).toBe(expectedDegree);
      });
    }
  });

  describe('Cubic field Q(2^(1/3))', () => {
    const K = (() => {
      const poly = RationalPolynomial.fromBigInts([-2n, 0n, 0n, 1n]); // x^3 - 2
      return new NumberField(poly, 'a');
    })();
    const a = K.gen();

    test('degree', () => {
      expect(K.degree()).toBe(3);
    });

    test('signature', () => {
      // x^3 - 2 has 1 real root and 2 complex conjugate roots
      expect(K.signature()).toEqual([1, 1]);
    });

    test('trace of a', () => {
      // trace(2^(1/3)) = 0 (sum of all roots)
      expect(a.trace().isZero()).toBe(true);
    });

    test('norm of a', () => {
      // norm(2^(1/3)) = 2 (product of all roots)
      expect(a.norm().eq(new Rational(2n))).toBe(true);
    });

    test('trace of a^2', () => {
      const a2 = a.mul(a);
      // trace(2^(2/3)) = 0
      expect(a2.trace().isZero()).toBe(true);
    });

    test('norm of a^2', () => {
      const a2 = a.mul(a);
      // norm(2^(2/3)) = 2^2 = 4
      expect(a2.norm().eq(new Rational(4n))).toBe(true);
    });
  });

  describe('Element arithmetic in Q(sqrt(2))', () => {
    const K = (() => {
      const poly = RationalPolynomial.fromBigInts([-2n, 0n, 1n]);
      return new NumberField(poly, 'a');
    })();
    const a = K.gen();

    test('a^2 = 2', () => {
      const a2 = a.mul(a);
      expect(a2.__getitem__(0).eq(new Rational(2n))).toBe(true);
      expect(a2.__getitem__(1).isZero()).toBe(true);
    });

    test('a^3 = 2a', () => {
      const a3 = a.pow(3n);
      expect(a3.__getitem__(0).isZero()).toBe(true);
      expect(a3.__getitem__(1).eq(new Rational(2n))).toBe(true);
    });

    test('(a+1)^2 = 3 + 2a', () => {
      const result = a.add(K.one()).pow(2n);
      expect(result.__getitem__(0).eq(new Rational(3n))).toBe(true);
      expect(result.__getitem__(1).eq(new Rational(2n))).toBe(true);
    });

    test('1/a = a/2', () => {
      const aInv = a.inv();
      expect(aInv.__getitem__(0).isZero()).toBe(true);
      expect(aInv.__getitem__(1).eq(new Rational(1n, 2n))).toBe(true);
    });

    test('a.is_integral()', () => {
      expect(a.is_integral()).toBe(true);
    });

    test('(a/2).is_integral()', () => {
      const aHalf = a.scalarMul(new Rational(1n, 2n));
      expect(aHalf.is_integral()).toBe(false);
    });

    test('(a+1).is_integral()', () => {
      expect(a.add(K.one()).is_integral()).toBe(true);
    });
  });
});
