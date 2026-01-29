/**
 * Unit tests for sage/crypto/boolean_function
 *
 * Tests for Boolean function implementation including Walsh-Hadamard transform,
 * nonlinearity, algebraic properties, and cryptographic measures.
 */
import { describe, expect, test } from 'bun:test';
import {
  BooleanFunction,
  hammingWeight,
  randomBooleanFunction,
  fromANF,
  createBentFunction,
  createAffineFunction,
  verifyParseval,
} from './boolean_function.js';

describe('BooleanFunction construction', () => {
  test('creates zero function from integer', () => {
    const f = new BooleanFunction(3);
    expect(f.nvariables()).toBe(3);
    expect(f.length()).toBe(8);

    // All values should be 0
    for (let i = 0; i < 8; i++) {
      expect(f.call(i)).toBe(0);
    }
  });

  test('creates from truth table array', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    expect(f.nvariables()).toBe(2);
    expect(f.call(0)).toBe(0);
    expect(f.call(1)).toBe(1);
    expect(f.call(2)).toBe(1);
    expect(f.call(3)).toBe(0);
  });

  test('creates from boolean array', () => {
    const f = new BooleanFunction([false, true, true, false, true, false, false, true]);
    expect(f.nvariables()).toBe(3);
    expect(f.call(0)).toBe(0);
    expect(f.call(1)).toBe(1);
    expect(f.call(7)).toBe(1);
  });

  test('creates from hex string', () => {
    const f = new BooleanFunction('6');
    // '6' = 0110 in binary (LSB first: 0, 1, 1, 0)
    expect(f.nvariables()).toBe(2);
    expect(f.truthTable('int')).toEqual([0, 1, 1, 0]);
  });

  test('creates from longer hex string', () => {
    const f = new BooleanFunction('0f');
    // '0f' = 0000 1111 (16 bits -> 4 variables)
    // LSB first: 0,0,0,0, 1,1,1,1
    expect(f.nvariables()).toBe(3);
    expect(f.truthTable('int')).toEqual([0, 0, 0, 0, 1, 1, 1, 1]);
  });

  test('copy constructor', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    const g = new BooleanFunction(f);
    expect(g.nvariables()).toBe(f.nvariables());
    expect(g.truthTable('int')).toEqual(f.truthTable('int'));

    // Modifying g should not affect f
    g.set(0, 1);
    expect(f.call(0)).toBe(0);
    expect(g.call(0)).toBe(1);
  });

  test('throws on non-power-of-2 array length', () => {
    expect(() => new BooleanFunction([0, 1, 1])).toThrow(
      'the length of the truth table must be a power of 2'
    );
  });

  test('throws on empty array', () => {
    expect(() => new BooleanFunction([])).toThrow(
      'the length of the truth table must be a power of 2'
    );
  });

  test('throws on invalid hex string length', () => {
    // 'xyz' has length 3, which is not a power of 2
    expect(() => new BooleanFunction('xyz')).toThrow('the length of the truth table must be a power of 2');
  });

  test('throws on invalid hex character', () => {
    // 'xx' has length 2 (power of 2), but contains invalid hex characters
    expect(() => new BooleanFunction('xx')).toThrow('invalid hexadecimal character');
  });
});

describe('BooleanFunction evaluation', () => {
  test('call with integer input', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    expect(f.call(0)).toBe(0);
    expect(f.call(1)).toBe(1);
    expect(f.call(2)).toBe(1);
    expect(f.call(3)).toBe(0);
    expect(f.call(4)).toBe(1);
    expect(f.call(5)).toBe(0);
    expect(f.call(6)).toBe(0);
    expect(f.call(7)).toBe(1);
  });

  test('call with bit array input (LSB first)', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    // [1, 0, 0] represents 1 (binary 001)
    expect(f.call([1, 0, 0])).toBe(1);
    // [0, 1, 0] represents 2 (binary 010)
    expect(f.call([0, 1, 0])).toBe(1);
    // [1, 1, 1] represents 7 (binary 111)
    expect(f.call([1, 1, 1])).toBe(1);
  });

  test('throws on out-of-range input', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    expect(() => f.call(4)).toThrow('index out of bound');
    expect(() => f.call(-1)).toThrow('index out of bound');
  });

  test('throws on wrong number of bits', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    expect(() => f.call([1, 0])).toThrow('expected 3 input bits, got 2');
  });

  test('get and set methods', () => {
    const f = new BooleanFunction(3);
    expect(f.get(0)).toBe(false);

    f.set(1, 1);
    f.set(3, true);
    expect(f.get(1)).toBe(true);
    expect(f.get(3)).toBe(true);
    expect(f.call(1)).toBe(1);
    expect(f.call(3)).toBe(1);
  });
});

describe('truth table formats', () => {
  test('bool format', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    expect(f.truthTable('bool')).toEqual([false, true, true, false]);
  });

  test('int format', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    expect(f.truthTable('int')).toEqual([0, 1, 1, 0]);
  });

  test('hex format', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    expect(f.truthTable('hex')).toBe('6');
  });

  test('hex format roundtrip', () => {
    const original = '0113077c165e76a8';
    const f = new BooleanFunction(original);
    expect(f.truthTable('hex')).toBe(original);
  });

  test('throws on unknown format', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    // @ts-expect-error - testing invalid format
    expect(() => f.truthTable('oct')).toThrow('unknown output format');
  });
});

describe('Walsh-Hadamard transform', () => {
  test('basic Walsh transform', () => {
    // Test with SageMath reference: BooleanFunction([1,0,0,1]) -> (0, 0, 0, -4)
    const f = new BooleanFunction([1, 0, 0, 1]);
    const W = f.walshHadamardTransform();
    expect(W).toEqual([0, 0, 0, -4]);
  });

  test('Walsh transform of XOR function', () => {
    // XOR function: f(x) = x0 XOR x1, truth table [0, 1, 1, 0]
    // This is the complement of [1, 0, 0, 1], so Walsh values are negated
    const f = new BooleanFunction([0, 1, 1, 0]);
    const W = f.walshHadamardTransform();
    expect(W).toEqual([0, 0, 0, 4]);
  });

  test('Walsh transform of zero function', () => {
    const f = new BooleanFunction(3);
    const W = f.walshHadamardTransform();
    // Zero function: all (-1)^0 = 1, so WHT is [8, 0, 0, ..., 0]
    expect(W[0]).toBe(8);
    for (let i = 1; i < 8; i++) {
      expect(W[i]).toBe(0);
    }
  });

  test('Walsh transform of all-ones function', () => {
    const f = new BooleanFunction([1, 1, 1, 1, 1, 1, 1, 1]);
    const W = f.walshHadamardTransform();
    // All-ones function: all (-1)^1 = -1, so WHT is [-8, 0, 0, ..., 0]
    expect(W[0]).toBe(-8);
    for (let i = 1; i < 8; i++) {
      expect(W[i]).toBe(0);
    }
  });

  test("Parseval's theorem holds", () => {
    // For any Boolean function, sum of W_f(a)^2 = 2^{2n}
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    const W = f.walshHadamardTransform();

    let sumSquares = 0;
    for (const w of W) {
      sumSquares += w * w;
    }

    expect(sumSquares).toBe(64); // 2^{2*3} = 64
  });

  test("Parseval's theorem for random functions", () => {
    for (let trial = 0; trial < 5; trial++) {
      const f = randomBooleanFunction(4);
      expect(verifyParseval(f.walshHadamardTransform())).toBe(true);
    }
  });

  test('absolute Walsh spectrum', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    const spectrum = f.absoluteWalshSpectrum();

    // Should map absolute values to their frequencies
    expect(spectrum.has(0) || spectrum.has(4) || spectrum.has(8)).toBe(true);

    // Sum of frequencies should equal 2^n
    let total = 0;
    for (const count of spectrum.values()) {
      total += count;
    }
    expect(total).toBe(8);
  });

  test('Walsh transform caching', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    const W1 = f.walshHadamardTransform();
    const W2 = f.walshHadamardTransform();

    // Should return same values
    expect(W1).toEqual(W2);
  });
});

describe('balanced functions', () => {
  test('balanced XOR function', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    expect(f.isBalanced()).toBe(true);
  });

  test('unbalanced function', () => {
    const f = new BooleanFunction([0, 0, 1, 1, 1, 1, 1, 1]);
    expect(f.isBalanced()).toBe(false);
  });

  test('zero function is not balanced', () => {
    const f = new BooleanFunction(3);
    expect(f.isBalanced()).toBe(false);
  });

  test('balanced is equivalent to W_f(0) = 0', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    const W = f.walshHadamardTransform();
    expect(f.isBalanced()).toBe(W[0] === 0);
  });
});

describe('nonlinearity', () => {
  test('nonlinearity of XOR function', () => {
    // XOR is an affine function, so nonlinearity = 0
    const xor = new BooleanFunction([0, 1, 1, 0]);
    expect(xor.nonlinearity()).toBe(0);
  });

  test('nonlinearity formula: NL = 2^{n-1} - max|W|/2', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    const W = f.walshHadamardTransform();
    const maxAbs = Math.max(...W.map(Math.abs));
    const expected = (4 - maxAbs / 2); // 2^{n-1} = 4 for n=3
    expect(f.nonlinearity()).toBe(expected);
  });

  test('nonlinearity caching', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    const nl1 = f.nonlinearity();
    const nl2 = f.nonlinearity();
    expect(nl1).toBe(nl2);
  });
});

describe('bent functions', () => {
  test('bent function creation', () => {
    const bent = createBentFunction(4);
    expect(bent.nvariables()).toBe(4);
    expect(bent.isBent()).toBe(true);
  });

  test('bent function has maximum nonlinearity', () => {
    const bent = createBentFunction(4);
    // Max nonlinearity for n=4 is 2^3 - 2^1 = 6
    expect(bent.nonlinearity()).toBe(6);
  });

  test('bent functions only exist for even n', () => {
    expect(() => createBentFunction(3)).toThrow('bent functions require an even number');
  });

  test('bent function from known hex string', () => {
    // Known bent function on 6 variables
    const f = new BooleanFunction('0113077C165E76A8');
    expect(f.nvariables()).toBe(6);
    expect(f.isBent()).toBe(true);
    // Max nonlinearity for n=6 is 2^5 - 2^2 = 28
    expect(f.nonlinearity()).toBe(28);
  });

  test('non-bent function correctly identified', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    expect(f.isBent()).toBe(false);
  });

  test('odd-variable functions are never bent', () => {
    const f = randomBooleanFunction(3);
    expect(f.isBent()).toBe(false);
  });
});

describe('correlation immunity', () => {
  test('correlation immunity of constant function', () => {
    const f = new BooleanFunction(3); // Zero function
    // Zero function has max correlation immunity
    expect(f.correlationImmunity()).toBe(2);
  });

  test('correlation immunity basic computation', () => {
    const f = new BooleanFunction('7969817CC5893BA6AC326E47619F5AD0');
    expect(f.correlationImmunity()).toBe(2);
  });

  test('resiliency order', () => {
    // Resiliency = balanced + correlation immune
    const f = new BooleanFunction('077CE5A2F8831A5DF8831A5D077CE5A26996699669699696669999665AA5A55A');
    const order = f.resiliencyOrder();
    expect(order).toBe(3);
  });

  test('unbalanced function has resiliency -1', () => {
    const f = new BooleanFunction([0, 0, 1, 1, 1, 1, 1, 1]);
    expect(f.resiliencyOrder()).toBe(-1);
  });
});

describe('autocorrelation', () => {
  test('autocorrelation computation', () => {
    const f = new BooleanFunction('03');
    const ac = f.autocorrelation();
    expect(ac[0]).toBe(8); // Delta(0) = 2^n
    expect(ac[1]).toBe(8);
  });

  test('autocorrelation Delta(0) = 2^n always', () => {
    const f = randomBooleanFunction(4);
    const ac = f.autocorrelation();
    expect(ac[0]).toBe(16); // 2^4
  });

  test('absolute autocorrelation spectrum', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    const spectrum = f.absoluteAutocorrelation();

    let total = 0;
    for (const count of spectrum.values()) {
      total += count;
    }
    expect(total).toBe(8);
  });

  test('absolute indicator', () => {
    const f = new BooleanFunction('7969817CC5893BA6AC326E47619F5AD0');
    const ai = f.absoluteIndicator();
    expect(ai).toBe(32);
  });

  test('sum of square indicator', () => {
    const f = new BooleanFunction('7969817CC5893BA6AC326E47619F5AD0');
    const ssi = f.sumOfSquareIndicator();
    expect(ssi).toBe(32768);
  });
});

describe('derivative', () => {
  test('derivative computation', () => {
    const f = new BooleanFunction([0, 1, 0, 1, 0, 1, 0, 1]);
    const df = f.derivative(1);
    // D_1 f(x) = f(x) + f(x XOR 1)
    // f is x0, so D_1 f = x0 + (x0 XOR 1) = x0 + x0 + 1 = 1
    const expected = [1, 1, 1, 1, 1, 1, 1, 1];
    expect(df.truthTable('int')).toEqual(expected);
  });

  test('derivative with bit array direction', () => {
    const f = new BooleanFunction([0, 1, 0, 1, 0, 1, 0, 1]);
    const df = f.derivative([1, 0, 0]); // Direction 1
    const expected = [1, 1, 1, 1, 1, 1, 1, 1];
    expect(df.truthTable('int')).toEqual(expected);
  });

  test('derivative of zero function is zero', () => {
    const f = new BooleanFunction(3);
    const df = f.derivative(1);
    expect(df.truthTable('int')).toEqual([0, 0, 0, 0, 0, 0, 0, 0]);
  });

  test('throws on invalid direction', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    expect(() => f.derivative(5)).toThrow('index out of bound');
    expect(() => f.derivative([1])).toThrow('expected 2 bits, got 1');
  });
});

describe('linear structures', () => {
  test('detect linear structure', () => {
    const f = new BooleanFunction([0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0]);
    expect(f.isLinearStructure(1)).toBe(true);
  });

  test('detect non-linear structure', () => {
    const f = new BooleanFunction([0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0]);
    expect(f.isLinearStructure(7)).toBe(false);
  });

  test('has linear structure', () => {
    const f = new BooleanFunction([0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0]);
    expect(f.hasLinearStructure()).toBe(true);
  });

  test('no linear structure', () => {
    const f = new BooleanFunction([0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1]);
    expect(f.hasLinearStructure()).toBe(false);
  });
});

describe('plateaued functions', () => {
  test('bent function is plateaued', () => {
    const bent = createBentFunction(4);
    expect(bent.isPlateaued()).toBe(true);
  });

  test('detect plateaued function', () => {
    const f = new BooleanFunction('0113077C165E76A8');
    expect(f.isPlateaued()).toBe(true);
  });
});

describe('algebraic normal form', () => {
  test('ANF coefficients of XOR', () => {
    // f = x0 XOR x1 has ANF x0 + x1
    const f = new BooleanFunction([0, 1, 1, 0]);
    const anf = f.algebraicNormalFormCoefficients();
    // Index 1 = x0, index 2 = x1
    expect(anf[0]).toBe(0); // constant term
    expect(anf[1]).toBe(1); // x0
    expect(anf[2]).toBe(1); // x1
    expect(anf[3]).toBe(0); // x0*x1
  });

  test('ANF coefficients roundtrip', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 1, 1]);
    const anf = f.algebraicNormalFormCoefficients();
    const g = fromANF(anf);
    expect(g.truthTable('int')).toEqual(f.truthTable('int'));
  });

  test('ANF string representation', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    const anfStr = f.algebraicNormalForm();
    expect(anfStr).toBe('x0 + x1');
  });

  test('ANF of constant 1', () => {
    const f = new BooleanFunction([1, 1]);
    expect(f.algebraicNormalForm()).toBe('1');
  });

  test('ANF of zero function', () => {
    const f = new BooleanFunction(3);
    expect(f.algebraicNormalForm()).toBe('0');
  });

  test('ANF with custom variable names', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    expect(f.algebraicNormalForm(['a', 'b'])).toBe('a + b');
  });

  test('ANF with product terms', () => {
    const f = new BooleanFunction([0, 0, 0, 1, 0, 0, 0, 1]);
    // This is x0*x1
    expect(f.algebraicNormalForm()).toContain('x0*x1');
  });
});

describe('algebraic degree', () => {
  test('algebraic degree of linear function', () => {
    const f = new BooleanFunction([0, 1, 1, 0]); // x0 + x1
    expect(f.algebraicDegree()).toBe(1);
  });

  test('algebraic degree of zero function', () => {
    const f = new BooleanFunction(3);
    expect(f.algebraicDegree()).toBe(-1);
  });

  test('algebraic degree of constant 1', () => {
    const f = new BooleanFunction([1, 1, 1, 1]);
    expect(f.algebraicDegree()).toBe(0);
  });

  test('algebraic degree of quadratic function', () => {
    const f = new BooleanFunction([0, 0, 0, 1, 0, 0, 0, 1]);
    expect(f.algebraicDegree()).toBe(2); // Contains x0*x1
  });

  test('algebraic degree of cubic function', () => {
    // Create a function with x0*x1*x2 term
    const f = new BooleanFunction([0, 0, 0, 0, 0, 0, 0, 1]);
    expect(f.algebraicDegree()).toBe(3);
  });
});

describe('algebraic immunity', () => {
  test('algebraic immunity of simple function', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    const ai = f.algebraicImmunity();
    expect(ai).toBeGreaterThanOrEqual(0);
    expect(ai).toBeLessThanOrEqual(3);
  });

  test('zero function has AI = 0', () => {
    const f = new BooleanFunction(2);
    expect(f.algebraicImmunity()).toBe(0);
  });
});

describe('Boolean function operations', () => {
  test('complement', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    const g = f.complement();
    expect(g.truthTable('int')).toEqual([1, 0, 0, 1]);
  });

  test('add (XOR)', () => {
    const f = new BooleanFunction([0, 1, 0, 1]);
    const g = new BooleanFunction([0, 0, 1, 1]);
    const h = f.add(g);
    expect(h.truthTable('int')).toEqual([0, 1, 1, 0]); // XOR
  });

  test('mul (AND)', () => {
    const f = new BooleanFunction([0, 1, 0, 1]);
    const g = new BooleanFunction([0, 0, 1, 1]);
    const h = f.mul(g);
    expect(h.truthTable('int')).toEqual([0, 0, 0, 1]); // AND
  });

  test('concatenate', () => {
    const f = new BooleanFunction([0, 1, 0, 1]);
    const g = new BooleanFunction([1, 0, 1, 0]);
    const h = f.concatenate(g);
    expect(h.nvariables()).toBe(3);
    expect(h.truthTable('int')).toEqual([0, 1, 0, 1, 1, 0, 1, 0]);
  });

  test('operations require same number of variables', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    const g = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);

    expect(() => f.add(g)).toThrow('same number of variables');
    expect(() => f.mul(g)).toThrow('same number of variables');
    expect(() => f.concatenate(g)).toThrow('same number of variables');
  });

  test('equals', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    const g = new BooleanFunction([0, 1, 1, 0]);
    const h = new BooleanFunction([0, 1, 0, 1]);

    expect(f.equals(g)).toBe(true);
    expect(f.equals(h)).toBe(false);
  });
});

describe('symmetric functions', () => {
  test('constant function is symmetric', () => {
    const f = new BooleanFunction(3);
    expect(f.isSymmetric()).toBe(true);
  });

  test('XOR is symmetric', () => {
    // XOR on 2 variables: f(x0, x1) = x0 + x1
    // Output depends only on Hamming weight: wt=0 -> 0, wt=1 -> 1, wt=2 -> 0
    // So XOR IS symmetric!
    const f = new BooleanFunction([0, 1, 1, 0]);
    expect(f.isSymmetric()).toBe(true);
  });

  test('non-symmetric function', () => {
    // f(x0, x1, x2) = x0 (depends on specific bit, not weight)
    const f = new BooleanFunction([0, 1, 0, 1, 0, 1, 0, 1]);
    expect(f.isSymmetric()).toBe(false);
  });

  test('majority function is symmetric', () => {
    // Majority on 3 bits: output 1 iff at least 2 inputs are 1
    // Weights: 0->0, 1->0, 2->1, 3->1
    const f = new BooleanFunction([0, 0, 0, 1, 0, 1, 1, 1]);
    expect(f.isSymmetric()).toBe(true);
  });
});

describe('helper functions', () => {
  test('hammingWeight', () => {
    expect(hammingWeight(0)).toBe(0);
    expect(hammingWeight(1)).toBe(1);
    expect(hammingWeight(2)).toBe(1);
    expect(hammingWeight(3)).toBe(2);
    expect(hammingWeight(7)).toBe(3);
    expect(hammingWeight(255)).toBe(8);
  });

  test('randomBooleanFunction', () => {
    const f = randomBooleanFunction(5);
    expect(f.nvariables()).toBe(5);
    expect(f.length()).toBe(32);

    // Random function should not typically be all zeros
    // (probability 2^{-32})
    const allZero = (f.truthTable('int') as number[]).every((x) => x === 0);
    // Very unlikely to be all zeros
    expect(allZero).toBe(false);
  });

  test('fromANF', () => {
    // ANF: x0 + x1 (indices 1 and 2)
    const f = fromANF([0, 1, 1, 0]);
    expect(f.truthTable('int')).toEqual([0, 1, 1, 0]);
  });

  test('createAffineFunction', () => {
    // f(x) = x0 + x2 + 1
    const f = createAffineFunction(0b101, 1, 3);
    expect(f.call(0)).toBe(1); // 0 + 0 + 1 = 1
    expect(f.call(1)).toBe(0); // 1 + 0 + 1 = 0
    expect(f.call(4)).toBe(0); // 0 + 1 + 1 = 0
    expect(f.call(5)).toBe(1); // 1 + 1 + 1 = 1
  });

  test('verifyParseval', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);
    expect(verifyParseval(f.walshHadamardTransform())).toBe(true);
  });
});

describe('string representation', () => {
  test('toString', () => {
    const f = new BooleanFunction(4);
    expect(f.toString()).toBe('Boolean function with 4 variables');
  });

  test('toString singular', () => {
    const f = new BooleanFunction(1);
    expect(f.toString()).toBe('Boolean function with 1 variable');
  });
});

describe('iteration', () => {
  test('iterate over truth table', () => {
    const f = new BooleanFunction([0, 1, 1, 0]);
    const values: boolean[] = [];
    for (const v of f) {
      values.push(v);
    }
    expect(values).toEqual([false, true, true, false]);
  });

  test('spread operator', () => {
    const f = new BooleanFunction([1, 0, 1, 0]);
    const values = [...f];
    expect(values).toEqual([true, false, true, false]);
  });
});

describe('cryptographic properties verification', () => {
  test('affine functions have nonlinearity 0', () => {
    // Test various affine functions
    for (let a = 0; a < 8; a++) {
      for (const b of [0, 1] as const) {
        const f = createAffineFunction(a, b, 3);
        expect(f.nonlinearity()).toBe(0);
      }
    }
  });

  test('bent function properties', () => {
    const bent = createBentFunction(6);

    // Bent functions are always plateaued
    expect(bent.isPlateaued()).toBe(true);

    // Bent functions are NOT balanced
    expect(bent.isBalanced()).toBe(false);

    // Max nonlinearity for n=6: 2^5 - 2^2 = 28
    expect(bent.nonlinearity()).toBe(28);
  });

  test('Walsh spectrum of bent function', () => {
    const bent = createBentFunction(4);
    const W = bent.walshHadamardTransform();

    // All Walsh values have absolute value 2^{n/2} = 4
    for (const w of W) {
      expect(Math.abs(w)).toBe(4);
    }
  });

  test('random function nonlinearity distribution', () => {
    // Random functions on 8 variables typically have nonlinearity
    // around 100-112 (out of max 120)
    const nonlinearities: number[] = [];
    for (let i = 0; i < 10; i++) {
      const f = randomBooleanFunction(8);
      nonlinearities.push(f.nonlinearity());
    }

    // Average should be reasonably high
    const avg = nonlinearities.reduce((a, b) => a + b, 0) / nonlinearities.length;
    expect(avg).toBeGreaterThan(90);
    expect(avg).toBeLessThan(115);
  });
});

describe('edge cases', () => {
  test('1-variable function', () => {
    const f = new BooleanFunction([0, 1]);
    expect(f.nvariables()).toBe(1);
    // f(0)=0, f(1)=1 means values are [1, -1]
    // Walsh: W(0) = 1 + (-1) = 0, W(1) = 1 - (-1) = 2
    expect(f.walshHadamardTransform()).toEqual([0, 2]);
    expect(f.isBalanced()).toBe(true);
    expect(f.nonlinearity()).toBe(0); // Linear function
  });

  test('single value truth table', () => {
    const f0 = new BooleanFunction([0]);
    const f1 = new BooleanFunction([1]);

    expect(f0.nvariables()).toBe(0);
    expect(f0.call(0)).toBe(0);
    expect(f0.isBalanced()).toBe(false);

    expect(f1.nvariables()).toBe(0);
    expect(f1.call(0)).toBe(1);
    expect(f1.isBalanced()).toBe(false);
  });

  test('large function (8 variables)', () => {
    const f = randomBooleanFunction(8);
    expect(f.nvariables()).toBe(8);
    expect(f.length()).toBe(256);

    // Should be able to compute all properties
    const W = f.walshHadamardTransform();
    expect(W.length).toBe(256);
    expect(verifyParseval(W)).toBe(true);

    const nl = f.nonlinearity();
    expect(nl).toBeGreaterThanOrEqual(0);
    expect(nl).toBeLessThanOrEqual(120); // Max for 8 variables
  });

  test('cache clearing on modification', () => {
    const f = new BooleanFunction([0, 1, 1, 0, 1, 0, 0, 1]);

    // Compute Walsh transform (caches it)
    const W1 = f.walshHadamardTransform();

    // Modify function
    f.set(0, 1);

    // Walsh transform should be recomputed
    const W2 = f.walshHadamardTransform();
    expect(W1).not.toEqual(W2);
  });
});
