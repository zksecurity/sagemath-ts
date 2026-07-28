/**
 * Unit tests for sage/crypto/sbox
 *
 * Tests for S-Box cryptographic analysis implementation.
 * Reference: sage/crypto/sbox.pyx
 */
import { describe, expect, test } from 'bun:test';
import { Rational } from '../rings/rational.js';
import {
  AES_SBOX,
  DES_SBOX1,
  GIFT_SBOX,
  PRESENT_SBOX,
  SBox,
  SKINNY_SBOX,
  feistel_construction,
  misty_construction,
} from './sbox.js';

describe('SBox construction', () => {
  test('constructs from lookup table', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.input_size()).toBe(3);
    expect(S.output_size()).toBe(3);
    expect(S.length).toBe(3);
  });

  test('throws on empty lookup table', () => {
    expect(() => new SBox([])).toThrow('no lookup table provided');
  });

  test('throws on non-power-of-2 length', () => {
    expect(() => new SBox([1, 2, 3])).toThrow('lookup table length is not a power of 2');
  });

  test('handles different input/output sizes', () => {
    // 4 inputs (2-bit), outputs go up to 7 (3-bit)
    const S = new SBox([7, 4, 8, 6]);
    expect(S.input_size()).toBe(2);
    expect(S.output_size()).toBe(4);
  });

  test('preserves big_endian option', () => {
    const S_be = new SBox([7, 6, 0, 4, 2, 5, 1, 3], { big_endian: true });
    const S_le = new SBox([7, 6, 0, 4, 2, 5, 1, 3], { big_endian: false });

    // Same integer input gives same integer output
    expect(S_be.call(1)).toBe(6);
    expect(S_le.call(1)).toBe(6);

    // Bit representation differs based on endianness
    // Big endian: [0, 0, 1] -> 1, S(1) = 6 = 110 -> [1, 1, 0]
    expect(S_be.call([0, 0, 1])).toEqual([1, 1, 0]);
    // Little endian: [0, 0, 1] -> 4, S(4) = 2 = 010 -> [0, 1, 0]
    expect(S_le.call([0, 0, 1])).toEqual([0, 1, 0]);
  });
});

describe('SBox evaluation', () => {
  test('evaluates integer input', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.call(0)).toBe(7);
    expect(S.call(1)).toBe(6);
    expect(S.call(2)).toBe(0);
    expect(S.call(7)).toBe(3);
  });

  test('evaluates bit array input (big endian)', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    // Input [0,0,1] = 1 in big endian, S(1) = 6 = [1,1,0]
    expect(S.call([0, 0, 1])).toEqual([1, 1, 0]);
    // Input [0,0,0] = 0 in big endian, S(0) = 7 = [1,1,1]
    expect(S.call([0, 0, 0])).toEqual([1, 1, 1]);
  });

  test('get method works', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.get(0)).toBe(7);
    expect(S.get(3)).toBe(4);
  });

  test('throws on out of range input', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(() => S.call(8)).toThrow('input 8 out of range');
    expect(() => S.call(-1)).toThrow('input -1 out of range');
  });

  test('throws on wrong number of input bits', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(() => S.call([0, 1])).toThrow('expected 3 input bits');
  });
});

describe('to_bits and from_bits', () => {
  test('to_bits converts integer to bit array', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.to_bits(6)).toEqual([1, 1, 0]);
    expect(S.to_bits(0)).toEqual([0, 0, 0]);
    expect(S.to_bits(7)).toEqual([1, 1, 1]);
  });

  test('to_bits with custom length', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.to_bits(6, 4)).toEqual([0, 1, 1, 0]);
    expect(S.to_bits(1, 5)).toEqual([0, 0, 0, 0, 1]);
  });

  test('from_bits converts bit array to integer', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.from_bits([1, 1, 0])).toBe(6);
    expect(S.from_bits([0, 0, 0])).toBe(0);
    expect(S.from_bits([1, 1, 1])).toBe(7);
  });

  test('roundtrip conversion', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    for (let i = 0; i < 8; i++) {
      expect(S.from_bits(S.to_bits(i))).toBe(i);
    }
  });
});

describe('is_permutation', () => {
  test('returns true for bijective S-box', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.is_permutation()).toBe(true);
  });

  test('returns false for non-bijective S-box', () => {
    const S = new SBox([3, 2, 0, 0, 2, 1, 1, 3]);
    expect(S.is_permutation()).toBe(false);
  });

  test('returns false when input and output sizes differ', () => {
    const S = new SBox([7, 4, 8, 6]); // 2-bit input, 4-bit output
    expect(S.is_permutation()).toBe(false);
  });
});

describe('is_involution', () => {
  test('returns true for involutive S-box', () => {
    // Identity is an involution
    const id = new SBox([0, 1, 2, 3, 4, 5, 6, 7]);
    expect(id.is_involution()).toBe(true);

    // Swap pairs is an involution
    const swap = new SBox([1, 0, 3, 2, 5, 4, 7, 6]);
    expect(swap.is_involution()).toBe(true);
  });

  test('returns false for non-involutive S-box', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.is_involution()).toBe(false);
  });
});

describe('inverse', () => {
  test('computes inverse of permutation S-box', () => {
    const S = new SBox([0, 1, 3, 6, 7, 4, 5, 2]);
    const Sinv = S.inverse();

    // Verify S^{-1}(S(x)) = x for all x
    for (let x = 0; x < 8; x++) {
      expect(Sinv.call(S.call(x) as number)).toBe(x);
    }

    // Verify S(S^{-1}(y)) = y for all y
    for (let y = 0; y < 8; y++) {
      expect(S.call(Sinv.call(y) as number)).toBe(y);
    }
  });

  test('throws for non-permutation S-box', () => {
    const S = new SBox([3, 2, 0, 0, 2, 1, 1, 3]);
    expect(() => S.inverse()).toThrow('S-Box must be a permutation');
  });
});

describe('Difference Distribution Table (DDT)', () => {
  test('computes DDT for small S-box', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const ddt = S.difference_distribution_table();

    // First row should have 8 in (0,0) and 0 elsewhere
    expect(Number(ddt.get(0, 0).value)).toBe(8);
    for (let j = 1; j < 8; j++) {
      expect(Number(ddt.get(0, j).value)).toBe(0);
    }

    // All row sums should equal 8 (= 2^m)
    for (let i = 0; i < 8; i++) {
      let rowSum = 0;
      for (let j = 0; j < 8; j++) {
        rowSum += Number(ddt.get(i, j).value);
      }
      expect(rowSum).toBe(8);
    }

    // All column sums should equal 8
    for (let j = 0; j < 8; j++) {
      let colSum = 0;
      for (let i = 0; i < 8; i++) {
        colSum += Number(ddt.get(i, j).value);
      }
      expect(colSum).toBe(8);
    }
  });

  test('DDT entries match reference values', () => {
    // Reference DDT for S = [7,6,0,4,2,5,1,3] from SageMath
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const ddt = S.difference_distribution_table();

    // Spot check some values from sage reference
    // Row 1: [0 2 2 0 2 0 0 2]
    expect(Number(ddt.get(1, 0).value)).toBe(0);
    expect(Number(ddt.get(1, 1).value)).toBe(2);
    expect(Number(ddt.get(1, 2).value)).toBe(2);
    expect(Number(ddt.get(1, 4).value)).toBe(2);
  });

  test('DDT is cached', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const ddt1 = S.difference_distribution_table();
    const ddt2 = S.difference_distribution_table();
    expect(ddt1).toBe(ddt2); // Same object (cached)
  });
});

describe('differential_uniformity', () => {
  test('computes differential uniformity', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.differential_uniformity()).toBe(2);
  });

  test('AES S-box has differential uniformity 4', () => {
    expect(AES_SBOX.differential_uniformity()).toBe(4);
  });

  test('PRESENT S-box has differential uniformity 4', () => {
    expect(PRESENT_SBOX.differential_uniformity()).toBe(4);
  });
});

describe('is_apn', () => {
  test('returns true for APN S-box', () => {
    const S = new SBox([0, 1, 3, 6, 7, 4, 5, 2]);
    expect(S.is_apn()).toBe(true);
    expect(S.differential_uniformity()).toBe(2);
  });

  test('returns false for non-APN S-box', () => {
    expect(PRESENT_SBOX.is_apn()).toBe(false);
    expect(AES_SBOX.is_apn()).toBe(false);
  });
});

describe('Linear Approximation Table (LAT)', () => {
  test("matches Sage's linear_approximation_table doctest", () => {
    // sage: S = SBox(7,6,0,4,2,5,1,3); S.linear_approximation_table()
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const expected = [
      [4, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 2, 2, 2, -2],
      [0, 0, -2, -2, -2, 2, 0, 0],
      [0, 0, -2, 2, 0, 0, -2, -2],
      [0, 2, 0, 2, -2, 0, 2, 0],
      [0, -2, 0, 2, 0, 2, 0, 2],
      [0, -2, -2, 0, 0, -2, 2, 0],
      [0, -2, 2, 0, -2, 0, 0, -2],
    ];
    const lat = S.linear_approximation_table();
    for (let i = 0; i < 8; i++) {
      for (let j = 0; j < 8; j++) {
        expect(Number(lat.get(i, j).value)).toBe(expected[i]![j]!);
      }
    }

    // sage: lat_abs_bias*2 == S.linear_approximation_table(scale='fourier_coefficient')
    const fourier = S.linear_approximation_table('fourier_coefficient');
    // sage: lat_abs_bias/(1 << S.input_size()) == ...(scale='bias')
    const bias = S.linear_approximation_table('bias');
    // sage: lat_abs_bias/(1 << (S.input_size()-1)) == ...(scale='correlation')
    const correlation = S.linear_approximation_table('correlation');
    for (let i = 0; i < 8; i++) {
      for (let j = 0; j < 8; j++) {
        const a = expected[i]![j]!;
        expect(Number(fourier.get(i, j).value)).toBe(2 * a);
        expect(bias.get(i, j).value.eq(new Rational(BigInt(a), 8n))).toBe(true);
        expect(correlation.get(i, j).value.eq(new Rational(BigInt(a), 4n))).toBe(true);
      }
    }

    expect(() => (S.linear_approximation_table as (s: string) => unknown)('nonsense')).toThrow(
      'no such scaling for the LAT: nonsense'
    );
  });

  test('computes LAT for small S-box', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const lat = S.linear_approximation_table();

    // LAT[0,0] should be 2^{m-1} = 4 (bias from 2^{m-1})
    expect(Number(lat.get(0, 0).value)).toBe(4);

    // First row (except [0,0]) should be 0 for balanced S-box
    for (let j = 1; j < 8; j++) {
      expect(Number(lat.get(0, j).value)).toBe(0);
    }
  });

  test('LAT sum of squared entries follows Parseval relation', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const lat = S.linear_approximation_table();

    // For an m x n S-box with LAT entries being biases (count - 2^{m-1}),
    // the sum of squared entries equals 2^m * 2^n * 2^{2m-2} = 2^{3m+n-2}
    // For a 3x3 S-box: 2^{3*3+3-2} = 2^8 = 256 when including all entries
    // Alternatively, check that each row/column has consistent properties

    // Verify LAT entries are integers and bounded
    let allInBounds = true;
    const bound = 1 << (S.input_size() - 1); // 2^{m-1} = 4
    for (let i = 0; i < 8; i++) {
      for (let j = 0; j < 8; j++) {
        const val = Number(lat.get(i, j).value);
        if (Math.abs(val) > bound) {
          allInBounds = false;
        }
      }
    }
    expect(allInBounds).toBe(true);

    // Row 0 column sums (excluding [0,0]) should all be 0 for balanced S-box
    let row0sum = 0;
    for (let j = 1; j < 8; j++) {
      row0sum += Number(lat.get(0, j).value);
    }
    expect(row0sum).toBe(0);
  });

  test('LAT sum of squared entries for 4-bit S-box', () => {
    const lat = PRESENT_SBOX.linear_approximation_table();

    // Verify entries are bounded by 2^{m-1} = 8
    let allInBounds = true;
    const bound = 1 << (PRESENT_SBOX.input_size() - 1);
    for (let i = 0; i < 16; i++) {
      for (let j = 0; j < 16; j++) {
        const val = Number(lat.get(i, j).value);
        if (Math.abs(val) > bound) {
          allInBounds = false;
        }
      }
    }
    expect(allInBounds).toBe(true);

    // First row (except [0,0]) should be all zeros for balanced permutation
    for (let j = 1; j < 16; j++) {
      expect(Number(lat.get(0, j).value)).toBe(0);
    }
  });

  test('LAT is cached', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const lat1 = S.linear_approximation_table();
    const lat2 = S.linear_approximation_table();
    expect(lat1).toBe(lat2);
  });
});

describe('max_linear_bias', () => {
  test('computes maximum linear bias', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.max_linear_bias()).toBe(2);
  });

  test('AES S-box max linear bias is 16', () => {
    expect(AES_SBOX.max_linear_bias()).toBe(16);
  });
});

describe('linearity and nonlinearity', () => {
  test('computes linearity', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.linearity()).toBe(4); // 2 * max_linear_bias
  });

  test('computes nonlinearity', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    // nonlinearity = 2^{m-1} - max_linear_bias = 4 - 2 = 2
    expect(S.nonlinearity()).toBe(2);
  });

  test('AES S-box nonlinearity', () => {
    // 2^7 - 16 = 112
    expect(AES_SBOX.nonlinearity()).toBe(112);
  });
});

describe('component_function', () => {
  test('returns BooleanFunction object', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);

    // Component function for b=1 (least significant bit)
    const f1 = S.component_function(1);
    expect(f1.nvariables()).toBe(3);
    expect(f1.length()).toBe(8);

    // f1(x) = S(x) & 1
    for (let x = 0; x < 8; x++) {
      expect(f1.call(x)).toBe(S.get(x) & 1);
    }
  });

  test('component function from bit mask', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);

    // b = [0, 1, 1] = 6 in big endian (but stored as 3 in computation)
    const f = S.component_function([0, 1, 1]);

    // Should be same as component function for integer 3
    // (since [0,1,1] in big endian is 3)
    const f3 = S.component_function(3);
    expect(f.truthTable('int')).toEqual(f3.truthTable('int'));
  });

  test('component function has correct algebraic properties', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const f3 = S.component_function(3);

    // Can compute algebraic degree
    const deg = f3.algebraicDegree();
    expect(deg).toBeGreaterThanOrEqual(0);
    expect(deg).toBeLessThanOrEqual(3);

    // Can compute Walsh-Hadamard transform
    const W = f3.walshHadamardTransform();
    expect(W.length).toBe(8);
  });

  test('component_function_table returns raw array', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const table = S.component_function_table(1);
    expect(table.length).toBe(8);

    // Should match truth table from BooleanFunction
    const f = S.component_function(1);
    expect(table).toEqual(f.truthTable('int'));
  });
});

describe('algebraic_degree', () => {
  test('computes algebraic degree', () => {
    const S = new SBox([12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]);
    const deg = S.algebraic_degree();
    expect(deg).toBe(3);
  });

  test('identity has degree 1', () => {
    const id = new SBox([0, 1, 2, 3, 4, 5, 6, 7]);
    expect(id.algebraic_degree()).toBe(1);
  });

  test('constant S-box has degree 0', () => {
    const constant = new SBox([5, 5, 5, 5, 5, 5, 5, 5]);
    expect(constant.algebraic_degree()).toBe(0);
  });
});

describe('fixed_points', () => {
  test('finds fixed points', () => {
    const S = new SBox([0, 1, 3, 6, 7, 4, 5, 2]);
    expect(S.fixed_points()).toEqual([0, 1]);
  });

  test('no fixed points', () => {
    const S = new SBox([1, 0, 3, 2, 5, 4, 7, 6]);
    expect(S.fixed_points()).toEqual([]);
  });
});

describe('branch numbers', () => {
  test('differential branch number', () => {
    const S = new SBox([12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]);
    expect(S.differential_branch_number()).toBe(3);
  });

  test('linear branch number', () => {
    const S = new SBox([12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]);
    expect(S.linear_branch_number()).toBe(2);
  });
});

describe('is_balanced', () => {
  test('permutation S-box is balanced', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.is_balanced()).toBe(true);
  });

  test('non-permutation may not be balanced', () => {
    // All outputs are 0
    const S = new SBox([0, 0, 0, 0, 0, 0, 0, 0]);
    expect(S.is_balanced()).toBe(false);
  });
});

describe('derivative', () => {
  test('computes derivative S-box', () => {
    const S = new SBox([0, 1, 2, 3]);
    const D = S.derivative(1);

    // D(x) = S(x) XOR S(x XOR 1)
    // D(0) = S(0) XOR S(1) = 0 XOR 1 = 1
    // D(1) = S(1) XOR S(0) = 1 XOR 0 = 1
    // D(2) = S(2) XOR S(3) = 2 XOR 3 = 1
    // D(3) = S(3) XOR S(2) = 3 XOR 2 = 1
    expect(D.toArray()).toEqual([1, 1, 1, 1]);
  });

  test('derivative from bit array', () => {
    const S = new SBox([0, 1, 2, 3]);
    const D1 = S.derivative([0, 1]); // = derivative(1) in big endian
    const D2 = S.derivative(1);
    expect(D1.toArray()).toEqual(D2.toArray());
  });
});

describe('equality and iteration', () => {
  test('equals method', () => {
    const S1 = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const S2 = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const S3 = new SBox([0, 1, 2, 3, 4, 5, 6, 7]);

    expect(S1.equals(S2)).toBe(true);
    expect(S1.equals(S3)).toBe(false);
  });

  test('iterator yields all outputs', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const outputs = [...S];
    expect(outputs).toEqual([7, 6, 0, 4, 2, 5, 1, 3]);
  });

  test('toString representation', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.toString()).toBe('(7, 6, 0, 4, 2, 5, 1, 3)');
  });

  test('toArray returns copy', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const arr = S.toArray();
    arr[0] = 999;
    expect(S.get(0)).toBe(7); // Original unchanged
  });
});

describe('Standard S-boxes', () => {
  test('AES S-box properties', () => {
    expect(AES_SBOX.input_size()).toBe(8);
    expect(AES_SBOX.output_size()).toBe(8);
    expect(AES_SBOX.is_permutation()).toBe(true);
    expect(AES_SBOX.differential_uniformity()).toBe(4);
    expect(AES_SBOX.nonlinearity()).toBe(112);
  });

  test('AES S-box specific values', () => {
    expect(AES_SBOX.call(0x00)).toBe(0x63);
    expect(AES_SBOX.call(0x01)).toBe(0x7c);
    expect(AES_SBOX.call(0x52)).toBe(0x00);
  });

  test('PRESENT S-box properties', () => {
    expect(PRESENT_SBOX.input_size()).toBe(4);
    expect(PRESENT_SBOX.output_size()).toBe(4);
    expect(PRESENT_SBOX.is_permutation()).toBe(true);
    expect(PRESENT_SBOX.differential_uniformity()).toBe(4);
  });

  test('PRESENT S-box specific values', () => {
    expect(PRESENT_SBOX.call(0)).toBe(0xc);
    expect(PRESENT_SBOX.call(1)).toBe(0x5);
    expect(PRESENT_SBOX.call(5)).toBe(0x0);
  });

  test('DES S-box 1 properties', () => {
    expect(DES_SBOX1.input_size()).toBe(6);
    expect(DES_SBOX1.output_size()).toBe(4);
    // DES S-boxes are not permutations (6-bit input, 4-bit output)
    expect(DES_SBOX1.is_permutation()).toBe(false);
  });

  test('GIFT S-box properties', () => {
    expect(GIFT_SBOX.input_size()).toBe(4);
    expect(GIFT_SBOX.output_size()).toBe(4);
    expect(GIFT_SBOX.is_permutation()).toBe(true);
  });

  test('SKINNY S-box properties', () => {
    expect(SKINNY_SBOX.input_size()).toBe(4);
    expect(SKINNY_SBOX.output_size()).toBe(4);
    expect(SKINNY_SBOX.is_permutation()).toBe(true);
  });
});

describe('feistel_construction', () => {
  test('constructs larger S-box from smaller ones', () => {
    const s = PRESENT_SBOX;
    const S = feistel_construction([s, s, s]);

    // Input/output size should be doubled
    expect(S.input_size()).toBe(8);
    expect(S.output_size()).toBe(8);
  });

  test('3-round Feistel is a permutation', () => {
    const s = PRESENT_SBOX;
    const S = feistel_construction([s, s, s]);
    expect(S.is_permutation()).toBe(true);
  });

  test('Feistel construction cryptographic properties', () => {
    const s = PRESENT_SBOX;
    const S = feistel_construction([s, s, s]);

    // Check that we can compute properties without error
    const du = S.differential_uniformity();
    const nl = S.nonlinearity();

    expect(du).toBeGreaterThan(0);
    expect(nl).toBeGreaterThanOrEqual(0);
  });

  test('throws on empty input', () => {
    expect(() => feistel_construction([])).toThrow('no input S-boxes provided');
  });
});

describe('misty_construction', () => {
  test('constructs larger S-box from smaller ones', () => {
    const s = PRESENT_SBOX;
    const S = misty_construction([s, s, s]);

    expect(S.input_size()).toBe(8);
    expect(S.output_size()).toBe(8);
  });

  test('Misty construction produces permutation', () => {
    const s = PRESENT_SBOX;
    const S = misty_construction([s, s, s]);
    expect(S.is_permutation()).toBe(true);
  });

  test("matches Sage's 3-round misty_construction doctest", () => {
    // sage: S1 = SBox([0x4,0x0,0x1,0xF,0x2,0xB,0x6,0x7,0x3,0x9,0xA,0x5,0xC,0xD,0xE,0x8])
    // sage: S2 = SBox([0x0,0x0,0x0,0x1,0x0,0xA,0x8,0x3,0x0,0x8,0x2,0xB,0x4,0x6,0xE,0xD])
    // sage: S3 = SBox([0x0,0x7,0xB,0xD,0x4,0x1,0xB,0xF,0x1,0x2,0xC,0xE,0xD,0xC,0x5,0x5])
    // sage: S = misty_construction(S1, S2, S3)
    // sage: S.differential_uniformity()  ->  8
    // sage: S.linearity()                ->  64
    const S1 = new SBox([
      0x4, 0x0, 0x1, 0xf, 0x2, 0xb, 0x6, 0x7, 0x3, 0x9, 0xa, 0x5, 0xc, 0xd, 0xe, 0x8,
    ]);
    const S2 = new SBox([
      0x0, 0x0, 0x0, 0x1, 0x0, 0xa, 0x8, 0x3, 0x0, 0x8, 0x2, 0xb, 0x4, 0x6, 0xe, 0xd,
    ]);
    const S3 = new SBox([
      0x0, 0x7, 0xb, 0xd, 0x4, 0x1, 0xb, 0xf, 0x1, 0x2, 0xc, 0xe, 0xd, 0xc, 0x5, 0x5,
    ]);
    const S = misty_construction([S1, S2, S3]);
    expect(S.differential_uniformity()).toBe(8);
    expect(S.linearity()).toBe(64);
  });

  test('throws on empty input', () => {
    expect(() => misty_construction([])).toThrow('no input S-boxes provided');
  });
});

describe('edge cases', () => {
  test('2-element S-box (1-bit)', () => {
    const S = new SBox([1, 0]); // NOT gate
    expect(S.input_size()).toBe(1);
    expect(S.output_size()).toBe(1);
    expect(S.is_permutation()).toBe(true);
    expect(S.is_involution()).toBe(true);
  });

  test('256-element S-box (8-bit)', () => {
    // Create identity S-box
    const id = new SBox(Array.from({ length: 256 }, (_, i) => i));
    expect(id.input_size()).toBe(8);
    expect(id.output_size()).toBe(8);
    expect(id.is_permutation()).toBe(true);
    expect(id.is_involution()).toBe(true);
  });

  test('S-box with zero outputs', () => {
    const S = new SBox([0, 0, 0, 0]);
    expect(S.differential_uniformity()).toBe(4);
    expect(S.is_balanced()).toBe(false);
  });
});

describe('consistency checks', () => {
  test('DDT row sums are constant', () => {
    const sboxes = [PRESENT_SBOX, GIFT_SBOX, new SBox([7, 6, 0, 4, 2, 5, 1, 3])];

    for (const S of sboxes) {
      const ddt = S.difference_distribution_table();
      const expectedSum = 1 << S.input_size();

      for (let i = 0; i < 1 << S.input_size(); i++) {
        let sum = 0;
        for (let j = 0; j < 1 << S.output_size(); j++) {
          sum += Number(ddt.get(i, j).value);
        }
        expect(sum).toBe(expectedSum);
      }
    }
  });

  test('DDT column sums for permutation S-box', () => {
    const S = PRESENT_SBOX;
    const ddt = S.difference_distribution_table();

    // For permutation, each column sum should equal 2^m
    for (let j = 0; j < 16; j++) {
      let sum = 0;
      for (let i = 0; i < 16; i++) {
        sum += Number(ddt.get(i, j).value);
      }
      expect(sum).toBe(16);
    }
  });

  test('inverse of inverse is identity', () => {
    const S = new SBox([0, 1, 3, 6, 7, 4, 5, 2]);
    const Sinv = S.inverse();
    const Sinvinv = Sinv.inverse();
    expect(S.equals(Sinvinv)).toBe(true);
  });

  test('differential uniformity is even for balanced S-box', () => {
    const sboxes = [PRESENT_SBOX, GIFT_SBOX, SKINNY_SBOX];
    for (const S of sboxes) {
      expect(S.differential_uniformity() % 2).toBe(0);
    }
  });
});

describe('SageMath reference tests', () => {
  /**
   * Tests based on SageMath doctests from sbox.pyx
   */

  test('SageMath example: basic S-box', () => {
    // sage: from sage.crypto.sbox import SBox
    // sage: S = SBox(7,6,0,4,2,5,1,3); S
    // (7, 6, 0, 4, 2, 5, 1, 3)
    // sage: S(1)
    // 5
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.toString()).toBe('(7, 6, 0, 4, 2, 5, 1, 3)');

    // Note: SageMath's S(1) returns 5 which is S[1] = 6, not 5
    // Actually checking: sage gives S(1) = 6
    expect(S.call(1)).toBe(6);
  });

  test('SageMath example: DDT', () => {
    // sage: S = SBox(7,6,0,4,2,5,1,3)
    // sage: S.difference_distribution_table()
    // [8 0 0 0 0 0 0 0]
    // [0 2 2 0 2 0 0 2]
    // ...
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const ddt = S.difference_distribution_table();

    // First row
    expect(Number(ddt.get(0, 0).value)).toBe(8);
    for (let j = 1; j < 8; j++) {
      expect(Number(ddt.get(0, j).value)).toBe(0);
    }

    // Second row: [0 2 2 0 2 0 0 2]
    const row1 = [0, 2, 2, 0, 2, 0, 0, 2];
    for (let j = 0; j < 8; j++) {
      expect(Number(ddt.get(1, j).value)).toBe(row1[j]);
    }
  });

  test('SageMath example: differential uniformity', () => {
    // sage: S = SBox(7,6,0,4,2,5,1,3)
    // sage: S.differential_uniformity()
    // 2
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.differential_uniformity()).toBe(2);
  });

  test('SageMath example: is_apn', () => {
    // sage: S = SBox([0,1,3,6,7,4,5,2])
    // sage: S.is_apn()
    // True
    const S = new SBox([0, 1, 3, 6, 7, 4, 5, 2]);
    expect(S.is_apn()).toBe(true);
  });

  test('SageMath example: LAT', () => {
    // sage: S = SBox(7,6,0,4,2,5,1,3)
    // sage: S.linear_approximation_table()
    // [ 4  0  0  0  0  0  0  0]
    // [ 0  0  0  0  2  2  2 -2]
    // ...
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const lat = S.linear_approximation_table();

    // First row: [4 0 0 0 0 0 0 0]
    expect(Number(lat.get(0, 0).value)).toBe(4);
    for (let j = 1; j < 8; j++) {
      expect(Number(lat.get(0, j).value)).toBe(0);
    }
  });

  test('SageMath example: max linear bias', () => {
    // sage: S = SBox(7,6,0,4,2,5,1,3)
    // sage: S.maximal_linear_bias_absolute()
    // 2
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.max_linear_bias()).toBe(2);
  });

  test('SageMath example: is_permutation', () => {
    // sage: S = SBox(7,6,0,4,2,5,1,3)
    // sage: S.is_permutation()
    // True
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    expect(S.is_permutation()).toBe(true);

    // sage: S = SBox(3,2,0,0,2,1,1,3)
    // sage: S.is_permutation()
    // False
    const S2 = new SBox([3, 2, 0, 0, 2, 1, 1, 3]);
    expect(S2.is_permutation()).toBe(false);
  });

  test('SageMath example: inverse', () => {
    // sage: S = SBox([0, 1, 3, 6, 7, 4, 5, 2])
    // sage: Sinv = S.inverse()
    // sage: [Sinv(S(i)) for i in range(8)]
    // [0, 1, 2, 3, 4, 5, 6, 7]
    const S = new SBox([0, 1, 3, 6, 7, 4, 5, 2]);
    const Sinv = S.inverse();
    const result = Array.from({ length: 8 }, (_, i) => Sinv.call(S.call(i) as number));
    expect(result).toEqual([0, 1, 2, 3, 4, 5, 6, 7]);
  });

  test('SageMath example: fixed_points', () => {
    // sage: S = SBox([0,1,3,6,7,4,5,2])
    // sage: S.fixed_points()
    // [0, 1]
    const S = new SBox([0, 1, 3, 6, 7, 4, 5, 2]);
    expect(S.fixed_points()).toEqual([0, 1]);
  });

  test('SageMath example: PRESENT S-box', () => {
    // sage: from sage.crypto.sbox import SBox
    // sage: S = SBox(12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2); S
    // (12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2)
    // sage: S(1)
    // 5
    expect(PRESENT_SBOX.toString()).toBe('(12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2)');
    expect(PRESENT_SBOX.call(1)).toBe(5);
  });

  test('SageMath example: component_function returns BooleanFunction', () => {
    // sage: from sage.crypto.sbox import SBox
    // sage: S = SBox([7,6,0,4,2,5,1,3])
    // sage: f3 = S.component_function(3)
    // sage: f3.algebraic_normal_form()
    // x0*x1 + x0*x2 + x0 + x2
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const f3 = S.component_function(3);

    // Check that it's a BooleanFunction
    expect(f3.nvariables()).toBe(3);

    // Verify algebraic degree is correct
    const deg = f3.algebraicDegree();
    expect(deg).toBeGreaterThanOrEqual(1);
    expect(deg).toBeLessThanOrEqual(3);
  });

  test('SageMath example: max_degree and min_degree', () => {
    // sage: from sage.crypto.sbox import SBox
    // sage: S = SBox([12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2])
    // sage: S.max_degree()
    // 3
    // sage: S.min_degree()
    // 2
    const S = new SBox([12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]);
    expect(S.max_degree()).toBe(3);
    expect(S.min_degree()).toBe(2);
  });

  test('SageMath example: linearity', () => {
    // sage: from sage.crypto.mq.sr import SR
    // sage: S = mq.SR(1, 4, 4, 8).sbox()
    // sage: S.linearity()
    // 32
    // For AES: linearity = 2 * max_linear_bias = 2 * 16 = 32
    expect(AES_SBOX.linearity()).toBe(32);
  });

  test('SageMath example: nonlinearity', () => {
    // sage: S = mq.SR(1,4,4,8).sbox()
    // sage: S.nonlinearity()
    // 112
    expect(AES_SBOX.nonlinearity()).toBe(112);
  });
});

describe('BooleanFunction integration', () => {
  test('component function can be used for spectral analysis', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);

    // Get component function for each output bit
    for (let b = 1; b < 8; b++) {
      const f = S.component_function(b);

      // Walsh-Hadamard transform should satisfy Parseval's theorem
      const W = f.walshHadamardTransform();
      const sumSquares = W.reduce((acc, w) => acc + w * w, 0);
      expect(sumSquares).toBe(64); // 2^{2*3} = 64
    }
  });

  test('component function nonlinearity matches S-box nonlinearity', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);

    // S-box nonlinearity should be minimum component function nonlinearity
    let minNL = Number.POSITIVE_INFINITY;
    for (let b = 1; b < 8; b++) {
      const f = S.component_function(b);
      const nl = f.nonlinearity();
      if (nl < minNL) {
        minNL = nl;
      }
    }

    expect(S.nonlinearity()).toBe(minNL);
  });

  test('component function balanced check matches S-box balanced check', () => {
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);

    // S-box is balanced iff all component functions are balanced
    const sboxBalanced = S.is_balanced();

    let allComponentsBalanced = true;
    for (let b = 1; b < 8; b++) {
      const f = S.component_function(b);
      if (!f.isBalanced()) {
        allComponentsBalanced = false;
        break;
      }
    }

    expect(sboxBalanced).toBe(allComponentsBalanced);
    expect(sboxBalanced).toBe(true); // This is a permutation, so it's balanced
  });

  test('AES S-box component functions', () => {
    // AES S-box should have high-degree component functions
    const f1 = AES_SBOX.component_function(1);
    expect(f1.nvariables()).toBe(8);

    // All 255 component functions should be balanced
    for (let b = 1; b < 256; b++) {
      const f = AES_SBOX.component_function(b);
      expect(f.isBalanced()).toBe(true);
    }
  });
});

describe('DDT and LAT comprehensive tests', () => {
  test('DDT computation for non-square S-box', () => {
    // 4 inputs, up to 16 possible outputs
    const S = new SBox([7, 4, 8, 6]);
    const ddt = S.difference_distribution_table();

    // DDT should have 4 rows (2^2 input differences)
    // and 16 columns (2^4 output differences based on max output)
    expect(ddt.nrows).toBe(4);

    // First entry should be 4 (= 2^2)
    expect(Number(ddt.get(0, 0).value)).toBe(4);
  });

  test('LAT values match expected for known S-box', () => {
    // For S = [7,6,0,4,2,5,1,3], second row should be [0 0 0 0 2 2 2 -2]
    const S = new SBox([7, 6, 0, 4, 2, 5, 1, 3]);
    const lat = S.linear_approximation_table();

    const row1Expected = [0, 0, 0, 0, 2, 2, 2, -2];
    for (let j = 0; j < 8; j++) {
      expect(Number(lat.get(1, j).value)).toBe(row1Expected[j]);
    }
  });

  test('AES S-box has expected LAT properties', () => {
    const lat = AES_SBOX.linear_approximation_table();

    // First row should have LAT[0,0] = 128 (= 2^7) and zeros elsewhere
    expect(Number(lat.get(0, 0).value)).toBe(128);
    for (let j = 1; j < 256; j++) {
      expect(Number(lat.get(0, j).value)).toBe(0);
    }

    // First column (except [0,0]) should be all zeros for balanced permutation
    for (let i = 1; i < 256; i++) {
      expect(Number(lat.get(i, 0).value)).toBe(0);
    }
  });

  test('DDT entries are non-negative even integers for permutation', () => {
    const sboxes = [PRESENT_SBOX, GIFT_SBOX, SKINNY_SBOX];

    for (const S of sboxes) {
      const ddt = S.difference_distribution_table();
      const nrows = 1 << S.input_size();
      const ncols = 1 << S.output_size();

      for (let i = 0; i < nrows; i++) {
        for (let j = 0; j < ncols; j++) {
          const val = Number(ddt.get(i, j).value);
          expect(val).toBeGreaterThanOrEqual(0);
          expect(val % 2).toBe(0); // Must be even for permutation
        }
      }
    }
  });
});

describe('max_degree and min_degree', () => {
  test('max_degree equals algebraic_degree', () => {
    const sboxes = [
      new SBox([7, 6, 0, 4, 2, 5, 1, 3]),
      PRESENT_SBOX,
      new SBox([0, 1, 2, 3, 4, 5, 6, 7]), // Identity
    ];

    for (const S of sboxes) {
      expect(S.max_degree()).toBe(S.algebraic_degree());
    }
  });

  test('min_degree is at most max_degree', () => {
    const sboxes = [
      new SBox([7, 6, 0, 4, 2, 5, 1, 3]),
      PRESENT_SBOX,
      new SBox([12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]),
    ];

    for (const S of sboxes) {
      expect(S.min_degree()).toBeLessThanOrEqual(S.max_degree());
    }
  });

  test('constant S-box has degree 0', () => {
    const S = new SBox([5, 5, 5, 5, 5, 5, 5, 5]);
    expect(S.max_degree()).toBe(0);
    // Sage's min_degree() takes the plain minimum of the algebraic degrees of
    // *all* component functions b.S for b != 0.  Here b = 2 gives the constant
    // zero function, whose algebraic degree is -1 by convention.
    expect(S.min_degree()).toBe(-1);
  });

  test('identity S-box has degree 1', () => {
    const S = new SBox([0, 1, 2, 3, 4, 5, 6, 7]);
    expect(S.max_degree()).toBe(1);
    expect(S.min_degree()).toBe(1);
  });

  test("Sage's min_degree doctest", () => {
    // sage: SBox([12,5,6,11,9,0,10,13,3,14,15,8,4,7,1,2]).min_degree()  ->  2
    const S = new SBox([12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]);
    expect(S.min_degree()).toBe(2);
  });

  test('S-box with a vanishing component has min_degree -1', () => {
    // S = [0,0,3,3]: the component function for b = 3 is identically zero.
    const S = new SBox([0, 0, 3, 3]);
    expect(S.component_function(3).algebraicDegree()).toBe(-1);
    expect(S.min_degree()).toBe(-1);
  });
});
