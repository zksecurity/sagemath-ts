/**
 * Tests for Reed-Muller codes
 */

import { describe, expect, it } from 'vitest';
import { ValueError } from '../errors.js';
import { GF2, GF2Element } from '../rings/finite_rings/gf2.js';
import {
  BinaryReedMullerCode,
  FirstOrderReedMullerCode,
  PuncturedReedMullerCode,
  ReedMullerCode,
  RepetitionCode,
} from './reed_muller_code.js';

describe('ReedMullerCode', () => {
  describe('constructor', () => {
    it('should create a valid RM(1, 3) code', () => {
      const rm = new ReedMullerCode(1n, 3n);

      expect(rm.order()).toBe(1);
      expect(rm.num_variables()).toBe(3);
      expect(rm.number_of_variables()).toBe(3); // alias
    });

    it('should create a valid RM(2, 4) code', () => {
      const rm = new ReedMullerCode(2n, 4n);

      expect(rm.order()).toBe(2);
      expect(rm.num_variables()).toBe(4);
    });

    it('should throw for r > m', () => {
      expect(() => new ReedMullerCode(3n, 2n)).toThrow(ValueError);
    });

    it('should throw for negative r', () => {
      expect(() => new ReedMullerCode(-1n, 3n)).toThrow(ValueError);
    });

    it('should throw for negative m', () => {
      expect(() => new ReedMullerCode(0n, -1n)).toThrow(ValueError);
    });

    it('should throw for non-integer r', () => {
      expect(() => new ReedMullerCode(1.5 as any, 3n)).toThrow(TypeError);
    });

    it('should throw for non-integer m', () => {
      expect(() => new ReedMullerCode(1n, 3.5 as any)).toThrow(TypeError);
    });

    it('should accept r = m (full space)', () => {
      const rm = new ReedMullerCode(3n, 3n);
      expect(rm.order()).toBe(3);
      expect(rm.is_full_space()).toBe(true);
    });

    it('should accept r = 0 (repetition code)', () => {
      const rm = new ReedMullerCode(0n, 3n);
      expect(rm.order()).toBe(0);
      expect(rm.is_repetition_code()).toBe(true);
    });
  });

  describe('code parameters', () => {
    it('should compute RM(1, 3) parameters correctly: [8, 4, 4]', () => {
      const rm = new ReedMullerCode(1n, 3n);

      expect(rm.length()).toBe(8); // 2^3 = 8
      expect(rm.dimension()).toBe(4); // 1 + 3 = 4
      expect(rm.minimum_distance()).toBe(4); // 2^(3-1) = 4
      expect(rm.parameters()).toEqual([8, 4, 4]);
    });

    it('should compute RM(2, 4) parameters correctly: [16, 11, 4]', () => {
      const rm = new ReedMullerCode(2n, 4n);

      expect(rm.length()).toBe(16); // 2^4 = 16
      expect(rm.dimension()).toBe(11); // 1 + 4 + 6 = 11
      expect(rm.minimum_distance()).toBe(4); // 2^(4-2) = 4
      expect(rm.parameters()).toEqual([16, 11, 4]);
    });

    it('should compute RM(0, 3) (repetition) parameters: [8, 1, 8]', () => {
      const rm = new ReedMullerCode(0n, 3n);

      expect(rm.length()).toBe(8);
      expect(rm.dimension()).toBe(1);
      expect(rm.minimum_distance()).toBe(8);
      expect(rm.parameters()).toEqual([8, 1, 8]);
    });

    it('should compute RM(2, 3) (parity check) parameters: [8, 7, 2]', () => {
      const rm = new ReedMullerCode(2n, 3n);

      expect(rm.length()).toBe(8);
      expect(rm.dimension()).toBe(7); // 1 + 3 + 3 = 7
      expect(rm.minimum_distance()).toBe(2); // 2^(3-2) = 2
      expect(rm.parameters()).toEqual([8, 7, 2]);
      expect(rm.is_single_parity_check()).toBe(true);
    });

    it('should compute RM(3, 3) (full space) parameters: [8, 8, 1]', () => {
      const rm = new ReedMullerCode(3n, 3n);

      expect(rm.length()).toBe(8);
      expect(rm.dimension()).toBe(8); // 1 + 3 + 3 + 1 = 8
      expect(rm.minimum_distance()).toBe(1); // 2^(3-3) = 1
      expect(rm.parameters()).toEqual([8, 8, 1]);
      expect(rm.is_full_space()).toBe(true);
    });

    it('should compute RM(1, 4) parameters: [16, 5, 8]', () => {
      const rm = new ReedMullerCode(1n, 4n);

      expect(rm.length()).toBe(16);
      expect(rm.dimension()).toBe(5); // 1 + 4 = 5
      expect(rm.minimum_distance()).toBe(8); // 2^(4-1) = 8
      expect(rm.parameters()).toEqual([16, 5, 8]);
    });

    it('should compute decoding radius correctly', () => {
      const rm13 = new ReedMullerCode(1n, 3n);
      expect(rm13.decoding_radius()).toBe(1); // floor((4-1)/2) = 1

      const rm24 = new ReedMullerCode(2n, 4n);
      expect(rm24.decoding_radius()).toBe(1); // floor((4-1)/2) = 1

      const rm14 = new ReedMullerCode(1n, 4n);
      expect(rm14.decoding_radius()).toBe(3); // floor((8-1)/2) = 3
    });

    it('should compute code rate correctly', () => {
      const rm13 = new ReedMullerCode(1n, 3n);
      expect(rm13.rate()).toBe(0.5); // 4/8

      const rm24 = new ReedMullerCode(2n, 4n);
      expect(rm24.rate()).toBe(11 / 16);
    });
  });

  describe('generator matrix', () => {
    it('should have correct dimensions for RM(1, 3)', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const G = rm.generator_matrix();

      expect(G.nrows).toBe(4); // dimension
      expect(G.ncols).toBe(8); // length
    });

    it('should have correct dimensions for RM(2, 4)', () => {
      const rm = new ReedMullerCode(2n, 4n);
      const G = rm.generator_matrix();

      expect(G.nrows).toBe(11); // dimension
      expect(G.ncols).toBe(16); // length
    });

    it('should have first row all ones for RM(1, m)', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const G = rm.generator_matrix();

      // First monomial is constant (1), so first row should be all 1s
      for (let j = 0; j < 8; j++) {
        expect(G.get(0, j).value).toBe(1);
      }
    });

    it('should be cached', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const G1 = rm.generator_matrix();
      const G2 = rm.generator_matrix();

      expect(G1).toBe(G2); // Same object reference
    });

    it('rows should be linearly independent', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const G = rm.generator_matrix();

      // Check that rows are not duplicates (basic linear independence check)
      const rows = G.rows();
      const seenRows = new Set<string>();
      for (const row of rows) {
        const key = row.map((x) => x.value).join(',');
        expect(seenRows.has(key)).toBe(false);
        seenRows.add(key);
      }
    });
  });

  describe('encode', () => {
    it('should encode message of correct length', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const message = [1, 0, 1, 0];
      const codeword = rm.encode(message);

      expect(codeword.length).toBe(8);
    });

    it('should encode zero message to zero codeword', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const message = [0, 0, 0, 0];
      const codeword = rm.encode(message);

      expect(codeword.every((c) => c.value === 0)).toBe(true);
    });

    it('should encode all-ones message correctly for RM(0, m)', () => {
      const rm = new ReedMullerCode(0n, 3n);
      const message = [1]; // Single coefficient for constant term
      const codeword = rm.encode(message);

      // Repetition code: all 1s
      expect(codeword.every((c) => c.value === 1)).toBe(true);
    });

    it('should throw for wrong message length', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const message = [1, 0, 1]; // 3 elements instead of 4

      expect(() => rm.encode(message)).toThrow(ValueError);
    });

    it('should accept GF2Element inputs', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const message = [GF2.__call__(1), GF2.__call__(0), GF2.__call__(1), GF2.__call__(0)];
      const codeword = rm.encode(message);

      expect(codeword.length).toBe(8);
    });

    it('should produce codewords in the code', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const message = [1, 1, 0, 1];
      const codeword = rm.encode(message);

      // Verify by checking minimum weight
      // A non-zero codeword should have weight at least d
      const weight = ReedMullerCode.hamming_weight(codeword);
      if (weight > 0) {
        expect(weight).toBeGreaterThanOrEqual(rm.minimum_distance());
      }
    });
  });

  describe('decode', () => {
    it('should decode a valid codeword correctly', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const message = [1, 0, 1, 0];
      const codeword = rm.encode(message);
      const decoded = rm.decode(codeword);

      expect(decoded.length).toBe(4);
      for (let i = 0; i < 4; i++) {
        expect(decoded[i]!.value).toBe(message[i]);
      }
    });

    it('should correct one error in RM(1, 3)', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const message = [1, 0, 1, 1];
      const codeword = rm.encode(message);

      // Flip one bit
      const corrupted = codeword.map((x) => GF2.__call__(x.value));
      corrupted[3] = GF2.__call__(1 - corrupted[3]!.value);

      const decoded = rm.decode(corrupted);

      for (let i = 0; i < 4; i++) {
        expect(decoded[i]!.value).toBe(message[i]);
      }
    });

    it('should correct multiple errors up to decoding radius', () => {
      const rm = new ReedMullerCode(1n, 4n);
      const message = [1, 0, 1, 0, 1];
      const codeword = rm.encode(message);

      // RM(1, 4) has d = 8, so can correct up to 3 errors
      const corrupted = codeword.map((x) => GF2.__call__(x.value));
      corrupted[0] = GF2.__call__(1 - corrupted[0]!.value);
      corrupted[5] = GF2.__call__(1 - corrupted[5]!.value);
      corrupted[10] = GF2.__call__(1 - corrupted[10]!.value);

      const decoded = rm.decode(corrupted);

      for (let i = 0; i < 5; i++) {
        expect(decoded[i]!.value).toBe(message[i]);
      }
    });

    it('should throw for wrong received length', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const received = [1, 0, 1, 0, 1, 0, 1]; // 7 elements instead of 8

      expect(() => rm.decode(received)).toThrow(ValueError);
    });

    it('should handle zero vector correctly', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const zeroWord = [0, 0, 0, 0, 0, 0, 0, 0];
      const decoded = rm.decode(zeroWord);

      expect(decoded.every((c) => c.value === 0)).toBe(true);
    });
  });

  describe('encode/decode round-trip', () => {
    it('should round-trip all messages for RM(1, 2)', () => {
      const rm = new ReedMullerCode(1n, 2n);
      const k = rm.dimension(); // 3

      // Test all 2^k = 8 possible messages
      for (let i = 0; i < 1 << k; i++) {
        const message: number[] = [];
        for (let j = 0; j < k; j++) {
          message.push((i >> j) & 1);
        }

        const codeword = rm.encode(message);
        const decoded = rm.decode(codeword);

        for (let j = 0; j < k; j++) {
          expect(decoded[j]!.value).toBe(message[j]);
        }
      }
    });

    it('should round-trip with random messages for RM(2, 4)', () => {
      const rm = new ReedMullerCode(2n, 4n);
      const k = rm.dimension(); // 11

      // Test a few random messages
      const testMessages = [
        [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
        [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
        [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
      ];

      for (const message of testMessages) {
        const codeword = rm.encode(message);
        const decoded = rm.decode(codeword);

        for (let j = 0; j < k; j++) {
          expect(decoded[j]!.value).toBe(message[j]);
        }
      }
    });
  });

  describe('Plotkin decomposition', () => {
    it('should decompose a codeword correctly', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const message = [1, 0, 1, 0];
      const codeword = rm.encode(message);

      const [u, v] = rm.plotkin_decomposition(codeword);

      expect(u.length).toBe(4);
      expect(v.length).toBe(4);

      // Verify: codeword should equal (u, u+v)
      for (let i = 0; i < 4; i++) {
        expect(codeword[i]!.value).toBe(u[i]!.value);
        expect(codeword[i + 4]!.value).toBe(u[i]!.add(v[i]!).value);
      }
    });

    it('should work with compose and decompose being inverses', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const message = [1, 1, 0, 1];
      const codeword = rm.encode(message);

      const [u, v] = rm.plotkin_decomposition(codeword);
      const reconstructed = ReedMullerCode.plotkin_compose(u, v);

      for (let i = 0; i < 8; i++) {
        expect(reconstructed[i]!.value).toBe(codeword[i]!.value);
      }
    });

    it('should throw for m = 0', () => {
      const rm = new ReedMullerCode(0n, 0n);
      const codeword = [GF2.__call__(1)];

      expect(() => rm.plotkin_decomposition(codeword)).toThrow(ValueError);
    });

    it('should throw for wrong codeword length', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const wrongWord = [1, 0, 1, 0]; // 4 elements instead of 8

      expect(() => rm.plotkin_decomposition(wrongWord)).toThrow(ValueError);
    });
  });

  describe('special codes', () => {
    it('should correctly identify repetition code', () => {
      const rm = new ReedMullerCode(0n, 4n);

      expect(rm.is_repetition_code()).toBe(true);
      expect(rm.is_first_order()).toBe(false);
      expect(rm.is_single_parity_check()).toBe(false);
      expect(rm.is_full_space()).toBe(false);
    });

    it('should correctly identify first-order code', () => {
      const rm = new ReedMullerCode(1n, 4n);

      expect(rm.is_repetition_code()).toBe(false);
      expect(rm.is_first_order()).toBe(true);
      expect(rm.is_single_parity_check()).toBe(false);
      expect(rm.is_full_space()).toBe(false);
    });

    it('should correctly identify single parity check code', () => {
      const rm = new ReedMullerCode(3n, 4n);

      expect(rm.is_repetition_code()).toBe(false);
      expect(rm.is_first_order()).toBe(false);
      expect(rm.is_single_parity_check()).toBe(true);
      expect(rm.is_full_space()).toBe(false);
    });

    it('should correctly identify full space', () => {
      const rm = new ReedMullerCode(4n, 4n);

      expect(rm.is_repetition_code()).toBe(false);
      expect(rm.is_first_order()).toBe(false);
      expect(rm.is_single_parity_check()).toBe(false);
      expect(rm.is_full_space()).toBe(true);
    });
  });

  describe('dual code', () => {
    it('should compute dual correctly', () => {
      const rm = new ReedMullerCode(1n, 4n);
      const dual = rm.dual();

      // Dual of RM(1, 4) is RM(4-1-1, 4) = RM(2, 4)
      expect(dual.order()).toBe(2);
      expect(dual.num_variables()).toBe(4);
    });

    it('should satisfy duality properties', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const dual = rm.dual();

      // RM(1, 3) has k = 4, RM(1, 3)^perp = RM(1, 3) (self-dual!)
      // Actually dual of RM(1, 3) is RM(3-1-1, 3) = RM(1, 3)
      expect(dual.order()).toBe(1);
      expect(dual.num_variables()).toBe(3);
      expect(rm.dimension() + dual.dimension()).toBe(rm.length()); // k + k^perp = n
    });

    it('should throw for RM(m, m)', () => {
      const rm = new ReedMullerCode(3n, 3n);
      expect(() => rm.dual()).toThrow(ValueError);
    });
  });

  describe('contains', () => {
    it('should return true for valid codewords', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const message = [1, 0, 1, 0];
      const codeword = rm.encode(message);

      expect(rm.contains(codeword)).toBe(true);
    });

    it('should return false for non-codewords', () => {
      const rm = new ReedMullerCode(1n, 3n);

      // A random word is unlikely to be a codeword
      // unless it happens to have weight >= d
      const word = [1, 0, 0, 0, 0, 0, 0, 0];
      // This is not a codeword (weight 1 < minimum distance 4)
      expect(rm.contains(word)).toBe(false);
    });

    it('should return false for wrong length', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const word = [1, 0, 1, 0, 1];

      expect(rm.contains(word)).toBe(false);
    });
  });

  describe('Hamming weight and distance', () => {
    it('should compute Hamming weight correctly', () => {
      expect(ReedMullerCode.hamming_weight([1, 0, 1, 1, 0, 0, 1, 0])).toBe(4);
      expect(ReedMullerCode.hamming_weight([0, 0, 0, 0])).toBe(0);
      expect(ReedMullerCode.hamming_weight([1, 1, 1, 1])).toBe(4);
    });

    it('should compute Hamming distance correctly', () => {
      const a = [1, 0, 1, 0, 1, 0, 1, 0];
      const b = [1, 1, 1, 1, 0, 0, 0, 0];
      expect(ReedMullerCode.hamming_distance(a, b)).toBe(4);
    });

    it('should throw for different lengths in Hamming distance', () => {
      const a = [1, 0, 1];
      const b = [1, 0, 1, 0];
      expect(() => ReedMullerCode.hamming_distance(a, b)).toThrow(ValueError);
    });
  });

  describe('minimum distance verification', () => {
    it('should verify minimum distance for RM(1, 3) by enumeration', () => {
      const rm = new ReedMullerCode(1n, 3n);
      const k = rm.dimension();
      const expectedD = rm.minimum_distance();

      let minWeight = Number.POSITIVE_INFINITY;

      // Enumerate all non-zero codewords
      for (let i = 1; i < 1 << k; i++) {
        const message: number[] = [];
        for (let j = 0; j < k; j++) {
          message.push((i >> j) & 1);
        }

        const codeword = rm.encode(message);
        const weight = ReedMullerCode.hamming_weight(codeword);
        minWeight = Math.min(minWeight, weight);
      }

      expect(minWeight).toBe(expectedD);
    });

    it('should verify minimum distance for RM(0, 3) (repetition)', () => {
      const rm = new ReedMullerCode(0n, 3n);

      // Only two codewords: all zeros and all ones
      const zero = rm.encode([0]);
      const one = rm.encode([1]);

      expect(ReedMullerCode.hamming_weight(zero)).toBe(0);
      expect(ReedMullerCode.hamming_weight(one)).toBe(8);
      expect(rm.minimum_distance()).toBe(8);
    });
  });
});

describe('PuncturedReedMullerCode', () => {
  it('should create a punctured code', () => {
    const rm = new ReedMullerCode(1n, 3n);
    const punctured = rm.puncture([0, 1]);

    expect(punctured.length()).toBe(6);
    expect(punctured.dimension()).toBe(4);
  });

  it('should encode correctly', () => {
    const rm = new ReedMullerCode(1n, 3n);
    const punctured = rm.puncture([0, 1]);

    const message = [1, 0, 1, 0];
    const codeword = punctured.encode(message);

    expect(codeword.length).toBe(6);

    // Verify against full encoding
    const fullCodeword = rm.encode(message);
    let pIdx = 0;
    for (let i = 0; i < 8; i++) {
      if (i !== 0 && i !== 1) {
        expect(codeword[pIdx]!.value).toBe(fullCodeword[i]!.value);
        pIdx++;
      }
    }
  });

  it('should throw for out of bounds positions', () => {
    const rm = new ReedMullerCode(1n, 3n);
    expect(() => rm.puncture([0, 8])).toThrow(ValueError);
    expect(() => rm.puncture([-1])).toThrow(ValueError);
  });

  it('should return punctured positions', () => {
    const rm = new ReedMullerCode(1n, 3n);
    const punctured = rm.puncture([2, 5]);

    expect(punctured.punctured_positions()).toEqual([2, 5]);
  });

  it('should return base code', () => {
    const rm = new ReedMullerCode(1n, 3n);
    const punctured = rm.puncture([0]);

    expect(punctured.base_code()).toBe(rm);
  });
});

describe('Factory functions', () => {
  it('BinaryReedMullerCode should create correct code', () => {
    const rm = BinaryReedMullerCode(2n, 4n);

    expect(rm.order()).toBe(2);
    expect(rm.num_variables()).toBe(4);
    expect(rm).toBeInstanceOf(ReedMullerCode);
  });

  it('FirstOrderReedMullerCode should create RM(1, m)', () => {
    const rm = FirstOrderReedMullerCode(5n);

    expect(rm.order()).toBe(1);
    expect(rm.num_variables()).toBe(5);
    expect(rm.is_first_order()).toBe(true);
  });

  it('RepetitionCode should create RM(0, m)', () => {
    const rm = RepetitionCode(4n);

    expect(rm.order()).toBe(0);
    expect(rm.num_variables()).toBe(4);
    expect(rm.is_repetition_code()).toBe(true);
    expect(rm.dimension()).toBe(1);
    expect(rm.minimum_distance()).toBe(16);
  });
});

describe('edge cases', () => {
  it('should handle RM(0, 0) - single bit code', () => {
    const rm = new ReedMullerCode(0n, 0n);

    expect(rm.length()).toBe(1);
    expect(rm.dimension()).toBe(1);
    expect(rm.minimum_distance()).toBe(1);

    const codeword = rm.encode([1]);
    expect(codeword.length).toBe(1);
    expect(codeword[0]!.value).toBe(1);
  });

  it('should handle RM(1, 1)', () => {
    const rm = new ReedMullerCode(1n, 1n);

    expect(rm.length()).toBe(2);
    expect(rm.dimension()).toBe(2);
    expect(rm.minimum_distance()).toBe(1);
  });

  it('should handle larger codes RM(2, 5)', () => {
    const rm = new ReedMullerCode(2n, 5n);

    expect(rm.length()).toBe(32);
    expect(rm.dimension()).toBe(16); // 1 + 5 + 10 = 16
    expect(rm.minimum_distance()).toBe(8); // 2^(5-2) = 8
  });

  it('toString should return meaningful description', () => {
    const rm = new ReedMullerCode(1n, 3n);
    const str = rm.toString();

    expect(str).toContain('Reed-Muller');
    expect(str).toContain('1');
    expect(str).toContain('3');
    expect(str).toContain('[8, 4, 4]');
  });

  it('eq should compare codes correctly', () => {
    const rm1 = new ReedMullerCode(1n, 3n);
    const rm2 = new ReedMullerCode(1n, 3n);
    const rm3 = new ReedMullerCode(2n, 3n);

    expect(rm1.eq(rm2)).toBe(true);
    expect(rm1.eq(rm3)).toBe(false);
  });
});
