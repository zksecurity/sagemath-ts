/**
 * Tests for Reed-Solomon codes for ZK applications (FRI, STARKs)
 */

import { describe, expect, it } from 'vitest';
import {
  ReedSolomonCode,
  DecodingError,
  createClassicalReedSolomonCode,
  createFRIReedSolomonCode,
  type FieldElement,
  type FiniteField,
} from './reed_solomon.js';
import {
  FiniteFieldElement,
  FiniteFieldPrime,
} from '../rings/finite_rings/finite_field_prime.js';

// Helper to create GF(p) fields
function GF(p: bigint): FiniteFieldPrime {
  return new FiniteFieldPrime(p, true);
}

describe('ReedSolomonCode', () => {
  describe('constructor', () => {
    it('should create a valid RS code with default evaluation points', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8, 4);

      expect(rs.n).toBe(8);
      expect(rs.k).toBe(4);
      expect(rs.field).toBe(F);
      expect(rs.evaluation_points().length).toBe(8);
    });

    it('should create a valid RS code with custom evaluation points', () => {
      const F = GF(17n);
      const evalPoints = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
      const rs = new ReedSolomonCode(F, 4, 2, evalPoints);

      expect(rs.n).toBe(4);
      expect(rs.k).toBe(2);
      expect(rs.evaluation_points()).toEqual(evalPoints);
    });

    it('should throw for k > n', () => {
      const F = GF(17n);
      expect(() => new ReedSolomonCode(F, 4, 5)).toThrow(ValueError);
    });

    it('should throw for k < 1', () => {
      const F = GF(17n);
      expect(() => new ReedSolomonCode(F, 4, 0)).toThrow(ValueError);
    });

    it('should throw for n > field order', () => {
      const F = GF(7n);
      expect(() => new ReedSolomonCode(F, 10, 5)).toThrow(ValueError);
    });

    it('should throw for non-distinct evaluation points', () => {
      const F = GF(17n);
      const evalPoints = [F.__call__(1), F.__call__(2), F.__call__(1), F.__call__(4)];
      expect(() => new ReedSolomonCode(F, 4, 2, evalPoints)).toThrow(ValueError);
    });

    it('should throw for wrong number of evaluation points', () => {
      const F = GF(17n);
      const evalPoints = [F.__call__(1), F.__call__(2), F.__call__(3)];
      expect(() => new ReedSolomonCode(F, 4, 2, evalPoints)).toThrow(ValueError);
    });
  });

  describe('code parameters', () => {
    it('should compute minimum distance correctly', () => {
      const F = GF(17n);

      // RS(8, 4) has d = 8 - 4 + 1 = 5
      const rs84 = new ReedSolomonCode(F, 8, 4);
      expect(rs84.minimum_distance()).toBe(5);

      // RS(10, 3) has d = 10 - 3 + 1 = 8
      const rs103 = new ReedSolomonCode(F, 10, 3);
      expect(rs103.minimum_distance()).toBe(8);

      // RS(5, 5) has d = 5 - 5 + 1 = 1
      const rs55 = new ReedSolomonCode(F, 5, 5);
      expect(rs55.minimum_distance()).toBe(1);
    });

    it('should compute decoding radius correctly', () => {
      const F = GF(17n);

      // RS(8, 4) can correct floor((8-4)/2) = 2 errors
      const rs84 = new ReedSolomonCode(F, 8, 4);
      expect(rs84.decoding_radius()).toBe(2);

      // RS(10, 3) can correct floor((10-3)/2) = 3 errors
      const rs103 = new ReedSolomonCode(F, 10, 3);
      expect(rs103.decoding_radius()).toBe(3);

      // RS(7, 6) can correct floor((7-6)/2) = 0 errors
      const rs76 = new ReedSolomonCode(F, 7, 6);
      expect(rs76.decoding_radius()).toBe(0);
    });

    it('should compute rate correctly', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8, 4);
      expect(rs.rate()).toBe(0.5);
    });

    it('should compute redundancy correctly', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8, 4);
      expect(rs.redundancy()).toBe(4);
    });
  });

  describe('encode', () => {
    it('should encode a message correctly', () => {
      const F = GF(7n);
      const evalPoints = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
      const rs = new ReedSolomonCode(F, 4, 2, evalPoints);

      // Message: [a, b] represents polynomial f(x) = a + b*x
      // f(1) = a + b, f(2) = a + 2b, f(3) = a + 3b, f(4) = a + 4b
      const message = [F.__call__(3), F.__call__(2)]; // f(x) = 3 + 2x
      const codeword = rs.encode(message);

      expect(codeword.length).toBe(4);
      expect(codeword[0]!.eq(F.__call__(5))).toBe(true);  // 3 + 2*1 = 5
      expect(codeword[1]!.eq(F.__call__(0))).toBe(true);  // 3 + 2*2 = 7 ≡ 0 (mod 7)
      expect(codeword[2]!.eq(F.__call__(2))).toBe(true);  // 3 + 2*3 = 9 ≡ 2 (mod 7)
      expect(codeword[3]!.eq(F.__call__(4))).toBe(true);  // 3 + 2*4 = 11 ≡ 4 (mod 7)
    });

    it('should encode the zero message to the zero codeword', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8, 4);
      const message = [F.zero(), F.zero(), F.zero(), F.zero()];
      const codeword = rs.encode(message);

      expect(codeword.every(c => c.isZero())).toBe(true);
    });

    it('should throw for wrong message length', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8, 4);
      const message = [F.__call__(1), F.__call__(2), F.__call__(3)]; // length 3, not 4

      expect(() => rs.encode(message)).toThrow(ValueError);
    });
  });

  describe('encode/decode round-trip', () => {
    it('should decode a valid codeword correctly', () => {
      const F = GF(17n);
      const evalPoints: FiniteFieldElement[] = [];
      for (let i = 1; i <= 8; i++) {
        evalPoints.push(F.__call__(i));
      }
      const rs = new ReedSolomonCode(F, 8, 4, evalPoints);

      const message = [F.__call__(5), F.__call__(3), F.__call__(7), F.__call__(2)];
      const codeword = rs.encode(message);
      const decoded = rs.decode(codeword);

      expect(decoded.length).toBe(4);
      for (let i = 0; i < 4; i++) {
        expect(decoded[i]!.eq(message[i]!)).toBe(true);
      }
    });

    it('should decode with one error correctly', () => {
      const F = GF(17n);
      const evalPoints: FiniteFieldElement[] = [];
      for (let i = 1; i <= 8; i++) {
        evalPoints.push(F.__call__(i));
      }
      const rs = new ReedSolomonCode(F, 8, 4, evalPoints); // can correct 2 errors

      const message = [F.__call__(5), F.__call__(3), F.__call__(7), F.__call__(2)];
      const codeword = rs.encode(message);

      // Introduce one error
      const corrupted = [...codeword];
      corrupted[3] = F.__call__(99); // corrupt position 3

      const decoded = rs.decode(corrupted);

      expect(decoded.length).toBe(4);
      for (let i = 0; i < 4; i++) {
        expect(decoded[i]!.eq(message[i]!)).toBe(true);
      }
    });

    it('should decode with maximum correctable errors', () => {
      const F = GF(17n);
      const evalPoints: FiniteFieldElement[] = [];
      for (let i = 1; i <= 10; i++) {
        evalPoints.push(F.__call__(i));
      }
      const rs = new ReedSolomonCode(F, 10, 4, evalPoints); // can correct 3 errors

      const message = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
      const codeword = rs.encode(message);

      // Introduce 3 errors (the maximum)
      const corrupted = [...codeword];
      corrupted[0] = F.__call__(99);
      corrupted[5] = F.__call__(88);
      corrupted[9] = F.__call__(77);

      const decoded = rs.decode(corrupted);

      expect(decoded.length).toBe(4);
      for (let i = 0; i < 4; i++) {
        expect(decoded[i]!.eq(message[i]!)).toBe(true);
      }
    });

    it('should fail to decode with too many errors', () => {
      const F = GF(17n);
      const evalPoints: FiniteFieldElement[] = [];
      for (let i = 1; i <= 8; i++) {
        evalPoints.push(F.__call__(i));
      }
      const rs = new ReedSolomonCode(F, 8, 4, evalPoints); // can correct 2 errors

      const message = [F.__call__(5), F.__call__(3), F.__call__(7), F.__call__(2)];
      const codeword = rs.encode(message);

      // Introduce 3 errors (one too many)
      const corrupted = [...codeword];
      corrupted[0] = F.__call__(99);
      corrupted[3] = F.__call__(88);
      corrupted[7] = F.__call__(77);

      expect(() => rs.decode(corrupted)).toThrow(DecodingError);
    });
  });

  describe('syndrome', () => {
    it('should compute zero syndromes for valid codewords', () => {
      const F = GF(17n);
      const evalPoints: FiniteFieldElement[] = [];
      for (let i = 1; i <= 8; i++) {
        evalPoints.push(F.__call__(i));
      }
      const rs = new ReedSolomonCode(F, 8, 4, evalPoints);

      const message = [F.__call__(5), F.__call__(3), F.__call__(7), F.__call__(2)];
      const codeword = rs.encode(message);
      const syndromes = rs.syndrome(codeword);

      // For a valid codeword, all syndromes should be zero
      expect(syndromes.length).toBe(4); // n - k = 8 - 4 = 4
      expect(syndromes.every(s => s.isZero())).toBe(true);
    });

    it('should compute non-zero syndromes for corrupted words', () => {
      const F = GF(17n);
      const evalPoints: FiniteFieldElement[] = [];
      for (let i = 1; i <= 8; i++) {
        evalPoints.push(F.__call__(i));
      }
      const rs = new ReedSolomonCode(F, 8, 4, evalPoints);

      const message = [F.__call__(5), F.__call__(3), F.__call__(7), F.__call__(2)];
      const codeword = rs.encode(message);

      // Corrupt the codeword
      const corrupted = [...codeword];
      corrupted[2] = F.__call__(99);

      const syndromes = rs.syndrome(corrupted);

      // At least one syndrome should be non-zero
      expect(syndromes.some(s => !s.isZero())).toBe(true);
    });
  });

  describe('FRI operations', () => {
    describe('fold', () => {
      it('should fold a codeword to half length', () => {
        const F = GF(17n);
        const rs = new ReedSolomonCode(F, 8, 4);

        const message = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
        const codeword = rs.encode(message);
        const challenge = F.__call__(5);

        const folded = rs.fold(codeword, challenge);

        expect(folded.length).toBe(4);
      });

      it('should throw for odd-length codewords', () => {
        const F = GF(17n);
        const rs = new ReedSolomonCode(F, 7, 3);

        const message = [F.__call__(1), F.__call__(2), F.__call__(3)];
        const codeword = rs.encode(message);
        const challenge = F.__call__(5);

        expect(() => rs.fold(codeword, challenge)).toThrow(ValueError);
      });

      it('should produce consistent results with different challenges', () => {
        const F = GF(17n);
        const rs = new ReedSolomonCode(F, 8, 4);

        const message = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
        const codeword = rs.encode(message);

        const folded1 = rs.fold(codeword, F.__call__(3));
        const folded2 = rs.fold(codeword, F.__call__(7));

        // Different challenges should give different results
        expect(folded1[0]!.eq(folded2[0]!)).toBe(false);
      });
    });

    describe('query', () => {
      it('should return correct values at queried indices', () => {
        const F = GF(17n);
        const rs = new ReedSolomonCode(F, 8, 4);

        const message = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
        const codeword = rs.encode(message);

        const queries = rs.query(codeword, [0, 3, 7]);

        expect(queries.length).toBe(3);
        expect(queries[0]![0]).toBe(0);
        expect(queries[0]![1].eq(codeword[0]!)).toBe(true);
        expect(queries[1]![0]).toBe(3);
        expect(queries[1]![1].eq(codeword[3]!)).toBe(true);
        expect(queries[2]![0]).toBe(7);
        expect(queries[2]![1].eq(codeword[7]!)).toBe(true);
      });

      it('should throw for out-of-bounds indices', () => {
        const F = GF(17n);
        const rs = new ReedSolomonCode(F, 8, 4);

        const message = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
        const codeword = rs.encode(message);

        expect(() => rs.query(codeword, [0, 8])).toThrow(ValueError);
        expect(() => rs.query(codeword, [-1])).toThrow(ValueError);
      });
    });

    describe('isClose', () => {
      it('should return true for valid codewords', () => {
        const F = GF(17n);
        const evalPoints: FiniteFieldElement[] = [];
        for (let i = 1; i <= 8; i++) {
          evalPoints.push(F.__call__(i));
        }
        const rs = new ReedSolomonCode(F, 8, 4, evalPoints);

        const message = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
        const codeword = rs.encode(message);

        expect(rs.isClose(codeword, 0)).toBe(true);
        expect(rs.isClose(codeword, 2)).toBe(true);
      });

      it('should return true for words within error-correction capability', () => {
        const F = GF(17n);
        const evalPoints: FiniteFieldElement[] = [];
        for (let i = 1; i <= 8; i++) {
          evalPoints.push(F.__call__(i));
        }
        const rs = new ReedSolomonCode(F, 8, 4, evalPoints); // can correct 2 errors

        const message = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
        const codeword = rs.encode(message);

        const corrupted = [...codeword];
        corrupted[0] = F.__call__(99);
        corrupted[5] = F.__call__(88);

        expect(rs.isClose(corrupted, 2)).toBe(true);
      });

      it('should return false for words too far from any codeword', () => {
        const F = GF(17n);
        const evalPoints: FiniteFieldElement[] = [];
        for (let i = 1; i <= 8; i++) {
          evalPoints.push(F.__call__(i));
        }
        const rs = new ReedSolomonCode(F, 8, 4, evalPoints);

        const message = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
        const codeword = rs.encode(message);

        const corrupted = [...codeword];
        corrupted[0] = F.__call__(99);
        corrupted[1] = F.__call__(88);
        corrupted[2] = F.__call__(77);
        corrupted[3] = F.__call__(66);

        expect(rs.isClose(corrupted, 2)).toBe(false);
      });
    });
  });

  describe('hammingDistance', () => {
    it('should compute Hamming distance correctly', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8, 4);

      const a = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4),
                 F.__call__(5), F.__call__(6), F.__call__(7), F.__call__(8)];
      const b = [F.__call__(1), F.__call__(2), F.__call__(99), F.__call__(4),
                 F.__call__(5), F.__call__(88), F.__call__(7), F.__call__(77)];

      expect(rs.hammingDistance(a, b)).toBe(3);
    });

    it('should return 0 for identical words', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8, 4);

      const a = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4),
                 F.__call__(5), F.__call__(6), F.__call__(7), F.__call__(8)];

      expect(rs.hammingDistance(a, a)).toBe(0);
    });
  });
});

describe('createClassicalReedSolomonCode', () => {
  it('should create a code with consecutive powers of primitive root', () => {
    const F = GF(17n); // 17 - 1 = 16, divisible by 8
    const rs = createClassicalReedSolomonCode(F, 8, 4);

    expect(rs.n).toBe(8);
    expect(rs.k).toBe(4);

    // Verify evaluation points are powers of a primitive 8th root of unity
    const evalPoints = rs.evaluation_points();
    expect(evalPoints.length).toBe(8);

    // Check omega^8 = 1 (using the ratio between consecutive points)
    const omega = evalPoints[1]!.mul(evalPoints[0]!.inv());
    expect(omega.pow(8).eq(F.one())).toBe(true);
  });

  it('should throw for invalid length', () => {
    const F = GF(17n); // 17 - 1 = 16, not divisible by 7
    expect(() => createClassicalReedSolomonCode(F, 7, 3)).toThrow(ValueError);
  });
});

describe('createFRIReedSolomonCode', () => {
  it('should create a code suitable for FRI', () => {
    const F = GF(17n); // 17 - 1 = 16 = 2^4
    const rs = createFRIReedSolomonCode(F, 8, 4);

    expect(rs.n).toBe(8);
    expect(rs.k).toBe(4);
    expect(rs.rate()).toBe(0.5);
  });

  it('should throw for non-power-of-2 length', () => {
    const F = GF(17n);
    expect(() => createFRIReedSolomonCode(F, 6, 3)).toThrow(ValueError);
  });
});

describe('edge cases and special scenarios', () => {
  it('should handle rate-1 code (n = k)', () => {
    const F = GF(17n);
    const evalPoints: FiniteFieldElement[] = [];
    for (let i = 1; i <= 4; i++) {
      evalPoints.push(F.__call__(i));
    }
    const rs = new ReedSolomonCode(F, 4, 4, evalPoints);

    expect(rs.minimum_distance()).toBe(1);
    expect(rs.decoding_radius()).toBe(0);

    const message = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
    const codeword = rs.encode(message);
    const decoded = rs.decode(codeword);

    for (let i = 0; i < 4; i++) {
      expect(decoded[i]!.eq(message[i]!)).toBe(true);
    }
  });

  it('should handle single-symbol message (k = 1)', () => {
    const F = GF(17n);
    const evalPoints: FiniteFieldElement[] = [];
    for (let i = 1; i <= 8; i++) {
      evalPoints.push(F.__call__(i));
    }
    const rs = new ReedSolomonCode(F, 8, 1, evalPoints);

    expect(rs.minimum_distance()).toBe(8);
    expect(rs.decoding_radius()).toBe(3);

    const message = [F.__call__(5)];
    const codeword = rs.encode(message);

    // All codeword values should be 5 (constant polynomial)
    expect(codeword.every(c => c.eq(F.__call__(5)))).toBe(true);
  });

  it('should handle large field', () => {
    const F = GF(65537n); // Fermat prime
    const evalPoints: FiniteFieldElement[] = [];
    for (let i = 1; i <= 32; i++) {
      evalPoints.push(F.__call__(i));
    }
    const rs = new ReedSolomonCode(F, 32, 16, evalPoints);

    const message: FiniteFieldElement[] = [];
    for (let i = 0; i < 16; i++) {
      message.push(F.__call__(i * 1000 + 1));
    }

    const codeword = rs.encode(message);
    const decoded = rs.decode(codeword);

    for (let i = 0; i < 16; i++) {
      expect(decoded[i]!.eq(message[i]!)).toBe(true);
    }
  });
});

// Import ValueError for test assertions
import { ValueError } from '../errors.js';
