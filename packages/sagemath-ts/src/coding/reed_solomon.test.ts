/**
 * Tests for Reed-Solomon codes for ZK applications (FRI, STARKs)
 */

import { describe, expect, it } from 'vitest';
import {
  type FiniteFieldElement,
  FiniteFieldPrime,
} from '../rings/finite_rings/finite_field_prime.js';
import {
  DecodingError,
  type FieldElement,
  type FiniteField,
  ReedSolomonCode,
  createClassicalReedSolomonCode,
  createFRIReedSolomonCode,
} from './reed_solomon.js';

// Helper to create GF(p) fields
function GF(p: bigint): FiniteFieldPrime {
  return new FiniteFieldPrime(p, true);
}

describe('ReedSolomonCode', () => {
  describe('constructor', () => {
    it('should create a valid RS code with default evaluation points', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8n, 4n);

      expect(rs.n).toBe(8);
      expect(rs.k).toBe(4);
      expect(rs.field).toBe(F);
      expect(rs.evaluation_points().length).toBe(8);
    });

    it('should create a valid RS code with custom evaluation points', () => {
      const F = GF(17n);
      const evalPoints = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
      const rs = new ReedSolomonCode(F, 4n, 2n, evalPoints);

      expect(rs.n).toBe(4);
      expect(rs.k).toBe(2);
      expect(rs.evaluation_points()).toEqual(evalPoints);
    });

    it('should throw for k > n', () => {
      const F = GF(17n);
      expect(() => new ReedSolomonCode(F, 4n, 5n)).toThrow(ValueError);
    });

    it('should throw for k < 1', () => {
      const F = GF(17n);
      expect(() => new ReedSolomonCode(F, 4n, 0n)).toThrow(ValueError);
    });

    it('should throw for n > field order', () => {
      const F = GF(7n);
      expect(() => new ReedSolomonCode(F, 10n, 5n)).toThrow(ValueError);
    });

    it('should throw for non-distinct evaluation points', () => {
      const F = GF(17n);
      const evalPoints = [F.__call__(1), F.__call__(2), F.__call__(1), F.__call__(4)];
      expect(() => new ReedSolomonCode(F, 4n, 2n, evalPoints)).toThrow(ValueError);
    });

    it('should throw for wrong number of evaluation points', () => {
      const F = GF(17n);
      const evalPoints = [F.__call__(1), F.__call__(2), F.__call__(3)];
      expect(() => new ReedSolomonCode(F, 4n, 2n, evalPoints)).toThrow(ValueError);
    });
  });

  describe('code parameters', () => {
    it('should compute minimum distance correctly', () => {
      const F = GF(17n);

      // RS(8, 4) has d = 8 - 4 + 1 = 5
      const rs84 = new ReedSolomonCode(F, 8n, 4n);
      expect(rs84.minimum_distance()).toBe(5);

      // RS(10, 3) has d = 10 - 3 + 1 = 8
      const rs103 = new ReedSolomonCode(F, 10n, 3n);
      expect(rs103.minimum_distance()).toBe(8);

      // RS(5, 5) has d = 5 - 5 + 1 = 1
      const rs55 = new ReedSolomonCode(F, 5n, 5n);
      expect(rs55.minimum_distance()).toBe(1);
    });

    it('should compute decoding radius correctly', () => {
      const F = GF(17n);

      // RS(8, 4) can correct floor((8-4)/2) = 2 errors
      const rs84 = new ReedSolomonCode(F, 8n, 4n);
      expect(rs84.decoding_radius()).toBe(2);

      // RS(10, 3) can correct floor((10-3)/2) = 3 errors
      const rs103 = new ReedSolomonCode(F, 10n, 3n);
      expect(rs103.decoding_radius()).toBe(3);

      // RS(7, 6) can correct floor((7-6)/2) = 0 errors
      const rs76 = new ReedSolomonCode(F, 7n, 6n);
      expect(rs76.decoding_radius()).toBe(0);
    });

    it('should compute rate correctly', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8n, 4n);
      expect(rs.rate()).toBe(0.5);
    });

    it('should compute redundancy correctly', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8n, 4n);
      expect(rs.redundancy()).toBe(4);
    });
  });

  describe('encode', () => {
    it('should encode a message correctly', () => {
      const F = GF(7n);
      const evalPoints = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
      const rs = new ReedSolomonCode(F, 4n, 2n, evalPoints);

      // Message: [a, b] represents polynomial f(x) = a + b*x
      // f(1) = a + b, f(2) = a + 2b, f(3) = a + 3b, f(4) = a + 4b
      const message = [F.__call__(3), F.__call__(2)]; // f(x) = 3 + 2x
      const codeword = rs.encode(message);

      expect(codeword.length).toBe(4);
      expect(codeword[0]!.eq(F.__call__(5))).toBe(true); // 3 + 2*1 = 5
      expect(codeword[1]!.eq(F.__call__(0))).toBe(true); // 3 + 2*2 = 7 ≡ 0 (mod 7)
      expect(codeword[2]!.eq(F.__call__(2))).toBe(true); // 3 + 2*3 = 9 ≡ 2 (mod 7)
      expect(codeword[3]!.eq(F.__call__(4))).toBe(true); // 3 + 2*4 = 11 ≡ 4 (mod 7)
    });

    it('should encode the zero message to the zero codeword', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8n, 4n);
      const message = [F.zero(), F.zero(), F.zero(), F.zero()];
      const codeword = rs.encode(message);

      expect(codeword.every((c) => c.isZero())).toBe(true);
    });

    it('should throw for wrong message length', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8n, 4n);
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
      const rs = new ReedSolomonCode(F, 8n, 4n, evalPoints);

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
      const rs = new ReedSolomonCode(F, 8n, 4n, evalPoints); // can correct 2 errors

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
      const rs = new ReedSolomonCode(F, 10n, 4n, evalPoints); // can correct 3 errors

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
      const rs = new ReedSolomonCode(F, 8n, 4n, evalPoints); // can correct 2 errors

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
      const rs = new ReedSolomonCode(F, 8n, 4n, evalPoints);

      const message = [F.__call__(5), F.__call__(3), F.__call__(7), F.__call__(2)];
      const codeword = rs.encode(message);
      const syndromes = rs.syndrome(codeword);

      // For a valid codeword, all syndromes should be zero
      expect(syndromes.length).toBe(4); // n - k = 8 - 4 = 4
      expect(syndromes.every((s) => s.isZero())).toBe(true);
    });

    it('should compute non-zero syndromes for corrupted words', () => {
      const F = GF(17n);
      const evalPoints: FiniteFieldElement[] = [];
      for (let i = 1; i <= 8; i++) {
        evalPoints.push(F.__call__(i));
      }
      const rs = new ReedSolomonCode(F, 8n, 4n, evalPoints);

      const message = [F.__call__(5), F.__call__(3), F.__call__(7), F.__call__(2)];
      const codeword = rs.encode(message);

      // Corrupt the codeword
      const corrupted = [...codeword];
      corrupted[2] = F.__call__(99);

      const syndromes = rs.syndrome(corrupted);

      // At least one syndrome should be non-zero
      expect(syndromes.some((s) => !s.isZero())).toBe(true);
    });
  });

  describe('FRI operations', () => {
    describe('fold', () => {
      it('should fold a codeword to half length', () => {
        const F = GF(17n);
        const rs = new ReedSolomonCode(F, 8n, 4n);

        const message = [F.__call__(1), F.__call__(2), F.__call__(3), F.__call__(4)];
        const codeword = rs.encode(message);
        const challenge = F.__call__(5);

        const folded = rs.fold(codeword, challenge);

        expect(folded.length).toBe(4);
      });

      it('should throw for odd-length codewords', () => {
        const F = GF(17n);
        const rs = new ReedSolomonCode(F, 7n, 3n);

        const message = [F.__call__(1), F.__call__(2), F.__call__(3)];
        const codeword = rs.encode(message);
        const challenge = F.__call__(5);

        expect(() => rs.fold(codeword, challenge)).toThrow(ValueError);
      });

      it('should produce consistent results with different challenges', () => {
        const F = GF(17n);
        const rs = new ReedSolomonCode(F, 8n, 4n);

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
        const rs = new ReedSolomonCode(F, 8n, 4n);

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
        const rs = new ReedSolomonCode(F, 8n, 4n);

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
        const rs = new ReedSolomonCode(F, 8n, 4n, evalPoints);

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
        const rs = new ReedSolomonCode(F, 8n, 4n, evalPoints); // can correct 2 errors

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
        const rs = new ReedSolomonCode(F, 8n, 4n, evalPoints);

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
      const rs = new ReedSolomonCode(F, 8n, 4n);

      const a = [
        F.__call__(1),
        F.__call__(2),
        F.__call__(3),
        F.__call__(4),
        F.__call__(5),
        F.__call__(6),
        F.__call__(7),
        F.__call__(8),
      ];
      const b = [
        F.__call__(1),
        F.__call__(2),
        F.__call__(99),
        F.__call__(4),
        F.__call__(5),
        F.__call__(88),
        F.__call__(7),
        F.__call__(77),
      ];

      expect(rs.hammingDistance(a, b)).toBe(3);
    });

    it('should return 0 for identical words', () => {
      const F = GF(17n);
      const rs = new ReedSolomonCode(F, 8n, 4n);

      const a = [
        F.__call__(1),
        F.__call__(2),
        F.__call__(3),
        F.__call__(4),
        F.__call__(5),
        F.__call__(6),
        F.__call__(7),
        F.__call__(8),
      ];

      expect(rs.hammingDistance(a, a)).toBe(0);
    });
  });
});

describe('createClassicalReedSolomonCode', () => {
  it('should create a code with consecutive powers of primitive root', () => {
    const F = GF(17n); // 17 - 1 = 16, divisible by 8
    const rs = createClassicalReedSolomonCode(F, 8n, 4n);

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
    expect(() => createClassicalReedSolomonCode(F, 7n, 3n)).toThrow(ValueError);
  });
});

describe('createFRIReedSolomonCode', () => {
  it('should create a code suitable for FRI', () => {
    const F = GF(17n); // 17 - 1 = 16 = 2^4
    const rs = createFRIReedSolomonCode(F, 8n, 4n);

    expect(rs.n).toBe(8);
    expect(rs.k).toBe(4);
    expect(rs.rate()).toBe(0.5);
  });

  it('should throw for non-power-of-2 length', () => {
    const F = GF(17n);
    expect(() => createFRIReedSolomonCode(F, 6n, 3n)).toThrow(ValueError);
  });
});

describe('edge cases and special scenarios', () => {
  it('should handle rate-1 code (n = k)', () => {
    const F = GF(17n);
    const evalPoints: FiniteFieldElement[] = [];
    for (let i = 1; i <= 4; i++) {
      evalPoints.push(F.__call__(i));
    }
    const rs = new ReedSolomonCode(F, 4n, 4n, evalPoints);

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
    const rs = new ReedSolomonCode(F, 8n, 1n, evalPoints);

    expect(rs.minimum_distance()).toBe(8);
    expect(rs.decoding_radius()).toBe(3);

    const message = [F.__call__(5)];
    const codeword = rs.encode(message);

    // All codeword values should be 5 (constant polynomial)
    expect(codeword.every((c) => c.eq(F.__call__(5)))).toBe(true);
  });

  it('should handle large field', () => {
    const F = GF(65537n); // Fermat prime
    const evalPoints: FiniteFieldElement[] = [];
    for (let i = 1; i <= 32; i++) {
      evalPoints.push(F.__call__(i));
    }
    const rs = new ReedSolomonCode(F, 32n, 16n, evalPoints);

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

describe('key equation decoding chain (SageMath GRSKeyEquationSyndromeDecoder)', () => {
  // Sage:
  //   F = GF(11); n, k = 10, 5
  //   C = codes.GeneralizedReedSolomonCode(F.list()[1:n+1], k)
  //   D = codes.decoders.GRSKeyEquationSyndromeDecoder(C)
  const F = GF(11n);
  const evalPoints = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10].map((v) => F.__call__(v));
  const code = () => new ReedSolomonCode(F, 10n, 5n, evalPoints);

  it('should reproduce parity_column_multipliers', () => {
    // sage: C.parity_column_multipliers()
    // [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    const eta = code().parity_column_multipliers();
    expect(eta.map((e) => e.toString())).toEqual([
      '10',
      '9',
      '8',
      '7',
      '6',
      '5',
      '4',
      '3',
      '2',
      '1',
    ]);
  });

  it('should reproduce the _syndrome doctest', () => {
    // sage: r = vector(F, (8, 2, 6, 10, 6, 10, 7, 6, 7, 2))
    // sage: D._syndrome(r)
    // [1, 10, 1, 10, 1]
    const r = [8, 2, 6, 10, 6, 10, 7, 6, 7, 2].map((v) => F.__call__(v));
    const S = code().syndrome(r);
    expect(S.map((e) => e.toString())).toEqual(['1', '10', '1', '10', '1']);
  });

  it('should reproduce the _forney_formula doctest', () => {
    // sage: R.<x> = F[]
    // sage: evaluator, locator = R(10), R([10, 10])
    // sage: D._forney_formula(evaluator, locator)
    // (0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
    const C = code();
    const R = new PolynomialRing<FiniteFieldElement>(F, 'x');
    const evaluator = new Polynomial<FiniteFieldElement>([F.__call__(10)], R);
    const locator = new Polynomial<FiniteFieldElement>([F.__call__(10), F.__call__(10)], R);

    const errors = C.forney_algorithm(locator, evaluator);
    const vector = Array.from({ length: 10 }, (_, i) =>
      errors.has(i) ? errors.get(i)!.toString() : '0'
    );
    expect(vector).toEqual(['0', '0', '0', '0', '0', '0', '0', '0', '0', '1']);
  });

  it('should locate and evaluate a single error exactly', () => {
    // A [10, 4] code over GF(11): put a known error at position 2 and check
    // that the locator/evaluator/Forney chain reports position 2 (not the
    // reflected position 8) with the exact error value.
    const C = new ReedSolomonCode(F, 10n, 4n, evalPoints);
    const message = [F.__call__(1), F.__call__(4), F.__call__(7), F.__call__(10)];
    const codeword = C.encode(message);

    const received = [...codeword];
    received[2] = received[2]!.add(F.__call__(5));

    const S = C.syndrome(received);
    const locator = C.error_locator(S);
    const evaluator = C.error_evaluator(S, locator);
    const errors = C.forney_algorithm(locator, evaluator);

    expect([...errors.keys()]).toEqual([2]);
    expect(errors.get(2)!.toString()).toBe('5');
  });

  it('should correct every error pattern of weight <= t through the chain', () => {
    const cases: Array<[bigint, number, number]> = [
      [11n, 10, 4],
      [13n, 12, 6],
      [17n, 16, 10],
    ];

    for (const [p, n, k] of cases) {
      const K = GF(p);
      const points = Array.from({ length: n }, (_, i) => K.__call__(i + 1));
      const C = new ReedSolomonCode(K, BigInt(n), BigInt(k), points);
      const t = C.decoding_radius();

      const message = Array.from({ length: k }, (_, i) =>
        K.__call__(Number((BigInt(i) * 3n + 1n) % p))
      );
      const codeword = C.encode(message);

      const subsets: number[][] = [];
      const cur: number[] = [];
      const rec = (start: number, size: number) => {
        if (cur.length === size) {
          subsets.push([...cur]);
          return;
        }
        for (let i = start; i < n; i++) {
          cur.push(i);
          rec(i + 1, size);
          cur.pop();
        }
      };
      for (let w = 1; w <= t; w++) rec(0, w);

      for (const positions of subsets) {
        for (let v = 1n; v < p; v++) {
          const received = [...codeword];
          for (const pos of positions) received[pos] = received[pos]!.add(K.__call__(Number(v)));

          const S = C.syndrome(received);
          const locator = C.error_locator(S);
          const evaluator = C.error_evaluator(S, locator);
          const errors = C.forney_algorithm(locator, evaluator);

          expect([...errors.keys()].sort((a, b) => a - b)).toEqual(positions);
          for (const pos of positions) {
            expect(errors.get(pos)!.toString()).toBe(String(v));
          }
        }
      }
    }
  }, 60_000);
});

describe('FRI fold domain', () => {
  it('should fold twice over the squared domain', () => {
    // GF(17) has a primitive 4th root of unity; the domain is [1, 13, 16, 4].
    // f(x) = 1 + 2x + 3x^2 + 4x^3 splits as f0(y) = 1 + 3y, f1(y) = 2 + 4y.
    // Folding with challenge 5 gives g(y) = f0 + 5*f1 = 11 + 6y on the squared
    // domain [1, 16]; folding g with challenge 7 gives 11 + 7*6 = 2, constant.
    const F = GF(17n);
    const rs = createClassicalReedSolomonCode(F, 4n, 2n);
    const domain = rs.evaluation_points();
    expect(domain.map((e) => e.toString())).toEqual(['1', '13', '16', '4']);

    const f = (x: FiniteFieldElement) =>
      F.__call__(1)
        .add(F.__call__(2).mul(x))
        .add(F.__call__(3).mul(x.mul(x)))
        .add(F.__call__(4).mul(x.mul(x).mul(x)));
    const codeword = domain.map((x) => f(x));

    const first = rs.fold(codeword, F.__call__(5));
    const squared = rs.fold_domain(2);
    expect(squared.map((e) => e.toString())).toEqual(['1', '16']);
    const expectedFirst = squared.map((y) => F.__call__(11).add(F.__call__(6).mul(y)));
    expect(first.map((e) => e.toString())).toEqual(expectedFirst.map((e) => e.toString()));

    // Second fold must use the squared domain, not the top-level one.
    const second = rs.fold(first, F.__call__(7));
    expect(second.map((e) => e.toString())).toEqual(['2']);
  });

  it('should accept an explicit domain', () => {
    const F = GF(17n);
    const rs = createClassicalReedSolomonCode(F, 4n, 2n);
    const codeword = [F.__call__(3), F.__call__(1), F.__call__(4), F.__call__(1)];
    const auto = rs.fold(codeword, F.__call__(5));
    const explicit = rs.fold(codeword, F.__call__(5), [...rs.evaluation_points()]);
    expect(auto.map((e) => e.toString())).toEqual(explicit.map((e) => e.toString()));
  });
});

// Import ValueError for test assertions
import { ValueError } from '../errors.js';
import { Polynomial } from '../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
