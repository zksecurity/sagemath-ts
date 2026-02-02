/**
 * Tests for BCH (Bose-Chaudhuri-Hocquenghem) codes
 */

import { describe, expect, it } from 'vitest';
import { ValueError } from '../errors.js';
import {
  FiniteFieldElement,
  FiniteFieldExtension,
  PrimeField,
  type PrimeFieldElement,
} from '../rings/finite_rings/finite_field_extension.js';
import {
  BCHCode,
  DecodingError,
  createBCHCode,
  createPrimitiveBCHCode,
  isReedSolomonCode,
} from './bch_code.js';

describe('BCHCode', () => {
  describe('constructor', () => {
    it('should create a valid BCH(15, 7) code over GF(2)', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 7n, F2);

      expect(bch.length()).toBe(15);
      expect(bch.designed_distance()).toBe(7);
    });

    it('should create a valid BCH(15, 5) code over GF(2)', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);

      expect(bch.length()).toBe(15);
      expect(bch.designed_distance()).toBe(5);
    });

    it('should create a valid BCH(31, 7) code over GF(2)', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(31n, 7n, F2);

      expect(bch.length()).toBe(31);
      expect(bch.designed_distance()).toBe(7);
    });

    it('should throw for delta < 1', () => {
      const F2 = new PrimeField(2);
      expect(() => new BCHCode(15n, 0n, F2)).toThrow(ValueError);
    });

    it('should throw for delta > n', () => {
      const F2 = new PrimeField(2);
      expect(() => new BCHCode(15n, 16n, F2)).toThrow(ValueError);
    });

    it('should accept custom offset (non-narrow-sense)', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2, 0n); // offset = 0

      expect(bch.offset).toBe(0);
      expect(bch.length()).toBe(15);
    });
  });

  describe('defining_set', () => {
    it('should compute correct defining set for BCH(15, 5)', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);

      const D = bch.defining_set();

      // For delta=5, b=1, l=1: D = {1, 2, 3, 4}
      expect(D.length).toBe(4);
      expect(D).toContain(1);
      expect(D).toContain(2);
      expect(D).toContain(3);
      expect(D).toContain(4);
    });

    it('should compute correct defining set for BCH(15, 7)', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 7n, F2);

      const D = bch.defining_set();

      // For delta=7, b=1, l=1: D = {1, 2, 3, 4, 5, 6}
      expect(D.length).toBe(6);
      for (let i = 1; i <= 6; i++) {
        expect(D).toContain(i);
      }
    });

    it('should handle delta=1 (trivial code)', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 1n, F2);

      const D = bch.defining_set();

      // For delta=1, defining set is empty
      expect(D.length).toBe(0);
    });
  });

  describe('generator_polynomial', () => {
    it('should compute generator polynomial for BCH(15, 5)', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);

      const g = bch.generator_polynomial();

      // BCH(15, 5) over GF(2) has generator of degree 8
      // g(x) = LCM(m_1(x), m_3(x)) where m_i are minimal polynomials
      // m_1(x) = x^4 + x + 1 (primitive polynomial)
      // m_3(x) = x^4 + x^3 + x^2 + x + 1
      // g(x) has degree 8
      expect(g.degree()).toBe(8);
    });

    it('should compute generator polynomial for BCH(15, 7)', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 7n, F2);

      const g = bch.generator_polynomial();

      // BCH(15, 7) over GF(2) - classic example
      // g(x) = x^10 + x^8 + x^5 + x^4 + x^2 + x + 1
      // degree should be 10
      expect(g.degree()).toBe(10);
    });

    it('should return 1 for delta=1', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 1n, F2);

      const g = bch.generator_polynomial();

      expect(g.degree()).toBe(0);
      expect(g.leading_coefficient().eq(1)).toBe(true);
    });

    it('should divide x^n - 1', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);

      const g = bch.generator_polynomial();
      const R = g.parent;

      // x^15 - 1 (in GF(2), this is x^15 + 1)
      const xn = R.monomial(15);
      const one = R.one();
      const xnMinusOne = xn.sub(one);

      // Check that g divides x^n - 1
      const [_q, remainder] = xnMinusOne.quo_rem(g);
      expect(remainder.isZero()).toBe(true);
    });
  });

  describe('dimension', () => {
    it('should compute correct dimension for BCH(15, 5)', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);

      // k = n - deg(g) = 15 - 8 = 7
      expect(bch.dimension()).toBe(7);
    });

    it('should compute correct dimension for BCH(15, 7)', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 7n, F2);

      // k = n - deg(g) = 15 - 10 = 5
      expect(bch.dimension()).toBe(5);
    });

    it('should return n for delta=1', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 1n, F2);

      // Trivial code: k = n
      expect(bch.dimension()).toBe(15);
    });
  });

  describe('minimum_distance', () => {
    it('should return at least the designed distance', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 7n, F2);

      expect(bch.minimum_distance()).toBeGreaterThanOrEqual(7);
    });
  });

  describe('decoding_radius', () => {
    it('should compute correct decoding radius', () => {
      const F2 = new PrimeField(2);

      // BCH(15, 7): t = floor((7-1)/2) = 3
      const bch7 = new BCHCode(15n, 7n, F2);
      expect(bch7.decoding_radius()).toBe(3);

      // BCH(15, 5): t = floor((5-1)/2) = 2
      const bch5 = new BCHCode(15n, 5n, F2);
      expect(bch5.decoding_radius()).toBe(2);

      // BCH(15, 3): t = floor((3-1)/2) = 1
      const bch3 = new BCHCode(15n, 3n, F2);
      expect(bch3.decoding_radius()).toBe(1);
    });
  });

  describe('encode', () => {
    it('should encode a message to correct length', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__(i % 2));
      }

      const codeword = bch.encode(message);

      expect(codeword.length).toBe(15);
    });

    it('should produce valid codewords', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__(i % 2));
      }

      const codeword = bch.encode(message);

      // Verify it's a valid codeword (zero syndrome)
      expect(bch.is_codeword(codeword)).toBe(true);
    });

    it('should throw for wrong message length', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);

      const shortMessage = [F2.__call__(1), F2.__call__(0), F2.__call__(1)];

      expect(() => bch.encode(shortMessage)).toThrow(ValueError);
    });

    it('should encode zero message to zero codeword', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.zero());
      }

      const codeword = bch.encode(message);

      expect(codeword.every((c) => c.isZero())).toBe(true);
    });
  });

  describe('syndrome', () => {
    it('should compute zero syndromes for valid codewords', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__(i % 2));
      }

      const codeword = bch.encode(message);
      const syndromes = bch.syndrome(codeword);

      // All syndromes should be zero for a valid codeword
      expect(syndromes.every((s) => s.isZero())).toBe(true);
    });

    it('should compute non-zero syndromes for corrupted words', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__(i % 2));
      }

      const codeword = bch.encode(message);

      // Introduce an error
      const corrupted = [...codeword];
      corrupted[0] = corrupted[0]!.add(F2.one());

      const syndromes = bch.syndrome(corrupted);

      // At least one syndrome should be non-zero
      expect(syndromes.some((s) => !s.isZero())).toBe(true);
    });

    it('should return delta-1 syndromes', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 7n, F2);
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__(1));
      }

      const codeword = bch.encode(message);
      const syndromes = bch.syndrome(codeword);

      expect(syndromes.length).toBe(6); // delta - 1 = 7 - 1 = 6
    });
  });

  describe('decode', () => {
    it('should decode a valid codeword correctly', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__(i % 2));
      }

      const codeword = bch.encode(message);
      const decoded = bch.decode(codeword);

      // Decoded should equal original codeword
      for (let i = 0; i < 15; i++) {
        expect(decoded[i]!.eq(codeword[i]!)).toBe(true);
      }
    });

    it('should correct a single error', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2); // can correct 2 errors
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__(i % 2));
      }

      const codeword = bch.encode(message);

      // Introduce one error
      const corrupted = [...codeword];
      corrupted[3] = corrupted[3]!.add(F2.one());

      const decoded = bch.decode(corrupted);

      // Decoded should equal original codeword
      for (let i = 0; i < 15; i++) {
        expect(decoded[i]!.eq(codeword[i]!)).toBe(true);
      }
    });

    it('should correct up to floor((delta-1)/2) errors', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2); // t = 2
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__(1));
      }

      const codeword = bch.encode(message);

      // Introduce 2 errors (maximum correctable)
      const corrupted = [...codeword];
      corrupted[0] = corrupted[0]!.add(F2.one());
      corrupted[7] = corrupted[7]!.add(F2.one());

      const decoded = bch.decode(corrupted);

      // Decoded should equal original codeword
      for (let i = 0; i < 15; i++) {
        expect(decoded[i]!.eq(codeword[i]!)).toBe(true);
      }
    });
  });

  describe('encode/decode round-trip', () => {
    it('should successfully round-trip with no errors', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__((i * i + 1) % 2));
      }

      const codeword = bch.encode(message);
      const decodedCodeword = bch.decode(codeword);
      const decodedMessage = bch.decode_to_message(codeword);

      // Codeword should be unchanged
      for (let i = 0; i < 15; i++) {
        expect(decodedCodeword[i]!.eq(codeword[i]!)).toBe(true);
      }

      // Message should be recovered
      for (let i = 0; i < k; i++) {
        expect(decodedMessage[i]!.eq(message[i]!)).toBe(true);
      }
    });

    it('should successfully round-trip with correctable errors', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 7n, F2); // t = 3
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__(i % 2));
      }

      const codeword = bch.encode(message);

      // Introduce 3 errors
      const corrupted = [...codeword];
      corrupted[1] = corrupted[1]!.add(F2.one());
      corrupted[5] = corrupted[5]!.add(F2.one());
      corrupted[12] = corrupted[12]!.add(F2.one());

      const decoded = bch.decode(corrupted);

      // Should recover original codeword
      for (let i = 0; i < 15; i++) {
        expect(decoded[i]!.eq(codeword[i]!)).toBe(true);
      }
    });
  });

  describe('is_codeword', () => {
    it('should return true for valid codewords', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__(i % 2));
      }

      const codeword = bch.encode(message);

      expect(bch.is_codeword(codeword)).toBe(true);
    });

    it('should return false for corrupted words', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__(i % 2));
      }

      const codeword = bch.encode(message);

      // Introduce an error
      const corrupted = [...codeword];
      corrupted[5] = corrupted[5]!.add(F2.one());

      expect(bch.is_codeword(corrupted)).toBe(false);
    });

    it('should return false for wrong length', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 5n, F2);

      const shortWord: PrimeFieldElement[] = [];
      for (let i = 0; i < 10; i++) {
        shortWord.push(F2.zero());
      }

      expect(bch.is_codeword(shortWord)).toBe(false);
    });
  });

  describe('BCH(31, 7) over GF(2)', () => {
    it('should have correct parameters', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(31n, 7n, F2);

      expect(bch.length()).toBe(31);
      expect(bch.designed_distance()).toBe(7);
      expect(bch.decoding_radius()).toBe(3);
    });

    it('should produce valid codewords', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(31n, 7n, F2);
      const k = bch.dimension();

      const message: PrimeFieldElement[] = [];
      for (let i = 0; i < k; i++) {
        message.push(F2.__call__(i % 2));
      }

      const codeword = bch.encode(message);

      expect(codeword.length).toBe(31);
      expect(bch.is_codeword(codeword)).toBe(true);
    });
  });

  describe('toString', () => {
    it('should return descriptive string', () => {
      const F2 = new PrimeField(2);
      const bch = new BCHCode(15n, 7n, F2);

      const str = bch.toString();

      expect(str).toContain('15');
      expect(str).toContain('5'); // dimension
      expect(str).toContain('BCH');
      expect(str).toContain('7'); // designed distance
    });
  });
});

describe('createBCHCode', () => {
  it('should create a narrow-sense BCH code by default', () => {
    const bch = createBCHCode(15n, 5n, 2n);

    expect(bch.offset).toBe(1);
    expect(bch.length()).toBe(15);
    expect(bch.designed_distance()).toBe(5);
  });

  it('should create a non-narrow-sense BCH code when specified', () => {
    const bch = createBCHCode(15n, 5n, 2n, false);

    expect(bch.offset).toBe(0);
  });
});

describe('createPrimitiveBCHCode', () => {
  it('should create a primitive BCH code with correct length', () => {
    // GF(2), m=4: n = 2^4 - 1 = 15
    const bch = createPrimitiveBCHCode(4n, 5n, 2n);

    expect(bch.length()).toBe(15);
    expect(bch.designed_distance()).toBe(5);
  });

  it('should create a primitive BCH code for m=5', () => {
    // GF(2), m=5: n = 2^5 - 1 = 31
    const bch = createPrimitiveBCHCode(5n, 7n, 2n);

    expect(bch.length()).toBe(31);
    expect(bch.designed_distance()).toBe(7);
  });
});

describe('isReedSolomonCode', () => {
  it('should identify Reed-Solomon codes', () => {
    // Reed-Solomon: n = q - 1 with narrow-sense parameters
    // GF(16), n = 15
    const F16 = new FiniteFieldExtension(2, 4);
    const bch = new BCHCode(15n, 5n, F16, 1n, 1n);

    expect(isReedSolomonCode(bch)).toBe(true);
  });

  it('should not identify non-RS BCH codes', () => {
    // BCH over GF(2) with n = 15 is not RS (would need GF(16))
    const F2 = new PrimeField(2);
    const bch = new BCHCode(15n, 5n, F2);

    expect(isReedSolomonCode(bch)).toBe(false);
  });
});
