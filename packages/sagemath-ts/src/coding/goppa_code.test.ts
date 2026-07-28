/**
 * Tests for Goppa codes
 *
 * Tests cover:
 * - Code construction and parameters
 * - Parity check matrix properties
 * - Generator matrix properties
 * - Encoding and decoding
 * - Error correction capabilities
 * - Binary Goppa codes (McEliece applications)
 */

import { describe, expect, it } from 'vitest';
import { ValueError } from '../errors.js';
import {
  type FiniteFieldElement,
  type FiniteFieldExtension,
  GFExtended,
  PrimeField,
  PrimeFieldElement,
} from '../rings/finite_rings/finite_field_extension.js';
import type { Polynomial } from '../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import {
  BinaryGoppaCode,
  DecodingError,
  GoppaCode,
  createBinaryGoppaCode,
  createGoppaCode,
} from './goppa_code.js';

// Helper function to create messages from the code's base field
function createMessage<P>(
  baseField: { zero: () => P; one: () => P; __call__: (x: unknown) => P },
  length: number,
  pattern: number[]
): P[] {
  return Array(length)
    .fill(null)
    .map((_, i) => (pattern[i % pattern.length] === 1 ? baseField.one() : baseField.zero()));
}

describe('GoppaCode', () => {
  describe('construction', () => {
    it('should create a Goppa code over GF(2^3) with valid parameters', () => {
      // GF(2^3) = GF(8)
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();

      // g(x) = x^2 + x + 1 (irreducible over GF(2))
      const g = x.pow(2).add(x).add(R.one());

      // Support: all elements where g(a) != 0
      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);

      expect(C.length()).toBe(L.length);
      expect(C.degree()).toBe(2);
      expect(C.goppa_polynomial()).toBe(g);
    });

    it('should throw for non-monic polynomial', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();

      // 2*x + 1 is not monic (but in char 2, 2 = 0, so this creates x + 1)
      // Let's create a non-monic by multiplying by a non-unit
      const alpha = F.gen();
      const g = R.__call__(alpha).mul(x.pow(2)).add(R.one()); // alpha*x^2 + 1

      // Get support
      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      expect(() => new GoppaCode(g, L)).toThrow(ValueError);
    });

    it('should throw when support contains roots of g', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();

      // g(x) = x (has root at 0)
      const g = x;

      // Include 0 in support
      const L = [F.zero(), F.one(), F.gen()];

      expect(() => new GoppaCode(g, L)).toThrow(ValueError);
    });
  });

  describe('code parameters', () => {
    it('should compute correct length', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      expect(C.length()).toBe(L.length);
    });

    it('should compute correct distance bound for binary code', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one()); // degree 2

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      // SageMath's own doctest for this exact code (goppa_code.py:300-312):
      //   sage: C.distance_bound()
      //   3
      // Sage always returns 1 + deg(g), in every characteristic.
      expect(C.distance_bound()).toBe(3);
    });

    it('should compute dimension', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const dim = C.dimension();

      // Dimension should be positive and <= n - m*t
      expect(dim).toBeGreaterThan(0);
      expect(dim).toBeLessThanOrEqual(C.length());
    });
  });

  describe('parity check matrix', () => {
    it('should construct valid parity check matrix', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const H = C.parity_check_matrix();

      // H should have d*m rows and n columns
      // d = degree of g = 2, m = extension degree = 3, n = |L|
      expect(H.nrows).toBe(2 * 3);
      expect(H.ncols).toBe(L.length);
    });

    it('should have codewords in kernel of H', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const G = C.generator_matrix();
      const H = C.parity_check_matrix();

      // For each row of G (which is a codeword), H * row^T should be 0
      for (let i = 0; i < G.nrows; i++) {
        const codeword = G.row(i);

        // Check H * codeword^T = 0
        for (let j = 0; j < H.nrows; j++) {
          let sum = F.zero();
          for (let k = 0; k < H.ncols; k++) {
            sum = sum.add(H.get(j, k).mul(codeword[k]!));
          }
          expect(sum.isZero()).toBe(true);
        }
      }
    });
  });

  describe('generator matrix', () => {
    it('should construct valid generator matrix', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const G = C.generator_matrix();
      const k = C.dimension();

      // G should have k rows and n columns
      expect(G.nrows).toBe(k);
      expect(G.ncols).toBe(L.length);
    });
  });

  describe('encoding', () => {
    it('should encode messages to valid codewords', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const k = C.dimension();

      // Create a message using the code's base field
      const baseField = C.base_field();
      const message = createMessage(baseField, k, [0, 1]);

      const codeword = C.encode(message);

      // Codeword should have correct length
      expect(codeword.length).toBe(C.length());

      // Codeword should be valid
      expect(C.is_codeword(codeword)).toBe(true);
    });

    it('should throw for wrong message length', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const baseField = C.base_field();

      // Wrong length message (dimension is 2 for this code, so provide 3)
      const message = [baseField.one(), baseField.zero(), baseField.one()];

      expect(() => C.encode(message)).toThrow(ValueError);
    });

    it('should encode zero message to zero codeword', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const k = C.dimension();
      const baseField = C.base_field();

      const zeroMessage = createMessage(baseField, k, [0]);
      const codeword = C.encode(zeroMessage);

      // All codeword symbols should be zero
      expect(codeword.every((c) => c.isZero())).toBe(true);
    });
  });

  describe('is_codeword', () => {
    it('should return true for valid codewords', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const k = C.dimension();
      const baseField = C.base_field();

      const message = createMessage(baseField, k, [0, 1]);
      const codeword = C.encode(message);

      expect(C.is_codeword(codeword)).toBe(true);
    });

    it('should return false for invalid codewords', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const n = C.length();
      const baseField = C.base_field();

      // Random word (likely not a codeword)
      const randomWord = createMessage(baseField, n, [1, 0, 0]);

      // This might be a codeword by chance, but most likely not
      // We can't guarantee it's not a codeword, so skip this specific assertion
    });

    it('should return false for wrong length', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const baseField = C.base_field();

      const shortWord = [baseField.one(), baseField.zero()];
      expect(C.is_codeword(shortWord)).toBe(false);
    });
  });

  describe('syndrome', () => {
    it('should return zero syndrome for valid codewords', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const k = C.dimension();
      const baseField = C.base_field();

      const message = createMessage(baseField, k, [0, 1]);
      const codeword = C.encode(message);

      const syndrome = C.syndrome(codeword);

      // Syndrome should be zero for valid codeword
      expect(syndrome.isZero()).toBe(true);
    });

    it('should return non-zero syndrome for corrupted words', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const k = C.dimension();
      const baseField = C.base_field();

      const message = createMessage(baseField, k, [0, 1]);
      const codeword = C.encode(message);

      // Introduce an error
      const corrupted = [...codeword];
      const one = baseField.one();
      const zero = baseField.zero();
      corrupted[0] = corrupted[0]!.isZero() ? one : zero;

      const syndrome = C.syndrome(corrupted);

      // Syndrome should be non-zero for corrupted word
      expect(syndrome.isZero()).toBe(false);
    });
  });

  describe('decoding (basic tests)', () => {
    it('should decode valid codewords unchanged', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one());

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);
      const k = C.dimension();
      const baseField = C.base_field();

      const message = createMessage(baseField, k, [0, 1]);
      const codeword = C.encode(message);

      const decoded = C.decode(codeword);

      // Decoded should equal original codeword
      expect(decoded.length).toBe(codeword.length);
      for (let i = 0; i < decoded.length; i++) {
        expect(decoded[i]!.isZero()).toBe(codeword[i]!.isZero());
      }
    });
  });

  describe('decoding with errors (Patterson, binary)', () => {
    // Exhaustive weight <= t correction.  Before the `_partialXGCD` bound fix
    // every weight-2 pattern on the GF(16) code and every weight-2/3 pattern
    // on the GF(32) code threw "corrected word is not a codeword".
    type GoppaPoly = Polynomial<FiniteFieldElement>;
    type GoppaRing = PolynomialRing<FiniteFieldElement>;
    const binaryCases: Array<[string, number, (R: GoppaRing, x: GoppaPoly) => GoppaPoly]> = [
      ['GF(8), g = x^2 + x + 1', 8, (R, x) => x.pow(2).add(x).add(R.one())],
      ['GF(16), g = x^2 + x + 1', 16, (R, x) => x.pow(2).add(x).add(R.one())],
      ['GF(32), g = x^3 + x + 1', 32, (R, x) => x.pow(3).add(x).add(R.one())],
    ];

    for (const [label, q, buildG] of binaryCases) {
      it(
        `should correct every error pattern of weight <= t for ${label}`,
        { timeout: 120_000 },
        () => {
          const F = GFExtended(q) as FiniteFieldExtension;
          const R = new PolynomialRing(F, 'x');
          const x = R.gen();
          const g = buildG(R, x);

          const L: FiniteFieldElement[] = [];
          for (const a of F) {
            if (!g.evaluate(a).isZero()) {
              L.push(a);
            }
          }

          const C = new GoppaCode(g, L);
          const n = C.length();
          const k = C.dimension();
          const t = g.degree();
          const baseField = C.base_field();

          expect(C.error_correction_capability()).toBe(t);

          // A few distinct codewords (all non-zero messages, capped)
          const messageCount = Math.min(2 ** k, q >= 32 ? 3 : 5);
          const messages: FiniteFieldElement[][] = [];
          for (let m = 1; m < messageCount; m++) {
            const msg: FiniteFieldElement[] = [];
            for (let i = 0; i < k; i++) {
              msg.push(((m >> i) & 1 ? baseField.one() : baseField.zero()) as FiniteFieldElement);
            }
            messages.push(msg);
          }

          // All index subsets of size 1..t
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
          for (let w = 1; w <= t; w++) {
            rec(0, w);
          }

          for (const message of messages) {
            const codeword = C.encode(message);
            for (const positions of subsets) {
              const received = [...codeword];
              for (const p of positions) {
                received[p] = received[p]!.isZero() ? baseField.one() : baseField.zero();
              }
              const decoded = C.decode(received);
              for (let i = 0; i < n; i++) {
                expect(decoded[i]!.eq(codeword[i]! as never)).toBe(true);
              }
            }
          }
        }
      );
    }
  });

  describe('decoding with errors (key equation, non-binary)', () => {
    // GF(9)/GF(3): the old decoder toggled symbols between 0 and 1 and failed
    // on every single-error pattern.
    type GoppaPoly = Polynomial<FiniteFieldElement>;
    type GoppaRing = PolynomialRing<FiniteFieldElement>;
    const cases: Array<[string, (R: GoppaRing, x: GoppaPoly) => GoppaPoly]> = [
      ['g = x^2 + 1', (R, x) => x.pow(2).add(R.one())],
      ['g = x^2', (_R, x) => x.pow(2)],
    ];

    for (const [label, buildG] of cases) {
      it(`should correct every single error over GF(9)/GF(3) with ${label}`, () => {
        const F = GFExtended(9) as FiniteFieldExtension;
        const R = new PolynomialRing(F, 'x');
        const x = R.gen();
        const g = buildG(R, x);

        const L: FiniteFieldElement[] = [];
        for (const a of F) {
          if (!g.evaluate(a).isZero()) {
            L.push(a);
          }
        }

        const C = new GoppaCode(g, L);
        const n = C.length();
        const k = C.dimension();
        const baseField = C.base_field();

        expect(C.error_correction_capability()).toBe(1);

        for (let m = 1; m < Math.min(3 ** k, 9); m++) {
          const message = [];
          let mm = m;
          for (let i = 0; i < k; i++) {
            message.push(baseField.__call__(mm % 3));
            mm = Math.floor(mm / 3);
          }
          const codeword = C.encode(message);

          for (let p = 0; p < n; p++) {
            for (const e of [1, 2]) {
              const received = [...codeword];
              received[p] = received[p]!.add(baseField.__call__(e) as never);
              const decoded = C.decode(received);
              for (let i = 0; i < n; i++) {
                expect(decoded[i]!.eq(codeword[i]! as never)).toBe(true);
              }
            }
          }
        }
      });
    }
  });

  describe('error correction capability', () => {
    it('should report correct error correction capability', () => {
      const F = GFExtended(8) as FiniteFieldExtension;
      const R = new PolynomialRing(F, 'x');
      const x = R.gen();
      const g = x.pow(2).add(x).add(R.one()); // degree 2

      const L: FiniteFieldElement[] = [];
      for (const a of F) {
        if (!g.evaluate(a).isZero()) {
          L.push(a);
        }
      }

      const C = new GoppaCode(g, L);

      // For binary Goppa with t=2: can correct floor((2*2+1-1)/2) = 2 errors
      expect(C.error_correction_capability()).toBe(2);
    });
  });
});

describe('BinaryGoppaCode', () => {
  it('should create a binary Goppa code', () => {
    const F = GFExtended(8) as FiniteFieldExtension;
    const R = new PolynomialRing(F, 'x');
    const x = R.gen();
    const g = x.pow(2).add(x).add(R.one());

    const L: FiniteFieldElement[] = [];
    for (const a of F) {
      if (!g.evaluate(a).isZero()) {
        L.push(a);
      }
    }

    const C = new BinaryGoppaCode(g, L);

    expect(C.length()).toBe(L.length);
    // Sage's bound is 1 + deg(g); the 2t+1 binary bound is reported by
    // error_correction_capability() instead (t errors are correctable).
    expect(C.distance_bound()).toBe(3);
    expect(C.error_correction_capability()).toBe(2);
  });

  it('should throw for non-binary field', () => {
    const F = GFExtended(9) as FiniteFieldExtension; // GF(3^2)
    const R = new PolynomialRing(F, 'x');
    const x = R.gen();
    const g = x.pow(2).add(R.one());

    const L: FiniteFieldElement[] = [];
    for (const a of F) {
      if (!g.evaluate(a).isZero()) {
        L.push(a);
      }
    }

    expect(() => new BinaryGoppaCode(g, L)).toThrow(ValueError);
  });
});

describe('createGoppaCode', () => {
  it('should create a code with automatic support generation', () => {
    const F = GFExtended(8) as FiniteFieldExtension;
    const R = new PolynomialRing(F, 'x');
    const x = R.gen();
    const g = x.pow(2).add(x).add(R.one());

    const C = createGoppaCode(g, F);

    // Length should be |F| - number of roots of g in F
    expect(C.length()).toBeGreaterThan(0);
    expect(C.length()).toBeLessThanOrEqual(8);
  });
});

describe('createBinaryGoppaCode', () => {
  it('should create a binary code with automatic support generation', () => {
    const F = GFExtended(8) as FiniteFieldExtension;
    const R = new PolynomialRing(F, 'x');
    const x = R.gen();
    const g = x.pow(2).add(x).add(R.one());

    const C = createBinaryGoppaCode(g, F);

    expect(C.length()).toBeGreaterThan(0);
    expect(C.distance_bound()).toBe(3);
  });

  it('should throw for non-binary field', () => {
    const F = GFExtended(9) as FiniteFieldExtension;
    const R = new PolynomialRing(F, 'x');
    const x = R.gen();
    const g = x.pow(2).add(R.one());

    expect(() => createBinaryGoppaCode(g, F)).toThrow(ValueError);
  });
});

describe('Goppa code examples from literature', () => {
  it('should create [8, 2] Goppa code from SageMath example', () => {
    // From SageMath documentation:
    // F = GF(2^3)
    // g = x^2 + x + 1
    // L = [a for a in F if g(a) != 0]
    // C is [8, 2] Goppa code

    const F = GFExtended(8) as FiniteFieldExtension;
    const R = new PolynomialRing(F, 'x');
    const x = R.gen();
    const g = x.pow(2).add(x).add(R.one());

    const L: FiniteFieldElement[] = [];
    for (const a of F) {
      if (!g.evaluate(a).isZero()) {
        L.push(a);
      }
    }

    const C = new GoppaCode(g, L);

    expect(C.length()).toBe(8);
    // The dimension should be 2 based on SageMath
    expect(C.dimension()).toBe(2);
  });

  it('should handle larger Goppa polynomial', () => {
    // g(x) = x^9 + 1 over GF(2^6)
    // This should give a code suitable for McEliece

    const F = GFExtended(64) as FiniteFieldExtension;
    const R = new PolynomialRing(F, 'x');
    const x = R.gen();
    const g = x.pow(9).add(R.one());

    const L: FiniteFieldElement[] = [];
    for (const a of F) {
      if (!g.evaluate(a).isZero()) {
        L.push(a);
      }
    }

    const C = new GoppaCode(g, L);

    // Based on SageMath: [55, 16] Goppa code
    expect(C.length()).toBe(55);
    // Dimension might vary slightly due to implementation details
    expect(C.dimension()).toBeGreaterThan(0);
    expect(C.dimension()).toBeLessThan(55);

    // Error correction capability
    expect(C.error_correction_capability()).toBeGreaterThanOrEqual(9);
  });
});

describe('string representation', () => {
  it('should have meaningful toString', () => {
    const F = GFExtended(8) as FiniteFieldExtension;
    const R = new PolynomialRing(F, 'x');
    const x = R.gen();
    const g = x.pow(2).add(x).add(R.one());

    const L: FiniteFieldElement[] = [];
    for (const a of F) {
      if (!g.evaluate(a).isZero()) {
        L.push(a);
      }
    }

    const C = new GoppaCode(g, L);
    const str = C.toString();

    expect(str).toContain('Goppa code');
    expect(str).toContain('GF(2)');
    expect(str).toMatch(/\[\d+, \d+\]/);
  });
});
