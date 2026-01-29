/**
 * Unit tests for FFT/NTT convolution module
 */
import { describe, expect, test } from 'bun:test';
import {
  FiniteFieldPrime,
  FiniteFieldElement,
} from '../finite_rings/finite_field_prime.js';
import { PolynomialRingConstructor } from './polynomial_ring.js';
import { Polynomial } from './polynomial_element.js';
import {
  ntt,
  intt,
  ntt_multiply,
  fft,
  ifft,
  fft_multiply,
  evaluate_on_domain,
  interpolate_from_domain,
  find_primitive_root,
  max_fft_size,
  dft_naive,
  idft_naive,
  convolve_naive,
  supports_ntt_size,
  find_ntt_prime,
  NTT_FRIENDLY_PRIMES,
} from './convolution.js';

/**
 * Create a finite field GF(p) for testing.
 */
function GF(p: bigint): FiniteFieldPrime {
  return new FiniteFieldPrime(p);
}

// ============================================================================
// Test find_primitive_root
// ============================================================================

describe('find_primitive_root', () => {
  test('finds 4th root of unity in GF(17)', () => {
    const F = GF(17n);
    const omega = find_primitive_root(4, F);

    // omega^4 should be 1
    expect(omega.pow(4).isOne()).toBe(true);

    // omega^2 should not be 1 (primitive check)
    expect(omega.pow(2).isOne()).toBe(false);

    // omega^1 should not be 1
    expect(omega.pow(1).isOne()).toBe(false);
  });

  test('finds 8th root of unity in GF(17)', () => {
    const F = GF(17n);
    const omega = find_primitive_root(8, F);

    // omega^8 should be 1
    expect(omega.pow(8).isOne()).toBe(true);

    // omega^4 should not be 1
    expect(omega.pow(4).isOne()).toBe(false);
  });

  test('finds 16th root of unity in GF(17)', () => {
    const F = GF(17n);
    // 17 - 1 = 16, so 16th roots exist
    const omega = find_primitive_root(16, F);

    expect(omega.pow(16).isOne()).toBe(true);
    expect(omega.pow(8).isOne()).toBe(false);
    expect(omega.pow(4).isOne()).toBe(false);
    expect(omega.pow(2).isOne()).toBe(false);
    expect(omega.pow(1).isOne()).toBe(false);
  });

  test('finds 2nd root of unity in GF(17)', () => {
    const F = GF(17n);
    const omega = find_primitive_root(2, F);

    // omega^2 should be 1
    expect(omega.pow(2).isOne()).toBe(true);

    // omega should be -1 (i.e., p - 1 = 16)
    expect(omega.value).toBe(16n);
  });

  test('throws for non-divisor of p-1', () => {
    const F = GF(17n);
    // 17 - 1 = 16, and 3 does not divide 16
    expect(() => find_primitive_root(3, F)).toThrow();
  });

  test('works for larger field GF(97)', () => {
    const F = GF(97n);
    // 97 - 1 = 96 = 2^5 * 3
    // So 32 divides 96

    const omega = find_primitive_root(32, F);
    expect(omega.pow(32).isOne()).toBe(true);
    expect(omega.pow(16).isOne()).toBe(false);
  });
});

// ============================================================================
// Test max_fft_size
// ============================================================================

describe('max_fft_size', () => {
  test('GF(17) supports up to 2^4 = 16', () => {
    const F = GF(17n);
    expect(max_fft_size(F)).toBe(4);
  });

  test('GF(97) supports up to 2^5 = 32', () => {
    const F = GF(97n);
    // 97 - 1 = 96 = 2^5 * 3
    expect(max_fft_size(F)).toBe(5);
  });

  test('GF(257) supports up to 2^8 = 256', () => {
    const F = GF(257n);
    // 257 - 1 = 256 = 2^8
    expect(max_fft_size(F)).toBe(8);
  });

  test('GF(3) supports up to 2^1 = 2', () => {
    const F = GF(3n);
    // 3 - 1 = 2 = 2^1
    expect(max_fft_size(F)).toBe(1);
  });
});

// ============================================================================
// Test NTT vs naive DFT
// ============================================================================

describe('NTT vs naive DFT', () => {
  test('NTT matches DFT for size 4 in GF(17)', () => {
    const F = GF(17n);
    const omega = find_primitive_root(4, F);
    const coeffs = [F.__call__(1n), F.__call__(2n), F.__call__(3n), F.__call__(4n)];

    const nttResult = ntt(coeffs, omega, F);
    const dftResult = dft_naive(coeffs, omega, F);

    for (let i = 0; i < 4; i++) {
      expect(nttResult[i]!.eq(dftResult[i]!)).toBe(true);
    }
  });

  test('NTT matches DFT for size 8 in GF(17)', () => {
    const F = GF(17n);
    const omega = find_primitive_root(8, F);
    const coeffs = [1n, 2n, 3n, 4n, 5n, 6n, 7n, 8n].map((v) => F.__call__(v));

    const nttResult = ntt(coeffs, omega, F);
    const dftResult = dft_naive(coeffs, omega, F);

    for (let i = 0; i < 8; i++) {
      expect(nttResult[i]!.eq(dftResult[i]!)).toBe(true);
    }
  });

  test('NTT matches DFT for size 16 in GF(17)', () => {
    const F = GF(17n);
    const omega = find_primitive_root(16, F);
    const coeffs: FiniteFieldElement[] = [];
    for (let i = 0; i < 16; i++) {
      coeffs.push(F.__call__(BigInt(i + 1)));
    }

    const nttResult = ntt(coeffs, omega, F);
    const dftResult = dft_naive(coeffs, omega, F);

    for (let i = 0; i < 16; i++) {
      expect(nttResult[i]!.eq(dftResult[i]!)).toBe(true);
    }
  });

  test('NTT matches DFT with random coefficients', () => {
    const F = GF(97n);
    const omega = find_primitive_root(8, F);

    // Random-ish coefficients
    const coeffs = [23n, 45n, 67n, 12n, 89n, 34n, 56n, 78n].map((v) => F.__call__(v));

    const nttResult = ntt(coeffs, omega, F);
    const dftResult = dft_naive(coeffs, omega, F);

    for (let i = 0; i < 8; i++) {
      expect(nttResult[i]!.eq(dftResult[i]!)).toBe(true);
    }
  });
});

// ============================================================================
// Test INTT (inverse NTT)
// ============================================================================

describe('INTT (inverse NTT)', () => {
  test('INTT(NTT(x)) = x for size 4', () => {
    const F = GF(17n);
    const omega = find_primitive_root(4, F);
    const original = [F.__call__(1n), F.__call__(2n), F.__call__(3n), F.__call__(4n)];

    const transformed = ntt(original, omega, F);
    const recovered = intt(transformed, omega, F);

    for (let i = 0; i < 4; i++) {
      expect(recovered[i]!.eq(original[i]!)).toBe(true);
    }
  });

  test('INTT(NTT(x)) = x for size 8', () => {
    const F = GF(17n);
    const omega = find_primitive_root(8, F);
    const original = [1n, 2n, 3n, 4n, 5n, 6n, 7n, 8n].map((v) => F.__call__(v));

    const transformed = ntt(original, omega, F);
    const recovered = intt(transformed, omega, F);

    for (let i = 0; i < 8; i++) {
      expect(recovered[i]!.eq(original[i]!)).toBe(true);
    }
  });

  test('INTT(NTT(x)) = x for size 16', () => {
    const F = GF(17n);
    const omega = find_primitive_root(16, F);
    const original: FiniteFieldElement[] = [];
    for (let i = 0; i < 16; i++) {
      original.push(F.__call__(BigInt((i * 7 + 3) % 17)));
    }

    const transformed = ntt(original, omega, F);
    const recovered = intt(transformed, omega, F);

    for (let i = 0; i < 16; i++) {
      expect(recovered[i]!.eq(original[i]!)).toBe(true);
    }
  });

  test('INTT matches naive IDFT', () => {
    const F = GF(17n);
    const omega = find_primitive_root(8, F);
    const values = [10n, 5n, 8n, 3n, 12n, 7n, 2n, 15n].map((v) => F.__call__(v));

    const inttResult = intt(values, omega, F);
    const idftResult = idft_naive(values, omega, F);

    for (let i = 0; i < 8; i++) {
      expect(inttResult[i]!.eq(idftResult[i]!)).toBe(true);
    }
  });
});

// ============================================================================
// Test NTT-based polynomial multiplication
// ============================================================================

describe('NTT-based polynomial multiplication', () => {
  test('multiply simple polynomials in GF(17)', () => {
    const F = GF(17n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    // f = x + 1
    const f = x.add(R.one());

    // g = x + 2
    const g = x.add(R.__call__(F.__call__(2n)));

    // Expected: (x + 1)(x + 2) = x^2 + 3x + 2
    const expected = x.mul(x).add(x.scalar_mul(F.__call__(3n))).add(R.__call__(F.__call__(2n)));

    const product = ntt_multiply(f, g, F);

    expect(product.eq(expected)).toBe(true);
  });

  test('multiply quadratic polynomials', () => {
    const F = GF(17n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    // f = x^2 + x + 1
    const f = x.mul(x).add(x).add(R.one());

    // g = x^2 - 1 = x^2 + 16 (mod 17)
    const g = x.mul(x).add(R.__call__(F.__call__(16n)));

    // Expected: (x^2 + x + 1)(x^2 + 16)
    // = x^4 + 16x^2 + x^3 + 16x + x^2 + 16
    // = x^4 + x^3 + (16 + 1)x^2 + 16x + 16
    // = x^4 + x^3 + 0*x^2 + 16x + 16
    const expected = f.mul(g);
    const product = ntt_multiply(f, g, F);

    expect(product.eq(expected)).toBe(true);
  });

  test('multiply by constant', () => {
    const F = GF(17n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    // f = x^2 + 3x + 5
    const f = x.mul(x).add(x.scalar_mul(F.__call__(3n))).add(R.__call__(F.__call__(5n)));

    // g = 4
    const g = R.__call__(F.__call__(4n));

    const expected = f.mul(g);
    const product = ntt_multiply(f, g, F);

    expect(product.eq(expected)).toBe(true);
  });

  test('multiply by zero gives zero', () => {
    const F = GF(17n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    const f = x.mul(x).add(x);
    const g = R.zero();

    const product = ntt_multiply(f, g, F);

    expect(product.isZero()).toBe(true);
  });

  test('NTT multiplication matches schoolbook for larger polynomials', () => {
    const F = GF(97n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    // f = 5x^4 + 3x^3 + 2x^2 + 7x + 11
    let f = R.zero();
    f = f.add(x.pow(4).scalar_mul(F.__call__(5n)));
    f = f.add(x.pow(3).scalar_mul(F.__call__(3n)));
    f = f.add(x.pow(2).scalar_mul(F.__call__(2n)));
    f = f.add(x.scalar_mul(F.__call__(7n)));
    f = f.add(R.__call__(F.__call__(11n)));

    // g = 6x^3 + 4x^2 + 8x + 9
    let g = R.zero();
    g = g.add(x.pow(3).scalar_mul(F.__call__(6n)));
    g = g.add(x.pow(2).scalar_mul(F.__call__(4n)));
    g = g.add(x.scalar_mul(F.__call__(8n)));
    g = g.add(R.__call__(F.__call__(9n)));

    const expected = f.mul(g); // Schoolbook multiplication
    const product = ntt_multiply(f, g, F);

    expect(product.eq(expected)).toBe(true);
  });
});

// ============================================================================
// Test naive convolution (for verification)
// ============================================================================

describe('Naive convolution', () => {
  test('convolve_naive computes correct product', () => {
    const F = GF(17n);

    // f = [1, 2, 3] represents 1 + 2x + 3x^2
    // g = [4, 5] represents 4 + 5x
    // Product: (1 + 2x + 3x^2)(4 + 5x)
    //        = 4 + 5x + 8x + 10x^2 + 12x^2 + 15x^3
    //        = 4 + 13x + 22x^2 + 15x^3
    //        = 4 + 13x + 5x^2 + 15x^3 (mod 17)

    const fCoeffs = [1n, 2n, 3n].map((v) => F.__call__(v));
    const gCoeffs = [4n, 5n].map((v) => F.__call__(v));

    const result = convolve_naive(fCoeffs, gCoeffs, F);

    expect(result[0]!.eq(F.__call__(4n))).toBe(true);
    expect(result[1]!.eq(F.__call__(13n))).toBe(true);
    expect(result[2]!.eq(F.__call__(5n))).toBe(true); // 22 mod 17 = 5
    expect(result[3]!.eq(F.__call__(15n))).toBe(true);
  });

  test('NTT multiplication matches naive convolution', () => {
    const F = GF(17n);

    const fCoeffs = [1n, 2n, 3n, 4n].map((v) => F.__call__(v));
    const gCoeffs = [5n, 6n, 7n, 8n].map((v) => F.__call__(v));

    // Naive result
    const naiveResult = convolve_naive(fCoeffs, gCoeffs, F);

    // NTT result (need to pad to power of 2)
    const n = 8; // Next power of 2 >= 4 + 4 - 1 = 7
    const omega = find_primitive_root(n, F);

    const fPadded: FiniteFieldElement[] = [];
    const gPadded: FiniteFieldElement[] = [];
    for (let i = 0; i < n; i++) {
      fPadded.push(i < fCoeffs.length ? fCoeffs[i]! : F.zero());
      gPadded.push(i < gCoeffs.length ? gCoeffs[i]! : F.zero());
    }

    const fNTT = ntt(fPadded, omega, F);
    const gNTT = ntt(gPadded, omega, F);

    const productNTT: FiniteFieldElement[] = [];
    for (let i = 0; i < n; i++) {
      productNTT.push(fNTT[i]!.mul(gNTT[i]!));
    }

    const nttResult = intt(productNTT, omega, F);

    // Compare (only first naiveResult.length coefficients matter)
    for (let i = 0; i < naiveResult.length; i++) {
      expect(nttResult[i]!.eq(naiveResult[i]!)).toBe(true);
    }
  });
});

// ============================================================================
// Test evaluate_on_domain and interpolate_from_domain
// ============================================================================

describe('Domain evaluation and interpolation', () => {
  test('evaluate on domain matches polynomial evaluation', () => {
    const F = GF(17n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    // f = x^2 + 3x + 2
    const f = x.mul(x).add(x.scalar_mul(F.__call__(3n))).add(R.__call__(F.__call__(2n)));

    const omega = find_primitive_root(4, F);
    const evals = evaluate_on_domain(f, omega, F, 4);

    // Verify each evaluation
    for (let i = 0; i < 4; i++) {
      const point = omega.pow(i);
      const expected = f.evaluate(point);
      expect(evals[i]!.eq(expected)).toBe(true);
    }
  });

  test('interpolate recovers original polynomial', () => {
    const F = GF(17n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    // f = x^3 + 2x^2 + 5x + 7
    let f = R.zero();
    f = f.add(x.pow(3));
    f = f.add(x.pow(2).scalar_mul(F.__call__(2n)));
    f = f.add(x.scalar_mul(F.__call__(5n)));
    f = f.add(R.__call__(F.__call__(7n)));

    const omega = find_primitive_root(4, F);
    const evals = evaluate_on_domain(f, omega, F, 4);
    const recovered = interpolate_from_domain(evals, omega, F, R);

    expect(recovered.eq(f)).toBe(true);
  });

  test('interpolate with larger domain', () => {
    const F = GF(17n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    // f = x^2 + 1
    const f = x.mul(x).add(R.one());

    // Use domain size 8 (larger than needed)
    const omega = find_primitive_root(8, F);
    const evals = evaluate_on_domain(f, omega, F, 8);
    const recovered = interpolate_from_domain(evals, omega, F, R);

    // The recovered polynomial should equal f (high-degree coeffs are 0)
    for (let i = 0; i <= f.degree(); i++) {
      expect(recovered.getCoeff(i).eq(f.getCoeff(i))).toBe(true);
    }
    // Higher degree coefficients should be zero
    for (let i = f.degree() + 1; i < 8; i++) {
      expect(recovered.getCoeff(i).isZero()).toBe(true);
    }
  });

  test('roundtrip with random polynomial', () => {
    const F = GF(97n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    // f = 23x^5 + 45x^4 + 12x^3 + 67x^2 + 89x + 34
    let f = R.zero();
    const coefficients = [34n, 89n, 67n, 12n, 45n, 23n];
    for (let i = 0; i < coefficients.length; i++) {
      f = f.add(x.pow(i).scalar_mul(F.__call__(coefficients[i]!)));
    }

    const omega = find_primitive_root(8, F);
    const evals = evaluate_on_domain(f, omega, F, 8);
    const recovered = interpolate_from_domain(evals, omega, F, R);

    // Check each coefficient
    for (let i = 0; i <= f.degree(); i++) {
      expect(recovered.getCoeff(i).eq(f.getCoeff(i))).toBe(true);
    }
  });
});

// ============================================================================
// Test fft_multiply (the convenience wrapper)
// ============================================================================

describe('fft_multiply convenience wrapper', () => {
  test('multiplies polynomials over finite field', () => {
    const F = GF(17n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    const f = x.add(R.one()); // x + 1
    const g = x.sub(R.one()); // x - 1 = x + 16

    // (x + 1)(x - 1) = x^2 - 1 = x^2 + 16
    const expected = f.mul(g);
    const product = fft_multiply(f, g);

    expect(product.eq(expected)).toBe(true);
  });
});

// ============================================================================
// Test NTT-friendly prime utilities
// ============================================================================

describe('NTT-friendly prime utilities', () => {
  test('supports_ntt_size correctly identifies support', () => {
    // GF(17): 17 - 1 = 16 = 2^4
    expect(supports_ntt_size(17n, 1)).toBe(true);
    expect(supports_ntt_size(17n, 2)).toBe(true);
    expect(supports_ntt_size(17n, 4)).toBe(true);
    expect(supports_ntt_size(17n, 8)).toBe(true);
    expect(supports_ntt_size(17n, 16)).toBe(true);
    expect(supports_ntt_size(17n, 32)).toBe(false);

    // GF(97): 97 - 1 = 96 = 2^5 * 3
    expect(supports_ntt_size(97n, 32)).toBe(true);
    expect(supports_ntt_size(97n, 64)).toBe(false);
  });

  test('find_ntt_prime finds valid primes', () => {
    const p = find_ntt_prime(256, 16);

    // Check that p is of the form k * 256 + 1
    expect((p - 1n) % 256n).toBe(0n);

    // Check that p has at least 15 bits (16 - 1 due to our search)
    expect(p >= (1n << 15n)).toBe(true);
  });

  test('well-known NTT primes are valid', () => {
    // P_30 = 3 * 2^30 + 1
    expect((NTT_FRIENDLY_PRIMES.P_30 - 1n) % (1n << 30n)).toBe(0n);

    // P_25 = 5 * 2^25 + 1
    expect((NTT_FRIENDLY_PRIMES.P_25 - 1n) % (1n << 25n)).toBe(0n);

    // P_SMALL = 17, supports up to 16
    expect((NTT_FRIENDLY_PRIMES.P_SMALL - 1n) % 16n).toBe(0n);
  });
});

// ============================================================================
// Edge cases
// ============================================================================

describe('Edge cases', () => {
  test('NTT of empty array', () => {
    const F = GF(17n);
    const omega = find_primitive_root(4, F);
    const result = ntt([], omega, F);
    expect(result.length).toBe(0);
  });

  test('NTT of single element', () => {
    const F = GF(17n);
    const omega = find_primitive_root(1, F); // omega = 1 for n = 1
    const coeffs = [F.__call__(5n)];
    const result = ntt(coeffs, omega, F);
    expect(result.length).toBe(1);
    expect(result[0]!.eq(coeffs[0]!)).toBe(true);
  });

  test('multiply identity polynomial', () => {
    const F = GF(17n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    const f = x.pow(3).add(x);
    const one = R.one();

    const product = ntt_multiply(f, one, F);
    expect(product.eq(f)).toBe(true);
  });

  test('multiply same polynomial (squaring)', () => {
    const F = GF(17n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    const f = x.add(R.__call__(F.__call__(2n))); // x + 2

    // (x + 2)^2 = x^2 + 4x + 4
    const expected = f.mul(f);
    const product = ntt_multiply(f, f, F);

    expect(product.eq(expected)).toBe(true);
  });

  test('throws for non-power-of-2 length', () => {
    const F = GF(17n);
    const omega = find_primitive_root(4, F);
    const coeffs = [F.__call__(1n), F.__call__(2n), F.__call__(3n)]; // length 3

    expect(() => ntt(coeffs, omega, F)).toThrow();
  });
});

// ============================================================================
// Performance-related tests (correctness with larger sizes)
// ============================================================================

describe('Larger NTT sizes', () => {
  test('NTT roundtrip for size 32 in GF(97)', () => {
    const F = GF(97n);
    const omega = find_primitive_root(32, F);

    const original: FiniteFieldElement[] = [];
    for (let i = 0; i < 32; i++) {
      original.push(F.__call__(BigInt((i * 13 + 7) % 97)));
    }

    const transformed = ntt(original, omega, F);
    const recovered = intt(transformed, omega, F);

    for (let i = 0; i < 32; i++) {
      expect(recovered[i]!.eq(original[i]!)).toBe(true);
    }
  });

  test('NTT roundtrip for size 256 in GF(257)', () => {
    const F = GF(257n);
    const omega = find_primitive_root(256, F);

    const original: FiniteFieldElement[] = [];
    for (let i = 0; i < 256; i++) {
      original.push(F.__call__(BigInt((i * 37 + 11) % 257)));
    }

    const transformed = ntt(original, omega, F);
    const recovered = intt(transformed, omega, F);

    for (let i = 0; i < 256; i++) {
      expect(recovered[i]!.eq(original[i]!)).toBe(true);
    }
  });

  test('polynomial multiplication with degree > 10', () => {
    const F = GF(97n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    // Build polynomials of degree 7
    let f = R.zero();
    let g = R.zero();

    for (let i = 0; i < 8; i++) {
      f = f.add(x.pow(i).scalar_mul(F.__call__(BigInt(i + 1))));
      g = g.add(x.pow(i).scalar_mul(F.__call__(BigInt(8 - i))));
    }

    const expected = f.mul(g);
    const product = ntt_multiply(f, g, F);

    expect(product.eq(expected)).toBe(true);
  });
});

// ============================================================================
// Specific value tests (manually verified)
// ============================================================================

describe('Specific value tests', () => {
  test('NTT of [1, 1, 1, 1] with omega=4 in GF(17)', () => {
    const F = GF(17n);
    // Find a 4th root of unity
    // In GF(17), the multiplicative group has order 16
    // Generator is 3 (3^16 = 1 mod 17)
    // 4th root of unity: g^(16/4) = g^4 = 3^4 = 81 mod 17 = 13
    // But we use find_primitive_root which may give a different one
    const omega = find_primitive_root(4, F);

    const coeffs = [F.__call__(1n), F.__call__(1n), F.__call__(1n), F.__call__(1n)];

    // The DFT of [1,1,1,1] at roots of unity:
    // y_k = sum_{j=0}^3 omega^{jk} for k = 0,1,2,3
    // y_0 = 1 + 1 + 1 + 1 = 4
    // y_1 = 1 + omega + omega^2 + omega^3 = 0 (since omega^4 = 1)
    // y_2 = 1 + omega^2 + omega^4 + omega^6 = 1 + omega^2 + 1 + omega^2 = 2 + 2*omega^2
    //     If omega^2 = -1 (which is 16 in GF(17)), then y_2 = 2 + 2*16 = 2 + 32 = 34 mod 17 = 0
    // y_3 = 1 + omega^3 + omega^6 + omega^9 = 1 + omega^3 + omega^2 + omega = 0

    const result = ntt(coeffs, omega, F);

    // Sum of all ones = 4
    expect(result[0]!.eq(F.__call__(4n))).toBe(true);

    // For other indices, they should sum to 0 due to root of unity property
    expect(result[1]!.eq(F.zero())).toBe(true);
    expect(result[2]!.eq(F.zero())).toBe(true);
    expect(result[3]!.eq(F.zero())).toBe(true);
  });

  test('DFT at omega=1 gives n times the sum', () => {
    // When omega = 1, DFT(a)[k] = sum of all a[j] for all k
    // This is a degenerate case but worth testing

    const F = GF(17n);
    const one = F.one();

    // Using omega = 1 (not primitive, but valid as a root of x^n - 1)
    const coeffs = [F.__call__(1n), F.__call__(2n), F.__call__(3n), F.__call__(4n)];

    // Manual DFT with omega = 1
    const sum = F.__call__(1n + 2n + 3n + 4n); // = 10
    const expected = [sum, sum, sum, sum];

    const result = dft_naive(coeffs, one, F);

    for (let i = 0; i < 4; i++) {
      expect(result[i]!.eq(expected[i]!)).toBe(true);
    }
  });
});

// Fix the typo in one test (double underscore)
// This was an intentional variation but let's keep the test correct

describe('Additional correctness checks', () => {
  test('polynomial multiplication preserves degree bound', () => {
    const F = GF(17n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    // f has degree 3, g has degree 2
    // product should have degree 5
    const f = x.pow(3).add(x);
    const g = x.pow(2).add(R.one());

    const product = ntt_multiply(f, g, F);

    expect(product.degree()).toBe(5);
  });

  test('multiplication is commutative', () => {
    const F = GF(17n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    const f = x.pow(2).add(x.scalar_mul(F.__call__(3n)));
    const g = x.pow(3).add(R.__call__(F.__call__(7n)));

    const fg = ntt_multiply(f, g, F);
    const gf = ntt_multiply(g, f, F);

    expect(fg.eq(gf)).toBe(true);
  });

  test('multiplication is associative', () => {
    const F = GF(97n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    const f = x.add(R.one());
    const g = x.add(R.__call__(F.__call__(2n)));
    const h = x.add(R.__call__(F.__call__(3n)));

    // (f * g) * h
    const fg = ntt_multiply(f, g, F);
    const fg_h = ntt_multiply(fg, h, F);

    // f * (g * h)
    const gh = ntt_multiply(g, h, F);
    const f_gh = ntt_multiply(f, gh, F);

    expect(fg_h.eq(f_gh)).toBe(true);
  });

  test('multiplication distributes over addition', () => {
    const F = GF(97n);
    const [R, x] = PolynomialRingConstructor(F, 'x');

    const f = x.add(R.one());
    const g = x.pow(2);
    const h = x.add(R.__call__(F.__call__(5n)));

    // f * (g + h)
    const gh = g.add(h);
    const f_gh = ntt_multiply(f, gh, F);

    // f * g + f * h
    const fg = ntt_multiply(f, g, F);
    const fh = ntt_multiply(f, h, F);
    const fg_plus_fh = fg.add(fh);

    expect(f_gh.eq(fg_plus_fh)).toBe(true);
  });
});
