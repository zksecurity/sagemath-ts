/**
 * Unit tests for FFT/NTT convolution module
 */
import { describe, expect, test } from 'bun:test';
import { type FiniteFieldElement, FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import {
  NTT_FRIENDLY_PRIMES,
  _convolution_fft,
  _convolution_naive,
  _negaconvolution_fft,
  _negaconvolution_naive,
  convolution,
  convolve_naive,
  dft_naive,
  evaluate_on_domain,
  fft,
  fft_multiply,
  find_ntt_prime,
  find_primitive_root,
  idft_naive,
  ifft,
  interpolate_from_domain,
  intt,
  max_fft_size,
  ntt,
  ntt_multiply,
  supports_ntt_size,
} from './convolution.js';
import { Polynomial } from './polynomial_element.js';
import { PolynomialRingConstructor } from './polynomial_ring.js';

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
    const expected = x
      .mul(x)
      .add(x.scalar_mul(F.__call__(3n)))
      .add(R.__call__(F.__call__(2n)));

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
    const f = x
      .mul(x)
      .add(x.scalar_mul(F.__call__(3n)))
      .add(R.__call__(F.__call__(5n)));

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
    const f = x
      .mul(x)
      .add(x.scalar_mul(F.__call__(3n)))
      .add(R.__call__(F.__call__(2n)));

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
    expect(p >= 1n << 15n).toBe(true);
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

// ============================================================================
// Large fields: the 32-bit shift boundary (M16, L22)
// ============================================================================

describe('max_fft_size beyond the int32 boundary', () => {
  test('Goldilocks 2^64 - 2^32 + 1 supports 2^32', () => {
    const goldilocks = GF(2n ** 64n - 2n ** 32n + 1n);
    expect(max_fft_size(goldilocks)).toBe(32);
    // 1 << 32 overflows a 32-bit shift; the comparison must be done in bigint
    expect(supports_ntt_size(goldilocks.characteristic, 1 << 20)).toBe(true);
  });

  test('a prime with v_2(p-1) = 31', () => {
    // 35 * 2^31 + 1 = 75161927681 is prime
    const F = GF(35n * 2n ** 31n + 1n);
    expect(max_fft_size(F)).toBe(31);
    expect(supports_ntt_size(F.characteristic, 1 << 20)).toBe(true);
  });

  test('ntt_multiply works over Goldilocks (M16)', () => {
    const goldilocks = GF(2n ** 64n - 2n ** 32n + 1n);
    const [R, x] = PolynomialRingConstructor(goldilocks, 'x');
    const f = x.pow(2).add(x).add(R.one());
    const g = x.add(R.__call__(goldilocks.__call__(2n)));
    const product = ntt_multiply(f, g, goldilocks);
    expect(product.eq(f.mul(g))).toBe(true);
    expect(product.toString()).toBe('x^3 + 3*x^2 + 3*x + 2');
  });

  test('ntt_multiply over a 2^31-friendly prime', () => {
    const F = GF(35n * 2n ** 31n + 1n);
    const [R, x] = PolynomialRingConstructor(F, 'x');
    const f = x.pow(3).sub(R.one());
    const g = x.pow(2).add(x).add(R.one());
    expect(ntt_multiply(f, g, F).eq(f.mul(g))).toBe(true);
  });
});

// ============================================================================
// fft / ifft (L22: these were imported but never exercised)
// ============================================================================

describe('fft and ifft', () => {
  test('fft matches the naive DFT in GF(17)', () => {
    const F = GF(17n);
    const omega = find_primitive_root(4, F);
    const coeffs = [F.__call__(1n), F.__call__(2n), F.__call__(3n), F.__call__(4n)];
    const fast = fft(coeffs, omega);
    const slow = dft_naive(coeffs, omega, F);
    expect(fast.length).toBe(slow.length);
    for (let i = 0; i < slow.length; i++) {
      expect(fast[i]!.eq(slow[i]!)).toBe(true);
    }
  });

  test('ifft inverts fft in GF(257)', () => {
    const F = GF(257n);
    const n = 8;
    const omega = find_primitive_root(n, F);
    const coeffs: FiniteFieldElement[] = [];
    for (let i = 0; i < n; i++) coeffs.push(F.__call__(BigInt(3 * i + 1)));
    const transformed = fft(coeffs, omega);
    const back = ifft(transformed, omega);
    const nInv = F.__call__(BigInt(n)).inv();
    for (let i = 0; i < n; i++) {
      expect(back[i]!.mul(nInv).eq(coeffs[i]!)).toBe(true);
    }
  });

  test('fft of a single element is the identity', () => {
    const F = GF(17n);
    const c = [F.__call__(5n)];
    expect(fft(c, F.one())[0]!.eq(c[0]!)).toBe(true);
    expect(ifft(c, F.one())[0]!.eq(c[0]!)).toBe(true);
  });
});

// ============================================================================
// Generic convolution (port of sage/rings/polynomial/convolution.py, M18)
// ============================================================================

describe('convolution (sage.rings.polynomial.convolution)', () => {
  test("Sage's doctests", () => {
    // sage: convolution([1, 2, 3, 4, 5], [6, 7]) == [6, 19, 32, 45, 58, 35]
    expect(convolution([1n, 2n, 3n, 4n, 5n], [6n, 7n])).toEqual([6n, 19n, 32n, 45n, 58n, 35n]);
    // sage: convolution([1, 2, 3], [4, 5, 6, 7]) == [4, 13, 28, 34, 32, 21]
    expect(convolution([1n, 2n, 3n], [4n, 5n, 6n, 7n])).toEqual([4n, 13n, 28n, 34n, 32n, 21n]);
    // sage: _convolution_naive([4, 5, 6, 7], [1, 2, 3]) == [4, 13, 28, 34, 32, 21]
    expect(_convolution_naive([4n, 5n, 6n, 7n], [1n, 2n, 3n])).toEqual([
      4n,
      13n,
      28n,
      34n,
      32n,
      21n,
    ]);
    // sage: _negaconvolution_naive([1, 2, 3], [3, 4, 5]) == [-19, -5, 22]
    expect(_negaconvolution_naive([1n, 2n, 3n], [3n, 4n, 5n])).toEqual([-19n, -5n, 22n]);
    // sage: _convolution_fft([1, 2, 3], [4, 5, 6]) == [4, 13, 28, 27, 18]
    expect(_convolution_fft([1n, 2n, 3n], [4n, 5n, 6n])).toEqual([4n, 13n, 28n, 27n, 18n]);
    // sage: _negaconvolution_fft(range(8), range(5, 13), 3)
    expect(
      _negaconvolution_fft([0n, 1n, 2n, 3n, 4n, 5n, 6n, 7n], [5n, 6n, 7n, 8n, 9n, 10n, 11n, 12n], 3)
    ).toEqual([-224n, -234n, -224n, -192n, -136n, -54n, 56n, 196n]);
  });

  test('fft and naive agree (Sage TESTS blocks)', () => {
    let state = 12345n;
    const rnd = (m: bigint) => {
      state = (state * 6364136223846793005n + 1442695040888963407n) % (1n << 64n);
      return (state >> 17n) % m;
    };

    for (let n = 3; n <= 8; n++) {
      const L1 = Array.from({ length: 1 << n }, () => rnd(100n));
      const L2 = Array.from({ length: 1 << n }, () => rnd(100n));
      expect(_negaconvolution_fft(L1, L2, n)).toEqual(_negaconvolution_naive(L1, L2));
    }

    for (let len1 = 4; len1 < 20; len1++) {
      for (let len2 = 4; len2 < 20; len2++) {
        const L1 = Array.from({ length: len1 }, () => rnd(100n));
        const L2 = Array.from({ length: len2 }, () => rnd(100n));
        expect(_convolution_fft(L1, L2)).toEqual(_convolution_naive(L1, L2));
      }
    }
  });

  test('long inputs take the FFT path and stay correct', () => {
    // sage: L3 = convolution(L1, L2); L3[2000] == sum(L1[i]*L2[2000-i] ...)
    const L1 = Array.from({ length: 1000 }, (_, i) => BigInt(i));
    const L2 = Array.from({ length: 3756 }, (_, i) => BigInt((i * 7) % 47));
    const L3 = convolution(L1, L2);
    expect(L3.length).toBe(1000 + 3756 - 1);
    let expected = 0n;
    for (let i = 0; i < 1000; i++) expected += L1[i]! * L2[2000 - i]!;
    expect(L3[2000]).toBe(expected);
  });

  test('convolution over a finite field agrees with polynomial multiplication', () => {
    const F = GF(257n);
    const L1 = Array.from({ length: 120 }, (_, i) => F.__call__(BigInt((7 * i + 1) % 257)));
    const L2 = Array.from({ length: 130 }, (_, i) => F.__call__(BigInt((11 * i + 5) % 257)));
    const conv = convolution(L1, L2);
    const [R] = PolynomialRingConstructor(F, 'x');
    const product = new Polynomial(L1, R).mul(new Polynomial(L2, R));
    for (let i = 0; i < conv.length; i++) {
      expect(conv[i]!.eq(product.getCoeff(i))).toBe(true);
    }
  });

  test('empty input raises ValueError', () => {
    expect(() => convolution([], [1n])).toThrow('cannot compute convolution of empty lists');
    expect(() => convolution([1n], [])).toThrow('cannot compute convolution of empty lists');
  });
});
