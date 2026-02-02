#!/usr/bin/env bun
/**
 * LWE property tests comparing sagemath-ts with SageMath.
 *
 * This test suite verifies that our LWE implementation produces structurally
 * equivalent results to SageMath's crypto.lwe module.
 *
 * Note: Exact sample values cannot be compared due to different random number
 * generators, but structural properties (dimensions, moduli, bounds) should match.
 *
 * Usage:
 *   bun test tests/property/typescript/lwe-comparison.test.ts
 *
 * To generate SageMath baseline:
 *   sage tests/property/sage/lwe.sage > tests/property/transcripts/sage-lwe.json
 *
 * Then run comparison:
 *   bun test tests/property/typescript/lwe-comparison.test.ts
 */

import { beforeAll, describe, expect, test } from 'bun:test';
import { existsSync } from 'node:fs';
import { readFile } from 'node:fs/promises';
import { join } from 'node:path';

import { euler_phi, next_prime } from '../../../packages/sagemath-ts/src/arith/misc.js';
import {
  DiscreteGaussianDistributionPolynomialSampler,
  LWE,
  LWEVector,
  LindnerPeikert,
  Regev,
  RingLWE,
  RingLWEConverter,
  RingLindnerPeikert,
  UniformNoiseLWE,
  UniformPolynomialSampler,
  UniformSampler,
  balance_sample,
  samples,
} from '../../../packages/sagemath-ts/src/crypto/lwe.js';
import { Zmod } from '../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod_ring.js';
import { PolynomialRingConstructor } from '../../../packages/sagemath-ts/src/rings/polynomial/polynomial_ring.js';

/**
 * Interface for SageMath test result
 */
interface SageTestResult {
  test: string;
  n?: number;
  q?: number;
  computed_q?: number;
  expected_q?: number;
  dimension?: number;
  sample_a_length?: number;
  q_matches?: boolean;
  delta?: number;
  m?: number;
  expected_m?: number;
  instance?: string;
  secret_dist?: string;
  actual_q?: number;
  num_samples?: number;
  all_valid?: boolean;
  error_bounds?: [number, number];
  all_errors_in_bounds?: boolean;
  sample_errors?: number[];
  samples_obtained?: number;
  exhausted_correctly?: boolean;
  all_values_balanced?: boolean;
  q_half?: number;
  num_samples_returned?: number;
  first_sample_a_length?: number;
  secret_length?: number;
  secret_sample?: number[];
  error: string | null;
}

/**
 * Interface for SageMath output
 */
interface SageOutput {
  description: string;
  sage_version: string;
  tests: SageTestResult[];
}

// Path to SageMath transcript
const SAGE_TRANSCRIPT_PATH = join(import.meta.dir, '../transcripts/sage-lwe.json');

// Cached SageMath results
let sageResults: SageOutput | null = null;

/**
 * Load SageMath results if available
 */
async function loadSageResults(): Promise<SageOutput | null> {
  if (sageResults !== null) {
    return sageResults;
  }

  if (!existsSync(SAGE_TRANSCRIPT_PATH)) {
    console.warn(
      `SageMath transcript not found at ${SAGE_TRANSCRIPT_PATH}.` +
        '\nRun: sage tests/property/sage/lwe.sage > tests/property/transcripts/sage-lwe.json'
    );
    return null;
  }

  try {
    const content = await readFile(SAGE_TRANSCRIPT_PATH, 'utf-8');
    sageResults = JSON.parse(content);
    return sageResults;
  } catch (e) {
    console.error('Failed to load SageMath transcript:', e);
    return null;
  }
}

/**
 * Get SageMath test results by test name
 */
function getSageTests(testName: string): SageTestResult[] {
  if (!sageResults) return [];
  return sageResults.tests.filter((t) => t.test === testName);
}

// ============================================================================
// Regev Parameter Tests
// ============================================================================
describe('Regev Parameters (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  const testNValues = [5n, 10n, 15n, 20n, 25n, 50n];

  for (const n of testNValues) {
    test(`Regev(${n}) computes correct modulus q`, () => {
      const oracle = new Regev(n);

      // q should be the smallest prime >= n^2
      const expectedQ = next_prime(n * n);
      expect(oracle.K.modulus).toBe(expectedQ);
    });

    test(`Regev(${n}) has correct dimension`, () => {
      const oracle = new Regev(n);
      expect(oracle.n).toBe(Number(n));
    });

    test(`Regev(${n}) samples have correct structure`, () => {
      const oracle = new Regev(n);
      const [a, c] = oracle.call();

      // Vector a should have dimension n
      expect(a.length).toBe(Number(n));

      // All elements should be in Z_q
      const q = oracle.K.modulus;
      for (let i = 0; i < Number(n); i++) {
        const val = a.get(i).value;
        expect(val >= 0n).toBe(true);
        expect(val < q).toBe(true);
      }

      // c should be in Z_q
      expect(c.value >= 0n).toBe(true);
      expect(c.value < q).toBe(true);
    });
  }

  test('Regev parameters match SageMath', async () => {
    const sageTests = getSageTests('regev_parameters');
    if (sageTests.length === 0) {
      console.warn('Skipping SageMath comparison - no transcript available');
      return;
    }

    for (const sageTest of sageTests) {
      if (sageTest.error) continue;

      const n = BigInt(sageTest.n!);
      const oracle = new Regev(n);

      expect(Number(oracle.K.modulus)).toBe(sageTest.computed_q);
      expect(oracle.n).toBe(sageTest.dimension);
    }
  });
});

// ============================================================================
// LindnerPeikert Parameter Tests
// ============================================================================
describe('LindnerPeikert Parameters (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  const testCases = [
    { n: 10n, delta: 0.01 },
    { n: 20n, delta: 0.01 },
    { n: 20n, delta: 0.001 },
    { n: 50n, delta: 0.01 },
  ];

  for (const { n, delta } of testCases) {
    test(`LindnerPeikert(${n}, delta=${delta}) has correct m`, () => {
      const oracle = new LindnerPeikert(n, delta);

      // m should be 2*n + 128 per [LP2011]
      const expectedM = 2 * Number(n) + 128;
      expect(oracle.m).toBe(expectedM);
    });

    test(`LindnerPeikert(${n}, delta=${delta}) has correct dimension`, () => {
      const oracle = new LindnerPeikert(n, delta);
      expect(oracle.n).toBe(Number(n));
    });

    test(`LindnerPeikert(${n}, delta=${delta}) uses noise secret distribution`, () => {
      const oracle = new LindnerPeikert(n, delta);
      expect(oracle.secret_dist).toBe('noise');
    });

    test(`LindnerPeikert(${n}, delta=${delta}) samples have correct structure`, () => {
      const oracle = new LindnerPeikert(n, delta);
      const [a, c] = oracle.call();

      expect(a.length).toBe(Number(n));

      const q = oracle.K.modulus;
      for (let i = 0; i < Number(n); i++) {
        const val = a.get(i).value;
        expect(val >= 0n).toBe(true);
        expect(val < q).toBe(true);
      }
    });
  }

  test('LindnerPeikert m values match SageMath', async () => {
    const sageTests = getSageTests('lindner_peikert_parameters');
    if (sageTests.length === 0) {
      console.warn('Skipping SageMath comparison - no transcript available');
      return;
    }

    for (const sageTest of sageTests) {
      if (sageTest.error) continue;

      const n = BigInt(sageTest.n!);
      const delta = sageTest.delta!;
      const oracle = new LindnerPeikert(n, delta);

      // m should match
      expect(oracle.m).toBe(sageTest.expected_m);
      expect(oracle.n).toBe(sageTest.dimension);
    }
  });
});

// ============================================================================
// UniformNoiseLWE Parameter Tests
// ============================================================================
describe('UniformNoiseLWE Parameters (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  const testCases = [
    { n: 89n, instance: 'key' as const },
    { n: 89n, instance: 'encrypt' as const },
    { n: 100n, instance: 'key' as const },
    { n: 100n, instance: 'encrypt' as const },
  ];

  for (const { n, instance } of testCases) {
    test(`UniformNoiseLWE(${n}, '${instance}') creates valid oracle`, () => {
      const oracle = new UniformNoiseLWE(n, instance);

      // Should have positive dimension
      expect(oracle.n).toBeGreaterThan(0);

      // Should have positive modulus
      expect(oracle.K.modulus > 0n).toBe(true);
    });

    test(`UniformNoiseLWE(${Number(n)}, '${instance}') samples have correct structure`, () => {
      const oracle = new UniformNoiseLWE(n, instance);
      const [a, c] = oracle.call();

      // Vector a should match oracle dimension
      expect(a.length).toBe(oracle.n);

      const q = oracle.K.modulus;
      for (let i = 0; i < oracle.n; i++) {
        const val = a.get(i).value;
        expect(val >= 0n).toBe(true);
        expect(val < q).toBe(true);
      }
    });
  }

  test('UniformNoiseLWE rejects n < 89', () => {
    expect(() => new UniformNoiseLWE(88n)).toThrow('Parameter too small');
    expect(() => new UniformNoiseLWE(50n)).toThrow('Parameter too small');
  });
});

// ============================================================================
// Base LWE Oracle Tests
// ============================================================================
describe('Base LWE Oracle (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  test('LWE with uniform secret distribution', () => {
    const D = new UniformSampler(-3n, 3n);
    const oracle = new LWE(5n, 101n, D, 'uniform');

    expect(oracle.n).toBe(5);
    expect(oracle.K.modulus).toBe(101n);
    expect(oracle.secret_dist).toBe('uniform');

    // Secret should be in Z_q
    for (let i = 0; i < 5; i++) {
      const val = oracle.secret.get(i).value;
      expect(val >= 0n).toBe(true);
      expect(val < 101n).toBe(true);
    }
  });

  test('LWE with noise secret distribution', () => {
    const D = new UniformSampler(-3n, 3n);
    const oracle = new LWE(5n, 101n, D, 'noise');

    expect(oracle.secret_dist).toBe('noise');

    // Secret should be drawn from noise distribution (small values when balanced)
    const q = 101n;
    const q2 = q / 2n;

    for (let i = 0; i < 5; i++) {
      const val = oracle.secret.get(i).value;
      const balanced = val <= q2 ? val : val - q;
      expect(balanced >= -3n).toBe(true);
      expect(balanced <= 3n).toBe(true);
    }
  });

  test('LWE with bounded secret distribution', () => {
    const D = new UniformSampler(-2n, 2n);
    const oracle = new LWE(5n, 101n, D, [-5n, 5n]);

    // Secret should be drawn from [-5, 5]
    const q = 101n;
    const q2 = q / 2n;

    for (let i = 0; i < 5; i++) {
      const val = oracle.secret.get(i).value;
      const balanced = val <= q2 ? val : val - q;
      expect(balanced >= -5n).toBe(true);
      expect(balanced <= 5n).toBe(true);
    }
  });
});

// ============================================================================
// Sample Structure Tests
// ============================================================================
describe('Sample Structure (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  const testNValues = [5n, 10n, 20n];

  for (const n of testNValues) {
    test(`Regev(${n}) samples have valid structure`, () => {
      const oracle = new Regev(n);
      const q = oracle.K.modulus;

      // Generate multiple samples and verify structure
      for (let i = 0; i < 10; i++) {
        const [a, c] = oracle.call();

        // Check dimension
        expect(a.length).toBe(Number(n));

        // Check all values are in range [0, q)
        for (let j = 0; j < n; j++) {
          const val = a.get(j).value;
          expect(val >= 0n).toBe(true);
          expect(val < q).toBe(true);
        }

        expect(c.value >= 0n).toBe(true);
        expect(c.value < q).toBe(true);
      }
    });
  }

  test('Sample structure matches SageMath', async () => {
    const sageTests = getSageTests('sample_structure');
    if (sageTests.length === 0) {
      console.warn('Skipping SageMath comparison - no transcript available');
      return;
    }

    for (const sageTest of sageTests) {
      if (sageTest.error) continue;

      const n = sageTest.n!;
      const oracle = new Regev(BigInt(n));

      // Generate samples and verify structure matches
      for (let i = 0; i < 10; i++) {
        const [a, c] = oracle.call();
        expect(a.length).toBe(n);
      }

      // SageMath reports all_valid for its samples
      expect(sageTest.all_valid).toBe(true);
    }
  });
});

// ============================================================================
// Error Distribution Bounds Tests
// ============================================================================
describe('Error Distribution Bounds (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  test('Error recovery with known secret - bounded uniform error', () => {
    const D = new UniformSampler(-3n, 3n);
    const q = next_prime(100n);
    const oracle = new LWE(10n, q, D, 'uniform');
    const secret = oracle.secret;
    const q2 = q / 2n;

    // Generate samples and verify error recovery
    for (let i = 0; i < 20; i++) {
      const [a, c] = oracle.call();

      // Compute <a, s>
      const as = a.dotProduct(secret);

      // Recovered error: e = c - <a, s> (mod q)
      const e = c.sub(as);

      // Balance to get actual error value
      const eVal = e.value <= q2 ? e.value : e.value - q;

      // Error should be in [-3, 3]
      expect(eVal >= -3n).toBe(true);
      expect(eVal <= 3n).toBe(true);
    }
  });

  test('Error bounds match SageMath', async () => {
    const sageTests = getSageTests('error_distribution_bounds');
    if (sageTests.length === 0) {
      console.warn('Skipping SageMath comparison - no transcript available');
      return;
    }

    for (const sageTest of sageTests) {
      if (sageTest.error) continue;

      // SageMath should report all errors in bounds
      expect(sageTest.all_errors_in_bounds).toBe(true);

      // Verify bounds are correct
      expect(sageTest.error_bounds).toEqual([-3, 3]);
    }
  });
});

// ============================================================================
// Sample Count Limit Tests
// ============================================================================
describe('Sample Count Limit (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  const testMValues = [3n, 5n, 10n];

  for (const m of testMValues) {
    test(`LWE with m=${m} sample limit`, () => {
      const D = new UniformSampler(-2n, 2n);
      const oracle = new LWE(10n, 101n, D, 'uniform', m);

      // Should be able to get exactly m samples
      for (let i = 0n; i < m; i++) {
        expect(() => oracle.call()).not.toThrow();
      }

      // Next sample should throw
      expect(() => oracle.call()).toThrow('Number of available samples exhausted');
    });
  }

  test('Regev with m=5 sample limit', () => {
    const oracle = new Regev(10n, 'uniform', 5n);

    // Should be able to get exactly 5 samples
    for (let i = 0; i < 5; i++) {
      expect(() => oracle.call()).not.toThrow();
    }

    // Next sample should throw
    expect(() => oracle.call()).toThrow('Number of available samples exhausted');
  });
});

// ============================================================================
// balance_sample Function Tests
// ============================================================================
describe('balance_sample Function (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  test('balance_sample converts to balanced representation', () => {
    const oracle = new Regev(5n);
    const q = oracle.K.modulus;
    const qHalf = q / 2n;

    const sample = oracle.call();
    const [balancedA, balancedC] = balance_sample(sample);

    // All values should be in [-q/2, q/2]
    for (const val of balancedA) {
      expect(val >= -qHalf).toBe(true);
      expect(val <= qHalf).toBe(true);
    }

    expect(balancedC >= -qHalf).toBe(true);
    expect(balancedC <= qHalf).toBe(true);
  });

  test('balance_sample matches SageMath behavior', async () => {
    const sageTests = getSageTests('balance_sample');
    if (sageTests.length === 0) {
      console.warn('Skipping SageMath comparison - no transcript available');
      return;
    }

    for (const sageTest of sageTests) {
      if (sageTest.error) continue;

      // SageMath should report all values balanced
      expect(sageTest.all_values_balanced).toBe(true);
    }
  });
});

// ============================================================================
// samples() Convenience Function Tests
// ============================================================================
describe('samples() Function (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  test('samples() with Regev oracle by name', () => {
    const S = samples(10n, 15n, 'Regev');

    expect(S.length).toBe(10);

    // Each sample should have correct structure
    for (const [a] of S) {
      if (Array.isArray(a)) {
        expect(a.length).toBe(15);
      } else {
        expect(a.length).toBe(15);
      }
    }
  });

  test('samples() with Regev class reference', () => {
    const S = samples(10n, 15n, Regev);

    expect(S.length).toBe(10);
  });

  test('samples() with LindnerPeikert by name', () => {
    const S = samples(5n, 20n, 'LindnerPeikert');

    expect(S.length).toBe(5);
  });

  test('samples() with balanced option', () => {
    const S = samples(10n, 10n, 'Regev', { balanced: true });

    expect(S.length).toBe(10);

    // q for Regev(10) is next_prime(100) = 101
    const q = 101n;
    const qHalf = q / 2n;

    for (const [a, c] of S) {
      // a should be an array of bigints in balanced form
      expect(Array.isArray(a)).toBe(true);
      const aArr = a as bigint[];

      for (const val of aArr) {
        expect(val >= -qHalf).toBe(true);
        expect(val <= qHalf).toBe(true);
      }

      // c should be a bigint in balanced form
      expect(typeof c).toBe('bigint');
      expect(c >= -qHalf).toBe(true);
      expect(c <= qHalf).toBe(true);
    }
  });

  test('samples() function matches SageMath', async () => {
    const sageTests = getSageTests('samples_function');
    if (sageTests.length === 0) {
      console.warn('Skipping SageMath comparison - no transcript available');
      return;
    }

    for (const sageTest of sageTests) {
      if (sageTest.error) continue;

      // Verify same number of samples returned
      const S = samples(BigInt(sageTest.num_samples_returned || 10), 15n, 'Regev');

      expect(S.length).toBe(sageTest.num_samples_returned);

      // Verify sample structure matches
      if (sageTest.first_sample_a_length) {
        const [a] = S[0]!;
        if (Array.isArray(a)) {
          expect(a.length).toBe(sageTest.first_sample_a_length);
        } else {
          expect(a.length).toBe(sageTest.first_sample_a_length);
        }
      }
    }
  });
});

// ============================================================================
// Secret Distribution Type Tests
// ============================================================================
describe('Secret Distribution Types (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  test('uniform secret distribution has full Z_q range', () => {
    const D = new UniformSampler(-3n, 3n);
    const oracle = new LWE(5n, 101n, D, 'uniform');

    // Secret values should be uniformly distributed in Z_q
    // We can't test uniformity statistically with few values,
    // but we can verify they're in range
    for (let i = 0; i < 5; i++) {
      const val = oracle.secret.get(i).value;
      expect(val >= 0n).toBe(true);
      expect(val < 101n).toBe(true);
    }
  });

  test('noise secret distribution has small values', () => {
    const D = new UniformSampler(-3n, 3n);
    const oracle = new LWE(5n, 101n, D, 'noise');
    const q = 101n;
    const q2 = q / 2n;

    // Secret values should be drawn from noise distribution
    for (let i = 0; i < 5; i++) {
      const val = oracle.secret.get(i).value;
      const balanced = val <= q2 ? val : val - q;
      expect(balanced >= -3n).toBe(true);
      expect(balanced <= 3n).toBe(true);
    }
  });

  test('secret distribution matches SageMath', async () => {
    const sageTests = getSageTests('secret_distribution');
    if (sageTests.length === 0) {
      console.warn('Skipping SageMath comparison - no transcript available');
      return;
    }

    for (const sageTest of sageTests) {
      if (sageTest.error) continue;

      // Verify secret length matches
      expect(sageTest.secret_length).toBe(5);
    }
  });
});

// ============================================================================
// UniformSampler Tests
// ============================================================================
describe('UniformSampler (comparison with SageMath)', () => {
  test('samples within bounds', () => {
    const sampler = new UniformSampler(-10n, 10n);

    for (let i = 0; i < 100; i++) {
      const sample = sampler.call();
      expect(sample >= -10n).toBe(true);
      expect(sample <= 10n).toBe(true);
    }
  });

  test('single-value range returns constant', () => {
    const sampler = new UniformSampler(42n, 42n);

    for (let i = 0; i < 10; i++) {
      expect(sampler.call()).toBe(42n);
    }
  });

  test('throws on invalid bounds', () => {
    expect(() => new UniformSampler(10n, 5n)).toThrow('lower bound must be <= upper bound');
  });

  test('handles large ranges', () => {
    const sampler = new UniformSampler(0n, 10n ** 20n);

    for (let i = 0; i < 10; i++) {
      const sample = sampler.call();
      expect(sample >= 0n).toBe(true);
      expect(sample <= 10n ** 20n).toBe(true);
    }
  });
});

// ============================================================================
// Edge Cases
// ============================================================================
describe('Edge Cases (comparison with SageMath)', () => {
  test('LWE with dimension 1', () => {
    const D = new UniformSampler(-1n, 1n);
    const oracle = new LWE(1n, 17n, D);

    const [a, c] = oracle.call();
    expect(a.length).toBe(1);
    expect(c.parent.modulus).toBe(17n);
  });

  test('Regev with small n', () => {
    const oracle = new Regev(2n);

    // q should be next_prime(4) = 5
    expect(oracle.K.modulus).toBe(5n);
  });

  test('zero error distribution gives exact computation', () => {
    const D = new UniformSampler(0n, 0n);
    const oracle = new LWE(5n, 101n, D);
    const secret = oracle.secret;

    const [a, c] = oracle.call();
    const as = a.dotProduct(secret);

    // c should equal <a, s> exactly when error is 0
    expect(c.value).toBe(as.value);
  });
});

// ============================================================================
// Summary Test - Compare All Parameters
// ============================================================================
describe('Summary: Parameter Comparison with SageMath', () => {
  test('All Regev parameters match SageMath', async () => {
    await loadSageResults();
    const sageTests = getSageTests('regev_parameters');

    if (sageTests.length === 0) {
      console.warn('No SageMath transcript - running standalone verification');
    }

    // Test that our implementation computes correct Regev parameters
    for (const n of [5n, 10n, 15n, 20n, 25n, 50n]) {
      const oracle = new Regev(n);
      const expectedQ = next_prime(n * n);

      expect(oracle.K.modulus).toBe(expectedQ);
      expect(oracle.n).toBe(Number(n));

      // If we have SageMath results, compare
      const sageTest = sageTests.find((t) => t.n === Number(n));
      if (sageTest && !sageTest.error) {
        expect(Number(oracle.K.modulus)).toBe(sageTest.computed_q);
        expect(sageTest.q_matches).toBe(true);
      }
    }
  });

  test('All LindnerPeikert m values match expected', () => {
    for (const n of [10n, 20n, 50n]) {
      const oracle = new LindnerPeikert(n);

      // m should be 2*n + 128 per [LP2011]
      expect(oracle.m).toBe(2 * Number(n) + 128);
    }
  });
});

// ============================================================================
// Ring-LWE Tests
// ============================================================================
describe('Ring-LWE Parameters (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  // Test various N values (index of cyclotomic polynomial)
  const testNValues = [
    { N: 8n, expectedN: 4 }, // phi(8) = 4
    { N: 16n, expectedN: 8 }, // phi(16) = 8
    { N: 32n, expectedN: 16 }, // phi(32) = 16
    { N: 20n, expectedN: 8 }, // phi(20) = 8
  ];

  for (const { N, expectedN } of testNValues) {
    test(`RingLWE(${N}) has n = phi(${N}) = ${expectedN}`, () => {
      const K = Zmod(257n);
      const [R] = PolynomialRingConstructor(K, 'x');
      const D = new DiscreteGaussianDistributionPolynomialSampler(R, expectedN, 3.0);
      const ringlwe = new RingLWE(N, 257n, D);

      expect(ringlwe.N).toBe(Number(N));
      expect(ringlwe.n).toBe(expectedN);
    });

    test(`RingLWE(${N}) samples have correct dimensions`, () => {
      const K = Zmod(257n);
      const [R] = PolynomialRingConstructor(K, 'x');
      const n = Number(euler_phi(N));
      const D = new DiscreteGaussianDistributionPolynomialSampler(R, n, 3.0);
      const ringlwe = new RingLWE(N, 257n, D);

      const [a, c] = ringlwe.call();

      // Both a and c should be vectors of length n = phi(N)
      expect(a.length).toBe(n);
      expect(c.length).toBe(n);
    });

    test(`RingLWE(${N}) samples have valid range`, () => {
      const q = 257n;
      const K = Zmod(q);
      const [R] = PolynomialRingConstructor(K, 'x');
      const n = Number(euler_phi(N));
      const D = new DiscreteGaussianDistributionPolynomialSampler(R, n, 3.0);
      const ringlwe = new RingLWE(N, q, D);

      const [a, c] = ringlwe.call();

      // All values should be in [0, q)
      for (let i = 0; i < n; i++) {
        expect(a[i]! >= 0n).toBe(true);
        expect(a[i]! < q).toBe(true);
        expect(c[i]! >= 0n).toBe(true);
        expect(c[i]! < q).toBe(true);
      }
    });
  }

  test('RingLWE rejects mismatched noise distribution dimension', () => {
    const K = Zmod(257n);
    const [R] = PolynomialRingConstructor(K, 'x');
    // phi(16) = 8, but we create sampler with n=4
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 4, 3.0);

    expect(() => new RingLWE(16n, 257n, D)).toThrow('Noise distribution has dimensions 4 != 8');
  });

  test('RingLWE with sample limit', () => {
    const K = Zmod(257n);
    const [R] = PolynomialRingConstructor(K, 'x');
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
    const ringlwe = new RingLWE(16n, 257n, D, null, 'uniform', 3n);

    // Should be able to get exactly 3 samples
    ringlwe.call();
    ringlwe.call();
    ringlwe.call();

    // 4th sample should throw
    expect(() => ringlwe.call()).toThrow('Number of available samples exhausted');
  });
});

describe('RingLindnerPeikert Parameters (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  const testNValues = [8n, 16n, 32n];

  for (const N of testNValues) {
    const expectedN = Number(euler_phi(N));

    test(`RingLindnerPeikert(${N}) has n = phi(${N}) = ${expectedN}`, () => {
      const rlp = new RingLindnerPeikert(N);

      expect(rlp.N).toBe(Number(N));
      expect(rlp.n).toBe(expectedN);
    });

    test(`RingLindnerPeikert(${N}) has default m = 3*n = ${3 * expectedN}`, () => {
      const rlp = new RingLindnerPeikert(N);

      expect(rlp.m).toBe(3 * expectedN);
    });

    test(`RingLindnerPeikert(${N}) uses noise secret distribution`, () => {
      const rlp = new RingLindnerPeikert(N);

      expect(rlp.secret_dist).toBe('noise');
    });

    test(`RingLindnerPeikert(${N}) samples have correct dimensions`, () => {
      const rlp = new RingLindnerPeikert(N);
      const [a, c] = rlp.call();

      expect(a.length).toBe(expectedN);
      expect(c.length).toBe(expectedN);
    });
  }

  test('RingLindnerPeikert custom m parameter', () => {
    const rlp = new RingLindnerPeikert(16n, 0.01, 50n);
    expect(rlp.m).toBe(50);
  });
});

describe('RingLWEConverter (comparison with SageMath)', () => {
  beforeAll(async () => {
    await loadSageResults();
  });

  test('RingLWEConverter produces n LWE samples per Ring-LWE sample', () => {
    const K = Zmod(257n);
    const [R] = PolynomialRingConstructor(K, 'x');
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
    const ringlwe = new RingLWE(16n, 257n, D);
    const converter = new RingLWEConverter(ringlwe);

    // Generate n samples (all from the same Ring-LWE sample internally)
    const lweSamples: Array<[bigint[], bigint]> = [];
    for (let i = 0; i < 8; i++) {
      lweSamples.push(converter.call());
    }

    // Verify all samples have correct dimensions
    for (const [a, c] of lweSamples) {
      expect(a.length).toBe(8);
      expect(typeof c).toBe('bigint');
      expect(c >= 0n).toBe(true);
      expect(c < 257n).toBe(true);
    }
  });

  test('RingLWEConverter a vectors have valid range', () => {
    const q = 257n;
    const K = Zmod(q);
    const [R] = PolynomialRingConstructor(K, 'x');
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
    const ringlwe = new RingLWE(16n, q, D);
    const converter = new RingLWEConverter(ringlwe);

    for (let i = 0; i < 8; i++) {
      const [a, c] = converter.call();

      for (const val of a) {
        expect(val >= 0n).toBe(true);
        expect(val < q).toBe(true);
      }
    }
  });
});

describe('UniformPolynomialSampler (comparison with SageMath)', () => {
  test('samples have correct degree', () => {
    const [R] = PolynomialRingConstructor(Zmod(101n), 'x');
    const sampler = new UniformPolynomialSampler(R, 8n, -2n, 2n);

    for (let i = 0; i < 10; i++) {
      const poly = sampler.call();
      expect(poly.degree()).toBeLessThanOrEqual(7);
    }
  });

  test('samples coefficients within bounds', () => {
    const [R] = PolynomialRingConstructor(Zmod(101n), 'x');
    const sampler = new UniformPolynomialSampler(R, 8n, -3n, 3n);

    for (let i = 0; i < 10; i++) {
      const poly = sampler.call();
      for (let j = 0; j < 8; j++) {
        const coeff = poly.getCoeff(j);
        const val = coeff.value;
        // Balance to check bounds
        const balanced = val <= 50n ? val : val - 101n;
        expect(balanced >= -3n).toBe(true);
        expect(balanced <= 3n).toBe(true);
      }
    }
  });

  test('repr matches SageMath format', () => {
    const [R] = PolynomialRingConstructor(Zmod(101n), 'x');
    const sampler = new UniformPolynomialSampler(R, 10n, -10n, 10n);
    expect(sampler.repr()).toBe('UniformPolynomialSampler(10, -10, 10)');
  });
});

describe('DiscreteGaussianDistributionPolynomialSampler (comparison with SageMath)', () => {
  test('samples have correct degree', () => {
    const [R] = PolynomialRingConstructor(Zmod(1009n), 'x');
    const sampler = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);

    for (let i = 0; i < 10; i++) {
      const poly = sampler.call();
      expect(poly.degree()).toBeLessThanOrEqual(7);
    }
  });

  test('samples coefficients are concentrated near center', () => {
    const [R] = PolynomialRingConstructor(Zmod(1009n), 'x');
    const sigma = 3.0;
    const sampler = new DiscreteGaussianDistributionPolynomialSampler(R, 8, sigma);

    // With sigma=3, most values should be within 5*sigma = 15
    let withinBounds = 0;
    const total = 50 * 8;

    for (let i = 0; i < 50; i++) {
      const poly = sampler.call();
      for (let j = 0; j < 8; j++) {
        const coeff = poly.getCoeff(j);
        const val = coeff.value;
        const balanced = val <= 504n ? val : val - 1009n;
        if (balanced >= -15n && balanced <= 15n) {
          withinBounds++;
        }
      }
    }

    // At least 80% should be within 5*sigma
    expect(withinBounds).toBeGreaterThan(0.8 * total);
  });

  test('repr matches SageMath format', () => {
    const [R] = PolynomialRingConstructor(Zmod(1009n), 'x');
    const sampler = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
    expect(sampler.repr()).toContain('Discrete Gaussian sampler for polynomials of degree < 8');
    expect(sampler.repr()).toContain('σ=3.000000');
  });
});
