/**
 * Unit tests for sage/stats/distributions/discrete_gaussian_integer
 *
 * Tests for the Discrete Gaussian Distribution sampler over integers.
 */
import { describe, expect, test } from 'bun:test';
import {
  DiscreteGaussianDistributionIntegerSampler,
  DiscreteGaussianInteger,
  klDivergence,
  statisticalDistance,
} from './discrete_gaussian_integer.js';

describe('DiscreteGaussianDistributionIntegerSampler', () => {
  describe('construction', () => {
    test('creates sampler with default parameters', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });

      expect(D.sigma).toBe(3);
      expect(D.c).toBe(0);
      expect(D.tau).toBe(6);
    });

    test('creates sampler with custom center', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 5,
        c: 10n,
      });

      expect(D.sigma).toBe(5);
      expect(D.c).toBe(10);
    });

    test('creates sampler with custom tau', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3,
        tau: 4n,
      });

      expect(D.tau).toBe(4);
    });

    test('creates sampler with explicit algorithm', () => {
      const D1 = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3,
        algorithm: 'uniform+table',
      });
      const D2 = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3,
        algorithm: 'uniform+online',
      });

      expect(D1.algorithm).toBe('uniform+table');
      expect(D2.algorithm).toBe('uniform+online');
    });

    test('throws on invalid sigma', () => {
      expect(() => new DiscreteGaussianDistributionIntegerSampler({ sigma: 0 })).toThrow(
        'sigma must be > 0'
      );
      expect(() => new DiscreteGaussianDistributionIntegerSampler({ sigma: -1 })).toThrow(
        'sigma must be > 0'
      );
    });

    test('throws on invalid tau', () => {
      expect(() => new DiscreteGaussianDistributionIntegerSampler({ sigma: 3, tau: 0.5 })).toThrow(
        'tau must be >= 1'
      );
    });

    test('throws on missing sigma', () => {
      // @ts-expect-error - Testing invalid input
      expect(() => new DiscreteGaussianDistributionIntegerSampler({})).toThrow('sigma is required');
    });
  });

  describe('sampling range', () => {
    test('computes correct bounds', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3,
        c: 0n,
        tau: 2n,
      });

      // Range should be [-ceil(3*2), ceil(3*2)] = [-6, 6]
      expect(D.lowerBound).toBe(-6n);
      expect(D.upperBound).toBe(6n);
    });

    test('handles non-integer center', () => {
      // Using DiscreteGaussianOptionsInternal allows number for c (for GPV algorithm)
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3,
        c: 1.5, // Internal API accepts number for fractional centers
        tau: 2n,
      });

      // floor(1.5) = 1, halfWidth = ceil(3*2) = 6
      // Range should be [1-6, 1+6] = [-5, 7]
      expect(D.lowerBound).toBe(-5n);
      expect(D.upperBound).toBe(7n);
    });

    test('samples stay within bounds', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3,
        c: 5n,
        tau: 3n,
      });

      for (let i = 0; i < 100; i++) {
        const sample = D.sample();
        expect(sample >= D.lowerBound).toBe(true);
        expect(sample <= D.upperBound).toBe(true);
      }
    });
  });

  describe('sampling', () => {
    test('sample returns bigint', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });
      const sample = D.sample();
      expect(typeof sample).toBe('bigint');
    });

    test('call is alias for sample', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });
      const sample = D.call();
      expect(typeof sample).toBe('bigint');
    });

    test('samples generates array of correct length', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });
      const samples = D.samples(50);
      expect(samples.length).toBe(50);
      samples.forEach((s) => expect(typeof s).toBe('bigint'));
    });

    test('mean is approximately c', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 5,
        c: 7n,
      });

      const samples = D.samples(10000);
      const mean = Number(samples.reduce((a, b) => a + b, 0n)) / samples.length;

      // Mean should be close to c=7 (within 0.5)
      expect(Math.abs(mean - 7)).toBeLessThan(0.5);
    });

    test('variance is approximately sigma^2', () => {
      const sigma = 4;
      const D = new DiscreteGaussianDistributionIntegerSampler({ sigma });

      const samples = D.samples(10000);
      const sampleNums = samples.map((s) => Number(s));
      const mean = sampleNums.reduce((a, b) => a + b, 0) / sampleNums.length;
      const variance =
        sampleNums.reduce((a, b) => a + (b - mean) * (b - mean), 0) / sampleNums.length;

      // Variance should be close to sigma^2=16 (within 10%)
      expect(Math.abs(variance - sigma * sigma)).toBeLessThan(sigma * sigma * 0.1);
    });

    test('samples are concentrated near center', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3,
        c: 0n,
      });

      const samples = D.samples(1000);
      let withinOneSigma = 0;
      let withinTwoSigma = 0;

      for (const s of samples) {
        const sNum = Number(s);
        if (Math.abs(sNum) <= 3) withinOneSigma++;
        if (Math.abs(sNum) <= 6) withinTwoSigma++;
      }

      // About 68% should be within 1 sigma, 95% within 2 sigma
      expect(withinOneSigma / 1000).toBeGreaterThan(0.55);
      expect(withinTwoSigma / 1000).toBeGreaterThan(0.9);
    });

    test('different algorithms produce similar distributions', () => {
      const D1 = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 5,
        algorithm: 'uniform+table',
      });
      const D2 = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 5,
        algorithm: 'uniform+online',
      });

      const samples1 = D1.samples(5000);
      const samples2 = D2.samples(5000);

      // Compute sample means
      const mean1 = Number(samples1.reduce((a, b) => a + b, 0n)) / samples1.length;
      const mean2 = Number(samples2.reduce((a, b) => a + b, 0n)) / samples2.length;

      // Means should be close (within 0.5)
      expect(Math.abs(mean1 - mean2)).toBeLessThan(0.5);
    });
  });

  describe('probability computations', () => {
    test('rho is highest at center', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3,
        c: 5n,
      });

      const rho5 = D.rho(5);
      const rho4 = D.rho(4);
      const rho6 = D.rho(6);
      const rho0 = D.rho(0);

      expect(rho5).toBeGreaterThan(rho4);
      expect(rho5).toBeGreaterThan(rho6);
      expect(rho5).toBeGreaterThan(rho0);
    });

    test('rho is symmetric around center', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3,
        c: 0n,
      });

      expect(Math.abs(D.rho(2) - D.rho(-2))).toBeLessThan(1e-10);
      expect(Math.abs(D.rho(5) - D.rho(-5))).toBeLessThan(1e-10);
    });

    test('probability sums to 1 (approximately)', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3,
        tau: 10n, // Large tau to capture most of the distribution
      });

      let totalProb = 0;
      for (let x = D.lowerBound; x <= D.upperBound; x++) {
        totalProb += D.probability(x);
      }

      expect(Math.abs(totalProb - 1)).toBeLessThan(1e-10);
    });

    test('probability is 0 outside support', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3,
        tau: 2n,
      });

      expect(D.probability(D.lowerBound - 1n)).toBe(0);
      expect(D.probability(D.upperBound + 1n)).toBe(0);
    });
  });

  describe('statistical properties', () => {
    test('mean matches theoretical value', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 4,
        c: 3n,
      });

      const computedMean = D.mean();
      // Mean should be very close to c
      expect(Math.abs(computedMean - 3)).toBeLessThan(0.01);
    });

    test('variance matches theoretical value', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 5,
        tau: 10n,
      });

      const computedVariance = D.variance();
      // Variance should be close to sigma^2
      expect(Math.abs(computedVariance - 25)).toBeLessThan(1);
    });

    test('stddev matches sigma', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 7,
        tau: 10n,
      });

      const computedStddev = D.stddev();
      expect(Math.abs(computedStddev - 7)).toBeLessThan(0.5);
    });
  });

  describe('support', () => {
    test('support returns correct range', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 2,
        c: 0n,
        tau: 2n,
      });

      const support = D.support();
      expect(support.length).toBe(Number(D.upperBound - D.lowerBound + 1n));
      expect(support[0]).toBe(D.lowerBound);
      expect(support[support.length - 1]).toBe(D.upperBound);
    });
  });

  describe('repr and toString', () => {
    test('repr contains parameters', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3.5,
        c: 2n,
        tau: 5n,
      });

      const repr = D.repr();
      expect(repr).toContain('DiscreteGaussianDistributionIntegerSampler');
      expect(repr).toContain('sigma=3.5');
      expect(repr).toContain('c=2');
      expect(repr).toContain('tau=5');
    });

    test('toString equals repr', () => {
      const D = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });
      expect(D.toString()).toBe(D.repr());
    });
  });

  describe('withOptions', () => {
    test('creates modified copy', () => {
      const D1 = new DiscreteGaussianDistributionIntegerSampler({
        sigma: 3,
        c: 0n,
      });
      const D2 = D1.withOptions({ c: 5n });

      expect(D1.c).toBe(0);
      expect(D2.c).toBe(5);
      expect(D2.sigma).toBe(3); // Unchanged
    });
  });
});

describe('DiscreteGaussianInteger factory', () => {
  test('creates sampler with positional arguments', () => {
    const D = DiscreteGaussianInteger(3, 5n, 4n);

    expect(D.sigma).toBe(3);
    expect(D.c).toBe(5);
    expect(D.tau).toBe(4);
  });

  test('uses defaults for optional arguments', () => {
    const D = DiscreteGaussianInteger(3);

    expect(D.sigma).toBe(3);
    expect(D.c).toBe(0);
    expect(D.tau).toBe(6);
  });
});

describe('statisticalDistance', () => {
  test('returns 0 for identical distributions', () => {
    const D1 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });
    const D2 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });

    expect(statisticalDistance(D1, D2)).toBeLessThan(1e-10);
  });

  test('returns value in [0, 1]', () => {
    const D1 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 2 });
    const D2 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 5 });

    const dist = statisticalDistance(D1, D2);
    expect(dist).toBeGreaterThanOrEqual(0);
    expect(dist).toBeLessThanOrEqual(1);
  });

  test('increases with parameter differences', () => {
    const D1 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });
    const D2 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 4 });
    const D3 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 6 });

    const dist12 = statisticalDistance(D1, D2);
    const dist13 = statisticalDistance(D1, D3);

    expect(dist13).toBeGreaterThan(dist12);
  });
});

describe('klDivergence', () => {
  test('returns 0 for identical distributions', () => {
    const D1 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });
    const D2 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });

    expect(klDivergence(D1, D2)).toBeLessThan(1e-10);
  });

  test('is non-negative', () => {
    const D1 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });
    const D2 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 4 });

    expect(klDivergence(D1, D2)).toBeGreaterThanOrEqual(0);
  });

  test('is asymmetric', () => {
    const D1 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3, tau: 10 });
    const D2 = new DiscreteGaussianDistributionIntegerSampler({ sigma: 5, tau: 10 });

    const kl12 = klDivergence(D1, D2);
    const kl21 = klDivergence(D2, D1);

    // KL divergence is generally asymmetric
    expect(Math.abs(kl12 - kl21)).toBeGreaterThan(1e-10);
  });
});

describe('integration with LWE', () => {
  test('can be used as error distribution', () => {
    const D = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3.19 });

    // Generate errors like in LWE
    const errors: bigint[] = [];
    for (let i = 0; i < 100; i++) {
      errors.push(D.call());
    }

    // Errors should be small (within tau*sigma of center)
    const bound = BigInt(Math.ceil(D.tau * D.sigma));
    for (const e of errors) {
      expect(e >= -bound).toBe(true);
      expect(e <= bound).toBe(true);
    }
  });
});

describe('edge cases', () => {
  test('works with very small sigma', () => {
    const D = new DiscreteGaussianDistributionIntegerSampler({
      sigma: 0.5,
      tau: 6n,
    });

    // With small sigma, samples should be concentrated
    const samples = D.samples(100);
    const uniqueValues = new Set(samples.map((s) => s.toString()));
    expect(uniqueValues.size).toBeLessThanOrEqual(10);
  });

  test('works with large sigma', () => {
    const D = new DiscreteGaussianDistributionIntegerSampler({
      sigma: 100,
      tau: 3n,
    });

    // Should still sample correctly
    const sample = D.sample();
    expect(typeof sample).toBe('bigint');
    expect(sample >= D.lowerBound).toBe(true);
    expect(sample <= D.upperBound).toBe(true);
  });

  test('works with negative center', () => {
    const D = new DiscreteGaussianDistributionIntegerSampler({
      sigma: 3,
      c: -10n,
    });

    const samples = D.samples(1000);
    const mean = Number(samples.reduce((a, b) => a + b, 0n)) / samples.length;

    // Mean should be close to -10
    expect(Math.abs(mean - -10)).toBeLessThan(0.5);
  });

  test('works with fractional center', () => {
    // Using internal API that accepts number for fractional centers (for GPV algorithm)
    const D = new DiscreteGaussianDistributionIntegerSampler({
      sigma: 3,
      c: 2.7, // Internal API accepts number for fractional centers
    });

    const samples = D.samples(1000);
    const mean = Number(samples.reduce((a, b) => a + b, 0n)) / samples.length;

    // Mean should be close to 2.7
    expect(Math.abs(mean - 2.7)).toBeLessThan(0.5);
  });
});

describe('SageMath compatibility', () => {
  test('matches SageMath parameter naming', () => {
    // In SageMath: DiscreteGaussianDistributionIntegerSampler(sigma=3.0, c=0, tau=6)
    const D = new DiscreteGaussianDistributionIntegerSampler({
      sigma: 3.0,
      c: 0n,
      tau: 6n,
    });

    expect(D.sigma).toBe(3.0);
    expect(D.c).toBe(0);
    expect(D.tau).toBe(6);
  });

  test('supports algorithm parameter like SageMath', () => {
    // SageMath uses 'uniform+table' or 'uniform+online'
    const D1 = new DiscreteGaussianDistributionIntegerSampler({
      sigma: 3,
      algorithm: 'uniform+table',
    });
    const D2 = new DiscreteGaussianDistributionIntegerSampler({
      sigma: 3,
      algorithm: 'uniform+online',
    });

    expect(D1.algorithm).toBe('uniform+table');
    expect(D2.algorithm).toBe('uniform+online');
  });

  test('auto-selects algorithm based on range size', () => {
    // Small range should use table
    const D1 = new DiscreteGaussianDistributionIntegerSampler({
      sigma: 100,
      tau: 6n, // Range ~1200
    });
    expect(D1.algorithm).toBe('uniform+table');

    // Large range should use online
    const D2 = new DiscreteGaussianDistributionIntegerSampler({
      sigma: 200000,
      tau: 6n, // Range > 10^6
    });
    expect(D2.algorithm).toBe('uniform+online');
  });
});
