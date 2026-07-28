/**
 * Unit tests for sage/stats/distributions/discrete_gaussian_lattice
 *
 * Tests for the Discrete Gaussian Distribution sampler over lattices (GPV algorithm).
 */
import { describe, expect, test } from 'bun:test';
import { set_random_seed } from '../../misc/randstate.js';
import { Rational } from '../../rings/rational.js';
import {
  DiscreteGaussianDistributionLatticeSampler,
  DiscreteGaussianDistributionPolynomialSampler,
  DiscreteGaussianLattice,
  DiscreteGaussianPolynomial,
  samplePreimage,
  sampleShortVector,
} from './discrete_gaussian_lattice.js';

describe('DiscreteGaussianDistributionLatticeSampler', () => {
  describe('construction', () => {
    test('creates sampler with identity basis', () => {
      const basis = [
        [1, 0],
        [0, 1],
      ];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 3 });

      expect(D.sigma).toBe(3);
      expect(D.rank).toBe(2);
      expect(D.degree).toBe(2);
      expect(D.tau).toBe(6);
    });

    test('creates sampler with custom center', () => {
      const basis = [
        [1, 0],
        [0, 1],
      ];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, {
        sigma: 3,
        c: [5, 5],
      });

      expect(D.c).toEqual([5, 5]);
    });

    test('creates sampler with bigint basis', () => {
      const basis = [
        [1n, 0n, 0n],
        [0n, 1n, 0n],
        [0n, 0n, 1n],
      ];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 5 });

      expect(D.rank).toBe(3);
      expect(D.degree).toBe(3);
    });

    test('throws on empty basis', () => {
      expect(() => new DiscreteGaussianDistributionLatticeSampler([], { sigma: 3 })).toThrow(
        'basis must be a non-empty array'
      );
    });

    test('throws on invalid sigma', () => {
      const basis = [
        [1, 0],
        [0, 1],
      ];
      expect(() => new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 0 })).toThrow(
        'sigma must be > 0'
      );
      expect(() => new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: -1 })).toThrow(
        'sigma must be > 0'
      );
    });

    test('throws on mismatched basis dimensions', () => {
      const basis = [
        [1, 0],
        [0, 1, 0], // Wrong dimension
      ];
      expect(() => new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 3 })).toThrow(
        'all basis vectors must have the same dimension'
      );
    });

    test('throws on mismatched center dimension', () => {
      const basis = [
        [1, 0],
        [0, 1],
      ];
      expect(
        () => new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 3, c: [1, 2, 3] })
      ).toThrow('c must be a vector of dimension 2');
    });
  });

  describe('sampling from Z^n', () => {
    test('samples from Z^2', () => {
      const basis = [
        [1, 0],
        [0, 1],
      ];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 3 });

      const sample = D.sample();
      expect(sample.length).toBe(2);
      expect(typeof sample[0]).toBe('bigint');
      expect(typeof sample[1]).toBe('bigint');
    });

    test('samples are lattice vectors', () => {
      // Z^2 lattice - all integer vectors are valid
      const basis = [
        [1, 0],
        [0, 1],
      ];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 3 });

      for (let i = 0; i < 50; i++) {
        const sample = D.sample();
        // For Z^2, every integer vector is a lattice vector
        expect(typeof sample[0]).toBe('bigint');
        expect(typeof sample[1]).toBe('bigint');
      }
    });

    test('mean is approximately center', () => {
      const basis = [
        [1, 0],
        [0, 1],
      ];
      const center = [5, -3];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, {
        sigma: 5,
        c: center,
      });

      const samples = D.samples(2000);
      const mean = [0, 0];
      for (const s of samples) {
        mean[0] += Number(s[0]!);
        mean[1] += Number(s[1]!);
      }
      mean[0] /= samples.length;
      mean[1] /= samples.length;

      // Mean should be close to center (within 1)
      expect(Math.abs(mean[0]! - center[0]!)).toBeLessThan(1);
      expect(Math.abs(mean[1]! - center[1]!)).toBeLessThan(1);
    });
  });

  describe('sampling from non-trivial lattice', () => {
    test('samples from 2D lattice with basis [[2,0],[1,2]]', () => {
      const basis = [
        [2, 0],
        [1, 2],
      ];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 5 });

      for (let i = 0; i < 50; i++) {
        const sample = D.sample();
        // Sample should be a linear combination of basis vectors
        // sample = a * [2,0] + b * [1,2] = [2a+b, 2b]
        // So sample[1] should be even
        expect(sample[1]! % 2n).toBe(0n);
      }
    });

    test('samples from 3D lattice', () => {
      const basis = [
        [1, 0, 0],
        [0, 2, 0],
        [0, 0, 3],
      ];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 5 });

      for (let i = 0; i < 30; i++) {
        const sample = D.sample();
        expect(sample.length).toBe(3);
        // Second coordinate should be even (multiple of 2)
        expect(sample[1]! % 2n).toBe(0n);
        // Third coordinate should be multiple of 3
        expect(sample[2]! % 3n).toBe(0n);
      }
    });
  });

  describe('call method', () => {
    test('call is alias for sample', () => {
      const basis = [
        [1, 0],
        [0, 1],
      ];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 3 });

      const sample = D.call();
      expect(sample.length).toBe(2);
      expect(typeof sample[0]).toBe('bigint');
    });
  });

  describe('samples method', () => {
    test('generates correct number of samples', () => {
      const basis = [
        [1, 0],
        [0, 1],
      ];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 3 });

      const samples = D.samples(25);
      expect(samples.length).toBe(25);
    });
  });

  describe('smoothing parameter', () => {
    test('computes smoothing parameter', () => {
      const basis = [
        [1, 0],
        [0, 1],
      ];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 10 });

      const eta = D.smoothingParameter();
      expect(eta).toBeGreaterThan(0);
      expect(Number.isFinite(eta)).toBe(true);
    });

    test('isAboveSmoothingParameter', () => {
      const basis = [
        [1, 0],
        [0, 1],
      ];

      // Large sigma should be above smoothing parameter
      const D1 = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 100 });
      expect(D1.isAboveSmoothingParameter()).toBe(true);

      // Very small sigma might not be (depends on epsilon)
      const D2 = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 0.1 });
      // This might be false for small sigma
      expect(typeof D2.isAboveSmoothingParameter()).toBe('boolean');
    });
  });

  describe('repr and toString', () => {
    test('repr contains key information', () => {
      const basis = [
        [1, 0],
        [0, 1],
      ];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, {
        sigma: 5,
        c: [1, 2],
      });

      const repr = D.repr();
      expect(repr).toContain('DiscreteGaussianDistributionLatticeSampler');
      expect(repr).toContain('rank=2');
      expect(repr).toContain('sigma=5');
    });

    test('toString equals repr', () => {
      const basis = [
        [1, 0],
        [0, 1],
      ];
      const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 3 });
      expect(D.toString()).toBe(D.repr());
    });
  });
});

describe('DiscreteGaussianLattice factory', () => {
  test('creates sampler with positional arguments', () => {
    const basis = [
      [1, 0],
      [0, 1],
    ];
    const D = DiscreteGaussianLattice(basis, 5, [1, 2], 4n);

    expect(D.sigma).toBe(5);
    expect(D.c).toEqual([1, 2]);
    expect(D.tau).toBe(4);
  });

  test('uses defaults for optional arguments', () => {
    const basis = [
      [1, 0],
      [0, 1],
    ];
    const D = DiscreteGaussianLattice(basis, 5);

    expect(D.sigma).toBe(5);
    expect(D.c).toEqual([0, 0]);
    expect(D.tau).toBe(6);
  });
});

describe('sampleShortVector', () => {
  test('samples a lattice vector', () => {
    const basis = [
      [1, 0],
      [0, 1],
    ];
    const v = sampleShortVector(basis, 3);

    expect(v.length).toBe(2);
    expect(typeof v[0]).toBe('bigint');
    expect(typeof v[1]).toBe('bigint');
  });
});

describe('samplePreimage', () => {
  test('samples near target', () => {
    const basis = [
      [1, 0],
      [0, 1],
    ];
    const target = [10, 20];

    // With large sigma, samples should be near target
    const v = samplePreimage(basis, 5, target);

    expect(v.length).toBe(2);
    // Should be reasonably close to target (within few sigma)
    expect(Math.abs(Number(v[0]!) - target[0]!)).toBeLessThan(30);
    expect(Math.abs(Number(v[1]!) - target[1]!)).toBeLessThan(30);
  });
});

describe('DiscreteGaussianDistributionPolynomialSampler', () => {
  describe('construction', () => {
    test('creates polynomial sampler', () => {
      const D = new DiscreteGaussianDistributionPolynomialSampler(8, { sigma: 3 });

      expect(D.n).toBe(8);
      expect(D.sigma).toBe(3);
    });

    test('throws on invalid n', () => {
      expect(() => new DiscreteGaussianDistributionPolynomialSampler(0, { sigma: 3 })).toThrow(
        'n must be > 0'
      );
      expect(() => new DiscreteGaussianDistributionPolynomialSampler(-5, { sigma: 3 })).toThrow(
        'n must be > 0'
      );
    });
  });

  describe('sampling', () => {
    test('samples polynomial of correct degree', () => {
      const D = new DiscreteGaussianDistributionPolynomialSampler(8, { sigma: 3 });

      const poly = D.sample();
      expect(poly.length).toBe(8);
      poly.forEach((coeff) => expect(typeof coeff).toBe('bigint'));
    });

    test('call is alias for sample', () => {
      const D = new DiscreteGaussianDistributionPolynomialSampler(4, { sigma: 3 });

      const poly = D.call();
      expect(poly.length).toBe(4);
    });

    test('samples generates multiple polynomials', () => {
      const D = new DiscreteGaussianDistributionPolynomialSampler(8, { sigma: 3 });

      const polys = D.samples(10);
      expect(polys.length).toBe(10);
      polys.forEach((poly) => expect(poly.length).toBe(8));
    });

    test('coefficients follow discrete Gaussian', () => {
      const sigma = 5;
      const D = new DiscreteGaussianDistributionPolynomialSampler(4, { sigma });

      // Sample many polynomials and collect all coefficients
      const coeffs: number[] = [];
      for (let i = 0; i < 500; i++) {
        const poly = D.sample();
        poly.forEach((c) => coeffs.push(Number(c)));
      }

      // Compute variance
      const mean = coeffs.reduce((a, b) => a + b, 0) / coeffs.length;
      const variance = coeffs.reduce((a, b) => a + (b - mean) * (b - mean), 0) / coeffs.length;

      // Variance should be close to sigma^2
      expect(Math.abs(variance - sigma * sigma)).toBeLessThan(sigma * sigma * 0.15);
    });
  });

  describe('repr', () => {
    test('contains parameters', () => {
      const D = new DiscreteGaussianDistributionPolynomialSampler(8, {
        sigma: 3.5,
        c: 1n,
      });

      const repr = D.repr();
      expect(repr).toContain('DiscreteGaussianDistributionPolynomialSampler');
      expect(repr).toContain('n=8');
      expect(repr).toContain('sigma=3.5');
    });
  });
});

describe('DiscreteGaussianPolynomial factory', () => {
  test('creates sampler with positional arguments', () => {
    const D = DiscreteGaussianPolynomial(8, 5, 2n, 4n);

    expect(D.n).toBe(8);
    expect(D.sigma).toBe(5);
    expect(D.c).toBe(2);
  });

  test('uses defaults for optional arguments', () => {
    const D = DiscreteGaussianPolynomial(8, 5);

    expect(D.n).toBe(8);
    expect(D.sigma).toBe(5);
    expect(D.c).toBe(0);
  });
});

describe('exact (non-integral) lattices', () => {
  // Sage keeps the basis over an exact ring and computes
  // v = sum z_i B[i] exactly (discrete_gaussian_lattice.py:842-872).
  // Rounding the basis to integers changes both the lattice and the width.
  test('samples of (1/2 Z)^2 stay in (1/2 Z)^2', () => {
    set_random_seed(0);
    const D = new DiscreteGaussianDistributionLatticeSampler(
      [
        [0.5, 0],
        [0, 0.5],
      ],
      { sigma: 3 }
    );
    expect(D.isIntegral).toBe(false);

    const samples = D.samplesExact(2048);
    let sawHalf = false;
    for (const v of samples) {
      for (const x of v) {
        // 2*x must be an integer, i.e. x lies in (1/2)Z ...
        expect(x.mul(new Rational(2n)).isInteger()).toBe(true);
        // ... and it must not be forced into Z.
        if (!x.isInteger()) sawHalf = true;
      }
    }
    expect(sawHalf).toBe(true);
  });

  test('(1/2 Z)^2 sampler has the requested width, not twice it', () => {
    set_random_seed(1);
    const sigma = 3;
    const D = new DiscreteGaussianDistributionLatticeSampler(
      [
        [0.5, 0],
        [0, 0.5],
      ],
      { sigma }
    );
    const samples = D.samplesExact(4096);
    const xs = samples.map((v) => v[0]!.toNumber());
    const mean = xs.reduce((a, b) => a + b, 0) / xs.length;
    const stddev = Math.sqrt(xs.reduce((a, b) => a + (b - mean) * (b - mean), 0) / xs.length);
    expect(Math.abs(mean)).toBeLessThan(0.25);
    // Rounding the basis to the identity would give stddev ~= 6.
    expect(Math.abs(stddev - sigma)).toBeLessThan(0.3);
  });

  test('sample() refuses a non-integral lattice', () => {
    const D = new DiscreteGaussianDistributionLatticeSampler(
      [
        [0.5, 0],
        [0, 0.5],
      ],
      { sigma: 3 }
    );
    expect(() => D.sample()).toThrow('lattice basis is not integral');
  });

  test('accepts an exact rational basis', () => {
    const D = new DiscreteGaussianDistributionLatticeSampler(
      [
        [new Rational(1n, 3n), new Rational(0n)],
        [new Rational(0n), new Rational(1n, 3n)],
      ],
      { sigma: 2 }
    );
    expect(D.basisExact[0]![0]!.eq(new Rational(1n, 3n))).toBe(true);
    for (const x of D.sampleExact()) {
      expect(x.mul(new Rational(3n)).isInteger()).toBe(true);
    }
  });

  test('integral bases still return bigint vectors', () => {
    const D = new DiscreteGaussianDistributionLatticeSampler(
      [
        [1, 0],
        [0, 1],
      ],
      { sigma: 3 }
    );
    expect(D.isIntegral).toBe(true);
    for (const x of D.sample()) {
      expect(typeof x).toBe('bigint');
    }
  });
});

describe('SageMath doctests', () => {
  // discrete_gaussian_lattice.py:665-680
  test('ZZ^3, sigma=3, c=(1,0,0): norm(mean - c) < 0.25', () => {
    set_random_seed(0);
    const D = new DiscreteGaussianDistributionLatticeSampler(
      [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
      ],
      { sigma: 3, c: [1, 0, 0] }
    );
    const L = D.samples(4096);
    const mean = [0, 0, 0];
    for (const v of L) {
      for (let j = 0; j < 3; j++) mean[j]! += Number(v[j]!);
    }
    for (let j = 0; j < 3; j++) mean[j]! /= L.length;
    const dist = Math.hypot(mean[0]! - 1, mean[1]!, mean[2]!);
    expect(dist).toBeLessThan(0.25);
  });

  test('ZZ^3, sigma=3, c=(1/2,0,0): norm(mean - c) < 0.25', () => {
    set_random_seed(0);
    const D = new DiscreteGaussianDistributionLatticeSampler(
      [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
      ],
      { sigma: 3, c: [new Rational(1n, 2n), new Rational(0n), new Rational(0n)] }
    );
    const L = D.samples(4096);
    const mean = [0, 0, 0];
    for (const v of L) {
      for (let j = 0; j < 3; j++) mean[j]! += Number(v[j]!);
    }
    for (let j = 0; j < 3; j++) mean[j]! /= L.length;
    const dist = Math.hypot(mean[0]! - 0.5, mean[1]!, mean[2]!);
    expect(dist).toBeLessThan(0.25);
  });

  test('M = [[1,2],[0,1]], sigma=20: 0.9 < mean|x| / mean|y| < 1.1', () => {
    set_random_seed(0);
    const D = new DiscreteGaussianDistributionLatticeSampler(
      [
        [1, 2],
        [0, 1],
      ],
      { sigma: 20 }
    );
    const L = D.samples(4096);
    let ax = 0;
    let ay = 0;
    for (const v of L) {
      ax += Math.abs(Number(v[0]!));
      ay += Math.abs(Number(v[1]!));
    }
    const div = ax / ay;
    expect(div).toBeGreaterThan(0.9);
    expect(div).toBeLessThan(1.1);
  });
});

describe('GPV algorithm properties', () => {
  test('samples are integer vectors', () => {
    const basis = [
      [2, 1],
      [0, 3],
    ];
    const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 10 });

    for (let i = 0; i < 100; i++) {
      const sample = D.sample();
      sample.forEach((x) => {
        expect(typeof x).toBe('bigint');
      });
    }
  });

  test('larger sigma gives more spread', () => {
    const basis = [
      [1, 0],
      [0, 1],
    ];
    const D1 = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 2 });
    const D2 = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: 10 });

    // Compute sample variances
    const computeVariance = (D: DiscreteGaussianDistributionLatticeSampler) => {
      const samples = D.samples(500);
      const x = samples.map((s) => Number(s[0]!));
      const mean = x.reduce((a, b) => a + b, 0) / x.length;
      return x.reduce((a, b) => a + (b - mean) * (b - mean), 0) / x.length;
    };

    const var1 = computeVariance(D1);
    const var2 = computeVariance(D2);

    // Larger sigma should give larger variance
    expect(var2).toBeGreaterThan(var1);
  });
});

describe('integration with cryptographic applications', () => {
  test('can sample error polynomials for Ring-LWE', () => {
    // Typical Ring-LWE parameters
    const n = 256;
    const sigma = 3.19;

    const D = new DiscreteGaussianDistributionPolynomialSampler(n, { sigma });
    const error = D.sample();

    expect(error.length).toBe(n);
    // Errors should be small (bounded by tau*sigma)
    const bound = Math.ceil(6 * sigma);
    for (const e of error) {
      expect(Number(e) >= -bound).toBe(true);
      expect(Number(e) <= bound).toBe(true);
    }
  });

  test('can sample from q-ary lattice basis', () => {
    // Simple example: basis for a q-ary lattice
    const q = 97;
    const basis = [
      [q, 0],
      [0, q],
    ];

    const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma: q / 4 });
    const sample = D.sample();

    // Sample should be in qZ^2
    expect(sample[0]! % BigInt(q)).toBe(0n);
    expect(sample[1]! % BigInt(q)).toBe(0n);
  });
});
