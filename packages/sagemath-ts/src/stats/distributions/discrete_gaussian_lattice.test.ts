/**
 * Unit tests for sage/stats/distributions/discrete_gaussian_lattice
 *
 * Tests for the Discrete Gaussian Distribution sampler over lattices (GPV algorithm).
 */
import { describe, expect, test } from 'bun:test';
import { RuntimeError } from '../../errors.js';
import { set_random_seed } from '../../misc/randstate.js';
import { Rational } from '../../rings/rational.js';
import {
  DiscreteGaussianDistributionLatticeSampler,
  DiscreteGaussianDistributionPolynomialSampler,
  DiscreteGaussianLattice,
  DiscreteGaussianPolynomial,
  qfrep,
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

      // Sage exposes the Gaussian parameter as a method
      // (discrete_gaussian_lattice.py:723).
      expect(D.sigma()).toBe(3);
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

      expect(D.cNumeric()).toEqual([5, 5]);
      expect(D.c().map((x) => x.toString())).toEqual(['5', '5']);
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

      // Sage's __repr__ (discrete_gaussian_lattice.py:794-820).
      expect(D.repr()).toBe(
        'Discrete Gaussian sampler with Gaussian parameter \u03c3 = 5.00000000000000, ' +
          'c=(1, 2) over lattice with basis\n\n[1 0]\n[0 1]'
      );
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

    expect(D.sigma()).toBe(5);
    expect(D.cNumeric()).toEqual([1, 2]);
    expect(D.tau).toBe(4);
  });

  test('uses defaults for optional arguments', () => {
    const basis = [
      [1, 0],
      [0, 1],
    ];
    const D = DiscreteGaussianLattice(basis, 5);

    expect(D.sigma()).toBe(5);
    expect(D.cNumeric()).toEqual([0, 0]);
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

/**
 * The features the audit flagged as missing (M146): a matrix Gaussian
 * parameter, Peikert's rounding parameter `r`, Cholesky, offline samples,
 * `_call_non_spherical`, `set_c`/`c`/`sigma`/`f` and
 * `_normalisation_factor_zz`.
 *
 * Every expected value below is quoted from a doctest in
 * `reference/sage/src/sage/stats/distributions/discrete_gaussian_lattice.py`.
 * The locally installed SageMath (10.3) predates the non-spherical sampler, so
 * the vendored source is the oracle for those.
 */
const I = (n: number): number[][] =>
  Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));

describe('sigma(), c(), set_c() and f()', () => {
  test('sigma() returns the scalar or the covariance matrix', () => {
    // discrete_gaussian_lattice.py:723-736
    expect(new DiscreteGaussianDistributionLatticeSampler(I(3), { sigma: 3.0 }).sigma()).toBe(3.0);

    // py:520-528: sigma_basis=True means Sigma = S*S^T
    const S = [
      [2, 0, 0],
      [-1, 3, 0],
      [2, -1, 1],
    ];
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), {
      sigma: S,
      sigma_basis: true,
    });
    expect(D.sigma()).toEqual([
      [4, -2, 4],
      [-2, 10, -5],
      [4, -5, 6],
    ]);
    expect(D.is_spherical).toBe(false);
  });

  test('a scaled identity Sigma is treated as spherical', () => {
    // py:475-482: DGL(ZZ^3, sigma=Matrix(3, 3, 2)) prints sigma = 2.0
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), {
      sigma: [
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 2],
      ],
    });
    expect(D.is_spherical).toBe(true);
    expect(D.sigma()).toBe(2);
  });

  test("a non positive definite Sigma is rejected with Sage's message", () => {
    // py:511-518
    expect(
      () =>
        new DiscreteGaussianDistributionLatticeSampler(I(2), {
          sigma: [
            [0, 1],
            [1, 0],
          ],
        })
    ).toThrow(
      'Sigma(=[0.000000000000000  1.00000000000000]\n' +
        '[ 1.00000000000000 0.000000000000000]) is not positive definite'
    );
    // Sage raises RuntimeError here (discrete_gaussian_lattice.py:570).
    expect(
      () =>
        new DiscreteGaussianDistributionLatticeSampler(I(2), {
          sigma: [
            [0, 1],
            [1, 0],
          ],
        })
    ).toThrow(RuntimeError);
  });

  test('set_c changes the center and the repr', () => {
    // py:758-772
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), { sigma: 3.0, c: [1, 0, 0] });
    expect(D.c().map((x) => x.toString())).toEqual(['1', '0', '0']);
    D.set_c([2, 0, 0]);
    expect(D.repr()).toBe(
      'Discrete Gaussian sampler with Gaussian parameter σ = 3.00000000000000, ' +
        'c=(2, 0, 0) over lattice with basis\n\n[1 0 0]\n[0 1 0]\n[0 0 1]'
    );
  });

  test('f() matches the non-spherical doctest values', () => {
    // py:704-709
    const Sigma = [
      [5, -2, 4],
      [-2, 10, -5],
      [4, -5, 5],
    ];
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), { sigma: Sigma });
    expect(D.f([1, 0, 1])).toBeCloseTo(0.802518797962478, 14);
    expect(D.f([1, 0, 3])).toBeCloseTo(0.00562800641440405, 16);
  });

  test('f() is the spherical Gaussian when sigma is a scalar', () => {
    const D = new DiscreteGaussianDistributionLatticeSampler(I(2), { sigma: 3.0, c: [1, 1] });
    expect(D.f([1, 1])).toBe(1);
    expect(D.f([2, 1])).toBeCloseTo(Math.exp(-1 / 18), 15);
  });

  test('repr of a matrix Gaussian parameter', () => {
    // py:805-814
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), {
      sigma: [
        [10, -6, 1],
        [-6, 5, -1],
        [1, -1, 2],
      ],
    });
    expect(D.repr()).toBe(
      'Discrete Gaussian sampler with Gaussian parameter Σ =\n' +
        '[ 10.0000000000000 -6.00000000000000  1.00000000000000]\n' +
        '[-6.00000000000000  5.00000000000000 -1.00000000000000]\n' +
        '[ 1.00000000000000 -1.00000000000000  2.00000000000000], c=(0, 0, 0) ' +
        'over lattice with basis\n\n[1 0 0]\n[0 1 0]\n[0 0 1]'
    );
  });
});

describe('_maximal_r and the non-spherical sampler', () => {
  const Sigma = [
    [5, -2, 4],
    [-2, 10, -5],
    [4, -5, 5],
  ];

  test('_maximal_r matches the doctest and keeps Sigma - r^2 Q positive definite', () => {
    // py:366-375
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), {
      sigma: Sigma,
      c: [7, 2, 5],
    });
    const r = D._maximal_r();
    expect(r).toBeGreaterThan(0.58402);
    expect(r).toBeLessThan(0.58403);

    // assert all(e_val >= -1e-12 for e_val in (D.sigma() - r^2*D.Q).eigenvalues())
    const M = Sigma.map((row, i) => row.map((v, j) => v - r * r * (i === j ? 1 : 0)));
    // Positive semidefinite <=> a Cholesky factorisation of M + eps*I exists.
    const n = 3;
    const L = Array.from({ length: n }, () => new Array<number>(n).fill(0));
    let psd = true;
    for (let i = 0; i < n && psd; i++) {
      for (let j = 0; j <= i && psd; j++) {
        let acc = M[i]![j]! + (i === j ? 1e-9 : 0);
        for (let k = 0; k < j; k++) acc -= L[i]![k]! * L[j]![k]!;
        if (i === j) {
          if (acc <= 0) psd = false;
          else L[i]![j] = Math.sqrt(acc);
        } else {
          L[i]![j] = acc / L[j]![j]!;
        }
      }
    }
    expect(psd).toBe(true);
  });

  test('_randomise stays within a sane radius', () => {
    // py:400-405: all(D._randomise([0,0,0]).norm() <= 16 for _ in range(100))
    set_random_seed(0);
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), { sigma: Sigma });
    for (let i = 0; i < 100; i++) {
      const v = D._randomise([0, 0, 0]);
      const norm = Math.sqrt(v.reduce((a, x) => a + Number(x) * Number(x), 0));
      expect(norm).toBeLessThanOrEqual(16);
    }
  });

  test('add_offline_samples and _call_non_spherical produce the right mean', () => {
    // py:894-905: mean of 2^12 samples is within 0.25 of the centre.
    set_random_seed(0);
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), {
      sigma: Sigma,
      c: [0.5, 0, 0],
    });
    expect(D.is_spherical).toBe(false);
    D.add_offline_samples(2 ** 12);
    expect(D.offline_samples.length).toBe(2 ** 12);

    const N = 2 ** 12;
    const mean = [0, 0, 0];
    for (let i = 0; i < N; i++) {
      const v = D._call_non_spherical();
      for (let j = 0; j < 3; j++) mean[j]! += v[j]!.toNumber() / N;
    }
    const c = D.cNumeric();
    const dist = Math.sqrt(mean.reduce((a, m, j) => a + (m - c[j]!) ** 2, 0));
    expect(dist).toBeLessThan(0.25);
    // The offline pool was consumed by the sampling loop.
    expect(D.offline_samples.length).toBe(0);
  });

  test('the empirical covariance of the non-spherical sampler is Sigma', () => {
    set_random_seed(1);
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), { sigma: Sigma });
    const N = 20000;
    const cov = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ];
    for (let t = 0; t < N; t++) {
      const v = D.sample().map(Number);
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) cov[i]![j]! += (v[i]! * v[j]!) / N;
      }
    }
    // 5 standard errors of the sample covariance of a Gaussian with this
    // Sigma is comfortably below 0.6 for N = 20000.
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        expect(Math.abs(cov[i]![j]! - Sigma[i]![j]!)).toBeLessThan(0.6);
      }
    }
  });

  test('an r that is too large is rejected', () => {
    expect(
      () => new DiscreteGaussianDistributionLatticeSampler(I(3), { sigma: Sigma, r: 10 })
    ).toThrow('is not positive definite');
  });
});

describe('_normalisation_factor_zz', () => {
  /** Sage's `...` doctest ellipsis: the printed value starts with `prefix`. */
  const expectPrefix = (value: number, prefix: string): void => {
    expect(String(value).startsWith(prefix)).toBe(true);
  };

  test('matches the doctest values', () => {
    // py:207-246
    expectPrefix(
      new DiscreteGaussianDistributionLatticeSampler(I(3), {
        sigma: 1.0,
      })._normalisation_factor_zz(),
      '15.7496'
    );
    expectPrefix(
      new DiscreteGaussianDistributionLatticeSampler(I(8), {
        sigma: 0.5,
      })._normalisation_factor_zz(3),
      '3.1653'
    );
    expectPrefix(
      new DiscreteGaussianDistributionLatticeSampler(I(8), {
        sigma: 0.5,
      })._normalisation_factor_zz(),
      '6.8249'
    );
    // py:450-451
    expectPrefix(
      new DiscreteGaussianDistributionLatticeSampler(I(2), {
        sigma: 3.0,
      })._normalisation_factor_zz(),
      '56.5486677646'
    );
    // py:244-246 (Sage prints the 100-bit rounding 1558545456544038969634991553;
    // our reals are doubles, so we only claim 15 significant digits).
    const big = new DiscreteGaussianDistributionLatticeSampler(I(8), {
      sigma: 1000,
    })._normalisation_factor_zz();
    expect(Math.abs(big / 1.558545456544039e27 - 1)).toBeLessThan(1e-14);
  });

  test('the non-spherical branch matches the doctest', () => {
    // py:255-258 and py:487-501
    const Sigma = [
      [5, -2, 4],
      [-2, 10, -5],
      [4, -5, 5],
    ];
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), {
      sigma: Sigma,
      c: [7, 2, 5],
    });
    const nf = D._normalisation_factor_zz();
    expect(String(nf).startsWith('78.6804')).toBe(true);
    // "We can compute the expected number of samples before sampling a vector"
    expect(String(1 / (D.f([11, 4, 8]) / nf)).startsWith('2553.9461')).toBe(true);
  });

  test('raises exactly as Sage does', () => {
    // py:248-277
    expect(() =>
      new DiscreteGaussianDistributionLatticeSampler(
        [
          [1, 3, 0],
          [-2, 5, 1],
        ],
        { sigma: 3 }
      )._normalisation_factor_zz()
    ).toThrow('basis must be a square matrix');

    expect(() =>
      new DiscreteGaussianDistributionLatticeSampler(I(3), {
        sigma: 1,
        c: [0.5, 0, 0],
      })._normalisation_factor_zz()
    ).toThrow('center must be at zero and basis must be trivial');

    expect(() =>
      new DiscreteGaussianDistributionLatticeSampler(
        [
          [0.5, 0, 0],
          [0, 0.5, 0],
          [0, 0, 0.5],
        ],
        { sigma: 1 }
      )._normalisation_factor_zz()
    ).toThrow('lattice must be integral');

    const M = [
      [1, 3, 0],
      [-2, 5, 1],
      [3, -4, 2],
    ];
    expect(() =>
      new DiscreteGaussianDistributionLatticeSampler(M, { sigma: 1.7 })._normalisation_factor_zz()
    ).toThrow('center must be at zero and basis must be trivial');
  });

  test('is consistent with the empirical sampling frequency', () => {
    // py:213-236: m*f(v)/nf/counter[v] -> 1
    set_random_seed(0);
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), { sigma: 1.0 });
    const nf = D._normalisation_factor_zz();
    const m = 200000;
    const counter = new Map<string, number>();
    for (let i = 0; i < m; i++) {
      const key = D.sample().join(',');
      counter.set(key, (counter.get(key) ?? 0) + 1);
    }
    for (const v of [
      [0, 0, 0],
      [1, 0, 0],
      [-1, 1, 0],
      [1, 1, 1],
    ]) {
      const observed = counter.get(v.join(',')) ?? 0;
      const expected = (m * D.f(v)) / nf;
      expect(observed).toBeGreaterThan(0);
      // 5 Poisson standard deviations.
      expect(Math.abs(observed - expected)).toBeLessThan(5 * Math.sqrt(expected));
    }
  });
});

describe('_call_simple / _call dispatch', () => {
  test('a trivial basis with an integral centre uses _call_simple', () => {
    // py:618-624 sets up the D/VS fast path.
    set_random_seed(0);
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), { sigma: 3.0, c: [1, 0, 0] });
    const N = 4096;
    const mean = [0, 0, 0];
    for (let i = 0; i < N; i++) {
      const v = D._call_simple();
      for (let j = 0; j < 3; j++) mean[j]! += v[j]!.toNumber() / N;
    }
    // py:828-832: norm(mean - c) < 0.25
    const c = D.cNumeric();
    expect(Math.sqrt(mean.reduce((a, m, j) => a + (m - c[j]!) ** 2, 0))).toBeLessThan(0.25);
  });

  test('a fractional centre uses the general GPV loop', () => {
    // py:848-852
    set_random_seed(0);
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), {
      sigma: 3.0,
      c: [0.5, 0, 0],
    });
    const N = 4096;
    const mean = [0, 0, 0];
    for (let i = 0; i < N; i++) {
      const v = D.sampleExact();
      for (let j = 0; j < 3; j++) mean[j]! += v[j]!.toNumber() / N;
    }
    const c = D.cNumeric();
    expect(Math.sqrt(mean.reduce((a, m, j) => a + (m - c[j]!) ** 2, 0))).toBeLessThan(0.25);
  });

  test('a skew basis samples with balanced coordinates', () => {
    // py:680-686
    set_random_seed(0);
    const D = new DiscreteGaussianDistributionLatticeSampler(
      [
        [1, 2],
        [0, 1],
      ],
      { sigma: 20.0 }
    );
    const N = 4096;
    let sx = 0;
    let sy = 0;
    for (let i = 0; i < N; i++) {
      const v = D.sample();
      sx += Math.abs(Number(v[0]!)) / N;
      sy += Math.abs(Number(v[1]!)) / N;
    }
    expect(sx / sy).toBeGreaterThan(0.9);
    expect(sx / sy).toBeLessThan(1.1);
  });
});

describe('qfrep', () => {
  test('matches PARI qfrep(Q, B, 0)', () => {
    // Golden values produced by PARI 2.15 through SageMath 10.3:
    //   Q.__pari__().qfrep(B, 0)
    const I = (n: number): number[][] =>
      Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
    const cases: [number[][], number, number[]][] = [
      [I(3), 10, [3, 6, 4, 3, 12, 12, 0, 6, 15, 12]],
      [I(5), 15, [5, 20, 40, 45, 56, 120, 160, 100, 125, 280, 280, 200, 280, 400, 480]],
      [I(8), 5, [8, 56, 224, 568, 1008]],
      [
        [
          [2, 1],
          [1, 3],
        ],
        12,
        [0, 1, 2, 0, 0, 0, 2, 1, 0, 1, 0, 2],
      ],
      [
        [
          [1, 1],
          [1, 2],
        ],
        20,
        [2, 2, 0, 2, 4, 0, 0, 2, 2, 4, 0, 0, 4, 0, 0, 2, 4, 2, 0, 4],
      ],
      [
        [
          [4, 1, 0],
          [1, 5, 2],
          [0, 2, 6],
        ],
        25,
        [0, 0, 0, 1, 1, 1, 2, 0, 1, 2, 1, 0, 1, 0, 1, 1, 2, 2, 1, 2, 2, 2, 1, 1, 1],
      ],
    ];
    for (const [Q, b, want] of cases) {
      expect(qfrep(Q, b).map(Number)).toEqual(want);
    }
  });

  test('counts every vector exactly once up to sign, by brute force', () => {
    const Q = [
      [3, 1, 1],
      [1, 4, 2],
      [1, 2, 5],
    ];
    const bound = 20;
    const brute = new Array<number>(bound).fill(0);
    for (let a = -6; a <= 6; a++) {
      for (let b2 = -6; b2 <= 6; b2++) {
        for (let c2 = -6; c2 <= 6; c2++) {
          const x = [a, b2, c2];
          let m = 0;
          for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) m += x[i]! * Q[i]![j]! * x[j]!;
          }
          if (m >= 1 && m <= bound) brute[m - 1]!++;
        }
      }
    }
    expect(qfrep(Q, bound).map(Number)).toEqual(brute.map((v) => v / 2));
  });

  test('returns an empty list for a non-positive bound', () => {
    expect(qfrep([[1]], 0)).toEqual([]);
  });
});
