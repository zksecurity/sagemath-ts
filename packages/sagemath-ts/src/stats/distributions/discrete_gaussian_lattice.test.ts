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
  _iter_vectors,
  DiscreteGaussianDistributionLatticeSampler,
  DiscreteGaussianDistributionPolynomialSampler,
  DiscreteGaussianLattice,
  DiscreteGaussianPolynomial,
  _mp,
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

/** Gauss-Jordan inverse of a small matrix, for the brute-force oracles below. */
const ratInverseSmall = (A: number[][]): number[][] => {
  const n = A.length;
  const M = A.map((r, i) => [...r, ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))]);
  for (let col = 0; col < n; col++) {
    let pivot = col;
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(M[row]![col]!) > Math.abs(M[pivot]![col]!)) pivot = row;
    }
    if (M[pivot]![col] === 0) throw new Error('singular');
    if (pivot !== col) {
      const tmp = M[col]!;
      M[col] = M[pivot]!;
      M[pivot] = tmp;
    }
    const d = M[col]![col]!;
    M[col] = M[col]!.map((x) => x / d);
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const f = M[row]![col]!;
      if (f === 0) continue;
      M[row] = M[row]!.map((x, j) => x - f * M[col]![j]!);
    }
  }
  return M.map((r) => r.slice(n));
};

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
  test('matches the doctest values', () => {
    // py:207-246
    expect(
      new DiscreteGaussianDistributionLatticeSampler(I(3), {
        sigma: 1.0,
      })
        ._normalisation_factor_zz()
        .toString()
    ).toStartWith('15.7496');
    expect(
      new DiscreteGaussianDistributionLatticeSampler(I(8), {
        sigma: 0.5,
      })
        ._normalisation_factor_zz(3)
        .toString()
    ).toStartWith('3.1653');
    expect(
      new DiscreteGaussianDistributionLatticeSampler(I(8), {
        sigma: 0.5,
      })
        ._normalisation_factor_zz()
        .toString()
    ).toStartWith('6.8249');
    // py:450-451
    expect(
      new DiscreteGaussianDistributionLatticeSampler(I(2), {
        sigma: 3.0,
      })
        ._normalisation_factor_zz()
        .toString()
    ).toStartWith('56.5486677646');
    // py:244-246: `round(D._normalisation_factor_zz(prec=100))`.  This is a
    // 28-digit integer, i.e. 91 significant bits — the whole reason the sum is
    // built in RealField(prec) and not in doubles.
    const big = new DiscreteGaussianDistributionLatticeSampler(I(8), {
      sigma: 1000,
    })._normalisation_factor_zz(undefined, 100);
    expect(big.round()).toBe(1558545456544038969634991553n);
    expect(big.precision()).toBe(100);
  });

  /**
   * Every value below was produced by running the vendored
   * `_normalisation_factor_zz` (a verbatim transcription of
   * `reference/sage/src/sage/stats/distributions/discrete_gaussian_lattice.py:190-355`)
   * on the locally installed SageMath 10.3 — its own `RealField`, its own
   * symbolic `pi`/`exp`, its own `LLL` and its own PARI `qfrep`.  (10.3 ships an
   * older `_normalisation_factor_zz`, so the method itself could not be called
   * directly.)  They are compared character for character, i.e. to the full
   * `RealField(prec)` printing width.
   */
  test('matches SageMath character for character, including at prec > 53', () => {
    const DGL = DiscreteGaussianDistributionLatticeSampler;
    const skew = [
      [1, 0],
      [2, 1],
    ];
    const uni = [
      [1, 0, 0],
      [3, 1, 0],
      [-2, 5, 1],
    ];
    const cases: [string, string][] = [
      [new DGL(I(3), { sigma: 1.0 })._normalisation_factor_zz().toString(), '15.7496101985309'],
      [new DGL(I(8), { sigma: 0.5 })._normalisation_factor_zz(3).toString(), '3.16536453178580'],
      [new DGL(I(8), { sigma: 0.5 })._normalisation_factor_zz().toString(), '6.82492448921763'],
      [new DGL(I(2), { sigma: 3.0 })._normalisation_factor_zz().toString(), '56.5486677646163'],
      [
        new DGL(I(8), { sigma: 1000 })._normalisation_factor_zz(undefined, 100).toString(),
        '1.5585454565440389696349915528e27',
      ],
      [
        new DGL(I(8), { sigma: 0.5 })._normalisation_factor_zz(undefined, 100).toString(),
        '6.8249244892176317661265358111',
      ],
      [
        new DGL(I(3), { sigma: 1.0 })._normalisation_factor_zz(undefined, 100).toString(),
        '15.749610198530875444830259920',
      ],
      [
        new DGL(I(3), { sigma: 1.0 })._normalisation_factor_zz(undefined, 200).toString(),
        '15.749610198530875444830259919841133336026975226544294742292',
      ],
      [
        new DGL(I(5), { sigma: 0.8 })._normalisation_factor_zz(undefined, 120).toString(),
        '32.427522816195800227147476189804401',
      ],
      [new DGL(I(4), { sigma: 0.9 })._normalisation_factor_zz(10).toString(), '25.3334099910513'],
      [
        new DGL(I(2), { sigma: 3.0 })._normalisation_factor_zz(undefined, 100).toString(),
        '56.548667764616278292327580899',
      ],
      [
        new DGL(I(3), { sigma: 2.5 })._normalisation_factor_zz(undefined, 100).toString(),
        '246.08765540191280850454134366',
      ],
      [new DGL(I(5), { sigma: 1.5 })._normalisation_factor_zz().toString(), '751.460169579992'],
      [new DGL(skew, { sigma: 0.7 })._normalisation_factor_zz().toString(), '3.07953682400995'],
      [new DGL(skew, { sigma: 2.0 })._normalisation_factor_zz().toString(), '25.1327412287183'],
      [new DGL(uni, { sigma: 0.6 })._normalisation_factor_zz().toString(), '3.41868216714722'],
      [new DGL(uni, { sigma: 0.6 })._normalisation_factor_zz(8).toString(), '3.38945789055808'],
    ];
    for (const [got, want] of cases) {
      expect(got).toBe(want);
    }
  }, 180000);

  test('agrees with exhaustive lattice enumeration', () => {
    // sum_{x in L} exp(-|x|^2 / (2 sigma^2)), summed directly over every
    // lattice point in a large box, versus the Poisson-summation value the
    // method returns.
    const bruteForce = (basis: number[][], sigma: number, radius: number): number => {
      const n = basis.length;
      const d = basis[0]!.length;
      let total = 0;
      const coeff = new Array<number>(n).fill(0);
      const rec = (i: number): void => {
        if (i === n) {
          const pt = new Array<number>(d).fill(0);
          for (let k = 0; k < n; k++) {
            for (let j = 0; j < d; j++) pt[j]! += coeff[k]! * basis[k]![j]!;
          }
          let ns = 0;
          for (let j = 0; j < d; j++) ns += pt[j]! * pt[j]!;
          total += Math.exp(-ns / (2 * sigma * sigma));
          return;
        }
        for (let v = -radius; v <= radius; v++) {
          coeff[i] = v;
          rec(i + 1);
        }
      };
      rec(0);
      return total;
    };
    const cases: [number[][], number, number][] = [
      [I(1), 1.0, 60],
      [I(2), 1.0, 60],
      [I(3), 1.0, 40],
      [I(4), 1.0, 20],
      [I(1), 0.5, 40],
      [I(2), 0.5, 40],
      [I(3), 0.5, 30],
      [I(5), 0.5, 10],
      [I(1), 2.0, 80],
      [I(2), 2.0, 60],
      [I(3), 2.0, 30],
      [I(2), 1.3, 60],
      [I(2), 1.05, 60],
      [I(3), 1.1, 30],
      [I(1), 3.0, 120],
      [I(2), 0.9, 40],
    ];
    for (const [B, sigma, radius] of cases) {
      const nf = new DiscreteGaussianDistributionLatticeSampler(B, {
        sigma,
      })._normalisation_factor_zz();
      const bf = bruteForce(B, sigma, radius);
      expect(Math.abs(nf.toNumber() / bf - 1)).toBeLessThan(1e-12);
    }
  });

  test('agrees with exhaustive enumeration on a skew unimodular basis', () => {
    // These bases span Z^n, so the sum is the same as over Z^n; enumerating the
    // *coefficients* in a box would miss short vectors, so we enumerate ambient
    // integer points instead.
    const bruteAmbient = (basis: number[][], sigma: number, radius: number): number => {
      const n = basis.length;
      const inv = ratInverseSmall(basis);
      let total = 0;
      const x = new Array<number>(n).fill(0);
      const rec = (i: number): void => {
        if (i === n) {
          for (let j = 0; j < n; j++) {
            let s = 0;
            for (let k = 0; k < n; k++) s += x[k]! * inv[k]![j]!;
            if (Math.abs(s - Math.round(s)) > 1e-9) return;
          }
          let ns = 0;
          for (let j = 0; j < n; j++) ns += x[j]! * x[j]!;
          total += Math.exp(-ns / (2 * sigma * sigma));
          return;
        }
        for (let v = -radius; v <= radius; v++) {
          x[i] = v;
          rec(i + 1);
        }
      };
      rec(0);
      return total;
    };
    const cases: [number[][], number, number][] = [
      [
        [
          [1, 0],
          [2, 1],
        ],
        0.7,
        40,
      ],
      [
        [
          [1, 0],
          [2, 1],
        ],
        2.0,
        60,
      ],
      [
        [
          [1, 0, 0],
          [3, 1, 0],
          [-2, 5, 1],
        ],
        0.6,
        25,
      ],
      [
        [
          [1, 0, 0],
          [3, 1, 0],
          [-2, 5, 1],
        ],
        2.0,
        35,
      ],
    ];
    for (const [B, sigma, radius] of cases) {
      const nf = new DiscreteGaussianDistributionLatticeSampler(B, {
        sigma,
      })._normalisation_factor_zz();
      expect(Math.abs(nf.toNumber() / bruteAmbient(B, sigma, radius) - 1)).toBeLessThan(1e-12);
    }
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
    expect(nf.toString()).toStartWith('78.6804');
    // Full-width value from the vendored algorithm run on SageMath 10.3.
    expect(nf.toString()).toBe('78.6804676242625');
    // "We can compute the expected number of samples before sampling a vector"
    expect(String(1 / (D.f([11, 4, 8]) / nf.toNumber())).startsWith('2553.9461')).toBe(true);
  });

  test('the non-spherical branch matches SageMath on further inputs', () => {
    const cases: [number[][], number[][], number[], string][] = [
      [
        I(3),
        [
          [10, -6, 1],
          [-6, 5, -1],
          [1, -1, 2],
        ],
        [0, 0, 0],
        '78.6806810436643',
      ],
      [
        I(2),
        [
          [3, 1],
          [1, 4],
        ],
        [1, -2],
        '20.8389657583938',
      ],
    ];
    for (const [B, Sigma, c, want] of cases) {
      const D = new DiscreteGaussianDistributionLatticeSampler(B, { sigma: Sigma, c });
      expect(D.is_spherical).toBe(false);
      expect(D._normalisation_factor_zz().toString()).toBe(want);
    }
  });

  test('the non-spherical branch is exactly the finite sum upstream computes', () => {
    // py:295-313: sum over `u in [-BOUND, BOUND]^n` of `f((u + base) * B)`,
    // with `base = round(B.LLL().solve_left(c))` and BOUND = 10 for n = 3.
    const Sigma = [
      [5, -2, 4],
      [-2, 10, -5],
      [4, -5, 5],
    ];
    const c = [7, 2, 5];
    const D = new DiscreteGaussianDistributionLatticeSampler(I(3), { sigma: Sigma, c });
    const BOUND = 10;
    let total = 0;
    for (let a = -BOUND; a <= BOUND; a++) {
      for (let b = -BOUND; b <= BOUND; b++) {
        for (let d = -BOUND; d <= BOUND; d++) {
          // B is the identity here, so (u + base) * B = u + base.
          total += D.f([a + c[0]!, b + c[1]!, d + c[2]!]);
        }
      }
    }
    expect(Math.abs(D._normalisation_factor_zz().toNumber() / total - 1)).toBeLessThan(1e-14);
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
      const expected = (m * D.f(v)) / nf.toNumber();
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

  test('floors a non-integral bound, as PARI gfloor does', () => {
    // `_normalisation_factor_zz(tau=3)` on ZZ^8 with sigma = 0.5 asks PARI for
    // `qfrep(Q, 1.5, 0)`, i.e. a bound of 1.
    expect(qfrep(I(8), 1.5).map(Number)).toEqual([8]);
    expect(qfrep(I(8), 1.999).map(Number)).toEqual([8]);
    expect(qfrep(I(8), 2).map(Number)).toEqual([8, 56]);
  });

  test('rejects a non-integral Gram matrix, as PARI does', () => {
    expect(() =>
      qfrep(
        [
          [1.5, 0],
          [0, 1],
        ],
        4
      )
    ).toThrow('incorrect type in qfminim');
  });
});

/**
 * The `RealField(prec)` layer `_normalisation_factor_zz` builds its sum in.
 *
 * MPFR is not vendored under `reference/`, so the oracles here are external:
 * `mpmath` at 400 bits for the transcendental functions, exact `Rational`
 * arithmetic for the correct rounding of `+ - * /`, and
 * `reference/sage/src/sage/rings/real_mpfr.pyx:1897-2149` for the printing.
 */
describe('RealField layer', () => {
  const {
    rnPi,
    rnLog2,
    rnSqrt,
    rnExp,
    rnFromBigInt,
    rnFromRational,
    rnMul,
    rnDiv,
    rnAdd,
    rnSub,
    rnOne,
    rnZero,
    rnSetPrec,
    rnCmp,
    rnEq,
    rnScalePow2,
    exactFromDouble,
  } = _mp;

  test('pi, log 2, sqrt and exp match mpmath to 200 bits', () => {
    // mpmath, mp.prec = 400, printed with mp.nstr(v, 60); a 200-bit RealField
    // prints floor(199*log10(2)) = 59 digits.
    expect(rnPi(200).toString()).toBe(
      '3.1415926535897932384626433832795028841971693993751058209749'
    );
    expect(rnLog2(200).toString()).toBe(
      '0.69314718055994530941723212145817656807550013436025525412068'
    );
    expect(rnSqrt(rnFromBigInt(2n, 200), 200).toString()).toBe(
      '1.4142135623730950488016887242096980785696718753769480731767'
    );
    expect(rnExp(rnOne(200), 200).toString()).toBe(
      '2.7182818284590452353602874713526624977572470936999595749670'
    );
    expect(rnExp(rnFromBigInt(-2n, 200), 200).toString()).toBe(
      '0.13533528323661269189399949497248440340763154590957588146816'
    );
    expect(rnExp(rnFromBigInt(1000n, 200), 200).toString()).toBe(
      '1.9700711140170469938888793522433231253169379853238457899528e434'
    );
    expect(rnSqrt(rnMul(rnFromBigInt(2n, 200), rnPi(200), 200), 200).toString()).toBe(
      '2.5066282746310005024157652848110452530069867406099383166299'
    );
    // exp of a large negative argument, evaluated with the argument itself
    // given exactly (mpmath at 400 bits, 150 digits).
    expect(rnExp(rnFromRational(Rational.from('-12345.678'), 700), 500).toString()).toStartWith(
      '2.18861436617239831946402056504514994069104544836380259974004317381196371090648392372027889334554500675072042063720865640876364356173687623476685948502e-5362'
    );
  });

  test('+ - * / are correctly rounded (round to nearest, ties to even)', () => {
    // Compare against the exact rational result, rounded once.
    const toRat = (x: { s: number; m: bigint; e: number }): Rational => {
      const mag = new Rational(x.m, 1n);
      const scale =
        x.e >= 0 ? new Rational(1n << BigInt(x.e), 1n) : new Rational(1n, 1n << BigInt(-x.e));
      const v = mag.mul(scale);
      return x.s < 0 ? v.neg() : v;
    };
    let seed = 12345;
    const rnd = (): number => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed;
    };
    let checks = 0;
    for (let t = 0; t < 4000; t++) {
      const p = 8 + (rnd() % 60);
      const a = new Rational(BigInt(rnd() - 1073741824), BigInt(1 + (rnd() % 100000)));
      const b = new Rational(BigInt(rnd() - 1073741824), BigInt(1 + (rnd() % 100000)));
      const A = rnFromRational(a, p);
      const B = rnFromRational(b, p);
      const ea = toRat(A);
      const eb = toRat(B);
      expect(rnEq(rnAdd(A, B, p), rnFromRational(ea.add(eb), p))).toBe(true);
      expect(rnEq(rnSub(A, B, p), rnFromRational(ea.sub(eb), p))).toBe(true);
      expect(rnEq(rnMul(A, B, p), rnFromRational(ea.mul(eb), p))).toBe(true);
      checks += 3;
      if (!eb.isZero()) {
        expect(rnEq(rnDiv(A, B, p), rnFromRational(ea.div(eb), p))).toBe(true);
        checks += 1;
      }
      // Ordering agrees with the exact rational ordering.
      expect(rnCmp(A, B)).toBe(ea.cmp(eb));
    }
    expect(checks).toBeGreaterThan(15000);
  });

  test('addition of a far smaller operand rounds back to the larger one', () => {
    // The shortcut that keeps `1 + exp(-pi^2 * 2e6)` cheap must agree with the
    // exact computation right at the boundary.
    const p = 53;
    const one = rnOne(p);
    const pow2 = (k: number) => rnScalePow2(rnOne(p), k);
    // Numbers in [1, 2) have ulp 2^-52, so half an ulp above 1 is 2^-53.
    expect(rnEq(rnAdd(one, pow2(-54), p), one)).toBe(true); // below half an ulp
    expect(rnEq(rnAdd(one, pow2(-53), p), one)).toBe(true); // exact tie -> to even
    expect(rnEq(rnAdd(one, pow2(-52), p), one)).toBe(false); // representable
    // 1.5 * 2^-52 is above half an ulp and must round up.
    expect(rnEq(rnAdd(one, rnAdd(pow2(-53), pow2(-54), p), p), one)).toBe(false);
    // Numbers in [1/2, 1) have ulp 2^-53, so half a gap below 1 is 2^-54.
    expect(rnEq(rnSub(one, pow2(-55), p), one)).toBe(true); // below half a gap
    expect(rnEq(rnSub(one, pow2(-54), p), one)).toBe(true); // exact tie -> to even
    expect(rnEq(rnSub(one, pow2(-53), p), one)).toBe(false); // representable
    // A truly astronomical gap: this is the case the shortcut exists for.
    expect(rnEq(rnAdd(one, pow2(-28500000), p), one)).toBe(true);
    expect(rnEq(rnSub(one, pow2(-28500000), p), one)).toBe(true);
  });

  test('printing follows real_mpfr.pyx (truncate=True)', () => {
    // digits = max(2, floor((prec - 1) * log10 2)); scientific iff |E-1| >= 6.
    expect(rnFromBigInt(3n, 53).toString()).toBe('3.00000000000000');
    expect(rnFromBigInt(5n, 53).toString()).toBe('5.00000000000000');
    expect(rnZero(53).toString()).toBe('0.000000000000000');
    expect(rnFromRational(new Rational(61n, 3n), 53).toString()).toBe('20.3333333333333');
    expect(rnFromRational(new Rational(1n, 10n), 53).toString()).toBe('0.100000000000000');
    expect(rnFromBigInt(-7n, 53).toString()).toBe('-7.00000000000000');
    // |E - 1| >= 6 switches to scientific notation.
    expect(rnFromBigInt(10n ** 6n, 53).toString()).toBe('1.00000000000000e6');
    expect(rnFromBigInt(10n ** 5n, 53).toString()).toBe('100000.000000000');
    expect(rnFromRational(new Rational(1n, 10n ** 6n), 53).toString()).toBe('1.00000000000000e-6');
    expect(rnPi(53).toString()).toBe('3.14159265358979');
  });

  test('the printed decimal is the correctly rounded one, at many precisions', () => {
    // Independent check: parse the printed string back as an exact rational and
    // verify (a) it has exactly `digits` significant digits and (b) it is the
    // nearest such decimal to the exact binary value.
    const parseSage = (s: string): { value: Rational; digits: number } => {
      let str = s;
      let sign = 1n;
      if (str.startsWith('-')) {
        sign = -1n;
        str = str.slice(1);
      }
      let exp10 = 0;
      const eIdx = str.indexOf('e');
      if (eIdx >= 0) {
        exp10 = Number.parseInt(str.slice(eIdx + 1), 10);
        str = str.slice(0, eIdx);
      }
      const dot = str.indexOf('.');
      const intPart = dot < 0 ? str : str.slice(0, dot);
      const fracPart = dot < 0 ? '' : str.slice(dot + 1);
      const digitsStr = (intPart + fracPart).replace(/^0+/, '') || '0';
      const scale = exp10 - fracPart.length;
      const mant = new Rational(sign * BigInt(intPart + fracPart), 1n);
      const p10 =
        scale >= 0
          ? new Rational(10n ** BigInt(scale), 1n)
          : new Rational(1n, 10n ** BigInt(-scale));
      return { value: mant.mul(p10), digits: digitsStr.length };
    };
    const exactOf = (x: { s: number; m: bigint; e: number }): Rational => {
      const scale =
        x.e >= 0 ? new Rational(1n << BigInt(x.e), 1n) : new Rational(1n, 1n << BigInt(-x.e));
      const v = new Rational(x.m, 1n).mul(scale);
      return x.s < 0 ? v.neg() : v;
    };
    let seed = 987654321;
    const rnd = (): number => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed;
    };
    for (let t = 0; t < 600; t++) {
      const p = 10 + (rnd() % 120);
      const num = BigInt(rnd() - 1073741824) * BigInt(1 + (rnd() % 1000));
      if (num === 0n) continue;
      const den = BigInt(1 + (rnd() % 1000000));
      const shift = (rnd() % 60) - 30;
      let q = new Rational(num, den);
      q =
        shift >= 0
          ? q.mul(new Rational(1n << BigInt(shift), 1n))
          : q.div(new Rational(1n << BigInt(-shift), 1n));
      const x = rnFromRational(q, p);
      const { value } = parseSage(x.toString());
      const wantDigits = Math.max(2, Math.floor((p - 1) * Math.LN2 * Math.LOG10E));
      // The printed decimal is within half an ulp (of its `wantDigits`-th
      // significant digit) of the exact binary value, i.e. it is the correctly
      // rounded `wantDigits`-digit decimal.
      const exact = exactOf(x);
      const err = value.sub(exact).abs();
      const ten = new Rational(10n, 1n);
      const magnitude = exact.abs();
      let ulp = new Rational(1n, 1n);
      while (ulp.mul(ten).cmp(magnitude) <= 0) ulp = ulp.mul(ten);
      while (ulp.cmp(magnitude) > 0) ulp = ulp.div(ten);
      for (let i = 1; i < wantDigits; i++) ulp = ulp.div(ten);
      expect(err.mul(new Rational(2n, 1n)).cmp(ulp)).toBeLessThanOrEqual(0);
    }
  });

  test('printing agrees with V8 toPrecision(15) at 53 bits', () => {
    // An entirely independent decimal conversion (V8's dtoa) for the default
    // RealField(53), whose truncate=True width is exactly 15 digits.
    const parse = (s: string): Rational => {
      let str = s;
      let sign = 1n;
      if (str.startsWith('-')) {
        sign = -1n;
        str = str.slice(1);
      }
      let exp10 = 0;
      const eIdx = str.indexOf('e');
      if (eIdx >= 0) {
        exp10 = Number.parseInt(str.slice(eIdx + 1), 10);
        str = str.slice(0, eIdx);
      }
      const dot = str.indexOf('.');
      const intPart = dot < 0 ? str : str.slice(0, dot);
      const fracPart = dot < 0 ? '' : str.slice(dot + 1);
      const scale = exp10 - fracPart.length;
      const mant = new Rational(sign * BigInt(intPart + fracPart || '0'), 1n);
      const p10 =
        scale >= 0
          ? new Rational(10n ** BigInt(scale), 1n)
          : new Rational(1n, 10n ** BigInt(-scale));
      return mant.mul(p10);
    };
    let seed = 24681357;
    const rnd = (): number => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed;
    };
    let compared = 0;
    for (let t = 0; t < 3000; t++) {
      const d = ((rnd() - 1073741824) / (1 + (rnd() % 100000))) * 2 ** ((rnd() % 60) - 30);
      if (d === 0 || !Number.isFinite(d)) continue;
      const mine = rnFromRational(exactFromDouble(d), 53).toString();
      const js = d.toPrecision(15);
      expect(parse(mine).eq(parse(js))).toBe(true);
      compared++;
    }
    expect(compared).toBeGreaterThan(2900);
  });

  test('round() and floor() follow real_mpfr.pyx:3034 and :3058', () => {
    // RR(0.49).round() == 0, RR(0.5).round() == 1, RR(-0.5).round() == -1
    expect(rnFromRational(new Rational(49n, 100n), 53).round()).toBe(0n);
    expect(rnFromRational(new Rational(1n, 2n), 53).round()).toBe(1n);
    expect(rnFromRational(new Rational(-49n, 100n), 53).round()).toBe(0n);
    expect(rnFromRational(new Rational(-1n, 2n), 53).round()).toBe(-1n);
    expect(rnFromRational(new Rational(5n, 2n), 53).round()).toBe(3n);
    expect(rnFromRational(new Rational(299n, 100n), 53).floor()).toBe(2n);
    expect(rnFromRational(new Rational(-299n, 100n), 53).floor()).toBe(-3n);
    expect(rnFromBigInt(7n, 53).floor()).toBe(7n);
  });

  test('exactFromDouble reproduces the exact value of a RealField(53) element', () => {
    expect(exactFromDouble(0.5).eq(new Rational(1n, 2n))).toBe(true);
    expect(exactFromDouble(1000).eq(new Rational(1000n, 1n))).toBe(true);
    expect(exactFromDouble(0).isZero()).toBe(true);
    // 1.7 is not a dyadic rational: the double is the exact value below.
    expect(exactFromDouble(1.7).toString()).toBe('7656119366529843/4503599627370496');
    // Round-tripping through a 53-bit RealField must be the identity (going
    // through Rational.toNumber() would overflow for the extreme exponents).
    for (const x of [0.1, -3.25, 1e-300, 1e300, 12345.6789, 2 ** -1074, 2 ** 1023]) {
      expect(rnFromRational(exactFromDouble(x), 53).toNumber()).toBe(x);
    }
  });

  test('setPrec is exact when widening and correctly rounded when narrowing', () => {
    const pi200 = rnPi(200);
    expect(rnEq(rnSetPrec(rnSetPrec(pi200, 400), 200), pi200)).toBe(true);
    expect(rnSetPrec(pi200, 53).toString()).toBe('3.14159265358979');
    expect(rnEq(rnSetPrec(rnZero(53), 100), rnZero(100))).toBe(true);
  });
});

describe('_iter_vectors (discrete_gaussian_lattice.py:78-118)', () => {
  test('the docstring example: coordinate 0 varies fastest', () => {
    // `list(_iter_vectors(2, -1, 2))`, verified against SageMath 10.3.
    expect([..._iter_vectors(2n, -1n, 2n)]).toEqual([
      [-1n, -1n],
      [0n, -1n],
      [1n, -1n],
      [-1n, 0n],
      [0n, 0n],
      [1n, 0n],
      [-1n, 1n],
      [0n, 1n],
      [1n, 1n],
    ]);
  });

  test('length 1 and the full count for length 3', () => {
    expect([..._iter_vectors(1n, 0n, 3n)]).toEqual([[0n], [1n], [2n]]);
    expect([..._iter_vectors(3n, -1n, 2n)].length).toBe(27);
  });

  test('rejects lower >= upper and n <= 0', () => {
    // The vendored reference (`discrete_gaussian_lattice.py:100,102`) spells
    // these lower-case; installed SageMath 10.3 capitalises them.  Following
    // the vendored source per CLAUDE.md.
    expect(() => [..._iter_vectors(2n, 2n, 2n)]).toThrow('expected lower < upper, but got 2 >= 2');
    expect(() => [..._iter_vectors(0n, -1n, 2n)]).toThrow('expected n>0 but got 0 <= 0');
  });
});
