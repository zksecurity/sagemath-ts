/**
 * Unit tests for sage/crypto/lwe
 *
 * Tests for Learning with Errors (LWE) oracle implementation.
 */
import { describe, expect, test } from 'bun:test';
import { Zmod } from '../rings/finite_rings/integer_mod_ring.js';
import { PolynomialRingConstructor } from '../rings/polynomial/polynomial_ring.js';
import {
  DiscreteGaussianDistributionPolynomialSampler,
  LWE,
  LWEVector,
  LindnerPeikert,
  Regev,
  RingLWE,
  RingLWEConverter,
  RingLindnerPeikert,
  UniformPolynomialSampler,
  UniformSampler,
  balance_sample,
  samples,
} from './lwe.js';

describe('UniformSampler', () => {
  test('samples within bounds', () => {
    const sampler = new UniformSampler(-5n, 5n);

    for (let i = 0; i < 100; i++) {
      const sample = sampler.call();
      expect(sample >= -5n).toBe(true);
      expect(sample <= 5n).toBe(true);
    }
  });

  test('constructs with bigint bounds', () => {
    const sampler = new UniformSampler(-100n, 100n);
    expect(sampler.lower_bound).toBe(-100n);
    expect(sampler.upper_bound).toBe(100n);
  });

  test('throws if lower > upper', () => {
    expect(() => new UniformSampler(10n, 5n)).toThrow('lower bound must be <= upper bound');
  });

  test('repr', () => {
    const sampler = new UniformSampler(-2n, 2n);
    expect(sampler.repr()).toBe('UniformSampler(-2, 2)');
  });

  test('samples from single-element range', () => {
    const sampler = new UniformSampler(42n, 42n);
    for (let i = 0; i < 10; i++) {
      expect(sampler.call()).toBe(42n);
    }
  });

  test('samples from non-negative range', () => {
    const sampler = new UniformSampler(0n, 10n);
    for (let i = 0; i < 50; i++) {
      const sample = sampler.call();
      expect(sample >= 0n).toBe(true);
      expect(sample <= 10n).toBe(true);
    }
  });
});

describe('LWEVector', () => {
  test('dot product computation', () => {
    const ring = Zmod(17n);
    const v1 = new LWEVector([ring.__call__(3n), ring.__call__(5n)], ring);
    const v2 = new LWEVector([ring.__call__(2n), ring.__call__(4n)], ring);

    // 3*2 + 5*4 = 6 + 20 = 26 = 9 mod 17
    const dot = v1.dotProduct(v2);
    expect(dot.value).toBe(9n);
  });

  test('toTuple returns bigint array', () => {
    const ring = Zmod(11n);
    const v = new LWEVector([ring.__call__(3n), ring.__call__(7n), ring.__call__(2n)], ring);
    expect(v.toTuple()).toEqual([3n, 7n, 2n]);
  });

  test('length property', () => {
    const ring = Zmod(13n);
    const v = new LWEVector(
      [ring.__call__(1n), ring.__call__(2n), ring.__call__(3n), ring.__call__(4n)],
      ring
    );
    expect(v.length).toBe(4);
  });

  test('throws on mismatched lengths for dot product', () => {
    const ring = Zmod(17n);
    const v1 = new LWEVector([ring.__call__(1n), ring.__call__(2n)], ring);
    const v2 = new LWEVector([ring.__call__(1n), ring.__call__(2n), ring.__call__(3n)], ring);

    expect(() => v1.dotProduct(v2)).toThrow('Vectors must have the same length');
  });
});

describe('LWE', () => {
  test('generates samples of correct structure', () => {
    const D = new UniformSampler(-3n, 3n);
    const lwe = new LWE(5n, 101n, D);

    const [a, c] = lwe.call();

    expect(a.length).toBe(5);
    expect(a.ring.modulus).toBe(101n);
    expect(c.parent.modulus).toBe(101n);
  });

  test('knowing secret allows recovery of error', () => {
    const D = new UniformSampler(-2n, 2n);
    const lwe = new LWE(10n, 1009n, D);
    const secret = lwe.secret;

    // Generate a sample and verify error recovery
    const [a, c] = lwe.call();

    // c = <a, s> + e, so e = c - <a, s>
    const as = a.dotProduct(secret);
    const recovered_e = c.sub(as);

    // The error should be small (within the distribution bounds)
    // Balance to get the actual error value
    const q = 1009n;
    const q2 = q / 2n;
    const e_val = recovered_e.value <= q2 ? recovered_e.value : recovered_e.value - q;

    expect(e_val >= -2n).toBe(true);
    expect(e_val <= 2n).toBe(true);
  });

  test('multiple samples verify error recovery', () => {
    const D = new UniformSampler(-5n, 5n);
    const lwe = new LWE(8n, 509n, D);
    const secret = lwe.secret;
    const q = 509n;
    const q2 = q / 2n;

    for (let i = 0; i < 20; i++) {
      const [a, c] = lwe.call();
      const as = a.dotProduct(secret);
      const recovered_e = c.sub(as);

      const e_val = recovered_e.value <= q2 ? recovered_e.value : recovered_e.value - q;
      expect(e_val >= -5n).toBe(true);
      expect(e_val <= 5n).toBe(true);
    }
  });

  test('uniform secret distribution', () => {
    const D = new UniformSampler(-2n, 2n);
    const lwe = new LWE(5n, 101n, D, 'uniform');

    // Secret should be a vector of length n with elements in Z_q
    expect(lwe.secret.length).toBe(5);
    for (let i = 0; i < 5; i++) {
      const val = lwe.secret.get(i).value;
      expect(val >= 0n).toBe(true);
      expect(val < 101n).toBe(true);
    }
  });

  test('noise secret distribution', () => {
    const D = new UniformSampler(-3n, 3n);
    const lwe = new LWE(5n, 101n, D, 'noise');

    // Secret should be drawn from the noise distribution
    // Since D is UniformSampler(-3, 3), secret entries should be small when balanced
    const q = 101n;
    const q2 = q / 2n;

    for (let i = 0; i < 5; i++) {
      const val = lwe.secret.get(i).value;
      const balanced = val <= q2 ? val : val - q;
      expect(balanced >= -3n).toBe(true);
      expect(balanced <= 3n).toBe(true);
    }
  });

  test('bounded secret distribution', () => {
    const D = new UniformSampler(-2n, 2n);
    const lwe = new LWE(5n, 101n, D, [-5n, 5n]);

    // Secret entries should be drawn uniformly from [-5, 5]
    const q = 101n;
    const q2 = q / 2n;

    for (let i = 0; i < 5; i++) {
      const val = lwe.secret.get(i).value;
      const balanced = val <= q2 ? val : val - q;
      expect(balanced >= -5n).toBe(true);
      expect(balanced <= 5n).toBe(true);
    }
  });

  test('sample count limit (m parameter)', () => {
    const D = new UniformSampler(-2n, 2n);
    const lwe = new LWE(5n, 101n, D, 'uniform', 3n);

    // Should be able to get 3 samples
    lwe.call();
    lwe.call();
    lwe.call();

    // 4th sample should throw
    expect(() => lwe.call()).toThrow('Number of available samples exhausted');
  });

  test('repr string', () => {
    const D = new UniformSampler(-2n, 2n);
    const lwe = new LWE(10n, 101n, D, 'uniform', 100n);

    const repr = lwe.repr();
    expect(repr).toContain('LWE');
    expect(repr).toContain('10');
    expect(repr).toContain('101');
    expect(repr).toContain('uniform');
  });

  test('samples method generates multiple samples', () => {
    const D = new UniformSampler(-2n, 2n);
    const lwe = new LWE(5n, 101n, D);

    const sampleList = lwe.samples(10n);
    expect(sampleList.length).toBe(10);

    for (const [a, c] of sampleList) {
      expect(a.length).toBe(5);
      expect(a.ring.modulus).toBe(101n);
      expect(c.parent.modulus).toBe(101n);
    }
  });
});

describe('Regev', () => {
  test('creates oracle with Regev parameters', () => {
    const regev = new Regev(20n);

    // q should be the smallest prime >= n^2 = 400
    expect(regev.K.modulus).toBe(401n);
    expect(regev.n).toBe(20);
  });

  test('generates valid samples', () => {
    const regev = new Regev(15n);
    const [a, c] = regev.call();

    expect(a.length).toBe(15);
    // q should be smallest prime >= 225
    expect(a.ring.modulus).toBe(227n);
  });

  test('error recovery works with Regev parameters', () => {
    const regev = new Regev(10n);
    const secret = regev.secret;
    const q = regev.K.modulus;
    const q2 = q / 2n;

    // Generate samples and verify errors are small
    for (let i = 0; i < 10; i++) {
      const [a, c] = regev.call();
      const as = a.dotProduct(secret);
      const recovered_e = c.sub(as);
      const e_val = recovered_e.value <= q2 ? recovered_e.value : recovered_e.value - q;

      // Error should be bounded by the distribution
      // For Regev, sigma ~= q/(sqrt(n)*log(n)^2*sqrt(2*pi))
      // With our uniform approximation, errors should be reasonably small
      expect(e_val >= -q / 4n).toBe(true);
      expect(e_val <= q / 4n).toBe(true);
    }
  });

  test('supports noise secret distribution', () => {
    const regev = new Regev(10n, 'noise');
    expect(regev.secret_dist).toBe('noise');
  });

  test('supports sample limit', () => {
    const regev = new Regev(10n, 'uniform', 5n);
    expect(regev.m).toBe(5);

    // Should be able to get 5 samples
    for (let i = 0; i < 5; i++) {
      regev.call();
    }

    expect(() => regev.call()).toThrow('Number of available samples exhausted');
  });
});

describe('LindnerPeikert', () => {
  test('creates oracle with LP parameters', () => {
    const lp = new LindnerPeikert(20n);

    // With LP parameters, m should be 2*n + 128 = 168
    expect(lp.m).toBe(168);
    expect(lp.n).toBe(20);

    // Secret distribution should be 'noise' per [LP2011]
    expect(lp.secret_dist).toBe('noise');
  });

  test('generates valid samples', () => {
    const lp = new LindnerPeikert(15n);
    const [a, c] = lp.call();

    expect(a.length).toBe(15);
    expect(c.parent.modulus).toBe(lp.K.modulus);
  });

  test('custom delta parameter', () => {
    const lp1 = new LindnerPeikert(20n, 0.01);
    const lp2 = new LindnerPeikert(20n, 0.001);

    // Different delta values might result in different q values
    // Both should still create valid oracles
    expect(lp1.n).toBe(20);
    expect(lp2.n).toBe(20);
  });

  test('custom m parameter', () => {
    const lp = new LindnerPeikert(20n, 0.01, 50n);
    expect(lp.m).toBe(50);
  });
});

describe('samples function', () => {
  test('generates samples using string oracle name', () => {
    const S = samples(10n, 15n, 'Regev');
    expect(S.length).toBe(10);
  });

  test('generates samples using class reference', () => {
    const S = samples(10n, 15n, Regev);
    expect(S.length).toBe(10);
  });

  test('generates samples using instance', () => {
    const regev = new Regev(15n);
    const S = samples(10n, 15n, regev);
    expect(S.length).toBe(10);
  });

  test('throws on dimension mismatch with instance', () => {
    const regev = new Regev(15n);
    expect(() => samples(10n, 20n, regev)).toThrow(
      'Passed LWE instance has n=15, but n=20 was passed to this function'
    );
  });

  test('balanced samples have correct range', () => {
    const S = samples(20n, 10n, 'Regev', { balanced: true });

    for (const [a, c] of S) {
      // a should be an array of bigints in balanced form
      expect(Array.isArray(a)).toBe(true);
      const aArr = a as bigint[];

      // All values should be balanced (between -q/2 and q/2)
      // For Regev with n=10, q = next_prime(100) = 101
      for (const val of aArr) {
        expect(val >= -51n).toBe(true);
        expect(val <= 50n).toBe(true);
      }

      // c should be a bigint in balanced form
      expect(typeof c).toBe('bigint');
      expect(c >= -51n).toBe(true);
      expect(c <= 50n).toBe(true);
    }
  });

  test('LindnerPeikert oracle by name', () => {
    const S = samples(5n, 20n, 'LindnerPeikert');
    expect(S.length).toBe(5);
  });

  test('seed makes the output reproducible', () => {
    // Sage: `if seed is not None: set_random_seed(seed)`
    const show = (S: ReturnType<typeof samples>) =>
      S.map(([a, c]) => {
        const av = a instanceof LWEVector ? a.toTuple() : (a as bigint[]);
        const cv = Array.isArray(c) ? c.map(String).join(',') : String(c);
        return `${av.map(String).join(',')}|${cv}`;
      }).join(';');

    const s1 = samples(3n, 20n, 'Regev', { seed: 1337 });
    const s2 = samples(3n, 20n, 'Regev', { seed: 1337 });
    expect(show(s1)).toBe(show(s2));

    const s3 = samples(3n, 20n, 'Regev', { seed: 1338 });
    expect(show(s3)).not.toBe(show(s1));
  });

  test('m is forwarded to the oracle constructor', () => {
    // Sage: `lwe = lwe(n, m=m, **kwds)`.  LindnerPeikert's default m is
    // 2*n+128 = 168, so without forwarding, asking for 200 samples would
    // exhaust the oracle.
    const S = samples(200n, 20n, 'LindnerPeikert');
    expect(S.length).toBe(200);

    // An explicitly constructed instance keeps its own (default) limit.
    const lp = new LindnerPeikert(20n);
    expect(lp.m).toBe(168);
    expect(() => samples(200n, 20n, lp)).toThrow('Number of available samples exhausted.');
  });

  test('kwds are passed through to the oracle constructor', () => {
    const S = samples(4n, 10n, 'Regev', { secret_dist: 'noise' });
    expect(S.length).toBe(4);
  });

  test('accepts Ring-LWE oracles, as Sage does', () => {
    // sage: D = DiscreteGaussianDistributionPolynomialSampler(ZZ['x'], 8, 5)
    // sage: rlwe = RingLWE(20, 257, D)
    // sage: for s in samples(10, 8, rlwe): ...
    const K = Zmod(257n);
    const [R] = PolynomialRingConstructor(K, 'x');
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 5);
    const rlwe = new RingLWE(20n, 257n, D);
    const S = samples(10n, 8n, rlwe);
    expect(S.length).toBe(10);
    for (const [a, c] of S) {
      expect((a as bigint[]).length).toBe(8);
      expect((c as bigint[]).length).toBe(8);
    }

    // balance_sample handles the vector right-hand side too.
    for (const s of S) {
      const [ba, bc] = balance_sample(s, 257n);
      expect(ba.every((v) => v >= -128n && v <= 128n)).toBe(true);
      expect((bc as bigint[]).every((v) => v >= -128n && v <= 128n)).toBe(true);
    }

    // and Ring-LWE converters
    const conv = new RingLWEConverter(rlwe);
    expect(samples(4n, 8n, conv).length).toBe(4);
  });
});

describe('balance_sample', () => {
  test('balances positive values', () => {
    const ring = Zmod(101n);
    const a = new LWEVector([ring.__call__(80n), ring.__call__(95n), ring.__call__(100n)], ring);
    const c = ring.__call__(90n);

    const [balA, balC] = balance_sample([a, c], 101n);

    // 80 > 50, so balanced is 80 - 101 = -21
    expect(balA[0]).toBe(-21n);
    // 95 > 50, so balanced is 95 - 101 = -6
    expect(balA[1]).toBe(-6n);
    // 100 > 50, so balanced is 100 - 101 = -1
    expect(balA[2]).toBe(-1n);
    // 90 > 50, so balanced is 90 - 101 = -11
    expect(balC).toBe(-11n);
  });

  test('keeps small positive values unchanged', () => {
    const ring = Zmod(101n);
    const a = new LWEVector([ring.__call__(10n), ring.__call__(25n), ring.__call__(50n)], ring);
    const c = ring.__call__(40n);

    const [balA, balC] = balance_sample([a, c], 101n);

    expect(balA[0]).toBe(10n);
    expect(balA[1]).toBe(25n);
    expect(balA[2]).toBe(50n);
    expect(balC).toBe(40n);
  });

  test('infers q from vector ring', () => {
    const ring = Zmod(67n);
    const a = new LWEVector([ring.__call__(60n), ring.__call__(5n)], ring);
    const c = ring.__call__(35n);

    const [balA, balC] = balance_sample([a, c]);

    // 60 > 33, so balanced is 60 - 67 = -7
    expect(balA[0]).toBe(-7n);
    // 5 <= 33, stays as 5
    expect(balA[1]).toBe(5n);
    // 35 > 33, so balanced is 35 - 67 = -32
    expect(balC).toBe(-32n);
  });
});

describe('LWE security verification', () => {
  test('samples appear uniformly random without secret', () => {
    // Without the secret, samples should look random
    // This is a statistical test - we check that the values appear distributed
    const D = new UniformSampler(-2n, 2n);
    const lwe = new LWE(10n, 101n, D);

    const cValues: bigint[] = [];
    for (let i = 0; i < 100; i++) {
      const [_a, c] = lwe.call();
      cValues.push(c.value);
    }

    // Check that c values span a reasonable range (not concentrated)
    const uniqueValues = new Set(cValues);
    // With 100 samples over Z_101, we expect some unique values
    // Note: This is inherently statistical and depends on the random seed
    // With small error distribution [-2,2] and n=10, we expect reasonable spread
    expect(uniqueValues.size).toBeGreaterThan(5);
  });

  test('knowing secret reveals error structure', () => {
    const D = new UniformSampler(-3n, 3n);
    const lwe = new LWE(8n, 509n, D);
    const secret = lwe.secret;
    const q = 509n;
    const q2 = q / 2n;

    const errors: bigint[] = [];
    for (let i = 0; i < 50; i++) {
      const [a, c] = lwe.call();
      const as = a.dotProduct(secret);
      const recovered_e = c.sub(as);
      const e_val = recovered_e.value <= q2 ? recovered_e.value : recovered_e.value - q;
      errors.push(e_val);
    }

    // All errors should be within the distribution bounds
    for (const e of errors) {
      expect(e >= -3n).toBe(true);
      expect(e <= 3n).toBe(true);
    }

    // Errors should have some variety (not all the same)
    const uniqueErrors = new Set(errors);
    expect(uniqueErrors.size).toBeGreaterThan(1);
  });
});

describe('edge cases', () => {
  test('LWE with dimension 1', () => {
    const D = new UniformSampler(-1n, 1n);
    const lwe = new LWE(1n, 17n, D);

    const [a, c] = lwe.call();
    expect(a.length).toBe(1);
    expect(c.parent.modulus).toBe(17n);
  });

  test('LWE with small modulus', () => {
    const D = new UniformSampler(0n, 0n);
    const lwe = new LWE(3n, 5n, D);

    const [a, c] = lwe.call();
    expect(a.length).toBe(3);
    expect(a.ring.modulus).toBe(5n);
  });

  test('Regev with small n', () => {
    const regev = new Regev(2n);
    // q should be next_prime(4) = 5
    expect(regev.K.modulus).toBe(5n);
  });

  test('zero error distribution', () => {
    // With zero error, c = <a, s> exactly
    const D = new UniformSampler(0n, 0n);
    const lwe = new LWE(5n, 101n, D);
    const secret = lwe.secret;

    const [a, c] = lwe.call();
    const as = a.dotProduct(secret);

    // c should equal <a, s> exactly
    expect(c.value).toBe(as.value);
  });
});

describe('UniformPolynomialSampler', () => {
  test('samples polynomial with correct degree', () => {
    const [R] = PolynomialRingConstructor(Zmod(101n), 'x');
    const sampler = new UniformPolynomialSampler(R, 8n, -2n, 2n);

    const poly = sampler.call();

    // Degree should be at most n-1 = 7
    expect(poly.degree()).toBeLessThanOrEqual(7);
  });

  test('samples coefficients within bounds', () => {
    const [R] = PolynomialRingConstructor(Zmod(101n), 'x');
    const sampler = new UniformPolynomialSampler(R, 8n, -3n, 3n);

    for (let i = 0; i < 10; i++) {
      const poly = sampler.call();
      for (let j = 0; j < 8; j++) {
        const coeff = poly.getCoeff(j);
        // In Z_101, values -3..3 map to 0,1,2,3 and 98,99,100
        const val = coeff.value;
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

  test('throws on invalid bounds', () => {
    const [R] = PolynomialRingConstructor(Zmod(101n), 'x');
    expect(() => new UniformPolynomialSampler(R, 8n, 10n, 5n)).toThrow(
      'lower bound must be <= upper bound'
    );
  });
});

describe('DiscreteGaussianDistributionPolynomialSampler', () => {
  test('samples polynomial with correct degree', () => {
    const [R] = PolynomialRingConstructor(Zmod(101n), 'x');
    const sampler = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);

    const poly = sampler.call();

    // Degree should be at most n-1 = 7
    expect(poly.degree()).toBeLessThanOrEqual(7);
  });

  test('repr matches SageMath format', () => {
    const [R] = PolynomialRingConstructor(Zmod(101n), 'x');
    const sampler = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
    expect(sampler.repr()).toContain('Discrete Gaussian sampler for polynomials of degree < 8');
    expect(sampler.repr()).toContain('σ=3.000000');
  });

  test('samples coefficients are small (near center)', () => {
    const [R] = PolynomialRingConstructor(Zmod(1009n), 'x');
    const sampler = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);

    // Collect samples and verify they're mostly small
    let smallCount = 0;
    const total = 50;

    for (let i = 0; i < total; i++) {
      const poly = sampler.call();
      for (let j = 0; j < 8; j++) {
        const coeff = poly.getCoeff(j);
        const val = coeff.value;
        // Balance to get signed value
        const balanced = val <= 504n ? val : val - 1009n;
        // With sigma=3, most values should be within 3*sigma = 9
        if (balanced >= -15n && balanced <= 15n) {
          smallCount++;
        }
      }
    }

    // At least 80% of coefficients should be within 5*sigma
    expect(smallCount).toBeGreaterThan(0.8 * total * 8);
  });
});

describe('RingLWE', () => {
  test('creates oracle with correct parameters', () => {
    const K = Zmod(257n);
    const [R] = PolynomialRingConstructor(K, 'x');
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
    const ringlwe = new RingLWE(16n, 257n, D);

    expect(ringlwe.N).toBe(16);
    expect(ringlwe.n).toBe(8); // phi(16) = 8
    expect(ringlwe.q).toBe(257n);
  });

  test('generates samples with correct dimensions', () => {
    const K = Zmod(257n);
    const [R] = PolynomialRingConstructor(K, 'x');
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
    const ringlwe = new RingLWE(16n, 257n, D);

    const [a, c] = ringlwe.call();

    expect(a.length).toBe(8);
    expect(c.length).toBe(8);
  });

  test('samples are in valid range', () => {
    const K = Zmod(257n);
    const [R] = PolynomialRingConstructor(K, 'x');
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
    const ringlwe = new RingLWE(16n, 257n, D);

    const [a, c] = ringlwe.call();

    for (let i = 0; i < 8; i++) {
      expect(a[i]! >= 0n).toBe(true);
      expect(a[i]! < 257n).toBe(true);
      expect(c[i]! >= 0n).toBe(true);
      expect(c[i]! < 257n).toBe(true);
    }
  });

  test('noise distribution dimension check', () => {
    const K = Zmod(257n);
    const [R] = PolynomialRingConstructor(K, 'x');
    // phi(16) = 8, but we create sampler with n=4 (wrong dimension)
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 4, 3.0);

    expect(() => new RingLWE(16n, 257n, D)).toThrow('Noise distribution has dimensions 4 != 8');
  });

  test('sample count limit', () => {
    const K = Zmod(257n);
    const [R] = PolynomialRingConstructor(K, 'x');
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
    const ringlwe = new RingLWE(16n, 257n, D, null, 'uniform', 3n);

    // Should be able to get 3 samples
    ringlwe.call();
    ringlwe.call();
    ringlwe.call();

    // 4th sample should throw
    expect(() => ringlwe.call()).toThrow('Number of available samples exhausted');
  });

  test('different N values give correct n = phi(N)', () => {
    // Test various N values
    const testCases = [
      { N: 8n, expectedN: 4 }, // phi(8) = 4
      { N: 16n, expectedN: 8 }, // phi(16) = 8
      { N: 32n, expectedN: 16 }, // phi(32) = 16
      { N: 20n, expectedN: 8 }, // phi(20) = 8
    ];

    for (const { N, expectedN } of testCases) {
      const K = Zmod(257n);
      const [R] = PolynomialRingConstructor(K, 'x');
      const D = new DiscreteGaussianDistributionPolynomialSampler(R, expectedN, 3.0);
      const ringlwe = new RingLWE(N, 257n, D);
      expect(ringlwe.n).toBe(expectedN);
    }
  });
});

describe('RingLindnerPeikert', () => {
  test('creates oracle with correct parameters', () => {
    const rlp = new RingLindnerPeikert(16n);

    expect(rlp.N).toBe(16);
    expect(rlp.n).toBe(8); // phi(16) = 8
    // Default m = 3 * n = 24
    expect(rlp.m).toBe(24);
    // Secret distribution should be 'noise'
    expect(rlp.secret_dist).toBe('noise');
  });

  test('generates valid samples', () => {
    const rlp = new RingLindnerPeikert(16n);
    const [a, c] = rlp.call();

    expect(a.length).toBe(8);
    expect(c.length).toBe(8);
  });

  test('custom m parameter', () => {
    const rlp = new RingLindnerPeikert(16n, 0.01, 50n);
    expect(rlp.m).toBe(50);
  });

  test('repr contains RingLWE info', () => {
    const rlp = new RingLindnerPeikert(16n);
    const repr = rlp.repr();
    expect(repr).toContain('RingLWE');
    expect(repr).toContain('16');
    expect(repr).toContain('noise');
  });
});

describe('RingLWEConverter', () => {
  test('converts Ring-LWE to LWE samples', () => {
    const K = Zmod(257n);
    const [R] = PolynomialRingConstructor(K, 'x');
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
    const ringlwe = new RingLWE(16n, 257n, D);
    const converter = new RingLWEConverter(ringlwe);

    const [a, c] = converter.call();

    // a should be a vector of length n
    expect(a.length).toBe(8);
    // c should be a scalar
    expect(typeof c).toBe('bigint');
    expect(c >= 0n).toBe(true);
    expect(c < 257n).toBe(true);
  });

  test('generates n LWE samples per Ring-LWE sample', () => {
    const K = Zmod(257n);
    const [R] = PolynomialRingConstructor(K, 'x');
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
    const ringlwe = new RingLWE(16n, 257n, D);
    const converter = new RingLWEConverter(ringlwe);

    // Get n samples (should all come from the same Ring-LWE sample)
    const samples: Array<[bigint[], bigint]> = [];
    for (let i = 0; i < 8; i++) {
      samples.push(converter.call());
    }

    // All samples should have correct dimensions
    for (const [a, c] of samples) {
      expect(a.length).toBe(8);
      expect(typeof c).toBe('bigint');
    }
  });

  test('repr includes Ring-LWE oracle info', () => {
    const K = Zmod(257n);
    const [R] = PolynomialRingConstructor(K, 'x');
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
    const ringlwe = new RingLWE(16n, 257n, D);
    const converter = new RingLWEConverter(ringlwe);

    const repr = converter.repr();
    expect(repr).toContain('RingLWEConverter');
    expect(repr).toContain('RingLWE');
  });

  test('successive samples are x^i * a in R_q, not a signed rotation', () => {
    // Sage: r = vector((x**(self._i % self.n) * R_q(a.list())).list()), ...
    // For N = 20, Phi_20 = x^8 - x^6 + x^4 - x^2 + 1, so multiplication by x
    // is NOT the negacyclic shift used when Phi_N = x^n + 1.
    const q = 257n;
    const K = Zmod(q);
    const [R] = PolynomialRingConstructor(K, 'x');
    const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 5);
    const ringlwe = new RingLWE(20n, q, D);
    expect(String(ringlwe.poly)).toBe('x^8 + 256*x^6 + x^4 + 256*x^2 + 1'); // x^8-x^6+x^4-x^2+1
    const converter = new RingLWEConverter(ringlwe);

    // Phi_20 coefficients, ascending, reduced mod q
    const phi = [1n, 0n, -1n, 0n, 1n, 0n, -1n, 0n, 1n].map((v) => ((v % q) + q) % q);
    const mulX = (c: bigint[]): bigint[] => {
      const lead = c[7]!;
      const out = new Array<bigint>(8);
      out[0] = (((-lead * phi[0]!) % q) + q) % q;
      for (let i = 1; i < 8; i++) out[i] = (((c[i - 1]! - lead * phi[i]!) % q) + q) % q;
      return out;
    };

    const [a0] = converter.call();
    let expected = mulX(a0);
    for (let i = 1; i < 8; i++) {
      const [ai] = converter.call();
      expect(ai).toEqual(expected);
      expected = mulX(expected);
    }

    // The negacyclic shift differs from the true product (given a nonzero
    // leading coefficient, which the x^2 and x^6 terms of Phi_20 pick up).
    const negacyclic = [((-a0[7]! % q) + q) % q, ...a0.slice(0, 7)];
    expect(negacyclic).not.toEqual(mulX(a0));
  });
});
