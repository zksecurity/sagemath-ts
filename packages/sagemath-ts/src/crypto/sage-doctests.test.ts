/**
 * SageMath Doctests for crypto module (LWE, lattice)
 *
 * This file contains tests ported from SageMath doctests.
 * Reference files:
 * - sage/crypto/lwe.py
 * - sage/crypto/lattice.py
 *
 * Tests are organized by the original doctest source.
 */
import { describe, expect, test } from 'bun:test';
import {
  UniformSampler,
  LWE,
  LWEVector,
  Regev,
  LindnerPeikert,
  UniformNoiseLWE,
  samples,
  balance_sample,
} from './lwe.js';
import { Zmod } from '../rings/finite_rings/integer_mod_ring.js';
import { next_prime } from '../arith/misc.js';

/**
 * Tests ported from sage/crypto/lwe.py
 */
describe('lwe.py doctests', () => {
  /**
   * From lwe.py module docstring:
   * sage: from sage.crypto.lwe import samples
   * sage: S = samples(30, 20, 'Regev')
   * sage: len(S)
   * 30
   */
  test('lwe.py - samples function returns correct number of samples', () => {
    const S = samples(30, 20, 'Regev');
    expect(S.length).toBe(30);
  });

  /**
   * From lwe.py module docstring:
   * sage: S[0][0].parent(), S[0][1].parent()
   * (Vector space of dimension 20 over Ring of integers modulo 401,
   *  Ring of integers modulo 401)
   */
  test('lwe.py - Regev sample dimension and modulus', () => {
    const S = samples(30, 20, 'Regev');
    const [a, _c] = S[0]!;

    // For Regev with n=20, q = next_prime(400) = 401
    if (a instanceof LWEVector) {
      expect(a.length).toBe(20);
      expect(a.ring.modulus).toBe(401n);
    } else {
      // Balanced case - array of bigints
      expect(a.length).toBe(20);
    }
  });

  /**
   * From lwe.py module docstring:
   * sage: from sage.crypto.lwe import samples, LindnerPeikert
   * sage: S = samples(30, 20, LindnerPeikert)
   * sage: len(S)
   * 30
   */
  test('lwe.py - samples with LindnerPeikert class', () => {
    const S = samples(30, 20, LindnerPeikert);
    expect(S.length).toBe(30);
  });

  /**
   * From lwe.py module docstring:
   * sage: from sage.crypto.lwe import LindnerPeikert
   * sage: lwe = LindnerPeikert(20)
   * sage: S = samples(30, 20, lwe)
   * sage: len(S)
   * 30
   */
  test('lwe.py - samples with LindnerPeikert instance', () => {
    const lwe = new LindnerPeikert(20);
    const S = samples(30, 20, lwe);
    expect(S.length).toBe(30);
  });

  /**
   * From lwe.py module docstring - balanced representation:
   * sage: from sage.crypto.lwe import samples
   * sage: for s in samples(30, 20, 'Regev', balanced=True):
   * ....:     s1 = list(s[0]) + [s[1]]
   * ....:     assert all(-401//2 <= b <= 401//2 for b in s1)
   */
  test('lwe.py - balanced samples have values in [-q//2, q//2]', () => {
    const S = samples(30, 20, 'Regev', { balanced: true });
    const q = 401n;
    const halfQ = q / 2n;

    for (const s of S) {
      const [a, c] = s;
      // a should be bigint[] when balanced
      const aArr = a as bigint[];
      const s1 = [...aArr, c as bigint];

      for (const b of s1) {
        expect(b >= -halfQ).toBe(true);
        expect(b <= halfQ).toBe(true);
      }
    }
  });

  /**
   * From UniformSampler class:
   * sage: from sage.crypto.lwe import UniformSampler
   * sage: sampler = UniformSampler(-2, 2); sampler
   * UniformSampler(-2, 2)
   */
  test('lwe.py - UniformSampler repr', () => {
    const sampler = new UniformSampler(-2, 2);
    expect(sampler.repr()).toBe('UniformSampler(-2, 2)');
  });

  /**
   * From UniformSampler class:
   * sage: sampler() in range(-2, 3)
   * True
   */
  test('lwe.py - UniformSampler sample in range', () => {
    const sampler = new UniformSampler(-2, 2);
    for (let i = 0; i < 100; i++) {
      const sample = sampler.call();
      expect(sample >= -2n).toBe(true);
      expect(sample <= 2n).toBe(true);
    }
  });

  /**
   * From UniformSampler.__call__ docstring:
   * sage: from sage.crypto.lwe import UniformSampler
   * sage: sampler = UniformSampler(-12, 12)
   * sage: sampler() in range(-12, 13)
   * True
   */
  test('lwe.py - UniformSampler(-12, 12) sample in range', () => {
    const sampler = new UniformSampler(-12, 12);
    for (let i = 0; i < 100; i++) {
      const sample = sampler.call();
      expect(sample >= -12n).toBe(true);
      expect(sample <= 12n).toBe(true);
    }
  });

  /**
   * From LWE class:
   * sage: from sage.crypto.lwe import LWE
   * sage: lwe = LWE(n=20, q=next_prime(400), D=D); lwe
   * LWE(20, 401, ...)
   */
  test('lwe.py - LWE construction with n=20, q=401', () => {
    const q = next_prime(400n);
    expect(q).toBe(401n);

    const D = new UniformSampler(-3, 3);
    const lwe = new LWE(20, q, D);

    expect(lwe.n).toBe(20);
    expect(lwe.K.modulus).toBe(401n);
  });

  /**
   * From LWE.__call__ docstring:
   * sage: LWE(10, 401, DiscreteGaussianDistributionIntegerSampler(3))()[0].parent()
   * Vector space of dimension 10 over Ring of integers modulo 401
   *
   * We use UniformSampler instead of DiscreteGaussian
   */
  test('lwe.py - LWE(10, 401) sample dimension', () => {
    const D = new UniformSampler(-3, 3);
    const lwe = new LWE(10, 401n, D);
    const [a, c] = lwe.call();

    expect(a.length).toBe(10);
    expect(a.ring.modulus).toBe(401n);
    expect(c.parent.modulus).toBe(401n);
  });

  /**
   * From LWE class - sample limit (m parameter):
   * sage: lwe = LWE(n=20, q=next_prime(400), D=D, m=30)
   * sage: _ = [lwe() for _ in range(30)]
   * sage: lwe()  # 31
   * Traceback (most recent call last):
   * ...
   * IndexError: Number of available samples exhausted.
   */
  test('lwe.py - LWE sample limit exhaustion', () => {
    const D = new UniformSampler(-3, 3);
    const lwe = new LWE(20, 401n, D, 'uniform', 30);

    // Get 30 samples
    for (let i = 0; i < 30; i++) {
      lwe.call();
    }

    // 31st should throw
    expect(() => lwe.call()).toThrow('Number of available samples exhausted');
  });

  /**
   * From LWE class - secret distribution as tuple:
   * sage: lwe = LWE(n=20, q=next_prime(400), D=D, secret_dist=(-3, 3)); lwe
   * LWE(20, 401, ..., (-3, 3), None)
   */
  test('lwe.py - LWE with bounded secret distribution', () => {
    const D = new UniformSampler(-3, 3);
    const lwe = new LWE(20, 401n, D, [-3, 3]);

    const repr = lwe.repr();
    expect(repr).toContain('(-3, 3)');
    expect(lwe.n).toBe(20);

    // Secret should be bounded
    const q = 401n;
    const q2 = q / 2n;
    for (let i = 0; i < lwe.n; i++) {
      const val = lwe.secret.get(i).value;
      const balanced = val <= q2 ? val : val - q;
      expect(balanced >= -3n).toBe(true);
      expect(balanced <= 3n).toBe(true);
    }
  });

  /**
   * From Regev class:
   * sage: from sage.crypto.lwe import Regev
   * sage: Regev(n=20)
   * LWE(20, 401, Discrete Gaussian sampler..., 'uniform', None)
   */
  test('lwe.py - Regev(20) parameters', () => {
    const regev = new Regev(20);
    // q = next_prime(n^2) = next_prime(400) = 401
    expect(regev.K.modulus).toBe(401n);
    expect(regev.n).toBe(20);
    expect(regev.secret_dist).toBe('uniform');
    expect(regev.m).toBe(null);
  });

  /**
   * From LindnerPeikert class:
   * sage: from sage.crypto.lwe import LindnerPeikert
   * sage: LindnerPeikert(n=20)
   * LWE(20, 2053, ..., 'noise', 168)
   */
  test('lwe.py - LindnerPeikert(20) parameters', () => {
    const lp = new LindnerPeikert(20);
    expect(lp.n).toBe(20);
    expect(lp.secret_dist).toBe('noise');
    // m = 2*n + 128 = 168
    expect(lp.m).toBe(168);
  });

  /**
   * From UniformNoiseLWE class:
   * sage: from sage.crypto.lwe import UniformNoiseLWE
   * sage: UniformNoiseLWE(89)
   * LWE(89, 64311834871, UniformSampler(0, 6577), 'noise', 131)
   */
  test('lwe.py - UniformNoiseLWE(89) key instance', () => {
    const lwe = new UniformNoiseLWE(89);
    expect(lwe.n).toBe(89);
    expect(lwe.secret_dist).toBe('noise');
    // m should be n1 for 'key' instance
    expect(lwe.m).toBe(131);
  });

  /**
   * From UniformNoiseLWE class:
   * sage: UniformNoiseLWE(89, instance='encrypt')
   * LWE(131, 64311834871, UniformSampler(0, 11109), 'noise', 181)
   */
  test('lwe.py - UniformNoiseLWE(89, encrypt) instance', () => {
    const lwe = new UniformNoiseLWE(89, 'encrypt');
    // For 'encrypt', n = n1
    expect(lwe.n).toBe(131);
    expect(lwe.secret_dist).toBe('noise');
    // m should be n2 + l for 'encrypt' instance
    expect(lwe.m).toBe(181);
  });

  /**
   * From UniformNoiseLWE class - parameter too small:
   * sage: UniformNoiseLWE(88)  # n < 89
   * Traceback (most recent call last):
   * ...
   * TypeError: Parameter too small
   */
  test('lwe.py - UniformNoiseLWE parameter too small', () => {
    expect(() => new UniformNoiseLWE(88)).toThrow('Parameter too small');
  });

  /**
   * From samples function:
   * sage: from sage.crypto.lwe import samples, Regev
   * sage: samples(2, 20, Regev, seed=1337)
   * [((199, 388, ...), 15), ((365, 227, ...), 143)]
   *
   * We can't test exact values without seeded random, but we can verify structure
   */
  test('lwe.py - samples function structure', () => {
    const S = samples(2, 20, Regev);
    expect(S.length).toBe(2);

    for (const [a, _c] of S) {
      if (a instanceof LWEVector) {
        expect(a.length).toBe(20);
      } else {
        expect(a.length).toBe(20);
      }
    }
  });

  /**
   * From samples function with balanced:
   * sage: samples(2, 20, Regev, balanced=True, seed=1337)
   * [((199, -13, -64, ...), 15), ((-36, -174, ...), 143)]
   */
  test('lwe.py - samples balanced returns bigint arrays', () => {
    const S = samples(2, 20, Regev, { balanced: true });

    for (const [a, c] of S) {
      // Should be bigint[] when balanced
      expect(Array.isArray(a)).toBe(true);
      expect(typeof c).toBe('bigint');
    }
  });

  /**
   * From balance_sample function:
   * sage: from sage.crypto.lwe import balance_sample, samples, Regev
   * sage: for s in samples(10, 5, Regev):
   * ....:     b = balance_sample(s)
   * ....:     assert all(-29//2 <= c <= 29//2 for c in b[0])
   * ....:     assert -29//2 <= b[1] <= 29//2
   *
   * For n=5, q = next_prime(25) = 29
   */
  test('lwe.py - balance_sample bounds check', () => {
    const regev = new Regev(5);
    const q = regev.K.modulus;
    const halfQ = q / 2n;

    for (let i = 0; i < 10; i++) {
      const sample = regev.call();
      const [bA, bC] = balance_sample(sample, q);

      // All values in a should be balanced
      for (const c of bA) {
        expect(c >= -halfQ).toBe(true);
        expect(c <= halfQ).toBe(true);
      }

      // c should be balanced
      expect(bC >= -halfQ).toBe(true);
      expect(bC <= halfQ).toBe(true);
    }
  });

  /**
   * From balance_sample function:
   * sage: assert all(s[0][j] == b[0][j] % 29 for j in range(5))
   * sage: assert s[1] == b[1] % 29
   */
  test('lwe.py - balance_sample mod equivalence', () => {
    const regev = new Regev(5);
    const q = regev.K.modulus;

    for (let i = 0; i < 10; i++) {
      const sample = regev.call();
      const [a, c] = sample;
      const [bA, bC] = balance_sample(sample, q);

      // Original and balanced should be congruent mod q
      for (let j = 0; j < 5; j++) {
        const original = a.get(j).value;
        const balanced = bA[j]!;
        // balanced mod q should equal original
        const balancedMod = ((balanced % q) + q) % q;
        expect(balancedMod).toBe(original);
      }

      const originalC = c.value;
      const balancedCMod = ((bC % q) + q) % q;
      expect(balancedCMod).toBe(originalC);
    }
  });
});

/**
 * Tests ported from sage/crypto/lattice.py
 *
 * Note: gen_lattice is not yet implemented (throws NotImplementedError),
 * so we test the expected interface and error handling.
 */
describe('lattice.py doctests', () => {
  /**
   * From gen_lattice docstring:
   * sage: sage.crypto.gen_lattice(m=10, seed=42)
   * [11  0  0  0  0  0  0  0  0  0]
   * [ 0 11  0  0  0  0  0  0  0  0]
   * ...
   *
   * Since gen_lattice is not implemented, we test the import
   */
  test('lattice.py - gen_lattice function exists', async () => {
    const { gen_lattice } = await import('./lattice.js');
    expect(typeof gen_lattice).toBe('function');
  });

  /**
   * From gen_lattice docstring - expected interface test
   * These tests document the expected API even though not implemented
   */
  test('lattice.py - gen_lattice throws NotImplementedError', async () => {
    const { gen_lattice } = await import('./lattice.js');
    expect(() => gen_lattice({ m: 10, seed: 42 })).toThrow('SAGE_NOT_IMPLEMENTED');
  });

  test('lattice.py - gen_lattice types interface exists', async () => {
    const lattice = await import('./lattice.js');
    expect(lattice.gen_lattice).toBeDefined();
    // Type exports would be checked at compile time
  });
});

/**
 * Additional doctest-style tests for LWE error recovery
 * These verify the core LWE functionality
 */
describe('lwe.py - LWE correctness tests', () => {
  /**
   * From LWE class - error recovery test:
   * The error e = c - <a, s> should match the noise distribution
   */
  test('lwe.py - error recovery with known secret', () => {
    const D = new UniformSampler(-3, 3);
    const lwe = new LWE(20, 401n, D);
    const secret = lwe.secret;
    const q = 401n;
    const q2 = q / 2n;

    // Generate samples and verify error recovery
    for (let i = 0; i < 20; i++) {
      const [a, c] = lwe.call();

      // c = <a, s> + e, so e = c - <a, s>
      const as = a.dotProduct(secret);
      const recovered_e = c.sub(as);

      // Balance the error
      const e_val = recovered_e.value <= q2 ? recovered_e.value : recovered_e.value - q;

      // Error should be within the distribution bounds
      expect(e_val >= -3n).toBe(true);
      expect(e_val <= 3n).toBe(true);
    }
  });

  /**
   * Test noise secret distribution produces small secrets
   */
  test('lwe.py - noise secret distribution', () => {
    const D = new UniformSampler(-5, 5);
    const lwe = new LWE(10, 101n, D, 'noise');
    const q = 101n;
    const q2 = q / 2n;

    for (let i = 0; i < lwe.n; i++) {
      const val = lwe.secret.get(i).value;
      const balanced = val <= q2 ? val : val - q;
      expect(balanced >= -5n).toBe(true);
      expect(balanced <= 5n).toBe(true);
    }
  });

  /**
   * Test that samples function with instance validates dimension
   */
  test('lwe.py - samples dimension validation', () => {
    const lwe = new Regev(15);
    expect(() => samples(10, 20, lwe)).toThrow(
      'Passed LWE instance has n=15, but n=20 was passed to this function'
    );
  });
});

/**
 * LWEVector doctest-style tests
 */
describe('lwe.py - LWEVector tests', () => {
  /**
   * Test vector dot product
   */
  test('lwe.py - LWEVector dot product', () => {
    const ring = Zmod(101n);
    const v1 = new LWEVector(
      [ring.__call__(3n), ring.__call__(5n), ring.__call__(7n)],
      ring
    );
    const v2 = new LWEVector(
      [ring.__call__(2n), ring.__call__(4n), ring.__call__(6n)],
      ring
    );

    // 3*2 + 5*4 + 7*6 = 6 + 20 + 42 = 68 mod 101 = 68
    const dot = v1.dotProduct(v2);
    expect(dot.value).toBe(68n);
  });

  /**
   * Test vector string representation
   */
  test('lwe.py - LWEVector toString', () => {
    const ring = Zmod(17n);
    const v = new LWEVector(
      [ring.__call__(3n), ring.__call__(5n), ring.__call__(7n)],
      ring
    );
    expect(v.toString()).toBe('(3, 5, 7)');
  });
});

/**
 * Edge case tests from doctests
 */
describe('lwe.py - edge cases', () => {
  /**
   * Test UniformSampler with same bounds
   */
  test('lwe.py - UniformSampler single value', () => {
    const sampler = new UniformSampler(5, 5);
    for (let i = 0; i < 10; i++) {
      expect(sampler.call()).toBe(5n);
    }
  });

  /**
   * Test LWE with m=0 (no samples allowed)
   */
  test('lwe.py - LWE with m=0', () => {
    const D = new UniformSampler(-1, 1);
    const lwe = new LWE(5, 17n, D, 'uniform', 0);
    expect(() => lwe.call()).toThrow('Number of available samples exhausted');
  });

  /**
   * Test Regev with very small n
   */
  test('lwe.py - Regev with n=3', () => {
    const regev = new Regev(3);
    // q = next_prime(9) = 11
    expect(regev.K.modulus).toBe(11n);
  });

  /**
   * Test samples with string 'LindnerPeikert'
   */
  test('lwe.py - samples with LindnerPeikert string', () => {
    const S = samples(5, 20, 'LindnerPeikert');
    expect(S.length).toBe(5);
    const [a, _c] = S[0]!;
    if (a instanceof LWEVector) {
      expect(a.length).toBe(20);
    }
  });

  /**
   * Test unknown oracle name
   */
  test('lwe.py - samples with unknown oracle', () => {
    expect(() => samples(5, 20, 'UnknownOracle' as 'Regev')).toThrow('Unknown LWE oracle');
  });
});

/**
 * Statistical tests (not exact doctests, but verify statistical properties)
 */
describe('lwe.py - statistical properties', () => {
  /**
   * Verify that uniform sampling produces roughly uniform distribution
   */
  test('lwe.py - UniformSampler distribution', () => {
    const sampler = new UniformSampler(-5, 5);
    const counts = new Map<bigint, number>();

    for (let i = 0; i < 1000; i++) {
      const val = sampler.call();
      counts.set(val, (counts.get(val) || 0) + 1);
    }

    // Should have all values from -5 to 5
    expect(counts.size).toBe(11);

    // Each value should appear roughly 1000/11 ~ 91 times
    // Allow for significant variance
    for (let v = -5n; v <= 5n; v++) {
      const count = counts.get(v) || 0;
      expect(count).toBeGreaterThan(30);
      expect(count).toBeLessThan(200);
    }
  });

  /**
   * Verify LWE c values appear random (without knowing secret)
   */
  test('lwe.py - LWE sample randomness', () => {
    const D = new UniformSampler(-2, 2);
    const lwe = new LWE(10, 101n, D);

    const cValues = new Set<bigint>();
    for (let i = 0; i < 100; i++) {
      const [_a, c] = lwe.call();
      cValues.add(c.value);
    }

    // Should have many distinct c values
    expect(cValues.size).toBeGreaterThan(50);
  });
});
