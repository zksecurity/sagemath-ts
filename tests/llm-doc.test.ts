/**
 * Executable check for the examples in `LLM.md`.
 *
 * `LLM.md` is the file agents and downstream tooling read to learn this API, and it is
 * the one document with no other test touching it. Every code sample there is mirrored
 * here so a signature change breaks the build instead of quietly misleading a reader.
 *
 * If a case here fails, fix `LLM.md` in the same commit.
 */

import { describe, expect, test } from 'bun:test';
import {
  CRT_list,
  DiscreteGaussianDistributionIntegerSampler,
  DiscreteGaussianInteger,
  EllipticCurve,
  GF,
  Integer,
  IntegerLattice,
  IntegerMatrixFromEntries,
  LWE,
  Mod,
  Rational,
  Regev,
  ZZ,
  Zmod,
  bsgs,
  createClassicalReedSolomonCode,
  crt,
  divisors,
  embedding_degree,
  euler_phi,
  factor,
  gcd,
  identity_integer_matrix,
  inverse_mod,
  is_prime,
  is_prime_power,
  is_square,
  isqrt,
  lcm,
  lllReduce,
  next_prime,
  power_mod,
  previous_prime,
  prime_factors,
  prime_range,
  sigma,
  squarefree_part,
  toBigInt,
  toRational,
  valuation,
  xgcd,
  zero_integer_matrix,
} from 'sagemath-ts';
import { crypto, arith, coding, groups, modules, schemes, stats } from 'sagemath-ts';

describe('LLM.md — importing', () => {
  test('the bare specifier resolves and the subpath entries exist', async () => {
    const sub = await import('sagemath-ts/arith');
    expect(sub.gcd(12n, 8n)).toBe(4n);
    const rings = await import('sagemath-ts/rings/finite_rings');
    expect(rings.GF(7n).__call__(3n).toString()).toBe('3');
  });
});

describe('LLM.md — the four rules', () => {
  test('JavaScript numbers are rejected for integer arguments', () => {
    // @ts-expect-error: documented as a TypeError, not a silent coercion.
    expect(() => gcd(12, 8)).toThrow(/JavaScript numbers are not accepted/);
  });

  test('strings are not coerced', () => {
    // @ts-expect-error: documented as a TypeError.
    expect(() => toBigInt('123')).toThrow(TypeError);
  });

  test('rings are not callable; __call__ is', () => {
    const F = GF(13n);
    expect(typeof (F as unknown as () => void)).toBe('object');
    expect(F.__call__(5n).pow(12n).toString()).toBe('1');
  });

  test('ring constructors and element arithmetic do accept number', () => {
    expect(GF(7).__call__(10).toString()).toBe('3');
    expect(Zmod(7).__call__(10).toString()).toBe('3');
    expect(Mod(10, 7).toString()).toBe('3');
    expect(GF(7n).__call__(3n).add(4).toString()).toBe('0');
  });

  test('matrices expose determinant(), not det()', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n],
      [3n, 4n],
    ]);
    expect((A as unknown as { det?: unknown }).det).toBeUndefined();
    expect(A.determinant().toString()).toBe('-2');
  });
});

describe('LLM.md — arithmetic', () => {
  test('the arithmetic block', () => {
    expect(gcd(12n, 8n)).toBe(4n);
    expect(lcm(12n, 8n)).toBe(24n);
    expect(xgcd(15n, 6n)).toEqual([3n, 1n, -2n]);
    expect(factor(60n)).toEqual([
      [2n, 2n],
      [3n, 1n],
      [5n, 1n],
    ]);
    expect(is_prime(97n)).toBe(true);
    expect(power_mod(2n, 100n, 1000000007n)).toBe(976371285n);
    expect(inverse_mod(3n, 7n)).toBe(5n);
    expect(euler_phi(60n)).toBe(16n);
    expect(crt(2n, 3n, 5n, 7n)).toBe(17n);
    expect(CRT_list([2n, 3n], [5n, 7n])).toBe(17n);
  });

  test('the primes block', () => {
    expect(next_prime(10n)).toBe(11n);
    expect(previous_prime(10n)).toBe(7n);
    expect(prime_range(10n, 30n)).toEqual([11n, 13n, 17n, 19n, 23n, 29n]);
    expect(prime_factors(60n)).toEqual([2n, 3n, 5n]);
    expect(is_prime_power(8n)).toBe(true);
    expect(is_prime_power(8n, true)).toEqual([2n, 3n]);
  });

  test('the divisibility block', () => {
    expect(divisors(12n)).toEqual([1n, 2n, 3n, 4n, 6n, 12n]);
    expect(sigma(12n, 0n)).toBe(6n);
    expect(sigma(12n, 1n)).toBe(28n);
    expect(valuation(12n, 2n)).toBe(2n);
    expect(is_square(16n)).toBe(true);
    expect(is_square(15n, true)).toEqual([false, null]);
    expect(isqrt(17n)).toBe(4n);
    expect(squarefree_part(12n)).toBe(3n);
  });
});

describe('LLM.md — factorization', () => {
  test('factor is the full PARI cascade, not trial division', () => {
    expect(factor(12345678901234567890n)).toEqual([
      [2n, 1n],
      [3n, 2n],
      [5n, 1n],
      [101n, 1n],
      [3541n, 1n],
      [3607n, 1n],
      [3803n, 1n],
      [27961n, 1n],
    ]);
  });

  test('a 128-bit semiprime stays well inside the documented budget', () => {
    const p = next_prime((1n << 64n) + 12345n);
    const q = next_prime((1n << 64n) + 98765n);
    const started = performance.now();
    expect(factor(p * q)).toEqual([
      [p, 1n],
      [q, 1n],
    ]);
    // Documented as ~56 ms; the assertion is loose enough to survive slow CI but tight
    // enough to catch a regression to trial division, which could never finish.
    expect(performance.now() - started).toBeLessThan(10_000);
  });

  test('factor(0) throws, as in SageMath', () => {
    expect(() => factor(0n)).toThrow();
  });
});

describe('LLM.md — rings and fields', () => {
  test('the rings block', () => {
    expect(ZZ.__call__(42n)).toBe(42n);
    expect(new Integer(42n).factor()).toEqual([
      [2n, 1n],
      [3n, 1n],
      [7n, 1n],
    ]);
    expect(new Rational(3n, 4n).add(new Rational(1n, 2n)).toString()).toBe('5/4');

    const Z7 = Zmod(7n);
    expect(Z7.__call__(3n).mul(Z7.__call__(5n)).toString()).toBe('1');
    expect(Mod(10n, 7n).toString()).toBe('3');

    const F13 = GF(13n);
    expect(F13.__call__(5n).pow(12n).toString()).toBe('1');
  });
});

describe('LLM.md — matrices', () => {
  test('the integer matrix block', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n],
      [3n, 4n],
    ]);
    const B = IntegerMatrixFromEntries([
      [5n, 6n],
      [7n, 8n],
    ]);

    expect(A.mul(B).toString()).toBe(
      IntegerMatrixFromEntries([
        [19n, 22n],
        [43n, 50n],
      ]).toString()
    );
    expect(A.determinant().toString()).toBe('-2');
    expect(A.rank()).toBe(2);
    expect(identity_integer_matrix(3).toString()).toBe(
      IntegerMatrixFromEntries([
        [1n, 0n, 0n],
        [0n, 1n, 0n],
        [0n, 0n, 1n],
      ]).toString()
    );
    expect(zero_integer_matrix(2, 3).toString()).toBe(
      IntegerMatrixFromEntries([
        [0n, 0n, 0n],
        [0n, 0n, 0n],
      ]).toString()
    );
    expect(A.hermite_form()).toBeDefined();
    expect(A.smith_form()).toBeDefined();
    expect(A.elementary_divisors()).toBeDefined();
    expect(A.right_kernel_matrix()).toBeDefined();
  });

  test('IntegerMatrix has no inverse(), as documented', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);
    expect((A as unknown as { inverse?: unknown }).inverse).toBeUndefined();
  });
});

describe('LLM.md — elliptic curves', () => {
  const F = GF(101n);
  const E = EllipticCurve(F, [1n, 2n]);

  test('the curve block', () => {
    expect(E.order()).toBe(100n);
    expect(E.discriminant().toString()).toBe('26');
    expect(E.j_invariant().toString()).toBe('4');

    const P = E.point(1n, 2n);
    expect(P.order()).toBe(4n);
    expect(P.mul(4n).isZero()).toBe(true);

    expect(E.is_on_curve(F.__call__(1n), F.__call__(2n))).toBe(true);
    expect(E.is_on_curve(F.__call__(3n), F.__call__(6n))).toBe(false);
  });

  test('lift_x throws when no point has that x-coordinate', () => {
    expect(() => E.lift_x(F.__call__(3n))).toThrow(/No point with x-coordinate/);
  });

  test('a random point is on the curve and killed by the group order', () => {
    const R = E.random_point();
    expect(R.mul(E.order()).isZero()).toBe(true);
  });

  test('embedding_degree', () => {
    expect(embedding_degree(E, 5n)).toBe(1n);
  });
});

describe('LLM.md — lattices', () => {
  test('IntegerLattice and lllReduce take matrices or bigint[][], not vectors', () => {
    expect(
      IntegerLattice([
        [1n, 0n, 0n],
        [0n, 1n, 0n],
        [0n, 0n, 1n],
      ])
    ).toBeDefined();
    const reduced = lllReduce(
      IntegerMatrixFromEntries([
        [1n, 2n, 3n],
        [4n, 5n, 6n],
        [7n, 8n, 10n],
      ])
    );
    expect(reduced.toString()).toBe(
      IntegerMatrixFromEntries([
        [0n, 0n, 1n],
        [-1n, 1n, 0n],
        [2n, 1n, 0n],
      ]).toString()
    );
  });
});

describe('LLM.md — LWE and discrete Gaussians', () => {
  test('LWE oracles are sampled with call(), positionally constructed', () => {
    const D = DiscreteGaussianInteger(3.2);
    expect(typeof D.call()).toBe('bigint');

    const lwe = new LWE(32n, 40961n, D);
    const [a, c] = lwe.call();
    expect(a.length).toBe(32);
    expect(c).toBeDefined();
    expect(lwe.samples(3n).length).toBe(3);
  });

  test('Regev is an LWE parameter set, not an encryption scheme', () => {
    const regev = new Regev(32n);
    const [a] = regev.call();
    expect(a.length).toBe(32);
    for (const absent of ['keygen', 'encrypt', 'decrypt']) {
      expect((regev as unknown as Record<string, unknown>)[absent]).toBeUndefined();
    }
  });

  test('the sampler class takes an options object', () => {
    const sampler = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3.2 });
    expect(typeof sampler.call()).toBe('bigint');
    expect(DiscreteGaussianInteger(3.2, 5n, 6n)).toBeDefined();
  });
});

describe('LLM.md — coding theory', () => {
  test('Reed-Solomon length must divide q-1', () => {
    const F = GF(929n);
    const RS = createClassicalReedSolomonCode(F, 928n, 4n);
    expect(RS.minimum_distance()).toBe(925);
  });
});

describe('LLM.md — groups', () => {
  test('bsgs is positional and defaults to the multiplicative operation', () => {
    const F = GF(101n);
    const g = F.__call__(2n);
    expect(bsgs(g, g.pow(37n), [0n, 100n])).toBe(37n);
  });
});

describe('LLM.md — coercion and namespaces', () => {
  test('the coercion block', () => {
    expect(toBigInt(42n)).toBe(42n);
    expect(toBigInt(new Integer(42n))).toBe(42n);
    expect(toRational(3n).toString()).toBe('3');
  });

  test('every documented namespace is exported', () => {
    expect(arith.gcd(12n, 8n)).toBe(4n);
    expect(typeof crypto.LWE).toBe('function');
    expect(typeof coding.ReedSolomonCode).toBe('function');
    expect(typeof groups.bsgs).toBe('function');
    expect(typeof modules.lllReduce).toBe('function');
    expect(typeof schemes.EllipticCurve).toBe('function');
    expect(typeof stats.DiscreteGaussianInteger).toBe('function');
  });
});
