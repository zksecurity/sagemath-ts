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
  OverflowError,
  QQ,
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
import { FpXQ_pow, FpXQ_inv, Fp_order, Z_isanypower } from '../packages/parigp-ts/src/index.js';

describe('LLM.md — importing', () => {
  test('the bare specifier resolves and the subpath entries exist', async () => {
    const sub = await import('sagemath-ts/arith');
    expect(sub.gcd(12n, 8n)).toBe(4n);
    const rings = await import('sagemath-ts/rings/finite_rings');
    expect(rings.GF(7n).__call__(3n).toString()).toBe('3');
  });

  test('the package root is a subset; the named gaps live on their subpaths', async () => {
    const root = new Set(Object.keys(await import('sagemath-ts')));

    // Documented as not-at-root. If one of these is later re-exported, update the
    // LLM.md table rather than deleting the case.
    for (const name of ['nth_prime', 'binomial', 'factorial', 'fibonacci', 'kronecker']) {
      expect(root.has(name)).toBe(false);
    }
    const ar = await import('sagemath-ts/arith');
    for (const name of ['nth_prime', 'binomial', 'factorial', 'fibonacci', 'kronecker']) {
      expect(name in ar).toBe(true);
    }

    for (const name of ['PolynomialRing', 'Polynomial', 'NumberField']) {
      expect(root.has(name)).toBe(false);
    }
    const rings = await import('sagemath-ts/rings');
    for (const name of ['PolynomialRing', 'Polynomial', 'NumberField']) {
      expect(name in rings).toBe(true);
    }

    // Every rings/polynomial export is absent from the root — the table says "all 82".
    const poly = await import('sagemath-ts/rings/polynomial');
    expect(Object.keys(poly).filter((n) => root.has(n))).toEqual([]);
  });

  test('exactly two names are bound to different implementations per entry point', async () => {
    const rootMod: Record<string, unknown> = await import('sagemath-ts');
    const subpaths = [
      'sagemath-ts/arith',
      'sagemath-ts/crypto',
      'sagemath-ts/rings',
      'sagemath-ts/rings/finite_rings',
      'sagemath-ts/rings/polynomial',
      'sagemath-ts/matrix',
      'sagemath-ts/schemes/elliptic_curves',
      'sagemath-ts/stats',
    ];
    const divergent = new Set<string>();
    for (const spec of subpaths) {
      const mod: Record<string, unknown> = await import(spec);
      for (const [name, value] of Object.entries(mod)) {
        if (name in rootMod && value !== rootMod[name]) divergent.add(name);
      }
    }
    expect([...divergent].sort()).toEqual(['order_from_multiple', 'rational_reconstruction']);
  });
});

describe('LLM.md — the four rules', () => {
  test('JavaScript numbers are rejected for integer arguments', () => {
    // @ts-expect-error: documented as a TypeError, not a silent coercion.
    expect(() => gcd(12, 8)).toThrow(/JavaScript numbers are not accepted/);
  });

  test('nth_prime accepts IntegerLike and first-prime numeric counts truncate', async () => {
    const { nth_prime, primes_first_n } = await import('sagemath-ts/arith');
    expect(nth_prime(5n)).toBe(11n);
    expect(nth_prime(new Integer(5n))).toBe(11n);
    expect(nth_prime(10000000n)).toBe(179424673n);
    expect(primes_first_n(2.5)).toEqual([2n, 3n]);
    expect(primes_first_n(0.5)).toEqual([]);
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

describe('LLM.md — the GF alias', () => {
  test('GF is GFExtended: prime powers, generator names and modulus options', async () => {
    const ext = await import(
      '../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js'
    );
    const ctor = await import(
      '../packages/sagemath-ts/src/rings/finite_rings/finite_field_constructor.js'
    );

    // The publicly reachable GF — root and subpath alike — is the extension alias.
    const sub = await import('sagemath-ts/rings/finite_rings');
    expect(GF).toBe(ext.GFExtended);
    expect(sub.GF).toBe(GF);
    expect(GF).not.toBe(ctor.GF);

    expect(GF(13n).constructor.name).toBe('PrimeField');
    expect(GF(4n).constructor.name).toBe('FiniteFieldExtension');
    // The prime-only constructor cannot do this, which is why the alias shadows it.
    expect(() => ctor.GF(4n)).toThrow();

    expect(GF(8n, 'b').constructor.name).toBe('FiniteFieldExtension');
    expect(
      GF(9n, { name: 'b', modulus: [1n, 0n, 1n] })
        .gen()
        .pow(2n)
        .toString()
    ).toBe('2');
    expect(
      GF(7n, 'a', { modulus: [-2n, 1n] })
        .gen()
        .toString()
    ).toBe('2');
    expect(GF(9n, { name: 'b', modulus: 'conway' }).gen().toString()).toBe('b');
    expect(() => GF(13n, { check: false } as unknown as string)).toThrow(
      "create_key_and_extra_args() got an unexpected keyword argument 'check'"
    );
  });
});

describe('LLM.md — rings and fields', () => {
  test('the rings block', () => {
    expect(ZZ.__call__(42n)).toBe(42n);
    expect(new Integer(2n).pow(-3n)).toBeInstanceOf(Rational);
    expect(String(new Integer(2n).pow(-3n))).toBe('1/8');
    expect(new Integer(0n).exact_log(2n)).toBe('-Infinity');
    expect(new Integer(0n).log(2n)).toBe('-Infinity');
    expect(new Integer(0n).log()).toBe('-Infinity');
    expect(String(new Integer(5n).__invert__())).toBe('1/5');
    expect(new Integer(0n).valuation(2n)).toBe('Infinity');
    expect(new Integer(0n).ord(2n)).toBe('Infinity');
    expect(new Integer(-1n).popcount()).toBe('Infinity');
    expect(new Integer(-1n).hamming_weight()).toBe('Infinity');
    expect(new Integer(5n).mod(-7n).value).toBe(5n);
    expect(new Integer(5n).__mod__(-7n).value).toBe(-2n);
    expect(new Integer(5n).mod(0n).value).toBe(5n);
    expect(String(new Integer(-1n).catalan_number())).toBe('-1/2');
    expect(new Integer(-179n).class_number()).toBe(5n);
    expect(ZZ.__call__()).toBe(0n);
    expect(ZZ.__call__(null)).toBe(0n);
    expect(ZZ.__call__('010')).toBe(10n);
    expect(ZZ.__call__('1_000')).toBe(1000n);
    expect(new Integer('ff', 16n).value).toBe(255n);
    expect(new Integer([1n, 0n, 1n], 2n).value).toBe(5n);
    expect(new Integer(new Rational(3n)).value).toBe(3n);
    expect(() => new Integer(new Rational(3n, 2n))).toThrow(
      'no conversion of this rational to integer'
    );
    expect(() => ZZ.__call__(new Rational(3n), 0n)).toThrow('_call_with_args not overridden');
    expect(new Integer(42n).factor()).toEqual([
      [2n, 1n],
      [3n, 1n],
      [7n, 1n],
    ]);
    expect(new Rational(3n, 4n).add(new Rational(1n, 2n)).toString()).toBe('5/4');
    expect(new Rational().toString()).toBe('0');
    expect(Rational.from(new Integer(7n)).toString()).toBe('7');
    expect(Rational.from().toString()).toBe('0');
    expect(Rational.from(null).toString()).toBe('0');
    expect(Rational.from(false).toString()).toBe('0');
    expect(QQ.__call__().toString()).toBe('0');
    expect(QQ.__call__(null).toString()).toBe('0');
    expect(() => QQ.__call__([null])).toThrow('unable to convert None');
    expect(QQ.__call__(true).toString()).toBe('1');
    expect(Rational.from(1 / 3).toString()).toBe('1/3');
    expect(Rational.fromString('010').toString()).toBe('8');
    expect(() => Rational.fromString('1.5')).toThrow('unable to convert');
    expect(Rational.from(1.5).toString()).toBe('3/2');
    expect(() => new Rational(1n, 0n)).toThrow('denominator must not be 0');
    expect(new Rational(1n, 2n ** 1074n).toNumber()).toBe(5e-324);
    expect(new Rational(0n).ord(2n)).toBe('Infinity');
    expect(new Rational(0n).gamma().toString()).toBe('Infinity');
    expect(new Rational(250n).roundToRational(-2).toString()).toBe('200');
    expect([...QQ.some_elements()]).toHaveLength(100);
    expect([...QQ.primes_of_bounded_norm_iter(new Rational(5n, 2n))]).toEqual([2n, 3n]);
    expect(() => QQ.gen(1)).toThrow('n must be 0');
    expect(() => QQ.__call__([new Rational(1n, 2n), 3n])).toThrow(
      'no conversion of this rational to integer'
    );
    expect(QQ.random_element()).toBeInstanceOf(Rational);
    const [root, exact] = new Integer(2n).nth_root(2147483647n, true);
    expect([root.value, exact]).toEqual([1n, false]);
    expect(() => new Integer(1n).nth_root(2147483648n)).toThrow(OverflowError);

    expect(Zmod()).toBe(ZZ);
    expect(Zmod(0n)).toBe(ZZ);
    expect(Zmod(-7n).modulus).toBe(7n);
    expect(Mod(10n, 0n)).toBe(10n);
    expect(Mod(10n, 12n, Zmod(7n)).value).toBe(3n);
    expect(Mod(new Rational(1n, 2n), 12n, GF(7n)).is_square()).toBe(true);
    const Z7 = Zmod(7n);
    expect(Z7.__call__(3n).mul(Z7.__call__(5n)).toString()).toBe('1');
    expect(Mod(10n, 7n).toString()).toBe('3');
    expect(Zmod(12n).__call__(new Rational(1n, 5n)).toString()).toBe('5');
    expect(Zmod(12n).__call__().toString()).toBe('0');
    expect(Mod(3n, 12n).add(Mod(2n, 30n)).modulus).toBe(6n);
    expect(Mod(3n, 12n).add(GF(2n).one()).toString()).toBe('0');

    const F13 = GF(13n);
    expect(F13.__call__(5n).pow(12n).toString()).toBe('1');
    expect(F13.__call__(5n).pow(new Integer(-2n)).toString()).toBe('12');
    expect(F13.__call__(new Integer(20n)).toString()).toBe('7');
    expect(F13.__call__(new Rational(1n, 2n)).toString()).toBe('7');
    expect(F13.__call__().toString()).toBe('0');
    expect(F13.__call__(1n).eq(new Integer(14n))).toBe(true);
    expect(Fp_order(2n, 12n, 13n)).toBe(12n);
    expect(GF(9n).__call__('a^2 + 1').toString()).toBe('a + 2');
    expect(GF(9n).__call__('1/a*a').toString()).toBe('1');
    expect(() => GF(9n).__call__('a^-1')).toThrow('not integral');
    expect(GF(7n).__call__(3n).sqrt().toString()).toBe('sqrt3');
    expect(GF(7n).__call__(3n).sqrt().pow(2n).toString()).toBe('3');
    expect(GF(7n).__call__(2n).sqrt({ all: true }).map(String)).toEqual(['3', '4']);
    expect(GF(7n).__call__(3n).sqrt({ extend: false, all: true })).toEqual([]);
    expect(Z_isanypower(676n)).toEqual([2, 26n]);
    expect(Z_isanypower(37n)).toEqual([0, 37n]);
    expect(Z_isanypower(-64n)).toEqual([3, -4n]);
    expect(Z_isanypower(-16n)).toEqual([0, -16n]);
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

  test('the generic matrix factory works over a field but not over ZZ', async () => {
    const { matrix } = await import('sagemath-ts/matrix');
    const F = GF(7n);
    const M = matrix(F, [
      [1n, 2n],
      [3n, 4n],
    ]);
    const I = matrix(F, [
      [1n, 0n],
      [0n, 1n],
    ]);
    expect(M.mul(I).toString()).toBe(M.toString());

    // ZZ's elements are bare bigints with no .mul, so the same call cannot multiply.
    expect(() =>
      matrix(ZZ, [
        [1n, 2n],
        [3n, 4n],
      ]).mul(
        matrix(ZZ, [
          [1n, 0n],
          [0n, 1n],
        ])
      )
    ).toThrow();
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

test('LLM.md — mixed finite-ring promotion and sequence multiplication', () => {
  expect(Mod(2n, 7n).mul('ab')).toBe('abab');
  expect(Mod(2n, 7n).mul([3n, 5n])).toEqual([3n, 5n, 3n, 5n]);
  expect(Mod(2n, 12n).add(GF(9n).gen()).toString()).toBe('a + 2');
  expect(GF(3n).one().eq(GF(9n).one())).toBe(true);
});

test('LLM.md — finite generator indices and integer conversion hooks', () => {
  expect(QQ.gen(0n).toString()).toBe('1');
  expect(ZZ.__call__(GF(7n).__call__(5n))).toBe(5n);
  expect(new Integer(GF(9n).__call__(2n)).value).toBe(2n);
  expect(GF(7n).gen(false).toString()).toBe('1');
  expect(Zmod(7n).gen(new Rational(0n)).toString()).toBe('1');
  expect(() => ZZ.__call__(GF(9n).gen())).toThrow('element is not in the prime field');
});

test('LLM.md — PARI quotient kernels and extension exponent validation', () => {
  expect(GF(9n).gen().pow(new Rational(-1n)).toString()).toBe('a + 2');
  expect(() => GF(9n).zero().pow(new Rational(-1n, 2n))).toThrow('');
  expect(() => GF(9n).gen().pow(new Rational(-1n, 2n))).toThrow(
    'no conversion of this rational to integer'
  );
  expect(FpXQ_pow([0n, 1n], -1n, [1n, 0n, 1n], 3n)).toEqual([0n, 2n]);
  expect(FpXQ_inv([1n, 1n], [1n, 0n, 1n], 9n, 3n)).toEqual([5n, 4n]);
});

test('LLM.md polynomial coefficient coercion and rational lifts', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const P = new PolynomialRing(GF(7n), 'x');
  expect(P.__call__([9n, 3n, 7n]).toString()).toBe('3*x + 2');
  expect(P.__call__().toString()).toBe('0');
  expect(P.gen(0n).toString()).toBe('x');
  expect(QQ.__call__(GF(7n).__call__(-2n)).toString()).toBe('5');
});

test('LLM.md — dense FLINT polynomial GCD', async () => {
  const { _fmpz_poly_gcd } = await import('../packages/flint-ts/src/index.js');
  expect(_fmpz_poly_gcd([-6n, 0n, 6n], [9n, -9n])).toEqual([-3n, 3n]);
});

test('LLM.md — FLINT monic normalization', async () => {
  const { _nmod_poly_make_monic } = await import('../packages/flint-ts/src/index.js');
  expect(_nmod_poly_make_monic([1n, 2n, 3n], 5n)).toEqual([2n, 4n, 1n]);
});

test('polynomial index conversion examples in LLM.md', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const indexExample = new PolynomialRing(QQ, 'x').__call__([1n, 2n, 3n]);
  expect(indexExample.getCoeff(new Integer(1n)).toString()).toBe('2');
  expect(indexExample.truncate(new Rational(3n, 2n)).coeffs.map(String)).toEqual(['1']);
  expect(indexExample.reverse(new Rational(3n, 2n)).coeffs.map(String)).toEqual(['2', '1']);
  expect(indexExample.shift(new Rational(-3n, 2n))?.coeffs.map(String)).toEqual(['2', '3']);
});

test('polynomial mixed-parent arithmetic and formatting examples in LLM.md', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const primePolynomial = new PolynomialRing(GF(7n), 'x').__call__([1n, 2n]);
  const E = GF(49n, 'a');
  const extensionPolynomial = new PolynomialRing(E, 'x').__call__([3n, 1n]);
  const sum = primePolynomial.add(extensionPolynomial);
  expect(sum.coeffs.map(String)).toEqual(['4', '3']);
  expect(sum.parent.base_ring === E).toBe(true);
  expect(new PolynomialRing(QQ, 'x').__call__([1n, -2n, -1n]).toString()).toBe('-x^2 - 2*x + 1');
});

test('FLINT quotient/remainder examples in LLM.md', async () => {
  const { _fmpz_poly_divrem, _nmod_poly_divrem } = await import(
    '../packages/flint-ts/src/index.js'
  );
  expect(_fmpz_poly_divrem([-7n, -5n], [3n])).toEqual([
    [-3n, -2n],
    [2n, 1n],
  ]);
  expect(_fmpz_poly_divrem([1n, 3n], [0n, 2n], true)).toBeNull();
  expect(_nmod_poly_divrem([1n, 0n, 1n], [1n, 1n], 7n)).toEqual([[6n, 1n], [2n]]);
});

import { _fmpq_poly_gcd } from '../packages/flint-ts/src/fmpq_poly/gcd.js';
test('LLM rational polynomial GCD numerator kernel', () => {
  expect(_fmpq_poly_gcd([-6n, -9n], [4n, 6n])).toEqual([[2n, 3n], 3n]);
});

import { _nmod_poly_xgcd } from '../packages/flint-ts/src/nmod_poly/xgcd.js';
test('LLM modular polynomial XGCD kernel', () => {
  expect(_nmod_poly_xgcd([2n, 2n], [], 7n)).toEqual([[1n, 1n], [4n], []]);
});

import { _fmpz_poly_resultant } from '../packages/flint-ts/src/fmpz_poly/resultant.js';
import { _nmod_poly_resultant } from '../packages/flint-ts/src/nmod_poly/resultant.js';
test('LLM native polynomial resultant kernels', () => {
  expect(_fmpz_poly_resultant([2n, 4n], [-6n, -9n])).toBe(-6n);
  expect(_nmod_poly_resultant([1n], [1n], 7n)).toBe(1n);
});

import { _fmpz_poly_xgcd, _fmpq_poly_xgcd } from '../packages/flint-ts/src/index.js';
test('LLM integer and rational polynomial XGCD kernels', () => {
  expect(_fmpz_poly_xgcd([2n, 1n], [3n, 0n, 0n, 1n])).toEqual([5n, [4n, -2n, 1n], [-1n]]);
  expect(_fmpq_poly_xgcd([2n, 2n], 6n, [], 1n)).toEqual([
    [[1n, 1n], 1n],
    [[3n], 1n],
    [[], 1n],
  ]);
});

test('LLM Sylvester matrices and explicit finite polynomial representatives', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const sylvesterRing = new PolynomialRing(QQ, 'x');
  expect(
    sylvesterRing
      .__call__([1n, 0n, 1n])
      .sylvester_matrix(new Rational(2n))
      .map((row) => row.map(String))
  ).toEqual([
    ['2', '0'],
    ['0', '2'],
  ]);
  expect(sylvesterRing.__call__(GF(49n, 'a').gen()).coeffs.map(String)).toEqual(['0', '1']);
});

test('LLM rational and real resultant APIs', async () => {
  const { RDF } = await import('sagemath-ts');
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const { _fmpq_poly_resultant } = await import('../packages/flint-ts/src/index.js');
  const { PariError, resultant } = await import('../packages/parigp-ts/src/index.js');
  const resultantRing = new PolynomialRing(QQ, 'x');
  expect(resultantRing.__call__([1n, 0n, 1n]).resultant(0.5).toString()).toBe('0.25');
  expect(RDF.__call__(1n << 1200n).toString()).toBe('+infinity');
  expect(_fmpq_poly_resultant([1n, 1n], 2n, [1n, -1n], 3n)).toEqual([1n, 3n]);
  expect(resultant([1, 2, 3, 4], [0.5])).toBe(0.125);
  expect(() => resultant([Infinity, 1], [1, 1])).toThrow(PariError);
});

test('LLM general real resultants and native inexact leading zero', async () => {
  const { PariError, resultant } = await import('../packages/parigp-ts/src/index.js');
  expect(resultant([1, 0, 1], [-2, 0, 0, 1])).toBe(5);
  expect(resultant([1, 0], [1, 1])).toBe(-1);
});

test('LLM polynomial derivative argument protocol', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const derivativeRing = new PolynomialRing(QQ, 'x');
  const derivativeExample = derivativeRing.__call__([-1n, 0n, new Rational(1n, 2n), 0n, -1n]);
  expect(derivativeExample.derivative(derivativeRing.gen(), 2n).coeffs.map(String)).toEqual([
    '1',
    '0',
    '-12',
  ]);
  expect(derivativeExample.diff([])).toBe(derivativeExample);
});

test('LLM derivative utilities and rational kernel', async () => {
  const { derivative_parse } = await import('sagemath-ts/misc/derivative');
  const { _fmpq_poly_derivative } = await import('../packages/flint-ts/src/index.js');
  expect(derivative_parse(['x', 2n, 2n])).toEqual(['x', 'x', null, null]);
  expect(_fmpq_poly_derivative([1n, 0n, 3n], 6n)).toEqual([[0n, 1n], 1n]);
});

test('LLM pseudo-division and constant-denominator fraction imports', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const { FractionField_generic } = await import('sagemath-ts/rings/fraction_field');
  const { FractionFieldElement } = await import('sagemath-ts/rings/fraction_field_element');
  const pseudoRing = new PolynomialRing(QQ, 'x');
  expect(
    pseudoRing
      .__call__(1n)
      .pseudo_quo_rem(pseudoRing.__call__([2n, 0n, 0n, 2n]))
      .map(String)
  ).toEqual(['0', '1/4']);
  const constantFractions = new FractionField_generic(pseudoRing);
  const value = constantFractions.__call__([1n, 2n], 3n).mul(constantFractions.__call__(2n, 3n));
  expect(value).toBeInstanceOf(FractionFieldElement);
  expect(value.toString()).toBe('4/9*x + 2/9');
});

test('LLM nonconstant fractions and native product kernels', async () => {
  const { FractionField } = await import('sagemath-ts/rings/fraction_field');
  const { FpTElement } = await import('sagemath-ts/rings/fraction_field_FpT');
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const { _nmod_poly_mul, _nmod_poly_pow } = await import('../packages/flint-ts/src/index.js');
  const fractions = FractionField(new PolynomialRing(QQ, 'x'));
  expect(fractions.__call__([1n, 2n, 1n], [1n, 1n]).toString()).toBe('x + 1');
  expect(fractions.__call__([1n, 1n]).inv().toString()).toBe('1/(x + 1)');
  expect(typeof FpTElement).toBe('function');
  expect(_nmod_poly_mul([1n, 1n], [6n, 1n], 7n)).toEqual([6n, 0n, 1n]);
  expect(_nmod_poly_pow([1n, 1n], 2n, 7n)).toEqual([1n, 2n, 1n]);
});

test('LLM polynomial coercion and fraction sections', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const sectionRing = new PolynomialRing(QQ, 'x');
  const sectionField = sectionRing.fraction_field();
  expect(sectionRing.has_coerce_map_from(new PolynomialRing(QQ, 'x'))).toBe(true);
  expect(sectionRing.has_coerce_map_from(new PolynomialRing(QQ, 'y'))).toBe(false);
  const sectionValue = sectionField.__call__([1n, 2n, 1n], [1n, 1n]);
  expect(sectionRing.__call__(sectionValue).toString()).toBe('x + 1');
  expect(sectionRing.__call__(sectionValue)).toBe(sectionValue.numerator());
});

test('LLM scalar hooks and partial polynomial sections', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const { FractionField_generic } = await import('sagemath-ts/rings/fraction_field');
  const { FractionFieldElement } = await import('sagemath-ts/rings/fraction_field_element');
  const { GF2 } = await import('../packages/sagemath-ts/src/rings/finite_rings/gf2.js');
  const scalarRing = new PolynomialRing(QQ, 'x');
  expect(QQ.__call__(scalarRing.__call__(new Rational(3n, 2n))).toString()).toBe('3/2');
  expect(scalarRing.__call__(3n)._integer_(ZZ)).toBe(3n);
  expect(scalarRing.__call__(3n)._rational_().toString()).toBe('3');
  expect(scalarRing.__call__(3n)._scalar_conversion(QQ).toString()).toBe('3');
  const scalarFraction = new FractionFieldElement(
    new FractionField_generic(scalarRing),
    [1n, 1n],
    [2n, 2n]
  );
  expect(scalarFraction._rational_().toString()).toBe('1/2');
  expect(scalarFraction._conversion(QQ).toString()).toBe('1/2');
  expect(GF2.__call__(scalarRing.__call__(3n))._integer_()).toBe(1n);
  expect(GF2.__call__(1n)._rational_().toString()).toBe('1');
});

test('polynomial integer powers and native backend exports', async () => {
  const { PolynomialRing, Polynomial } = await import('sagemath-ts/rings');
  const R = new PolynomialRing(QQ, 'x');
  const f = R.gen().add(R.one());
  const positive: typeof f = f.pow(3n);
  expect(positive.toString()).toBe('x^3 + 3*x^2 + 3*x + 1');
  const reciprocal = f.pow(-2n);
  expect(reciprocal.toString()).toBe('1/(x^2 + 2*x + 1)');
  const exponent: bigint = -2n;
  const dynamic = f.pow(exponent);
  if (!(dynamic instanceof Polynomial))
    expect(dynamic.denominator().toString()).toBe('x^2 + 2*x + 1');
  else throw new Error('negative power should be a fraction');
  // A literal union containing a negative value must retain the fraction alternative.
  const mixed = f.pow(-1n as -1n | 1n);
  if (!(mixed instanceof Polynomial)) expect(mixed.numerator().toString()).toBe('1');
  const { _fmpz_poly_pow, _fmpq_poly_pow } = await import('../packages/flint-ts/src/index.js');
  const { GF2X, GF2X_power, ZZ_pX_power, ZZ_pEX_power } = await import(
    '../packages/ntl-ts/src/index.js'
  );
  expect(_fmpz_poly_pow([1n, 1n], 3n)).toEqual([1n, 3n, 3n, 1n]);
  expect(_fmpq_poly_pow([2n], 2n, 3n)).toEqual([[8n], 8n]);
  expect(GF2X_power(new GF2X([1, 1]), 3n).rep()).toBe(15n);
  expect(ZZ_pX_power([1n, 1n], 3n, 7n)).toEqual([1n, 3n, 3n, 1n]);
  expect(ZZ_pEX_power([[0n, 1n]], 2n, [1n, 1n, 1n], 7n)).toEqual([[6n, 6n]]);
});

test('polynomial rational roots and truncated-series API examples', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const R = new PolynomialRing(QQ, 'x'),
    f = R.__call__([1n, 1n]);
  expect(R.__call__([1n, 2n, 1n]).pow(new Rational(1n, 2n)).toString()).toBe('x + 1');
  expect(R.__call__([1n, 2n, 1n]).nth_root(2n).toString()).toBe('x + 1');
  expect(f._nth_root_series(2n, 3n).toString()).toBe('-1/8*x^2 + 1/2*x + 1');
  expect(f.inverse_series_trunc(4n).toString()).toBe('-x^3 + x^2 - x + 1');
  expect(f.power_trunc(5n, 3n).toString()).toBe('10*x^2 + 5*x + 1');
  expect(f._power_trunc(5n, 3n).toString()).toBe('10*x^2 + 5*x + 1');
  expect(f._mul_trunc_(f, 2n).toString()).toBe('2*x + 1');
  expect(f.multiplication_trunc(f, 2n).toString()).toBe('2*x + 1');
  const { generic_power_trunc } = await import(
    '../packages/sagemath-ts/src/rings/polynomial/polynomial_element.js'
  );
  expect(generic_power_trunc(f, 5n, 3).toString()).toBe('10*x^2 + 5*x + 1');
});

test('FLINT and NTL truncated-series export examples', async () => {
  const F = await import('../packages/flint-ts/src/index.js');
  expect(F._fmpz_poly_mullow([1n, 1n], [1n, 1n], 2)).toEqual([1n, 2n]);
  expect(F._fmpq_poly_mullow([1n, 1n], 2n, [1n, 1n], 3n, 2)).toEqual([[1n, 2n], 6n]);
  expect(F._nmod_poly_mullow([1n, 1n], [6n, 1n], 2, 7n)).toEqual([6n]);
  expect(F._fmpz_poly_pow_trunc([1n, 1n], 5n, 3)).toEqual([1n, 5n, 10n]);
  expect(F._nmod_poly_pow_trunc([1n, 1n], 5n, 3, 7n)).toEqual([1n, 5n, 3n]);
  expect(F._fmpz_poly_inv_series([1n, 1n], 4)).toEqual([1n, -1n, 1n, -1n]);
  expect(F._fmpq_poly_inv_series_newton([2n, 1n], 1n, 3)).toEqual([[4n, -2n, 1n], 8n]);
  const { ZZ_pEX_InvTrunc } = await import('../packages/ntl-ts/src/index.js');
  expect(ZZ_pEX_InvTrunc([[1n], [1n]], 4, [1n, 1n, 1n], 7n)).toEqual([[1n], [6n], [1n], [6n]]);
});

test('truncated multiplication selects canonical scalar parents and preserves signed zero', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const { RDF } = await import('../packages/sagemath-ts/src/rings/real_double.js');
  const R = new PolynomialRing(QQ, 'x'),
    f = R.__call__([1n, 2n, 1n]);
  const h = f.multiplication_trunc(0.5, 3n);
  expect(h.toString()).toBe('0.5*x^2 + x + 0.5');
  expect(h.parent.base_ring).toBe(RDF);
  expect(f._mul_trunc_(null, 3n).toString()).toBe('0');
  expect(() => f.multiplication_trunc('2', 3n)).toThrow('no common canonical parent');
  const S = new PolynomialRing(RDF, 'x'),
    g = S.__call__([RDF.__call__(-0), RDF.__call__(1)]);
  expect(Object.is(g._mul_trunc_(S.one(), 3n).getCoeff(0).value, -0)).toBe(true);
});

test('LLM full exact polynomial product adapters', async () => {
  const F = await import('../packages/flint-ts/src/index.js');
  const N = await import('../packages/ntl-ts/src/index.js');
  expect(F._fmpz_poly_mul([1n, 2n], [3n, 4n])).toEqual([3n, 10n, 8n]);
  expect(F._fmpq_poly_mul([2n, 2n], 2n, [1n, 1n], 1n)).toEqual([[2n, 4n, 2n], 2n]);
  const a = [2n, 2n];
  expect(F._fmpq_poly_mul(a, 2n, a, 2n)).toEqual([[4n, 8n, 4n], 4n]);
  expect(N.ZZ_pX_mul([1n, 2n], [3n, 4n], 7n)).toEqual([3n, 3n, 1n]);
  expect(N.ZZ_pEX_mul([[1n], [0n, 1n]], [[1n], [0n, 1n]], [1n, 0n, 1n], 7n)).toEqual([
    [1n],
    [0n, 2n],
    [6n],
  ]);
});

test('LLM modular polynomial and RDF integer powers', async () => {
  const { RDF, GF, Integer } = await import('sagemath-ts');
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const R = new PolynomialRing(GF(7n), 'x'),
    x = R.gen();
  expect(x.add(R.one()).pow(5n, x.pow(2n).add(R.one())).toString()).toBe('3*x + 3');
  expect(RDF.__call__(2).pow(-3n).value).toBe(0.125);
  expect(RDF.__call__(2).pow(new Integer(-3n)).value).toBe(0.125);
  expect(RDF.__call__(4).inv().value).toBe(0.25);
  const F = await import('../packages/flint-ts/src/index.js');
  const m = [1n, 0n, 1n],
    inverse = F._nmod_poly_inv_series_newton(m.slice().reverse(), 3, 7n);
  expect(inverse).toEqual([1n, 0n, 6n]);
  expect(F._nmod_poly_powmod_ui_binexp([1n, 1n], 5n, m, 7n)).toEqual([3n, 3n]);
  expect(F._nmod_poly_powmod_fmpz_binexp_preinv([1n, 1n], 5n, m, inverse, 7n)).toEqual([3n, 3n]);
  expect(F._nmod_poly_powmod_x_fmpz_preinv(5n, m, inverse, 7n)).toEqual([0n, 1n]);
  const { _nmod_poly_preinv_remainder } = await import(
    '../packages/flint-ts/src/nmod_poly/powmod_binexp_preinv.js'
  );
  expect(_nmod_poly_preinv_remainder([1n, 2n, 1n], m, inverse, 7n)).toEqual([0n, 2n]);
  const N = await import('../packages/ntl-ts/src/index.js');
  expect(N.ZZ_pEX_PowerXMod(5n, [[1n], [], [1n]], [1n, 0n, 1n], 7n)).toEqual([[], [1n]]);
  expect(N.ZZ_pEX_PowerMod([[1n], [1n]], 5n, [[1n], [], [1n]], [1n, 0n, 1n], 7n)).toEqual([
    [3n],
    [3n],
  ]);
  expect(N.ZZ_pEX_XGCD([[1n], [1n]], [[1n], [], [1n]], [1n, 0n, 1n], 7n)).toEqual([
    [[1n]],
    [[4n], [3n]],
    [[4n]],
  ]);
  const G = await import('../packages/gsl-ts/src/index.js');
  expect(G.gsl_pow_int(2, -3)).toBe(0.125);
  expect(G.gsl_sf_log(1)).toBe(0);
  expect(G.gsl_sf_exp(0)).toBe(1);
});

describe('LLM.md — polynomial evaluation and composition', () => {
  test('public dispatch, cached generators and compiled evaluation', async () => {
    const { GF, RDF } = await import('sagemath-ts');
    const { PolynomialRing, CompiledPolynomialFunction } = await import(
      'sagemath-ts/rings/polynomial'
    );
    const R = new PolynomialRing(GF(7n), 'x'),
      x = R.gen(),
      f = x.pow(2n).add(R.one());
    expect(f.evaluate(3n).toString()).toBe('3');
    expect(f.compose(x.add(R.one())).toString()).toBe('x^2 + 2*x + 2');
    expect(x === R.gen()).toBe(true);
    expect(x.is_gen()).toBe(true);
    const c = new CompiledPolynomialFunction([1, 0, 0, 0, 1].map((v) => RDF.__call__(v)));
    expect(c.eval(RDF.__call__(2)).value).toBe(17);
    expect(c.__call__(RDF.__call__(2)).value).toBe(17);
    expect(c.toString()).toBe('CompiledPolynomialFunction((a4*((x)^2)^2+a0))');
  });
  test('native array adapters and unreduced rational evaluation', async () => {
    const F = await import('../packages/flint-ts/src/index.js');
    const N = await import('../packages/ntl-ts/src/index.js');
    expect(F._fmpz_poly_evaluate_fmpz([1n, 2n, 3n], 2n)).toBe(17n);
    expect(F._fmpz_poly_evaluate_fmpq([1n, 2n, 3n], 1n, 2n)).toEqual([11n, 4n]);
    expect(F._fmpq_poly_evaluate_fmpz([1n, 2n, 3n], 3n, 2n)).toEqual([17n, 3n]);
    expect(F._fmpq_poly_evaluate_fmpq([1n, 2n, 3n], 3n, 1n, 2n)).toEqual([11n, 12n]);
    expect(F._nmod_poly_evaluate_nmod([1n, 2n, 3n], 2n, 7n)).toBe(3n);
    expect(F._fmpz_poly_compose([1n, 0n, 1n], [1n, 1n])).toEqual([2n, 2n, 1n]);
    expect(F._fmpq_poly_compose([1n, 0n, 1n], 2n, [1n, 1n], 3n)).toEqual([[10n, 2n, 1n], 18n]);
    expect(F._nmod_poly_compose([1n, 0n, 1n], [1n, 1n], 7n)).toEqual([2n, 2n, 1n]);
    expect(F._fmpz_poly_taylor_shift([1n, 0n, 1n], 2n)).toEqual([5n, 4n, 1n]);
    expect(N.ZZ_pX_evaluate([1n, 2n, 3n], 2n, 7n)).toBe(3n);
    expect(N.ZZ_pEX_eval([[1n], [0n, 1n]], [2n], [1n, 0n, 1n], 7n)).toEqual([1n, 2n]);
  });
});

test('LLM.md — polynomial matrix evaluation and native real conversion', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings/polynomial');
  const { integer_to_real_double_dense } = await import('sagemath-ts/matrix');
  const { fmpz_get_d } = await import('../packages/flint-ts/src/index.js');
  const matrixPolynomial = new PolynomialRing(QQ, 'x').__call__([1n, 0n, 1n]);
  const matrixPoint = IntegerMatrixFromEntries([
    [1n, 2n],
    [3n, 4n],
  ]);
  const matrixValue = matrixPolynomial.evaluate(matrixPoint);
  expect(matrixValue.get(0, 0).toString()).toBe('8');
  expect(matrixValue.get(1, 1).toString()).toBe('23');
  const realMatrix = integer_to_real_double_dense(IntegerMatrixFromEntries([[9007199254740995n]]));
  expect(realMatrix.get(0, 0).value).toBe(9007199254740994);
  expect(fmpz_get_d(9007199254740995n)).toBe(9007199254740994);
});

test('LLM.md — specialized binary matrices and packed M4RI adapters', async () => {
  const { Matrix_mod2_dense } = await import('sagemath-ts/matrix');
  const { GF } = await import('sagemath-ts');
  const { PolynomialRing } = await import('sagemath-ts/rings/polynomial');
  const {
    mzd_init,
    mzd_add,
    mzd_mul,
    mzd_mul_m4rm,
    mzd_mul_naive,
    mzd_make_table,
    mzd_init_window,
  } = await import('../packages/m4ri-ts/src/index.js');
  const binaryPoint = new Matrix_mod2_dense(2, 2, [
      [1, 1],
      [0, 1],
    ]),
    binarySquare = binaryPoint.mul(binaryPoint);
  expect(binarySquare.row(0)).toEqual([1, 0]);
  expect(binarySquare.row(1)).toEqual([0, 1]);
  expect(binaryPoint._multiply_m4rm(binaryPoint, 0n).row(0)).toEqual([1, 0]);
  expect(binaryPoint._multiply_strassen(binaryPoint, 64n).row(1)).toEqual([0, 1]);
  const binaryPolynomial = new PolynomialRing(GF(2n), 'x').__call__([1n, 1n]),
    binaryValue = binaryPolynomial.evaluate([binaryPoint]);
  expect((binaryValue as Matrix_mod2_dense).row(0)).toEqual([0, 1]);
  const a = mzd_init(2, 2, [3n, 2n]);
  for (const multiply of [mzd_mul, mzd_mul_m4rm, mzd_mul_naive])
    expect(multiply(a, a).rows).toEqual([1n, 2n]);
  expect(mzd_add(a, a).rows).toEqual([0n, 0n]);
  expect(mzd_make_table(a, 0, 2)).toEqual({ rows: [0n, 3n, 1n, 2n], lookup: [0, 1, 3, 2] });
  expect(mzd_init_window(mzd_init(1, 65, [(1n << 64n) | 1n]), 0, 64, 1, 65).rows).toEqual([1n]);
});

test('LLM.md — binary matrix indices and slices', async () => {
  const { Matrix_mod2_dense } = await import('sagemath-ts/matrix');
  const { mzd_init, mzd_submatrix } = await import('../packages/m4ri-ts/src/index.js');
  const indexed = new Matrix_mod2_dense(2, 3, [
    [0, 1, 0],
    [1, 0, 1],
  ]);
  expect(indexed.get(-1n, new Integer(-1n))).toBe(1);
  expect(indexed.row(new Rational(1n))).toEqual([1, 0, 1]);
  expect(indexed.submatrix(0n, 1n, -2n, -1n).row(1n)).toEqual([0, 1]);
  expect(indexed.get(new Rational(1n, 2n), 1n)).toBe(1);
  expect(() => indexed.row(new Rational(1n, 2n))).toThrow(
    'unable to convert rational 1/2 to an integer'
  );
  indexed.set(-1n, -1n, 0);
  expect(indexed.get(1n, 2n)).toBe(0);
  expect(mzd_submatrix(mzd_init(1, 5, [22n]), 0, 1, 1, 4).rows).toEqual([3n]);
});

test('LLM.md — exact binary density and native real literals', async () => {
  const { Matrix_mod2_dense } = await import('sagemath-ts/matrix');
  const { create_RealNumber } = await import('sagemath-ts/rings');
  const densityMatrix = new Matrix_mod2_dense(1, 3, [[1, 0, 0]]);
  expect(densityMatrix.density().toString()).toBe('1/3');
  expect(densityMatrix.density(true).toString()).toBe('0.333333333333333');
  const exactLiteral = create_RealNumber('9007199254740993');
  expect(exactLiteral.exact_rational()).toEqual([9007199254740993n, 1n]);
  expect(exactLiteral.toNumber()).toBe(9007199254740992);
  expect(exactLiteral.numerical_approx(128).exact_rational()).toEqual([9007199254740993n, 1n]);
});

test('LLM.md — binary elimination cache and immutability', async () => {
  const { Matrix_mod2_dense } = await import('sagemath-ts/matrix');
  const eliminationMatrix = new Matrix_mod2_dense(2, 2, [
    [1, 1],
    [0, 1],
  ]);
  const cachedForm = eliminationMatrix.echelon_form('m4ri', 0, false, { k: 1n });
  expect(cachedForm.row(0)).toEqual([1, 1]);
  expect(cachedForm.is_immutable()).toBe(true);
  expect(eliminationMatrix.echelon_form('bogus')).toBe(cachedForm);
  expect(eliminationMatrix.rank('m4ri')).toBe(2);
  expect(eliminationMatrix.rank('bogus')).toBe(2);
  eliminationMatrix.set(0, 0, 0);
  expect(eliminationMatrix.rank()).toBe(1);
  expect(cachedForm.copy().is_mutable()).toBe(true);
});

test('LLM.md — standalone binary factorization parameters', async () => {
  const { Matrix_mod2_dense, ple, pluq } = await import('sagemath-ts/matrix');
  const factorMatrix = new Matrix_mod2_dense(4, 4, [
    [0, 1, 0, 1],
    [0, 1, 1, 1],
    [0, 0, 0, 1],
    [0, 1, 1, 0],
  ]);
  expect(ple(factorMatrix, 'russian', new Rational(3n, 2n))[0].row(0)).toEqual([1, 0, 0, 1]);
  expect(pluq(factorMatrix, 'naive', 0n)[0].row(0)).toEqual([1, 0, 1, 0]);
});

test('LLM.md — binary swaps, permutations and inverse', async () => {
  const { Matrix_mod2_dense } = await import('sagemath-ts/matrix');
  const swapExample = new Matrix_mod2_dense(3, 3, [
    [1, 1, 0],
    [0, 1, 1],
    [0, 0, 1],
  ]);
  swapExample.swap_rows(0n, 1n);
  swapExample.permute_rows([2, 1]);
  swapExample.permute_columns([1, 2, 3, 4]);
  expect(swapExample.inverse().row(0)).toEqual([1, 1, 1]);
});

test('LLM.md — binary linear systems and kernels', async () => {
  const { Matrix_mod2_dense } = await import('sagemath-ts/matrix');
  const system = new Matrix_mod2_dense(2, 3, [
    [1, 0, 1],
    [0, 1, 1],
  ]);
  const rhs = new Matrix_mod2_dense(2, 1, [[1], [0]]);
  expect(system.solve_right(rhs).list()).toEqual([1, 0, 0]);
  expect(system.right_kernel_matrix({ basis: 'pivot', algorithm: 'pluq' }).row(0)).toEqual([
    1, 1, 1,
  ]);
});

test('LLM.md — binary column ownership', async () => {
  const { Matrix_mod2_dense } = await import('sagemath-ts/matrix');
  const columnExample = new Matrix_mod2_dense(2, 2, [
    [1, 0],
    [1, 1],
  ]);
  const cachedColumns = columnExample.columns(false);
  expect(columnExample.columns()[0]).toBe(cachedColumns[0]);
  const editableColumn = cachedColumns[0]!.slice();
  editableColumn[0] = 0;
  expect(columnExample.columns()[0]).toEqual([1, 1]);
  expect(columnExample.transpose().row(0)).toEqual([1, 1]);
});

test('LLM.md — binary formatting and subdivisions', async () => {
  const { Matrix_mod2_dense } = await import('sagemath-ts/matrix');
  const formatted = new Matrix_mod2_dense(2, 3, [
    [0, 1, 0],
    [1, 0, 1],
  ]);
  const mapping = { 0: 'zero', 1: 'x' };
  expect(formatted.str(mapping, '.', 'one')).toBe('[  . one   .]\n[one   . one]');
  expect(mapping).toEqual({ 0: '.', 1: 'one' });
  formatted.subdivide([1n], [1n]);
  expect(formatted.augment(formatted, true).subdivisions()).toEqual([[1n], [1n, 3n, 4n]]);
  expect(formatted.transpose().get_subdivisions()).toEqual([[1n], [1n]]);
});

test('LLM.md — binary construction and cached rows', async () => {
  const { Matrix_mod2_dense } = await import('sagemath-ts/matrix');
  const { Rational } = await import('sagemath-ts');
  const constructed = new Matrix_mod2_dense(2n, 2n, [1n, 0n, 0n, 1n]);
  constructed.set(0, 1, new Rational(5n, 3n));
  expect(constructed.list()).toEqual([1, 1, 0, 1]);
  const cachedRow = constructed.row(0, true);
  expect(constructed.row(0, true)).toBe(cachedRow);
  expect(constructed.row(0)).not.toBe(cachedRow);
  expect(() => {
    cachedRow[0] = 0;
  }).toThrow('vector is immutable; please change a copy instead (use copy())');
  constructed.set(0, 0, 0n);
  expect(constructed.row(0, true)).not.toBe(cachedRow);
});

test('LLM.md — seeded binary matrices', async () => {
  const { random_matrix_gf2 } = await import('sagemath-ts/matrix');
  const { set_random_seed } = await import('sagemath-ts/misc/randstate');
  set_random_seed(42n);
  const firstRandom = random_matrix_gf2(2n, 65n);
  set_random_seed(42n);
  expect(firstRandom.eq(random_matrix_gf2(2n, 65n))).toBe(true);
  expect(random_matrix_gf2(2, 2, 1).list()).toEqual([1, 1, 1, 1]);
});

test('LLM.md — binary density inputs and PNG errors', async () => {
  const { Matrix_mod2_dense, random_matrix_gf2, to_png_data } = await import('sagemath-ts/matrix');
  const densityInput = new Matrix_mod2_dense(1, 2, [[1, 0]]);
  densityInput.randomize('0');
  expect(densityInput.list()).toEqual([1, 0]);
  expect(random_matrix_gf2(1, 2, new TextEncoder().encode('1')).list()).toEqual([1, 1]);
  expect(() => densityInput.randomize(null)).toThrow(
    "float() argument must be a string or a real number, not 'NoneType'"
  );
  expect(() => to_png_data(new Matrix_mod2_dense(0, 0))).toThrow(
    'cannot write image with dimensions 0 x 0'
  );
});

test('LLM.md — lattice Gaussian defaults and deferred centers', async () => {
  const { DiscreteGaussianDistributionLatticeSampler: DGL, DiscreteGaussianLattice } = await import(
    'sagemath-ts/stats'
  );
  const gaussian = new DGL([[1]]);
  expect(gaussian.sigma()).toBe(1);
  expect(DiscreteGaussianLattice([[1]]).sigma()).toBe(1);
  gaussian.set_c(null);
  expect(gaussian.c()).toBeNull();
  expect(gaussian.cNumeric()).toBeNull();
  expect(new DGL([])._call()).toBe(0);
  expect(new DGL([]).sampleExact()).toEqual([]);
});

test('LLM.md — Gaussian parameter strings and native double conversion', async () => {
  const {
    DiscreteGaussianDistributionIntegerSampler: DGI,
    DiscreteGaussianDistributionLatticeSampler: DGL,
  } = await import('sagemath-ts/stats');
  const { mpfr_init2, mpfr_set_d, mpfr_get_d } = await import('@sagemath-ts/mpfr-ts');
  expect(new DGI({ sigma: 1 / 128, c: -0 }).repr()).toBe(
    'Discrete Gaussian sampler over the Integers with sigma = 0.007812 and c = -0.000000'
  );
  expect(new DGL([[2]], { sigma: 1e6 }).repr().split(', c=')[0]).toBe(
    'Discrete Gaussian sampler with Gaussian parameter σ = 1.00000000000000e6'
  );
  const rounded = mpfr_init2(2);
  expect(mpfr_set_d(rounded, 1.25)).toBe(-1);
  expect(mpfr_get_d(rounded)).toBe(1);
});

test('LLM.md — rounded Gaussian centers retain exact integer bounds', async () => {
  const { DiscreteGaussianDistributionIntegerSampler: DGI } = await import('sagemath-ts/stats');
  const shifted = new DGI({ sigma: 2, c: 9007199254740993n, tau: 1 });
  expect(shifted.c).toBe(9007199254740992);
  expect(shifted.lowerBound).toBe(9007199254740990n);
  expect(shifted.upperBound).toBe(9007199254740994n);
});

test('LLM.md — native real numeric construction and MPFR integer setter', async () => {
  const { RealField } = await import('sagemath-ts/rings');
  const { mpfr_init2, mpfr_set_z, mpfr_get_d } = await import('@sagemath-ts/mpfr-ts');
  expect(new RealField(2).__call__(1.25).exact_rational()).toEqual([1n, 1n]);
  expect(new RealField(100).__call__(9007199254740993n).exact_rational()).toEqual([
    9007199254740993n,
    1n,
  ]);
  const integer = mpfr_init2(2);
  expect(mpfr_set_z(integer, 7n)).toBe(1);
  expect(mpfr_get_d(integer)).toBe(8);
});

test('LLM.md — native real predicates and integer rounding', async () => {
  const { RealField } = await import('sagemath-ts/rings');
  const { mpfr_init2, mpfr_set_d, mpfr_floor, mpfr_get_z } = await import('@sagemath-ts/mpfr-ts');
  const R = new RealField(128);
  expect(R.__call__('1e400').is_infinity()).toBe(false);
  expect(R.__call__('-1e-400').sign()).toBe(-1);
  expect(R.__call__('9007199254740993.25').floor()).toBe(9007199254740993n);
  const source = mpfr_init2(53),
    rounded = mpfr_init2(2);
  mpfr_set_d(source, 1.5);
  expect(mpfr_floor(rounded, source)).toBe(-2);
  expect(mpfr_get_z(source)).toEqual([2n, 2]);
});

test('LLM.md — native fractional parts and real comparison coercion', async () => {
  const { RealField } = await import('sagemath-ts/rings');
  const { mpfr_init2, mpfr_set_str, mpfr_frac, mpfr_get_d } = await import('@sagemath-ts/mpfr-ts');
  const R = new RealField(128);
  expect(R.__call__('9007199254740993.25').frac().exact_rational()).toEqual([1n, 4n]);
  const nearOne = R.__call__('1.00000000000000000001');
  expect(nearOne.equals(R.__call__(1))).toBe(false);
  expect(nearOne.equals(1)).toBe(true);
  const source = mpfr_init2(17),
    fraction = mpfr_init2(1);
  mpfr_set_str(source, '1.75');
  expect(mpfr_frac(fraction, source)).toBe(2);
  expect(mpfr_get_d(fraction)).toBe(1);
});

test('LLM.md — Gaussian copies delegate parameter coercion', async () => {
  const { DiscreteGaussianDistributionIntegerSampler: DGI } = await import('sagemath-ts/stats');
  const original = new DGI({ sigma: 2 });
  expect(original.withOptions({ c: 0.5 }).c).toBe(0.5);
  expect(original.withOptions({ c: 9007199254740993n }).c).toBe(9007199254740992);
  expect(original.c).toBe(0);
});

test('LLM.md — prime iterator bounds and successive differences', async () => {
  const { primes, differences, subfactorial } = await import('sagemath-ts/arith');
  expect([...primes(13n)]).toEqual([2n, 3n, 5n, 7n, 11n]);
  expect([...primes()]).toEqual([]);
  expect(primes(10n, Infinity).next().value).toBe(11n);
  expect(differences([1n, 4n, 9n], 3n)).toEqual([]);
  expect(subfactorial(8n)).toBe(14833n);
  expect(subfactorial(new Integer(8n))).toBe(14833n);
});

test('LLM.md — Dedekind algorithm selection', async () => {
  const { dedekind_sum } = await import('sagemath-ts/arith');
  expect(dedekind_sum(5n, -7n)).toEqual({ numerator: 0n, denominator: 1n });
  expect(dedekind_sum(5n, -7n, 'pari')).toEqual({ numerator: 1n, denominator: 14n });
});

test('LLM.md — integer rounding dispatch', async () => {
  const { integer_floor, integer_ceil, integer_trunc } = await import('sagemath-ts/arith');
  expect(integer_floor(new Rational(-5n, 2n))).toBe(-3n);
  expect(integer_ceil('2.5')).toBe(3n);
  expect(integer_trunc(new Rational(-5n, 2n))).toBe(-2n);
});

test('LLM.md — algebraic dependency parent dispatch and bounds', async () => {
  const { algebraic_dependency } = await import('sagemath-ts/arith');
  const { RealField } = await import('sagemath-ts/rings');
  const R = new RealField(128);
  expect(algebraic_dependency(R.__call__('1.5'), 2n)).toEqual([-3n, 2n]);
  expect(algebraic_dependency(5n, 0n, { height_bound: 5n })).toBeNull();
  expect(algebraic_dependency(Math.sqrt(2), 3n)).toEqual([-2n, 0n, 1n]);
});

test('LLM.md — native real addition and multiplication', async () => {
  const { mpfr_init2, mpfr_set_str, mpfr_add, mpfr_mul, mpfr_get_d } = await import(
    '@sagemath-ts/mpfr-ts'
  );
  const a = mpfr_init2(100),
    b = mpfr_init2(100),
    result = mpfr_init2(2);
  mpfr_set_str(a, '1.5');
  mpfr_set_str(b, '1.5');
  expect(mpfr_mul(result, a, b)).toBe(-1);
  expect(mpfr_get_d(result)).toBe(2);
  expect(mpfr_add(result, a, b)).toBe(0);
  expect(mpfr_get_d(result)).toBe(3);
});

test('LLM.md — IntegerLike prime-power traversal and ranges', async () => {
  const { Integer } = await import('sagemath-ts');
  const { next_prime_power, is_pseudoprime_power, prime_powers, eratosthenes } = await import(
    'sagemath-ts/arith'
  );
  expect(next_prime_power(new Integer(255n))).toBe(256n);
  expect(is_pseudoprime_power(new Integer(1000003n ** 6n), true)).toEqual([1000003n, 6n]);
  expect(prime_powers(new Integer(8n), new Integer(17n))).toEqual([8n, 9n, 11n, 13n, 16n]);
  expect(eratosthenes(new Integer(11n))).toEqual([2n, 3n, 5n, 7n, 11n]);
});

test('LLM.md — rational Hilbert symbols and wrapped conductor inputs', async () => {
  const { Integer, Rational } = await import('sagemath-ts');
  const { hilbert_symbol, hilbert_conductor, hilbert_conductor_inverse } = await import(
    'sagemath-ts/arith'
  );
  expect(hilbert_symbol(new Rational(1n, 2n), new Rational(1n, 3n), 2n, 'all')).toBe(-1n);
  expect(hilbert_conductor(new Integer(-3n), new Integer(-17n))).toBe(17n);
  expect(hilbert_conductor_inverse(new Integer(30n))).toEqual([-3n, -10n]);
});

test('LLM.md — arithmetic scalar edges and subset order', async () => {
  const { Integer } = await import('sagemath-ts');
  const {
    continuant,
    quadratic_residues,
    fundamental_discriminant,
    squarefree_divisors,
    odd_part,
    prime_to_m_part,
  } = await import('sagemath-ts/arith');
  expect(continuant([new Integer(7n), 2n], -1n)).toBe(7n);
  expect(quadratic_residues(new Integer(-8n))).toEqual([0n, 1n, 4n]);
  expect(fundamental_discriminant(0n)).toBe(0n);
  expect([...squarefree_divisors(30n)]).toEqual([1n, 2n, 3n, 6n, 5n, 10n, 15n, 30n]);
  expect(odd_part(new Integer(0n))).toBe(0n);
  expect(prime_to_m_part(new Integer(-240n), new Integer(-3n))).toBe(-80n);
});

test('LLM bounded arithmetic factories', async () => {
  const { get_gcd, get_inverse_mod, gcd } = await import('sagemath-ts/arith');
  expect(get_gcd(new Integer(4000n))(new Integer(-18n), 24n)).toBe(6n);
  expect(get_inverse_mod(6000n)(2n, -3n)).toBe(2n);
  expect(get_inverse_mod(600000n)(2n, -3n)).toBe(-4n);
  expect(get_gcd(4000000000n) === gcd).toBe(true);
});

test('LLM CRT list and basis dispatch', async () => {
  const { crt, CRT_list, CRT_basis } = await import('sagemath-ts/arith');
  expect(crt([2n, 3n], [5n, 7n])).toBe(17n);
  expect(CRT_list([10n], [3n])).toBe(10n);
  const residue = Mod(2n, 3n);
  expect(CRT_list([residue]) === residue).toBe(true);
  expect(CRT_list([residue, Mod(3n, 5n)]).modulus).toBe(15n);
  expect(CRT_basis([-3n])).toEqual([-2n]);
  expect(CRT_basis([7n, 6n, 10n], false)).toEqual([[120n, 120n, -140n, 21n], false]);
});

test('LLM product trees and factor bases', async () => {
  const { ProductTree, prod_with_derivative } = await import('sagemath-ts/rings');
  const { smooth_part, coprime_part } = await import('sagemath-ts/arith');
  const tree = new ProductTree([2n, 3n, 5n]);
  expect(tree.root()).toBe(30n);
  expect(tree.remainders(17n)).toEqual([1n, 2n, 2n]);
  expect(tree.interpolation([1n, 2n, 2n])).toBe(17n);
  expect(smooth_part(240n, tree)).toEqual([
    [2n, 4n],
    [3n, 1n],
    [5n, 1n],
  ]);
  expect(coprime_part(-126n, [6n, 14n])).toBe(-2n);
  expect(
    prod_with_derivative([
      [2n, 1n],
      [3n, 1n],
    ])
  ).toEqual([6n, 5n]);
});

test('LLM signed extended LCM and factor-based scalars', async () => {
  const { xlcm, dedekind_psi, carmichael_lambda } = await import('sagemath-ts/arith');
  expect(xlcm(new Integer(-12n), 18n)).toEqual([-36n, -4n, 9n]);
  expect(xlcm(12n, 0n)).toEqual([0n, 1n, 0n]);
  expect(dedekind_psi(new Integer(-6n))).toBe(-12n);
  expect(carmichael_lambda(new Integer(16n))).toBe(4n);
});

test('LLM native totient and divisor count', async () => {
  const { eulerphi, numdiv } = await import('@sagemath-ts/parigp-ts');
  expect(eulerphi(-12n)).toBe(4n);
  expect(eulerphi(0n)).toBe(2n);
  expect(numdiv(-12n)).toBe(6n);
});

test('LLM square decompositions and count coercion', async () => {
  const { two_squares, three_squares, sum_of_k_squares } = await import('sagemath-ts/arith');
  const { is_sum_of_two_squares_pyx } = await import('sagemath-ts/rings');
  expect(two_squares(new Integer(106n))).toEqual([5n, 9n]);
  expect(three_squares(new Integer(0n))).toEqual([0n, 0n, 0n]);
  expect(sum_of_k_squares(2n, 9634n)).toEqual([15n, 97n]);
  expect(sum_of_k_squares(2.9, 9634n)).toEqual([15n, 97n]);
  expect(is_sum_of_two_squares_pyx(21n)).toBe(false);
  expect(() => two_squares(21n)).toThrow('21 is not a sum of 2 squares');
});

test('LLM maximal-quotient reconstruction source behavior', async () => {
  const { mqrr_rational_reconstruction } = await import('sagemath-ts/arith');
  expect(mqrr_rational_reconstruction(new Integer(0n), 100n, 50n)).toEqual([0n, 1n]);
  expect(mqrr_rational_reconstruction(15n, 31n, 2n)).toEqual([15n, 1n]);
  expect(mqrr_rational_reconstruction(16n, 31n, 2n)).toBeNull();
});

test('LLM complex display records', async () => {
  const { sort_complex_numbers_for_display } = await import('sagemath-ts/arith');
  const nonreal = [{ re: -1, im: 1e-15 }, 'nonreal'] as const;
  const real = [{ re: 2, im: 0 }, 'real'] as const;
  const result = sort_complex_numbers_for_display([nonreal, real]);
  expect(result[0]).toBe(real);
  expect(result[1]).toBe(nonreal);
});

test('LLM delegated valuations and bit metadata', async () => {
  const { valuation } = await import('sagemath-ts');
  expect(valuation(0n, 0n)).toBe('Infinity');
  expect(valuation({_integer_: () => 40n}, 2n)).toBe(3n);
  const object = {_integer_: () => 0n, valuation: (_p: unknown) => [3n]};
  const result: bigint[] = valuation(object, 2n);
  expect(result).toEqual([3n]);
  expect(valuation(new Rational(5n, 9n), 3n)).toBe(-2n);
  const [exponent, unit] = new Rational(-4n, 17n).val_unit(new Integer(2n));
  expect(exponent).toBe(2n);
  expect(unit.eq(new Rational(-1n, 17n))).toBe(true);
  const integer = new Integer((1n << 262144n) - 1n);
  expect(integer.nbits()).toBe(262144n);
  expect(integer.bit_length()).toBe(262144n);
  expect(new Integer(integer).nbits()).toBe(262144n);
});

test('LLM Integer native factorial', () => {
  expect(new Integer(25n).factorial().value).toBe(15511210043330985984000000n);
  expect(() => new Integer(-1n).factorial()).toThrow('factorial only defined for nonnegative integers');
  expect(() => new Integer(1n << 64n).factorial()).toThrow('argument too large for factorial');
});

test('LLM PARI factorial real return', async () => {
  const { factorial } = await import('sagemath-ts/arith');
  const f = factorial(5n, 'pari');
  expect([f.s, f.e, f.m, f.p]).toEqual([1, 6n, 17293822569102704640n, 64]);
  expect(factorial(1n << 54n, 'pari').e).toBe(946788236117799971n);
  expect(factorial(5n)).toBe(120n);
  expect(() => factorial(1n << 63n, 'pari')).toThrow('Python int too large to convert to C long');
  expect(() => factorial(-1n, 'pari')).toThrow('factorial -- must be nonnegative');
});


test('LLM PARI unsigned word division', async () => {
  const { divru, real_1 } = await import('@sagemath-ts/parigp-ts');
  expect(divru(real_1(64), (1n << 64n) - 1n)).toEqual({
    s: 1, e: -64, p: 64, m: 9223372036854775809n,
  });
});


test('LLM PARI allocated zero precision', async () => {
  const { itor, rtor, real_0_bit } = await import('@sagemath-ts/parigp-ts');
  expect(itor(0n, 128)).toEqual({ s: 0, e: -128, m: 0n, p: 128 });
  expect(rtor(real_0_bit(-179), 64)).toEqual({ s: 0, e: -179, m: 0n, p: 64 });
});


test('LLM PARI product identity dispatch', async () => {
  const { mulrr } = await import('@sagemath-ts/parigp-ts');
  const x = { s: 1 as const, e: 0, p: 832, m: (5n << 829n) - 3n };
  expect(mulrr(x, x).m).toBe(mulrr(x, { ...x }).m + 1n);
});


test('LLM PARI real and integer square roots', async () => {
  const { itor, sqrtr, sqrti } = await import('@sagemath-ts/parigp-ts');
  expect(sqrtr(itor(-4n, 64))).toEqual({ re: 0n, im: { s: 1, e: 1, p: 64, m: 9223372036854775808n } });
  expect(sqrti(-9n)).toBe(3n);
});


test('LLM PARI real integer conversion error', async () => {
  const { gcvtoi, itor } = await import('@sagemath-ts/parigp-ts');
  expect(gcvtoi(itor(1n, 64))).toEqual([1n, -63]);
});


test('LLM native truncation and error-estimating conversion have distinct precision rules', async () => {
  const { gcvtoi, itor, truncr, PariError } = await import('@sagemath-ts/parigp-ts');
  const x = itor(1n << 64n, 64);
  expect(gcvtoi(x)).toEqual([18446744073709551616n, 1]);
  expect(() => truncr(x)).toThrow(PariError);
  expect(() => truncr(x)).toThrow('precision too low in truncr (precision loss in truncation)');
});


test('LLM public exponential preserves native exponent type and supports double conversion', async () => {
  const { itor, mpexp, rtodbl } = await import('@sagemath-ts/parigp-ts');
  const x = mpexp(itor(0n, 64));
  expect(x).toEqual({s:1,e:0n,m:9223372036854775808n,p:64});
  expect(rtodbl(x)).toBe(1);
});


test('LLM rational function field renaming maps and factory identity', async () => {
  const { FunctionField } = await import('sagemath-ts/rings/function_field');
  type ConstantField = import('sagemath-ts/rings/function_field').ConstantField<import('sagemath-ts/rings/function_field').ConstantFieldElement>;
  const K = FunctionField(QQ as unknown as ConstantField, 'x');
  const [L, fromL, toL] = K.change_variable_name('y');
  expect(toL(K.gen()).toString()).toBe('y');
  expect(fromL(L.gen()).toString()).toBe('x');
  expect(L.change_variable_name('x')[0]).toBe(K);
});


test('LLM function field underlying fraction and string construction', async () => {
  const {FunctionField} = await import('sagemath-ts/rings/function_field');
  type ConstantField = import('sagemath-ts/rings/function_field').ConstantField<import('sagemath-ts/rings/function_field').ConstantFieldElement>;
  const K=FunctionField(QQ as unknown as ConstantField,'x');
  const underlying=K.field(),reciprocal=K.__call__('1/x');
  expect(reciprocal.element().parent===underlying).toBe(true);
  expect(K.__call__(reciprocal.element()).element()===reciprocal.element()).toBe(true);
  expect(underlying.__call__('x==x','1/x').toString()).toBe('x');
  expect(reciprocal.pow(1n)===reciprocal).toBe(true);
});


test('LLM order basis, modular inverse and monoid identity', async () => {
  const {FunctionField} = await import('sagemath-ts/rings/function_field');
  type ConstantField = import('sagemath-ts/rings/function_field').ConstantField<import('sagemath-ts/rings/function_field').ConstantFieldElement>;
  const K=FunctionField(QQ as unknown as ConstantField,'x'),finiteOrder=K.maximal_order();
  expect(K.maximal_order_infinite().basis()).toEqual([1n]);
  expect(K.__call__(2n).inverse_mod(finiteOrder.ideal(K.zero())).toString()).toBe('1/2');
  const ideal=finiteOrder.ideal(K.gen());
  expect(finiteOrder.ideal_monoid()===finiteOrder.ideal_monoid()).toBe(true);
  expect(finiteOrder.ideal_monoid().__call__(ideal)===ideal).toBe(true);
});


test('LLM native FLINT word arithmetic', async () => {
  // Resolve the public package name from its declared workspace consumer.
  const {createRequire}=await import('node:module');
  const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
  const {n_sqrtmod,n_jacobi,n_is_square,n_preinvert_limb,n_powmod2_ui_preinv}=await import(require.resolve('@sagemath-ts/flint-ts'));
  expect(n_sqrtmod(100n,101n)).toBe(10n);
  expect(n_jacobi(-1n,7n)).toBe(-1);
  expect(n_is_square(4294967295n**2n)).toBe(true);
  expect(n_powmod2_ui_preinv(100n,3n,101n,n_preinvert_limb(101n))).toBe(100n);
});

test('LLM native polynomial roots, series and NTL squarefree kernels', async()=>{
  const {createRequire}=await import('node:module');
  const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
  const {_nmod_poly_sqrt,_gr_poly_sqrt_series_newton}=await import(require.resolve('@sagemath-ts/flint-ts'));
  const {ZZX_GCD,ZZX_SquareFreeDecomp}=await import(require.resolve('@sagemath-ts/ntl-ts'));
  expect(_nmod_poly_sqrt([4n,8n,4n],17n)).toEqual([2n,2n]);
  expect(_gr_poly_sqrt_series_newton([1n,2n,1n],8,2,17n)).toEqual([1n,1n]);
  expect(ZZX_GCD([-2n,0n,2n],[2n,-4n,2n])).toEqual([-2n,2n]);
  expect(ZZX_SquareFreeDecomp([1n,4n,4n])).toEqual([[[1n,2n],2]]);
});

test('LLM polynomial and underlying fraction square roots',async()=>{
  const {FunctionField}=await import('sagemath-ts/rings/function_field');
  const {FractionFieldElement}=await import('sagemath-ts/rings/fraction_field_element');
  type ConstantField=import('sagemath-ts/rings/function_field').ConstantField<import('sagemath-ts/rings/function_field').ConstantFieldElement>;
  const K=FunctionField(QQ as unknown as ConstantField,'x'),F=K.field(),R=K._ring;
  const raw=new FractionFieldElement(F,R.gen(),R.gen(),{reduce:false});
  expect(raw.is_square()).toBe(true);
  expect(String(raw.sqrt())).toBe('1');
  expect(String((R.__call__([1n,2n,1n]).is_square(true))[1])).toBe('x + 1');
  const zero=K.zero();
  expect(zero.sqrt()===zero).toBe(false);
  expect((F.gen().sqrt(false,true) as unknown[])).toEqual([]);
});

test('LLM squarefree decomposition retains native units and factor order',async()=>{
  const {createRequire}=await import('node:module');
  const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
  const {nmod_poly_factor_squarefree}=await import(require.resolve('@sagemath-ts/flint-ts'));
  const {PolynomialRing}=await import('sagemath-ts/rings/polynomial');
  expect(nmod_poly_factor_squarefree([0n,1n,0n,0n,0n,0n,1n],5n)).toEqual([[[0n,1n],1],[[1n,1n],5]]);
  const Qx=new PolynomialRing(QQ,'x');
  expect(Qx.__call__([0n,2n]).squarefree_decomposition().map(([f,e])=>[f.toString(),e])).toEqual([['2*x',1]]);
});

test('LLM binary coefficient roots and cached elements',async()=>{
  const {GF2,GF2Element}=await import('sagemath-ts/rings/finite_rings');
  const bit=GF2.__call__(1n);
  expect(bit.is_square()).toBe(true);
  expect(bit.sqrt()).toBe(bit);
  expect(bit.sqrt({all:true})).toEqual([bit]);
  expect(bit.add(1n)).toBe(GF2.zero());
  const rawZero=new GF2Element(0n,GF2);
  expect(rawZero.neg()).toBe(rawZero);
  expect(rawZero.sqrt()).toBe(GF2.zero());
});


test('LLM PARI multiplicity-indexed squarefree components',async()=>{
  const {createRequire}=await import('node:module');
  const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
  const {FpX_factor_squarefree}=await import(require.resolve('@sagemath-ts/parigp-ts'));
  expect(FpX_factor_squarefree([2n,3n,1n],101n)).toEqual([[2n,3n,1n]]);
  expect(FpX_factor_squarefree([1n,2n,1n],101n)).toEqual([[1n],[1n,1n]]);
  expect(FpX_factor_squarefree([2n,2n],101n)).toEqual([[2n,2n]]);
  expect(FpX_factor_squarefree([2n],101n)).toEqual([]);
  expect(FpX_factor_squarefree([2n],(1n<<127n)-1n)).toEqual([[2n]]);
});


test('LLM polynomial numerator parents and LCM backends',async()=>{
  const {PolynomialRing}=await import('sagemath-ts/rings/polynomial');
  const T=new PolynomialRing(QQ,'t'),f=T.__call__([QQ.__call__(1n).div(QQ.__call__(3n)),QQ.__call__(1n).div(QQ.__call__(2n))]);
  expect(f.denominator()).toBe(6n);
  expect(String(f.numerator())).toBe('3*t + 2');
  expect(String(f.numerator().parent.base_ring)).toBe('Integer Ring');
  expect(f.numerator().parent).toBe(f.numerator().parent);
  expect(String(f.lcm(T.__call__([1n,1n])))).toBe('t^2 + 5/3*t + 2/3');
  const {createRequire}=await import('node:module');
  const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
  const {_fmpz_poly_lcm,_fmpq_poly_lcm,fmpq_poly_get_numerator,fmpq_poly_get_denominator}=await import(require.resolve('@sagemath-ts/flint-ts'));
  expect(_fmpz_poly_lcm([-2n,2n],[3n,3n])).toEqual([-6n,0n,6n]);
  expect(_fmpq_poly_lcm([1n,2n],[1n,2n])).toEqual([[1n,2n],2n]);
  expect(fmpq_poly_get_numerator([-1n,0n,3n],2n)).toEqual([-1n,0n,3n]);
  expect(fmpq_poly_get_denominator([-1n,0n,3n],2n)).toBe(2n);
});


test('LLM divisor formatter, parent coercion and cached spaces',async()=>{
  const {FunctionField,DivisorGroup}=await import('sagemath-ts/rings/function_field');
  type CF=import('sagemath-ts/rings/function_field').ConstantField<import('sagemath-ts/rings/function_field').ConstantFieldElement>;
  const K=FunctionField(QQ as unknown as CF,'x'),DG=K.divisor_group(),P=K.maximal_order().ideal(K.gen()).place();
  expect(new DivisorGroup(K)).toBe(DG);
  expect(String(DG.__call__(P))).toBe('Place (x)');
  const D=P.divisor(0n);
  expect(D._format(v=>`<${v}>`,' @ ',' / ')).toBe('<0> @ <Place (x)>');
  expect(D.add(DG.zero()).support()).toHaveLength(0);
  expect(DG.zero().function_space()).toBe(DG.zero().function_space());
});


test('LLM ideal monoid identities and scalar membership',async()=>{
  const {FunctionField}=await import('sagemath-ts/rings/function_field');
  type CF=import('sagemath-ts/rings/function_field').ConstantField<import('sagemath-ts/rings/function_field').ConstantFieldElement>;
  const K=FunctionField(QQ as unknown as CF,'x'),O=K.maximal_order(),I=O.ideal(K.gen());
  expect(I.parent()).toBe(O.ideal_monoid());
  expect(I.pow(1n)).toBe(I);
  expect(I.pow(0n)).toBe(I.parent().one());
  expect(O.ideal(K.one()).contains(1n)).toBe(true);
});

test('LLM polynomial enumeration and rational cardinality examples', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings/polynomial');
  const { GF } = await import('sagemath-ts/rings/finite_rings');
  const { QQ } = await import('sagemath-ts/rings');
  const E = new PolynomialRing(GF(3n), 'x');
  expect([...E.monics({ max_degree: 1 })].map(String)).toEqual(['1', 'x', 'x + 1', 'x + 2']);
  expect([...E.polynomials({ max_degree: -1 })].map(String)).toEqual(['0']);
  expect(QQ.cardinality()).toBe('Infinity');
});

test('LLM immutable function-field matrix and mutable copy', async () => {
  const { QQ } = await import('sagemath-ts/rings');
  const { FunctionField } = await import('sagemath-ts/rings/function_field');
  const F = FunctionField(QQ as never, 'x');
  const M = F.gen().matrix();
  expect(Object.isFrozen(M) && Object.isFrozen(M[0])).toBe(true);
  const C = M.map(row => row.slice());
  C[0]![0] = F.one();
  expect(String(M[0]![0])).toBe('x');
});

test('LLM function-field element predicate parent protocol', async () => {
  const { FunctionField, is_FunctionFieldElement } = await import('sagemath-ts/rings/function_field');
  const { QQ } = await import('sagemath-ts/rings');
  const K = FunctionField(QQ as never, 'x');
  expect(is_FunctionFieldElement({ parent: () => K })).toBe(true);
  expect(is_FunctionFieldElement(1n)).toBe(false);
});


test('LLM matrix base validation and parent cache', async () => {
  const { RationalFunctionField } = await import('sagemath-ts/rings/function_field');
  const { QQ } = await import('sagemath-ts/rings');
  const K = new RationalFunctionField(QQ as never, 'documentedmatrixbase');
  expect(() => K.gen().matrix(QQ)).toThrow('base must be the rational function field itself');
  K.one().norm();
  expect(String(K.gen().matrix(QQ)[0]![0])).toBe('documentedmatrixbase');
  expect(() => K.gen().matrix([])).toThrow("unhashable type: 'list'");
});


test('LLM variable-map domain conversion and polynomial dictionaries', async () => {
  const { FunctionField } = await import('sagemath-ts/rings/function_field');
  const { PolynomialRing } = await import('sagemath-ts/rings/polynomial');
  const { QQ } = await import('sagemath-ts/rings');
  const K = FunctionField(QQ as never, 'x');
  const [L, , toL] = K.change_variable_name('y');
  expect(toL(1n).parent).toBe(L);
  expect(String(new PolynomialRing(QQ, 'x').__call__({ 2: 3n, 0: -2n }))).toBe('3*x^2 - 2');
});

test('LLM FLINT deterministic word random state', async () => {
  const { createRequire } = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const { flint_rand_init, n_randlimb } = await import(require.resolve('@sagemath-ts/flint-ts'));
  expect(n_randlimb(flint_rand_init())).toBe(15737102703946861599n);
});

test('LLM FLINT word-prime factorization', async () => {
  const { createRequire } = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const { nmod_poly_factor } = await import(require.resolve('@sagemath-ts/flint-ts'));
  expect(nmod_poly_factor([4n,0n,0n,0n,1n],5n)).toEqual([1n,[[[4n,1n],1],[[2n,1n],1],[[3n,1n],1],[[1n,1n],1]]]);
});

test('LLM PARI native random-state save and restore', async () => {
  const { setrand, getrand, pari_rand } = await import('@sagemath-ts/parigp-ts');
  setrand(1n);
  const saved = getrand();
  expect(pari_rand()).toBe(13282407956253574712n);
  setrand(saved);
  expect(pari_rand()).toBe(13282407956253574712n);
});

test('LLM PARI binary polynomial factors and kernel mutation', async () => {
  const { F2x_mul, F2x_factor, F2m_ker_sp } = await import('@sagemath-ts/parigp-ts');
  expect(F2x_mul(7n, 3n)).toBe(9n);
  expect(F2x_factor(63n)).toEqual([[3n,1],[7n,2]]);
  const columns = [3n,3n];
  expect(F2m_ker_sp(columns,2)).toEqual([3n]);
  expect(columns).toEqual([3n,1n]);
});


test('LLM PARI quotient power vectors preserve initial entries', async () => {
  const { FpXQ_powers, FpXQ_autpowers } = await import('@sagemath-ts/parigp-ts');
  expect(FpXQ_powers([1n, 1n], 2, [1n], 3n)).toEqual([[1n], [1n, 1n], []]);
  expect(FpXQ_autpowers([1n, 1n], 1, [1n], 3n)).toEqual([[], [0n, 1n], [1n, 1n]]);
});


test('LLM PARI column matrix products', async () => {
  const { ZM_mul, FpM_mul, F2m_mul } = await import('@sagemath-ts/parigp-ts');
  const A = [[], [0n, 1n, 3n], [0n, 2n, 4n]];
  expect(ZM_mul(A, A)).toEqual([[], [0n, 7n, 15n], [0n, 10n, 22n]]);
  expect(FpM_mul(A, A, 3n)).toEqual([[], [0n, 1n, 0n], [0n, 1n, 1n]]);
  expect(F2m_mul([3n, 2n], [3n, 1n], 2)).toEqual([1n, 3n]);
});


test('LLM PARI modular composition and automorphism trace', async () => {
  const { brent_kung_optpow, Flxq_powers, FpXQ_auttrace } = await import('@sagemath-ts/parigp-ts');
  expect(brent_kung_optpow(10, 2, 1)).toBe(5);
  expect(Flxq_powers([1n, 1n], 3, [1n, 0n, 1n], 5n)).toEqual([[1n], [1n, 1n], [0n, 2n], [3n, 2n]]);
  expect(FpXQ_auttrace([[0n, 1n], [1n, 1n]], 3n, [1n, 0n, 1n], 5n)).toEqual([[0n, 1n], [3n, 3n]]);
});


test('LLM PARI polynomial products and squares', async () => {
  const { ZX_sqr, Flx_mul, Flx_sqr, FpX_sqr } = await import('@sagemath-ts/parigp-ts');
  expect(ZX_sqr([1n, -1n, 1n])).toEqual([1n, -2n, 3n, -2n, 1n]);
  expect(Flx_mul([1n, 2n], [3n, 4n], 5n)).toEqual([3n, 0n, 3n]);
  expect(Flx_sqr([1n, 2n], 5n)).toEqual([1n, 4n, 4n]);
  expect(FpX_sqr([-1n, 1n], 3n)).toEqual([1n, 1n, 1n]);
});

test('LLM PARI polynomial division and reciprocal inverses', async () => {
  const { FpX_divrem, FpX_rem, Flx_divrem, Flx_rem, FpX_invBarrett, Flx_invBarrett } = await import('@sagemath-ts/parigp-ts');
  expect(FpX_divrem([1n, 2n, 3n], [1n, 1n], 5n)).toEqual([[4n, 3n], [2n]]);
  expect(FpX_rem([1n, 2n, 3n], [1n, 1n], 5n)).toEqual([2n]);
  expect(Flx_divrem([1n, 2n, 3n], [1n, 1n], 5n)).toEqual([[4n, 3n], [2n]]);
  expect(Flx_rem([1n, 2n], [2n], 8n)).toEqual([]);
  expect(FpX_invBarrett([1n, 2n, 3n, 1n], 5n)).toEqual([1n, 2n]);
  expect(Flx_invBarrett([1n, 2n, 3n, 1n], 5n)).toEqual([1n, 2n]);
});

test('LLM PARI unscaled GCD and row polynomial matrices', async () => {
  const { FpX_gcd, FpX_extgcd, Flx_gcd, Flx_extgcd, FpX_halfgcd, FpX_halfgcd_all, Flx_halfgcd, Flx_halfgcd_all, gen_pow_i } = await import('@sagemath-ts/parigp-ts');
  const a = [1n, 2n, 3n], b = [1n, 1n], M = [[[], [1n]], [[1n], [1n, 2n]]];
  expect(FpX_gcd(a,b,5n)).toEqual([2n]); expect(Flx_gcd(a,b,5n)).toEqual([2n]);
  expect(FpX_extgcd(a,b,5n)).toEqual([[2n],[1n],[1n,2n]]);
  expect(Flx_extgcd(a,b,5n)).toEqual([[2n],[1n],[1n,2n]]);
  expect(FpX_halfgcd(a,b,5n)).toEqual(M);expect(Flx_halfgcd(a,b,5n)).toEqual(M);
  expect(FpX_halfgcd_all(a,b,5n)).toEqual([M,[1n,1n],[2n]]);
  expect(Flx_halfgcd_all(a,b,5n)).toEqual([M,[1n,1n],[2n]]);
  expect(gen_pow_i(3n,3n,x=>x*x,(x,y)=>x*y)).toBe(27n);
});


test('LLM PARI Shoup minimal polynomials', async () => {
  const { FpXQ_minpoly, Flxq_minpoly, setrand, getrand } = await import('@sagemath-ts/parigp-ts');
  const saved = getrand();
  try {
    setrand(1n);
    expect(FpXQ_minpoly([0n, 1n], [1n, 0n, 1n], 7n)).toEqual([1n, 0n, 1n]);
    const consumed = getrand();
    setrand(1n);
    expect(Flxq_minpoly([0n, 1n], [1n, 0n, 1n], 7n)).toEqual([1n, 0n, 1n]);
    expect(getrand()).toBe(consumed);
  } finally { setrand(saved); }
});

test('LLM PARI distinct-degree factors and counts', async () => {
  const { FpX_ddf, Flx_ddf, FpX_nbfact, Flx_nbfact_by_degree } = await import('@sagemath-ts/parigp-ts');
  const f = [5n, 2n, 6n, 0n, 0n, 2n, 6n, 3n, 1n];
  expect(FpX_ddf(f, 7n).get(1)).toEqual([2n, 4n, 5n]);
  expect(Flx_ddf(f, 7n)[0]).toEqual([[2n, 4n, 5n], 1]);
  expect(FpX_nbfact(f, 7n)).toBe(4);
  expect(Flx_nbfact_by_degree(f, 7n)).toEqual({ D: [0, 2, 1, 0, 1, 0, 0, 0, 0], nb: 4 });
});

test('PARI modular square-root and fused-power examples', async () => {
  const { Fp_sqrt, Fp_sqrt_i, Fl_sqrt, gen_pow_fold } = await import('@sagemath-ts/parigp-ts');
  expect(Fp_sqrt(2n, 17n)).toBe(6n);
  expect(Fp_sqrt(3n, 17n)).toBeNull();
  expect(Fp_sqrt_i(2n, null, 17n)).toBe(6n);
  expect(Fl_sqrt(2n, 17n)).toBe(6n);
  expect(gen_pow_fold(3n, -5n, x => x * x, x => 3n * x * x)).toBe(243n);
});

test('LLM PARI full factors and monic normalization', async () => {
  const { FpX_factor, Flx_factor, FpX_normalize, Flx_normalize, getrand, setrand } =
    await import('@sagemath-ts/parigp-ts');
  const saved = getrand();
  try {
    setrand(1n);
    expect(FpX_factor([1n, 0n, 0n, 0n, 1n], 5n)).toEqual([
      [[2n, 0n, 1n], 1], [[3n, 0n, 1n], 1],
    ]);
    expect(Flx_factor([2n, 4n, 2n], 5n)).toEqual([[[1n, 1n], 2]]);
    expect(FpX_factor([], 7n)).toEqual([[[], 1]]);
    expect(FpX_normalize([2n, 4n, 2n], 5n)).toEqual([1n, 2n, 1n]);
    expect(Flx_normalize([2n, 4n, 2n], 5n)).toEqual([1n, 2n, 1n]);
  } finally {
    setrand(saved);
  }
});

test('LLM PARI polynomial roots and splitting predicates', async () => {
  const { FpX_roots, Flx_roots, Flx_nbroots, FpX_is_totally_split, Flx_is_totally_split } =
    await import('@sagemath-ts/parigp-ts');
  expect(FpX_roots([1n, 7n], 7n)).toEqual([]);
  expect(Flx_roots([6n, 0n, 1n], 7n)).toEqual([1n, 6n]);
  expect(Flx_nbroots([1n, 2n, 1n], 7n)).toBe(1);
  expect(Flx_is_totally_split([1n, 2n, 1n], 7n)).toBe(false);
  expect(FpX_is_totally_split([7n], 7n)).toBe(true);
  const p = (1n << 64n) - 59n;
  expect(FpX_roots([2n, p - 1n, p - 2n, 1n], p)).toEqual([p - 1n, 1n, 2n]);
});

test('LLM PARI extension-polynomial coefficient representations', async () => {
  const { FpXQX_sqr, FpXQX_red, FpXQX_normalize, FlxqX_mul, FlxqX_normalize,
    F2xqX_sqr, F2xqX_normalize } = await import('@sagemath-ts/parigp-ts');
  const T = [2n, 0n, 1n];
  expect(FpXQX_sqr([4n, 4n], T, 5n)).toEqual([16n, 32n, 16n]);
  expect(FpXQX_red([16n, 32n, 16n], T, 5n)).toEqual([1n, 2n, 1n]);
  expect(FpXQX_normalize([[1n, 0n, 1n], [1n]], T, 5n)).toEqual([[1n, 0n, 1n], 1n]);
  expect(FlxqX_normalize([[1n, 0n, 1n], [1n]], T, 5n)).toEqual([[4n], [1n]]);
  expect(FlxqX_mul([[1n, 1n], [1n]], [[1n, 1n], [1n]], T, 5n)).toEqual([[4n, 2n], [2n, 2n], [1n]]);
  expect(F2xqX_sqr([3n, 1n], 7n)).toEqual([2n, 0n, 1n]);
  expect(F2xqX_normalize([7n, 1n], 7n)).toEqual([0n, 1n]);
});

test('LLM PARI extension-polynomial division and cached reductions', async () => {
  const {FpXQX_divrem,FpXQX_rem,FlxqX_divrem,FlxqX_get_red,F2xqX_divrem} =
    await import('@sagemath-ts/parigp-ts');
  const T=[3n,0n,1n],x=[[1n],[],[1n]],S=[[1n],[1n]];
  expect(FpXQX_divrem(x,S,T,17n)).toEqual([[16n,1n],[2n]]);
  expect(FlxqX_divrem(x,S,T,17n)).toEqual([[[16n],[1n]],[[2n]]]);
  expect(FpXQX_rem([[1n]],S,T,17n)).toEqual([[1n]]);
  expect(F2xqX_divrem([1n,0n,1n],[1n,1n],7n)).toEqual([[1n,1n],[]]);
  const divisor=Array.from({length:16},()=>[1n]),cached=FlxqX_get_red(divisor,T,17n);
  expect(FlxqX_divrem(divisor,cached,T,17n)).toEqual([[[1n]],[]]);
});

test('LLM PARI extension GCD and native raw binary packing', async () => {
  const {FpXQX_gcd,FpXQX_extgcd,FlxqX_extgcd,F2xqX_halfgcd,F2xqX_mul}=
    await import('@sagemath-ts/parigp-ts');
  const T=[3n,0n,1n];
  expect(FpXQX_gcd([2n,2n],[4n,4n],T,17n)).toEqual([4n,4n]);
  expect(FpXQX_extgcd([1n,0n,1n],[1n,1n],T,17n)).toEqual([[2n],[1n],[1n,16n]]);
  expect(FlxqX_extgcd([[1n],[],[1n]],[[1n],[1n]],T,17n)).toEqual([[[2n]],[[1n]],[[1n],[16n]]]);
  expect(F2xqX_halfgcd([1n],[8n,1n],7n)).toEqual([[[1n],[]],[[1n,1n],[1n]]]);
  expect(F2xqX_mul([32n,1n],[1n],7n)).toEqual([]);
});

test('LLM PARI extension quotient inverses and power tables', async () => {
  const {FpXQXQ_inv,FpXQXQ_invsafe,FpXQXQ_powers,FlxqXQ_powu,F2xqXQ_inv,F2xqXQ_powers}=
    await import('@sagemath-ts/parigp-ts');
  const T=[3n,0n,1n],S=[1n,0n,1n];
  expect(FpXQXQ_inv([1n,1n],S,T,17n)).toEqual([9n,8n]);
  expect(FpXQXQ_invsafe(S,S,T,17n)).toBeNull();
  expect(FpXQXQ_powers([1n,1n],3,S,T,17n)).toEqual([[1n],[1n,1n],[0n,2n],[15n,2n]]);
  expect(FlxqXQ_powu([[1n],[1n]],3n,[[1n],[],[1n]],T,17n)).toEqual([[15n],[2n]]);
  expect(F2xqXQ_inv([2n,1n],[3n,1n,1n],7n)).toEqual([2n,3n]);
  expect(F2xqXQ_powers([2n,1n],3,[3n,1n,1n],7n)).toEqual([[1n],[2n,1n],[0n,1n],[3n,3n]]);
});

test('LLM PARI extension composition and packed matrices',async()=>{
 const api=await import('@sagemath-ts/parigp-ts');
 expect(api.FpXQX_FpXQXQ_eval([1n,1n,1n],[1n,1n],[1n,0n,1n],[3n,0n,1n],17n)).toEqual([2n,3n]);
 expect(api.FlxqX_FlxqXQ_eval([[1n],[1n],[1n]],[[1n],[1n]],[[1n],[],[1n]],[3n,0n,1n],17n)).toEqual([[2n],[3n]]);
 expect(api.F2xqX_F2xqXQ_eval([1n,2n,1n],[2n,1n],[3n,1n,1n],7n)).toEqual([2n,3n]);
 const A=[[],[[],[1n],[2n]],[[],[3n],[4n]]];
 expect(api.FlxqM_mul(A,A,[3n,0n,1n],17n)).toEqual([[],[[],[7n],[10n]],[[],[15n],[5n]]]);
});

test('LLM native power-table canonical storage and reciprocal errors',async()=>{
 const api=await import('@sagemath-ts/parigp-ts');
 expect(api.FpXQ_powers([18n,0n],1,[3n,0n,1n],17n)).toEqual([[1n],[18n]]);
 const T=[1n,...Array(35).fill(0n),2n];
 expect(api.FpX_FpXQV_eval([0n,0n],[],T,4n)).toEqual([]);
 expect(()=>api.FpXQ_powers([1n],0,T,4n)).toThrow('impossible inverse in Fp_inv: Mod(2, 4).');
});

test('LLM native automorphism storage and cache preparation', async () => {
  const api = await import('@sagemath-ts/parigp-ts');
  expect(api.FpXQ_autpowers([18n, 0n], 1, [3n, 0n, 1n], 17n)).toEqual([[], [0n, 1n], [18n]]);
  const T = [1n, ...Array(35).fill(0n), 2n];
  expect(() => api.FpXQ_auttrace([[1n], [1n]], 1n, T, 4n)).toThrow(
    'impossible inverse in Fp_inv: Mod(2, 4).'
  );
});

test('LLM coefficient substitution and binary scalar composition', async () => {
  const api = await import('@sagemath-ts/parigp-ts');
  expect(api.FpXY_FpXQ_evalx([18n, [18n, 0n]], [1n], [3n, 0n, 1n], 17n)).toEqual([18n, [1n]]);
  expect(api.FlxY_Flxq_evalx([[1n], [1n, 1n]], [2n], [3n, 0n, 1n], 17n)).toEqual([[1n], [3n]]);
  expect(api.F2xY_F2xq_evalx([3n, 5n], 2n, 7n)).toEqual([3n, 2n]);
  expect(api.F2x_F2xq_eval((1n << 4096n) - 1n, 2n, 7n)).toBe(1n);
});


test('LLM native extension automorphism tuple layouts and count casts', async () => {
  const { FpXQXQ_auttrace, FpXQXQ_autsum, F2xqXQ_auttrace } = await import('@sagemath-ts/parigp-ts');
  const T = [1n, 0n, 1n], S = [1n, 0n, 1n], B = [0n, 1n];
  expect(FpXQXQ_auttrace([B, [2n, 1n]], 3n, S, T, 17n)).toEqual([[0n,1n],[6n,3n]]);
  expect(FpXQXQ_autsum([[0n,1n], B, [2n,1n]], 3n, S, T, 17n)).toEqual([[0n,1n],[0n,1n],[2n,11n]]);
  expect(FpXQXQ_auttrace([B, [1n]], -1n, S, T, 17n)).toEqual([[0n,1n],[]]);
  expect(F2xqXQ_auttrace([2n, [0n,1n], [1n,1n]], 3n, [1n,0n,1n], 7n)).toEqual([2n,[0n,1n],[1n,1n]]);
});


test('LLM extension projections and truncated products', async () => {
  const api = await import('@sagemath-ts/parigp-ts');
  const T = [1n, 0n, 1n], saved = api.getrand();
  try {
    api.setrand(1n);
    expect(api.random_FpXQX(2, T, 17n)).toEqual([[13n,11n],[11n,3n]]);
    api.setrand(1n);
    expect(api.random_FlxqX(2, T, 17n)).toEqual([[13n,11n],[11n,3n]]);
    expect(api.FpXQX_dotproduct([1n,1n], [1n,-1n], T, 17n)).toBe(0n);
    expect(api.FlxqX_dotproduct([[1n],[2n]], [[3n],[4n]], T, 17n)).toEqual([11n]);
    expect(api.FpXQXn_sqr([18n,19n], 2, T, 17n)).toEqual([324n,684n]);
    expect(api.FpXQXn_mul([18n,19n], [18n,19n], 2, T, 17n)).toEqual([1n,4n]);
    expect(api.FlxqXn_mul([[1n],[2n]], [[1n],[2n]], 2, T, 17n)).toEqual([[1n],[4n]]);
    expect(api.FlxqXn_sqr([[1n],[2n]], 2, T, 17n)).toEqual([[1n],[4n]]);
  } finally { api.setrand(saved); }
});


test('LLM extension minimal polynomials', async () => {
  const api = await import('@sagemath-ts/parigp-ts'), saved = api.getrand();
  try {
    api.setrand(1n);
    expect(api.FpXQXQ_minpoly([0n,1n], [1n,0n,1n], [3n,0n,1n], 17n)).toEqual([[1n],[],1n]);
    api.setrand(1n);
    expect(api.FlxqXQ_minpoly([[],[1n]], [[1n],[],[1n]], [3n,0n,1n], 17n)).toEqual([[1n],[],[1n]]);
    api.setrand(1n);
    expect(api.FpXQXQ_minpoly([[1n,1n]], [1n,0n,1n], [3n,0n,1n], 17n)).toEqual([[16n,16n],1n]);
    api.setrand(1n);
    expect(api.FpXQXQ_minpoly([0n,1n], [1n,0n,0n,0n,1n], [1n,1n,1n], 2n)).toEqual([[1n],[],[],[],1n]);
  } finally { api.setrand(saved); }
});


test('LLM exact quotient cache and inverse errors', async () => {
  const { FpXQ_pow, FpXQ_inv, PariError } = await import('@sagemath-ts/parigp-ts');
  const T = [1n, ...Array<bigint>(90).fill(0n), 2n];
  expect(FpXQ_pow([1n], 0n, T, 4n)).toEqual([1n]);
  for (const [run, message] of [
    [() => FpXQ_pow([1n], 2n, T, 4n), 'impossible inverse in Fl_inv: Mod(2, 4).'],
    [() => FpXQ_inv([], [1n,0n,1n], 17n), 'impossible inverse in FpXQ_inv: 0.'],
    [() => FpXQ_pow([], -2n, [1n,0n,1n], 17n), 'impossible inverse in Flxq_inv: 0.'],
  ] as const) {
    expect(run).toThrow(PariError);
    expect(run).toThrow(message);
  }
});


test('LLM extension Frobenius APIs', async () => {
  const api = await import('@sagemath-ts/parigp-ts');
  const T = [3n,0n,1n], S = [1n,0n,1n], W = [[1n],[],[1n]];
  expect(api.Flx_Frobenius(T,17n)).toEqual([0n,16n]);
  expect(api.FpXQX_Frobenius(S,T,17n)).toEqual([0n,1n]);
  expect(api.FlxqX_Frobenius(W,T,17n)).toEqual([[],[1n]]);
  expect(api.F2xqX_Frobenius([1n,1n,1n],7n)).toEqual([0n,1n]);
  expect(api.FpXQXQ_halfFrobenius([1n,1n],S,T,17n)).toEqual([1n]);
  expect(api.FlxqXQ_halfFrobenius([[1n],[1n]],W,T,17n)).toEqual([[1n]]);
});


test('LLM prime-polynomial signed and zero modulus boundaries',async()=>{
  const api=await import('@sagemath-ts/parigp-ts');
  expect(api.FpX_red([-1n,18n,0n],-17n)).toEqual([16n,1n]);
  expect(api.FpX_add([0n,1n],[0n,-1n],0n)).toEqual([]);
  expect(api.FpX_Fp_mul([1n,1n,1n],0n,0n)).toEqual([]);
  expect(api.FpX_Fp_mul([0n,0n],1n,0n)).toEqual([]);
  expect(()=>api.FpX_red([1n],0n)).toThrow(api.PariError);
  expect(()=>api.FpX_red([1n],0n)).toThrow('impossible inverse in dvmdii: 0.');
});


test('LLM extension root counts, squarefree checks and derivatives',async()=>{
  const r=await import('@sagemath-ts/parigp-ts');
  const T=[1n,0n,1n],S=[1n,0n,1n],W=[[1n],[],[1n]];
  expect(r.FpXQX_split_part(S,T,3n)).toEqual([1n,0n,1n]);
  expect(r.FpXQX_nbroots(S,T,3n)).toBe(2);
  expect(r.FlxqX_nbroots(W,T,3n)).toBe(2);
  expect(r.F2xqX_nbroots([1n,1n,1n],7n)).toBe(2);
  expect(r.FqX_nbroots(S,null,3n)).toBe(0);
  expect(r.FqX_nbroots(S,T,3n)).toBe(2);
  expect(r.FpXQX_is_squarefree([1n,-2n,1n],T,3n)).toBe(false);
  expect(r.FlxqX_is_squarefree(W,T,3n)).toBe(true);
  expect(r.FpXX_deriv([1n,[2n,3n],4n],17n)).toEqual([[2n,3n],8n]);
  expect(r.FlxX_deriv([[1n],[2n,3n],[4n]],17n)).toEqual([[2n,3n],[8n]]);
});


test('LLM signed composition and matrix boundaries', async () => {
  const { FpXQ_powers, FpX_FpXQ_eval, FpM_mul } = await import('@sagemath-ts/parigp-ts');
  expect(FpXQ_powers([-1n], 2, [1n,1n], -2n)).toEqual([[1n],[-1n],[1n]]);
  expect(FpXQ_powers([-1n], 3, [1n,1n], -2n)).toEqual([[1n],[1n],[1n],[1n]]);
  expect(FpX_FpXQ_eval([1n,2n,3n], [1n], [1n,1n], -17n)).toEqual([6n]);
  expect(FpM_mul([[]], [[]], 0n)).toEqual([[]]);
  expect(FpM_mul([[],[0n,1n],[0n,1n]], [[],[0n,1n,-1n]], 0n)).toEqual([[],[0n,0n]]);
});


test('LLM polynomial observation boundaries', async () => {
  const { FpX_deriv, FpX_eval, FpX_center, FpX_div_by_X_x } = await import('@sagemath-ts/parigp-ts');
  expect(FpX_deriv([1n,0n,0n], 0n)).toEqual([]);
  expect(FpX_eval([-17n,1n], 17n, 0n)).toBe(0n);
  expect(FpX_eval([1n,2n,3n], -1n, -17n)).toBe(2n);
  expect(FpX_center([1n,16n], 17n, 8n)).toEqual([1n,-1n]);
  expect(FpX_div_by_X_x([1n,1n], 17n, 0n)).toEqual([1n]);
});


test('LLM exported scalar arithmetic boundaries', async () => {
  const { Fp_inv, Fp_mul, Fp_red, Fp_div, Fp_sqr } = await import('@sagemath-ts/parigp-ts');
  expect(Fp_inv(-1n,17n)).toBe(16n);
  expect(Fp_mul(-1n,1n,17n)).toBe(16n);
  expect(Fp_red(0n,0n)).toBe(0n);
  expect(Fp_div(0n,17n,17n)).toBe(0n);
  expect(Fp_div(1n,2n,(1n<<64n)+13n)).toBe(9223372036854775815n);
  expect(()=>Fp_sqr(0n,0n)).toThrow('impossible inverse in dvmdii: 0.');
});


test('LLM PARI scalar power shortcuts', async () => {
  const { Fp_pow } = await import('@sagemath-ts/parigp-ts');
  expect(Fp_pow(17n,0n,17n)).toBe(0n);
  expect(Fp_pow(0n,0n,0n)).toBe(0n);
  expect(Fp_pow(1n,0n,0n)).toBe(1n);
  expect(Fp_pow(-1n,-1n,17n)).toBe(16n);
});


test('LLM PARI scalar order boundaries', async () => {
  const { Fp_order, znorder, xgcd } = await import('@sagemath-ts/parigp-ts');
  expect(Fp_order(2n,0n,17n)).toBe(8n);
  expect(Fp_order(2n,1n<<64n,17n)).toBe(8n);
  expect(znorder(2n,-17n,0n)).toBe(8n);
  expect(xgcd(0n,0n)).toEqual([0n,0n,0n]);
});


test('LLM PARI rational polynomial conversion', async () => {
  const { QPoly_normalize, QPoly_to_fractions, QPoly_to_FpX } = await import('@sagemath-ts/parigp-ts');
  expect(QPoly_normalize({num:[2n,4n,0n],den:-6n})).toEqual({num:[-1n,-2n],den:3n});
  expect(QPoly_to_fractions({num:[2n,4n,0n],den:-6n})).toEqual([[-1n,3n],[-2n,3n]]);
  expect(QPoly_to_FpX({num:[2n,4n,6n],den:2n},2n)).toEqual([1n,0n,1n]);
});


test('LLM PARI Galois integer helpers', async () => {
  const { ulcm, factoru_small, logint, Forprime } = await import('@sagemath-ts/parigp-ts');
  expect(ulcm(3,2**63)).toBe(2**63);
  expect(factoru_small(0)).toEqual([[0,1]]);
  expect(logint(3n**8192n,3n)).toBe(8192);
  expect(new Forprime(499979).next()).toBe(500009);
});


test('LLM PARI integer-polynomial denominator refinement', async () => {
  const { ZX_neg, ZX_is_squarefree, indexpartial } = await import('@sagemath-ts/parigp-ts');
  expect(ZX_neg([1n, 0n])).toEqual([-1n]);
  expect(ZX_is_squarefree([2n])).toBe(true);
  expect(indexpartial([-2n, 0n, 0n, 0n, 1n])).toBe(8n);
  expect(indexpartial([1n, -2n, 1n])).toBe(0n);
});


test('LLM PARI native Hensel normalization and Bezout lifting', async () => {
  const { ZpX_liftfact, bezout_lift_fact, ZpX_ZpXQ_liftroot } = await import('@sagemath-ts/parigp-ts');
  expect(ZpX_liftfact([4n,1n],[[],[0n,1n]],2n,1)).toEqual([[],[4n,1n]]);
  const f=[-9n,-2n,-1n], factors=[[],[0n,1n],[2n,1n]];
  expect(ZpX_liftfact(f,factors,3n,3)).toEqual([[],[18n,1n],[11n,1n]]);
  expect(bezout_lift_fact(f,factors,3n,3)).toEqual([[],[10n,23n],[18n,4n]]);
  expect(ZpX_ZpXQ_liftroot([-1n,0n,1n],[4n],[1n,0n,1n],3n,1)).toEqual([4n]);
});


test('LLM PARI inverse Vandermonde storage', async () => {
  const { FpV_invVandermonde } = await import('@sagemath-ts/parigp-ts');
  expect(FpV_invVandermonde([0n,0n,1n],1n,5n)).toEqual([[0n,0n,0n],[0n,1n,0n],[0n,4n,1n]]);
});


test('LLM PARI permutation word exponents and prefix comparison', async () => {
  const { perm_powu, perm_cycles, cyc_pow, vecsmall_lexcmp } = await import('@sagemath-ts/parigp-ts');
  const p=[0,2,3,1];
  expect(perm_powu(p,-1)).toEqual([0,2,3,1]);
  expect(cyc_pow(perm_cycles(p),-1)).toEqual([[],[0,1,3,2]]);
  expect(vecsmall_lexcmp([0],[0,1,2,3])).toBe(-1);
});


test('LLM PARI symmetric-polynomial integer combination and scalar zero', async () => {
  const { sympol_eval } = await import('@sagemath-ts/parigp-ts');
  const orbits = [[], [0n, 2n], [0n, 3n]];
  expect(sympol_eval({ v: [0, 4], w: [0, 1] }, orbits, 5n)).toEqual([0n,8n,12n]);
  expect(sympol_eval({ v: [0, -1], w: [0, 1] }, orbits, 5n)).toEqual([0n,-2n,-3n]);
  expect(sympol_eval({ v: [0, 0], w: [0, 1] }, orbits, 5n)).toBe(0n);
});


test('LLM PARI unit-subgroup ordering and unsigned order argument', async () => {
  const { listznstarelts } = await import('@sagemath-ts/parigp-ts');
  expect(listznstarelts(8, 2)).toEqual([[0,1],[0,1,3],[0,1,7],[0,1,5]]);
  expect(listznstarelts(9, -1)).toEqual([[0,1],[0,1,4,7]]);
});


test('LLM PARI subgroup source-subpath bound modes', async () => {
  const { subgrouplist } = await import('@sagemath-ts/parigp-ts/src/subgroup.ts');
  expect(subgrouplist([4n])).toEqual([[[4n]],[[2n]],[[1n]]]);
  expect(subgrouplist([4n], 2n)).toEqual([[[2n]],[[1n]]]);
  expect(subgrouplist([4n], [2n])).toEqual([[[2n]]]);
});


test('LLM PARI conjugate identity and linear shortcuts', async () => {
  const { galoisconj4 } = await import('@sagemath-ts/parigp-ts');
  expect(galoisconj4([-2n, 1n])).toEqual([{ num:[0n,1n],den:1n }]);
  expect(galoisconj4([-2n, 0n, 0n, 1n])).toEqual([{ num:[0n,1n],den:1n }]);
  expect(galoisconj4([-1n, 1n])).toEqual([{ num:[1n],den:1n }]);
});


test('LLM PARI symmetric-search source-subpath helper', async () => {
  const { fixedfieldsurmer } = await import('@sagemath-ts/parigp-ts/src/galconj.ts');
  expect(fixedfieldsurmer(17n, [[], [0n,0n,1n]], [0,1])).toEqual({v:[0,1],w:[0,1]});
});


test('LLM PARI remaining Galois source-subpath helpers', async () => {
  const { eulerphiu, radicalu, findpsi, inittest, galois_test_perm, galoisfindgroups } =
    await import('@sagemath-ts/parigp-ts/src/galconj.ts');
  expect(eulerphiu(0)).toBe(2);
  expect(eulerphiu(-12)).toBe(4);
  expect(radicalu(12)).toBe(6);
  const out={Tmod:[] as bigint[][],psi:[] as number[],p:0};
  expect(findpsi(-4n,2,[1n,0n,1n],{num:[0n,-1n],den:1n},2,out)).toBe(3);
  expect(out.psi).toEqual([0,1]);
  const td=inittest([0n,1n,2n],[[],[0n,1n,0n],[0n,0n,1n]],10n,1n<<96n);
  expect(galois_test_perm(td,[0,2,1])).toBe(true);
  expect(galoisfindgroups([[0,1,3],[0,1,2]],[0,1],2)).toEqual([[0,1,3]]);
});


test('LLM PARI matrix-kernel source-subpath helpers', async () => {
  const { FpM_ker, lcmBig } = await import('@sagemath-ts/parigp-ts/src/galconj.ts');
  const { FpM_ker: kernel } = await import('@sagemath-ts/parigp-ts/src/alglin1.ts');
  const { Flm_ker } = await import('@sagemath-ts/parigp-ts/src/Flv.ts');
  const { F3m_ker } = await import('@sagemath-ts/parigp-ts/src/F3v.ts');
  expect(FpM_ker([[1n,1n],[0n,0n]],2,2,17n)).toEqual([[16n,1n]]);
  expect(lcmBig(-6n,8n)).toBe(24n);
  expect(kernel([[],[0n,1n],[0n,1n]],17n)).toEqual([[],[0n,16n,1n]]);
  expect(Flm_ker([[],[0n,1n],[0n,1n]],17n)).toEqual([[],[0n,16n,1n]]);
  expect(F3m_ker([1n,1n],1)).toEqual([6n]);
});


test('LLM PARI exact Pollard-Brent budget example', async () => {
  const {Z_pollardbrent}=await import('@sagemath-ts/parigp-ts/src/ifactor.ts');
  expect(Z_pollardbrent(1022117n,2**26,0)).toEqual([1009n,1013n]);
});

test('LLM PARI power-stage source-subpath helpers', async () => {
  const { is_357_power, is_pth_power, forprime } = await import('../packages/parigp-ts/src/ifactor.js');
  expect(is_357_power(32n, 7)).toEqual([5, 2n, 2]);
  expect(is_357_power(3n ** 21n, 7)).toEqual([3, 2187n, 5]);
  const primes = forprime(11, 100);
  expect(is_pth_power(103n, { next: () => primes.next().value ?? null }, 1)).toEqual([0, 103n]);
});

test('LLM internal PARI ECM source example', async () => {
  const { ECM_loop } = await import('../packages/parigp-ts/src/_ecm.js');
  expect(ECM_loop(4295229443n, 8, 1, 142)).toBe(65537n);
});

test('LLM PARI signed power and MPQS helper examples', async () => {
  const { is_kth_power } = await import('../packages/parigp-ts/src/ifactor.js');
  const { mpqsInternals } = await import('../packages/parigp-ts/src/mpqs.js');
  expect(is_kth_power(-16n, 3)).toBeNull();
  expect(() => is_kth_power(-8n, 3)).toThrow('sorry, sqrtnr for x < 0 is not yet implemented.');
  expect(is_kth_power(2n, 17886697)).toBeNull();
  expect(mpqsInternals.Fl_sqrt(4, 5)).toBe(2);
});

test('LLM PARI sparse binary kernel example', async () => {
  const { F2Ms_ker } = await import('../packages/parigp-ts/src/F2v.js');
  const { getrand, setrand } = await import('../packages/parigp-ts/src/random.js');
  const saved = getrand();
  try { setrand(1n); expect(F2Ms_ker([[], [], []], 641)).toEqual([4n, 6n, 5n]); }
  finally { setrand(saved); }
});

test('LLM PARI MPQS scalar helper examples', async () => {
  const { mpqsInternals: P } = await import('../packages/parigp-ts/src/mpqs.js');
  expect(P.krouu(1, 2 ** 32)).toBe(1);
  expect(P.kroiu(-100n, 3)).toBe(-1);
  const h = P.newHandle(); h.N = 23n; h.size_of_FB = 6;
  P.mpqs_FB_ctor(h).p.set([0, 0, 2, 3, 5, 7, 11, 13]);
  expect(P.mpqs_factorback(h, [-1048574])).toBe(12n);
  expect(P.mpqs_check_rel(h, { Y: 5n, relp: [1048578] }, 1, 0)).toBeUndefined();
});

test('LLM PARI MPQS inverse example', async () => {
  const { mpqsInternals: P } = await import('../packages/parigp-ts/src/mpqs.js');
  const { PariError } = await import('../packages/parigp-ts/src/errors.js');
  expect(P.Fl_inv(3, 97)).toBe(65);
  expect(P.Fl_inv(0, 1)).toBe(0);
  expect(() => P.Fl_inv(0, 2)).toThrow(PariError);
  expect(() => P.Fl_inv(0, 2)).toThrow('impossible inverse in Fl_inv: Mod(0, 2).');
});

test('LLM PARI MPQS relation hash example', async () => {
  const { relationHash } = await import('../packages/parigp-ts/src/_mpqs_hash.js');
  const { mpqs } = await import('../packages/parigp-ts/src/mpqs.js');
  const { getrand, setrand } = await import('../packages/parigp-ts/src/random.js');
  expect(relationHash({ Y: 1n, relp: [1048578] })).toBe(6171988346495030103n);
  const saved = getrand();
  try { setrand(1n); expect(mpqs(1000509270188293n)).toEqual([[10005089n, 1n], [100000037n, 1n]]); }
  finally { setrand(saved); }
});

test('LLM PARI MPQS warning example', async () => {
  const { mpqsInternals: P } = await import('../packages/parigp-ts/src/mpqs.js');
  const { RelationTable } = await import('../packages/parigp-ts/src/_mpqs_hash.js');
  const h = P.newHandle(); h.N = 143n; h.size_of_FB = 6; h.debug = true;
  P.mpqs_FB_ctor(h).p.set([0, 0, 2, 3, 5, 7, 11, 13]);
  const rows = new RelationTable(20); rows.add({ Y: 10n, relp: [] });
  const warnings: unknown[] = [], saved = console.warn;
  console.warn = message => { warnings.push(message); };
  try {
    expect(P.mpqs_solve_linear_system(h, rows)).toEqual([[13n, 1n], [11n, 1n]]);
    expect(warnings).toEqual(['MPQS: wrong relation found after Gauss']);
  } finally { console.warn = saved; }
});

test('LLM: seeded Pollard lambda with a signed custom hash', async () => {
  const { groups, Mod } = await import('sagemath-ts');
  const { set_random_seed } = await import('sagemath-ts/misc/randstate');
  set_random_seed(0);
  expect(groups.discrete_log_lambda(
    Mod(300n, 17n), Mod(1n, 17n), [0n, 256n], '+',
    undefined, undefined, undefined, v => -(v.value ** 2n + 1n)
  )).toBe(232n);
});

test('LLM: Pollard rho prime-order subgroup', async () => {
  const { groups, Mod } = await import('sagemath-ts');
  const { set_random_seed } = await import('sagemath-ts/misc/randstate');
  set_random_seed(0);
  const base = Mod(4n, 1019n);
  expect(groups.discrete_log_rho(base.pow(250n), base, 509n, '*')).toBe(250n);
});

test('LLM: empty group factor lists request automatic factorization', async () => {
  const { groups, Mod } = await import('sagemath-ts');
  expect(groups.order_from_multiple(Mod(1n, 12n), 12n, [], '+')).toBe(12n);
});

test('LLM: group iterator advances before returning', async () => {
  const { groups, Mod } = await import('sagemath-ts');
  let calls = 0;
  const iterator = groups.multiples(Mod(2n, 11n), 2n, Mod(5n, 11n), false, 'other',
    (x, y) => { calls++; return x.add(y); });
  const first = iterator.next();
  expect(first.done).toBe(false);
  if (first.done) throw new Error('expected a value');
  expect(first.value.value).toBe(5n);
  expect(calls).toBe(1);
});

test('LLM: number-field scalar actions and additive BSGS', async () => {
  const { groups } = await import('sagemath-ts');
  const { QuadraticField, Rational } = await import('sagemath-ts/rings');
  const a = QuadraticField.create(2n, 'a').gen();
  expect(a.mul(new Rational(7n, 3n)).list().map(String)).toEqual(['0', '7/3']);
  expect(groups.bsgs(a, a.mul(31n), [0n, 64n], '+')).toBe(31n);
});

test('LLM: number-field coercion and rational scalar arithmetic', async () => {
  const { QuadraticField, Integer, Rational } = await import('sagemath-ts/rings');
  const K = QuadraticField.create(2n, 'a');
  expect(K.gen().add(new Integer(3n)).list().map(String)).toEqual(['3', '1']);
  expect(K.gen().div(3n).list().map(String)).toEqual(['0', '1/3']);
  expect(K.__call__(1.5).eq(new Rational(3n, 2n))).toBe(true);
});

test('LLM: number-field coefficient indexing', async () => {
  const { QuadraticField } = await import('sagemath-ts/rings');
  const a = QuadraticField.create(2n, 'a').gen().add(2n);
  expect(a.__getitem__(-2).toString()).toBe('2');
  expect(a.__getitem__(-1n).toString()).toBe('1');
});

test('LLM: integer-rational arithmetic and generic inversion', async () => {
  const { Integer, Rational, groups } = await import('sagemath-ts');
  expect(new Integer(2n).add(new Rational(1n, 3n)).toString()).toBe('7/3');
  expect(new Integer(2n).eq(new Rational(2n))).toBe(true);
  expect(groups.bsgs(new Integer(2n), new Integer(8n), [0n, 10n], '*')).toBe(3n);
  const ops = groups.parseGroupOps('*', undefined, undefined, undefined, new Integer(2n));
  expect(ops.inverse(new Integer(2n)).toString()).toBe('1/2');
});

test('LLM: rational equality with real and exact integer operands', async () => {
  const { Rational } = await import('sagemath-ts');
  expect(new Rational(1n, 10n).eq(0.1)).toBe(true);
  const q = new Rational(2n ** 100n + 1n);
  expect(q.eq(2 ** 100)).toBe(true);
  expect(q.eq(2n ** 100n)).toBe(false);
});


test('LLM: RDF scalar arithmetic and binary group multiples', async () => {
  const { RDF, Rational } = await import('sagemath-ts/rings');
  const { groups } = await import('sagemath-ts');
  expect(RDF.__call__(1.5).add(2n).value).toBe(3.5);
  expect(RDF.__call__(1.5).mul(new Rational(1n, 3n)).value).toBe(0.5);
  expect(Object.is(groups.multiple(RDF.__call__(-0), 2n, '*').value, -0)).toBe(true);
});


test('LLM: Pollard lambda single-value bounds retain Sage failure', async () => {
  const { groups, Mod } = await import('sagemath-ts');
  expect(() => groups.discrete_log_lambda(Mod(3n, 11n), Mod(1n, 11n), [3n, 3n], '+'))
    .toThrow('integer modulo by zero');
});


test('LLM: number-field defining equation, reduction and Integer powers', async () => {
  const { NumberField, NumberFieldElement, RationalPolynomial, Rational, Integer } = await import('sagemath-ts/rings');
  const f = RationalPolynomial.fromBigInts([-1n, 2n]);
  const K = new NumberField(f, 'a');
  expect(K.gen().list().map(String)).toEqual(['1/2']);
  expect(K.polynomial()).toBe(f);
  const a = new NumberFieldElement(K, [1n, 2n, 3n].map(c => new Rational(c)));
  expect(a.list().map(String)).toEqual(['11/4']);
  expect(a.pow(new Integer(1n))).toBe(a);
});

test('LLM.md — maximal-order different and general ideal powers', async () => {
  const { QuadraticField, Integer } = await import('sagemath-ts/rings');
  const K = QuadraticField.create(-23n, 'a');
  expect(K.different().norm().toString()).toBe('23');
  expect(K.different()).toBe(K.different());
  const I = K.ideal(2n, K.gen().add(1n));
  expect(I.pow(new Integer(1n))).toBe(I);
  expect(I.mul(I.inverse()).eq(K.ideal(1n))).toBe(true);
});

test('LLM.md — ideal bases and field-valued two-generator tuples', async () => {
  const { QuadraticField, Rational } = await import('sagemath-ts/rings');
  const K = QuadraticField.create(-5n, 'a');
  const I = K.ideal(new Rational(3n, 2n));
  const [a, b] = I.gens_two();
  expect(a.toString()).toBe('3/2');
  expect(b.is_zero()).toBe(true);
  expect(a.parent()).toBe(K);
  expect(I.free_module()).toBe(I.free_module());
  expect(K.ideal(0n).integral_basis().length).toBe(0);
});

test('LLM.md — coprime ideal numerators and denominators', async () => {
  const { QuadraticField, Rational } = await import('sagemath-ts/rings');
  const K = QuadraticField.create(-1n, 'i');
  const I = K.ideal(K.gen().mul(4n).add(3n).div(5n));
  expect(I.numerator().norm().toString()).toBe('5');
  expect(I.denominator().norm().toString()).toBe('5');
  expect(I.numerator().div(I.denominator()).eq(I)).toBe(true);
  const half = K.ideal(new Rational(1n, 2n));
  expect(half.intersection(half).eq(half)).toBe(true);
  expect(K.ideal().is_zero()).toBe(true);
});

test('LLM.md — internal PARI fast LLL stage', async () => {
  const { fplll_fast } = await import('@sagemath-ts/parigp-ts/src/lll.js');
  const [status, B, U] = fplll_fast([[1n, 1n], [0n, 2n]]);
  expect(status).toBe(0);
  expect(B).toEqual([[1n, 1n], [-1n, 1n]]);
  expect(U).toEqual([[1n, 0n], [-1n, 1n]]);
});

test('LLM.md — internal PARI DPE reduction stage', async () => {
  const { fplll_dpe } = await import('@sagemath-ts/parigp-ts/src/lll.js');
  const [status, G, B, U] = fplll_dpe([[1n, 1n], [0n, 2n]]);
  expect(status).toBe(0);
  expect(G).toEqual([[2n, 0n], [0n, 2n]]);
  expect(B).toEqual([[1n, 1n], [-1n, 1n]]);
  expect(U).toEqual([[1n, 0n], [-1n, 1n]]);
});

test('LLM.md — internal real LLL stages use bit precision', async () => {
  const { fplll, fplll_heuristic } = await import('@sagemath-ts/parigp-ts/src/lll.js');
  const B = [[1n, 1n], [0n, 2n]];
  expect(fplll_heuristic(B, 0.99, 0.51, false, true, 128, 256)[1]).toEqual([[1n, 1n], [-1n, 1n]]);
  expect(fplll(B, null, 0.99, 0.51, false, true, true, 128)[4]?.map(r => r.p)).toEqual([128, 128]);
});


test('LLM.md — internal Householder QR column adapters', async () => {
  const { QR_init, R_from_QR, gaussred_from_QR } = await import('../packages/parigp-ts/src/bibli1.js');
  const I = [[1n, 0n], [0n, 1n]];
  expect(QR_init(I)[1]).toEqual([1n, 1n]);
  expect(R_from_QR(I)?.[1]?.[1]).toBe(1n);
  expect(gaussred_from_QR(I)?.[0]).toEqual([1n, 0n]);
});


test('LLM.md — real Cholesky and upper triangular inverse', async () => {
  const { RgM_Cholesky } = await import('../packages/parigp-ts/src/alglin2.js');
  const { RgM_inv_upper } = await import('../packages/parigp-ts/src/alglin1.js');
  const { itor } = await import('../packages/parigp-ts/src/qfb.js');
  const D = [[4n, 0n], [0n, 9n]].map(c => c.map(x => itor(x, 128)));
  expect(RgM_Cholesky(D)?.[0]?.[0]).toEqual(itor(2n, 128));
  expect((RgM_inv_upper(D)[0]?.[0] as ReturnType<typeof itor>).e).toBe(-2);
});


test('LLM.md — integer/real matrix products and rescaling', async () => {
  const { RgV_dotsquare, RgV_dotproduct, gram_matrix, RgM_mul } = await import('../packages/parigp-ts/src/RgV.js');
  const { RgM_rescale_to_int } = await import('../packages/parigp-ts/src/polarit2.js');
  const { itor } = await import('../packages/parigp-ts/src/qfb.js');
  expect(RgV_dotsquare([1n, 2n])).toBe(5n);
  expect(RgV_dotproduct([1n, 2n], [3n, 4n])).toBe(11n);
  expect(gram_matrix([[1n, 2n], [3n, 4n]])).toEqual([[5n, 11n], [11n, 25n]]);
  expect(RgM_mul([[1n, 2n], [3n, 4n]], [[1n, 0n], [0n, 1n]])).toEqual([[1n, 2n], [3n, 4n]]);
  expect(RgM_rescale_to_int([[3n, itor(16n, 64)]])).toEqual([[2n, 8n]]);
});


test('LLM.md — adaptive Gram-Schmidt precision and bounds', async () => {
  const { drop, potential, spread, condition_bound, GS_extraprec, gramschmidt_upper,
    gramschmidt_dynprec, RgM_Cholesky_dynprec } = await import('../packages/parigp-ts/src/lll.js');
  const I = [[1n, 0n], [0n, 1n]];
  expect([drop(I), potential(I), spread(I), condition_bound(I), GS_extraprec(I)]).toEqual([0n, 0n, 0n, 0n, 4n]);
  expect(gramschmidt_upper(I)[0]?.[0]?.p).toBe(64);
  expect((gramschmidt_dynprec(I)[0]?.[0] as { p: number }).p).toBe(64);
  expect((RgM_Cholesky_dynprec(I)[0]?.[0] as { p: number }).p).toBe(64);
});


test('LLM.md — modular pivots and rational matrix solve adapters', async () => {
  const { Flm_pivots, Flm_gauss } = await import('../packages/parigp-ts/src/Flv.js');
  const { ZM_pivots, ZM_rank, ZM_gauss } = await import('../packages/parigp-ts/src/alglin1.js');
  const A = [[], [0n, 1n, 0n], [0n, 0n, 2n]], B = [[], [0n, 3n, 3n]];
  expect(Flm_pivots(A, 101n)).toEqual([[0, 1, 2], 0]);
  expect(Flm_gauss(A, B, 101n)).toEqual([[], [0n, 3n, 52n]]);
  expect(ZM_pivots(A)).toEqual([[0, 1, 2], 0]);
  expect(ZM_rank(A)).toBe(2);
  expect(ZM_gauss(A, B)).toEqual([[[], [0n, 6n, 3n]], 2n]);
});

test('LLM.md — adaptive LLL and FLATTER adapters', async () => {
  const { ZM_lll, lllfp, flat, ZM_flatter, ZM_flatter_rank, ZM_flattergram, LLL_ALL } = await import('../packages/parigp-ts/src/lll.js');
  const { itor, mkqfb, redimagsl2 } = await import('../packages/parigp-ts/src/qfb.js');
  const I = [[1n, 0n], [0n, 1n]], B = [[1n, 0n], [3n, 1n]], U = [[1n, 0n], [-3n, 1n]];
  expect(ZM_lll(B)).toEqual(U);
  expect(ZM_lll(B, .99, LLL_ALL)).toEqual([[], U]);
  expect(lllfp(B.map(c => c.map(x => itor(x, 64))))).toEqual(U);
  expect(flat(I)).toEqual([I, I, 0n, 0n]);
  expect(ZM_flatter(I)).toBeNull();
  expect(ZM_flatter_rank(I, 2)).toBeNull();
  expect(ZM_flattergram(I)).toEqual(I);
  expect(redimagsl2(mkqfb(1n, 6n, 10n, -4n)).U).toEqual([[1n, -3n], [0n, 1n]]);
});


test('LLM.md — ideal list coercion and fractional identity', async () => {
  const { QuadraticField, Rational } = await import('sagemath-ts/rings');
  const K = QuadraticField.create(-1n, 'i');
  const I = K.ideal([new Rational(1n, 2n)]);
  expect(K.fractional_ideal(I)).toBe(I);
  expect(
    I.intersection([new Rational(3n, 2n), 3n])
      .norm()
      .toString()
  ).toBe('9/4');
  expect(I.intersection([]).is_zero()).toBe(true);
  expect(K.ideal().gens()).toHaveLength(1);
});


test('LLM.md — field coefficient vectors and ideal operator coercion', async () => {
  const { NumberField, RationalPolynomial, Rational } = await import('sagemath-ts/rings');
  const K = new NumberField(RationalPolynomial.fromBigInts([3n, 0n, 1n]), 'a');
  expect(K.__call__([1n, new Rational(1n, 2n)]).toString()).toBe('1/2*a + 1');
  expect(K._pari_integral_basis().map(String)).toEqual(['1', '1/2*a - 1/2']);
  const I = K.ideal(2n);
  expect(I.contains(2n)).toBe(true);
  expect(I.contains(new Rational(1n, 2n))).toBe(false);
  expect(I.divides(6n)).toBe(true);
  expect(I.is_coprime(3n)).toBe(true);
  expect(I.add([3n]).norm().toString()).toBe('1');
  expect(I.div(new Rational(1n, 2n)).norm().toString()).toBe('16');
});


test('LLM.md — native Gram image and Hermite kernel columns', async () => {
  const { lllgramint, ZM_hnflll } = await import('@sagemath-ts/parigp-ts');
  expect(
    lllgramint([
      [1n, 1n],
      [1n, 1n],
    ])
  ).toEqual([[1n, 0n]]);
  expect(lllgramint([[0n]])).toEqual([]);
  expect(ZM_hnflll([[], [0n, 0n, 0n], [0n, 0n, 0n]], true, true)).toEqual({
    H: [[]],
    B: [[], [0n, 1n, 0n], [0n, 0n, 1n]],
  });
});


test('LLM.md — direct base ideals and prime factory classes',async()=>{
  const {NumberField,RationalPolynomial,NumberFieldIdeal,NumberFieldFractionalIdeal}=await import('sagemath-ts/rings');
  const K=new NumberField(RationalPolynomial.fromBigInts([-2n,0n,1n]),'a');
  const base=new NumberFieldIdeal(K,[K.gen()]);
  expect(()=>base.inverse()).toThrow("bad operand type for unary ~: 'NumberFieldIdeal'");
  const P=K.decomposition(2n)[0]![0];
  expect(P instanceof NumberFieldFractionalIdeal).toBe(true);
  expect(P.ramification_index()).toBe(2n);
  expect(P.residue_class_degree()).toBe(1n);
  expect(P.inverse().mul(P).norm().toString()).toBe('1');
});


test('LLM.md — ideal valuation receiver, fractions and infinity', async () => {
  const { NumberField, RationalPolynomial } = await import('sagemath-ts/rings');
  const VK = new NumberField(RationalPolynomial.fromBigInts([-2n, 0n, 1n]), 'a');
  const VP = VK.decomposition(2n)[0]![0];
  const VI = VK.ideal(VK.gen().add(1n).div(VK.__call__(2n)));
  expect(VI.valuation(VP)).toBe(-2n);
  expect(VK.ideal(0n).valuation(VP)).toBe('Infinity');
  expect(VK.ideal(8n).valuation(VP)).toBe(6n);
  expect(VP.prime_below()).toBe(2n);
});


test('LLM.md — fractional ideal factorization, order and cache',async()=>{
  const {NumberField,RationalPolynomial,Rational}=await import('sagemath-ts/rings');
  const FK=new NumberField(RationalPolynomial.fromBigInts([1n,0n,1n]),'i');
  const FI=FK.ideal(FK.gen().mul(4n).add(3n).div(FK.__call__(5n)));
  expect(FI.norm().toString()).toBe('1');
  expect(FI.factor().map(([P,e])=>[P.norm().toString(),e])).toEqual([['5',1n],['5',-1n]]);
  expect(FK.factor(new Rational(1n,2n)).map(([,e])=>e)).toEqual([-2n]);
  expect(FI.factor()).toBe(FK.factor(FI));
});


test('LLM.md — exact absolute trace/norm and QQ polynomial adapters',async()=>{
 const {NumberField,RationalPolynomial}=await import('sagemath-ts/rings');
 const TQ=new NumberField(RationalPolynomial.fromBigInts([1n,-1n,2n]),'a');
 expect(TQ.gen().trace().toString()).toBe('1/2');
 expect(TQ.gen().norm().toString()).toBe('1/2');
 expect(TQ.gen().add(1n).trace().toString()).toBe('5/2');
 expect(TQ.gen().add(1n).norm().toString()).toBe('2');
 const {RgXQ_trace,RgXQ_norm}=await import('@sagemath-ts/parigp-ts/src/RgX.js');
 expect(RgXQ_trace([[1n,1n],2n],[[-2n,0n,0n,1n],1n])).toEqual([3n,2n]);
 expect(RgXQ_norm([[1n,1n],2n],[[-2n,0n,0n,1n],1n])).toEqual([3n,8n]);
});


test('LLM.md — trace/norm base fields and result parents',async()=>{
 const {QQ,NumberField,RationalPolynomial}=await import('sagemath-ts/rings');
 const TB=new NumberField(RationalPolynomial.fromBigInts([-2n,0n,0n,1n]),'a'),tb=TB.gen().add(1n);
 expect(tb.trace(QQ).toString()).toBe('3');expect(tb.norm(QQ).toString()).toBe('3');
 expect(tb.trace(TB)).toBe(tb);expect(tb.norm(TB)).toBe(tb);
 const T1=new NumberField(RationalPolynomial.fromBigInts([-3n,1n]),'b');
 expect(tb.trace(T1).parent()).toBe(T1);expect(tb.norm(T1).toString()).toBe('3');
});


test('LLM.md — adaptive LLL squared norm output', async () => {
 const { ZM_lll_norms, LLL_INPLACE } = await import('@sagemath-ts/parigp-ts/src/lll.js');
 const { mpreal_to_frac } = await import('@sagemath-ts/parigp-ts/src/qfb.js');
 const [NB, NN] = ZM_lll_norms([[1n, 0n], [0n, 1n]], 0.99, LLL_INPLACE);
 expect(NB).toEqual([[1n, 0n], [0n, 1n]]);
 expect(NN!.map(mpreal_to_frac)).toEqual([[9223372036854775808n, 9223372036854775808n], [9223372036854775808n, 9223372036854775808n]]);
 expect(ZM_lll_norms([[1n]])[1]).toBeNull();
});


test('LLM.md — permuted Hermite and factor progress kernels',async()=>{
 const {ZM_hnfperm,hnfperm,ZM_hnf_knapsack}=await import('@sagemath-ts/parigp-ts/src/hnf_snf.js');
 const {LLL_check_progress}=await import('@sagemath-ts/parigp-ts/src/QX_factor.js');
 const {itor}=await import('@sagemath-ts/parigp-ts/src/qfb.js');
 const A=[[1n,0n],[2n,1n]],all=[[[1n,0n],[0n,1n]],[[-2n,1n],[1n,0n]],[2,1]];
 expect(hnfperm(A)).toEqual(all);expect(ZM_hnfperm(A,true,true)).toEqual(all);
 expect(ZM_hnf_knapsack(A)).toEqual([[0n,1n],[1n,0n]]);
 expect(LLL_check_progress(itor(3n,64),2,[[1n,0n,0n],[0n,1n,0n],[0n,0n,2n]],true)).toEqual([[1n,0n],[0n,1n]]);
});
test('LLM.md — integer polynomial bounds and sparse evaluation',async()=>{
 const {Mignotte_bound,Beauzamy_bound,factor_bound,root_bound}=await import('@sagemath-ts/parigp-ts/src/QX_factor.js');
 const {vecbinomial}=await import('@sagemath-ts/parigp-ts/src/bibli2.js');
 const {ZX_Z_eval}=await import('@sagemath-ts/parigp-ts/src/ZX.js');
 const {ceil_safe}=await import('@sagemath-ts/parigp-ts/src/gen3.js');
 const P=[-2n,0n,1n];expect(ceil_safe(Mignotte_bound(P))).toBe(4n);expect(ceil_safe(Beauzamy_bound(P))).toBe(4n);
 expect(factor_bound(P)).toBe(4n);expect(root_bound(P)).toBe(2n);
 expect(vecbinomial(5)).toEqual([1n,5n,10n,10n,5n,1n]);expect(ZX_Z_eval([-2n,0n,0n,1n],3n)).toBe(25n);
});
test('LLM.md — native conservative ceiling, half powers and mixed comparison',async()=>{
 const {ceil_safe}=await import('@sagemath-ts/parigp-ts/src/gen3.js');
 const {powruhalf}=await import('@sagemath-ts/parigp-ts/src/trans1.js');
 const {itor}=await import('@sagemath-ts/parigp-ts/src/qfb.js');
 const {cmpir,cmpri}=await import('@sagemath-ts/parigp-ts/src/kernel/none/level1.js');
 expect(ceil_safe(2n)).toBe(2n);expect(ceil_safe(itor(2n,64))).toBe(3n);expect(ceil_safe([3n,2n])).toBe(2n);
 expect(powruhalf(itor(-4n,64),1n)).toEqual({re:0n,im:{s:1,e:1,m:9223372036854775808n,p:64}});
 expect(cmpir((1n<<64n)+1n,itor((1n<<64n)+1n,64))).toBe(0);expect(cmpri(itor(1n,64),2n)).toBe(-1);
});


test('LLM.md — native bounded recombination and exact quotient',async()=>{
 const {ZX_divides_i,ZX_divides,cmbf_precs,cmbf}=await import('@sagemath-ts/parigp-ts/src/QX_factor.js');
 expect(ZX_divides_i([2n,2n],[1n,1n],0n)).toEqual([2n]);
 expect(ZX_divides([-1n,0n,1n],[1n,1n])).toEqual([-1n,1n]);
 expect(cmbf_precs(5n,100n,100n)).toEqual([1,16,3,152587890625n,125n]);
 expect(cmbf([-1n,0n,1n],[[1n,1n],[124n,1n]],3n,5n,3,1,1)).toEqual([[[1n,1n],[-1n,1n]],[[[1n,1n]],[[124n,1n]]],0,true]);
});
test('LLM.md — modular factor products, subset cap and centered ties',async()=>{
 const {cmbf_maxK}=await import('@sagemath-ts/parigp-ts/src/nffactor.js');
 const {FpXV_prod}=await import('@sagemath-ts/parigp-ts/src/FpX.js');
 const {centermodii}=await import('@sagemath-ts/parigp-ts/src/polarit2.js');
 expect(cmbf_maxK(11)).toBe(3);
 expect(FpXV_prod([[1n,1n],[4n,1n]],5n)).toEqual([4n,0n,1n]);
 expect(FpXV_prod([],5n)).toBe(1n);
 expect(centermodii(-5n,10n,5n)).toBe(-5n);
});


test('LLM.md — general Hermite forms and transformation removal',async()=>{
 const {ZM_hnf,ZM_hnfall,ZM_hnfall_i,hnfall}=await import('@sagemath-ts/parigp-ts/src/hnf_snf.js');
 const A=[[1n,2n],[2n,4n],[1n,0n]],H=[[1n,0n],[0n,2n]];
 expect(ZM_hnf(A)).toEqual(H);
 expect(hnfall(A)).toEqual([H,[[-2n,1n,0n],[0n,0n,1n],[1n,0n,-1n]]]);
 expect(ZM_hnfall_i(A,true,2)).toEqual([H,[[0n,0n,1n],[1n,0n,-1n]]]);
 expect(ZM_hnfall(A,false,0)[0]).toEqual([[0n,0n],[1n,0n],[0n,2n]]);
});


test('LLM.md — native integer factors, multiplicities and scaled GCD quotient',async()=>{
 const {ZX_factor,QX_factor,ZX_is_irred,ZX_squff,ZX_gcd_all}=await import('@sagemath-ts/parigp-ts/src/QX_factor.js');
 expect(ZX_factor([-1n,0n,0n,0n,1n])).toEqual([[[-1n,1n],1],[[1n,1n],1],[[1n,0n,1n],1]]);
 expect(QX_factor([[-7n,0n,7n],3n])).toEqual([[[-1n,1n],1],[[1n,1n],1]]);
 expect(ZX_is_irred([-2n,0n,1n])).toBe(true);
 expect(ZX_squff([-1n,3n,-3n,1n])).toEqual([[[-1n,1n]],[3]]);
 expect(ZX_gcd_all([1n,3n,2n],[-1n,-1n,2n])).toEqual([[1n,2n],[2n,2n]]);
});
test('LLM.md — native distinct-factor order and knapsack products',async()=>{
 const {pick_prime,ZX_DDF,ZX_DDF_max,combine_factors,chk_factors,chk_factors_get}=await import('@sagemath-ts/parigp-ts/src/QX_factor.js');
 const P=[-1n,0n,0n,0n,1n];expect(pick_prime(P,0)).toBe(3n);
 expect(ZX_DDF_max(P,0)).toEqual([[1n,1n],[-1n,1n],[1n,0n,1n]]);
 expect(ZX_DDF(P)).toEqual([[1n,0n,1n],[1n,1n],[-1n,1n]]);
 expect(combine_factors([-1n,0n,1n],[[1n,1n],[2n,1n]],3n,1)).toEqual([[1n,1n],[-1n,1n]]);
 expect(chk_factors([-1n,0n,1n],[[1n,0n],[0n,1n]],3n,[[1n,1n],[124n,1n]],125n)).toEqual([[-1n,1n],[1n,1n]]);
 expect(chk_factors_get(null,[[1n,1n],[4n,1n]],[1n,0n],null,5n)).toEqual([1n,1n]);
});
test('LLM.md — exact and modular Newton sums with cached prefixes',async()=>{
 const {polsym,polsym_gen}=await import('@sagemath-ts/parigp-ts/src/polarit2.js');
 expect(polsym([1n,0n,2n],4)).toEqual([[2n,1n],[0n,1n],[-1n,1n],[0n,1n],[1n,2n]]);
 expect(polsym_gen([-2n,0n,1n],[2n,0n],4,null,5n)).toEqual([2n,0n,-1n,0n,-2n]);
});
test('LLM.md — degree-200 number-field constructor uses native factorization',async()=>{
 const {NumberField,RationalPolynomial}=await import('sagemath-ts/rings');
 const f=RationalPolynomial.fromBigInts([-2n,...Array<bigint>(199).fill(0n),1n]);
 expect(f.isIrreducible()).toBe(true);expect(new NumberField(f,'a').degree()).toBe(200);
});


test('LLM.md — NTL cached traces and chopped combinations', async () => {
 const {createRequire}=await import('node:module');
 const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {ComputeTrace,ChopTraces,DenseChopTraces}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
 expect(ComputeTrace([0n,0n,91n],[-2n,0n,1n],2,125n)).toEqual([0n,4n,91n]);
 expect(ChopTraces([0n,0n,99n],[4n,8n],2,[5n,5n],5n,125n,2n)).toEqual([2n,1n,99n]);
 expect(DenseChopTraces([0n,0n,99n],[4n,8n],2,2,5n,5n,125n,2n,[[1n,2n],[-1n,1n]])).toEqual([-1n,0n,99n]);
});
test('LLM.md — NTL trace precisions and factor lattice', async () => {
 const {createRequire}=await import('node:module');
 const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {Compute_pb,Compute_pdelta,Compute_pb_eff,d1_val,BuildReductionMatrix}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
 expect(Compute_pb([],[],3,1,2n,4)).toEqual([[3],[27n]]);
 expect(Compute_pdelta(1,3n,3,4)).toEqual([3,27n]);
 expect(Compute_pb_eff(3,2,1n,4,1)).toEqual([4,81n]);
 expect(d1_val(8,16,8)).toBe(5);
 expect(BuildReductionMatrix(2,2,9n,[[1n,2n],[3n,4n]],[[1n,0n],[1n,-1n]])).toEqual([[[2n,0n,1n,2n],[2n,-2n,-2n,-2n],[0n,0n,9n,0n],[0n,0n,0n,9n]],2]);
});


test('LLM.md — exact NTL integer products', async () => {
 const {createRequire}=await import('node:module');
 const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {mul,sqr}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZX1.js'));
 expect(mul([-2n,0n,3n],[5n,-7n,1n])).toEqual([-10n,14n,13n,-21n,3n]);
 expect(sqr([-2n,0n,3n])).toEqual([4n,0n,-12n,0n,9n]);
});
test('LLM.md — NTL multifactor Hensel lifting', async () => {
 const {createRequire}=await import('node:module');
 const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {MultiLift}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
 expect(MultiLift([[1n,1n],[2n,1n]],[7n,3n,1n],1,5n)).toEqual([[1n,1n],[2n,1n]]);
 expect(MultiLift([[1n,1n],[2n,1n]],[7n,3n,1n],4,5n)).toEqual([[281n,1n],[347n,1n]]);
});


test('LLM.md — NTL exact lattice bases and solutions',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {LLL_plus,image,LatticeSolve}=await import(require.resolve('@sagemath-ts/ntl-ts/src/LLL.js'));
 expect(LLL_plus([[1n,1n,1n],[-1n,0n,2n],[3n,5n,6n]],3,4,true)).toEqual([3,[1n,1n,2n,9n],[[0n,1n,0n],[1n,0n,1n],[-1n,0n,2n]],[[-4n,-1n,1n],[5n,1n,-1n],[0n,1n,0n]]]);
 expect(image([[2n,0n],[0n,3n],[1n,1n]])[2]).toEqual([[0n,0n],[1n,0n],[0n,1n]]);
 expect(LatticeSolve([[2n,0n],[0n,3n],[1n,1n]],[1n,2n])).toEqual([1,[1n,1n,-1n]]);
 expect(LatticeSolve([[2n,0n],[0n,2n]],[1n,0n],0,[19n])).toEqual([0,[19n]]);
});
test('LLM.md — NTL extended GCD zero convention',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {_ntl_gexteucl}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lip.js'));
 expect(_ntl_gexteucl(0n,0n)).toEqual([1n,0n,0n]);expect(_ntl_gexteucl(3n,2n)).toEqual([1n,-1n,1n]);
});


test('LLM.md — NTL root bounds and factor lattice rows',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {PolyEval,RootBound,CutAway}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
 expect(PolyEval([-2n,0n,1n],3n)).toBe(7n);expect(RootBound([-2n,0n,1n])).toBe(2n);expect(RootBound([-2n,3n])).toBe(3n);expect(RootBound([-7n])).toBe(7n);
 expect(CutAway([1n,1n,1n,100n],[[1n,0n,0n],[0n,1n,0n],[0n,0n,10n]],1,2,1)).toEqual([[1n,0n],[0n,1n]]);
});
test('LLM.md — NTL additional factor lifting',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {AdditionalLifting}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
 expect(AdditionalLifting(5n,1,[[1n,1n],[2n,1n]],5n,4,[7n,3n,1n],false)).toEqual([625n,4,[[281n,1n],[347n,1n]]]);
});


test('LLM.md — NTL factor degree patterns',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {CalcPossibleDegrees,NumFactors,unpack,SubPattern}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
 expect(CalcPossibleDegrees([0,2,1])).toBe(31n);expect(NumFactors([0,2,1])).toBe(3);expect(unpack(-13n,4)).toEqual([1,0,1,1,0]);expect(SubPattern([0,2,1],[0,1,0])).toEqual([0,1,1]);
});
test('LLM.md — NTL factor products and divisibility cache',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {mul,InvMul,ConstTermTest,BalCopy}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
 const W=[[1n,1n],[2n,1n],[3n,1n]];expect(mul(W,101n)).toEqual([6n,11n,6n,1n]);expect(InvMul(W,[1],101n)).toEqual([3n,4n,1n]);
 expect(ConstTermTest(W,[1,2],12n,1n,[7n,8n,9n],0,101n)).toEqual([1,[2n,6n,9n],1]);expect(BalCopy([8n,9n,16n],16n)).toEqual([8n,-7n]);
});


test('LLM.md — NTL word polynomial powers and inverse gcd',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {zz_pXModulus,PowerXMod,PowerXPlusAMod,PowerMod,InvModStatus}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));
 const F=new zz_pXModulus([1n,0n,1n],5n);expect(PowerXMod(3n,F)).toEqual([0n,4n]);expect(PowerXMod(-1n,F)).toEqual([0n,4n]);expect(PowerXPlusAMod(1n,3n,F)).toEqual([3n,2n]);expect(PowerMod([1n,1n],3n,F)).toEqual([3n,2n]);expect(InvModStatus([2n,1n],F.f,5n)).toEqual([1,[2n,1n]]);
});
test('LLM.md — NTL word polynomial reduction and rebuild',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {zz_pXModulus,rem,MulMod,SqrMod,build}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));
 const F=new zz_pXModulus([1n,0n,1n],5n);expect(rem([1n,2n,3n,4n],F)).toEqual([3n,3n]);expect(MulMod([1n,1n],[2n,1n],F)).toEqual([1n,3n]);expect(SqrMod([1n,1n],F)).toEqual([0n,2n]);build(F,[-1n,1n]);expect(rem([1n,2n,3n],F)).toEqual([1n]);
});


test('LLM.md — NTL word matrix products',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {mul}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_lzz_p.js'));const A=[[1n,2n],[3n,4n]];
 expect(mul(A,[[5n,6n],[7n,8n]],101n)).toEqual([[19n,22n],[43n,50n]]);expect(mul(A,A,5n)).toEqual([[2n,0n],[0n,2n]]);
});
test('LLM.md — NTL word matrix zero inner dimension',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {mul}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_lzz_p.js'));
 expect(mul([[],[],[]],[],17n,4,0)).toEqual([[0n,0n,0n,0n],[0n,0n,0n,0n],[0n,0n,0n,0n]]);
});


test('LLM.md — NTL block composition and shared arguments',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {zz_pXModulus}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));const {zz_pXNewArgument,build,CompMod,Comp2Mod,Comp3Mod}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX1.js'));
 const F=new zz_pXModulus([1n,0n,1n],5n),H=new zz_pXNewArgument();build(H,[1n,1n],F,2);expect(H.mat).toEqual([[1n,0n],[1n,1n]]);expect(H.poly).toEqual([0n,2n]);expect(CompMod([1n,2n,3n],H,F)).toEqual([3n,3n]);expect(Comp2Mod([1n,2n,3n],[0n,1n],[1n,1n],F)).toEqual([[3n,3n],[1n,1n]]);expect(Comp3Mod([1n,2n,3n],[0n,1n],[1n],[1n,1n],F)).toEqual([[3n,3n],[1n,1n],[1n]]);
});
test('LLM.md — NTL composition cache reduction and word GCD',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {zz_pXModulus}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));const {zz_pXNewArgument,build,reduce,CompMod,GCD}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX1.js'));
 const F=new zz_pXModulus([1n,0n,1n],5n),H=new zz_pXNewArgument();build(H,[1n,1n],F,2);const G=new zz_pXModulus([-1n,1n],5n);reduce(H,G);expect(H.mat).toEqual([[1n],[2n]]);expect(H.poly).toEqual([4n]);expect(CompMod([1n,2n,3n],H,G)).toEqual([2n]);expect(GCD([-2n,0n,2n],[-3n,3n],5n)).toEqual([4n,1n]);
});
test('LLM.md — NTL composition trace and Frobenius iteration',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {zz_pXModulus}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));const {TraceMap,PowerCompose}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js'));
 const F=new zz_pXModulus([1n,0n,1n],3n),b=[0n,2n];expect(TraceMap([1n,1n],2,F,b)).toEqual([2n]);expect(PowerCompose(b,2,F)).toEqual([0n,1n]);
});


test('LLM.md — NTL distinct-degree groups and Frobenius',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {NewDDF,SFCanZass1}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js'));
 const f=[-2n,2n,-1n,1n],[groups,h]=SFCanZass1(f,5n);expect(groups).toEqual([[[4n,1n],1],[[2n,0n,1n],2]]);expect(h).toEqual([3n,4n,4n]);expect(NewDDF(f,h,5n)).toEqual(groups);
});
test('LLM.md — NTL distinct-degree constant and linear cases',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {NewDDF}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js'));
 expect(NewDDF([1n],[],5n)).toEqual([]);expect(NewDDF([-1n,1n],[1n],5n)).toEqual([[[4n,1n],1]]);
});


test('LLM.md — NTL sequence minimum polynomial and unused suffix',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {MinPolySeq}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX1.js'));const a=[0n,1n,1n,2n,3n,5n];expect(MinPolySeq(a,3,101n)).toEqual([100n,100n,1n]);expect(MinPolySeq([...a,999n,-1n],3,101n)).toEqual([100n,100n,1n]);
});
test('LLM.md — NTL zero and impulse sequence polynomials',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {MinPolySeq}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX1.js'));expect(MinPolySeq([0n,0n,0n,0n],2,101n)).toEqual([1n]);expect(MinPolySeq([1n,0n,0n,0n],2,101n)).toEqual([0n,1n]);expect(MinPolySeq([],0,101n)).toEqual([1n]);
});


test('LLM.md — NTL prime iterator and reset',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {PrimeSeq}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const primes=new PrimeSeq();expect(Array.from({length:6},()=>primes.next())).toEqual([2n,3n,5n,7n,11n,13n]);primes.reset(100n);expect(primes.next()).toBe(101n);expect(primes.next()).toBe(103n);
});
test('LLM.md — NTL prime iterator exhaustion and restart',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {PrimeSeq}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const primes=new PrimeSeq();primes.reset(1073676289n);expect(primes.next()).toBe(0n);expect(primes.next()).toBe(0n);primes.reset(-5n);expect(primes.next()).toBe(2n);
});


test('LLM.md — NTL cached multiplier and transposed product',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {zz_pXModulus,zz_pXMultiplier,MulMod}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));const {UpdateMap}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX1.js'));
 const F=new zz_pXModulus([1n,0n,1n],5n),B=new zz_pXMultiplier([1n,1n],F);expect(B.val()).toEqual([1n,1n]);expect(B.UseFFT).toBe(0);expect(MulMod([1n,2n],B,F)).toEqual([4n,3n]);expect(UpdateMap([1n,2n],B,F)).toEqual([3n,1n]);
});
test('LLM.md — NTL powers projected through a shared argument',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {zz_pXModulus}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));const {zz_pXNewArgument,build,ProjectPowers}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX1.js'));
 const F=new zz_pXModulus([1n,0n,1n],5n),a=[1n,2n],H=new zz_pXNewArgument();expect(ProjectPowers(a,5,[0n,1n],F)).toEqual([1n,2n,4n,3n,1n]);build(H,[0n,1n],F,2);expect(ProjectPowers(a,5,H,F)).toEqual([1n,2n,4n,3n,1n]);expect(ProjectPowers(a,0,[0n,1n],F)).toEqual([]);
});


test('LLM.md — NTL SHA and generic ChaCha bytes',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {sha256,RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 const hex=(a:Uint8Array)=>Array.from(a,x=>x.toString(16).padStart(2,'0')).join('');
 expect(hex(sha256(new Uint8Array()))).toBe('e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855');
 expect(hex(new RandomStream(new Uint8Array(32)).get(16))).toBe('76b8e0ada0f13d90405d6ae55386bd28');
});
test('LLM.md — NTL stream copies and nonce reset',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 const stream=new RandomStream(new Uint8Array(32));stream.get(17);const copy=new RandomStream(stream);
 expect(Array.from(stream.get(5))).toEqual([210,25,184,160,141]);expect(Array.from(copy.get(5))).toEqual([210,25,184,160,141]);stream.set_nonce(0n);expect(Array.from(stream.get(4))).toEqual([118,184,224,173]);
});


test('LLM.md — NTL bounded-sampling overload byte consumption',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {RandomStream,RandomBnd}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 const a=new RandomStream(new Uint8Array(32)),b=new RandomStream(new Uint8Array(32));
 expect(RandomBnd(256n,a)).toBe(118n);expect(RandomBnd(256n,b,{word:true})).toBe(118n);expect(Array.from(a.get(4))).toEqual([224,173,160,241]);expect(Array.from(b.get(4))).toEqual([184,224,173,160]);
});
test('LLM.md — NTL exact-length overload byte consumption',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {RandomStream,RandomLen_long,RandomLen_ZZ}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 const a=new RandomStream(new Uint8Array(32)),b=new RandomStream(new Uint8Array(32));
 expect(RandomLen_long(9,a)).toBe(374n);expect(RandomLen_ZZ(9,b)).toBe(374n);expect(Array.from(a.get(4))).toEqual([184,224,173,160]);expect(Array.from(b.get(4))).toEqual([224,173,160,241]);
});


test('LLM.md — NTL probable and certified element minimum polynomials',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const {zz_pXModulus}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));const {ProbMinPolyMod,MinPolyMod}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX1.js'));
 const F=new zz_pXModulus([1n,1n,1n],2n),a=new RandomStream(new Uint8Array(32)),b=new RandomStream(new Uint8Array(32));
 expect(ProbMinPolyMod([0n,1n],F,a)).toEqual([1n]);expect(MinPolyMod([0n,1n],F,b)).toEqual([1n,1n,1n]);expect(Array.from(a.get(4))).toEqual([224,173,160,241]);expect(Array.from(b.get(4))).toEqual([160,241,61,144]);
});
test('LLM.md — NTL deterministic quotient projection',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {zz_pXModulus}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));const {DoMinPolyMod,IrredPolyMod}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX1.js'));
 const F=new zz_pXModulus([1n,1n,1n],2n);expect(DoMinPolyMod([0n,1n],F,2,[1n])).toEqual([1n,1n,1n]);expect(IrredPolyMod([0n,1n],F)).toEqual([1n,1n,1n]);
});


test('LLM.md — NTL monic product from roots',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {BuildFromRoots}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX1.js'));expect(BuildFromRoots([1n,2n,0n],5n)).toEqual([0n,2n,2n,1n]);
});
test('LLM.md — NTL normalized random polynomial and consumed bytes',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const {random}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));
 const stream=new RandomStream(new Uint8Array(32));expect(random(3,5n,stream)).toEqual([]);expect(Array.from(stream.get(4))).toEqual([241,61,144,64]);
});


test('LLM.md — NTL ordered roots and one-root search',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const {FindRoots,FindRoot}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js'));
 expect(FindRoots([0n,1n,1n],2n,new RandomStream(new Uint8Array(32)))).toEqual([1n,0n]);expect(FindRoot([0n,1n,1n],2n,new RandomStream(new Uint8Array(32)))).toBe(1n);
});
test('LLM.md — NTL equal-degree and square-free factor recovery',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const {EDF,SFCanZass}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js'));
 const f=[2n,1n,0n,1n,1n];expect(EDF(f,[0n,0n,0n,1n],2,3n,new RandomStream(new Uint8Array(32)))).toEqual([[1n,0n,1n],[2n,1n,1n]]);expect(SFCanZass(f,3n,new RandomStream(new Uint8Array(32)))).toEqual([[1n,0n,1n],[2n,1n,1n]]);
});


test('LLM.md — NTL requested context changes construction crossover',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {zz_pXModulus}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));
 const F=new zz_pXModulus([1n,...Array<bigint>(63).fill(0n),3n],9n,{maxroot:60});expect(F.n).toBe(64);expect(F.f[64]).toBe(3n);
});
test('LLM.md — NTL factorization forwards requested context',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const {SFCanZass}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pXFactoring.js'));
 const stream=new RandomStream(new Uint8Array(32));expect(SFCanZass([2n,1n,0n,1n,1n],3n,stream,{maxroot:60})).toEqual([[1n,0n,1n],[2n,1n,1n]]);expect(Array.from(stream.get(4))).toEqual([93,106,229,83]);
});


test('LLM.md — NTL small-prime factors and explicit selected context',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {LocalInfoT,SmallPrimeFactorization}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 const info=new LocalInfoT(),stream=new RandomStream(new Uint8Array(32));expect(SmallPrimeFactorization(info,[-1n,0n,1n],stream)).toEqual([[1n,1n],[2n,1n]]);expect(info.context).toEqual({p:3n,maxroot:2});expect(info.NumPrimes).toBe(6);expect(Array.from(stream.get(4))).toEqual([184,224,173,160]);
});
test('LLM.md — NTL local vector storage survives shrink and growth',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {LocalInfoT,SmallPrimeFactorization}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 const info=new LocalInfoT();info.p=[101n,103n,107n];info.p=[2n];info.pattern=[[9],[8],[7]];info.pattern=[[]];info.s.reset(3n);
 expect(SmallPrimeFactorization(info,[-2n,0n,1n],new RandomStream(new Uint8Array(32)),{InitNumPrimes:3})).toBeNull();expect(info.p).toEqual([3n,103n,107n]);expect(info.pattern).toEqual([[0,0,1],[8],[7]]);
});


test('LLM NTL retained factor-information cache update',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {LocalInfoT,SmallPrimeFactorization,UpdateLocalInfo}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
 const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 const info=new LocalInfoT();SmallPrimeFactorization(info,[-1n,0n,1n],new RandomStream(new Uint8Array(32)));
 expect(UpdateLocalInfo(info,[[1n,1n]],[[-1n,1n]],[1n,1n],1,125n)).toEqual([2n]);expect(info.NumFactors).toBe(1);expect(info.PossibleDegrees).toBe(3n);expect(info.context).toEqual({p:3n,maxroot:2});
});
test('LLM NTL retained factor-information additional prime',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {LocalInfoT,UpdateLocalInfo}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
 const info=new LocalInfoT();info.n=2;info.PossibleDegrees=7n;
 expect(UpdateLocalInfo(info,[],[],[-1n,0n,1n],1,125n,{van_Hoeij:0,MaxNumPrimes:3})).toBeNull();expect(info.p).toEqual([3n]);expect(info.NumPrimes).toBe(1);expect(info.context).toBeNull();
});


test('LLM NTL basic integer-factor cardinality search',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {LocalInfoT,CardinalitySearch}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
 const info=new LocalInfoT();info.n=4;info.PossibleDegrees=31n;const W=[[-1n,1n],[1n,1n],[2n,1n],[3n,1n]];
 expect(CardinalitySearch([],[-6n,-5n,5n,5n,1n],W,info,1,100,2n**127n-1n)).toEqual([[[-1n,1n],[1n,1n],[2n,1n]],[3n,1n],[[3n,1n]]]);
});
test('LLM NTL optimized integer-factor cardinality search',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {LocalInfoT,CardinalitySearch1}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
 const info=new LocalInfoT();info.n=4;info.PossibleDegrees=31n;const W=[[-1n,1n],[1n,1n],[2n,1n],[3n,1n]];
 expect(CardinalitySearch1([],[-6n,-5n,5n,5n,1n],W,info,2,100,2n**127n-1n)).toEqual([[[-1n,0n,1n]],[6n,5n,1n],[[2n,1n],[3n,1n]]]);
});
test('LLM NTL integer-factor recovery order',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {LocalInfoT,FindTrueFactors}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
 const info=new LocalInfoT();info.n=4;info.PossibleDegrees=31n;const W=[[-1n,1n],[1n,1n],[2n,1n],[3n,1n]];
 expect(FindTrueFactors([-6n,-5n,5n,5n,1n],W,2n**127n-1n,info,100)).toEqual([[-1n,1n],[1n,1n],[2n,1n],[3n,1n]]);
});


test('LLM NTL probable-prime bounds and consumed bytes',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {ComputePrimeBound,ErrBoundTest,ProbPrime,RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 expect(ComputePrimeBound(2305)).toBe(23050n);expect(ErrBoundTest(128,8,80)).toBe(0);
 const stream=new RandomStream(new Uint8Array(32));expect(ProbPrime(17n,stream,{NumTrials:2})).toBe(1);expect(Array.from(stream.get(4))).toEqual([64,93,106,229]);
});
test('LLM NTL sequential prime generators',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {RandomPrime,OldRandomPrime,RandomPrime_long,RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const fresh=()=>new RandomStream(new Uint8Array(32));
 expect(RandomPrime(32,fresh())).toBe(2419978657n);expect(OldRandomPrime(32,fresh())).toBe(2419978657n);expect(RandomPrime_long(32,fresh())).toBe(2915387179n);
});
test('LLM NTL error-controlled and seeded-block prime generators',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {GenPrime,GenPrime_long,RandomPrime,RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const fresh=()=>new RandomStream(new Uint8Array(32));
 expect(GenPrime(32,fresh())).toBe(2419978657n);expect(GenPrime_long(60,fresh())).toBe(1007410285329247537n);
 expect(RandomPrime(256,fresh(),{NumTrials:2})).toBe(99719095675712960086022384902837345215225263017301052027152607768417317746133n);
});


test('LLM NTL word elimination and left kernel',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {gauss,image,kernel}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_lzz_p.js'));const A=[[1n,2n],[2n,4n],[0n,1n]];
 expect(gauss(A,5n)).toEqual([2,[[1n,2n],[0n,1n],[0n,0n]]]);expect(image(A,5n)).toEqual([[1n,2n],[0n,1n]]);expect(kernel(A,5n)).toEqual([[3n,1n,0n]]);
});
test('LLM NTL word inverse and determinant',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {relaxed_inv,inv,relaxed_determinant,determinant}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_lzz_p.js'));const A=[[1n,2n],[3n,4n]];
 expect(relaxed_inv(A,5n)).toEqual([3n,[[3n,1n],[4n,2n]]]);expect(inv(A,5n)).toEqual([[3n,1n],[4n,2n]]);expect(relaxed_determinant(A,5n)).toBe(3n);expect(determinant(A,5n)).toBe(3n);
});
test('LLM NTL word solve orientations',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {relaxed_solve,solve}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_lzz_p.js'));const A=[[1n,2n],[3n,4n]],b=[1n,0n];
 expect(relaxed_solve(A,b,5n)).toEqual([3n,[3n,4n]]);expect(solve(A,b,5n)).toEqual([3n,[3n,4n]]);expect(relaxed_solve(A,b,5n,{left:true})).toEqual([3n,[3n,1n]]);
});
test('LLM NTL composite pivots and retained singular outputs',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {relaxed_inv,inv,relaxed_solve}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_lzz_p.js'));
 expect(relaxed_inv([[3n,1n],[1n,1n]],9n)).toEqual([2n,[[5n,4n],[4n,6n]]]);expect(()=>inv([[3n,1n],[1n,1n]],9n)).toThrow('InvMod: inverse undefined');
 expect(relaxed_inv([[3n]],9n,{previous:[[7n]]})).toEqual([0n,[[7n]]]);expect(relaxed_solve([[3n]],[1n],9n,{previous:[11n,12n]})).toEqual([0n,[2n,3n]]);
});


test('LLM NTL FFT-prime recognition and root cap',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {IsFFTPrime,CalcMaxRoot}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const stream=new RandomStream(new Uint8Array(32));
 expect(IsFFTPrime(97n,stream)).toEqual([1,46n]);expect(IsFFTPrime(21n,stream,123n)).toEqual([0,123n]);expect(CalcMaxRoot(97n)).toBe(5);
});
test('LLM NTL FFT root and inverse-power tables',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {InitFFTPrimeInfo}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const info=InitFFTPrimeInfo(17n,3n);
 expect(info.RootTable).toEqual([[1n,16n,13n,9n,3n],[1n,16n,4n,2n,6n]]);expect(info.TwoInvTable).toEqual([1n,9n,13n,15n,16n]);
});
test('LLM NTL FFT cached initialization and following bytes',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {FFTPrimeContext,UseFFTPrime,GetFFTPrime,GetFFTPrimeRecip}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const context=new FFTPrimeContext(),stream=new RandomStream(new Uint8Array(32));
 expect(context.length()).toBe(0);UseFFTPrime(0,context,stream);expect(context.length()).toBe(1);expect(GetFFTPrime(0,context)).toBe(882705526964617217n);expect(GetFFTPrimeRecip(0,context)).toBe(1.1328806373726086e-18);expect(context.get(0).RootTable[0][0]).toBe(1n);UseFFTPrime(0,context,stream);expect(Array.from(stream.get(4))).toEqual([152,186,151,124]);
});
test('LLM NTL FFT repeated-index candidate helper',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {FFTPrimeContext,NextFFTPrime,UseFFTPrime,GetFFTPrime}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const context=new FFTPrimeContext(),stream=new RandomStream(new Uint8Array(32));
 expect(NextFFTPrime(0,context,stream)[0]).toBe(882705526964617217n);expect(context.length()).toBe(0);expect(NextFFTPrime(0,context,stream)[0]).toBe(882705526964617217n);UseFFTPrime(0,context,stream);expect(GetFFTPrime(0,context)).toBe(882705526964617217n);
});

test('LLM NTL scalar CRT and signed interval',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {CRT,CRTInRange}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 expect(CRTInRange(-3n,6n)).toBe(0);expect(CRTInRange(3n,6n)).toBe(1);expect(CRT(1n,3n,0n,2n)).toEqual([1,-2n,6n]);expect(CRT(1n,3n,1n,2n)).toEqual([0,1n,6n]);
});
test('LLM NTL matrix CRT',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {CRT}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_ZZ.js'));
 expect(CRT([[1n,-1n]],3n,[[0n,1n]],2n)).toEqual([1,[[-2n,-1n]],6n]);
});
test('LLM NTL integer determinant bound',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {DetBound}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_ZZ.js'));
 expect(DetBound([[3n,4n],[0n,1n]])).toBe(3);
});
test('LLM NTL arbitrary-modulus determinant',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {determinant}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_ZZ_p.js'));const p=(1n<<127n)-1n;
 expect(determinant([[p+2n,-1n],[1n,p+3n]],p)).toBe(7n);
});

test('LLM NTL integer determinant and adjugate',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {determinant,inv}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_ZZ.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const context=new FFTPrimeContext(),stream=new RandomStream(new Uint8Array(32)),A=[[1n,2n],[3n,4n]];
 expect(determinant(A,context,stream)).toBe(-2n);expect(inv(A,context,stream,{status:true})).toEqual([-2n,[[4n,-2n],[-3n,1n]]]);
});
test('LLM NTL strict integral inverse with negative determinant',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {inv}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_ZZ.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const context=new FFTPrimeContext(),stream=new RandomStream(new Uint8Array(32));
 expect(inv([[1n,2n],[3n,5n]],context,stream)).toEqual([[-5n,2n],[3n,-1n]]);
});
test('LLM NTL retained singular integer inverse',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {inv}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_ZZ.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const context=new FFTPrimeContext(),stream=new RandomStream(new Uint8Array(32));
 expect(inv([[1n,2n],[2n,4n]],context,stream,{status:true,previous:[[7n]]})).toEqual([0n,[[7n]]]);
});
test('LLM NTL deterministic integer reconstruction path',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {determinant}=await import(require.resolve('@sagemath-ts/ntl-ts/src/mat_ZZ.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const context=new FFTPrimeContext(),stream=new RandomStream(new Uint8Array(32));
 expect(determinant([[1n,1n<<1024n],[0n,1n]],context,stream,{deterministic:true})).toBe(1n);expect(context.length()).toBe(18);
});

test('LLM NTL certified square elimination keeps signed scale',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {gauss}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const context=new FFTPrimeContext(),stream=new RandomStream(new Uint8Array(32));
 expect(gauss([[1n,2n],[3n,4n]],context,stream)).toEqual([-2n,[[-2n,0n],[0n,-2n]]]);
});
test('LLM NTL certified rectangular elimination',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {gauss}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const context=new FFTPrimeContext(),stream=new RandomStream(new Uint8Array(32));
 expect(gauss([[1n,2n,3n],[0n,1n,4n]],context,stream)).toEqual([1n,[[1n,0n,-5n],[0n,1n,4n]]]);
});


test('LLM NTL modular product and alias-square stream thresholds',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {mul}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ_pX.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 const context=new FFTPrimeContext(),stream=new RandomStream(new Uint8Array(32)),p=(1n<<127n)-1n,a=Array<bigint>(80).fill(1n);
 mul(a,a.slice(),p,{context,stream});expect(context.length()).toBe(0);const square=mul(a,a,p,{context,stream});expect([square.length,square[0],square[79],square[158]]).toEqual([159,1n,80n,1n]);expect(context.length()).toBe(5);
});
test('LLM NTL selected polynomial product shares FFT initialization',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));
 const {mul}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 const context=new FFTPrimeContext(),stream=new RandomStream(new Uint8Array(32)),W=Array.from({length:4},()=>Array<bigint>(101).fill(1n));const product=mul(W,(1n<<127n)-1n,[3,1,3,0],{context,stream});expect([product.length,product[0],product[400]]).toEqual([401,1n,1n]);expect(context.length()).toBe(5);
});
test('LLM NTL Schoenhage Strassen padding ratio',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {SSRatio}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZX1.js'));
 expect(SSRatio(1,1,1,1)).toBe(1.4);expect(SSRatio(0,127,199,127)).toBe(0);
});


test('LLM NTL integer polynomial zero-numerator dispatch',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {divide,PlainDivide,HomDivide}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZX1.js'));const options={previous:[17n,-2n]};
 expect(divide([],[1n,1n],options)).toEqual([0,[17n,-2n]]);expect(PlainDivide([],[1n,1n],options)).toEqual([0,[17n,-2n]]);expect(HomDivide([],[1n,1n],options)).toEqual([1,[]]);expect(divide([],[],options)).toEqual([1,[]]);
});
test('LLM NTL integer polynomial scalar division',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {divide}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZX1.js'));
 expect(divide([2n,-4n],-2n)).toEqual([1,[-1n,2n]]);expect(divide([1n,2n],2n,{previous:[7n]})).toEqual([0,[7n]]);
});
test('LLM NTL homomorphic polynomial division preserves shared prime state',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {HomDivide}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZX1.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const context=new FFTPrimeContext(),stream=new RandomStream(new Uint8Array(32));
 expect(HomDivide([1n,2n,1n],[1n,1n],{state:{context,stream}})).toEqual([1,[1n,1n]]);expect(context.length()).toBe(2);
});


test('LLM NTL GotThem accepts and preserves certified factor order',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {GotThem}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const context=new FFTPrimeContext(),stream=new RandomStream(new Uint8Array(32));
 const B=[[1n,1n,1n,1n,0n,0n,0n,0n],[0n,0n,0n,0n,1n,1n,1n,1n]],W=[...Array.from({length:4},()=>[1n,1n]),...Array.from({length:4},()=>[-1n,1n])],f=[1n,0n,-4n,0n,6n,0n,-4n,0n,1n];
 expect(GotThem([],B,W,f,16,65537n,context,stream)).toEqual([1,[[1n,4n,6n,4n,1n],[1n,-4n,6n,-4n,1n]]]);expect(context.length()).toBe(3);
});
test('LLM NTL GotThem rejects small factor groups without changing the prefix',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {GotThem}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 expect(GotThem([[7n]],[[1n,1n,1n]],[[1n,1n],[1n,1n],[1n,1n]],[1n,3n,3n,1n],16,65537n,new FFTPrimeContext(),new RandomStream(new Uint8Array(32)))).toEqual([0,[[7n]]]);
});


test('LLM NTL cold small-prime factor order and following bytes',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {LocalInfoT,SmallPrimeFactorization}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 const context=new FFTPrimeContext(),stream=new RandomStream(Uint8Array.from({length:32},(_,i)=>i));
 expect(SmallPrimeFactorization(new LocalInfoT(),[-1n,0n,1n],stream,{context})).toEqual([[2n,1n],[1n,1n]]);expect(context.length()).toBe(1);expect(Array.from(stream.get(4))).toEqual([246,94,222,93]);
});
test('LLM NTL failed word initialization retains generated FFT primes',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {zz_pXModulus}=await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
 const state={context:new FFTPrimeContext(),stream:new RandomStream(new Uint8Array(32))};expect(()=>new zz_pXModulus(null,3n,{maxroot:240,state})).toThrow('zz_pInit: too many primes');expect(state.context.length()).toBe(5);expect(new zz_pXModulus(null,3n,{state}).PrimeCnt).toBe(1);
});


test('LLM NTL local factor updates preserve shared initialization bytes',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {LocalInfoT,UpdateLocalInfo}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const state={context:new FFTPrimeContext(),stream:new RandomStream(new Uint8Array(32))},info=new LocalInfoT();info.n=2;info.PossibleDegrees=7n;
 expect(UpdateLocalInfo(info,[],[],[-1n,0n,1n],1,125n,{van_Hoeij:0,MaxNumPrimes:3,state})).toBeNull();expect(info.p).toEqual([3n]);expect(state.context.length()).toBe(1);expect(Array.from(state.stream.get(4))).toEqual([152,186,151,124]);
});
test('LLM NTL cardinality search shares exact-division state',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {LocalInfoT,CardinalitySearch}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const state={context:new FFTPrimeContext(),stream:new RandomStream(new Uint8Array(32))},info=new LocalInfoT();info.n=20;info.PossibleDegrees=(1n<<21n)-1n;
 const g=[1n,...Array<bigint>(9).fill(0n),1n],h=[2n,...Array<bigint>(9).fill(0n),1n],f=[2n,...Array<bigint>(9).fill(0n),3n,...Array<bigint>(9).fill(0n),1n];expect(CardinalitySearch([],f,[g,h],info,1,100,65537n,{state})).toEqual([[g],h,[h]]);expect(state.context.length()).toBe(2);
});
test('LLM NTL complementary factor products share FFT state',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {InvMul}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const state={context:new FFTPrimeContext(),stream:new RandomStream(new Uint8Array(32))},W=Array.from({length:4},()=>Array<bigint>(101).fill(1n));const f=InvMul(W,[],65537n,state);expect([f.length,f[0],f[400]]).toEqual([401,1n,1n]);expect(state.context.length()).toBe(2);
});


test('LLM NTL homomorphic integer product initializes shared primes',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {mul}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZX1.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const state={context:new FFTPrimeContext(),stream:new RandomStream(new Uint8Array(32))},a=Array<bigint>(150).fill(1n),c=mul(a,a,state);expect([c.length,c[0],c[149],c[298]]).toEqual([299,1n,150n,1n]);expect(state.context.length()).toBe(1);
});
test('LLM NTL SS integer square does not initialize FFT primes',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {sqr}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZX1.js'));const {FFTPrimeContext}=await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));const {RandomStream}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));const state={context:new FFTPrimeContext(),stream:new RandomStream(new Uint8Array(32))};expect(sqr(Array<bigint>(80).fill(1n<<1919n),state).length).toBe(159);expect(state.context.length()).toBe(0);
});
test('LLM NTL integer transform-choice thresholds',async()=>{
 const {createRequire}=await import('node:module');const require=createRequire(new URL('../packages/sagemath-ts/package.json',import.meta.url));const {ChooseSS}=await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZX1.js'));expect([ChooseSS(149,1,149,1),ChooseSS(79,1920,79,1920),ChooseSS(79,3400,8192,3400),ChooseSS(79,6000,16384,6000),ChooseSS(79,10000,65536,10000),ChooseSS(79,17500,65536,17500)]).toEqual([false,true,true,true,true,true]);
});


test('LLM ntl-truncated-transform-roundtrip', async () => {
  const { createRequire } = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const { RandomStream } = await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
  const { FFTPrimeContext, UseFFTPrime, FFTFwd_trunc, FFTRev1_trunc } = await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));
  const context=new FFTPrimeContext();UseFFTPrime(0,context,new RandomStream(new Uint8Array(32)));
  const a=[1n,2n,3n,4n],info=context.get(0),frequencies=FFTFwd_trunc(a,2,info,4,4);
  expect(FFTRev1_trunc(frequencies,2,info,4)).toEqual([1n,2n,3n,4n]);
});
test('LLM ntl-transform-prefix-rounding', async () => {
  const { createRequire } = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const { FFTRoundUp } = await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT_impl.js'));
  expect(FFTRoundUp(401,9)).toBe(416);
});
test('LLM ntl-word-failed-rebuild-remainder', async () => {
  const { createRequire } = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const { zz_pXModulus, build, rem } = await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));
  const good=[1n,...Array<bigint>(47).fill(0n),1n],F=new zz_pXModulus(good,6n);
  expect(()=>build(F,[...good.slice(0,48),2n])).toThrow('InvMod: inverse undefined');
  const a=[...Array<bigint>(94).fill(0n),1n];
  expect(rem(a,F).flatMap((x:bigint,i:number)=>x?[[i,x]]:[])).toEqual([[30,5n],[46,5n]]);
});


test('LLM ntl-word-exact-crt-centering', async () => {
  const { createRequire } = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const { zz_pXModulus, FromModularRep } = await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));
  expect(FromModularRep([[441352763482308608n]],new zz_pXModulus(null,6n))).toEqual([2n]);
});

test('LLM NTL transposed transforms', async () => {
  const { createRequire } = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const { RandomStream } = await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
  const { FFTPrimeContext, UseFFTPrime, FFTFwd_trans, FFTRev1_trans } = await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));
  const context = new FFTPrimeContext();
  UseFFTPrime(0, context, new RandomStream(new Uint8Array(32)));
  const info = context.get(0);
  expect(FFTFwd_trans([1n, 0n, 0n, 0n], 2, info)).toEqual([1n, 1n, 1n, 1n]);
  expect(FFTRev1_trans([1n, 1n, 1n, 1n], 2, info)).toEqual([1n, 0n, 0n, 0n]);
});
test('LLM NTL projection retains transforms after a plain rebuild', async () => {
  const { createRequire } = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const { zz_pXModulus, zz_pXMultiplier, build } = await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));
  const { UpdateMap } = await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX1.js'));
  const F = new zz_pXModulus([1n, ...Array<bigint>(47).fill(0n), 1n], 6n);
  const B = new zz_pXMultiplier([2n, ...Array<bigint>(46).fill(0n), 3n], F);
  build(F, [1n, 1n]);
  expect(UpdateMap([1n], B, F)).toEqual([2n]);
});

test('LLM ntl-big-quotient-arithmetic', async () => {
  const { createRequire } = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const { ZZ_pXModulus, ZZ_pXMultiplier, rem, MulMod, SqrMod } = await import(
    require.resolve('@sagemath-ts/ntl-ts/src/ZZ_pX.js')
  );
  const F = new ZZ_pXModulus([1n, 0n, 1n], 7n);
  const B = new ZZ_pXMultiplier([3n, 1n], F);
  expect(rem([0n, 0n, 0n, 1n], F)).toEqual([0n, 6n]);
  expect(MulMod([1n, 2n], B, F)).toEqual([1n]);
  expect(SqrMod([1n, 2n], F)).toEqual([4n, 4n]);
});

test('LLM ntl-big-quotient-crt', async () => {
  const { createRequire } = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const { ZZ_pXModulus, FromModularRep } = await import(
    require.resolve('@sagemath-ts/ntl-ts/src/ZZ_pX.js')
  );
  const F = new ZZ_pXModulus(null, 6n);
  expect(FromModularRep([[441352763482308608n]], F)).toEqual([3n]);
});

test('LLM ordinary word products and optional inverse context', async () => {
  const {createRequire} = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const {mul, sqr, InvMod, InvModStatus} = await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));
  expect(mul([1n,2n], [3n,4n], 5n)).toEqual([3n,0n,3n]);
  expect(sqr([1n,2n], 5n)).toEqual([1n,4n,4n]);
  expect(InvMod([1n,1n], [1n,0n,1n], 5n, {maxroot:0})).toEqual([3n,2n]);
  expect(InvModStatus([1n,1n], [1n,0n,1n], 5n, {maxroot:0})).toEqual([0,[3n,2n]]);
});
test('LLM word extended GCD coefficients and two zero polynomials', async () => {
  const {createRequire} = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const {XGCD} = await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX1.js'));
  expect(XGCD([1n,0n,1n], [1n,1n], 5n)).toEqual([[1n],[3n],[3n,2n]]);
  expect(XGCD([], [], 5n)).toEqual([[],[1n],[]]);
});

test('LLM Hensel lifting retains shared FFT state across precision stages', async () => {
  const {createRequire} = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const {MultiLift} = await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
  const {FFTPrimeContext} = await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));
  const {RandomStream} = await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
  const state = {context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32))};
  const a = [1n,...Array<bigint>(21).fill(0n),1n], b = [2n,...Array<bigint>(21).fill(0n),1n];
  const f = [2n,...Array<bigint>(21).fill(0n),3n,...Array<bigint>(21).fill(0n),1n];
  expect(MultiLift([a,b], f, 9, 17n, {state})).toEqual([a,b]);
  expect(state.context.length()).toBe(2);
});
test('LLM additional lifting initializes its degree-dependent context', async () => {
  const {createRequire} = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const {AdditionalLifting} = await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
  const {FFTPrimeContext} = await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));
  const {RandomStream} = await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
  const state = {context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32))};
  expect(AdditionalLifting(5n, 1, [[1n,1n],[2n,1n]], 5n, 4, [7n,3n,1n], false, state)).toEqual([625n,4,[[281n,1n],[347n,1n]]]);
  expect(state.context.length()).toBe(1);
});

test('LLM FFT-prime coefficient context', async () => {
  const {createRequire} = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const {zz_pXModulus,rem} = await import(require.resolve('@sagemath-ts/ntl-ts/src/lzz_pX.js'));
  const p = 882705526964617217n, F = new zz_pXModulus([1n,0n,1n], p, {fftPrime:0});
  expect([F.PrimeCnt,F.NumPrimes,F.MaxRoot]).toEqual([0,1,25]);
  expect(rem([0n,0n,1n], F)).toEqual([882705526964617216n]);
});
test('LLM integer GCD and squarefree shared state', async () => {
  const {createRequire} = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const {ZZX_GCD,ZZX_SquareFreeDecomp} = await import(require.resolve('@sagemath-ts/ntl-ts'));
  const {FFTPrimeContext} = await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));
  const {RandomStream} = await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
  const state = {context:new FFTPrimeContext(), stream:new RandomStream(new Uint8Array(32))};
  expect(ZZX_GCD([-2n,0n,2n], [2n,-4n,2n], state)).toEqual([-2n,2n]);
  expect(ZZX_SquareFreeDecomp([1n,4n,4n], state)).toEqual([[[1n,2n],2]]);
  expect(state.context.length()).toBe(2);
});

test('LLM complete NTL factorization preserves content, multiplicity and order', async () => {
  const {createRequire} = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const {ZZX_factor,ZZX_SFFactor} = await import(require.resolve('@sagemath-ts/ntl-ts'));
  const {FFTPrimeContext} = await import(require.resolve('@sagemath-ts/ntl-ts/src/FFT.js'));
  const {RandomStream} = await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZ.js'));
  const state = {context:new FFTPrimeContext(),stream:new RandomStream(new Uint8Array(32))};
  expect(ZZX_factor([-2n,-4n,-2n],0,{state})).toEqual([-2n,[[[1n,1n],2]]]);
  expect(state.context.length()).toBe(2);
  expect(ZZX_SFFactor([-1n,0n,0n,0n,1n])).toEqual([[1n,0n,1n],[1n,1n],[-1n,1n]]);
});
test('LLM NTL primitive driver and full word-prime precision helpers', async () => {
  const {createRequire} = await import('node:module');
  const require = createRequire(new URL('../packages/sagemath-ts/package.json', import.meta.url));
  const {ll_SFFactor,Compute_pb,Compute_pdelta,Compute_pb_eff} = await import(require.resolve('@sagemath-ts/ntl-ts/src/ZZXFactoring.js'));
  expect(ll_SFFactor([1n,3n,2n])).toEqual([[1n,1n],[1n,2n]]);
  const p = 882705526964617217n;
  expect(Compute_pb([],[],p,1,p*p,4)).toEqual([[3],[p**3n]]);
  expect(Compute_pdelta(0,1n,p,120)).toEqual([3,p**3n]);
  expect(Compute_pb_eff(p,1,p*p,4,120)).toEqual([5,p**5n]);
});
test('LLM number-field constructor uses the NTL degree range', async () => {
  const {RationalPolynomial,NumberField} = await import('sagemath-ts/rings');
  const f = RationalPolynomial.fromBigInts([-2n,...Array<bigint>(300).fill(0n),1n]);
  expect(f.isIrreducible()).toBe(true);
  expect(new NumberField(f,'a').degree()).toBe(301);
});

test('LLM PARI conjugate bound and word derivative helpers', async () => {
  const {numberofconjugates,Flx_deriv,Flx_is_squarefree} = await import('@sagemath-ts/parigp-ts');
  expect(numberofconjugates([-2n,0n,0n,1n])).toBe(1);
  expect(numberofconjugates([-2n,0n,0n,1n],-1n)).toBe(3);
  expect(Flx_deriv([1n,0n,1n],3n)).toEqual([0n,2n]);
  expect(Flx_is_squarefree([1n,0n,1n],2n)).toBe(false);
  expect(Flx_is_squarefree([1n],2n)).toBe(true);
});


test('LLM.md — PARI finite-field universal comparison', async () => {
  const { cmp_universal } = await import('@sagemath-ts/parigp-ts/src/gen2.js');
  const { PariType } = await import('@sagemath-ts/parigp-ts/src/types.js');
  const ff = { type: PariType.t_FFELT as const, p: 7n, degree: 2,
    definingPoly: [3n, 6n, 1n] };
  expect(cmp_universal({ ...ff, value: [1n, 5n] }, { ...ff, value: [2n, 1n] })).toBe(-1);
  expect(cmp_universal({ ...ff, value: [2n, 1n] }, { ...ff, value: [1n, 5n] })).toBe(1);
  expect(cmp_universal({ ...ff, value: [0n] }, { ...ff, value: [] })).toBe(0);
});


test('LLM.md — polynomial distinct-root options', async () => {
  const { PolynomialRing } = await import('sagemath-ts/rings');
  const rootRing = new PolynomialRing(GF(7n), 'x');
  const rootX = rootRing.gen();
  const rootPolynomial = rootX.pow(2n).mul(rootX.sub(rootRing.one()));
  expect(rootPolynomial.roots().map(([r, m]) => [String(r), m])).toEqual([['1', 1], ['0', 2]]);
  expect(rootPolynomial.roots({ multiplicities: false }).map(String)).toEqual(['0', '1']);
});

test('LLM.md — modular roots and cached field parents', async () => {
  const {PolynomialRing} = await import('sagemath-ts/rings');
  const K = Zmod(7n), R = new PolynomialRing(K,'x');
  expect(String(K.field())).toBe('Finite Field of size 7');
  expect(K.field()).toBe(K.field());
  expect(K.factored_order()).toEqual([[7n,1n]]);
  const f = R.__call__([0n,-1n,1n]);
  expect(f.roots().map(([r,m])=>[String(r),m,r.parent===K.field()])).toEqual([['0',1,true],['1',1,true]]);
  expect(f.roots({multiplicities:false}).map(r=>[String(r),r.parent===K])).toEqual([['0',true],['1',true]]);
  const S=new PolynomialRing(Zmod(8n),'x');
  expect(S.__call__([-1n,0n,1n]).roots({multiplicities:false}).map(String)).toEqual(['1','5','3','7']);
  const T=new PolynomialRing(Zmod(4n),'x');
  expect(T.__call__([2n,2n]).roots({multiplicities:false}).map(String)).toEqual(['0','2']);
});

test('LLM.md — extension-field Weierstrass isomorphisms', async () => {
  const { GFpn } = await import('sagemath-ts/rings/finite_rings');
  const { EllipticCurveGeneric, _isomorphisms, baseWI } = await import('sagemath-ts/schemes/elliptic_curves');
  const K = GFpn(2n, 2, [1, 1], 'a');
  const E = new EllipticCurveGeneric(K, [K.zero(), K.zero(), K.one(), K.zero(), K.zero()]);
  const automorphisms = [..._isomorphisms(E, E)];
  expect(automorphisms.length).toBe(24);
  expect(automorphisms[0]!.map(String)).toEqual(['1', '0', '0', '0']);
  expect(automorphisms.every(t => new baseWI(...t).call(E.a_invariants()).every((v, i) => v.eq(E.a_invariants()[i]!)))).toBe(true);
});

test('LLM.md — generic-curve isomorphism ordering and comparison', async () => {
  const { GFpn } = await import('sagemath-ts/rings/finite_rings');
  const { EllipticCurveGeneric, WeierstrassIsomorphism } = await import('sagemath-ts/schemes/elliptic_curves');
  const K = GFpn(2n, 2, [1, 1], 'a');
  const C = new EllipticCurveGeneric(K, [K.zero(), K.zero(), K.one(), K.zero(), K.zero()]);
  const autos = C.automorphisms();
  expect(autos.slice(0, 2).map(t => t.map(String))).toEqual([['1', '0', '0', '0'], ['1', '0', '0', '1']]);
  expect(C.is_isomorphic(C)).toBe(true);
  expect(C.isomorphism_to(C).map(String)).toEqual(['1', '0', '0', '0']);
  const a = new WeierstrassIsomorphism(C, autos[0], C), b = new WeierstrassIsomorphism(C, autos[1], C);
  expect(WeierstrassIsomorphism._comparison_impl(a, b, 'lt')).toBe(true);
  expect(WeierstrassIsomorphism._comparison_impl(null, a, 'eq')).toBe(null);
});


test('documented extension scalar roots preserve native order', async () => {
  const { GFpn } = await import('sagemath-ts/rings/finite_rings');
  const { setrand, getrand } = await import('../packages/parigp-ts/src/index.js');
  const state = getrand();
  try {
    const F = GFpn(3n, 2, [1, 0], 'a');
    setrand(17n);
    expect(F.__call__(2n).sqrt({ all: true }).map(String)).toEqual(['a', '2*a']);
  } finally {
    setrand(state);
  }
});


test('documented native extension square predicate and dependency norm',async()=>{
  const {GFpn}=await import('sagemath-ts/rings/finite_rings');
  const {FpX_resultant,FF_norm,FF_issquare,PariType,getrand}=await import('../packages/parigp-ts/src/index.js');
  const F=GFpn(3n,2,[1,0],'a'),state=getrand();
  expect(F.gen().is_square()).toBe(true);
  expect(F.gen().add(F.one()).is_square()).toBe(false);
  expect(FpX_resultant([1n,0n,1n],[0n,1n],3n)).toBe(1n);
  const a={type:PariType.t_FFELT as PariType.t_FFELT,p:3n,degree:2,value:[1n,1n],definingPoly:[1n,0n,1n]};
  expect(FF_norm(a)).toBe(2n);expect(FF_issquare(a)).toBe(false);
  expect(getrand()).toBe(state);
});


test('LLM.md extension curve coordinate lifting', async () => {
  const { GFpn } = await import('sagemath-ts/rings/finite_rings');
  const { EllipticCurveGeneric } = await import('sagemath-ts/schemes/elliptic_curves');
  const K = GFpn(3n, 2, [1, 0], 'a');
  const E = new EllipticCurveGeneric(K, [K.zero(), K.zero(), K.zero(), K.one(), K.zero()]);
  expect(E.is_x_coord(K.one())).toBe(true);
  expect(E.lift_x(K.one(), true).map(P => String(P.y()))).toEqual(['a', '2*a']);
});


test('documented optimized curve coordinate adapters', () => {
  const F = GF(3n);
  const E = EllipticCurve(F, [1n, 0n]);
  expect(E.is_x_coord('0')).toBe(true);
  expect(E.lift_x(0n).curve).toBe(E);
  expect(E.lift_x(1n, true)).toEqual([]);
  expect(E.lift_x(1n, true, true).map(String)).toEqual(['(1 : y : 1)', '(1 : 2*y : 1)']);
  expect(() => E.lift_x(QQ.zero())).toThrow('Unable to construct a point');
  expect(() => E.lift_x('0')).toThrow("'str' object has no attribute 'parent'");
});

test('documented generic y-extension and named polynomial APIs', async () => {
  const { EllipticCurveGeneric } = await import('sagemath-ts/schemes/elliptic_curves');
  const { RationalPolynomial, NumberField } = await import('sagemath-ts/rings');
  const F = GF(3n);
  const E = new EllipticCurveGeneric(F, [F.zero(), F.zero(), F.zero(), F.one(), F.one()]);
  expect(E.lift_x(F.__call__(2n), true)).toEqual([]);
  const points = E.lift_x(F.__call__(2n), true, true);
  expect(points.map(P => String(P.y()))).toEqual(['y', '2*y']);
  expect(String(points[0]!.curve.base_ring)).toBe('Finite Field in y of size 3^2');
  expect(points[0]!.curve.is_on_curve(points[0]!.x(), points[0]!.y())).toBe(true);
  const f = new RationalPolynomial([new Rational(-2n), Rational.zero(), Rational.one()], 'y');
  expect(String(f)).toBe('y^2 - 2');
  expect(String(f.derivative())).toBe('2*y');
  expect(f.scale(Rational.zero()).variableName).toBe('y');
  const K = new NumberField(f, 'a');
  expect(String(K)).toBe('Number Field in a with defining polynomial y^2 - 2');
  expect(K.zero().isZero()).toBe(true);
  expect(K.one().isZero()).toBe(false);
});


test('documented finite-extension and native trace APIs', async () => {
  const { GFpn } = await import('sagemath-ts/rings/finite_rings');
  const { FF_trace, FpXQ_trace, Flxq_trace, F2xq_trace, PariType } = await import('@sagemath-ts/parigp-ts');
  const F = GFpn(3n, 2, [1, 0], 'a');
  expect(F.gen().add(F.one()).trace().toString()).toBe('2');
  expect(F.gen().trace().parent).toBe(F.baseField);
  expect(FpXQ_trace([1n, 1n], [1n, 0n, 1n], 3n)).toBe(2n);
  expect(Flxq_trace([1n, 1n], [1n, 0n, 1n], 3n)).toBe(2n);
  expect(F2xq_trace(3n, 11n)).toBe(1n);
  expect(FF_trace({type: PariType.t_FFELT, p: 3n, degree: 2, value: [1n, 1n], definingPoly: [1n, 0n, 1n]})).toBe(2n);
});


test('documented finite-extension characteristic polynomial and native resultant APIs', async () => {
  const { GFpn } = await import('sagemath-ts/rings/finite_rings');
  const { FpX_FpXY_resultant, FpV_polint } = await import('@sagemath-ts/parigp-ts');
  const F = GFpn(3n, 2, [1, 0], 'a');
  expect(F.gen().charpoly('y').toString()).toBe('y^2 + 1');
  expect(F.gen().norm().toString()).toBe('1');
  expect(F.gen().norm().parent).toBe(F.baseField);
  expect(F.one().charpoly().toString()).toBe('x^2 + x + 1');
  expect(
    FpX_FpXY_resultant(
      [1n, 0n, 2n],
      [
        [1n, 1n],
        [0n, 1n],
      ],
      5n
    )
  ).toEqual([2n, 4n, 3n]);
  expect(FpV_polint([0n, 1n, 2n], [1n, 3n, 7n], 17n)).toEqual([1n, 1n, 1n]);
});


test('documented torsion group-relation dependency and additive order', async () => {
  const { Mod, groups } = await import('sagemath-ts');
  expect(Mod(2n, 20n).additive_order()).toBe(10n);
  expect(Mod(2n, 90308402384902n).additive_order()).toBe(45154201192451n);
  expect(groups.linear_relation(Mod(2n, 20n), Mod(4n, 20n))).toEqual([2n, 1n]);
  expect(groups.linear_relation(Mod(0n, 20n), Mod(2n, 20n))).toEqual([10n, 0n]);
});


test('LLM.md — cached finite-curve points and group', async () => {
  const { EllipticCurve, GF } = await import('sagemath-ts');
  const E = EllipticCurve(GF(5n), [1n, 3n]);
  expect(E.points().map(String)).toEqual(['(0 : 1 : 0)', '(1 : 0 : 1)', '(4 : 1 : 1)', '(4 : 4 : 1)']);
  expect(E.points()).toBe(E.points());
  expect(E.abelian_group().invariants).toEqual([4n]);
  expect(E.abelian_group()).toBe(E.abelian_group());
});


test('LLM.md — default five-coefficient model preservation', async () => {
  const { EllipticCurve, GF } = await import('sagemath-ts');
  const K = GF(2n);
  const E = EllipticCurve(K, [1n, 0n, 0n, 0n, 1n]);
  expect(E.a_invariants().map(String)).toEqual(['1', '0', '0', '0', '1']);
  expect(E.point([K.zero(), K.one()]).mul(2n).toString()).toBe('(0 : 1 : 0)');
});


test('LLM.md — general-model PARI scalar multiplication', async () => {
  const { EllipticCurve, GF, Integer } = await import('sagemath-ts');
  const { ellmul, ell_to_a4a6_bc } = await import('@sagemath-ts/parigp-ts');
  const K = GF(11n);
  const E = EllipticCurve(K, [1n, 0n, 1n, 0n, 0n]);
  expect(E.pari_curve()).toBe(E.__pari__());
  const model = E.pari_curve();
  if ('field' in model) throw new Error('expected a prime-field model');
  expect(ell_to_a4a6_bc(model, 11n)).toEqual([5n, 6n, [6n, 3n, 3n, 9n]]);
  expect(ellmul(model, {isInfinity:false,x:0n,y:0n}, 2n)).toEqual({isInfinity:false,x:0n,y:10n});
  expect(E.point([K.zero(),K.zero()]).mul(new Integer(2n)).toString()).toBe('(0 : 10 : 1)');
});


test('LLM.md — general-model PARI orders and bounded hybrid factorization', async () => {
  const { ellcard, ellorder } = await import('@sagemath-ts/parigp-ts');
  const { factor_trial_division } = await import('sagemath-ts/rings');
  const K = GF(11n);
  const E = EllipticCurve(K, [1n, 0n, 1n, 0n, 0n]);
  const model = E.pari_curve();
  if ('field' in model) throw new Error('expected a prime-field model');
  expect(ellcard(model)).toBe(6n);
  expect(ellorder(model, {isInfinity: false, x: 0n, y: 0n})).toBe(3n);
  expect(E.point([K.zero(), K.zero()]).order({algorithm: 'pari'})).toBe(3n);
  expect(new Integer(143n).factor({limit: 9n})).toEqual([[143n, 1n]]);
  expect(new Integer(143n).factor({limit: new Integer(12n)})).toEqual([[11n, 1n], [13n, 1n]]);
  expect(factor_trial_division(-143n, 12n)).toEqual([[-1n, 1n], [11n, 1n], [13n, 1n]]);
  const C = EllipticCurve(K, [0n, 0n, 0n, 1n, 1n]);
  expect(C.cardinality_pari()).toBe(14n);
  expect(C.cardinality()).toBe(14n);
  expect(C.order()).toBe(14n);
  expect(C.point([K.zero(), K.one()]).order({algorithm: 'hybrid'})).toBe(7n);
});

test('LLM.md — optimized finite-point algorithms and shared cache', () => {
  const K = GF(11n);
  const E = EllipticCurve(K, [1n, 1n]);
  expect(E.point(0n, 1n).order({algorithm: 'generic_small'})).toBe(7n);
  expect(E.point(0n, 1n).additive_order({algorithm: 'hybrid'})).toBe(7n);
  expect(E.point(0n, 1n).order({algorithm: 'pari'})).toBe(7n);
  expect(E._order).toBe(14n);
});


test('LLM.md — optimized scalar coercion and cached PARI representation', () => {
  const E = EllipticCurve(GF(11n), [1n, 1n]);
  const P = E.point(0n, 1n);
  expect(P.order()).toBe(7n);
  const Q = P.mul(new Integer(2n));
  expect(Q.toString()).toBe('(3 : 3 : 1)');
  expect((Q as unknown as {_order: bigint | null})._order).toBe(7n);
  expect(P.mul(2).toString()).toBe('(3 : 3 : 1)');
  expect(E.pari_curve()).toBe(E.__pari__());
  expect(E.pari_curve()).toBe(E.toPari());
});


test('LLM.md — native binary elliptic dependency kernels', async () => {
  const { F2xq_inv, F2xq_invsafe, F2xqE_dbl, F2xqE_mul } = await import('@sagemath-ts/parigp-ts');
  const P = {isInfinity: false as const, x: 1n, y: 2n};
  expect(F2xq_inv(2n, 7n)).toBe(3n);
  expect(F2xq_invsafe(3n, 5n)).toBeNull();
  expect(F2xqE_dbl(P, 1n, 7n)).toEqual({isInfinity: false, x: 0n, y: 1n});
  expect(F2xqE_mul(P, 4n, 1n, 7n)).toEqual({isInfinity: true});
});


test('LLM.md — binary model initialization and original coordinates', async () => {
  const { PariType, ellinit_Fq, ellmul, FF_ellmul } = await import('@sagemath-ts/parigp-ts');
  const field = {type: PariType.t_FFELT as const, p: 2n, degree: 2,
    definingPoly: [1n, 1n, 1n], value: [0n, 1n]};
  const E = ellinit_Fq([1n, 1n, 0n, 0n, 1n], field)!;
  const P = {isInfinity: false as const, x: 1n, y: field};
  const Q = ellmul(E, P, 2n);
  expect(Q.isInfinity).toBe(false);
  if (Q.isInfinity) throw new Error('expected an affine point');
  expect([Q.x.value, Q.y.value]).toEqual([[], [1n]]);
  expect(FF_ellmul(E, P, 4n)).toEqual({isInfinity: true});
  expect(ellinit_Fq([0n, 0n, 0n, 1n, 1n], field)).toBeNull();
});

test('LLM.md — binary Sage scalar dispatch', async () => {
  const { FiniteFieldExtension } = await import('sagemath-ts/rings/finite_rings');
  const K = new FiniteFieldExtension(2n, 2, [1, 1], 'a');
  const E = EllipticCurve(K, [1n, 1n, 0n, 0n, 1n]);
  const P = E.point([K.one(), K.gen()]);
  expect(P.mul(2n).toString()).toBe('(0 : 1 : 1)');
  expect(P.mul(4n).is_zero()).toBe(true);
  expect(E.pari_curve()).toBe(E.__pari__());
});


test('LLM.md — odd-extension elliptic kernels', async () => {
  const { FlxqE_dbl, FlxqE_mul, FpXQE_mul } = await import('@sagemath-ts/parigp-ts');
  const T = [1n, 0n, 1n];
  const P = {isInfinity: false as const, x: [1n], y: [1n]};
  expect(FlxqE_dbl(P, [[1n]], T, 3n)).toEqual({isInfinity: false, x: [1n], y: [2n]});
  expect(FlxqE_mul(P, 3n, [[1n]], T, 3n)).toEqual({isInfinity: true});
  expect(FpXQE_mul({isInfinity: false, x: [1n], y: [2n]}, 2n, [2n], T, 7n))
    .toEqual({isInfinity: false, x: [], y: [1n]});
});

test('LLM.md — characteristic-three and odd-extension scalar callers', async () => {
  const { FiniteFieldExtension } = await import('sagemath-ts/rings/finite_rings');
  const K = GF(3n);
  const E = EllipticCurve(K, [1n, 0n, 0n, 1n, 1n]);
  expect(E.point([K.zero(), K.one()]).mul(2n).toString()).toBe('(0 : 2 : 1)');
  const L = new FiniteFieldExtension(3n, 2, [1, 0], 'a');
  const C = EllipticCurve(L, [0n, 0n, 0n, 1n, 1n]);
  const P = C.point([L.gen(), L.one()]);
  expect(P.mul(2n).toString()).toBe('(a + 1 : 0 : 1)');
  expect(P.mul(4n).is_zero()).toBe(true);
  expect(C.pari_curve()).toBe(C.__pari__());
  const model = C.pari_curve();
  if (!('field' in model) || !('oddModel' in model)) throw new Error('expected odd field model');
  expect(model.oddModel[0]).toEqual([1n]);
});

test('LLM.md — finite-field j-invariant and short initialization', async () => {
  const { PariType, ellinit_Fq } = await import('@sagemath-ts/parigp-ts');
  const field = {type: PariType.t_FFELT as const, p: 3n, degree: 2,
    definingPoly: [1n, 0n, 1n], value: [0n, 1n]};
  const E = ellinit_Fq([field], field)!;
  expect(E.a2.value).toEqual([0n, 1n]);
  expect(E.a6.value).toEqual([1n]);
  expect(E.j.value).toEqual([0n, 1n]);
  expect(ellinit_Fq([1n, 1n], field)!.a4.value).toEqual([1n]);
  expect(ellinit_Fq([0n, 0n], field)).toBeNull();
});

test('LLM.md — small-prime initialization and field accessors', async () => {
  const { ellinit, ellj, elldisc, ellcoeffs, ellisnonsingular, ellmul } = await import('@sagemath-ts/parigp-ts');
  const E = ellinit([1n, 0n, 0n, 1n, 1n], 3n)!;
  expect(ellj(E).value).toEqual([2n]);
  expect(elldisc(E).value).toEqual([2n]);
  expect(ellcoeffs(E)[0].value).toEqual([1n]);
  expect(ellisnonsingular(E)).toBe(true);
  const Q = ellmul(E, {isInfinity: false, x: 0n, y: 1n}, 2n);
  expect(Q.isInfinity).toBe(false);
  if (Q.isInfinity) throw new Error('expected affine point');
  expect(Q.y.value).toEqual([2n]);
  expect(ellinit([0n, 0n], 3n)).toBeNull();
});

test('LLM.md — supplied-bound finite-field point orders', async () => {
  const { ellinit, FF_ellorder, gen_order } = await import('@sagemath-ts/parigp-ts');
  const E = ellinit([1n, 0n, 0n, 1n, 1n], 3n)!;
  const P = {isInfinity: false as const, x: 0n, y: 1n};
  expect(FF_ellorder(E, P, 12n)).toBe(3n);
  expect(FF_ellorder(E, P, [[2n, 2n], [3n, 1n]])).toBe(3n);
  expect(FF_ellorder(E, P, [12n, [[2n, 2n], [3n, 1n]]])).toBe(3n);
  expect(gen_order(4n, 12n, (x, n) => x * n % 12n, x => x === 0n)).toBe(3n);
});

test('LLM.md — base-field extension traces and cardinality', async () => {
  const { elltrace_extension, Fp_ffellcard } = await import('@sagemath-ts/parigp-ts');
  expect(elltrace_extension(3n, 2, 7n)).toBe(-5n);
  expect(elltrace_extension(3n, 0, 7n)).toBe(2n);
  expect(Fp_ffellcard(1n, 1n, 49n, 2, 7n)).toBe(55n);
});

test('LLM.md — word-field point kernels', async () => {
  const { Fle_mulu, Fle_order, Flj_dbl_pre } = await import('@sagemath-ts/parigp-ts');
  const P = {isInfinity: false as const, x: 0n, y: 1n};
  expect(Fle_order(P, 5n, 1n, 7n)).toBe(5n);
  expect(Fle_mulu(P, 5n, 1n, 7n).isInfinity).toBe(true);
  expect(Flj_dbl_pre({X: 3n, Y: 4n, Z: 0n}, 1n, 7n, 0n).X).toBe(3n);
});

test('LLM.md — direct extension square roots', async () => {
  const { F2xq_sqrt, Flxq_sqrt, Fl2_sqrt_pre } = await import('@sagemath-ts/parigp-ts');
  expect(F2xq_sqrt(2n, 11n)).toBe(6n);
  expect(Flxq_sqrt([2n], [1n,0n,1n], 3n)).toEqual([0n,1n]);
  expect(Fl2_sqrt_pre([2n,0n], 2n, 3n, 0n)).toEqual([0n,1n]);
});

test('LLM.md — ternary supersingular count corrects the upstream twist sign', async () => {
  const { F3xq_ellcardj, setrand } = await import('@sagemath-ts/parigp-ts');
  setrand(4n);
  expect(F3xq_ellcardj([0n,0n,2n], [2n,1n,2n], [1n,2n,0n,1n], 27n, 3)).toBe(19n);
});

test('LLM.md — constant-j extension counts', async () => {
  const { Flxq_ellcardj, FpXQ_ellcardj } = await import('@sagemath-ts/parigp-ts');
  for (const count of [Flxq_ellcardj, FpXQ_ellcardj])
    expect(count([], [1n,1n], 0n, [2n,4n,1n], 25n, 5n, 2)).toBe(21n);
});


test('LLM.md — extension-field Shanks count', async () => {
  const { Flxq_ellcard_Shanks, setrand } = await import('@sagemath-ts/parigp-ts');
  setrand(1n);
  expect(Flxq_ellcard_Shanks([1n,1n], [2n,1n], 289n, [3n,16n,1n], 17n)).toBe(306n);
});


test('LLM.md — polynomial p-adic counting dependencies', async () => {
  const { ZpXQ_inv, ZpXQ_invlift, ZpXQ_div, Flx_Teichmuller } = await import('@sagemath-ts/parigp-ts');
  const T = [1n,0n,1n];
  expect(ZpXQ_inv([1n,1n], T, 3n, 5)).toEqual([122n,121n]);
  expect(ZpXQ_invlift([1n,1n], [2n,1n], T, 3n, 5)).toEqual([122n,121n]);
  expect(ZpXQ_div([1n], [1n,1n], T, 243n, 3n, 5)).toEqual([122n,121n]);
  expect(Flx_Teichmuller([2n,4n,1n], 5n, 3)).toEqual([57n,89n,1n]);
});


test('LLM.md — scalar lifts and cyclotomic counting dependencies', async () => {
  const { Zp_sqrtnlift, Zp_sqrtlift, ZpXQ_frob_cyc, ZpXQ_norm_pcyc, ZpXQ_sqrtnorm_pcyc } = await import('@sagemath-ts/parigp-ts');
  expect(Zp_sqrtnlift(13n,3n,2n,5n,3)).toBe(67n);
  expect(Zp_sqrtlift(4n,2n,5n,3)).toBe(2n);
  expect(Zp_sqrtnlift(783n,1n,1n,2n,13)).toBe(783n);
  const T=[1n,1n,1n];
  expect(ZpXQ_frob_cyc([0n,1n],T,125n,5n)).toEqual([124n,124n]);
  expect(ZpXQ_norm_pcyc([2n,1n],T,125n,5n)).toBe(3n);
  expect(ZpXQ_sqrtnorm_pcyc([3n,3n],T,125n,5n,3)).toBe(122n);
});


test('LLM.md — p-adic exponential and quotient logarithm', async () => {
  const { Zp_inv, Zp_invlift, Zp_div, Zp_exp, ZpXQ_log } = await import('@sagemath-ts/parigp-ts');
  expect(Zp_inv(2n, 5n, 3)).toBe(63n);
  expect(Zp_invlift(2n, 3n, 5n, 3)).toBe(63n);
  expect(Zp_div(7n, 2n, 5n, 1)).toBe(3n);
  expect(Zp_exp(3n, 3n, 5)).toBe(229n);
  expect(ZpXQ_log([4n,3n], [1n,0n,1n], 3n, 5)).toEqual([66n,12n]);
});
