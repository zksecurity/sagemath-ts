/** Rational arithmetic compared with Sage; numeric results use exact binary64 bits. */
import { Matrix } from '../../../../packages/sagemath-ts/src/matrix/matrix_generic.js';
import { norm } from '../../../../packages/sagemath-ts/src/matrix/matrix_operations.js';
import { Integer } from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';
import { QQ, RationalField } from '../../../../packages/sagemath-ts/src/rings/rational_field.js';
const unary = [
  'abs',
  'neg',
  'inv',
  'toString',
  'floor',
  'ceil',
  'round',
  'trunc',
  'isZero',
  'isOne',
  'isInteger',
  'isPositive',
  'isNegative',
  'asIntegerRatio',
  'height',
  'continued_fraction_list',
  'continued_fraction',
  'is_square',
  'period',
  'real',
  'imag',
  'conjugate',
  'norm',
  'relative_norm',
  'absolute_norm',
  'trace',
  'support',
  'is_unit',
  'is_integral',
  'additive_order',
  'multiplicative_order',
  'simplify',
  'list',
  'signAsInt',
];
const binary = [
  'add',
  'sub',
  'mul',
  'div',
  'eq',
  'lt',
  'le',
  'gt',
  'ge',
  'cmp',
  'rational_gcd',
  'rational_lcm',
  'content',
];
const integer = [
  'pow',
  'valuation',
  'ord',
  'is_nth_power',
  'nth_root',
  'val_unit',
  'padic_valuation',
  'denominator_valuation',
  'numerator_valuation',
  'str',
];
const set = ['prime_to_S_part', 'is_S_unit', 'is_S_integral'];
function normalize(value: unknown): unknown {
  if (typeof value === 'boolean') return value;
  if (Array.isArray(value)) return value.map(normalize);
  return String(value);
}
function compare(method: string, n: bigint, d: bigint, args: unknown[]): string {
  try {
    const x = new Rational(n, d);
    const fn = (x as unknown as Record<string, (...args: unknown[]) => unknown>)[method]!;
    return JSON.stringify({ value: normalize(fn.apply(x, args)) });
  } catch (error) {
    return JSON.stringify({ error: error instanceof Error ? error.name : String(error) });
  }
}
export const functions = Object.fromEntries([
  ...unary.map((method) => [
    `rational_${method}`,
    (n: bigint, d: bigint) => compare(method, n, d, []),
  ]),
  ...binary.map((method) => [
    `rational_${method}`,
    (n: bigint, d: bigint, a: bigint, b: bigint) => compare(method, n, d, [new Rational(a, b)]),
  ]),
  ...integer.map((method) => [
    `rational_${method}`,
    (n: bigint, d: bigint, a: bigint) => compare(method, n, d, [method === 'str' ? Number(a) : a]),
  ]),
  ...set.map((method) => [
    `rational_${method}`,
    (n: bigint, d: bigint, a: bigint[]) => compare(method, n, d, [a]),
  ]),
]);
for (const method of ['toNumber', 'n', 'numerical_approx', 'global_height', 'matrix_norm']) {
  functions[`rational_${method}`] = (n: bigint, d: bigint) => {
    try {
      const x = new Rational(n, d);
      let value: string;
      if (method === 'global_height') value = String(BigInt(Math.round(x.global_height() * 1e8)));
      else {
        const view = new DataView(new ArrayBuffer(8));
        view.setFloat64(
          0,
          method === 'matrix_norm'
            ? norm(new Matrix(QQ, 1, 1, [[x]]), 2)
            : x[method as 'toNumber' | 'n' | 'numerical_approx']()
        );
        value = view.getBigUint64(0).toString(16).padStart(16, '0');
      }
      return JSON.stringify({ value });
    } catch (error) {
      return JSON.stringify({ error: error instanceof Error ? error.name : String(error) });
    }
  };
}
functions.rational_sqrt = (n: bigint, d: bigint) => compare('sqrt', n, d, [{ extend: false }]);
functions.rational_gamma = (n: bigint, d: bigint) => compare('gamma', n, d, []);

for (const method of [
  'from',
  'qq_from',
  'fromTuple',
  'zero',
  'one',
  'constructor',
  'minpoly',
  'charpoly',
  'factorial',
  'roundToRational',
  'ndigits',
  'nbits',
  'local_height',
  'round_mode',
  'continued_fraction_hj',
]) {
  functions[`rational_${method}`] = (n: bigint, d: bigint, arg = 0n) => {
    try {
      let value: unknown;
      if (method === 'from') value = Rational.from(Number(n) / Number(d));
      else if (method === 'qq_from') value = QQ.__call__(Number(n) / Number(d));
      else if (method === 'fromTuple') value = Rational.fromTuple([n, d]);
      else if (method === 'constructor') value = new Rational(n, d);
      else if (method === 'zero') value = Rational.zero();
      else if (method === 'one') value = Rational.one();
      else {
        const x = new Rational(n, d);
        if (method === 'minpoly' || method === 'charpoly') value = x[method]().coeffs;
        else if (method === 'roundToRational') value = x.roundToRational(Number(arg));
        else if (method === 'ndigits') value = x.ndigits(arg);
        else if (method === 'nbits') value = x.nbits();
        else if (method === 'factorial') value = x.factorial();
        else if (method === 'local_height') value = BigInt(Math.round(x.local_height(arg) * 1e8));
        else if (method === 'round_mode')
          value = x.round((['even', 'away', 'toward', 'up', 'down', 'odd'] as const)[Number(arg)]);
        else if (method === 'continued_fraction_hj') value = x.continued_fraction_list('hj');
      }
      return JSON.stringify({ value: normalize(value) });
    } catch (error) {
      return JSON.stringify({ error: error instanceof Error ? error.name : String(error) });
    }
  };
}
functions.rational_fromString = (codes: bigint[]) => {
  try {
    return JSON.stringify({
      value: String(Rational.fromString(String.fromCodePoint(...codes.map(Number)))),
    });
  } catch (error) {
    return JSON.stringify({ error: error instanceof Error ? error.name : String(error) });
  }
};

for (const method of [
  'numerator',
  'numer',
  'denominator',
  'denom',
  'sign',
  'from_rational',
  'from_bigint',
]) {
  functions[`rational_${method}`] = (n: bigint, d: bigint) => {
    const x = new Rational(n, d);
    return JSON.stringify({
      value: String(
        method === 'from_rational'
          ? Rational.from(x)
          : method === 'from_bigint'
            ? Rational.from(n)
            : x[method as 'numerator']
      ),
    });
  };
}
for (const method of ['qq_fromString', 'from_string']) {
  functions[`rational_${method}`] = (codes: bigint[]) => {
    try {
      const str = String.fromCodePoint(...codes.map(Number));
      return JSON.stringify({
        value: String(method === 'qq_fromString' ? QQ.__call__(str) : Rational.from(str)),
      });
    } catch (error) {
      return JSON.stringify({ error: error instanceof Error ? error.name : String(error) });
    }
  };
}

functions.rational_valuation_error = (n: bigint, d: bigint, p: bigint, method: bigint) => {
  const names = [
    'valuation',
    'ord',
    'val_unit',
    'numerator_valuation',
    'denominator_valuation',
    'padic_valuation',
  ] as const;
  try {
    return JSON.stringify({ value: normalize(new Rational(n, d)[names[Number(method)]!](p)) });
  } catch (error) {
    return JSON.stringify({
      error: error instanceof Error ? error.name : String(error),
      message: error instanceof Error ? error.message : String(error),
    });
  }
};

function qqResult(run: () => unknown): string {
  try {
    return JSON.stringify({ value: normalize(run()) });
  } catch (error) {
    return JSON.stringify({ error: error instanceof Error ? error.name : String(error) });
  }
}
const qqMethods = [
  'getInstance',
  'zero',
  'one',
  'is_field',
  'characteristic',
  'is_ring',
  'is_integral_domain',
  'is_prime_field',
  'is_absolute',
  'is_finite',
  'degree',
  'absolute_degree',
  'ngens',
  'gens',
  'an_element',
  'some_elements',
  'discriminant',
  'absolute_discriminant',
  'relative_discriminant',
  'class_number',
  'signature',
  'order',
  'maximal_order',
  'ring_of_integers',
  'number_field',
  'power_basis',
  'toString',
  '_latex_',
];
for (const method of qqMethods) {
  functions[`qq_${method}`] = () =>
    qqResult(() => {
      if (method === 'getInstance') return RationalField.getInstance();
      if (method === 'some_elements') return [...QQ.some_elements()];
      return (QQ as unknown as Record<string, () => unknown>)[method]!();
    });
}
for (const method of [
  'gen',
  'zeta',
  'primes_of_bounded_norm_iter',
  'iterator',
  'range_by_height',
]) {
  functions[`qq_${method}`] = (a: bigint, b?: bigint) =>
    qqResult(() => {
      if (method === 'gen' || method === 'zeta') return QQ[method](Number(a));
      if (method === 'primes_of_bounded_norm_iter') return [...QQ.primes_of_bounded_norm_iter(a)];
      if (method === 'range_by_height') return [...QQ.range_by_height(a, b)];
      const it = QQ[Symbol.iterator]();
      return Array.from({ length: Number(a) }, () => it.next().value);
    });
}
for (const method of ['selmer_generators', 'selmer_group_iterator']) {
  functions[`qq_${method}`] = (S: bigint[], m: bigint) =>
    qqResult(() =>
      method === 'selmer_generators'
        ? QQ.selmer_generators(S, m, { orders: true })
        : [...QQ.selmer_group_iterator(S, m)]
    );
}
functions.qq_quadratic_defect = (n: bigint, d: bigint, p: bigint) =>
  qqResult(() => QQ.quadratic_defect(new Rational(n, d), p));
functions.qq_pair = (n: bigint, d: bigint, a: bigint, b: bigint, mode: bigint) =>
  qqResult(() => {
    const x = mode === 0n ? new Rational(n, d) : Number(n) / Number(d);
    const y = mode === 0n ? new Rational(a, b) : Number(a) / Number(b);
    return QQ.__call__([x, y]);
  });
functions.qq_prime_fraction_bound = (n: bigint, d: bigint) =>
  qqResult(() => [...QQ.primes_of_bounded_norm_iter(new Rational(n, d))]);
functions.qq___call__ = (n: bigint) => qqResult(() => QQ.__call__(n));

functions.rational_root_constructed = (base: bigint, k: bigint, offset: bigint) => {
  const x = new Rational(base ** k + offset, 3n ** k);
  const exact = x.is_nth_power(k);
  return JSON.stringify({ value: [exact, exact ? String(x.nth_root(k)) : 'none'] });
};

functions.rational_factory_scalar = (mode: bigint) =>
  qqResult(() => {
    const cases = [
      () => new Rational(),
      () => Rational.from(new Integer(7n)),
      () => Rational.from(true),
      () => Rational.from(null),
      () => Rational.from(),
      () => QQ.__call__(null),
      () => QQ.__call__(true),
      () => QQ.__call__([null, 1n]),
      () => QQ.__call__([1n, null]),
      () => QQ.__call__([true, 2n]),
      () => QQ.__call__([false, 2n]),
      () => QQ.__call__([1n, false]),
      () => QQ.__call__([null]),
      () => Rational.from(false),
      () => Rational.from(new Integer(-3n)),
    ];
    return cases[Number(mode)]!();
  });

functions.qq_tuple_string = (codes: bigint[], den: bigint) =>
  qqResult(() => QQ.__call__([String.fromCodePoint(...codes.map(Number)), den]));


functions.rational_eq_float = (n: bigint, d: bigint, bits: bigint) => {
  const bytes = new DataView(new ArrayBuffer(8));
  bytes.setBigUint64(0, bits, false);
  return new Rational(n, d).eq(bytes.getFloat64(0, false));
};
