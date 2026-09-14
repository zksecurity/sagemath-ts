import { integer_rational } from '../integer_rational.js';
import { Mod } from '../../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod.js';
/** Integer methods compared directly with sage.rings.integer and its GMP/PARI delegates. */
import {
  Integer,
  IntegerRing,
  ZZ,
} from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';

const unary = [
  'abs',
  'sign',
  'is_prime',
  'is_unit',
  'isqrt',
  'is_square',
  'nbits',
  'prime_divisors',
  'divisors',
  'squarefree_part',
  'next_prime',
  'next_prime_power',
  'is_prime_power',
  'is_perfect_power',
  'is_irreducible',
  'is_pseudoprime',
  'is_squarefree',
  'is_discriminant',
  'is_fundamental_discriminant',
  'factorial',
  'radical',
  'popcount',
  'multiplicative_order',
  'sqrtrem',
  'bits',
  'bit_length',
  'numerator',
  'denominator',
  'floor',
  'ceil',
  'round',
  'trunc',
  '__invert__',
  'previous_prime',
  'class_number',
  'factor',
  'euler_phi',
  'moebius',
  'number_of_divisors',
  'hamming_weight',
  'square',
  'cube',
  'content',
  'primitive_part',
  'continued_fraction',
  'frac',
  'is_power_of_two',
  'is_even',
  'is_odd',
  'bell_number',
  'catalan_number',
  'fibonacci',
  'lucas_number',
  'number_of_partitions',
  'primitive_root',
  'nth_prime',
  'prime_pi',
  'carmichael_lambda',
  'neg',
  'isZero',
  'valueOf',
  'toString',
];
const binary = [
  'gcd',
  'lcm',
  'xgcd',
  'mod',
  'quo_rem',
  'ndigits',
  'valuation',
  'nth_root',
  'exact_log',
  'prime_to_m_part',
  'jacobi',
  'kronecker',
  'binomial',
  'digits',
  '__mod__',
  '__floordiv__',
  'divides',
  'inverse_mod',
  'ord',
  'trial_division',
  'sigma',
  'legendre_symbol',
  'is_coprime',
  'add',
  'sub',
  'mul',
  'div',
  'pow',
  'eq',
  'lt',
  'le',
  'gt',
  'ge',
];

function normalize(value: unknown): unknown {
  if (typeof value === 'boolean') return value;
  if (Array.isArray(value)) return value.map(normalize);
  return String(value);
}
function compare(method: string, n: bigint, args: bigint[]): string {
  try {
    const x = new Integer(n);
    const fn = (x as unknown as Record<string, (...args: bigint[]) => unknown>)[method]!;
    return JSON.stringify({ value: normalize(fn.apply(x, args)) });
  } catch (error) {
    return JSON.stringify({ error: error instanceof Error ? error.name : String(error) });
  }
}
export const functions = Object.fromEntries([
  ...unary.map((method) => [`integer_${method}`, (n: bigint) => compare(method, n, [])]),
  ...binary.map((method) => [
    `integer_${method}`,
    (n: bigint, a: bigint) => compare(method, n, [a]),
  ]),
]);

for (const method of [
  'digit_sum',
  'sqrt_mod',
  'nth_root_mod',
  'is_primitive_root',
  'powermod',
  'log',
  'real_log',
  'global_height',
  'is_quadratic_residue',
  'bit',
  'core',
  'is_strong_pseudoprime',
]) {
  functions[`integer_${method}`] = (n: bigint, ...args: bigint[]) => {
    try {
      const x = new Integer(n);
      const fn = (x as unknown as Record<string, (...args: bigint[]) => unknown>)[method]!;
      let value: unknown;
      if (method === 'nth_root_mod') {
        try {
          const root = fn.apply(x, args) as Integer;
          value =
            root.powermod(args[0]!, args[1]!).value === ((n % args[1]!) + args[1]!) % args[1]!;
          // A bogus root must differ from the oracle even when no roots exist.
          if (!value) value = 'invalid root';
        } catch (error) {
          if (
            error instanceof Error &&
            error.name === 'ValueError' &&
            error.message === 'no n-th root'
          )
            value = false;
          else throw error;
        }
      } else {
        value = fn.apply(x, args);
        if (value === null) value = 'None';
        if (method === 'real_log' || method === 'global_height')
          value = BigInt(Math.round((value as number) * 1e8));
      }
      return JSON.stringify({ value: normalize(value) });
    } catch (error) {
      return JSON.stringify({ error: error instanceof Error ? error.name : String(error) });
    }
  };
}

functions.integer_root_constructed = (base: bigint, k: bigint, offset: bigint, sign: bigint) => {
  try {
    return JSON.stringify({
      value: normalize(new Integer(sign * (base ** k + offset)).nth_root(k, true)),
    });
  } catch (error) {
    return JSON.stringify({ error: error instanceof Error ? error.name : String(error) });
  }
};
functions.integer_root_truncated = (n: bigint, k: bigint) => {
  try {
    return JSON.stringify({ value: normalize(new Integer(n).nth_root(k, true)) });
  } catch (error) {
    return JSON.stringify({ error: error instanceof Error ? error.name : String(error) });
  }
};

function integerConversion(value: unknown, base: bigint, wrapper: bigint): string {
  try {
    const result = wrapper
      ? new (Integer as unknown as new (value: unknown, base: bigint) => Integer)(value, base)
      : (ZZ.__call__ as unknown as (value: unknown, base: bigint) => bigint)(value, base);
    return JSON.stringify({ value: String(result) });
  } catch (error) {
    return JSON.stringify({
      error: error instanceof Error ? error.name : String(error),
      message: error instanceof Error ? error.message : String(error),
    });
  }
}
functions.zz_string = (codes: bigint[], base: bigint, wrapper: bigint) =>
  integerConversion(String.fromCodePoint(...codes.map(Number)), base, wrapper);
functions.zz_number = (n: bigint, d: bigint, wrapper: bigint) =>
  integerConversion(Number(n) / Number(d), 0n, wrapper);
functions.zz_rational = (n: bigint, d: bigint, wrapper: bigint) =>
  integerConversion(new Rational(n, d), 0n, wrapper);
functions.zz_wrapped_integer = (n: bigint, wrapper: bigint) =>
  integerConversion(new Integer(n), 0n, wrapper);
functions.zz_boolean = (n: bigint, wrapper: bigint) => integerConversion(n !== 0n, 0n, wrapper);
functions.zz_empty = (mode: bigint, wrapper: bigint) =>
  integerConversion(mode === 0n ? undefined : null, 0n, wrapper);
functions.zz_digits = (digits: bigint[], base: bigint, wrapper: bigint) =>
  integerConversion(digits, base, wrapper);
for (const method of [
  'getInstance',
  'zero',
  'one',
  'characteristic',
  'is_field',
  'is_ring',
  'is_integral_domain',
  'toString',
]) {
  functions[`zz_${method}`] = () => {
    const value =
      method === 'getInstance'
        ? IntegerRing.getInstance()
        : (ZZ as unknown as Record<string, () => unknown>)[method]!();
    return JSON.stringify({ value: normalize(value) });
  };
}

functions.zz_rational_no_base = (n: bigint, d: bigint, wrapper: bigint) => {
  try {
    const x = new Rational(n, d);
    return JSON.stringify({ value: String(wrapper ? new Integer(x) : ZZ.__call__(x)) });
  } catch (error) {
    return JSON.stringify({
      error: error instanceof Error ? error.name : String(error),
      message: error instanceof Error ? error.message : String(error),
    });
  }
};

functions.zz_large_string = (n: bigint, base: bigint, wrapper: bigint) =>
  integerConversion(n.toString(Number(base)), base, wrapper);
functions.zz_modular = (n: bigint, m: bigint, wrapper: bigint) => {
  try {
    const x = Mod(n, m);
    return JSON.stringify({ value: String(wrapper ? new Integer(x) : ZZ.__call__(x)) });
  } catch (error) {
    return JSON.stringify({
      error: error instanceof Error ? error.name : String(error),
      message: error instanceof Error ? error.message : String(error),
    });
  }
};
functions.zz_hook = (n: bigint, wrapper: bigint) => {
  const x = {
    _integer_: (parent: IntegerRing) => {
      if (parent !== ZZ) throw new Error('wrong parent');
      return new Integer(n);
    },
  };
  return JSON.stringify({ value: String(wrapper ? new Integer(x) : ZZ.__call__(x)) });
};

functions.integer_rational = integer_rational;
