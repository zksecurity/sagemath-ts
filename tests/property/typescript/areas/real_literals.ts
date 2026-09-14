import {
  mpfr_init2,
  mpfr_set_str,
  mpfr_set_d,
  mpfr_set_z,
  mpfr_get_d,
  mpfr_get_str,
} from '../../../../packages/mpfr-ts/src/index.js';
const TEXTS = [
  '',
  'x',
  ' ',
  ' +.5',
  '1.2junk',
  '1e',
  '1e+',
  '-nan',
  '@NaN@',
  'nan(foo)',
  '@Inf@',
  'Inf',
  '-Infinity',
  '+0',
  '-0',
  '1.',
  ' .1 ',
  '1_0',
  '1.2e3 ',
  '0e999999999999999999999',
  '1e4096',
  '1e-4096',
  '@nan',
  'nan@',
  'infinityx',
  '-@NaN@',
];
function conversionText(kind: bigint, a: bigint, b: bigint): string {
  if (kind === 0n) return String(a);
  if (kind === 1n) return `${a}e${b}`;
  if (kind === 2n) {
    const n = Number(a) / Number(b);
    return Number.isInteger(n) ? n.toFixed(1) : String(n);
  }
  if (kind === 3n) return TEXTS[Number(a)]!;
  if (kind === 4n)
    return (a < 0n ? '-' : '') + '0.' + '0'.repeat(Number(b)) + String(a < 0n ? -a : a);
  return b >= 0n ? String(a * 2n ** b) : `${a * 5n ** (-b)}e${b}`;
}
function doubleBits(x: number): string {
  if (Number.isNaN(x)) return 'NaN';
  const b = new ArrayBuffer(8),
    v = new DataView(b);
  v.setFloat64(0, x, false);
  return v.getBigUint64(0, false).toString(16).padStart(16, '0');
}
export const functions: Record<string, (...args: any[]) => unknown> = {
  mpfr_conversion: (kind: bigint, a: bigint, b: bigint, precision: bigint, digits: bigint) => {
    const text = conversionText(kind, a, b);
    const mantissa = text.split(/[eE]/)[0]!.replace(/^[-0.]+/, '');
    const sigfigs = mantissa.length - Number(mantissa.includes('.'));
    const p = precision
      ? Number(precision)
      : Math.max(53, Math.floor(3.321928094887363 * sigfigs) + 1);
    const x = mpfr_init2(p),
      status = mpfr_set_str(x, text),
      [formatted, exponent] = mpfr_get_str(10, Number(digits), x);
    return JSON.stringify([
      status,
      x.kind,
      x.sign,
      String(x.mantissa),
      x.exponent,
      p,
      doubleBits(mpfr_get_d(x)),
      formatted,
      exponent,
    ]);
  },
};

import {
  RealField,
  RealLiteral,
  create_RealNumber,
  type RealNumber,
} from '../../../../packages/sagemath-ts/src/rings/real_mpfr.js';
const OBJECTS = [
  null,
  true,
  false,
  NaN,
  Infinity,
  -Infinity,
  -0,
  1.0,
  1e-5,
  1e-4,
  1e16,
  1e20,
  1.133759543500045e153,
];
export function realLiteralFrame(x: RealNumber): unknown {
  let rational: unknown;
  try {
    rational = x.exact_rational().map(String);
  } catch (e) {
    rational = [(e as Error).name, (e as Error).message];
  }
  const [s, m, e] = x.sign_mantissa_exponent();
  return [
    x.toString(),
    x.parent().toString(),
    x.precision(),
    rational,
    [s, String(m), String(e)],
    doubleBits(x.toNumber()),
    x.constructor.name,
    x instanceof RealLiteral ? x.literal : null,
    x instanceof RealLiteral ? x.base : null,
  ];
}
functions.real_literal_conversion = (
  kind: bigint,
  a: bigint,
  b: bigint,
  precision: bigint,
  pad: bigint,
  min: bigint,
  action: bigint,
  arg: bigint
) => {
  const value =
    kind === 6n
      ? OBJECTS[Number(a)]
      : kind === 2n
        ? Number(a) / Number(b)
        : conversionText(kind, a, b);
  let x: RealNumber = precision
    ? new RealLiteral(new RealField(Number(precision)), String(value), 10)
    : create_RealNumber(value, { pad: Number(pad), min_prec: Number(min) });
  if (action === 1n) x = x.neg();
  else if (action === 2n) x = (x as RealLiteral).numerical_approx(Number(arg));
  else if (action === 3n) x = (x as RealLiteral).numerical_approx(undefined, Number(arg));
  else if (action === 4n) x = x.abs();
  else if (action === 5n) x = (x as RealLiteral).numerical_approx(Number(arg)).neg();
  return JSON.stringify(realLiteralFrame(x));
};

const SEQUENCES = [
  ['-1', '', 'nan', '-nan', '@nan', '-inf', 'nan(foo)', '1e+', ''],
  ['1.234567890123456789', '1.2junk', '-0', 'nan', '@NaN@', '-0.25', ''],
  ['-1', 'nan', '-nan', '+nan', '-inf', 'nan(payload)', '0', '-NaN'],
];
functions.mpfr_sequence = (kind: bigint, precision: bigint, digits: bigint) => {
  const x = mpfr_init2(Number(precision));
  return JSON.stringify(
    SEQUENCES[Number(kind)]!.map((text) => {
      const status = mpfr_set_str(x, text),
        [formatted, exponent] = mpfr_get_str(10, Number(digits), x);
      return [
        status,
        x.kind,
        x.sign,
        String(x.mantissa),
        x.exponent,
        x.precision,
        doubleBits(mpfr_get_d(x)),
        formatted,
        exponent,
      ];
    })
  );
};

functions.real_literal_format = (
  a: bigint,
  exponent: bigint,
  precision: bigint,
  digits: bigint,
  flags: bigint
) => {
  const x = new RealLiteral(new RealField(Number(precision)), `${a}e${exponent}`);
  return x.str({
    digits: Number(digits),
    truncate: Boolean(flags & 1n),
    skip_zeroes: Boolean(flags & 2n),
    no_sci: [undefined, false, true, 2][Number(flags >> 2n)] as undefined | boolean | 2,
    e: 'E',
  });
};

functions.mpfr_double_conversion = (raw: bigint, precision: bigint, digits: bigint) => {
  const buffer = new DataView(new ArrayBuffer(8));
  buffer.setBigUint64(0, raw, false);
  const x = mpfr_init2(Number(precision)),
    status = mpfr_set_d(x, buffer.getFloat64(0, false));
  const [formatted, exponent] = mpfr_get_str(10, Number(digits), x);
  return JSON.stringify([
    status,
    x.kind,
    x.sign,
    String(x.mantissa),
    x.exponent,
    x.precision,
    doubleBits(mpfr_get_d(x)),
    formatted,
    exponent,
  ]);
};

const mpfrDoubleSequences = [
  [13830554455654793216n, 9221120237041090560n],
  [9223372036854775808n, 18444492273895866368n, 9221120237041090560n],
  [18442240474082181120n, 9221120237041090560n, 0n],
  [4608308318706860032n, 13832806255468478464n, 1n, 9223372036854775809n],
  [9218868437227405311n, 18442240474082181119n, 4503599627370496n],
  [9221120237041090560n, 18442240474082181120n, 9218868437227405312n, 9223372036854775808n],
  [4607182418800017409n, 4607182418800017407n, 0n],
  [9218868437227405313n, 18442240474082181121n, 13830554455654793216n],
];
functions.mpfr_double_sequence = (kind: bigint, precision: bigint, digits: bigint) => {
  const buffer = new DataView(new ArrayBuffer(8)),
    x = mpfr_init2(Number(precision));
  return JSON.stringify(
    mpfrDoubleSequences[Number(kind)]!.map((raw) => {
      buffer.setBigUint64(0, raw, false);
      const status = mpfr_set_d(x, buffer.getFloat64(0, false)),
        [formatted, exponent] = mpfr_get_str(10, Number(digits), x);
      return [
        status,
        x.kind,
        x.sign,
        String(x.mantissa),
        x.exponent,
        x.precision,
        doubleBits(mpfr_get_d(x)),
        formatted,
        exponent,
      ];
    })
  );
};

import { repr53, fixed53 } from '../../../../packages/sagemath-ts/src/types/real_format.js';
functions.real_format53 = (raw: bigint) => {
  const buffer = new DataView(new ArrayBuffer(8));
  buffer.setBigUint64(0, raw, false);
  const value = buffer.getFloat64(0, false);
  return JSON.stringify([repr53(value), fixed53(value)]);
};

functions.real_numeric_constructor = (
  kind: bigint,
  value: bigint,
  precision: bigint,
  action: bigint
) => {
  const bits = new DataView(new ArrayBuffer(8));
  bits.setBigUint64(0, value & ((1n << 64n) - 1n), false);
  const source =
    kind === 0n ? bits.getFloat64(0, false) : kind === 1n ? value : conversionText(3n, value, 0n);
  let x = new RealField(Number(precision)).__call__(source);
  if (action === 1n) x = x.neg();
  else if (action === 2n) x = x.abs();
  return JSON.stringify([realLiteralFrame(x), x.str(), x.str({ digits: 6 })]);
};

functions.mpfr_integer_conversion = (value: bigint, precision: bigint, digits: bigint) => {
  const x = mpfr_init2(Number(precision)),
    status = mpfr_set_z(x, value);
  const [formatted, exponent] = mpfr_get_str(10, Number(digits), x);
  return JSON.stringify([
    status,
    x.kind,
    x.sign,
    String(x.mantissa),
    x.exponent,
    Number(precision),
    doubleBits(mpfr_get_d(x)),
    formatted,
    exponent,
  ]);
};

import {
  mpfr_sgn,
  mpfr_nan_p,
  mpfr_inf_p,
  mpfr_number_p,
  mpfr_integer_p,
  mpfr_cmp_si,
  mpfr_get_z,
  mpfr_rint,
  mpfr_roundeven,
  mpfr_round,
  mpfr_trunc,
  mpfr_ceil,
  mpfr_floor,
} from '../../../../packages/mpfr-ts/src/index.js';
const REAL_OBSERVER_TEXTS = [
  '0',
  '-0',
  'NaN',
  '+infinity',
  '-infinity',
  '1',
  '-1',
  '1.00000000000000000001',
  '-1.00000000000000000001',
  '9007199254740993',
  '1e400',
  '-1e400',
  '1e-400',
  '-1e-400',
  '0.5',
  '-0.5',
  '-NaN',
  '-@NaN@',
];
const REAL_OBSERVER_OPS = [
  'sign',
  'is_NaN',
  'is_positive_infinity',
  'is_negative_infinity',
  'is_infinity',
  'is_integer',
  'is_square',
  'multiplicative_order',
  'floor',
  'ceil',
  'round',
  'trunc',
] as const;
functions.real_observer = (kind: bigint, a: bigint, b: bigint, precision: bigint, op: bigint) => {
  const text = kind === 0n ? REAL_OBSERVER_TEXTS[Number(a)]! : `${a}e${b}`;
  const x = new RealField(Number(precision)).__call__(text);
  const value = x[REAL_OBSERVER_OPS[Number(op)]!]();
  return typeof value === 'boolean'
    ? value
      ? 'True'
      : 'False'
    : value === Infinity
      ? '+Infinity'
      : String(value);
};
functions.mpfr_observer = (
  kind: bigint,
  a: bigint,
  b: bigint,
  precision: bigint,
  target: bigint,
  operation: bigint,
  mode: bigint,
  integer: bigint
) => {
  const text = kind === 0n ? REAL_OBSERVER_TEXTS[Number(a)]! : `${a}e${b}`;
  const x = mpfr_init2(Number(precision));
  mpfr_set_str(x, text);
  if (operation === 0n)
    return JSON.stringify([
      mpfr_nan_p(x),
      mpfr_inf_p(x),
      mpfr_number_p(x),
      mpfr_integer_p(x),
      mpfr_sgn(x),
      mpfr_cmp_si(x, integer),
    ]);
  const modes = ['RNDN', 'RNDZ', 'RNDU', 'RNDD', 'RNDA', 'RNDNA'] as const;
  if (operation === 7n) {
    const [value, status] = mpfr_get_z(x, modes[Number(mode)] as 'RNDN');
    return JSON.stringify([String(value), status]);
  }
  const result = operation === 8n ? x : mpfr_init2(Number(target));
  const status =
    operation === 1n || operation === 8n
      ? mpfr_rint(result, x, modes[Number(mode)]!)
      : [mpfr_roundeven, mpfr_round, mpfr_trunc, mpfr_ceil, mpfr_floor][Number(operation) - 2]!(
          result,
          x
        );
  const [formatted, exponent] = mpfr_get_str(10, 0, result);
  return JSON.stringify([
    status,
    result.kind,
    result.sign,
    String(result.mantissa),
    result.exponent,
    result.precision,
    doubleBits(mpfr_get_d(result)),
    formatted,
    exponent,
  ]);
};

functions.mpfr_rint_extreme = (
  sign: bigint,
  exponent: bigint,
  precision: bigint,
  target: bigint,
  mode: bigint
) => {
  const x = mpfr_init2(Number(precision)),
    result = mpfr_init2(Number(target));
  mpfr_set_z(x, sign * 3n);
  x.exponent = Number(exponent);
  const modes = ['RNDN', 'RNDZ', 'RNDU', 'RNDD', 'RNDA', 'RNDNA'] as const;
  const status = mpfr_rint(result, x, modes[Number(mode)]!);
  return JSON.stringify([
    status,
    result.kind,
    result.sign,
    String(result.mantissa),
    result.exponent,
    result.precision,
  ]);
};

import { mpfr_set, mpfr_frac, mpfr_cmp } from '../../../../packages/mpfr-ts/src/index.js';
functions.real_fraction = (
  kind: bigint,
  a: bigint,
  b: bigint,
  precision: bigint,
  literal: bigint
) => {
  const text = kind === 0n ? REAL_OBSERVER_TEXTS[Number(a)]! : `${a}e${b}`;
  const R = new RealField(Number(precision));
  const x = literal ? new RealLiteral(R, text) : R.__call__(text);
  return JSON.stringify(realLiteralFrame(x.frac()));
};
functions.real_compare = (
  kind: bigint,
  a: bigint,
  b: bigint,
  precision: bigint,
  otherKind: bigint,
  otherA: bigint,
  otherB: bigint,
  otherPrecision: bigint,
  flags: bigint
) => {
  const left = kind === 0n ? REAL_OBSERVER_TEXTS[Number(a)]! : `${a}e${b}`;
  const right = otherKind === 0n ? REAL_OBSERVER_TEXTS[Number(otherA)]! : `${otherA}e${otherB}`;
  const R = new RealField(Number(precision)),
    S = new RealField(Number(otherPrecision));
  const x = flags & 1n ? new RealLiteral(R, left) : R.__call__(left);
  const y =
    flags & 4n
      ? Number(right.replace('infinity', 'Infinity'))
      : flags & 2n
        ? new RealLiteral(S, right)
        : S.__call__(right);
  return JSON.stringify([x.cmp(y), x.equals(y)]);
};
functions.mpfr_copy_fraction = (
  kind: bigint,
  a: bigint,
  b: bigint,
  precision: bigint,
  target: bigint,
  alias: bigint
) => {
  const text = kind === 0n ? REAL_OBSERVER_TEXTS[Number(a)]! : `${a}e${b}`;
  const x = mpfr_init2(Number(precision));
  mpfr_set_str(x, text);
  const copy = alias ? x : mpfr_init2(Number(target));
  const status = mpfr_set(copy, x);
  const frame = (value: typeof x, inexact: number) => {
    const [digits, exponent] = mpfr_get_str(10, 0, value);
    return [
      inexact,
      value.kind,
      value.sign,
      String(value.mantissa),
      value.exponent,
      value.precision,
      doubleBits(mpfr_get_d(value)),
      digits,
      exponent,
    ];
  };
  const copied = frame(copy, status),
    comparison = mpfr_cmp(x, copy);
  const fraction = alias ? x : mpfr_init2(Number(target));
  const fractional = frame(fraction, mpfr_frac(fraction, x));
  return JSON.stringify([copied, comparison, fractional]);
};

import { mpfr_add, mpfr_mul } from '../../../../packages/mpfr-ts/src/index.js';
const BINARY_TEXTS = [
  '0',
  '-0',
  'nan',
  '-nan',
  'inf',
  '-inf',
  '1',
  '-1',
  '1.5',
  '-1.5',
  '1.25',
  '-1.25',
  '1.75',
  '-1.75',
  '0.1',
  '-0.1',
  '1e1000',
  '-1e1000',
  '1e-1000',
  '-1e-1000',
  '1.00000000000000000000000000000000000001',
  '-1.00000000000000000000000000000000000001',
];
functions.mpfr_binary = (
  left: bigint,
  right: bigint,
  lp: bigint,
  rp: bigint,
  target: bigint,
  operation: bigint,
  alias: bigint,
  le: bigint,
  re: bigint
) => {
  const x = mpfr_init2(Number(lp)),
    y = mpfr_init2(Number(rp)),
    z = mpfr_init2(Number(target));
  mpfr_set_str(x, BINARY_TEXTS[Number(left)]!);
  mpfr_set_str(y, BINARY_TEXTS[Number(right)]!);
  if (le && x.kind === 'finite') x.exponent = Number(le);
  if (re && y.kind === 'finite') y.exponent = Number(re);
  z.kind = 'zero';
  z.sign = -1;
  const result = alias === 1n ? x : alias === 2n ? y : z;
  const status = operation === 0n ? mpfr_add(result, x, y) : mpfr_mul(result, x, y);
  return JSON.stringify([
    status,
    result.kind,
    result.sign,
    String(result.mantissa),
    result.exponent,
    result.precision,
  ]);
};

functions.mpfr_binary_halfway = (
  precision: bigint,
  extra: bigint,
  parity: bigint,
  sign: bigint,
  otherSign: bigint,
  gap: bigint,
  swapped: bigint,
  alias: bigint
) => {
  const mantissa = ((1n << (precision - 1n)) + parity) * (1n << extra) + (1n << (extra - 1n));
  let x = mpfr_init2(Number(precision + extra)),
    y = mpfr_init2(Number(precision + 3n));
  mpfr_set_str(x, String(sign * mantissa));
  mpfr_set_str(y, String(otherSign));
  x.exponent = Number(precision + extra);
  y.exponent = x.exponent - Number(gap);
  if (swapped) [x, y] = [y, x];
  const result = alias === 1n ? x : alias === 2n ? y : mpfr_init2(Number(precision));
  const status = mpfr_add(result, x, y);
  return JSON.stringify([
    status,
    result.kind,
    result.sign,
    String(result.mantissa),
    result.exponent,
    result.precision,
  ]);
};
