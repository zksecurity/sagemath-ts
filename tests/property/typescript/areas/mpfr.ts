/**
 * sagemath-ts side of the `mpfr` property-test area.
 *
 * Covers `sage/rings/real_mpfr` (`RealField` / `RealNumber`) and
 * `sage/rings/complex_mpfr` (`ComplexField` / `ComplexNumber`).
 *
 * Cases: tests/property/cases/mpfr.cases.json
 * SageMath counterpart: tests/property/python/areas/mpfr.py
 *
 * The operand tables and the two formatting regimes are documented in the
 * SageMath counterpart; both files must list the tables in the SAME ORDER,
 * since the cases JSON addresses operands by index.
 */

import {
  ComplexField,
  ComplexNumber,
} from '../../../../packages/sagemath-ts/src/rings/complex_mpfr.js';
import {
  RealField,
  RealNumber,
  RoundingMode,
} from '../../../../packages/sagemath-ts/src/rings/real_mpfr.js';

// ---------------------------------------------------------------------------
// Operand tables.  KEEP IN SYNC (same order!) with python/areas/mpfr.py.
// ---------------------------------------------------------------------------

/** Real operands as decimal strings (see the SageMath counterpart). */
const REAL_VALUES: string[] = [
  '0', // 0
  '-0', // 1   negative zero
  '1', // 2
  '-1', // 3
  '0.5', // 4
  '-0.5', // 5
  '2', // 6
  '-2', // 7
  '3', // 8
  '0.1', // 9   not exactly representable
  '0.3333333333333333', // 10
  '9007199254740992', // 11  2^53
  '9007199254740993', // 12  2^53 + 1  (rounds to 2^53)
  '9007199254740991', // 13  2^53 - 1
  '4503599627370496.5', // 14  2^52 + 1/2, an exact tie
  '0.49999999999999994', // 15  the double just below 1/2
  '-0.49999999999999994', // 16
  '1e-320', // 17  subnormal
  '5e-324', // 18  smallest subnormal
  '1.7976931348623157e308', // 19  largest finite double
  '3.141592653589793', // 20  pi
  '2.718281828459045', // 21  e
  '1e-16', // 22
  '100', // 23
  '1e100', // 24
  '-1e100', // 25
  '0.0001', // 26
  '1e16', // 27
  '1e17', // 28
  '1.5', // 29
  '-1.5', // 30
  '2.5', // 31  tie for round()
  '-2.5', // 32  tie for round()
  '1e-8', // 33
  '1000000', // 34
  '6.02e23', // 35
  'NaN', // 36
  '+infinity', // 37
  '-infinity', // 38
  '-8', // 39
  '0.75', // 40
  '10', // 41
  '-0.75', // 42
  '1e-300', // 43
  '1e300', // 44
  '-3.5', // 45
  '20.333333333333332', // 46  61/3 at 53 bits
  '0.9999999999999999', // 47  the double just below 1
  '1.0000000000000002', // 48  the double just above 1
  '-0.9999999999999999', // 49
  '6', // 50
  '0.999', // 51
  '1e-5', // 52
  '-1e-5', // 53
  '12', // 54
];

/** Complex operands as (real, imaginary) decimal-string pairs. */
const COMPLEX_VALUES: Array<[string, string]> = [
  ['0', '0'], // 0
  ['1', '0'], // 1
  ['-1', '0'], // 2   on the sqrt/log branch cut, from above
  ['-1', '-0'], // 3   on the branch cut, from below (-0.0)
  ['0', '1'], // 4
  ['0', '-1'], // 5
  ['1', '1'], // 6
  ['1', '-1'], // 7
  ['-1', '1'], // 8
  ['-1', '-1'], // 9
  ['-4', '0'], // 10
  ['-4', '-0'], // 11
  ['2', '0'], // 12  arccos/arcsin branch cut, |z| > 1
  ['-2', '0'], // 13
  ['0.5', '0'], // 14
  ['-0.5', '0'], // 15
  ['2', '1'], // 16
  ['-3', '0.5'], // 17  sqrt "avoid_branch" path, im > 0
  ['-3', '-0.5'], // 18  sqrt "avoid_branch" path, im < 0
  ['-3', '4'], // 19  |im| > |re|, re < 0
  ['-3', '-4'], // 20
  ['1e-100', '1e-100'], // 21
  ['1', '1e-100'], // 22
  ['-1', '1e-100'], // 23  just above the branch cut
  ['-1', '-1e-100'], // 24  just below the branch cut
  ['0', '2'], // 25
  ['0', '-2'], // 26
  ['3', '-4'], // 27
  ['1e300', '1e300'], // 28  a^2 + b^2 overflows a double
  ['1e-200', '1e-200'], // 29  a^2 + b^2 underflows a double
  ['0.5', '0.5'], // 30
  ['2', '-0'], // 31
  ['0', '0.5'], // 32
  ['0', '-0.5'], // 33
  ['1.5', '0'], // 34
  ['-1.5', '0'], // 35
  ['0.0001', '-0.0001'], // 36
  ['27', '0'], // 37
  ['-27', '0'], // 38
  ['3.141592653589793', '0'], // 39
  ['0', '3.141592653589793'], // 40
  ['1e-8', '1'], // 41
  ['-1e-8', '1'], // 42
  ['100', '0.001'], // 43
  ['-100', '0.001'], // 44
  ['-100', '-0.001'], // 45
  ['NaN', '1'], // 46
  ['1', 'NaN'], // 47
];

/** Rounding-mode names, indexed to match SageMath's `rnd=` strings. */
const RND_MODES: RoundingMode[] = [
  RoundingMode.RNDN,
  RoundingMode.RNDZ,
  RoundingMode.RNDD,
  RoundingMode.RNDU,
  RoundingMode.RNDA,
];

/** Significant digits for the two formatting regimes. */
const EXACT = 17;
const DISPLAY = 15;

const R = new RealField(53);
const C = new ComplexField(53);

// ---------------------------------------------------------------------------
// Operand construction and formatting
// ---------------------------------------------------------------------------

/**
 * Parse one of the table's decimal strings the way `mpfr_set_str` does.
 *
 * `parseFloat` is correctly rounded and so agrees with MPFR on every finite
 * decimal, but it does not know SageMath's `+infinity` / `-infinity`
 * spellings, so those are mapped explicitly.
 */
function parseReal(s: string): number {
  if (s === 'NaN') return Number.NaN;
  if (s === '+infinity') return Number.POSITIVE_INFINITY;
  if (s === '-infinity') return Number.NEGATIVE_INFINITY;
  return Number.parseFloat(s);
}

function rv(i: bigint): RealNumber {
  return R.__call__(parseReal(REAL_VALUES[Number(i)]!));
}

function cv(i: bigint): ComplexNumber {
  const [re, im] = COMPLEX_VALUES[Number(i)]!;
  return C.__call__(parseReal(re), parseReal(im));
}

/**
 * Format a real value with `digits` significant decimal digits.
 *
 * Signed zeros, NaN and both infinities get their own literal spellings so
 * that a sign flip on a branch cut cannot hide inside "0".
 */
function fmt(x: number, digits: number): string {
  if (Number.isNaN(x)) return 'NaN';
  if (x === Number.POSITIVE_INFINITY) return '+infinity';
  if (x === Number.NEGATIVE_INFINITY) return '-infinity';
  if (x === 0) return Object.is(x, -0) ? '-0' : '0';
  const [mant, exp] = x.toExponential(digits - 1).split('e') as [string, string];
  return `${mant}e${Number(exp)}`;
}

function fmtc(z: ComplexNumber, digits: number): string {
  return `(${fmt(z.real(), digits)}, ${fmt(z.imag(), digits)})`;
}

/** Format a real *or* complex result (SageMath sometimes widens the codomain). */
function fmtv(v: RealNumber | ComplexNumber | number | unknown, digits: number): string {
  if (v instanceof ComplexNumber) return fmtc(v, digits);
  if (typeof v === 'number') return fmt(v, digits);
  if (v instanceof RealNumber) return fmt(v.toNumber(), digits);
  // Mirrors the Python side's `except (TypeError, ValueError): return str(v)`
  // fallback: SageMath sometimes leaves the field entirely (`CC(0).gamma()` is
  // the `UnsignedInfinityRing` generator, printed as `Infinity`).
  return String(v);
}

function fmtlist(zs: ComplexNumber[], digits: number): string {
  return `[${zs.map((z) => fmtc(z, digits)).join(', ')}]`;
}

/** Format a rational the way Python's `str(Rational)` does. */
function fmtrat([num, den]: [bigint, bigint]): string {
  return den === 1n ? `${num}` : `${num}/${den}`;
}

function fmtbool(b: boolean): string {
  return b ? 'True' : 'False';
}

/** Format a multiplicative/additive order (SageMath prints `+Infinity`). */
function fmtorder(n: number): string {
  return n === Number.POSITIVE_INFINITY ? '+Infinity' : String(n);
}

// ---------------------------------------------------------------------------
// Real field / real number
// ---------------------------------------------------------------------------

const rawFunctions: Record<string, (...args: never[]) => string> = {
  real_field_str: (rnd: bigint) => new RealField(53, false, RND_MODES[Number(rnd)]!).toString(),
  real_field_prec: (prec: bigint) => String(new RealField(Number(prec)).precision()),
  real_field_bad_prec: (prec: bigint) => String(new RealField(Number(prec)).precision()),
  real_field_rounding_mode: (rnd: bigint) =>
    RoundingMode[new RealField(53, false, RND_MODES[Number(rnd)]!).rounding_mode()]!,
  real_field_characteristic: () => String(new RealField(53).characteristic()),

  real_rnd_div: (rnd: bigint, i: bigint, j: bigint) => {
    const F = new RealField(53, false, RND_MODES[Number(rnd)]!);
    const a = F.__call__(parseReal(REAL_VALUES[Number(i)]!));
    const b = F.__call__(parseReal(REAL_VALUES[Number(j)]!));
    return fmt(a.div(b).toNumber(), EXACT);
  },
  real_rnd_sqrt: (rnd: bigint, i: bigint) => {
    const F = new RealField(53, false, RND_MODES[Number(rnd)]!);
    return fmtv(F.__call__(parseReal(REAL_VALUES[Number(i)]!)).sqrt(), EXACT);
  },

  real_str: (i: bigint) => rv(i).toString(),
  real_add: (i: bigint, j: bigint) => fmt(rv(i).add(rv(j)).toNumber(), EXACT),
  real_sub: (i: bigint, j: bigint) => fmt(rv(i).sub(rv(j)).toNumber(), EXACT),
  real_mul: (i: bigint, j: bigint) => fmt(rv(i).mul(rv(j)).toNumber(), EXACT),
  real_div: (i: bigint, j: bigint) => fmt(rv(i).div(rv(j)).toNumber(), EXACT),
  real_neg: (i: bigint) => fmt(rv(i).neg().toNumber(), EXACT),
  real_abs: (i: bigint) => fmt(rv(i).abs().toNumber(), EXACT),
  real_sign: (i: bigint) => String(rv(i).sign()),
  real_pow: (i: bigint, j: bigint) => fmtv(rv(i).pow(rv(j)), EXACT),
  real_sqrt: (i: bigint) => fmtv(rv(i).sqrt(), EXACT),
  real_cube_root: (i: bigint) => fmt(rv(i).cube_root().toNumber(), DISPLAY),
  real_nth_root: (i: bigint, n: bigint) => fmt(rv(i).nth_root(Number(n)).toNumber(), DISPLAY),

  real_floor: (i: bigint) => String(rv(i).floor()),
  real_ceil: (i: bigint) => String(rv(i).ceil()),
  real_round: (i: bigint) => String(rv(i).round()),
  real_trunc: (i: bigint) => String(rv(i).trunc()),
  real_frac: (i: bigint) => fmt(rv(i).frac().toNumber(), EXACT),

  real_exact_rational: (i: bigint) => fmtrat(rv(i).exact_rational()),
  real_simplest_rational: (i: bigint) => fmtrat(rv(i).simplest_rational()),
  real_nearby_rational_denom: (i: bigint, d: bigint) => fmtrat(rv(i).nearby_rational(undefined, d)),
  real_nearby_rational_error: (i: bigint, j: bigint) =>
    fmtrat(rv(i).nearby_rational(rv(j).toNumber(), undefined)),
  real_nearby_rational_both: (i: bigint) => fmtrat(rv(i).nearby_rational(0.1, 10n)),
  real_nearby_rational_none: (i: bigint) => fmtrat(rv(i).nearby_rational()),

  real_sign_mantissa_exponent: (i: bigint) => {
    const [s, m, e] = rv(i).sign_mantissa_exponent();
    return `(${s}, ${m}, ${e})`;
  },
  real_fp_rank: (i: bigint) => String(rv(i).fp_rank()),
  real_ulp: (i: bigint) => fmt(rv(i).ulp().toNumber(), EXACT),
  real_epsilon: (i: bigint) => fmt(rv(i).epsilon().toNumber(), EXACT),
  real_nextabove: (i: bigint) => fmt(rv(i).nextabove().toNumber(), EXACT),
  real_nextbelow: (i: bigint) => fmt(rv(i).nextbelow().toNumber(), EXACT),
  real_nexttoward: (i: bigint, j: bigint) => fmt(rv(i).nexttoward(rv(j)).toNumber(), EXACT),

  real_exp: (i: bigint) => fmt(rv(i).exp().toNumber(), DISPLAY),
  real_exp2: (i: bigint) => fmt(rv(i).exp2().toNumber(), DISPLAY),
  real_exp10: (i: bigint) => fmt(rv(i).exp10().toNumber(), DISPLAY),
  real_expm1: (i: bigint) => fmt(rv(i).expm1().toNumber(), DISPLAY),
  real_log: (i: bigint) => fmtv(rv(i).log(), DISPLAY),
  real_log_base: (i: bigint, b: bigint) => fmtv(rv(i).log(Number(b)), DISPLAY),
  real_log2: (i: bigint) => fmtv(rv(i).log2(), DISPLAY),
  real_log10: (i: bigint) => fmtv(rv(i).log10(), DISPLAY),
  real_log1p: (i: bigint) => fmtv(rv(i).log1p(), DISPLAY),

  real_sin: (i: bigint) => fmt(rv(i).sin().toNumber(), DISPLAY),
  real_cos: (i: bigint) => fmt(rv(i).cos().toNumber(), DISPLAY),
  real_tan: (i: bigint) => fmt(rv(i).tan().toNumber(), DISPLAY),
  real_arcsin: (i: bigint) => fmt(rv(i).arcsin().toNumber(), DISPLAY),
  real_arccos: (i: bigint) => fmt(rv(i).arccos().toNumber(), DISPLAY),
  real_arctan: (i: bigint) => fmt(rv(i).arctan().toNumber(), DISPLAY),
  real_sinh: (i: bigint) => fmt(rv(i).sinh().toNumber(), DISPLAY),
  real_cosh: (i: bigint) => fmt(rv(i).cosh().toNumber(), DISPLAY),
  real_tanh: (i: bigint) => fmt(rv(i).tanh().toNumber(), DISPLAY),
  real_arcsinh: (i: bigint) => fmt(rv(i).arcsinh().toNumber(), DISPLAY),
  real_arccosh: (i: bigint) => fmt(rv(i).arccosh().toNumber(), DISPLAY),
  real_arctanh: (i: bigint) => fmt(rv(i).arctanh().toNumber(), DISPLAY),
  real_cot: (i: bigint) => fmt(rv(i).cot().toNumber(), DISPLAY),
  real_sec: (i: bigint) => fmt(rv(i).sec().toNumber(), DISPLAY),
  real_csc: (i: bigint) => fmt(rv(i).csc().toNumber(), DISPLAY),
  real_coth: (i: bigint) => fmt(rv(i).coth().toNumber(), DISPLAY),
  real_sech: (i: bigint) => fmt(rv(i).sech().toNumber(), DISPLAY),
  real_csch: (i: bigint) => fmt(rv(i).csch().toNumber(), DISPLAY),

  real_gamma: (i: bigint) => fmt(rv(i).gamma().toNumber(), DISPLAY),
  real_log_gamma: (i: bigint) => fmtv(rv(i).log_gamma(), DISPLAY),
  real_zeta: (i: bigint) => fmt(rv(i).zeta().toNumber(), DISPLAY),
  real_erf: (i: bigint) => fmt(rv(i).erf().toNumber(), DISPLAY),
  real_erfc: (i: bigint) => fmt(rv(i).erfc().toNumber(), DISPLAY),
  real_eint: (i: bigint) => fmt(rv(i).eint().toNumber(), DISPLAY),
  real_agm: (i: bigint, j: bigint) => fmt(rv(i).agm(rv(j)).toNumber(), DISPLAY),
  real_j0: (i: bigint) => fmt(rv(i).j0().toNumber(), DISPLAY),
  real_j1: (i: bigint) => fmt(rv(i).j1().toNumber(), DISPLAY),
  real_jn: (i: bigint, n: bigint) => fmt(rv(i).jn(Number(n)).toNumber(), DISPLAY),
  real_y0: (i: bigint) => fmt(rv(i).y0().toNumber(), DISPLAY),
  real_y1: (i: bigint) => fmt(rv(i).y1().toNumber(), DISPLAY),
  real_yn: (i: bigint, n: bigint) => fmt(rv(i).yn(Number(n)).toNumber(), DISPLAY),

  real_is_NaN: (i: bigint) => fmtbool(rv(i).is_NaN()),
  real_is_infinity: (i: bigint) => fmtbool(rv(i).is_infinity()),
  real_is_integer: (i: bigint) => fmtbool(rv(i).is_integer()),
  real_is_square: (i: bigint) => fmtbool(rv(i).is_square()),
  real_multiplicative_order: (i: bigint) => fmtorder(rv(i).multiplicative_order()),

  // -------------------------------------------------------------------------
  // Complex field / complex number
  // -------------------------------------------------------------------------

  complex_field_str: (prec: bigint) => new ComplexField(Number(prec)).toString(),
  complex_field_prec: (prec: bigint) => String(new ComplexField(Number(prec)).prec()),
  complex_field_is_exact: () => fmtbool(new ComplexField(53).is_exact()),
  complex_field_characteristic: () => String(new ComplexField(53).characteristic()),
  complex_field_gen: () => fmtc(new ComplexField(53).gen(), EXACT),
  complex_field_ngens: () => String(new ComplexField(53).ngens()),
  complex_field_zeta: (n: bigint) => fmtc(new ComplexField(53).zeta(Number(n)), DISPLAY),

  complex_str: (i: bigint) => cv(i).toString(),
  complex_add: (i: bigint, j: bigint) => fmtc(cv(i).add(cv(j)), EXACT),
  complex_sub: (i: bigint, j: bigint) => fmtc(cv(i).sub(cv(j)), EXACT),
  complex_mul: (i: bigint, j: bigint) => fmtc(cv(i).mul(cv(j)), EXACT),
  complex_div: (i: bigint, j: bigint) => fmtc(cv(i).div(cv(j)), EXACT),
  complex_inv: (i: bigint) => fmtc(cv(i).inv(), EXACT),
  complex_neg: (i: bigint) => fmtc(cv(i).neg(), EXACT),
  complex_conjugate: (i: bigint) => fmtc(cv(i).conjugate(), EXACT),
  complex_abs: (i: bigint) => fmt(cv(i).abs(), EXACT),
  complex_norm: (i: bigint) => fmt(cv(i).norm(), EXACT),
  complex_argument: (i: bigint) => fmt(cv(i).argument(), DISPLAY),

  complex_sqrt: (i: bigint) => fmtc(cv(i).sqrt() as ComplexNumber, EXACT),
  complex_sqrt_all: (i: bigint) => fmtlist(cv(i).sqrt(true) as ComplexNumber[], EXACT),
  complex_nth_root: (i: bigint, n: bigint) =>
    fmtc(cv(i).nth_root(Number(n)) as ComplexNumber, DISPLAY),
  complex_nth_root_all: (i: bigint, n: bigint) =>
    fmtlist(cv(i).nth_root(Number(n), true) as ComplexNumber[], DISPLAY),

  complex_exp: (i: bigint) => fmtc(cv(i).exp(), DISPLAY),
  complex_log: (i: bigint) => fmtc(cv(i).log(), DISPLAY),
  complex_log_base: (i: bigint, b: bigint) => fmtc(cv(i).log(Number(b)), DISPLAY),
  complex_cos: (i: bigint) => fmtc(cv(i).cos(), DISPLAY),
  complex_sin: (i: bigint) => fmtc(cv(i).sin(), DISPLAY),
  complex_tan: (i: bigint) => fmtc(cv(i).tan(), DISPLAY),
  complex_cosh: (i: bigint) => fmtc(cv(i).cosh(), DISPLAY),
  complex_sinh: (i: bigint) => fmtc(cv(i).sinh(), DISPLAY),
  complex_tanh: (i: bigint) => fmtc(cv(i).tanh(), DISPLAY),
  complex_arccos: (i: bigint) => fmtc(cv(i).arccos(), DISPLAY),
  complex_arcsin: (i: bigint) => fmtc(cv(i).arcsin(), DISPLAY),
  complex_arctan: (i: bigint) => fmtc(cv(i).arctan(), DISPLAY),
  complex_arccosh: (i: bigint) => fmtc(cv(i).arccosh(), DISPLAY),
  complex_arcsinh: (i: bigint) => fmtc(cv(i).arcsinh(), DISPLAY),
  complex_arctanh: (i: bigint) => fmtc(cv(i).arctanh(), DISPLAY),
  complex_cot: (i: bigint) => fmtc(cv(i).cot(), DISPLAY),
  complex_sec: (i: bigint) => fmtc(cv(i).sec(), DISPLAY),
  complex_csc: (i: bigint) => fmtc(cv(i).csc(), DISPLAY),
  complex_coth: (i: bigint) => fmtc(cv(i).coth(), DISPLAY),
  complex_sech: (i: bigint) => fmtc(cv(i).sech(), DISPLAY),
  complex_csch: (i: bigint) => fmtc(cv(i).csch(), DISPLAY),

  complex_gamma: (i: bigint) => fmtv(cv(i).gamma(), DISPLAY),
  complex_zeta: (i: bigint) => fmtv(cv(i).zeta(), DISPLAY),
  complex_dilog: (i: bigint) => fmtv(cv(i).dilog(), DISPLAY),
  complex_eta: (i: bigint) => fmtv(cv(i).eta(), DISPLAY),
  complex_agm: (i: bigint, j: bigint) => fmtv(cv(i).agm(cv(j)), DISPLAY),

  complex_is_real: (i: bigint) => fmtbool(cv(i).is_real()),
  complex_is_imaginary: (i: bigint) => fmtbool(cv(i).is_imaginary()),
  complex_is_integer: (i: bigint) => fmtbool(cv(i).is_integer()),
  complex_is_square: (i: bigint) => fmtbool(cv(i).is_square()),
  complex_is_NaN: (i: bigint) => fmtbool(cv(i).is_NaN()),
  complex_multiplicative_order: (i: bigint) => fmtorder(cv(i).multiplicative_order()),
  complex_additive_order: (i: bigint) => fmtorder(cv(i).additive_order()),
  complex_algdep: (i: bigint, n: bigint) =>
    `[${cv(i)
      .algebraic_dependency(Number(n))
      .map((c) => c.toString())
      .join(', ')}]`,
};

/**
 * Render a thrown error as `ERROR: <message>` instead of propagating it.
 *
 * `compare.ts` scores "both sides raised" as a pass without looking at the
 * messages, so a throw has to become a *value* for the harness to check that
 * SageMath and the port reject the same inputs for the same reason.
 */
function guard(f: (...args: never[]) => string): (...args: never[]) => string {
  return (...args: never[]) => {
    try {
      return f(...args);
    } catch (e) {
      return `ERROR: ${e instanceof Error ? e.message : String(e)}`;
    }
  };
}

export const functions: Record<string, (...args: never[]) => string> = Object.fromEntries(
  Object.entries(rawFunctions).map(([name, f]) => [name, guard(f)])
);
