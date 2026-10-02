import { polynomialQuotient } from './_polynomial_quotient.js';
/** PARI basemath/Flx.c polynomial sampling; ascending coefficients.
 * @see Deviation: PARI word random-state adapter
 */
import { random_Fl } from './random.js';
export function random_Flx(length: number, p: bigint): bigint[] {
  if (!Number.isSafeInteger(length) || length < 0)
    throw new RangeError('length must be nonnegative');
  const f = Array.from({ length }, () => random_Fl(p));
  while (f.length && f[f.length - 1] === 0n) f.pop();
  return f;
}

import { type FpX } from './ffinit.js';
import { quotientPowers, evaluatePowers } from './_polynomial_composition.js';

/** PARI Flxq_powers: ascending polynomial arrays with native initial entries.
 * @see Deviation: PARI modular composition input contracts
 */
export function Flxq_powers(x: FpX, n: number, T: FpX, p: bigint): FpX[] {
  if (p < 1n || p >= 1n << 64n) throw new RangeError('modulus must be a positive word integer');
  return quotientPowers(x, n, T, p, undefined, true);
}
/** PARI Flx_FlxqV_eval; powers are the zero-indexed Flxq_powers table.
 * @see Deviation: PARI modular composition input contracts
 */
export function Flx_FlxqV_eval(Q: FpX, powers: FpX[], T: FpX, p: bigint): FpX {
  return evaluatePowers(Q, powers, T, p, true);
}
/** PARI Flx_Flxq_eval: sqrt(degree) baby steps and blocked matrix evaluation.
 * @see Deviation: PARI modular composition input contracts
 */
export function Flx_Flxq_eval(Q: FpX, x: FpX, T: FpX, p: bigint): FpX {
  Q=trimPolynomial(Q);
  if (!Q.length) return [];
  return Flx_FlxqV_eval(Q, Flxq_powers(x, Math.floor(Math.sqrt(Q.length - 1)), T, p), T, p);
}

import {
  PACKED_PRODUCT_BIT_BUDGET,
  splitProduct,
  splitSquare,
  packUnsigned,
  unpackUnsigned,
  trimPolynomial,
} from './_polynomial_packing.js';

/** Native maxbitcoeffpol packs only when it saves a complete 64-bit limb. */
function packingBits(p: bigint, n: number): number {
  const bound = (p - 1n) * (p - 1n) * BigInt(n);
  let bits = bound.toString(2).length;
  if (bits <= 16) {
    if (Math.ceil((n * bits) / 64) === Math.ceil(n / 4)) bits = 16;
  } else if (bits <= 32) {
    if (Math.ceil((n * bits) / 64) === Math.ceil(n / 2)) bits = 32;
  } else {
    const limbs = Math.ceil(bits / 64);
    if (Math.ceil((n * bits) / 64) === n * limbs) bits = limbs * 64;
  }
  return bits;
}
function wordPolynomial(x: readonly bigint[], p: bigint): bigint[] {
  if (p < 1n || p >= 1n << 64n) throw new RangeError('modulus must be a positive word integer');
  if (x.some((c) => c < 0n || c >= p))
    throw new RangeError('word polynomial coefficients must be reduced');
  return trimPolynomial(x);
}
/** PARI Flx_mul with fixed 64-bit GMP basecase/packed-integer cutoffs.
 * @see Deviation: PARI polynomial multiplication adapters
 */
export function Flx_mul(x: readonly bigint[], y: readonly bigint[], p: bigint): bigint[] {
  let a = wordPolynomial(x, p),
    b = wordPolynomial(y, p),
    v = 0;
  let va = 0,
    vb = 0;
  while (a[va] === 0n) va++;
  while (b[vb] === 0n) vb++;
  v = va + vb;
  a = a.slice(va);
  b = b.slice(vb);
  if (a.length < b.length) [a, b] = [b, a];
  if (!b.length) return [];
  let result: bigint[];
  if (b.length >= (p <= 3037000493n ? 30 : 8)) {
    const bits = packingBits(p, b.length);
    if ((a.length + b.length - 1) * bits > PACKED_PRODUCT_BIT_BUDGET)
      result = splitProduct(a, b, (u, v) => Flx_mul(u, v, p), p);
    else
      result = unpackUnsigned(
        packUnsigned(a, bits) * packUnsigned(b, bits),
        bits,
        a.length + b.length - 1
      ).map((c) => c % p);
  } else {
    result = new Array<bigint>(a.length + b.length - 1).fill(0n);
    for (let i = 0; i < a.length; i++)
      for (let j = 0; j < b.length; j++) result[i + j]! += a[i]! * b[j]!;
    result = result.map((c) => c % p);
  }
  return trimPolynomial([...new Array<bigint>(v).fill(0n), ...result]);
}
/** PARI Flx_sqr has its own native basecase/packed-integer cutoffs.
 * @see Deviation: PARI polynomial multiplication adapters
 */
export function Flx_sqr(x: readonly bigint[], p: bigint): bigint[] {
  let a = wordPolynomial(x, p),
    v = 0;
  let va = 0;
  while (a[va] === 0n) va++;
  v = 2 * va;
  a = a.slice(va);
  if (!a.length) return [];
  let result: bigint[];
  if (a.length >= (p <= 3037000493n ? 37 : 14)) {
    const bits = packingBits(p, a.length);
    if ((2 * a.length - 1) * bits > PACKED_PRODUCT_BIT_BUDGET)
      result = splitSquare(a, (u) => Flx_sqr(u, p), p);
    else {
      const z = packUnsigned(a, bits);
      result = unpackUnsigned(z * z, bits, 2 * a.length - 1).map((c) => c % p);
    }
  } else {
    result = new Array<bigint>(2 * a.length - 1).fill(0n);
    for (let i = 0; i < a.length; i++) {
      result[2 * i]! += a[i]! * a[i]!;
      for (let j = 0; j < i; j++) result[i + j]! += 2n * a[i]! * a[j]!;
    }
    result = result.map((c) => c % p);
  }
  return trimPolynomial([...new Array<bigint>(v).fill(0n), ...result]);
}

import { PariError } from './errors.js';
import {
  residue,
  inverseCoefficient,
  divisionBasecase,
  inverseBarrett,
  divisionBarrett,
} from './_polynomial_division.js';

/** PARI reciprocal inverse with fixed 64-bit GMP basecase/Newton cutoffs.
 * @see Deviation: PARI polynomial division adapters
 */
export function Flx_invBarrett(T: bigint[], p: bigint): bigint[] {
  const t = wordPolynomial(T, p);
  return inverseBarrett(t, p, true, t.length < (p <= 3037000493n ? 200 : 22), (a, b) =>
    Flx_mul(a, b, p)
  );
}
function wordDivisionBasecase(
  a: bigint[],
  b: bigint[],
  p: bigint,
  remainder: boolean
): [bigint[], bigint[]] {
  if (!b.length) throw new PariError('impossible inverse in Flx_divrem: Vecsmall([0]).');
  if (b.length === 1) {
    if (remainder) return [[], []];
    const inv = b[0] === 1n ? 1n : inverseCoefficient(b[0]!, p, true);
    return [trimPolynomial(a.map((c) => residue(c * inv, p))), []];
  }
  if (a.length < b.length) return [[], a.slice()];
  return divisionBasecase(a, b, p, true);
}
function wordDivision(
  x: bigint[],
  T: bigint[],
  p: bigint,
  remainder: boolean
): [bigint[], bigint[]] {
  const a = wordPolynomial(x, p),
    b = wordPolynomial(T, p),
    d = a.length - b.length;
  if (remainder && d < 0) return [[], a];
  const small = p <= 3037000493n,
    cutoff = remainder ? (small ? 159 : 89) : small ? 161 : 14;
  if (remainder && !b.length && d + 3 < cutoff) {
    if (p === 1n) throw new RangeError('remainder divisor must be nonzero');
    inverseCoefficient(0n, p, true);
  }
  if (d + 3 < cutoff) return wordDivisionBasecase(a, b, p, remainder);
  return divisionBarrett(
    a,
    b,
    Flx_invBarrett(b, p),
    p,
    (u, v) => Flx_mul(u, v, p),
    (u, v) => wordDivisionBasecase(u, v, p, remainder)
  );
}
/** PARI word polynomial quotient/remainder with native Barrett dispatch.
 * @see Deviation: PARI polynomial division adapters
 */
export function Flx_divrem(x: bigint[], T: bigint[], p: bigint): [bigint[], bigint[]] {
  return wordDivision(x, T, p, false);
}
/** PARI word remainder has its own cutoff and constant-divisor shortcut.
 * @see Deviation: PARI polynomial division adapters
 */
export function Flx_rem(x: bigint[], T: bigint[], p: bigint): bigint[] {
  return wordDivision(x, T, p, true)[1];
}

import { FpX_add, FpX_sub } from './ffinit.js';
import {
  type GcdArithmetic,
  type HalfGcdResult,
  type PolynomialMatrix,
  polynomialGcd,
  polynomialExtendedGcd,
  polynomialHalfGcd,
} from './_polynomial_gcd.js';
function wordGcdArithmetic(p: bigint): GcdArithmetic {
  const small = p <= 3037000493n;
  return {
    add: (a, b) => FpX_add(a, b, p),
    sub: (a, b) => FpX_sub(a, b, p),
    mul: (a, b) => Flx_mul(a, b, p),
    divrem: (a, b) => Flx_divrem(a, b, p),
    rem: (a, b) => Flx_rem(a, b, p),
    halfLimit: small ? 120 : 31,
    gcdLimit: small ? 426 : 937,
    extendedLimit: small ? 17 : 22,
    word: true,
  };
}
/** PARI Flx_gcd returns the unscaled Euclidean result.
 * @see Deviation: PARI polynomial GCD adapters
 */
export function Flx_gcd(x: bigint[], y: bigint[], p: bigint): bigint[] {
  return polynomialGcd(wordPolynomial(x, p), wordPolynomial(y, p), wordGcdArithmetic(p));
}
/** Internal native output-pointer selection, shared by quotient inversion. */
export function _Flx_extgcd(
  x: bigint[],
  y: bigint[],
  p: bigint,
  needU: boolean
): [bigint[], bigint[], bigint[]] {
  return polynomialExtendedGcd(
    wordPolynomial(x, p),
    wordPolynomial(y, p),
    wordGcdArithmetic(p),
    needU
  );
}
/** Native unscaled [gcd,u,v], with u*x+v*y=gcd.
 * @see Deviation: PARI polynomial GCD adapters
 */
export function Flx_extgcd(x: bigint[], y: bigint[], p: bigint): [bigint[], bigint[], bigint[]] {
  return _Flx_extgcd(x, y, p, true);
}
/** Native half-GCD matrix and transformed polynomials; matrix rows are zero-indexed.
 * @see Deviation: PARI polynomial GCD adapters
 */
export function Flx_halfgcd_all(x: bigint[], y: bigint[], p: bigint): HalfGcdResult {
  return polynomialHalfGcd(wordPolynomial(x, p), wordPolynomial(y, p), wordGcdArithmetic(p));
}
/** Native half-GCD transformation matrix.
 * @see Deviation: PARI polynomial GCD adapters
 */
export function Flx_halfgcd(x: bigint[], y: bigint[], p: bigint): PolynomialMatrix {
  return Flx_halfgcd_all(x, y, p)[0];
}

import { polynomialMinimalPolynomial } from './_polynomial_minpoly.js';
/** PARI word-field Shoup minimal polynomial, preserving native random-state use.
 * @see Deviation: PARI minimal-polynomial input contracts
 */
export function Flxq_minpoly(x: FpX, T: FpX, p: bigint): FpX {
  return polynomialMinimalPolynomial(wordPolynomial(x, p), wordPolynomial(T, p), p, true);
}

import { polynomialNormalize } from './_polynomial_normalize.js';
/** Native reduced-word normalization, including the zero-polynomial inverse error.
 * @see Deviation: PARI polynomial normalization adapters
 */
export function Flx_normalize(f: FpX, p: bigint): FpX {
  return polynomialNormalize(f, p, true);
}


import { extensionProjection } from './_extension_projection.js';
/** Native FlxqX_dotproduct; ascending coefficient arrays, without variable metadata.
 * @see Deviation: PARI extension projection adapters
 */
export function FlxqX_dotproduct(x: bigint[][], y: bigint[][], T: bigint[], p: bigint): bigint[] {
  return extensionProjection(1, 1, 0, T, p, x, y) as bigint[];
}

import { wordFrobenius } from './_extension_frobenius.js';
/** Native word Frobenius; p=2 bypasses reciprocal preparation.
 * @see Deviation: PARI extension Frobenius adapters
 */
export function Flx_Frobenius(T: FpX, p: bigint): FpX {
  const modulus = wordPolynomial(T, p);
  if (p === 1n) return [0n, 1n];
  return wordFrobenius(modulus, p);
}


/** Native Flx_deriv; coefficients must be reduced modulo the word prime.
 * @see Deviation: PARI conjugate-count bound and word derivative adapters
 */
export function Flx_deriv(z: bigint[], p: bigint): bigint[] {
  const f = wordPolynomial(z, p);
  const out: bigint[] = [];
  for (let i = 1; i < f.length; i++) out.push((BigInt(i) * f[i]!) % p);
  return trimPolynomial(out);
}

/** Native Flx_is_squarefree, including zero and constant polynomials.
 * @see Deviation: PARI conjugate-count bound and word derivative adapters
 */
export function Flx_is_squarefree(z: bigint[], p: bigint): boolean {
  return Flx_gcd(z, Flx_deriv(z, p), p).length === 1;
}

/** Native Euclidean/half-GCD resultant; word inputs have reduced coefficients.
 * @see Deviation: PARI finite-field norm and square-predicate adapters
 */
export function Flx_resultant(x: bigint[], y: bigint[], p: bigint): bigint {
  return polynomialResultant(wordPolynomial(x, p), wordPolynomial(y, p), wordGcdArithmetic(p), p);
}
/** Native resultant norm, including the leading-modulus correction.
 * @see Deviation: PARI finite-field norm and square-predicate adapters
 */
export function Flxq_norm(x: bigint[], T: bigint[], p: bigint): bigint {
  x = wordPolynomial(x, p);
  T = wordPolynomial(T, p);
  const value = Flx_resultant(T, x, p),
    leading = T[T.length - 1]!;
  return leading === 1n || !x.length
    ? value
    : Fp_div(value, Fp_powu(leading, BigInt(x.length - 1), p), p);
}
/** Native square predicate via norm and the quadratic symbol.
 * @see Deviation: PARI finite-field norm and square-predicate adapters
 */
export function Flxq_issquare(x: bigint[], T: bigint[], p: bigint): boolean {
  x = wordPolynomial(x, p);
  if (!x.length || p === 2n) return true;
  return kronecker(Flxq_norm(x, T, p), p) === 1;
}

import { polynomialResultant } from './_polynomial_gcd.js';
import { Fp_div, kronecker } from './ff.js';
import { Fp_powu } from './arith1.js';

/** PARI Flx.c:3848, with reduced word coefficients.
 * @see Deviation: PARI finite-field trace adapters
 */
export function Flxq_trace(x: bigint[], T: bigint[], p: bigint): bigint {
  T = trimPolynomial(T);
  const n = T.length - 2, ctx = polynomialQuotient(T, p, true);
  const z = ctx.reduce(ctx.multiply(x, Flx_deriv(T, p)));
  return z.length - 1 < n ? 0n : Fp_div(z[n]!, T[n + 1]!, p);
}

/** PARI Flx.c:2849: word product/remainder-tree interpolation.
 * @see Deviation: PARI bivariate polynomial storage
 */
export function Flv_polint(x: bigint[], y: bigint[], p: bigint): bigint[] {
  return _polint_tree(x, y, p, true);
}
import { _polint_tree } from './FpX.js';
/** PARI Flx.c:3861: resultant of T(Y) and X-x(Y).
 * @see Deviation: PARI bivariate polynomial storage
 */
export function Flxq_charpoly(x: bigint[], T: bigint[], p: bigint): bigint[] {
  const Q = x.map((c) => trimPolynomial([c === 0n ? 0n : p - c]));
  Q[0] = [x[0] ? p - x[0] : 0n, 1n];
  return Flx_FlxY_resultant(T, Q, p);
}
import { Flx_FlxY_resultant } from './polarit3.js';

import { wordSquareRoot } from './_modular_sqrt.js';
import { FpX_Fp_mul } from './ffinit.js';
import { polynomialQuotientInverse } from './_polynomial_quotient_power.js';
import { gen_pow_fold, gen_powu_i } from './bb_group.js';
import { brent_kung_optpow } from './RgX.js';
import { F2xq_sqrt } from './F2x.js';

/** Flx.c:4690. Quadratic field x+y*sqrt(D), for odd p and nonsquare D.
 * @see Deviation: PARI native extension square-root adapters
 */
export function Fl2_sqrt_pre(z: readonly [bigint, bigint], D: bigint, p: bigint, _pi: bigint): [bigint, bigint] | null {
  const red = (a: bigint) => ((a % p) + p) % p;
  const halve = (a: bigint) => (a & 1n ? a + p : a) / 2n;
  let q = p - 1n; while (!(q & 1n)) q >>= 1n;
  const generator = Fp_powu(D, q, p);
  const sqrt = (a: bigint) => wordSquareRoot(a, p, generator);
  const [a, b] = z;
  if (b === 0n) return kronecker(a, p) === 1 ? [sqrt(a)!, 0n] : [0n, sqrt(Fp_div(a, D, p))!];
  const s = sqrt(red(a * a - D * b * b));
  if (s === null) return null;
  let as2 = halve(red(a + s));
  if (kronecker(as2, p) === -1) as2 = red(as2 - s);
  const u = sqrt(as2)!;
  return [u, Fp_div(b, 2n * u % p, p)];
}

/** Static Flx.c:3523 sum of iterated automorphism products, i>0. */
function wordSumAutSum(a: bigint[], xi: bigint[], i: number, T: bigint[], p: bigint): bigint[] {
  const ctx = polynomialQuotient(T, p, true), d = T.length - 1;
  const mul = (a: bigint[], b: bigint[]) => ctx.reduce(ctx.multiply(a, b));
  const zeta = Flx_Flxq_eval(a, xi, T, p);
  const xp = ctx.powers(xi, brent_kung_optpow(d - 1, 2 * (BigInt(i).toString(2).replaceAll('0', '').length - 1), 1));
  type Triple = [bigint[], bigint[], bigint[]];
  const square = ([x, z, v]: Triple): Triple => {
    const powers = ctx.powers(x, brent_kung_optpow(d - 1, 3, 1));
    return [ctx.evaluate(x, powers), mul(z, ctx.evaluate(z, powers)), FpX_add(v, mul(z, ctx.evaluate(v, powers)), p)];
  };
  const fused = (v: Triple): Triple => {
    const s = square(v), z = mul(zeta, ctx.evaluate(s[1], xp));
    return [ctx.evaluate(s[0], xp), z, FpX_add(s[2], z, p)];
  };
  const result = gen_pow_fold<Triple>([xi, zeta, zeta], BigInt(i), square, fused);
  return mul(a, FpX_add([1n], result[2], p));
}

/** Flx.c:3715: native quadratic, odd-constant and random trace branches.
 * Inputs are reduced; T is irreducible. pi is a retained machine optimization.
 * @see Deviation: PARI native extension square-root adapters
 */
export function Flxq_sqrt_pre(z: bigint[], T: bigint[], p: bigint, pi: bigint): bigint[] | null {
  const d = T.length - 1;
  if (p === 2n) {
    const pack = (v: bigint[]) => v.reduce((a, b, i) => a | (b << BigInt(i)), 0n);
    const r = F2xq_sqrt(pack(z), pack(T));
    return Array.from({length:r === 0n ? 0 : r.toString(2).length}, (_, i) => (r >> BigInt(i)) & 1n);
  }
  if (d === 2) {
    const [c, b, a] = T as [bigint, bigint, bigint], y = z[1] ?? 0n;
    if (a === 1n && b === 0n) {
      const r = Fl2_sqrt_pre([z[0] ?? 0n, y], c === 0n ? 0n : p - c, p, pi);
      return r === null ? null : trimPolynomial(r);
    }
    const b2 = (b & 1n ? b + p : b) / 2n, t = Fp_div(b2, a, p);
    const red = (a: bigint) => ((a % p) + p) % p;
    const D = red(b2 * b2 - a * c), x = red((z[0] ?? 0n) - y * t);
    const r = Fl2_sqrt_pre([x, y], D, p, pi);
    return r === null ? null : trimPolynomial([red(r[0] + r[1] * t), r[1]]);
  }
  if (z.length <= 1 && d % 2 === 1) {
    const r = wordSquareRoot(z[0] ?? 0n, p);
    return r === null ? null : r === 0n ? [] : [r];
  }
  if (!z.length) return [];
  const ctx = polynomialQuotient(T, p, true);
  const mul = (a: bigint[], b: bigint[]) => ctx.reduce(ctx.multiply(a, b));
  const sqr = (a: bigint[]) => ctx.reduce(Flx_sqr(a, p));
  const xi = Flx_Frobenius(T, p);
  let c: bigint[], b: bigint[], newZ: bigint[];
  do {
    do c = random_Flx(d, p); while (!c.length);
    newZ = mul(z, sqr(c));
    const alpha = gen_powu_i(newZ, p >> 1n, sqr, mul);
    b = FpX_add(wordSumAutSum(alpha, xi, d - 2, T, p), [1n], p);
  } while (!b.length);
  const x = mul(newZ, sqr(b));
  if (x.length > 1) return null;
  const beta = wordSquareRoot(x[0] ?? 0n, p);
  if (beta === null) return null;
  return FpX_Fp_mul(polynomialQuotientInverse(mul(b, c), T, p, true), beta, p);
}
/** Flx.c:3782; BigInt reduction requires no preinverse. */
export function Flxq_sqrt(z: bigint[], T: bigint[], p: bigint): bigint[] | null {
  return Flxq_sqrt_pre(z, T, p, 0n);
}

/** Flx.c:3807: inverse Frobenius by residue-class splitting and a dot product.
 * sqx contains powers 0..p-1 of X^(1/p) modulo T. The native reciprocal pi
 * is accepted as a cache hint; bigint reductions need no reciprocal word.
 * @see Deviation: PARI Teichmuller element lifting adapters
 */
export function Flxq_lroot_fast_pre(
  a: bigint[], sqx: bigint[][], T: bigint[], p: bigint, _pi: bigint
): bigint[] {
  const k = Number(p), parts = Array.from({ length: k }, () => [] as bigint[]);
  for (let i = 0; i < a.length; i++) parts[i % k]!.push(a[i]!);
  let sum: bigint[] = [];
  for (let i = 0; i < k; i++) sum = FpX_add(sum, Flx_mul(trimPolynomial(parts[i]!), sqx[i]!, p), p);
  return Flx_rem(sum, T, p);
}
