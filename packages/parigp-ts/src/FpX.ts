import { polynomialQuotient } from './_polynomial_quotient.js';
/** PARI basemath/FpX.c polynomial sampling; ascending coefficients.
 * @see Deviation: PARI word random-state adapter
 */
import { randomi } from './random.js';
export function random_FpX(length: number, p: bigint): bigint[] {
  if (!Number.isSafeInteger(length) || length < 0)
    throw new RangeError('length must be nonnegative');
  const f = Array.from({ length }, () => randomi(p));
  while (f.length && f[f.length - 1] === 0n) f.pop();
  return f;
}

import { FpX_red, FpX_add, type FpX } from './ffinit.js';
import { Flxq_powers, Flx_Flxq_eval } from './Flx.js';
import { quotientPowers, evaluatePowers } from './_polynomial_composition.js';
import { brent_kung_optpow } from './RgX.js';
import { gen_powu_i } from './bb_group.js';

/** PARI FpXQ_powers delegates word tables longer than three entries to Flx.
 * @see Deviation: PARI quotient power count adapter
 * @see Deviation: PARI signed composition and matrix boundaries
 */
export function FpXQ_powers(x: FpX, n: number, T: FpX, p: bigint): FpX[] {
  const pp = p < 0n ? -p : p;
  if (n > 2 && pp > 0n && pp < 1n << 64n) return Flxq_powers(FpX_red(x, pp), n, FpX_red(T, pp), pp);
  return quotientPowers(x, n, T, p);
}
/** PARI FpX_FpXQV_eval; powers are the zero-indexed FpXQ_powers table.
 * @see Deviation: PARI modular composition input contracts
 * @see Deviation: PARI signed composition and matrix boundaries
 */
export function FpX_FpXQV_eval(Q: FpX, powers: FpX[], T: FpX, p: bigint): FpX {
  return evaluatePowers(Q, powers, T, p, false);
}
/** PARI FpX_FpXQ_eval delegates word primes and uses Brent–Kung for large primes.
 * @see Deviation: PARI modular composition input contracts
 * @see Deviation: PARI signed composition and matrix boundaries
 */
export function FpX_FpXQ_eval(Q: FpX, x: FpX, T: FpX, p: bigint): FpX {
  Q=trimPolynomial(Q);
  if (!Q.length) return [];
  const pp = p < 0n ? -p : p;
  if (pp > 0n && pp < 1n << 64n) return Flx_Flxq_eval(FpX_red(Q, pp), FpX_red(x, pp), FpX_red(T, pp), pp);
  return FpX_FpXQV_eval(Q, FpXQ_powers(x, Math.floor(Math.sqrt(Q.length - 1)), T, p), T, p);
}

/** PARI FpXQ_auttrace: power the (automorphism, additive trace) pair.
 * @see Deviation: PARI modular composition input contracts
 * @see Deviation: PARI signed composition and matrix boundaries
 */
export function FpXQ_auttrace(x: [FpX, FpX], n: bigint, T: FpX, p: bigint): [FpX, FpX] {
  if (n <= 0n || n >= 1n << 64n)
    throw new RangeError('gen_powu_i requires a nonzero unsigned word');
  x = [trimPolynomial(x[0]), trimPolynomial(x[1])];
  T = trimPolynomial(T);
  const ctx = polynomialQuotient(T, p, false);
  const multiply = (a: [FpX, FpX], b: [FpX, FpX]): [FpX, FpX] => {
    const d = brent_kung_optpow(Math.max(b[0].length, b[1].length) - 1, 2, 1),
      v = ctx.powers(a[0], d);
    return [ctx.evaluate(b[0], v), FpX_add(a[1], ctx.evaluate(b[1], v), p)];
  };
  const result = gen_powu_i(x, n, (a) => multiply(a, a), multiply);
  return [result[0].slice(), result[1].slice()];
}

import { Flx_mul, Flx_sqr } from './Flx.js';
import { ZX_mul, ZX_sqr, ZX_deriv } from './ZX.js';

/** Native FpX multiplication delegates to word or integer polynomial products.
 * @see Deviation: PARI polynomial multiplication adapters
 * @see Deviation: PARI signed polynomial word dispatch
 */
export function FpX_mul(a: FpX, b: FpX, p: bigint): FpX {
  const pp = p < 0n ? -p : p;
  if (pp > 0n && pp < 1n << 64n) return Flx_mul(FpX_red(a, pp), FpX_red(b, pp), pp);
  return FpX_red(ZX_mul(a, b), p);
}
/** Native FpX squaring preserves the backend's separate square dispatch.
 * @see Deviation: PARI polynomial multiplication adapters
 * @see Deviation: PARI signed polynomial word dispatch
 */
export function FpX_sqr(a: FpX, p: bigint): FpX {
  const pp = p < 0n ? -p : p;
  if (pp > 0n && pp < 1n << 64n) return Flx_sqr(FpX_red(a, pp), pp);
  return FpX_red(ZX_sqr(a), p);
}

import { PariError } from './errors.js';
import { trimPolynomial } from './_polynomial_packing.js';
import {
  residue,
  inverseCoefficient,
  divisionBasecase,
  inverseBarrett,
  divisionBarrett,
} from './_polynomial_division.js';
import { Flx_divrem, Flx_rem } from './Flx.js';

/** PARI FpX_invBarrett: reciprocal inverse through basecase or Newton products.
 * @see Deviation: PARI polynomial division adapters
 * @see Deviation: PARI signed polynomial word dispatch
 */
export function FpX_invBarrett(T: FpX, p: bigint): FpX {
  const t = trimPolynomial(T);
  return inverseBarrett(t, p, false, t.length + 2 <= 111, (a, b) => FpX_mul(a, b, p));
}
function fieldDivisionBasecase(a: FpX, b: FpX, p: bigint, remainder: boolean): [FpX, FpX] {
  if (!b.length) throw new PariError('impossible inverse in FpX_divrem: 0.');
  if (a.length < b.length) return [[], FpX_red(a, p)];
  if (b.length === 1) {
    if (remainder) return [[], []];
    const inv = b[0] === 1n ? 1n : inverseCoefficient(b[0]!, p, false);
    return [
      FpX_red(
        a.map((c) => c * inv),
        p
      ),
      [],
    ];
  }
  const pp = p < 0n ? -p : p;
  if (pp > 0n && pp < 1n << 64n) {
    const x = FpX_red(a, pp),
      T = FpX_red(b, pp);
    return remainder ? [[], Flx_rem(x, T, pp)] : Flx_divrem(x, T, pp);
  }
  return divisionBasecase(a, b, p, false);
}
function fieldDivision(x: FpX, T: FpX, p: bigint, remainder: boolean): [FpX, FpX] {
  const a = trimPolynomial(x),
    b = trimPolynomial(T),
    d = a.length - b.length;
  if (remainder && d < 0) return [[], FpX_red(a, p)];
  if (d + 3 < (remainder ? 111 : 113)) return fieldDivisionBasecase(a, b, p, remainder);
  const pp = p < 0n ? -p : p;
  if (pp > 0n && pp < 1n << 64n) {
    const aa = FpX_red(a, pp),
      bb = FpX_red(b, pp);
    return remainder ? [[], Flx_rem(aa, bb, pp)] : Flx_divrem(aa, bb, pp);
  }
  return divisionBarrett(
    a,
    b,
    FpX_invBarrett(b, p),
    p,
    (u, v) => FpX_mul(u, v, p),
    (u, v) => fieldDivisionBasecase(u, v, p, remainder)
  );
}
/** PARI FpX_divrem preserves native word conversion and reciprocal dispatch.
 * @see Deviation: PARI polynomial division adapters
 * @see Deviation: PARI signed polynomial word dispatch
 */
export function FpX_divrem(x: FpX, T: FpX, p: bigint): [FpX, FpX] {
  return fieldDivision(x, T, p, false);
}
/** PARI FpX_rem preserves smaller-dividend reduction and constant shortcuts.
 * @see Deviation: PARI polynomial division adapters
 * @see Deviation: PARI signed polynomial word dispatch
 */
export function FpX_rem(x: FpX, T: FpX, p: bigint): FpX {
  return fieldDivision(x, T, p, true)[1];
}

import { FpX_sub } from './ffinit.js';
import { Flx_gcd, _Flx_extgcd, Flx_halfgcd_all } from './Flx.js';
import {
  type GcdArithmetic,
  type HalfGcdResult,
  type PolynomialMatrix,
  polynomialGcd,
  polynomialExtendedGcd,
  polynomialHalfGcd,
} from './_polynomial_gcd.js';
function fieldGcdArithmetic(p: bigint): GcdArithmetic {
  return {
    add: (a, b) => FpX_add(a, b, p),
    sub: (a, b) => FpX_sub(a, b, p),
    mul: (a, b) => FpX_mul(a, b, p),
    divrem: (a, b) => FpX_divrem(a, b, p),
    rem: (a, b) => FpX_rem(a, b, p),
    halfLimit: 45,
    gcdLimit: 194,
    extendedLimit: 19,
    word: false,
  };
}
/** PARI gcd preserves the final Euclidean coefficient scale.
 * @see Deviation: PARI polynomial GCD adapters
 * @see Deviation: PARI signed polynomial word dispatch
 */
export function FpX_gcd(x: FpX, y: FpX, p: bigint): FpX {
  const a = FpX_red(x, p),
    b = FpX_red(y, p);
  const pp = p < 0n ? -p : p;
  return pp > 0n && pp < 1n << 64n ? Flx_gcd(a, b, pp) : polynomialGcd(a, b, fieldGcdArithmetic(p));
}
/** Internal native output-pointer selection, shared by quotient inversion. */
export function _FpX_extgcd(x: FpX, y: FpX, p: bigint, needU: boolean): [FpX, FpX, FpX] {
  const a = FpX_red(x, p),
    b = FpX_red(y, p);
  const pp = p < 0n ? -p : p;
  return pp > 0n && pp < 1n << 64n
    ? _Flx_extgcd(a, b, pp, needU)
    : polynomialExtendedGcd(a, b, fieldGcdArithmetic(p), needU);
}
/** Native unscaled [gcd,u,v], with u*x+v*y=gcd.
 * @see Deviation: PARI polynomial GCD adapters
 * @see Deviation: PARI signed polynomial word dispatch
 */
export function FpX_extgcd(x: FpX, y: FpX, p: bigint): [FpX, FpX, FpX] {
  return _FpX_extgcd(x, y, p, true);
}
/** Native half-GCD matrix and transformed polynomials; matrix rows are zero-indexed.
 * @see Deviation: PARI polynomial GCD adapters
 * @see Deviation: PARI signed polynomial word dispatch
 */
export function FpX_halfgcd_all(x: FpX, y: FpX, p: bigint): HalfGcdResult {
  const pp = p < 0n ? -p : p;
  if (pp > 0n && pp < 1n << 64n) return Flx_halfgcd_all(FpX_red(x, pp), FpX_red(y, pp), pp);
  return polynomialHalfGcd(trimPolynomial(x), trimPolynomial(y), fieldGcdArithmetic(p));
}
/** Native half-GCD transformation matrix.
 * @see Deviation: PARI polynomial GCD adapters
 * @see Deviation: PARI signed polynomial word dispatch
 */
export function FpX_halfgcd(x: FpX, y: FpX, p: bigint): PolynomialMatrix {
  return FpX_halfgcd_all(x, y, p)[0];
}

import { polynomialMinimalPolynomial } from './_polynomial_minpoly.js';
import { Flxq_minpoly } from './Flx.js';
/** PARI Shoup minimal polynomial with native word conversion and RNG consumption.
 * @see Deviation: PARI minimal-polynomial input contracts
 */
export function FpXQ_minpoly(x: FpX, T: FpX, p: bigint): FpX {
  if (p > 0n && p < 1n << 64n) return Flxq_minpoly(FpX_red(x, p), FpX_red(T, p), p);
  return polynomialMinimalPolynomial(trimPolynomial(x), trimPolynomial(T), p, false);
}

import { polynomialNormalize } from './_polynomial_normalize.js';
/** Native field normalization preserves a monic input's unreduced coefficients.
 * @see Deviation: PARI polynomial normalization adapters
 */
export function FpX_normalize(f: FpX, p: bigint): FpX {
  return polynomialNormalize(f, p, false);
}

/** Fp_powu for the gaps (at least two) in native sparse polynomial evaluation.
 * Generic exponent two squares the original integer before reduction.
 * @see Deviation: PARI polynomial observation boundaries
 */
function FpX_eval_power(x: bigint, n: number, p: bigint): bigint {
  const pp = p < 0n ? -p : p;
  if (n === 2 && (pp === 0n || pp >= 1n << 64n)) return residue(x * x, p);
  const base = residue(x, p);
  if (pp >= 1n << 64n && base === 1n) return 1n;
  return gen_powu_i(
    base,
    BigInt(n),
    (a) => residue(a * a, p),
    (a, b) => residue(a * b, p)
  );
}

/** Native sparse Horner evaluation, including constant and zero-argument shortcuts.
 * @see Deviation: PARI polynomial observation boundaries
 */
export function FpX_eval(f: FpX, x: bigint, p: bigint): bigint {
  f = trimPolynomial(f);
  if (f.length <= 1 || x === 0n) return f.length ? residue(f[0]!, p) : 0n;
  let value = f[f.length - 1]!;
  for (let i = f.length - 2; i >= 0; ) {
    let j = i;
    while (f[j] === 0n) {
      if (j === 0) {
        const power = i === j ? x : FpX_eval_power(x, i - j + 1, p);
        return residue(value * power, p);
      }
      j--;
    }
    const power = i === j ? x : FpX_eval_power(x, i - j + 1, p);
    value = residue(f[j]! + value * power, p);
    i = j - 1;
  }
  return residue(value, p);
}

/** PARI FpX.c:1193. */
export function FpX_deriv(f: FpX, p: bigint): FpX {
  return FpX_red(ZX_deriv(f), p);
}

/** PARI FpX_div_by_X_x with an omitted remainder output pointer.
 * @see Deviation: PARI polynomial observation boundaries
 */
export function FpX_div_by_X_x(T: FpX, a: bigint, p: bigint): FpX {
  T = trimPolynomial(T);
  const n = T.length - 1;
  if (n <= 0) return [];
  const q: FpX = new Array(n);
  q[n - 1] = T[n]!;
  for (let i = n - 2; i >= 0; i--) q[i] = residue(T[i + 1]! + a * q[i + 1]!, p);
  return q;
}

import { Fp_inv, Fp_mul, Fp_add, Fp_neg } from './ff.js';
import { FpX_Fp_mul } from './ffinit.js';
import { producttree_scheme } from './bb_group.js';
/** Native Montgomery batch inverse (FpX.c:1417), for a nonempty vector. */
export function FpV_inv(x: bigint[], p: bigint): bigint[] {
  if (!x.length) throw new RangeError('FpV_inv requires a nonempty vector');
  const y = [x[0]!];
  for (let i = 1; i < x.length; i++) y[i] = Fp_mul(y[i - 1]!, x[i]!, p);
  let u = Fp_inv(y.at(-1)!, p);
  for (let i = x.length - 1; i > 0; i--) {
    y[i] = Fp_mul(u, y[i - 1]!, p);
    u = Fp_mul(u, x[i]!, p);
  }
  y[0] = u;
  return y;
}
/** Native product-tree leaves have degree one or two, preserving scalar phases. */
function FpV_producttree(x: bigint[], p: bigint): FpX[][] {
  let offset = 0;
  const leaves = producttree_scheme(x.length).map((size) => {
    const a = x[offset++]!;
    if (size === 1) return [Fp_neg(a, p), 1n];
    const b = x[offset++]!;
    return [Fp_mul(a, b, p), Fp_neg(Fp_add(a, b, p), p), 1n];
  });
  const tree = [leaves];
  while (tree.at(-1)!.length > 1) {
    const last = tree.at(-1)!,
      next: FpX[] = [];
    for (let i = 0; i < last.length; i += 2) next.push(FpX_mul(last[i]!, last[i + 1]!, p));
    tree.push(next);
  }
  return tree;
}
/** Remainder tree followed by sparse evaluation at each leaf's original points. */
function FpX_FpV_multieval_tree(P: FpX, x: bigint[], tree: FpX[][], p: bigint): bigint[] {
  let level = [P];
  for (let depth = tree.length - 2; depth >= 0; depth--) {
    const divisors = tree[depth]!,
      next: FpX[] = [];
    for (let i = 0; i < divisors.length; i++) next.push(FpX_rem(level[i >> 1]!, divisors[i]!, p));
    level = next;
  }
  const result: bigint[] = [];
  let offset = 0;
  for (let i = 0; i < tree[0]!.length; i++) {
    const degree = tree[0]![i]!.length - 1;
    for (let k = 0; k < degree; k++) result.push(FpX_eval(level[i]!, x[offset++]!, p));
  }
  return result;
}
/** PARI FpX.c:1865, returning native column-major matrix storage.
 * The points vector is nonempty; galconj retains its existing empty-matrix adapter.
 * @see Deviation: PARI Vandermonde interpolation adapters
 */
export function FpV_invVandermonde(L: bigint[], den: bigint, p: bigint): bigint[][] {
  if (!L.length) throw new RangeError('FpV_invVandermonde requires a nonempty vector');
  const tree = FpV_producttree(L, p),
    T = tree.at(-1)![0]!;
  const inverse = FpV_inv(FpX_FpV_multieval_tree(FpX_deriv(T, p), L, tree, p), p);
  const R = inverse.map((v) => Fp_mul(v, den, p));
  return L.map((a, i) => {
    const P = FpX_Fp_mul(FpX_div_by_X_x(T, a, p), R[i]!, p);
    return Array.from({ length: L.length }, (_, j) => P[j] ?? 0n);
  });
}

import { gen_product } from './bb_group.js';
/** Native balanced modular polynomial product; the empty vector returns integer 1.
 * @see Deviation: PARI bounded factor recombination adapters
 */
export function FpXV_prod(values: bigint[][], p: bigint): bigint[] | bigint {
  const result = gen_product(values, (a, b) => FpX_mul(a, b, p));
  return typeof result === 'bigint' ? result : result.slice();
}

/** Native Euclidean/half-GCD resultant; word inputs have reduced coefficients.
 * @see Deviation: PARI finite-field norm and square-predicate adapters
 */
export function FpX_resultant(x: bigint[], y: bigint[], p: bigint): bigint {
  x = trimPolynomial(x);
  y = trimPolynomial(y);
  if (!x.length || !y.length) return 0n;
  if (p > 0n && p < 1n << 64n) return Flx_resultant(FpX_red(x, p), FpX_red(y, p), p);
  return polynomialResultant(x, y, fieldGcdArithmetic(p), p);
}
/** Native resultant norm, including the leading-modulus correction.
 * @see Deviation: PARI finite-field norm and square-predicate adapters
 */
export function FpXQ_norm(x: bigint[], T: bigint[], p: bigint): bigint {
  x = trimPolynomial(x);
  T = trimPolynomial(T);
  const value = FpX_resultant(T, x, p),
    leading = T[T.length - 1]!;
  return leading === 1n || !x.length
    ? value
    : Fp_div(value, Fp_powu(leading, BigInt(x.length - 1), p), p);
}
/** Native square predicate via norm and the quadratic symbol.
 * @see Deviation: PARI finite-field norm and square-predicate adapters
 */
export function FpXQ_issquare(x: bigint[], T: bigint[], p: bigint): boolean {
  x = trimPolynomial(x);
  T = trimPolynomial(T);
  if (!x.length || p === 2n || p === -2n) return true;
  if (x.length === 1) return (T.length - 1) % 2 === 0 || Fp_issquare(x[0]!, p);
  return kronecker(FpXQ_norm(x, T, p), p) !== -1;
}

import { polynomialResultant } from './_polynomial_gcd.js';
import { Flx_resultant } from './Flx.js';
import { Fp_div, Fp_issquare, kronecker } from './ff.js';
import { Fp_powu } from './arith1.js';

/** PARI FpX.c:3028, derivative product in the quotient ring.
 * Inputs are reduced coefficient vectors and a nonconstant separable modulus.
 * @see Deviation: PARI finite-field trace adapters
 */
export function FpXQ_trace(x: bigint[], T: bigint[], p: bigint): bigint {
  T = trimPolynomial(T);
  const derivative = FpX_deriv(T, p), n = derivative.length - 1;
  const ctx = polynomialQuotient(T, p, false);
  const z = ctx.reduce(ctx.multiply(x, derivative));
  return z.length - 1 < n ? 0n : Fp_div(z[n]!, T[n + 1]!, p);
}
