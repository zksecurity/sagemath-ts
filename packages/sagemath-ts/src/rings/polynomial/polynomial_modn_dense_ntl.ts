/**
 * @module sage/rings/polynomial/polynomial_modn_dense_ntl
 * @description Coppersmith's method for small roots of univariate polynomials
 * over `Zmod(N)` (`small_roots`), following Alexander May's PhD thesis.
 *
 * Sage defines `small_roots` as a module-level function here and binds it as a
 * method of the NTL and FLINT `Zmod(N)[x]` element classes; this port exposes
 * the function and `Polynomial.small_roots`.
 */

import {
  mpfr_ceil,
  mpfr_cmp,
  mpfr_get_z,
  mpfr_init2,
  mpfr_mul,
  mpfr_pow,
  mpfr_set_d,
  mpfr_set_z,
  type mpfr_t,
} from '@sagemath-ts/mpfr-ts';
import { ArithmeticError, ValueError, ZeroDivisionError } from '../../errors.js';
import { IntegerMatrix, LLL } from '../../matrix/matrix_integer.js';
import { type IntegerLike, toBigInt } from '../../types/coercion.js';
import { python_int_hash, python_set } from '../../types/python_set.js';
import { Integer } from '../integer_ring.js';
import type { CoefficientRing, Polynomial, RingElement } from './polynomial_element.js';
import { PolynomialRing } from './polynomial_ring.js';

/** A residue-class coefficient: `Zmod(N)` and prime-field elements. */
export type ModnCoefficient = RingElement & { readonly value: bigint };

/**
 * Keyword arguments of `small_roots`. The LLL options are passed through to
 * `Matrix_integer_dense.LLL` exactly as Sage's `**kwds`.
 */
export interface SmallRootsOptions {
  /** Absolute bound for the root (default: `ceil(1/2 N^(beta^2/delta - epsilon))`). */
  X?: IntegerLike | null;
  /** Find roots modulo a factor `b >= N^beta` of `N` (default 1.0, so `b = N`). */
  beta?: number;
  /** Lattice parameter (default `beta / 8`). */
  epsilon?: number | null;
  delta?: number;
  eta?: number;
  algorithm?: string;
  fp?: string;
  prec?: number;
  early_red?: boolean;
  use_givens?: boolean;
  use_siegel?: boolean;
}

/** `ZZ` as a coefficient ring of `Integer` elements, for `f.change_ring(ZZ)`. */
const ZZ_COEFFICIENTS = {
  zero: () => new Integer(0n),
  one: () => new Integer(1n),
  __call__: (x: unknown): Integer => (x instanceof Integer ? x : new Integer(x as bigint)),
  is_field: () => false,
  is_integral_domain: () => true,
  characteristic: () => 0n,
  toString: () => 'Integer Ring',
} as unknown as CoefficientRing<Integer & RingElement>;

/** An element of Sage's `RR` (53-bit MPFR, nearest-even). */
function RR(): mpfr_t {
  return mpfr_init2(53);
}

function RR_from_integer(n: bigint): mpfr_t {
  const x = RR();
  mpfr_set_z(x, n);
  return x;
}

function RR_from_double(d: number): mpfr_t {
  const x = RR();
  mpfr_set_d(x, d);
  return x;
}

/** `RealNumber.ceil()` (real_mpfr.pyx), including its singular-value error. */
function RR_ceil(x: mpfr_t): bigint {
  if (x.kind === 'nan' || x.kind === 'inf')
    throw new ValueError('Calling ceil() on infinity or NaN');
  const rounded = RR();
  mpfr_ceil(rounded, x);
  return mpfr_get_z(rounded, 'RNDZ')[0];
}

function RR_floor(x: number): bigint {
  if (!Number.isFinite(x)) throw new ValueError('Calling floor() on infinity or NaN');
  return BigInt(Math.floor(x));
}

/** `Integer ** RealNumber`: the integer is coerced into `RR` first. */
function RR_integer_pow(n: bigint, y: number): mpfr_t {
  const z = RR();
  mpfr_pow(z, RR_from_integer(n), RR_from_double(y));
  return z;
}

/** Python's `max(a, b)`: `b` only when `b > a`, so a NaN first argument wins. */
function python_max(a: number, b: number): number {
  return b > a ? b : a;
}

function abs(n: bigint): bigint {
  return n < 0n ? -n : n;
}

function gcd(a: bigint, b: bigint): bigint {
  a = abs(a);
  b = abs(b);
  while (b !== 0n) [a, b] = [b, a % b];
  return a;
}

function polyMul(a: bigint[], b: bigint[]): bigint[] {
  const c: bigint[] = new Array(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++)
    for (let j = 0; j < b.length; j++) c[i + j] = c[i + j]! + a[i]! * b[j]!;
  return c;
}

function characteristic(self: Polynomial<ModnCoefficient>): bigint {
  const ring = self.parent.base_ring as { characteristic?: bigint | (() => bigint | number) };
  const c = ring.characteristic;
  return typeof c === 'function' ? BigInt(c.call(ring)) : BigInt(c!);
}

/**
 * Return the small roots of this polynomial modulo some factor `b` of
 * `N = self.base_ring().characteristic()` with `b >= N^beta`, using
 * Coppersmith's method.
 *
 * Small means `|x| < X`, where `X` is either given or chosen as
 * `ceil(1/2 N^(beta^2/delta - epsilon))` so that the algorithm runs in
 * polynomial time. Roots larger than `X` may also be returned. The real
 * parameters are computed in Sage's 53-bit `RR`, so `X` matches Sage exactly.
 *
 * @param self - A monic polynomial over `Zmod(N)` or a prime field
 * @param options - `X`, `beta`, `epsilon` and the LLL keyword arguments
 * @returns The roots, as elements of the base ring, in the iteration order of
 *   the Python set Sage builds
 * @throws ArithmeticError if the polynomial is not monic
 * @throws ValueError if `beta` is not in `(0, 1]`
 * @see Reference: sage/rings/polynomial/polynomial_modn_dense_ntl.pyx:473 (small_roots)
 * @see Deviation: Coppersmith small roots
 */
export function small_roots<C extends ModnCoefficient>(
  self: Polynomial<C>,
  options: SmallRootsOptions = {}
): C[] {
  const N = characteristic(self as unknown as Polynomial<ModnCoefficient>);

  if (!self.is_monic()) throw new ArithmeticError('Polynomial must be monic.');

  const beta = options.beta ?? 1.0;
  if (beta <= 0.0 || beta > 1.0) throw new ValueError('0.0 < beta <= 1.0 not satisfied.');

  const f = self.coeffs.map((c) => c.value);
  const delta = f.length - 1;
  const epsilon = options.epsilon ?? beta / 8;

  const m = RR_ceil(
    RR_from_double(python_max((beta * beta) / (delta * epsilon), (7 * beta) / delta))
  );
  const t = RR_floor(Number(BigInt(delta) * m) * (1 / beta - 1));

  let X: bigint;
  if (options.X === undefined || options.X === null) {
    X = RR_ceil(
      mpfr_mul_new(RR_from_double(0.5), RR_integer_pow(N, (beta * beta) / delta - epsilon))
    );
  } else {
    X = toBigInt(options.X);
  }

  // g_{i,j} = x^j N^(m-i) f^i for i < m, j < delta; then h_i = x^i f^m for i < t.
  const mm = Number(m);
  const tt = Number(t);
  const fPowers: bigint[][] = [[1n]];
  for (let i = 1; i <= mm; i++) fPowers.push(polyMul(fPowers[i - 1]!, f));
  const g: bigint[][] = [];
  for (let i = 0; i < mm; i++) {
    const scale = N ** (m - BigInt(i));
    for (let j = 0; j < delta; j++)
      g.push([...new Array<bigint>(j).fill(0n), ...fPowers[i]!.map((c) => c * scale)]);
  }
  for (let i = 0; i < tt; i++) g.push([...new Array<bigint>(i).fill(0n), ...fPowers[mm]!]);

  const ncols = delta * mm + Math.max(delta, tt);
  const Xpow: bigint[] = [1n];
  for (let j = 1; j < ncols; j++) Xpow.push(Xpow[j - 1]! * X);
  const rows = g.map((gi) => {
    const row: bigint[] = new Array(ncols).fill(0n);
    for (let j = 0; j < gi.length; j++) row[j] = gi[j]! * Xpow[j]!;
    return row;
  });

  const B = LLL(
    new IntegerMatrix(rows.length, ncols, rows),
    options.delta,
    options.eta,
    options.algorithm,
    options.fp,
    options.prec,
    options.early_red,
    options.use_givens,
    options.use_siegel
  ) as IntegerMatrix;

  const reduced: bigint[] = [];
  for (let i = 0; i < ncols; i++) {
    if (Xpow[i] === 0n) throw new ZeroDivisionError('Integer division by zero');
    reduced.push(floorDiv(B.get(0, i).value, Xpow[i]!));
  }
  const R = new PolynomialRing(ZZ_COEFFICIENTS, 'x')
    .__call__(reduced)
    .roots()
    .map(([r]) => r.value);

  const ZmodN = self.parent.base_ring;
  const roots = python_set(
    R.filter((r) => abs(r) <= X).map((r) => ZmodN.__call__(r)),
    (r) => python_int_hash(r.value),
    (a, b) => a.value === b.value
  );
  const Nbeta = RR_integer_pow(N, beta);
  return roots.filter(
    (root) => mpfr_cmp(RR_from_integer(gcd(N, evaluate(f, root.value, N))), Nbeta) >= 0
  );
}

function mpfr_mul_new(a: mpfr_t, b: mpfr_t): mpfr_t {
  const z = RR();
  mpfr_mul(z, a, b);
  return z;
}

function floorDiv(a: bigint, b: bigint): bigint {
  const q = a / b;
  return a % b !== 0n && a < 0n !== b < 0n ? q - 1n : q;
}

/** `ZZ(self(root))`: the lift of the polynomial's value at `root`. */
function evaluate(f: bigint[], x: bigint, N: bigint): bigint {
  let acc = 0n;
  for (let i = f.length - 1; i >= 0; i--) acc = (acc * x + f[i]!) % N;
  return acc;
}
