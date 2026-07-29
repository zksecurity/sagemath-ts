/**
 * @module sage/algebras/quatalg/quaternion_algebra
 * @description Quaternion algebras, their orders and their fractional ideals
 *
 * Port of: sage/algebras/quatalg/quaternion_algebra.py
 * Reference: reference/sage/src/sage/algebras/quatalg/quaternion_algebra.py
 *
 * SCOPE OF THIS PORT
 * ------------------
 * Only quaternion algebras over `QQ` are supported -- exactly the case in which
 * SageMath itself supports orders, ideals, maximal orders and ideal arithmetic
 * (`QuaternionOrder`/`QuaternionFractionalIdeal_rational` raise
 * `NotImplementedError` over other base rings).  Constructing an algebra over a
 * finite field, a number field or a function field raises `NotImplementedError`
 * naming what is missing.
 *
 * @see Deviation: Quaternion algebras only over QQ
 * @see Deviation: `quadratic_form()` returns a Gram matrix
 * @see Deviation: `free_module()` returns a `ZZLattice`
 * @see Deviation: infinite places are represented by the string `'infinity'`
 */

import { lllgramint, qf_ZM_apply, qfrep } from '@sagemath-ts/parigp-ts';
import {
  factor,
  gcd,
  hilbert_conductor_inverse,
  hilbert_symbol as hilbert_symbol_ZZ,
  valuation as integer_valuation,
  is_prime,
  kronecker_symbol,
  prime_divisors,
  sqrt_mod,
} from '../../arith/misc.js';
import { NotImplementedError, RuntimeError, ValueError } from '../../errors.js';
import {
  IntegerMatrix,
  hermite_normal_form,
  kernel_matrix,
  left_kernel_matrix,
} from '../../matrix/matrix_integer.js';
import { Integer, ZZ } from '../../rings/integer_ring.js';
import {
  PowerSeriesElement as _PowerSeriesElement,
  PowerSeriesRing as _PowerSeriesRing,
} from '../../rings/power_series_ring.js';
import { Rational } from '../../rings/rational.js';
import { QQ } from '../../rings/rational_field.js';
import {
  integral_matrix_and_denom_from_rational_quaternions,
  rational_matrix_from_rational_quaternions,
  rational_quaternions_from_integral_matrix_and_denom,
} from './quaternion_algebra_cython.js';
import {
  MatrixQQ,
  QuaternionAlgebraElement_rational_field,
  type QuaternionLike,
  type RationalMatrix,
  determinantQQ as determinant,
  inverseQQ as matrix_inverse,
  solveLeftQQ as solve_left,
} from './quaternion_algebra_element.js';

type Quat = QuaternionAlgebraElement_rational_field;

// ===========================================================================
// Small exact-linear-algebra helpers over QQ / ZZ
// ===========================================================================

/**
 * A power series over `ZZ`, structurally typed.
 *
 * As for {@link RationalMatrix}, `Integer` does not satisfy the repo's
 * `RingElement` constraint, so `PowerSeriesRing<Integer>` cannot be written
 * down under `strict`; the runtime objects are genuine `PowerSeriesElement`s.
 *
 * @see Deviation: `Matrix<Rational>` façade
 */
export interface IntegerPowerSeries {
  prec(): number;
  list(): Integer[];
  toString(): string;
}

const PowerSeriesRingZZ = _PowerSeriesRing as unknown as new (
  ring: unknown,
  name?: string,
  default_prec?: number
) => unknown;

const PowerSeriesElementZZ = _PowerSeriesElement as unknown as new (
  parent: unknown,
  coefficients: Integer[],
  prec: number
) => IntegerPowerSeries;

/**
 * `ZZ` as a coefficient ring whose elements are `Integer` objects (the plain
 * `ZZ` of `sage/rings/integer_ring` works with `bigint`s, which do not carry
 * the ring-element methods the power series code calls).
 */
const ZZ_as_coefficient_ring = {
  zero: () => new Integer(0n),
  one: () => new Integer(1n),
  __call__: (x: unknown): Integer => (x instanceof Integer ? x : new Integer(x as bigint)),
  is_field: () => false,
  toString: () => 'Integer Ring',
};

function lcmBigInt(a: bigint, b: bigint): bigint {
  if (a === 0n || b === 0n) return 0n;
  return (a / gcd(a, b)) * b;
}

function absBigInt(a: bigint): bigint {
  return a < 0n ? -a : a;
}

/** Build a `RationalMatrix` from a list of rows. */
function matrix_QQ(rows: Rational[][]): RationalMatrix {
  const m = rows.length;
  const n = m === 0 ? 0 : (rows[0] as Rational[]).length;
  return new MatrixQQ(
    QQ,
    m,
    n,
    rows.map((r) => r.slice())
  );
}

/** The rows of a rational matrix. */
function rows_QQ(M: RationalMatrix): Rational[][] {
  const out: Rational[][] = [];
  for (let i = 0; i < M.nrows; i++) out.push(M.row(i));
  return out;
}

/**
 * Sage's `Matrix_rational_dense._clear_denom`: return `(A, d)` with `A` integral,
 * `d` the lcm of all entry denominators and `A = d * self`.
 */
function clear_denom(rows: Rational[][]): [IntegerMatrix, bigint] {
  let d = 1n;
  for (const row of rows) {
    for (const e of row) d = lcmBigInt(d, e.denominator);
  }
  const m = rows.length;
  const n = m === 0 ? 0 : (rows[0] as Rational[]).length;
  const A = new IntegerMatrix(m, n);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      const e = (rows[i] as Rational[])[j] as Rational;
      A.set(i, j, e.numerator * (d / e.denominator));
    }
  }
  return [A, d];
}

/** Hermite normal form without zero rows (Sage's `_hnf_pari(0, include_zero_rows=False)`). */
function hnf(Z: IntegerMatrix): IntegerMatrix {
  return hermite_normal_form(Z, 'default', undefined, false) as IntegerMatrix;
}

function integer_rows(A: IntegerMatrix): bigint[][] {
  const out: bigint[][] = [];
  for (let i = 0; i < A.nrows; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < A.ncols; j++) row.push(A.get(i, j).value);
    out.push(row);
  }
  return out;
}

function integer_matrix_from_rows(rows: bigint[][]): IntegerMatrix {
  const m = rows.length;
  const n = m === 0 ? 0 : (rows[0] as bigint[]).length;
  const A = new IntegerMatrix(m, n);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) A.set(i, j, (rows[i] as bigint[])[j] as bigint);
  }
  return A;
}

/**
 * A `ZZ`-lattice inside `QQ^n`, given by an echelonized (Hermite) basis.
 *
 * This stands in for `sage.modules.free_module.FreeModule_submodule_with_basis`
 * over `ZZ` embedded in `QQ^n`; `sage.modules` is not part of this port's
 * scope, so `free_module()` returns this object instead.
 *
 * @see Deviation: `free_module()` returns a `ZZLattice`
 */
export class ZZLattice {
  /** Echelon (Hermite) basis rows. */
  readonly rows: Rational[][];
  readonly degree: number;

  private constructor(rows: Rational[][], degree: number) {
    this.rows = rows;
    this.degree = degree;
  }

  /** The `ZZ`-span of the given rational vectors. */
  static span(gens: Rational[][], degree?: number): ZZLattice {
    const deg = degree ?? (gens.length > 0 ? (gens[0] as Rational[]).length : 0);
    if (gens.length === 0) return new ZZLattice([], deg);
    const [A, d] = clear_denom(gens);
    const H = hnf(A);
    const rows: Rational[][] = integer_rows(H).map((r) => r.map((e) => new Rational(e, d)));
    return new ZZLattice(rows, deg);
  }

  /** Rank of the lattice. */
  get rank(): number {
    return this.rows.length;
  }

  /** Basis of the lattice (echelonized). */
  basis(): Rational[][] {
    return this.rows.map((r) => r.slice());
  }

  /** Echelon basis matrix. */
  basis_matrix(): RationalMatrix {
    return matrix_QQ(this.rows);
  }

  /** Equality of lattices. */
  eq(other: ZZLattice): boolean {
    if (this.rank !== other.rank || this.degree !== other.degree) return false;
    for (let i = 0; i < this.rank; i++) {
      for (let j = 0; j < this.degree; j++) {
        if (
          !((this.rows[i] as Rational[])[j] as Rational).eq(
            (other.rows[i] as Rational[])[j] as Rational
          )
        ) {
          return false;
        }
      }
    }
    return true;
  }

  /** Whether every element of this lattice lies in `other`. */
  is_submodule(other: ZZLattice): boolean {
    for (const v of this.rows) {
      if (!other.contains(v)) return false;
    }
    return true;
  }

  /** Whether the vector `v` lies in this lattice. */
  contains(v: Rational[]): boolean {
    if (this.rank === 0) return v.every((e) => e.isZero());
    const B = this.basis_matrix();
    let x: RationalMatrix;
    try {
      x = solve_left(B, matrix_QQ([v]));
    } catch {
      return false;
    }
    // check the solution really is a solution and is integral
    const prod = x.mul(B);
    for (let j = 0; j < this.degree; j++) {
      if (!prod.get(0, j).eq(v[j] as Rational)) return false;
    }
    for (let j = 0; j < x.ncols; j++) {
      if (x.get(0, j).denominator !== 1n) return false;
    }
    return true;
  }

  /** Intersection with another lattice of the same degree. */
  intersection(other: ZZLattice): ZZLattice {
    if (this.rank === 0 || other.rank === 0) return new ZZLattice([], this.degree);
    // Clear denominators of both bases simultaneously
    let d = 1n;
    for (const row of this.rows.concat(other.rows)) {
      for (const e of row) d = lcmBigInt(d, e.denominator);
    }
    const scale = (rows: Rational[][]): bigint[][] =>
      rows.map((r) => r.map((e) => e.numerator * (d / e.denominator)));
    const A = scale(this.rows);
    const B = scale(other.rows);
    const S = integer_matrix_from_rows(A.concat(B));
    const K = left_kernel_matrix(S); // rows (x, y) with x*A + y*B = 0
    const gens: Rational[][] = [];
    for (const k of integer_rows(K)) {
      const v: Rational[] = new Array(this.degree).fill(null).map(() => Rational.zero());
      for (let i = 0; i < A.length; i++) {
        const c = k[i] as bigint;
        if (c === 0n) continue;
        for (let j = 0; j < this.degree; j++) {
          v[j] = (v[j] as Rational).add(new Rational(c * ((A[i] as bigint[])[j] as bigint), d));
        }
      }
      gens.push(v);
    }
    return ZZLattice.span(gens, this.degree);
  }

  /** Sum (as `ZZ`-modules) with another lattice. */
  add(other: ZZLattice): ZZLattice {
    return ZZLattice.span(this.rows.concat(other.rows), this.degree);
  }

  /** Scale the lattice by a rational number. */
  scale(c: Rational | bigint): ZZLattice {
    const r = c instanceof Rational ? c : new Rational(c);
    return ZZLattice.span(
      this.rows.map((row) => row.map((e) => e.mul(r))),
      this.degree
    );
  }

  /**
   * The index `[other : self]` of this lattice in `other`, both of full rank.
   *
   * @see Reference: sage/modules/free_module.py:index_in
   */
  index_in(other: ZZLattice): Rational {
    if (this.rank !== this.degree || other.rank !== other.degree) {
      throw new ValueError('index_in is only implemented for full rank lattices');
    }
    const d1 = determinant(this.basis_matrix());
    const d2 = determinant(other.basis_matrix());
    const q = d1.div(d2);
    return q.abs();
  }

  toString(): string {
    return `Free module of degree ${this.degree} and rank ${this.rank} over Integer Ring`;
  }
}

// ===========================================================================
// Rational Hilbert symbol
// ===========================================================================

/**
 * `hilbert_symbol(a, b, p)` for rational `a`, `b`.
 *
 * Mirrors SageMath's normalisation `a -> a.numerator()*a.denominator()`
 * (`sage/arith/misc.py:4985`) and then delegates to the integral Hilbert
 * symbol of `sage.arith.misc`.
 *
 * @see Reference: sage/arith/misc.py:4922 (hilbert_symbol)
 */
export function hilbert_symbol_QQ(a: Rational, b: Rational, p: bigint): bigint {
  const A = a.numerator * a.denominator;
  const B = b.numerator * b.denominator;
  if (p === -1n) {
    return A < 0n && B < 0n ? -1n : 1n;
  }
  return hilbert_symbol_ZZ(A, B, p);
}

// ===========================================================================
// Constructor
// ===========================================================================

const _quaternion_algebra_cache = new Map<string, QuaternionAlgebra_ab>();

/** The marker used for the (unique) infinite place of `QQ`. */
export const INFINITE_PLACE_QQ = 'infinity';

function normalize_names(names: string | readonly string[] | undefined): [string, string, string] {
  if (names === undefined) return ['i', 'j', 'k'];
  const list = typeof names === 'string' ? names.split(',').map((s) => s.trim()) : names.slice();
  if (list.length !== 3) {
    throw new ValueError('the number of names must equal the number of generators');
  }
  return [list[0] as string, list[1] as string, list[2] as string];
}

/**
 * Construct a quaternion algebra.
 *
 * There are four input formats, of which this port supports those over `QQ`:
 *
 * - `QuaternionAlgebra(a, b)` -- the algebra `(a, b)_QQ`
 * - `QuaternionAlgebra(QQ, a, b)` -- ditto
 * - `QuaternionAlgebra(D)` -- for a squarefree `D >= 1`, the rational
 *   quaternion algebra of discriminant `D`
 * - `QuaternionAlgebra(QQ, primes, inv_archimedean)` -- the rational quaternion
 *   algebra ramified exactly at the given places
 *
 * @example
 * ```typescript
 * QuaternionAlgebra(-2n, -3n).toString()
 * // 'Quaternion Algebra (-2, -3) with base ring Rational Field'
 * QuaternionAlgebra(15n).invariants()  // [-3, 5]
 * ```
 *
 * @see Reference: quaternion_algebra.py:97 (QuaternionAlgebraFactory), :323
 *   (create_key), :441 (create_object)
 */
export function QuaternionAlgebra(
  arg0: unknown,
  arg1?: unknown,
  arg2?: unknown,
  names: string | readonly string[] = 'i,j,k'
): QuaternionAlgebra_ab {
  let a: Rational;
  let b: Rational;

  const isQQ = (x: unknown): boolean => x === QQ;

  if (arg1 === undefined && arg2 === undefined) {
    // QuaternionAlgebra(D)
    const D = toBigIntStrict(arg0);
    const [ai, bi] = hilbert_conductor_inverse(D);
    a = new Rational(ai);
    b = new Rational(bi);
  } else if (arg2 === undefined) {
    // QuaternionAlgebra(a, b)
    a = coerce_to_QQ(arg0);
    b = coerce_to_QQ(arg1);
  } else if (Array.isArray(arg1) && Array.isArray(arg2)) {
    // QuaternionAlgebra(K, primes, inv_archimedean)
    if (!isQQ(arg0)) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: quaternion algebra construction via ramification over a number field'
      );
    }
    const invs = (arg2 as unknown[]).map((x) => coerce_to_QQ(x));
    for (const inv of invs) {
      if (!(inv.isZero() || inv.eq(new Rational(1n, 2n)))) {
        throw new ValueError(
          'list of local invariants specifying ramification should contain only 0 and 1/2'
        );
      }
    }
    const primeSet = new Set<bigint>();
    for (const p of arg1 as unknown[]) {
      const q = toBigIntStrict(p);
      if (!is_prime(q)) {
        throw new ValueError(
          'quaternion algebra constructor requires a list of primes specifying the ramification'
        );
      }
      primeSet.add(q);
    }
    const primes = [...primeSet];
    // Reference: ``if len(arg2) > 1 or (len(arg2) == 1 and
    //   is_odd(len(primes) + 2*arg2[0]))`` (quaternion_algebra.py:379)
    if (invs.length > 1) {
      throw new ValueError(
        'quaternion algebra over the rationals must have an even number of ramified places'
      );
    }
    if (invs.length === 1) {
      const parity = (invs[0] as Rational)
        .mul(new Rational(2n))
        .add(new Rational(BigInt(primes.length)));
      if (parity.denominator !== 1n || parity.numerator % 2n !== 0n) {
        throw new ValueError(
          'quaternion algebra over the rationals must have an even number of ramified places'
        );
      }
    }
    let D = 1n;
    for (const p of primes) D *= p;
    const [ai, bi] = hilbert_conductor_inverse(D);
    a = new Rational(ai);
    b = new Rational(bi);
  } else {
    // QuaternionAlgebra(K, a, b)
    if (!isQQ(arg0)) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: quaternion algebras over base rings other than QQ'
      );
    }
    a = coerce_to_QQ(arg1);
    b = coerce_to_QQ(arg2);
  }

  if (a.isZero() || b.isZero()) {
    throw new ValueError(
      `defining elements of quaternion algebra (${a.toString()}, ${b.toString()}) are not invertible in Rational Field`
    );
  }

  const nm = normalize_names(names);
  const key = `${a.toString()}|${b.toString()}|${nm.join(',')}`;
  const cached = _quaternion_algebra_cache.get(key);
  if (cached) return cached;
  const A = new QuaternionAlgebra_ab(QQ, a, b, nm);
  _quaternion_algebra_cache.set(key, A);
  return A;
}

function coerce_to_QQ(x: unknown): Rational {
  if (x instanceof Rational) return x;
  if (x instanceof Integer) return QQ.__call__(x.value);
  if (typeof x === 'bigint' || typeof x === 'string') return QQ.__call__(x);
  if (typeof x === 'number') {
    if (!Number.isInteger(x)) {
      // SageMath would build the algebra over RR here.
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: quaternion algebras over base rings other than QQ (inexact invariant)'
      );
    }
    return QQ.__call__(BigInt(x));
  }
  if (x === null || x === undefined || typeof x !== 'object') {
    throw new ValueError('a and b must be elements of a ring with characteristic not 2');
  }
  // Anything else is an element of some other ring: SageMath would build the
  // quaternion algebra over that ring, which this port does not support.
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: quaternion algebras over base rings other than QQ'
  );
}

function toBigIntStrict(x: unknown): bigint {
  if (typeof x === 'bigint') return x;
  if (x instanceof Rational) {
    if (x.denominator !== 1n) throw new ValueError('argument must be an integer');
    return x.numerator;
  }
  if (typeof x === 'number' && Number.isInteger(x)) return BigInt(x);
  if (typeof x === 'object' && x !== null && 'value' in (x as Record<string, unknown>)) {
    return (x as { value: bigint }).value;
  }
  throw new ValueError('argument must be an integer');
}

// ===========================================================================
// The algebra
// ===========================================================================

/**
 * A quaternion algebra `(a, b)_QQ`: the `QQ`-algebra generated by `i`, `j`
 * subject to `i^2 = a`, `j^2 = b` and `j i = -i j`.
 *
 * @see Reference: quaternion_algebra.py:784 (QuaternionAlgebra_ab)
 */
export class QuaternionAlgebra_ab {
  readonly _a: Rational;
  readonly _b: Rational;
  private readonly _names: [string, string, string];
  private _gens: [Quat, Quat, Quat] | null = null;
  private _ramified_places_cache: [bigint[], string[]] | null = null;
  private _maximal_order_cache: Map<string, QuaternionOrder> = new Map();

  constructor(base_ring: typeof QQ, a: Rational, b: Rational, names: [string, string, string]) {
    if (base_ring !== QQ) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: quaternion algebras over base rings other than QQ'
      );
    }
    this._a = a;
    this._b = b;
    this._names = names;
  }

  /** The base ring, always `QQ`. */
  base_ring(): typeof QQ {
    return QQ;
  }

  /** The names of the three generators. */
  variable_names(): [string, string, string] {
    return this._names;
  }

  /**
   * Create an element of this algebra.
   *
   * @see Reference: quaternion_algebra.py:806 (__init__, element construction)
   */
  __call__(x: QuaternionLike): Quat {
    return new QuaternionAlgebraElement_rational_field(this, x);
  }

  /** The multiplicative unit. */
  one(): Quat {
    return this.__call__(1n);
  }

  /** The zero element. */
  zero(): Quat {
    return this.__call__(0n);
  }

  /**
   * The number of generators as a `K`-vector space, not including 1: always 3.
   *
   * @see Reference: quaternion_algebra.py:472 (ngens)
   */
  ngens(): number {
    return 3;
  }

  /**
   * The generators `i`, `j`, `k`.
   *
   * @see Reference: quaternion_algebra.py:1283 (gens)
   */
  gens(): [Quat, Quat, Quat] {
    if (this._gens === null) {
      this._gens = [
        this.__call__([0n, 1n, 0n, 0n]),
        this.__call__([0n, 0n, 1n, 0n]),
        this.__call__([0n, 0n, 0n, 1n]),
      ];
    }
    return this._gens;
  }

  /**
   * The `n`-th generator.
   *
   * @see Reference: quaternion_algebra.py:1260 (gen)
   */
  gen(n: number = 0): Quat {
    const g = this.gens();
    if (n < 0 || n > 2) throw new IndexError('n must be between 0 and 2, inclusive');
    return g[n] as Quat;
  }

  /**
   * The fixed basis `1, i, j, k` of this algebra.
   *
   * @see Reference: quaternion_algebra.py:491 (basis)
   */
  basis(): [Quat, Quat, Quat, Quat] {
    const [i, j, k] = this.gens();
    return [this.one(), i, j, k];
  }

  /**
   * The structural invariants `a`, `b`.
   *
   * @see Reference: quaternion_algebra.py:1199 (invariants)
   */
  invariants(): [Rational, Rational] {
    return [this._a, this._b];
  }

  /**
   * The Gram matrix `diag(2, -2a, -2b, 2ab)` of the reduced norm.
   *
   * @see Reference: quaternion_algebra.py:1316 (inner_product_matrix)
   */
  inner_product_matrix(): RationalMatrix {
    const a = this._a;
    const b = this._b;
    const two = new Rational(2n);
    const diag = [two, a.mul(two).neg(), b.mul(two).neg(), a.mul(b).mul(two)];
    const rows: Rational[][] = [];
    for (let i = 0; i < 4; i++) {
      const row: Rational[] = [];
      for (let j = 0; j < 4; j++) row.push(i === j ? (diag[i] as Rational) : Rational.zero());
      rows.push(row);
    }
    return matrix_QQ(rows);
  }

  /**
   * Whether elements are represented exactly: always true over `QQ`.
   *
   * @see Reference: quaternion_algebra.py:621 (is_exact)
   */
  is_exact(): boolean {
    return true;
  }

  /** Always false. @see Reference: quaternion_algebra.py:639 (is_field) */
  is_field(): boolean {
    return false;
  }

  /** Always false over `QQ`. @see Reference: quaternion_algebra.py:652 (is_finite) */
  is_finite(): boolean {
    return false;
  }

  /** Always false. @see Reference: quaternion_algebra.py:670 (is_integral_domain) */
  is_integral_domain(): boolean {
    return false;
  }

  /** Always true. @see Reference: quaternion_algebra.py:683 (is_noetherian) */
  is_noetherian(): boolean {
    return true;
  }

  /** Always false (quaternion algebras are noncommutative). */
  is_commutative(): boolean {
    return false;
  }

  /**
   * The number of elements: `+Infinity` over `QQ`.
   *
   * @see Reference: quaternion_algebra.py:696 (order)
   */
  order(): 'Infinity' {
    return 'Infinity';
  }

  /**
   * Return a random element of this quaternion algebra.
   *
   * @see Reference: quaternion_algebra.py:712 (random_element)
   */
  random_element(numBound?: bigint, denBound?: bigint): Quat {
    const c: Rational[] = [];
    for (let i = 0; i < 4; i++) {
      c.push(numBound === undefined ? QQ.random_element() : QQ.random_element(numBound, denBound));
    }
    return this.__call__(c);
  }

  /**
   * The free module `QQ^4` with the reduced-norm inner product; here just the
   * inner product matrix together with the rank.
   *
   * @see Reference: quaternion_algebra.py:747 (free_module)
   */
  free_module(): { rank: number; inner_product_matrix: RationalMatrix } {
    return { rank: 4, inner_product_matrix: this.inner_product_matrix() };
  }

  /** Alias for {@link free_module}. @see Reference: quaternion_algebra.py:767 */
  vector_space(): { rank: number; inner_product_matrix: RationalMatrix } {
    return this.free_module();
  }

  /**
   * Whether this quaternion algebra is definite, i.e. `a < 0` and `b < 0`.
   *
   * @see Reference: quaternion_algebra.py:1342 (is_definite)
   */
  is_definite(): boolean {
    const [a, b] = this.invariants();
    return a.lt(0n) && b.lt(0n);
  }

  /**
   * Over `QQ`, the same as {@link is_definite}.
   *
   * @see Reference: quaternion_algebra.py:1369 (is_totally_definite)
   */
  is_totally_definite(): boolean {
    return this.is_definite();
  }

  /**
   * The places of `QQ` at which this quaternion algebra ramifies.
   *
   * Returns `[finite, infinite]`; the (unique) infinite place of `QQ` is
   * represented by the string `'infinity'` because ring morphisms are not part
   * of this port.
   *
   * @see Reference: quaternion_algebra.py:1416 (ramified_places)
   * @see Deviation: infinite places are represented by the string `'infinity'`
   */
  ramified_places(inf: boolean = true): [bigint[], string[]] | bigint[] {
    if (this._ramified_places_cache === null) {
      const a = this._a;
      const b = this._b;
      const candidates = new Set<bigint>([2n]);
      for (const x of [a.numerator, a.denominator, b.numerator, b.denominator]) {
        for (const p of prime_divisors(absBigInt(x))) candidates.add(p);
      }
      const ram_fin = [...candidates]
        .filter((p) => hilbert_symbol_QQ(a, b, p) === -1n)
        .sort((x, y) => (x < y ? -1 : x > y ? 1 : 0));
      const ram_inf = this.is_definite() ? [INFINITE_PLACE_QQ] : [];
      this._ramified_places_cache = [ram_fin, ram_inf];
    }
    const [fin, infp] = this._ramified_places_cache;
    if (!inf) return fin.slice();
    return [fin.slice(), infp.slice()];
  }

  /**
   * The finite primes at which this quaternion algebra ramifies.
   *
   * @see Reference: quaternion_algebra.py:1554 (ramified_primes)
   */
  ramified_primes(): bigint[] {
    return this.ramified_places(false) as bigint[];
  }

  /**
   * The discriminant: the product of the ramified finite places.
   *
   * @see Reference: quaternion_algebra.py:1602 (discriminant)
   */
  discriminant(): bigint {
    let d = 1n;
    for (const p of this.ramified_primes()) d *= p;
    return d;
  }

  /**
   * Whether every nonzero element is invertible.
   *
   * @see Reference: quaternion_algebra.py:539 (is_division_algebra)
   */
  is_division_algebra(): boolean {
    const [fin, infp] = this.ramified_places(true) as [bigint[], string[]];
    return !(fin.length === 0 && infp.length === 0);
  }

  /**
   * Whether this algebra is isomorphic to the 2x2 matrix ring over `QQ`.
   *
   * @see Reference: quaternion_algebra.py:580 (is_matrix_ring)
   */
  is_matrix_ring(): boolean {
    const [fin, infp] = this.ramified_places(true) as [bigint[], string[]];
    return fin.length === 0 && infp.length === 0;
  }

  /**
   * Whether this quaternion algebra is isomorphic to `A`.
   *
   * @see Reference: quaternion_algebra.py:1655 (is_isomorphic)
   */
  is_isomorphic(A: QuaternionAlgebra_ab): boolean {
    if (!(A instanceof QuaternionAlgebra_ab)) {
      throw new TypeError('A must be a quaternion algebra of the form (a,b)_K');
    }
    const x = this.ramified_primes();
    const y = A.ramified_primes();
    return x.length === y.length && x.every((p, idx) => p === (y[idx] as bigint));
  }

  /** Equality: same invariants. @see Reference: quaternion_algebra.py:1217 (__eq__) */
  eq(other: unknown): boolean {
    if (!(other instanceof QuaternionAlgebra_ab)) return false;
    return this._a.eq(other._a) && this._b.eq(other._b);
  }

  /**
   * The order of this quaternion algebra with the given basis.
   *
   * @see Reference: quaternion_algebra.py:1738 (quaternion_order)
   */
  quaternion_order(basis: readonly QuaternionLike[], check: boolean = true): QuaternionOrder {
    return new QuaternionOrder(
      this,
      basis.map((x) => this.__call__(x)),
      check
    );
  }

  /**
   * The quaternion ideal with the given generators over `ZZ`.
   *
   * @see Reference: quaternion_algebra.py:1765 (ideal)
   */
  ideal(
    gens: readonly QuaternionLike[],
    options?: {
      left_order?: QuaternionOrder | null;
      right_order?: QuaternionOrder | null;
      check?: boolean;
    }
  ): QuaternionFractionalIdeal_rational {
    const g = gens.map((x) => this.__call__(x));
    return new QuaternionFractionalIdeal_rational(
      this,
      g,
      options?.left_order ?? null,
      options?.right_order ?? null,
      options?.check ?? true
    );
  }

  /**
   * Return a maximal order in this quaternion algebra.
   *
   * The algorithm is Voight's ([Voi2012]); if `take_shortcuts` is true and the
   * discriminant is prime with invariants of a nice form, Proposition 5.2 of
   * [Piz1980] is used instead.
   *
   * @see Reference: quaternion_algebra.py:866 (maximal_order)
   */
  maximal_order(options?: {
    take_shortcuts?: boolean;
    order_basis?: readonly QuaternionLike[] | null;
  }): QuaternionOrder {
    const take_shortcuts = options?.take_shortcuts ?? true;
    const order_basis = options?.order_basis ?? null;

    const cacheKey = `${take_shortcuts}|${
      order_basis === null ? '' : order_basis.map((x) => this.__call__(x).toString()).join(';')
    }`;
    const cached = this._maximal_order_cache.get(cacheKey);
    if (cached) return cached;

    const result = this._maximal_order(take_shortcuts, order_basis);
    this._maximal_order_cache.set(cacheKey, result);
    return result;
  }

  private _maximal_order(
    take_shortcuts: boolean,
    order_basis_in: readonly QuaternionLike[] | null
  ): QuaternionOrder {
    const d_A = this.discriminant();

    // --- Pizer's shortcut ------------------------------------------------
    const [a, b] = this.invariants();
    if (
      !order_basis_in &&
      take_shortcuts &&
      is_prime(d_A) &&
      a.denominator === 1n &&
      b.denominator === 1n
    ) {
      let ai = a.numerator;
      let bi = b.numerator;
      let [i, j, k] = this.gens();

      if (
        (ai !== -1n && bi === -1n) ||
        bi === -2n ||
        (ai !== -1n && ai !== -2n && mod(-ai, 8n) !== 1n)
      ) {
        [ai, bi] = [bi, ai];
        [i, j] = [j, i];
        k = i.mul(j);
      }

      let basis: Quat[] = [];
      const half = new Rational(1n, 2n);
      const quarter = new Rational(1n, 4n);
      if (ai === -1n && bi === -1n) {
        basis = [this.one().add(i).add(j).add(k).scalar_mul(half), i, j, k];
      } else if (ai === -1n && is_prime(-bi) && mod(-bi, 4n) === 3n) {
        basis = [this.one().add(j).scalar_mul(half), i.add(k).scalar_mul(half), j, k];
      } else if (ai === -2n && is_prime(-bi) && mod(-bi, 8n) === 5n) {
        basis = [
          this.one().add(j).add(k).scalar_mul(half),
          i.add(j.scalar_mul(2n)).add(k).scalar_mul(quarter),
          j,
          k,
        ];
      } else if (is_prime(-ai) && is_prime(-bi)) {
        const q = -bi;
        const p = -ai;
        if (mod(q, 4n) === 3n && kronecker_symbol(p, q) === -1n) {
          let t = 0n;
          while (mod(t * t * p + 1n, q) !== 0n) t += 1n;
          basis = [
            this.one().add(j).scalar_mul(half),
            i.add(k).scalar_mul(half),
            j.add(k.scalar_mul(t)).scalar_mul(new Rational(-1n, q)),
            k,
          ];
        }
      }

      if (basis.length > 0) {
        return this.quaternion_order(basis);
      }
    }

    // --- Voight's algorithm ---------------------------------------------
    let order_basis: Quat[];
    if (!order_basis_in) {
      const [i, j, k] = this.gens();
      order_basis = [this.one(), i, j, k];
    } else {
      order_basis = order_basis_in.map((x) => this.__call__(x));
    }

    let R: QuaternionOrder;
    let d_R: Rational;
    try {
      R = this.quaternion_order(order_basis);
      d_R = R.discriminant();
    } catch (e) {
      throw new ValueError(
        'order_basis is not a basis of an order of the given quaternion algebra'
      );
    }

    const basis = basis_for_quaternion_lattice(order_basis);
    const e_new_gens: Quat[] = [];

    if (d_R.denominator !== 1n) {
      throw new RuntimeError('discriminant of an order must be an integer');
    }
    for (const [p] of factor(d_R.numerator)) {
      let e: Quat[] = basis.slice();
      let disc = d_R;
      while (
        integer_valuation(disc.numerator, p) - integer_valuation(disc.denominator, p) >
        integer_valuation(d_A, p)
      ) {
        const f = normalize_basis_at_p(e.slice(), p);

        const A = matrix_QQ(e.map((g) => g.coefficient_tuple()));
        const rhs = matrix_QQ(f.map(([vec]) => vec.coefficient_tuple()));
        const x_rows = rows_QQ(solve_left(A, rhs, false));
        const denoms = x_rows.map((row) => {
          let d = 1n;
          for (const c of row) d = lcmBigInt(d, c.denominator);
          return d;
        });

        const e_n: Quat[] = [];
        for (let idx = 0; idx < 4; idx++) {
          const [vec, val] = f[idx] as [Quat, bigint];
          const v = floorDivBigInt(val, 2n);
          const scale = new Rational(denoms[idx] as bigint, 1n).div(new Rational(p).pow(v));
          e_n.push(vec.scalar_mul(scale));
        }

        // stable sort by valuation mod 2
        const withKey = e_n.map((vec, m) => ({
          vec,
          key: mod((f[m] as [Quat, bigint])[1], 2n),
        }));
        withKey.sort((u, v) => Number(u.key - v.key));
        let en = withKey.map((u) => u.vec);

        if (p !== 2n) {
          if (
            !val_eq(rational_valuation_or_null(constant_of(en[1] as Quat, en[1] as Quat), p), 0n)
          ) {
            if (
              val_eq(rational_valuation_or_null(constant_of(en[2] as Quat, en[2] as Quat), p), 0n)
            ) {
              [en[1], en[2]] = [en[2] as Quat, en[1] as Quat];
            } else {
              [en[1], en[3]] = [en[3] as Quat, en[1] as Quat];
            }
          }

          const aa = constant_of(en[1] as Quat, en[1] as Quat);
          const bb = constant_of(en[2] as Quat, en[2] as Quat);

          if (val_gt(rational_valuation_or_null(bb, p), 0n)) {
            const aInt = to_integer(aa);
            const aMod = mod(aInt, p);
            const root = sqrt_mod(aMod, p);
            if (root !== null) {
              let x = root;
              if (mod(x * x - aInt, p * p) === 0n) x = x + p;
              const g = this.__call__(x)
                .sub(en[1] as Quat)
                .mul(en[2] as Quat)
                .scalar_mul(new Rational(1n, p));
              en[2] = g;
              en[3] = (en[1] as Quat).mul(g);
            }
          }
        } else {
          const t = (en[1] as Quat).reduced_trace();
          const aa = (en[1] as Quat).reduced_norm().neg();
          const bb = constant_of(en[2] as Quat, en[2] as Quat);

          if (val_eq(rational_valuation_or_null(t, p), 0n)) {
            if (val_gt(rational_valuation_or_null(bb, p), 0n)) {
              const aInt = to_integer(aa);
              const tInt = to_integer(t);
              let x = aInt;
              if (mod(x * x - tInt * x + aInt, 4n) === 0n) x = x + 2n;
              const g = this.__call__(x)
                .sub(en[1] as Quat)
                .mul(en[2] as Quat)
                .scalar_mul(new Rational(1n, 2n));
              en[2] = g;
              en[3] = (en[1] as Quat).mul(g);
            }
          } else {
            const [y, z, w] = maxord_solve_aux_eq(to_integer(aa), to_integer(bb), 2n);
            const g = this.one()
              .add((en[1] as Quat).scalar_mul(y))
              .add((en[2] as Quat).scalar_mul(z))
              .add((en[1] as Quat).mul(en[2] as Quat).scalar_mul(w))
              .scalar_mul(new Rational(1n, 2n));
            const h = (en[1] as Quat)
              .scalar_mul(new Rational(z).mul(bb))
              .sub((en[2] as Quat).scalar_mul(new Rational(y).mul(aa)));
            const en1 = en.slice();
            en = [en1[0] as Quat, g, h, g.mul(h)];
            const check = Rational.one()
              .sub(aa.mul(new Rational(y * y)))
              .sub(bb.mul(new Rational(z * z)))
              .add(aa.mul(bb).mul(new Rational(w * w)));
            if (val_gt(rational_valuation_or_null(check, 2n), 2n)) {
              en = basis_for_quaternion_lattice(e.concat(en.slice(1)));
            }
          }
        }

        e = en;

        const L = e.map((x) => e.map((y) => x.pair(y)));
        const det = determinant(matrix_QQ(L));
        disc = det.sqrt() as Rational;
      }
      for (let idx = 1; idx < e.length; idx++) e_new_gens.push(e[idx] as Quat);
    }

    const e_new = basis_for_quaternion_lattice(basis.concat(e_new_gens));
    return this.quaternion_order(e_new);
  }

  /**
   * Return an order in this quaternion algebra with the given level.
   *
   * @see Reference: quaternion_algebra.py:1129 (order_with_level)
   */
  order_with_level(level: bigint): QuaternionOrder {
    if (this.ramified_primes().length > 1) {
      throw new NotImplementedError(
        'currently this algorithm only works when the quaternion algebra is only ramified at one finite prime'
      );
    }
    const lev = absBigInt(level);
    const N = this.discriminant();
    const N1 = gcd(lev, N);
    const M1 = lev / N1;

    let O = this.maximal_order();
    const B = O.basis();

    for (const [p, r] of factor(M1)) {
      // Reference: ``a = int(-p) // 2``, i.e. floor division (quaternion_algebra.py:1183)
      const aShift = floorDivBigInt(-p, 2n);
      // Mirrors the reference loop: ``x`` keeps the last candidate even when no
      // vector satisfies the Kronecker condition (quaternion_algebra.py:1184-1191)
      let x: Quat = this.zero();
      for (const v of gf_vectors(p, 4)) {
        let cand = this.zero();
        for (let idx = 0; idx < 4; idx++) {
          cand = cand.add((B[idx] as Quat).scalar_mul((v[idx] as bigint) + aShift));
        }
        x = cand;
        const t = cand.reduced_trace();
        const nrm = cand.reduced_norm();
        const D = t.mul(t).sub(nrm.mul(new Rational(4n)));
        if (kronecker_symbol(to_integer(D), p) === 1n) {
          break;
        }
      }
      // roots of X^2 - tr(x) X + n(x) over GF(p)
      const tr = to_integer(x.reduced_trace());
      const nr = to_integer(x.reduced_norm());
      const roots: bigint[] = [];
      for (let t = 0n; t < p; t++) {
        if (mod(t * t - tr * t + nr, p) === 0n) roots.push(t);
      }
      if (roots.length === 0) {
        // SageMath raises IndexError from ``roots()[0][0]`` here
        throw new IndexError('list index out of range');
      }
      const aRoot = roots[0] as bigint;
      const pr = pow_bigint(p, r);
      const xa = x.sub(this.__call__(aRoot));
      let xar = this.one();
      for (let e = 0n; e < r; e++) xar = xar.mul(xa);
      const I = O._left_ideal_basis([this.__call__(pr), xar]);
      O = O._right_order_from_ideal_basis(I);
    }
    return O;
  }

  /**
   * Return mod `p` splitting data `I, J, K` for this algebra at the unramified
   * odd prime `p`.
   *
   * @see Reference: quaternion_algebra.py:1794 (modp_splitting_data)
   */
  modp_splitting_data(p: bigint): [bigint[][], bigint[][], bigint[][]] {
    if (!is_prime(p)) {
      throw new ValueError(`p (=${p}) must be prime`);
    }
    if (p === 2n) {
      throw new NotImplementedError('p must be odd');
    }
    if (mod(this.discriminant(), p) === 0n) {
      throw new ValueError(`p (=${p}) must be an unramified prime`);
    }

    const [i, j] = this.gens();
    const i2 = to_Fp(i.mul(i).get(0), p);
    const j2 = to_Fp(j.mul(j).get(0), p);

    const I: bigint[][] = [
      [0n, i2],
      [1n, 0n],
    ];
    if (i2 === 0n) {
      throw new NotImplementedError(
        'algorithm for computing local splittings not implemented in general (currently require the first invariant to be coprime to p)'
      );
    }
    const i2inv = mod_inverse(i2, p);
    let a: bigint | null = null;
    let b = 0n;
    for (let bb = 0n; bb < p; bb++) {
      if (bb === 0n) continue;
      const c = mod(j2 + i2inv * bb * bb, p);
      const s = sqrt_mod(c, p);
      if (s !== null) {
        a = mod(-s, p);
        b = bb;
        break;
      }
    }

    if (a === null) {
      // fallback search over all 2x2 matrices, as in the reference implementation
      for (const J of all_2x2_matrices(p)) {
        const K = mat_mul(I, J, p);
        if (mat_eq(mat_mul(J, J, p), scalar_2x2(j2)) && mat_eq(K, mat_neg(mat_mul(J, I, p), p))) {
          return [I, J, K];
        }
      }
      throw new RuntimeError('bug in modp_splitting_data: no splitting found');
    }

    const J: bigint[][] = [
      [a, b],
      [mod((j2 - a * a) * mod_inverse(b, p), p), mod(-a, p)],
    ];
    const K = mat_mul(I, J, p);
    if (!mat_eq(K, mat_neg(mat_mul(J, I, p), p))) {
      throw new RuntimeError("bug in that I,J don't skew commute");
    }
    return [I, J, K];
  }

  /**
   * Return the map from the (`p`-integral) quaternion algebra to `2x2`
   * matrices over `GF(p)`.
   *
   * @see Reference: quaternion_algebra.py:1905 (modp_splitting_map)
   */
  modp_splitting_map(p: bigint): (q: Quat) => bigint[][] {
    const [I, J, K] = this.modp_splitting_data(p);
    return (q: Quat): bigint[][] => {
      const v = q.coefficient_tuple().map((c) => to_Fp(c, p));
      const out: bigint[][] = [
        [v[0] as bigint, 0n],
        [0n, v[0] as bigint],
      ];
      const addScaled = (M: bigint[][], c: bigint): void => {
        for (let r = 0; r < 2; r++) {
          for (let s = 0; s < 2; s++) {
            const cur = (out[r] as bigint[])[s] as bigint;
            const add = c * ((M[r] as bigint[])[s] as bigint);
            (out[r] as bigint[])[s] = mod(cur + add, p);
          }
        }
      };
      addScaled(I, v[1] as bigint);
      addScaled(J, v[2] as bigint);
      addScaled(K, v[3] as bigint);
      return out;
    };
  }

  /** @see Reference: quaternion_algebra.py:1296 (_repr_) */
  toString(): string {
    return `Quaternion Algebra (${this._a.toString()}, ${this._b.toString()}) with base ring Rational Field`;
  }
}

// ===========================================================================
// Orders
// ===========================================================================

/**
 * An order in a rational quaternion algebra.
 *
 * @see Reference: quaternion_algebra.py:1957 (QuaternionOrder)
 */
export class QuaternionOrder {
  private readonly __basis: Quat[];
  private readonly __quaternion_algebra: QuaternionAlgebra_ab;
  private _unit_ideal: QuaternionFractionalIdeal_rational | null = null;
  private _free_module: ZZLattice | null = null;

  /**
   * @param A - a quaternion algebra
   * @param basis - list of 4 integral quaternions in `A`
   * @param check - whether to do type and other consistency checks
   *
   * @see Reference: quaternion_algebra.py:1968 (__init__)
   */
  constructor(A: QuaternionAlgebra_ab, basis: readonly Quat[], check: boolean = true) {
    if (check) {
      if (!Array.isArray(basis)) {
        throw new TypeError('basis must be a list or tuple');
      }
      if (basis.length !== 4) {
        throw new ValueError('basis must have length 4');
      }
      const b = basis.map((x) => A.__call__(x));

      const M = matrix_QQ(b.map((x) => x.coefficient_tuple()));
      if (determinant(M).isZero()) {
        throw new ValueError('basis must have rank 4');
      }

      const one = matrix_QQ([[Rational.one(), Rational.zero(), Rational.zero(), Rational.zero()]]);
      const v = solve_left(M, one, false);
      for (let j = 0; j < 4; j++) {
        if (v.get(0, j).denominator !== 1n) {
          throw new ValueError('lattice must contain 1');
        }
      }

      const M1 = basis_for_quaternion_lattice(b);
      const prods: Quat[] = [];
      for (const x of b) for (const y of b) prods.push(x.mul(y));
      const M2 = basis_for_quaternion_lattice(b.concat(prods));
      if (M1.length !== M2.length || !M1.every((x, idx) => x.eq(M2[idx] as Quat))) {
        throw new ValueError('given lattice must be a ring');
      }
    }
    this.__basis = basis.map((x) => A.__call__(x));
    this.__quaternion_algebra = A;
  }

  /** The base ring of an order is `ZZ`. */
  base_ring(): typeof ZZ {
    return ZZ;
  }

  /**
   * Construct an element of this order, or throw if `x` is not in it.
   *
   * @see Reference: quaternion_algebra.py:2081 (_element_constructor_)
   */
  __call__(x: QuaternionLike): Quat {
    const y = this.quaternion_algebra().__call__(x);
    if (!this.unit_ideal().contains(y)) {
      throw new TypeError(`${y.toString()} does not lie in ${this.toString()}`);
    }
    return y;
  }

  /** The multiplicative unit. @see Reference: quaternion_algebra.py:2112 (one) */
  one(): Quat {
    return this.quaternion_algebra().one();
  }

  /** The zero element. */
  zero(): Quat {
    return this.quaternion_algebra().zero();
  }

  /** Generators of this order. @see Reference: quaternion_algebra.py:2123 (gens) */
  gens(): Quat[] {
    return this.__basis.slice();
  }

  /** Always 4. @see Reference: quaternion_algebra.py:2134 (ngens) */
  ngens(): number {
    return 4;
  }

  /** The `n`-th generator. @see Reference: quaternion_algebra.py:2145 (gen) */
  gen(n: number): Quat {
    const g = this.__basis[n];
    if (g === undefined) throw new IndexError('generator index out of range');
    return g;
  }

  /** The fixed basis of this order. @see Reference: quaternion_algebra.py:2236 (basis) */
  basis(): Quat[] {
    return this.__basis.slice();
  }

  /**
   * The ambient quaternion algebra.
   *
   * @see Reference: quaternion_algebra.py:2247 (quaternion_algebra)
   */
  quaternion_algebra(): QuaternionAlgebra_ab {
    return this.__quaternion_algebra;
  }

  /** @see Reference: quaternion_algebra.py:2258 (_repr_) */
  toString(): string {
    return `Order of ${this.quaternion_algebra().toString()} with basis (${this.basis()
      .map((x) => x.toString())
      .join(', ')})`;
  }

  /**
   * Return a random element of this order.
   *
   * @see Reference: quaternion_algebra.py:2271 (random_element)
   */
  random_element(x?: bigint, y?: bigint): Quat {
    let s = this.quaternion_algebra().zero();
    for (const b of this.basis()) {
      s = s.add(b.scalar_mul(ZZ.random_element(x, y)));
    }
    return s;
  }

  /**
   * Return the intersection of this order with `other`.
   *
   * @see Reference: quaternion_algebra.py:2294 (intersection)
   */
  intersection(other: QuaternionOrder): QuaternionOrder {
    if (!(other instanceof QuaternionOrder)) {
      throw new TypeError('other must be a QuaternionOrder');
    }
    const A = this.quaternion_algebra();
    if (!other.quaternion_algebra().eq(A)) {
      throw new ValueError('self and other must be in the same ambient quaternion algebra');
    }
    const B = this.free_module();
    const C = other.free_module();
    const inter = B.intersection(C);
    return new QuaternionOrder(
      A,
      inter.basis().map((v) => A.__call__(v))
    );
  }

  /**
   * The free `ZZ`-module corresponding to this order.
   *
   * @see Reference: quaternion_algebra.py:2342 (free_module)
   */
  free_module(): ZZLattice {
    if (this._free_module === null) {
      this._free_module = ZZLattice.span(
        this.basis().map((g) => g.coefficient_tuple()),
        4
      );
    }
    return this._free_module;
  }

  /**
   * The discriminant `sqrt(det(Tr(e_i * conj(e_j))))` of this order.
   *
   * @see Reference: quaternion_algebra.py:2366 (discriminant)
   */
  discriminant(): Rational {
    const e = this.basis();
    const L = e.map((x) => e.map((y) => x.pair(y)));
    return determinant(matrix_QQ(L)).sqrt() as Rational;
  }

  /**
   * Whether this order is maximal in the ambient quaternion algebra.
   *
   * @see Reference: quaternion_algebra.py:2390 (is_maximal)
   */
  is_maximal(): boolean {
    return this.discriminant().eq(new Rational(this.quaternion_algebra().discriminant()));
  }

  /**
   * Return a basis for the left ideal of this order with given generators.
   *
   * @see Reference: quaternion_algebra.py:2425 (_left_ideal_basis)
   */
  _left_ideal_basis(gens: readonly Quat[]): Quat[] {
    const prods: Quat[] = [];
    for (const b of this.basis()) for (const g of gens) prods.push(b.mul(g));
    return basis_for_quaternion_lattice(prods);
  }

  /**
   * Given a basis for a left ideal `I`, return the right order of elements `x`
   * with `I x` contained in `I`.
   *
   * @see Reference: quaternion_algebra.py:2445 (_right_order_from_ideal_basis)
   */
  _right_order_from_ideal_basis(basis: readonly Quat[]): QuaternionOrder {
    const B = this.basis();
    const Z = this.quaternion_algebra();

    const I = matrix_QQ(basis.map((f) => f.coefficient_tuple()));
    const psi = basis.map((f) => matrix_QQ(Z.basis().map((x) => f.mul(x).coefficient_tuple())));
    const psi_inv = psi.map((x) => matrix_inverse(x));
    const W = psi_inv.map((x) => I.mul(x));

    let X = ZZLattice.span(rows_QQ(matrix_QQ(B.map((b) => b.coefficient_tuple()))), 4);
    for (const A of W) {
      X = X.intersection(ZZLattice.span(rows_QQ(A), 4));
    }
    const C = X.basis().map((b) => Z.__call__(b));
    return Z.quaternion_order(C);
  }

  /**
   * The left ideal of this order generated by the given generators.
   *
   * @see Reference: quaternion_algebra.py:2492 (left_ideal)
   */
  left_ideal(
    gens: readonly QuaternionLike[] | QuaternionLike,
    check: boolean = true,
    options?: { is_basis?: boolean }
  ): QuaternionFractionalIdeal_rational {
    const A = this.quaternion_algebra();
    const is_basis = options?.is_basis ?? false;
    let basis: Quat[];
    if (is_basis) {
      basis = (gens as readonly QuaternionLike[]).map((x) => A.__call__(x));
    } else {
      const gl = (
        Array.isArray(gens) && !(gens instanceof QuaternionAlgebraElement_rational_field)
          ? (gens as readonly QuaternionLike[])
          : [gens as QuaternionLike]
      ).map((x) => A.__call__(x));
      const prods: Quat[] = [];
      for (const b of this.basis()) for (const g of gl) prods.push(b.mul(g));
      basis = basis_for_quaternion_lattice(prods);
      check = false;
    }
    return new QuaternionFractionalIdeal_rational(A, basis, this, null, check);
  }

  /**
   * The right ideal of this order generated by the given generators.
   *
   * @see Reference: quaternion_algebra.py:2535 (right_ideal)
   */
  right_ideal(
    gens: readonly QuaternionLike[] | QuaternionLike,
    check: boolean = true,
    options?: { is_basis?: boolean }
  ): QuaternionFractionalIdeal_rational {
    const A = this.quaternion_algebra();
    const is_basis = options?.is_basis ?? false;
    let basis: Quat[];
    if (is_basis) {
      basis = (gens as readonly QuaternionLike[]).map((x) => A.__call__(x));
    } else {
      const gl = (
        Array.isArray(gens) && !(gens instanceof QuaternionAlgebraElement_rational_field)
          ? (gens as readonly QuaternionLike[])
          : [gens as QuaternionLike]
      ).map((x) => A.__call__(x));
      const prods: Quat[] = [];
      for (const b of this.basis()) for (const g of gl) prods.push(g.mul(b));
      basis = basis_for_quaternion_lattice(prods);
      check = false;
    }
    return new QuaternionFractionalIdeal_rational(A, basis, null, this, check);
  }

  /**
   * The unit ideal of this order.
   *
   * @see Reference: quaternion_algebra.py:2579 (unit_ideal)
   */
  unit_ideal(): QuaternionFractionalIdeal_rational {
    if (this._unit_ideal === null) {
      this._unit_ideal = new QuaternionFractionalIdeal_rational(
        this.quaternion_algebra(),
        this.basis(),
        this,
        this,
        false
      );
    }
    return this._unit_ideal;
  }

  /**
   * The basis matrix of this order (for the basis returned by `basis()`).
   *
   * @see Reference: quaternion_algebra.py:2593 (basis_matrix)
   */
  basis_matrix(): RationalMatrix {
    return matrix_QQ(this.basis().map((x) => x.coefficient_tuple()));
  }

  /**
   * Every order equals its own unit ideal: multiplication is ideal
   * multiplication.
   *
   * @see Reference: quaternion_algebra.py:2629 (__mul__)
   */
  mul(
    other: QuaternionOrder | QuaternionFractionalIdeal_rational | Quat | Rational | bigint
  ): QuaternionFractionalIdeal_rational {
    return this.unit_ideal().mul(other);
  }

  /**
   * Left multiplication `other * self`.
   *
   * @see Reference: quaternion_algebra.py:2643 (__rmul__)
   */
  rmul(other: Quat | Rational | bigint): QuaternionFractionalIdeal_rational {
    return this.unit_ideal().scale(other, true);
  }

  /**
   * Ideal addition.
   *
   * @see Reference: quaternion_algebra.py:2646 (__add__)
   */
  add(
    other: QuaternionOrder | QuaternionFractionalIdeal_rational
  ): QuaternionFractionalIdeal_rational {
    return this.unit_ideal().add(other);
  }

  /** Equality of orders. @see Reference: quaternion_algebra.py:2169 (__richcmp__) */
  eq(other: unknown): boolean {
    if (!(other instanceof QuaternionOrder)) return false;
    return this.unit_ideal().free_module().eq(other.unit_ideal().free_module());
  }

  /** `self <= other` (containment). @see Reference: quaternion_algebra.py:2169 */
  le(other: QuaternionOrder): boolean {
    return this.unit_ideal().free_module().is_submodule(other.unit_ideal().free_module());
  }

  /** `self >= other`. */
  ge(other: QuaternionOrder): boolean {
    return other.le(this);
  }

  /** `self < other`. */
  lt(other: QuaternionOrder): boolean {
    return this.le(other) && !this.eq(other);
  }

  /** `self > other`. */
  gt(other: QuaternionOrder): boolean {
    return other.lt(this);
  }

  /**
   * The normalized quadratic form associated to this order, as its (integral,
   * primitive) Gram matrix.
   *
   * @see Reference: quaternion_algebra.py:2660 (quadratic_form)
   * @see Deviation: `quadratic_form()` returns a Gram matrix
   */
  quadratic_form(): IntegerMatrix {
    return this.unit_ideal().quadratic_form();
  }

  /**
   * The ternary quadratic form associated to this order, as its Gram matrix.
   *
   * This is the form obtained by restricting the pairing
   * `(x, y) = (conj(x)*y).reduced_trace()` to the trace-zero subspace `G` of
   * `ZZ + 2*self`.
   *
   * @see Reference: quaternion_algebra.py:2680 (ternary_quadratic_form)
   * @see Deviation: `quadratic_form()` returns a Gram matrix
   */
  ternary_quadratic_form(
    include_basis: boolean = false
  ): RationalMatrix | [RationalMatrix, Quat[]] {
    const Q = this.quaternion_algebra();
    const twoR = this.free_module().scale(2n);
    const Z = ZZLattice.span([Q.one().coefficient_tuple()], 4);
    const S = twoR.add(Z);

    // Intersect with the trace-zero submodule: {v in S : reduced_trace(v) = 0}
    const traces = S.basis().map((v) => Q.__call__(v).reduced_trace());
    let d = 1n;
    for (const t of traces) d = lcmBigInt(d, t.denominator);
    const T = integer_matrix_from_rows(traces.map((t) => [t.numerator * (d / t.denominator)]));
    const K = left_kernel_matrix(T);
    const gens: Rational[][] = [];
    for (const k of integer_rows(K)) {
      const v: Rational[] = [Rational.zero(), Rational.zero(), Rational.zero(), Rational.zero()];
      const sb = S.basis();
      for (let i = 0; i < sb.length; i++) {
        const c = k[i] as bigint;
        if (c === 0n) continue;
        for (let j = 0; j < 4; j++) {
          v[j] = (v[j] as Rational).add(
            ((sb[i] as Rational[])[j] as Rational).mul(new Rational(c))
          );
        }
      }
      gens.push(v);
    }
    const G = ZZLattice.span(gens, 4);
    const B = G.basis().map((a) => Q.__call__(a));
    const m = matrix_QQ(B.map((y) => B.map((x) => x.pair(y))));
    if (include_basis) return [m, B];
    return m;
  }

  /**
   * Compute an isomorphism from this order to `other`.
   *
   * Only implemented for maximal orders in definite rational quaternion
   * algebras.  With `conjugator: true` the quaternion `gamma` with
   * `other = gamma^-1 * self * gamma` is returned; otherwise the images of
   * `i, j, k` under the isomorphism together with a function applying it.
   *
   * @see Reference: quaternion_algebra.py:2746 (isomorphism_to)
   * @see Deviation: ring morphisms are returned as image data
   */
  isomorphism_to(
    other: QuaternionOrder,
    options?: { conjugator?: boolean; B?: number }
  ): Quat | { im_gens: [Quat, Quat, Quat]; apply: (x: Quat) => Quat } {
    const conjugator = options?.conjugator ?? false;
    const B = options?.B ?? 10;

    if (!(other instanceof QuaternionOrder)) {
      throw new TypeError('not a quaternion order');
    }
    const Q = this.quaternion_algebra();
    if (!other.quaternion_algebra().eq(Q)) {
      throw new TypeError('not an order in the same quaternion algebra');
    }
    if (!Q.is_definite()) {
      throw new NotImplementedError('only implemented for definite quaternion orders');
    }
    const dQ = new Rational(Q.discriminant());
    if (!(this.discriminant().eq(dQ) && other.discriminant().eq(dQ))) {
      throw new NotImplementedError('only implemented for maximal orders');
    }

    const t1 = this.unit_ideal().theta_series_vector(B);
    const t2 = other.unit_ideal().theta_series_vector(B);
    if (t1.length !== t2.length || !t1.every((c, idx) => c === (t2[idx] as bigint))) {
      throw new ValueError('quaternion orders not isomorphic');
    }

    const attempt = (O1: QuaternionOrder, O2: QuaternionOrder): Quat | null => {
      const N = O1.intersection(O2).free_module().index_in(O1.free_module());
      const I = O1.unit_ideal().mul(O2.unit_ideal()).scale(N, true);
      const gamma = I.minimal_element();
      if (!O1.unit_ideal().scale(gamma, false).eq(I)) return null;
      if (!O2.unit_ideal().scale(gamma, true).eq(I)) return null;
      return gamma;
    };

    const candidates: (Quat | null)[] = [null, ...Q.gens()];
    for (const alpha of candidates) {
      let other_conj = other;
      if (alpha !== null) {
        other_conj = Q.quaternion_order(
          other.basis().map((b) => alpha.mul(b).mul(alpha.inverse()))
        );
      }
      const g = attempt(this, other_conj);
      if (g !== null) {
        const gamma = alpha === null ? g : g.mul(alpha);
        if (conjugator) return gamma;
        const ginv = gamma.inverse();
        const ims = Q.gens().map((gen) => ginv.mul(gen).mul(gamma)) as [Quat, Quat, Quat];
        return {
          im_gens: ims,
          apply: (x: Quat): Quat => ginv.mul(x).mul(gamma),
        };
      }
    }

    const a = Q.invariants()[0].neg();
    const b = Q.invariants()[1].neg();
    const sqfree = [squarefree_part_QQ(a), squarefree_part_QQ(b), squarefree_part_QQ(a.mul(b))];
    const guaranteed = Q.ramified_primes().every((p) => sqfree.some((s) => s === p));
    if (guaranteed) {
      throw new ValueError('quaternion orders not isomorphic');
    }
    throw new NotImplementedError(
      'isomorphism_to was not able to recognize the given orders as isomorphic'
    );
  }
}

// ===========================================================================
// Fractional ideals
// ===========================================================================

/**
 * A fractional ideal in a rational quaternion algebra.
 *
 * @see Reference: quaternion_algebra.py:2980 (QuaternionFractionalIdeal_rational)
 */
export class QuaternionFractionalIdeal_rational {
  private readonly __Q: QuaternionAlgebra_ab;
  private readonly __gens: Quat[];
  private __left_order: QuaternionOrder | null;
  private __right_order: QuaternionOrder | null;
  private __basis_matrix: RationalMatrix | null = null;
  private __free_module: ZZLattice | null = null;
  private __gram_matrix: RationalMatrix | null = null;
  private __theta_series_vector: bigint[] | null = null;

  /**
   * @param Q - the ambient quaternion algebra
   * @param basis - tuple of length 4 of elements whose `ZZ`-span is the ideal
   * @param left_order - a quaternion order or `null`
   * @param right_order - a quaternion order or `null`
   * @param check - if `false`, do no type checking
   *
   * @see Reference: quaternion_algebra.py:2996 (__init__)
   */
  constructor(
    Q: QuaternionAlgebra_ab,
    basis: readonly Quat[],
    left_order: QuaternionOrder | null = null,
    right_order: QuaternionOrder | null = null,
    check: boolean = true
  ) {
    let b = basis.map((v) => Q.__call__(v));
    if (check) {
      if (left_order !== null && !(left_order instanceof QuaternionOrder)) {
        throw new TypeError('left_order must be a quaternion order or None');
      }
      if (right_order !== null && !(right_order instanceof QuaternionOrder)) {
        throw new TypeError('right_order must be a quaternion order or None');
      }
      if (!Array.isArray(basis)) {
        throw new TypeError('basis must be a list or tuple');
      }
      const L = ZZLattice.span(
        b.map((v) => v.coefficient_tuple()),
        4
      );
      b = L.basis().map((v) => Q.__call__(v));
      if (b.length !== 4) {
        throw new ValueError('fractional ideal must have rank 4');
      }
    }
    this.__Q = Q;
    this.__gens = b;
    this.__left_order = left_order;
    this.__right_order = right_order;
  }

  /** The ambient quaternion algebra. @see Reference: quaternion_algebra.py:3107 */
  quaternion_algebra(): QuaternionAlgebra_ab {
    return this.__Q;
  }

  /** Alias of {@link quaternion_algebra}. */
  ring(): QuaternionAlgebra_ab {
    return this.__Q;
  }

  /** The generators of this ideal. */
  gens(): Quat[] {
    return this.__gens.slice();
  }

  /** A basis for this fractional ideal. @see Reference: quaternion_algebra.py:3278 (basis) */
  basis(): Quat[] {
    return this.gens();
  }

  /** @see Reference: quaternion_algebra.py:3251 (__repr__) */
  toString(): string {
    return `Fractional ideal (${this.gens()
      .map((x) => x.toString())
      .join(', ')})`;
  }

  /**
   * Return a random element of this fractional ideal.
   *
   * @see Reference: quaternion_algebra.py:3265 (random_element)
   */
  random_element(x?: bigint, y?: bigint): Quat {
    let s = this.__Q.zero();
    for (const g of this.gens()) {
      s = s.add(g.scalar_mul(ZZ.random_element(x, y)));
    }
    return s;
  }

  /**
   * Scale this fractional ideal by `alpha`, on the right by default.
   *
   * @see Reference: quaternion_algebra.py:3032 (scale)
   */
  scale(
    alpha: Quat | Rational | bigint,
    left: boolean = false
  ): QuaternionFractionalIdeal_rational {
    const Q = this.quaternion_algebra();
    const a =
      alpha instanceof QuaternionAlgebraElement_rational_field ? alpha : Q.__call__(alpha as never);
    if (a.is_zero()) {
      throw new ValueError('the scaling factor must be nonzero');
    }
    const gens = left
      ? basis_for_quaternion_lattice(this.basis().map((b) => a.mul(b)))
      : basis_for_quaternion_lattice(this.basis().map((b) => b.mul(a)));
    const isScalar = a.is_constant();
    const left_order = isScalar || !left ? this.__left_order : null;
    const right_order = isScalar || left ? this.__right_order : null;
    return Q.ideal(gens, { check: false, left_order, right_order });
  }

  /**
   * Compute the left or right order of this ideal.
   *
   * @see Reference: quaternion_algebra.py:3122 (_compute_order)
   */
  _compute_order(side: 'left' | 'right' = 'left'): QuaternionOrder {
    let action: 'left' | 'right';
    if (side === 'left') action = 'right';
    else if (side === 'right') action = 'left';
    else throw new ValueError("side must be 'left' or 'right'");
    const Q = this.quaternion_algebra();
    const M = this.basis().map((b) => b.inverse().matrix(action));
    const B = this.basis_matrix();
    const invs = M.map((m) => B.mul(m));
    const inter = intersection_of_row_modules_over_ZZ(invs);
    const ISB = ZZLattice.span(rows_QQ(inter), 4)
      .basis()
      .map((v) => Q.__call__(v));
    return Q.quaternion_order(ISB);
  }

  /** The left order of this ideal. @see Reference: quaternion_algebra.py:3190 (left_order) */
  left_order(): QuaternionOrder {
    if (this.__left_order === null) {
      this.__left_order = this._compute_order('left');
    }
    return this.__left_order;
  }

  /** The right order of this ideal. @see Reference: quaternion_algebra.py:3215 (right_order) */
  right_order(): QuaternionOrder {
    if (this.__right_order === null) {
      this.__right_order = this._compute_order('right');
    }
    return this.__right_order;
  }

  /**
   * Basis matrix in Hermite normal form for this ideal.
   *
   * @see Reference: quaternion_algebra.py:3363 (basis_matrix)
   */
  basis_matrix(): RationalMatrix {
    if (this.__basis_matrix === null) {
      const B = rational_matrix_from_rational_quaternions(this.gens());
      const [C, d] = clear_denom(rows_QQ(B));
      const H = hnf(C);
      this.__basis_matrix = matrix_QQ(
        integer_rows(H).map((row) => row.map((e) => new Rational(e, d)))
      );
    }
    return this.__basis_matrix;
  }

  /**
   * The underlying free `ZZ`-module of this ideal.
   *
   * @see Reference: quaternion_algebra.py:3709 (free_module)
   */
  free_module(): ZZLattice {
    if (this.__free_module === null) {
      this.__free_module = ZZLattice.span(rows_QQ(this.basis_matrix()), 4);
    }
    return this.__free_module;
  }

  /** Equality of ideals. @see Reference: quaternion_algebra.py:3291 (_richcmp_) */
  eq(other: unknown): boolean {
    if (!(other instanceof QuaternionFractionalIdeal_rational)) return false;
    return this.free_module().eq(other.free_module());
  }

  /** `self <= other` (containment). */
  le(other: QuaternionFractionalIdeal_rational): boolean {
    return this.free_module().is_submodule(other.free_module());
  }

  /** `self >= other`. */
  ge(other: QuaternionFractionalIdeal_rational): boolean {
    return other.le(this);
  }

  /** `self < other`. */
  lt(other: QuaternionFractionalIdeal_rational): boolean {
    return this.le(other) && !this.eq(other);
  }

  /** `self > other`. */
  gt(other: QuaternionFractionalIdeal_rational): boolean {
    return other.lt(this);
  }

  /**
   * An LLL reduced basis of this ideal (definite algebras only).
   *
   * @see Reference: quaternion_algebra.py:3387 (reduced_basis)
   */
  reduced_basis(): Quat[] {
    if (!this.quaternion_algebra().is_definite()) {
      throw new TypeError('the quaternion algebra must be definite');
    }
    const [G] = clear_denom(rows_QQ(this.gram_matrix()));
    const u = lllgramint(integer_rows(G));
    if (u === null) {
      throw new ValueError(
        'qflllgram did not return a square matrix, perhaps the matrix is not positive definite'
      );
    }
    // PARI matrices are column-major: column l of u gives the l-th reduced vector
    const basis = this.basis();
    const out: Quat[] = [];
    for (let l = 0; l < u.length; l++) {
      let s = this.__Q.zero();
      const col = u[l] as bigint[];
      for (let i = 0; i < col.length; i++) {
        s = s.add((basis[i] as Quat).scalar_mul(col[i] as bigint));
      }
      out.push(s);
    }
    return out;
  }

  /**
   * The Gram matrix of this fractional ideal.
   *
   * @see Reference: quaternion_algebra.py:3560 (gram_matrix)
   */
  gram_matrix(): RationalMatrix {
    if (this.__gram_matrix === null) {
      const A = this.gens();
      const two = new Rational(2n);
      this.__gram_matrix = matrix_QQ(A.map((b) => A.map((a) => a.pair(b).mul(two))));
    }
    return this.__gram_matrix;
  }

  /**
   * The normalized quadratic form of this ideal, as its (integral, primitive)
   * Gram matrix.
   *
   * @see Reference: quaternion_algebra.py:3454 (quadratic_form)
   * @see Deviation: `quadratic_form()` returns a Gram matrix
   */
  quadratic_form(): IntegerMatrix {
    const [C] = clear_denom(rows_QQ(this.gram_matrix()));
    let g = 0n;
    for (const row of integer_rows(C)) {
      for (const e of row) g = gcd(g, e);
    }
    if (g !== 1n && g !== 0n) {
      return integer_matrix_from_rows(integer_rows(C).map((row) => row.map((e) => e / g)));
    }
    return C;
  }

  /**
   * The theta series coefficients of this ideal, as a vector of `B` integers.
   *
   * @see Reference: quaternion_algebra.py:3419 (theta_series_vector)
   */
  theta_series_vector(B: number): bigint[] {
    if (this.__theta_series_vector !== null && this.__theta_series_vector.length >= B) {
      return this.__theta_series_vector.slice(0, B);
    }
    const M = integer_rows(this.quadratic_form());
    const rep = B <= 1 ? [] : qfrep(M, BigInt(B - 1), 1);
    const v = [1n, ...rep.map((c) => 2n * c)].slice(0, B);
    this.__theta_series_vector = v;
    return v;
  }

  /**
   * The normalized theta series of this ideal, as a power series over `ZZ`.
   *
   * @see Reference: quaternion_algebra.py:3513 (theta_series)
   */
  theta_series(B: number, varName: string = 'q'): IntegerPowerSeries {
    const v = this.theta_series_vector(B);
    const R = new PowerSeriesRingZZ(ZZ_as_coefficient_ring, varName);
    return new PowerSeriesElementZZ(
      R,
      v.map((c) => new Integer(c)),
      B
    );
  }

  /**
   * An element of minimal norm in this ideal (definite algebras only).
   *
   * @see Reference: quaternion_algebra.py:3485 (minimal_element)
   */
  minimal_element(): Quat {
    if (!this.quaternion_algebra().is_definite()) {
      throw new ValueError('quaternion algebra must be definite');
    }
    // PARI's qfminim(q, NULL, NULL, 1) LLL-reduces the Gram matrix and returns
    // the basis vector of smallest diagonal entry (bibli1.c:1355-1365).
    const M = integer_rows(this.quadratic_form());
    const u = lllgramint(M);
    if (u === null) {
      throw new ValueError('qflllgram failed; the quadratic form must be positive definite');
    }
    const red = qf_ZM_apply(M, u);
    let t = 0;
    for (let i = 1; i < red.length; i++) {
      if (
        (red[i] as bigint[])[i] !== undefined &&
        ((red[i] as bigint[])[i] as bigint) < ((red[t] as bigint[])[t] as bigint)
      ) {
        t = i;
      }
    }
    const col = u[t] as bigint[];
    const basis = this.basis();
    let s = this.__Q.zero();
    for (let i = 0; i < col.length; i++) {
      s = s.add((basis[i] as Quat).scalar_mul(col[i] as bigint));
    }
    return s;
  }

  /**
   * The reduced norm of this fractional ideal.
   *
   * @see Reference: quaternion_algebra.py:3582 (norm)
   */
  norm(): Rational {
    const G = matrix_QQ(
      rows_QQ(this.gram_matrix()).map((r) => r.map((e) => e.div(new Rational(2n))))
    );
    let r = determinant(G).abs();
    if (!r.is_square()) throw new RuntimeError('first is bad!');
    r = r.sqrt() as Rational;
    const R = this.__left_order ?? this.__right_order ?? this.left_order();
    r = r.div(R.discriminant());
    if (!r.is_square()) throw new RuntimeError('second is bad!');
    return r.sqrt() as Rational;
  }

  /**
   * The ideal generated by the conjugates of the generators of this ideal.
   *
   * @see Reference: quaternion_algebra.py:3621 (conjugate)
   */
  conjugate(): QuaternionFractionalIdeal_rational {
    return this.quaternion_algebra().ideal(
      this.basis().map((b) => b.conjugate()),
      { left_order: this.__right_order, right_order: this.__left_order }
    );
  }

  /**
   * The product of this ideal with `right` (an ideal, an order, or a scalar).
   *
   * @see Reference: quaternion_algebra.py:3638 (__mul__)
   */
  mul(
    right: QuaternionFractionalIdeal_rational | QuaternionOrder | Quat | Rational | bigint
  ): QuaternionFractionalIdeal_rational {
    let r: QuaternionFractionalIdeal_rational;
    if (right instanceof QuaternionOrder) {
      r = right.unit_ideal();
    } else if (right instanceof QuaternionFractionalIdeal_rational) {
      r = right;
    } else {
      return this.scale(right as Quat | Rational | bigint, false);
    }
    const gens: Quat[] = [];
    for (const a of this.basis()) for (const b of r.basis()) gens.push(a.mul(b));
    const basis = basis_for_quaternion_lattice(gens);
    return this.quaternion_algebra().ideal(basis, { check: false });
  }

  /**
   * The sum of this ideal and `other`.
   *
   * @see Reference: quaternion_algebra.py:3669 (__add__)
   */
  add(
    other: QuaternionFractionalIdeal_rational | QuaternionOrder
  ): QuaternionFractionalIdeal_rational {
    const o = other instanceof QuaternionOrder ? other.unit_ideal() : other;
    if (!(o instanceof QuaternionFractionalIdeal_rational)) {
      throw new TypeError('can only add quaternion ideals');
    }
    return this.quaternion_algebra().ideal(this.basis().concat(o.basis()));
  }

  /**
   * The intersection of this ideal with `J`.
   *
   * @see Reference: quaternion_algebra.py:3769 (intersection)
   */
  intersection(J: QuaternionFractionalIdeal_rational): QuaternionFractionalIdeal_rational {
    const V = this.free_module().intersection(J.free_module());
    const [H, d] = clear_denom(V.basis());
    const A = this.quaternion_algebra();
    const gens = rational_quaternions_from_integral_matrix_and_denom(A, H, d);
    return A.ideal(gens);
  }

  /**
   * The product of this ideal and the conjugate of `J`.
   *
   * @see Reference: quaternion_algebra.py:3785 (multiply_by_conjugate)
   */
  multiply_by_conjugate(J: QuaternionFractionalIdeal_rational): QuaternionFractionalIdeal_rational {
    const Jbar = J.basis().map((b) => b.conjugate());
    const gens: Quat[] = [];
    for (const a of this.basis()) for (const b of Jbar) gens.push(a.mul(b));
    const basis = basis_for_quaternion_lattice(gens);
    return this.quaternion_algebra().ideal(basis, { check: false });
  }

  /**
   * Whether `x` lies in this ideal.
   *
   * @see Reference: quaternion_algebra.py:4149 (__contains__)
   */
  contains(x: QuaternionLike): boolean {
    let q: Quat;
    try {
      q = this.quaternion_algebra().__call__(x);
    } catch {
      return false;
    }
    return this.free_module().contains(q.coefficient_tuple());
  }

  /**
   * Whether this ideal is integral, i.e. contained in its left order.
   *
   * @see Reference: quaternion_algebra.py:4312 (is_integral)
   */
  is_integral(): boolean {
    if (this.__left_order !== null) {
      return this.free_module().is_submodule(this.left_order().free_module());
    }
    if (this.__right_order !== null) {
      return this.free_module().is_submodule(this.right_order().free_module());
    }
    const sq = this.mul(this);
    return sq.free_module().is_submodule(this.free_module());
  }

  /**
   * Decompose this integral ideal as `(primitive ideal, g)`.
   *
   * @see Reference: quaternion_algebra.py:4342 (primitive_decomposition)
   */
  primitive_decomposition(): [QuaternionFractionalIdeal_rational, bigint] {
    if (!this.is_integral()) {
      throw new ValueError('primitive ideals are defined only for integral ideals');
    }
    const I_basis = this.basis_matrix();
    const O_basis = this.left_order().basis_matrix();
    const M = solve_left(O_basis, I_basis, false);
    let g = 0n;
    for (const row of rows_QQ(M)) {
      for (const e of row) {
        if (e.denominator !== 1n) {
          throw new RuntimeError('ideal is not contained in its left order');
        }
        g = gcd(g, e.numerator);
      }
    }
    if (g === 1n) return [this, g];
    const J = this.scale(new Rational(1n, g));
    return [J, g];
  }

  /**
   * Whether this ideal is primitive.
   *
   * @see Reference: quaternion_algebra.py:4399 (is_primitive)
   */
  is_primitive(): boolean {
    const [, g] = this.primitive_decomposition();
    return g === 1n;
  }

  /**
   * Whether this ideal is principal (definite algebras only).
   *
   * @see Reference: quaternion_algebra.py:4108 (is_principal)
   */
  is_principal(certificate: false): boolean;
  is_principal(certificate: true): [boolean, Quat | null];
  is_principal(certificate?: boolean): boolean | [boolean, Quat | null];
  is_principal(certificate: boolean = false): boolean | [boolean, Quat | null] {
    if (!this.quaternion_algebra().is_definite()) {
      throw new NotImplementedError(
        'principality test not implemented in indefinite quaternion algebras'
      );
    }
    const c = this.theta_series_vector(2)[1] as bigint;
    if (!certificate) return c !== 0n;
    if (c === 0n) return [false, null];
    return [true, this.minimal_element()];
  }

  /**
   * Whether this ideal and `J` are equivalent as right ideals.
   *
   * @see Reference: quaternion_algebra.py:4032 (is_right_equivalent)
   */
  is_right_equivalent(
    J: QuaternionFractionalIdeal_rational,
    B: number = 10,
    certificate: boolean = false
  ): boolean | [boolean, Quat | null] {
    if (!(J instanceof QuaternionFractionalIdeal_rational)) {
      throw new TypeError('J must be a fractional ideal in a rational quaternion algebra');
    }
    if (!this.right_order().eq(J.right_order())) {
      throw new ValueError('self and J must be right ideals over the same order');
    }
    if (!this.quaternion_algebra().is_definite()) {
      throw new NotImplementedError(
        'equivalence test of ideals not implemented for indefinite quaternion algebras'
      );
    }
    if (B > 0) {
      const t1 = this.theta_series_vector(B);
      const t2 = J.theta_series_vector(B);
      if (t1.length !== t2.length || !t1.every((c, idx) => c === (t2[idx] as bigint))) {
        return certificate ? [false, null] : false;
      }
    }
    const IJbar = this.multiply_by_conjugate(J);
    const scaled = IJbar.scale(J.norm().inv(), true);
    return certificate ? scaled.is_principal(true) : scaled.is_principal(false);
  }

  /**
   * Whether this ideal and `J` are equivalent as left ideals.
   *
   * @see Reference: quaternion_algebra.py:3998 (is_left_equivalent)
   */
  is_left_equivalent(
    J: QuaternionFractionalIdeal_rational,
    B: number = 10,
    certificate: boolean = false
  ): boolean | [boolean, Quat | null] {
    if (certificate) {
      const [is_equiv, cert] = this.conjugate().is_right_equivalent(J.conjugate(), B, true) as [
        boolean,
        Quat | null,
      ];
      if (is_equiv && cert !== null) return [true, cert.conjugate()];
      return [false, null];
    }
    return this.conjugate().is_right_equivalent(J.conjugate(), B, false);
  }

  /**
   * The pushforward of this ideal through the ideal `J`.
   *
   * @see Reference: quaternion_algebra.py:3809 (pushforward)
   */
  pushforward(
    J: QuaternionFractionalIdeal_rational,
    side?: 'left' | 'right' | null
  ): QuaternionFractionalIdeal_rational {
    if (!(J instanceof QuaternionFractionalIdeal_rational)) {
      throw new TypeError('can only pushforward through a quaternion ideal');
    }
    if (side === 'left') {
      if (!this.left_order().eq(J.left_order())) {
        throw new ValueError('self and J must have the same left orders');
      }
      if (!this.is_integral() || !J.is_integral()) {
        throw new NotImplementedError(
          'quaternion ideal pushforward not implemented for non-integral ideals'
        );
      }
      const Jnorm = J.norm();
      if (!this.norm().rational_gcd(Jnorm).eq(1n)) {
        throw new ValueError('self and J must have coprime norms');
      }
      return J.conjugate().mul(this.intersection(J)).scale(Jnorm.inv(), true);
    }
    if (side === 'right') {
      if (!this.right_order().eq(J.right_order())) {
        throw new ValueError('self and J must have the same right orders');
      }
      return this.conjugate().pushforward(J.conjugate(), 'left').conjugate();
    }
    if (side === undefined || side === null) {
      const same_left = this.left_order().eq(J.left_order());
      const same_right = this.right_order().eq(J.right_order());
      if (!same_left && !same_right) {
        throw new ValueError('self and J must share a left or right order');
      }
      if (same_left && same_right) {
        throw new ValueError(
          'self and J have same left and right orders, side of pushforward must be specified'
        );
      }
      return same_left ? this.pushforward(J, 'left') : this.pushforward(J, 'right');
    }
    throw new ValueError('side must be "left", "right" or None');
  }

  /**
   * The pullback of this ideal through the ideal `J`.
   *
   * @see Reference: quaternion_algebra.py:3903 (pullback)
   */
  pullback(
    J: QuaternionFractionalIdeal_rational,
    side?: 'left' | 'right' | null
  ): QuaternionFractionalIdeal_rational {
    if (!(J instanceof QuaternionFractionalIdeal_rational)) {
      throw new TypeError('can only pullback through a quaternion ideal');
    }
    if (side === 'left') {
      if (!this.left_order().eq(J.right_order())) {
        throw new ValueError('left order of self should be right order of J');
      }
      if (!this.is_integral() || !J.is_integral()) {
        throw new NotImplementedError(
          'quaternion ideal pullback not implemented for non-integral ideals'
        );
      }
      const N = this.norm();
      if (!N.rational_gcd(J.norm()).eq(1n)) {
        throw new ValueError('self and J must have coprime norms');
      }
      return J.mul(this).add(J.left_order().unit_ideal().scale(N, true));
    }
    if (side === 'right') {
      if (!this.right_order().eq(J.left_order())) {
        throw new ValueError('right order of self should be left order of J');
      }
      return this.conjugate().pullback(J.conjugate(), 'left').conjugate();
    }
    if (side === undefined || side === null) {
      const is_left = this.left_order().eq(J.right_order());
      const is_right = this.right_order().eq(J.left_order());
      if (!is_left && !is_right) {
        throw new ValueError('left order of self must equal right order of J, or vice versa');
      }
      if (is_left && is_right) {
        throw new ValueError(
          'self and J have same left and right orders, side of pullback must be specified'
        );
      }
      return is_left ? this.pullback(J, 'left') : this.pullback(J, 'right');
    }
    throw new ValueError('side must be "left", "right" or None');
  }

  /**
   * The right subideals `J` of this ideal with `I/J` an `GF(p)`-vector space
   * of dimension 2.
   *
   * @see Reference: quaternion_algebra.py:4173 (cyclic_right_subideals)
   */
  cyclic_right_subideals(p: bigint, alpha?: Quat | null): QuaternionFractionalIdeal_rational[] {
    const R = this.right_order();
    const Q = this.quaternion_algebra();
    const basis = basis_for_quaternion_lattice(this.basis(), false);
    const f = Q.modp_splitting_map(p);

    let A: bigint[][];
    let scale: Rational = Rational.one();
    let IB: RationalMatrix;
    const spanRows = (mats: bigint[][][]): bigint[][] =>
      mats.map((m) => [
        (m[0] as bigint[])[0] as bigint,
        (m[0] as bigint[])[1] as bigint,
        (m[1] as bigint[])[0] as bigint,
        (m[1] as bigint[])[1] as bigint,
      ]);
    // The reference catches both ValueError (dependent vectors) and
    // ZeroDivisionError (a basis vector has a denominator divisible by p)
    // and then retries with a rescaled ideal (quaternion_algebra.py:4285-4299).
    let candidate: bigint[][] | null = null;
    try {
      candidate = spanRows(basis.map((a) => f(a)));
      if (!mod_matrix_is_invertible(candidate, p)) candidate = null;
    } catch {
      candidate = null;
    }
    if (candidate !== null) {
      A = candidate;
      IB = matrix_QQ(basis.map((b) => b.coefficient_tuple()));
    } else {
      const [Bm, d] = clear_denom(basis.map((b) => b.coefficient_tuple()));
      let g = 0n;
      for (const row of integer_rows(Bm)) for (const e of row) g = gcd(g, e);
      const IBrows = integer_rows(Bm).map((row) => row.map((e) => new Rational(e / g)));
      IB = matrix_QQ(IBrows);
      scale = new Rational(g, d);
      let cand2: bigint[][];
      try {
        cand2 = spanRows(IBrows.map((r) => f(Q.__call__(r))));
      } catch (e) {
        // Here we could replace the ideal by an *equivalent* ideal that works.
        // This is always possible.  However, upstream has not implemented that
        // algorithm yet either.
        throw new NotImplementedError(
          `general algorithm not implemented (${(e as Error).message})`
        );
      }
      if (!mod_matrix_is_invertible(cand2, p)) {
        throw new NotImplementedError(
          'general algorithm not implemented (the given basis vectors must be linearly independent)'
        );
      }
      A = cand2;
    }

    const Ai = mod_matrix_inverse(A, p);
    // Lift the mod-p inverse to QQ (Sage: Ai.change_ring(QQ)) and multiply
    const AiQ = matrix_QQ(Ai.map((row) => row.map((e) => new Rational(mod(e, p)))));
    const AiB = AiQ.mul(IB);
    const [AiBint] = clear_denom(rows_QQ(AiB));
    const AiBrows = integer_rows(AiBint);

    const pBrows0 = rows_QQ(IB).map((row) => row.map((e) => e.mul(new Rational(p))));
    const [pB, d2] = clear_denom(pBrows0);
    const pBrows = integer_rows(pB);

    const ans: QuaternionFractionalIdeal_rational[] = [];
    let lines: [bigint, bigint][];
    if (alpha === undefined || alpha === null) {
      lines = p1list(p);
    } else {
      let x = f(alpha);
      const a = f(alpha);
      lines = [];
      for (let n = 0n; n <= p; n++) {
        lines.push(
          p1_normalize(p, (x[0] as bigint[])[0] as bigint, (x[0] as bigint[])[1] as bigint)
        );
        x = mat_mul(x, a, p);
      }
    }

    for (const [u, v] of lines) {
      // z = matrix(QQ,2,4,[0,-v,0,u, -v,0,u,0]) * AiB
      const Z: bigint[][] = [
        [0n, mod(-v, p), 0n, u],
        [mod(-v, p), 0n, u, 0n],
      ];
      const z: bigint[][] = [];
      for (let r = 0; r < 2; r++) {
        const row: bigint[] = [0n, 0n, 0n, 0n];
        for (let c = 0; c < 4; c++) {
          let s = 0n;
          for (let t = 0; t < 4; t++) {
            s += ((Z[r] as bigint[])[t] as bigint) * ((AiBrows[t] as bigint[])[c] as bigint);
          }
          row[c] = s;
        }
        z.push(row);
      }
      const G = integer_matrix_from_rows(z.map((row) => row.map((e) => e * d2)).concat(pBrows));
      const H = hnf(G);
      let gens = rational_quaternions_from_integral_matrix_and_denom(Q, H, d2);
      if (!scale.eq(1n)) {
        gens = gens.map((gg) => gg.scalar_mul(scale));
      }
      ans.push(R.right_ideal(gens, false));
    }
    return ans;
  }
}

// ===========================================================================
// Utility functions
// ===========================================================================

/**
 * Return a basis for the `ZZ`-lattice in a quaternion algebra spanned by the
 * given generators.
 *
 * @param gens - list of elements of a single quaternion algebra
 * @param reverse - when computing the HNF do it on the basis `(k,j,i,1)`
 *   instead of `(1,i,j,k)`
 *
 * @example
 * ```typescript
 * // A.<i,j,k> = QuaternionAlgebra(-1,-7)
 * // basis_for_quaternion_lattice([i+j, i-j, 2*k, A(1/3)])
 * // [1/3, 2*i, i + j, 2*k]
 * ```
 *
 * @see Reference: quaternion_algebra.py:4425 (basis_for_quaternion_lattice)
 */
export function basis_for_quaternion_lattice(
  gens: readonly Quat[],
  reverse: boolean = true
): Quat[] {
  if (gens.length === 0) return [];
  const [Z, d] = integral_matrix_and_denom_from_rational_quaternions(gens, reverse);
  const H = hnf(Z);
  const A = (gens[0] as Quat).parent();
  return rational_quaternions_from_integral_matrix_and_denom(A, H, d, reverse);
}

/**
 * Intersect the `ZZ`-modules with basis matrices the full rank 4x4 `QQ`
 * matrices in the list `v`.
 *
 * @see Reference: quaternion_algebra.py:4457 (intersection_of_row_modules_over_ZZ)
 */
export function intersection_of_row_modules_over_ZZ(v: readonly RationalMatrix[]): RationalMatrix {
  if (v.length <= 0) {
    throw new ValueError('v must have positive length');
  }
  if (v.length === 1) {
    return v[0] as RationalMatrix;
  }
  if (v.length === 2) {
    const a = v[0] as RationalMatrix;
    const b = v[1] as RationalMatrix;
    const stacked = rows_QQ(a).concat(rows_QQ(b));
    const [s] = clear_denom(stacked);
    const K = kernel_matrix(s.transpose());
    const n = a.nrows;
    const Kfirst = matrix_QQ(
      integer_rows(K).map((row) => row.slice(0, n).map((e) => new Rational(e)))
    );
    return Kfirst.mul(a);
  }
  const w = intersection_of_row_modules_over_ZZ(v.slice(0, 2));
  return intersection_of_row_modules_over_ZZ([w, ...v.slice(2)]);
}

/**
 * Compute a basis of a `ZZ`-module that is normalized at `p`.
 *
 * @param e - basis of a `ZZ`-module (this array is modified)
 * @param p - prime at which the basis should be normalized
 * @param B - bilinear form (default: `x.pair(y)`)
 * @returns pairs `(basis element, valuation of its orthogonal summand)`
 *
 * @see Reference: quaternion_algebra.py:4501 (normalize_basis_at_p)
 */
export function normalize_basis_at_p(
  e: Quat[],
  p: bigint,
  B: (x: Quat, y: Quat) => Rational = (x, y) => x.pair(y)
): [Quat, bigint][] {
  const N = e.length;
  if (N === 0) return [];

  let min_m = 0;
  let min_n = 0;
  let min_v: bigint | null = null; // null = +Infinity

  for (let m = 0; m < N; m++) {
    for (let n = m; n < N; n++) {
      const v = rational_valuation_or_null(B(e[m] as Quat, e[n] as Quat), p);
      if (val_lt(v, min_v) || (val_eq(v, min_v) && min_m !== min_n && m === n)) {
        min_m = m;
        min_n = n;
        min_v = v;
      }
    }
  }
  if (min_v === null) {
    throw new ValueError('the bilinear form is identically zero on the given basis');
  }

  if (min_m === min_n || p !== 2n) {
    let f0: Quat;
    if (min_m === min_n) {
      f0 = e[min_m] as Quat;
    } else {
      f0 = (e[min_m] as Quat).add(e[min_n] as Quat);
    }

    [e[0], e[min_m]] = [e[min_m] as Quat, e[0] as Quat];

    const c = B(f0, f0);
    for (let l = 1; l < N; l++) {
      e[l] = (e[l] as Quat).sub(f0.scalar_mul(B(e[l] as Quat, f0).div(c)));
    }

    const f = normalize_basis_at_p(e.slice(1), p, B);
    f.unshift([f0, min_v - integer_valuation(p, 2n)]);
    return f;
  }

  // p = 2 and only off-diagonal entries have minimal valuation
  if (
    val_gt(
      rational_valuation_or_null(B(e[min_m] as Quat, e[min_m] as Quat), p),
      rational_valuation_or_null(B(e[min_n] as Quat, e[min_n] as Quat), p)
    )
  ) {
    [e[min_m], e[min_n]] = [e[min_n] as Quat, e[min_m] as Quat];
  }

  let f0 = (e[min_m] as Quat).scalar_mul(
    new Rational(pow_bigint(p, min_v)).div(B(e[min_m] as Quat, e[min_n] as Quat))
  );
  let f1 = e[min_n] as Quat;

  const v01 = rational_valuation_or_null(B(f0, f1), p);
  if (val_lt(v01 === null ? null : v01 + 1n, rational_valuation_or_null(B(f0, f0), p))) {
    const g = f0;
    f0 = f0.add(f1);
    f1 = g;
  }

  e[min_m] = e[0] as Quat;
  e[min_n] = e[1] as Quat;

  const B00 = B(f0, f0);
  const B11 = B(f1, f1);
  const B01 = B(f0, f1);
  const d = B00.mul(B11).sub(B01.mul(B01));
  const tu: [Rational, Rational][] = [];
  for (let l = 2; l < N; l++) {
    tu.push([
      B01.mul(B(f1, e[l] as Quat)).sub(B11.mul(B(f0, e[l] as Quat))),
      B01.mul(B(f0, e[l] as Quat)).sub(B00.mul(B(f1, e[l] as Quat))),
    ]);
  }
  for (let l = 2; l < N; l++) {
    const [t, u] = tu[l - 2] as [Rational, Rational];
    e[l] = (e[l] as Quat).add(f0.scalar_mul(t.div(d))).add(f1.scalar_mul(u.div(d)));
  }

  const f = normalize_basis_at_p(e.slice(2, N), p, B);
  return (
    [
      [f0, min_v],
      [f1, min_v],
    ] as [Quat, bigint][]
  ).concat(f);
}

/**
 * Given `a`, `b` and the even prime `p = 2`, find `(y, z, w)` with `y` a unit
 * mod `p^2` such that `1 - a y^2 - b z^2 + a b w^2 = 0 mod p^2`.
 *
 * @see Reference: quaternion_algebra.py:4638 (maxord_solve_aux_eq)
 */
export function maxord_solve_aux_eq(a: bigint, b: bigint, p: bigint): [bigint, bigint, bigint] {
  if (p !== 2n) {
    throw new NotImplementedError('algorithm only implemented over ZZ at the moment');
  }
  const v_a = a === 0n ? -1n : integer_valuation(a, p);
  const v_b = b === 0n ? -1n : integer_valuation(b, p);
  if (v_a !== 0n) {
    throw new RuntimeError('a must have v_p(a)=0');
  }
  if (v_b !== 0n && v_b !== 1n) {
    throw new RuntimeError('b must have v_p(b) in {0,1}');
  }
  const key = `${mod(a, 4n)},${mod(b, 4n)}`;
  const lut: Record<string, [bigint, bigint, bigint]> = {
    '1,1': [1n, 1n, 1n],
    '1,2': [1n, 0n, 0n],
    '1,3': [1n, 0n, 0n],
    '3,1': [1n, 1n, 1n],
    '3,2': [1n, 0n, 1n],
    '3,3': [1n, 1n, 1n],
  };
  const r = lut[key];
  if (r === undefined) {
    throw new ValueError(`no solution for (a, b) = (${a}, ${b}) mod 4`);
  }
  return r;
}

// ===========================================================================
// Small helpers
// ===========================================================================

/** Python-style modulo (result has the sign of the modulus). */
function mod(a: bigint, m: bigint): bigint {
  const r = a % m;
  return r < 0n === m < 0n || r === 0n ? r : r + m;
}

function pow_bigint(base: bigint, e: bigint): bigint {
  if (e < 0n) throw new ValueError('negative exponent');
  let r = 1n;
  let b = base;
  let n = e;
  while (n > 0n) {
    if (n & 1n) r *= b;
    b *= b;
    n >>= 1n;
  }
  return r;
}

/** Python's `//` for bigints (floor division). */
function floorDivBigInt(a: bigint, b: bigint): bigint {
  const q = a / b;
  if (a % b !== 0n && a < 0n !== b < 0n) return q - 1n;
  return q;
}

/** The `p`-adic valuation of a nonzero rational. */
function rational_valuation(r: Rational, p: bigint): bigint {
  if (r.isZero()) {
    throw new ValueError('valuation of zero is infinite');
  }
  return integer_valuation(r.numerator, p) - integer_valuation(r.denominator, p);
}

function rational_valuation_or_null(r: Rational, p: bigint): bigint | null {
  if (r.isZero()) return null;
  return integer_valuation(r.numerator, p) - integer_valuation(r.denominator, p);
}

function val_lt(a: bigint | null, b: bigint | null): boolean {
  if (a === null) return false;
  if (b === null) return true;
  return a < b;
}

function val_eq(a: bigint | null, b: bigint | null): boolean {
  if (a === null || b === null) return a === b;
  return a === b;
}

function val_gt(a: bigint | null, b: bigint | null): boolean {
  return val_lt(b, a);
}

/** The rational constant term of a product `x*y` known to be a scalar. */
function constant_of(x: Quat, y: Quat): Rational {
  const prod = x.mul(y);
  if (!prod.is_constant()) {
    throw new ValueError('element is not a scalar');
  }
  return prod.get(0);
}

/** Convert a rational known to be integral to a bigint. */
function to_integer(r: Rational): bigint {
  if (r.denominator !== 1n) {
    throw new ValueError(`no conversion of ${r.toString()} to integer`);
  }
  return r.numerator;
}

/** The squarefree part of a positive rational (used by `isomorphism_to`). */
function squarefree_part_QQ(r: Rational): bigint {
  const n = r.numerator * r.denominator;
  let s = 1n;
  const sign = n < 0n ? -1n : 1n;
  const m = absBigInt(n);
  for (const [p, e] of factor(m)) {
    if (e % 2n === 1n) s *= p;
  }
  void m;
  return sign * s;
}

/** The image of a rational in `GF(p)`. */
function to_Fp(r: Rational, p: bigint): bigint {
  const den = mod(r.denominator, p);
  if (den === 0n) {
    throw new ValueError('denominator is not invertible mod p');
  }
  return mod(mod(r.numerator, p) * mod_inverse(den, p), p);
}

function mod_inverse(a: bigint, p: bigint): bigint {
  let [old_r, r] = [mod(a, p), p];
  let [old_s, s] = [1n, 0n];
  while (r !== 0n) {
    const q = old_r / r;
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
  }
  if (old_r !== 1n && old_r !== -1n) {
    throw new ValueError('element is not invertible mod p');
  }
  return mod(old_s * old_r, p);
}

function mat_mul(A: bigint[][], B: bigint[][], p: bigint): bigint[][] {
  const n = A.length;
  const out: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < n; j++) {
      let s = 0n;
      for (let k = 0; k < n; k++) {
        s += ((A[i] as bigint[])[k] as bigint) * ((B[k] as bigint[])[j] as bigint);
      }
      row.push(mod(s, p));
    }
    out.push(row);
  }
  return out;
}

function mat_neg(A: bigint[][], p: bigint): bigint[][] {
  return A.map((row) => row.map((e) => mod(-e, p)));
}

function mat_eq(A: bigint[][], B: bigint[][]): boolean {
  return A.every((row, i) => row.every((e, j) => e === ((B[i] as bigint[])[j] as bigint)));
}

function scalar_2x2(c: bigint): bigint[][] {
  return [
    [c, 0n],
    [0n, c],
  ];
}

/**
 * Integer vectors of length `n` with entries in `[0, max_part]` summing to
 * `w`, in SageMath's `IntegerVectors` order (decreasing lexicographic).
 *
 * @see Reference: sage/combinat/integer_vector.py (IntegerVectors)
 */
function* integer_vectors(w: number, n: number, max_part: number): Generator<number[]> {
  if (n === 0) {
    if (w === 0) yield [];
    return;
  }
  const hi = Math.min(max_part, w);
  const lo = Math.max(0, w - max_part * (n - 1));
  for (let v = hi; v >= lo; v--) {
    for (const rest of integer_vectors(w - v, n - 1, max_part)) {
      yield [v, ...rest];
    }
  }
}

/**
 * All 2x2 matrices over `GF(p)`, in SageMath's `MatrixSpace` iteration order
 * (by weight, then by `IntegerVectors` order).
 *
 * @see Reference: sage/matrix/matrix_space.py:1548 (MatrixSpace.__iter__)
 */
function* all_2x2_matrices(p: bigint): Generator<bigint[][]> {
  const order = Number(p);
  for (let weight = 0; weight <= (order - 1) * 4; weight++) {
    for (const iv of integer_vectors(weight, 4, order - 1)) {
      yield [
        [BigInt(iv[0] as number), BigInt(iv[1] as number)],
        [BigInt(iv[2] as number), BigInt(iv[3] as number)],
      ];
    }
  }
}

/** All vectors of `GF(p)^n`, in SageMath's iteration order. */
function* gf_vectors(p: bigint, n: number): Generator<bigint[]> {
  const v = new Array<bigint>(n).fill(0n);
  for (;;) {
    yield v.slice();
    let i = 0;
    while (i < n) {
      v[i] = ((v[i] as bigint) + 1n) % p;
      if (v[i] !== 0n) break;
      i++;
    }
    if (i === n) return;
  }
}

/** Whether the 4x4 matrix over `GF(p)` is invertible. */
function mod_matrix_is_invertible(A: bigint[][], p: bigint): boolean {
  try {
    mod_matrix_inverse(A, p);
    return true;
  } catch {
    return false;
  }
}

/** Inverse of a square matrix over `GF(p)`. */
function mod_matrix_inverse(A: bigint[][], p: bigint): bigint[][] {
  const n = A.length;
  const M = A.map((row, i) =>
    row
      .map((e) => mod(e, p))
      .concat(new Array<bigint>(n).fill(0n).map((_, j) => (i === j ? 1n : 0n)))
  );
  for (let col = 0; col < n; col++) {
    let piv = -1;
    for (let r = col; r < n; r++) {
      if (((M[r] as bigint[])[col] as bigint) !== 0n) {
        piv = r;
        break;
      }
    }
    if (piv === -1) throw new ValueError('matrix is not invertible mod p');
    [M[col], M[piv]] = [M[piv] as bigint[], M[col] as bigint[]];
    const inv = mod_inverse((M[col] as bigint[])[col] as bigint, p);
    for (let j = 0; j < 2 * n; j++) {
      (M[col] as bigint[])[j] = mod(((M[col] as bigint[])[j] as bigint) * inv, p);
    }
    for (let r = 0; r < n; r++) {
      if (r === col) continue;
      const f = (M[r] as bigint[])[col] as bigint;
      if (f === 0n) continue;
      for (let j = 0; j < 2 * n; j++) {
        (M[r] as bigint[])[j] = mod(
          ((M[r] as bigint[])[j] as bigint) - f * ((M[col] as bigint[])[j] as bigint),
          p
        );
      }
    }
  }
  return M.map((row) => row.slice(n));
}

/**
 * The projective line `P^1(GF(p))` in SageMath's `P1List` order.
 *
 * SageMath uses `sage.modular.modsym.p1list.P1List`, which is not part of this
 * port; for prime `p` its list is `[(0,1), (1,0), (1,1), ..., (1,p-1)]`.
 *
 * @see Reference: sage/modular/modsym/p1list.pyx (P1List)
 * @see Deviation: P1List is inlined for prime level
 */
function p1list(p: bigint): [bigint, bigint][] {
  const out: [bigint, bigint][] = [[0n, 1n]];
  for (let v = 0n; v < p; v++) out.push([1n, v]);
  return out;
}

/** `P1List.normalize` for prime level `p`. */
function p1_normalize(p: bigint, u: bigint, v: bigint): [bigint, bigint] {
  const uu = mod(u, p);
  const vv = mod(v, p);
  if (uu === 0n) {
    if (vv === 0n) return [0n, 0n];
    return [0n, 1n];
  }
  return [1n, mod(vv * mod_inverse(uu, p), p)];
}

/** Thrown where SageMath raises `IndexError`. */
class IndexError extends ValueError {
  override name = 'IndexError';
}
