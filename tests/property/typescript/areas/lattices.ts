/**
 * sagemath-ts side of the `lattices` property-test area.
 *
 * Covers `src/modules/free_module.ts`, `free_module_element.ts`,
 * `free_module_integer.ts` and `bkz.ts`.
 *
 * Cases: tests/property/cases/lattices.cases.json
 * SageMath counterpart: tests/property/python/areas/lattices.py
 *
 * Every function returns an **already formatted string**, mirroring the
 * SageMath area module character for character, so the two transcripts can be
 * compared byte-for-byte without touching the shared runners.
 *
 * Matrices arrive flattened row-major together with the column count, because
 * the shared argument generators only produce integers and flat integer lists.
 */

import {
  IntegerMatrixFromEntries,
  hermite_normal_form,
} from '../../../../packages/sagemath-ts/src/matrix/matrix_integer.js';
import {
  FreeModule,
  type FreeModuleField,
  type FreeModuleGeneric,
  type FreeModulePID,
  VectorSpace,
} from '../../../../packages/sagemath-ts/src/modules/free_module.js';
import type { RingLike } from '../../../../packages/sagemath-ts/src/modules/free_module_element.js';
import { IntegerLattice } from '../../../../packages/sagemath-ts/src/modules/free_module_integer.js';
import { GF } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_constructor.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';

// ---------------------------------------------------------------------------
// Base rings
//
// The port models a base ring structurally; these are the same literals the
// unit tests in packages/sagemath-ts/src/modules/free_module.test.ts use.
// ---------------------------------------------------------------------------

const ZZ: RingLike = {
  zero: () => 0n,
  one: () => 1n,
  is_field: () => false,
  toString: () => 'Integer Ring',
};

const QQ: RingLike = {
  zero: () => Rational.zero(),
  one: () => Rational.one(),
  is_field: () => true,
  toString: () => 'Rational Field',
};

// ---------------------------------------------------------------------------
// Formatting helpers (mirrored in tests/property/python/areas/lattices.py)
// ---------------------------------------------------------------------------

const fmtBool = (b: boolean): string => (b ? 'True' : 'False');

const fmtList = (xs: readonly unknown[]): string => `[${xs.map((x) => String(x)).join(', ')}]`;

const fmtMat = (rows: readonly (readonly unknown[])[]): string =>
  `[${rows.map((row) => fmtList(row)).join(', ')}]`;

/** SageMath prints `+Infinity`; the port uses `Number.POSITIVE_INFINITY`.
 *  @see DEVIATIONS.md "Infinity Representation" */
const fmtScalar = (x: unknown): string => {
  if (typeof x === 'number' && x === Number.POSITIVE_INFINITY) {
    return '+Infinity';
  }
  if (typeof x === 'number' && x === Number.NEGATIVE_INFINITY) {
    return '-Infinity';
  }
  return String(x);
};

const fmtModule = (W: FreeModuleGeneric): string =>
  `rank=${W.rank()} basis=${fmtMat(W.basisMatrix() as unknown[][])}`;

const matRows = (M: {
  nrows: number;
  ncols: number;
  get: (i: number, j: number) => { value: bigint };
}): bigint[][] => {
  const out: bigint[][] = [];
  for (let i = 0; i < M.nrows; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < M.ncols; j++) {
      row.push(M.get(i, j).value);
    }
    out.push(row);
  }
  return out;
};

// ---------------------------------------------------------------------------
// Input reshaping
// ---------------------------------------------------------------------------

function reshape<T>(flat: readonly T[], ncols: bigint): T[][] {
  const n = Number(ncols);
  if (n <= 0) {
    throw new Error('ncols must be positive');
  }
  if (flat.length % n !== 0) {
    throw new Error(`flat length ${flat.length} is not a multiple of ncols ${n}`);
  }
  const rows: T[][] = [];
  for (let i = 0; i < flat.length / n; i++) {
    rows.push(flat.slice(i * n, (i + 1) * n) as T[]);
  }
  return rows;
}

/** `flat` holds (numerator, denominator) pairs, row-major. */
function reshapeQQ(flat: readonly bigint[], ncols: bigint): Rational[][] {
  if (flat.length % 2 !== 0) {
    throw new Error('rational payload must have even length');
  }
  const entries: Rational[] = [];
  for (let i = 0; i < flat.length / 2; i++) {
    entries.push(new Rational(flat[2 * i]!, flat[2 * i + 1]!));
  }
  return reshape(entries, ncols);
}

const zzModule = (ncols: bigint) => FreeModule(ZZ, Number(ncols));
const qqSpace = (ncols: bigint) => VectorSpace(QQ, Number(ncols));
const gfSpace = (ncols: bigint, p: bigint) =>
  VectorSpace(GF(p) as unknown as RingLike, Number(ncols));

const elements = (V: FreeModuleGeneric, rows: unknown[][]) => rows.map((r) => V.createElement(r));

// ---------------------------------------------------------------------------
// modules/free_module.ts -- spans, echelon bases, rank
// ---------------------------------------------------------------------------

function zz_span(flat: bigint[], ncols: bigint): string {
  const V = zzModule(ncols);
  return fmtModule(V.span(elements(V, reshape(flat, ncols))));
}

function zz_span_qq_entries(flat: bigint[], ncols: bigint): string {
  const V = zzModule(ncols);
  return fmtModule(V.span(elements(V, reshapeQQ(flat, ncols))));
}

function zz_span_of_basis(flat: bigint[], ncols: bigint): string {
  const V = zzModule(ncols) as FreeModulePID;
  const W = V.spanOfBasis(elements(V, reshape(flat, ncols)));
  return `basis=${fmtMat(W.basisMatrix() as unknown[][])} echelon=${fmtMat(
    W.echelonizedBasisMatrix()
  )}`;
}

function qq_span(flat: bigint[], ncols: bigint): string {
  const V = qqSpace(ncols);
  const rows = reshape(flat, ncols).map((r) => r.map((x) => new Rational(x, 1n)));
  const W = V.subspace(elements(V, rows));
  return `dim=${W.dimension()} basis=${fmtMat(W.basisMatrix() as unknown[][])}`;
}

function qq_span_qq_entries(flat: bigint[], ncols: bigint): string {
  const V = qqSpace(ncols);
  const W = V.subspace(elements(V, reshapeQQ(flat, ncols)));
  return `dim=${W.dimension()} basis=${fmtMat(W.basisMatrix() as unknown[][])}`;
}

function gf_span(flat: bigint[], ncols: bigint, p: bigint): string {
  const V = gfSpace(ncols, p);
  const W = V.subspace(elements(V, reshape(flat, ncols)));
  return `dim=${W.dimension()} basis=${fmtMat(W.basisMatrix() as unknown[][])}`;
}

function zz_intersection(a: bigint[], b: bigint[], ncols: bigint): string {
  const V = zzModule(ncols);
  const A = V.span(elements(V, reshape(a, ncols))) as FreeModulePID;
  const B = V.span(elements(V, reshape(b, ncols))) as FreeModulePID;
  return `${fmtModule(A.intersection(B))} | ${fmtModule(B.intersection(A))}`;
}

function qq_intersection(a: bigint[], b: bigint[], ncols: bigint): string {
  const V = qqSpace(ncols);
  const toRows = (flat: bigint[]) =>
    reshape(flat, ncols).map((r) => r.map((x) => new Rational(x, 1n)));
  const A = V.subspace(elements(V, toRows(a)));
  const B = V.subspace(elements(V, toRows(b)));
  return `${fmtModule(A.intersection(B))} | ${fmtModule(B.intersection(A))}`;
}

function gf_intersection(a: bigint[], b: bigint[], ncols: bigint, p: bigint): string {
  const V = gfSpace(ncols, p);
  const A = V.subspace(elements(V, reshape(a, ncols)));
  const B = V.subspace(elements(V, reshape(b, ncols)));
  return `${fmtModule(A.intersection(B))} | ${fmtModule(B.intersection(A))}`;
}

function zz_sum(a: bigint[], b: bigint[], ncols: bigint): string {
  const V = zzModule(ncols);
  const A = V.span(elements(V, reshape(a, ncols)));
  const B = V.span(elements(V, reshape(b, ncols)));
  return fmtModule(A.add(B));
}

function zz_discriminant(flat: bigint[], ncols: bigint): string {
  const V = zzModule(ncols);
  return String(V.span(elements(V, reshape(flat, ncols))).discriminant());
}

function zz_saturation(flat: bigint[], ncols: bigint): string {
  const V = zzModule(ncols);
  const W = V.span(elements(V, reshape(flat, ncols))) as FreeModulePID;
  return fmtModule(W.saturation());
}

function zz_index_in(a: bigint[], b: bigint[], ncols: bigint): string {
  const V = zzModule(ncols);
  const A = V.span(elements(V, reshape(a, ncols)));
  const B = V.span(elements(V, reshape(b, ncols)));
  return fmtScalar(A.indexIn(B));
}

function zz_index_in_ambient(flat: bigint[], ncols: bigint): string {
  const V = zzModule(ncols);
  const W = V.span(elements(V, reshape(flat, ncols)));
  return fmtScalar(W.indexIn(V));
}

function zz_coordinates(flat: bigint[], ncols: bigint, target: bigint[]): string {
  const V = zzModule(ncols);
  const W = V.span(elements(V, reshape(flat, ncols)));
  return fmtList(W.coordinates(V.createElement(target)));
}

function zz_coordinates_user_basis(flat: bigint[], ncols: bigint, target: bigint[]): string {
  const V = zzModule(ncols) as FreeModulePID;
  const W = V.spanOfBasis(elements(V, reshape(flat, ncols)));
  return fmtList(W.coordinates(V.createElement(target)));
}

function zz_is_submodule(a: bigint[], b: bigint[], ncols: bigint): string {
  const V = zzModule(ncols);
  const A = V.span(elements(V, reshape(a, ncols)));
  const B = V.span(elements(V, reshape(b, ncols)));
  return `${fmtBool(A.isSubmodule(B))} ${fmtBool(B.isSubmodule(A))}`;
}

function zz_module_eq(a: bigint[], b: bigint[], ncols: bigint): string {
  const V = zzModule(ncols);
  const A = V.span(elements(V, reshape(a, ncols)));
  const B = V.span(elements(V, reshape(b, ncols)));
  return fmtBool(A.equals(B));
}

function zz_quotient_invariants(flat: bigint[], ncols: bigint): string {
  const V = zzModule(ncols);
  const W = V.span(elements(V, reshape(flat, ncols)));
  return fmtList(V.quotientModule(W).invariants());
}

function qq_complement(flat: bigint[], ncols: bigint): string {
  const V = qqSpace(ncols);
  const rows = reshape(flat, ncols).map((r) => r.map((x) => new Rational(x, 1n)));
  const W = V.subspace(elements(V, rows)) as FreeModuleField;
  return fmtMat(W.complement().basisMatrix() as unknown[][]);
}

function gf_complement(flat: bigint[], ncols: bigint, p: bigint): string {
  const V = gfSpace(ncols, p);
  const W = V.subspace(elements(V, reshape(flat, ncols))) as FreeModuleField;
  return fmtMat(W.complement().basisMatrix() as unknown[][]);
}

function gf_cardinality(n: bigint, p: bigint): string {
  return fmtScalar(FreeModule(GF(p) as unknown as RingLike, Number(n)).cardinality());
}

// ---------------------------------------------------------------------------
// modules/free_module_element.ts
// ---------------------------------------------------------------------------

function vector_ops_zz(a: bigint[], b: bigint[]): string {
  const V = FreeModule(ZZ, a.length);
  const u = V.createElement(a);
  const v = V.createElement(b);
  return [
    `add=${fmtList(u.add(v).list())}`,
    `sub=${fmtList(u.sub(v).list())}`,
    `neg=${fmtList(u.neg().list())}`,
    `dot=${String(u.dotProduct(v))}`,
    `pairwise=${fmtList(u.pairwiseProduct(v).list())}`,
    `hw=${u.hammingWeight()}`,
    `support=${fmtList(u.support())}`,
    `norm2=${String(u.dotProduct(u))}`,
  ].join(' ');
}

function vector_cross_zz(a: bigint[], b: bigint[]): string {
  const V = FreeModule(ZZ, a.length);
  return fmtList(V.createElement(a).crossProduct(V.createElement(b)).list());
}

function vector_ops_gf(a: bigint[], b: bigint[], p: bigint): string {
  const V = gfSpace(BigInt(a.length), p);
  const u = V.createElement(a);
  const v = V.createElement(b);
  return [
    `add=${fmtList(u.add(v).list())}`,
    `dot=${String(u.dotProduct(v))}`,
    `pairwise=${fmtList(u.pairwiseProduct(v).list())}`,
    `is_zero=${fmtBool(u.isZero())}`,
    `hw=${u.hammingWeight()}`,
  ].join(' ');
}

function vector_scalar_zz(a: bigint[], c: bigint): string {
  const V = FreeModule(ZZ, a.length);
  return fmtList(V.createElement(a).mul(c).list());
}

// ---------------------------------------------------------------------------
// Exact LLL-reducedness predicate
//
// Deliberately written from the definition with exact rational arithmetic, so
// that it is an *independent* check on the reduction routines under test
// rather than a re-use of the port's own verifier.  delta / eta are the
// SageMath defaults 0.99 / 0.501 as exact rationals.
// ---------------------------------------------------------------------------

interface Frac {
  n: bigint;
  d: bigint;
}

function bgcd(a: bigint, b: bigint): bigint {
  let x = a < 0n ? -a : a;
  let y = b < 0n ? -b : b;
  while (y) {
    [x, y] = [y, x % y];
  }
  return x;
}

function frac(n: bigint, d: bigint = 1n): Frac {
  if (d === 0n) {
    throw new Error('zero denominator');
  }
  if (d < 0n) {
    n = -n;
    d = -d;
  }
  const g = bgcd(n, d) || 1n;
  return { n: n / g, d: d / g };
}

const fAdd = (a: Frac, b: Frac): Frac => frac(a.n * b.d + b.n * a.d, a.d * b.d);
const fSub = (a: Frac, b: Frac): Frac => frac(a.n * b.d - b.n * a.d, a.d * b.d);
const fMul = (a: Frac, b: Frac): Frac => frac(a.n * b.n, a.d * b.d);
const fDiv = (a: Frac, b: Frac): Frac => frac(a.n * b.d, a.d * b.n);
const fAbs = (a: Frac): Frac => ({ n: a.n < 0n ? -a.n : a.n, d: a.d });
/** sign of a - b */
const fCmp = (a: Frac, b: Frac): number => {
  const v = a.n * b.d - b.n * a.d;
  return v > 0n ? 1 : v < 0n ? -1 : 0;
};
const fIsZero = (a: Frac): boolean => a.n === 0n;

const DELTA: Frac = frac(99n, 100n);
const ETA: Frac = frac(501n, 1000n);

function isLLLReducedExact(rows: readonly bigint[][]): boolean {
  const n = rows.length;
  if (n <= 1) {
    return true;
  }
  const m = rows[0]!.length;
  // Gram-Schmidt over QQ: gs[i] = b_i - sum_{j<i} mu[i][j] * gs[j]
  const gs: Frac[][] = [];
  const mu: Frac[][] = [];
  const normsq: Frac[] = [];
  for (let i = 0; i < n; i++) {
    let cur: Frac[] = rows[i]!.map((x) => frac(x));
    mu.push([]);
    for (let j = 0; j < i; j++) {
      if (fIsZero(normsq[j]!)) {
        // Linearly dependent input: the definition does not apply.
        throw new Error('basis is linearly dependent');
      }
      let dot = frac(0n);
      for (let k = 0; k < m; k++) {
        dot = fAdd(dot, fMul(frac(rows[i]![k]!), gs[j]![k]!));
      }
      const coeff = fDiv(dot, normsq[j]!);
      mu[i]!.push(coeff);
      cur = cur.map((x, k) => fSub(x, fMul(coeff, gs[j]![k]!)));
    }
    gs.push(cur);
    let nn = frac(0n);
    for (let k = 0; k < m; k++) {
      nn = fAdd(nn, fMul(cur[k]!, cur[k]!));
    }
    normsq.push(nn);
  }
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < i; j++) {
      if (fCmp(fAbs(mu[i]![j]!), ETA) > 0) {
        return false;
      }
    }
  }
  for (let i = 1; i < n; i++) {
    const lhs = fMul(DELTA, normsq[i - 1]!);
    const rhs = fAdd(normsq[i]!, fMul(fMul(mu[i]![i - 1]!, mu[i]![i - 1]!), normsq[i - 1]!));
    if (fCmp(lhs, rhs) > 0) {
      return false;
    }
  }
  return true;
}

function hnfRows(rows: readonly bigint[][], ncols: bigint): bigint[][] {
  if (rows.length === 0) {
    return [];
  }
  const M = IntegerMatrixFromEntries(rows as bigint[][]);
  const H = hermite_normal_form(M, 'default', undefined, false) as Parameters<typeof matRows>[0];
  return matRows(H);
}

function latticeFingerprint(rows: readonly bigint[][], ncols: bigint): string {
  const hnf = hnfRows(rows, ncols);
  let first = 0n;
  if (rows.length > 0) {
    for (const x of rows[0]!) {
      first += x * x;
    }
  }
  return [
    `rows=${rows.length}`,
    `hnf=${fmtMat(hnf)}`,
    `reduced=${fmtBool(isLLLReducedExact(rows))}`,
    `first_norm2=${first}`,
  ].join(' ');
}

// ---------------------------------------------------------------------------
// modules/free_module_integer.ts + modules/bkz.ts
// ---------------------------------------------------------------------------

function lattice_rank_degree(flat: bigint[], ncols: bigint): string {
  const L = IntegerLattice(reshape(flat, ncols));
  return `rank=${L.rank()} degree=${L.degree()}`;
}

function lattice_basis(flat: bigint[], ncols: bigint): string {
  const L = IntegerLattice(reshape(flat, ncols));
  return fmtMat(matRows(L.basisMatrix()));
}

/**
 * Canonical fingerprint of the constructor's LLL-reduced basis.
 *
 * Unlike `lattice_basis` this does not pin the reduced basis itself: an
 * LLL-reduced basis is not unique (`v` and `-v`, or `v` and `v - b_0` when
 * `mu` is exactly 1/2, are equally valid), so only the lattice it spans (its
 * HNF), the reducedness predicate and the first squared norm are compared.
 */
function lattice_lll_fingerprint(flat: bigint[], ncols: bigint): string {
  const L = IntegerLattice(reshape(flat, ncols));
  return latticeFingerprint(matRows(L.basisMatrix()), ncols);
}

function lll_exact(flat: bigint[], ncols: bigint): string {
  const L = IntegerLattice(reshape(flat, ncols), { lllReduce: false });
  return fmtMat(matRows(L.LLL()));
}

function lll_invariants(flat: bigint[], ncols: bigint): string {
  const L = IntegerLattice(reshape(flat, ncols), { lllReduce: false });
  return latticeFingerprint(matRows(L.LLL()), ncols);
}

function bkz_exact(flat: bigint[], ncols: bigint, blockSize: bigint): string {
  const L = IntegerLattice(reshape(flat, ncols), { lllReduce: false });
  return fmtMat(matRows(L.BKZ({ blockSize: Number(blockSize) })));
}

function bkz_invariants(flat: bigint[], ncols: bigint, blockSize: bigint): string {
  const L = IntegerLattice(reshape(flat, ncols), { lllReduce: false });
  return latticeFingerprint(matRows(L.BKZ({ blockSize: Number(blockSize) })), ncols);
}

function hkz_first_norm2(flat: bigint[], ncols: bigint): string {
  const L = IntegerLattice(reshape(flat, ncols), { lllReduce: false });
  const R = matRows(L.HKZ());
  let s = 0n;
  for (const x of R[0]!) {
    s += x * x;
  }
  return String(s);
}

function shortest_vector_norm2(flat: bigint[], ncols: bigint): string {
  const L = IntegerLattice(reshape(flat, ncols));
  const v = L.shortestVector();
  let s = 0n;
  for (const x of v) {
    s += x * x;
  }
  return String(s);
}

function closest_vector_dist2(flat: bigint[], ncols: bigint, target: bigint[]): string {
  const L = IntegerLattice(reshape(flat, ncols));
  const v = L.closestVector(target);
  let s = 0n;
  for (let i = 0; i < v.length; i++) {
    const d = v[i]! - target[i]!;
    s += d * d;
  }
  return String(s);
}

function lattice_volume(flat: bigint[], ncols: bigint): string {
  return String(IntegerLattice(reshape(flat, ncols)).volume());
}

function lattice_discriminant(flat: bigint[], ncols: bigint): string {
  return String(IntegerLattice(reshape(flat, ncols)).discriminant());
}

function lattice_is_unimodular(flat: bigint[], ncols: bigint): string {
  return fmtBool(IntegerLattice(reshape(flat, ncols)).isUnimodular());
}

export const functions = {
  // free_module.ts
  zz_span,
  zz_span_qq_entries,
  zz_span_of_basis,
  qq_span,
  qq_span_qq_entries,
  gf_span,
  zz_intersection,
  qq_intersection,
  gf_intersection,
  zz_sum,
  zz_discriminant,
  zz_saturation,
  zz_index_in,
  zz_index_in_ambient,
  zz_coordinates,
  zz_coordinates_user_basis,
  zz_is_submodule,
  zz_module_eq,
  zz_quotient_invariants,
  qq_complement,
  gf_complement,
  gf_cardinality,
  // free_module_element.ts
  vector_ops_zz,
  vector_cross_zz,
  vector_ops_gf,
  vector_scalar_zz,
  // free_module_integer.ts / bkz.ts
  lattice_rank_degree,
  lattice_basis,
  lattice_lll_fingerprint,
  lll_exact,
  lll_invariants,
  bkz_exact,
  bkz_invariants,
  hkz_first_norm2,
  shortest_vector_norm2,
  closest_vector_dist2,
  lattice_volume,
  lattice_discriminant,
  lattice_is_unimodular,
};
