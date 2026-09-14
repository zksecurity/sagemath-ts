import { pari_lll_heuristic, pari_lll_proved, pari_lll_real_scalar } from '../pari_lll_real.js';
import { pari_lll_dpe, pari_lll_dpe_resource } from '../pari_lll_dpe.js';
import { pari_lll_fast } from '../pari_lll_fast.js';
import {pari_ideal_matrix} from '../pari_ideal_matrix.js';
/**
 * sagemath-ts side of the `matrix_extended` property-test area.
 *
 * Cases: tests/property/cases/matrix_extended.cases.json
 * SageMath counterpart: tests/property/python/areas/matrix_extended.py
 *
 * Covers the modules that an adversarial audit found to be under-tested:
 *
 *   - `src/matrix/matrix_operations.ts`      (charpoly, minpoly, right_kernel_matrix, ...)
 *   - `src/matrix/matrix_decompositions.ts`  (rref, pivots, hessenberg, smith, jordan, ...)
 *   - `src/matrix/matrix_integer.ts`         (HNF, SNF, LLL, saturation, Frobenius, ...)
 *   - `src/matrix/matrix_modn.ts`            (Matrix_modn_dense)
 *   - `src/matrix/matrix_mod2.ts`            (Matrix_mod2_dense)
 *
 * Historically wrong here, and therefore deliberately targeted below:
 *   * `charpoly` had a sign error         -> `gf_charpoly`, `qq_charpoly`, `zz_charpoly`
 *   * `right_kernel_matrix` echelonized the *transpose* (so it only ever agreed
 *     when the pivot columns happened to be `0..r-1`) -> every `*_rkm*` case uses
 *     matrices whose pivot columns have gaps
 *   * `minpoly` returned the minimal polynomial of `e0` -> `gf_minpoly`,
 *     `qq_minpoly` on matrices with repeated eigenvalues / non-cyclic vectors
 *   * `smith_form` looped                  -> `zz_snf`, `zz_elementary_divisors`
 *
 * Every expected value is produced by running real SageMath (`sage`); nothing in
 * here is hand-derived.  Both sides format their own results into plain strings
 * so that the runners' generic formatters never get a chance to disagree.
 */

import {
  elementary_divisors,
  hermite_form,
  hessenberg_form,
  jordan_form,
  pivots,
  rref,
  smith_form,
} from '../../../../packages/sagemath-ts/src/matrix/matrix_decompositions.js';
import {
  IntegerMatrix,
  LLL,
  elementary_divisors_integer,
  frobenius_form_integer,
  gcd_integer_matrix,
  height,
  hermite_normal_form,
  index_in_saturation,
  integer_valued_polynomials_generators,
  is_LLL_reduced,
  is_primitive,
  kernel_matrix,
  left_kernel_matrix as left_kernel_matrix_integer,
  null_ideal,
  p_minimal_polynomials,
  pivots_integer,
  rank_integer,
  saturation,
  smith_form_integer,
} from '../../../../packages/sagemath-ts/src/matrix/matrix_integer.js';
import { Matrix_mod2_dense } from '../../../../packages/sagemath-ts/src/matrix/matrix_mod2.js';
import { Matrix_modn_dense } from '../../../../packages/sagemath-ts/src/matrix/matrix_modn.js';
import {
  adjugate,
  change_ring,
  charpoly,
  determinant,
  eigenvalues,
  inverse,
  is_diagonalizable,
  left_kernel_matrix,
  minors,
  minpoly,
  permanent,
  rank,
  right_kernel_matrix,
  right_nullity,
  solve_right,
} from '../../../../packages/sagemath-ts/src/matrix/matrix_operations.js';
import { MatrixSpace } from '../../../../packages/sagemath-ts/src/matrix/matrix_space.js';
import { GF } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_constructor.js';
import { Zmod } from '../../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod_ring.js';
import { QQ } from '../../../../packages/sagemath-ts/src/rings/rational_field.js';

// ---------------------------------------------------------------------------
// Formatting helpers (area-private; mirrored exactly in the Python module)
// ---------------------------------------------------------------------------

/** `[[1, 2], [3, 4]]`; a matrix with no rows renders as `[]`. */
// biome-ignore lint/suspicious/noExplicitAny: uniform over four unrelated matrix classes
function fmtMat(M: any): string {
  const rows: string[] = [];
  for (let i = 0; i < M.nrows; i++) {
    const row: string[] = [];
    for (let j = 0; j < M.ncols; j++) {
      row.push(String(M.get(i, j)));
    }
    rows.push(`[${row.join(', ')}]`);
  }
  return `[${rows.join(', ')}]`;
}

/** `[1, 2, 3]` */
function fmtList(xs: readonly unknown[]): string {
  return `[${xs.map((x) => String(x)).join(', ')}]`;
}

/** Coefficient list of a polynomial, constant term first (`Polynomial.list()`). */
// biome-ignore lint/suspicious/noExplicitAny: Polynomial<R> for several R
function fmtPoly(p: any): string {
  const d: number = p.degree();
  if (d < 0) {
    return '[]';
  }
  const coeffs: string[] = [];
  for (let i = 0; i <= d; i++) {
    coeffs.push(String(p.getCoeff(i)));
  }
  return `[${coeffs.join(', ')}]`;
}

/** `True` / `False`, matching Python's `str(bool)`. */
function fmtBool(b: boolean): string {
  return b ? 'True' : 'False';
}

/**
 * Run `f`, returning `ErrorName: message` instead of propagating.
 *
 * Used only where SageMath's own exception text has been checked to be
 * byte-identical, so that "both sides refuse this input, for the same reason"
 * is a real assertion rather than the harness' vacuous both-errored pass.
 */
function guard(f: () => string): string {
  try {
    return f();
  } catch (e) {
    const err = e as Error;
    return `${err.name}: ${err.message}`;
  }
}

// ---------------------------------------------------------------------------
// Matrix builders
// ---------------------------------------------------------------------------

function reshape(entries: bigint[], nrows: number, ncols: number): bigint[][] {
  if (entries.length !== nrows * ncols) {
    throw new Error(`expected ${nrows * ncols} entries, got ${entries.length}`);
  }
  const rows: bigint[][] = [];
  for (let i = 0; i < nrows; i++) {
    rows.push(entries.slice(i * ncols, (i + 1) * ncols));
  }
  return rows;
}

// biome-ignore lint/suspicious/noExplicitAny: the ring/matrix generics are erased here on purpose
function gfMatrix(p: bigint, nrows: number, ncols: number, entries: bigint[]): any {
  const F: any = GF(p);
  const grid = reshape(entries, nrows, ncols).map((row) => row.map((v) => F.__call__(v)));
  return MatrixSpace(F, nrows, ncols).__call__(grid);
}

// biome-ignore lint/suspicious/noExplicitAny: see above
function qqMatrix(nrows: number, ncols: number, entries: bigint[]): any {
  const R: any = QQ;
  const grid = reshape(entries, nrows, ncols).map((row) => row.map((v) => R.__call__(v)));
  return MatrixSpace(R, nrows, ncols).__call__(grid);
}

function zzMatrix(nrows: number, ncols: number, entries: bigint[]): IntegerMatrix {
  return new IntegerMatrix(nrows, ncols, reshape(entries, nrows, ncols));
}

function modnMatrix(n: bigint, nrows: number, ncols: number, entries: bigint[]): Matrix_modn_dense {
  return new Matrix_modn_dense(nrows, ncols, n, reshape(entries, nrows, ncols));
}

function mod2Matrix(nrows: number, ncols: number, entries: bigint[]): Matrix_mod2_dense {
  const grid = reshape(entries, nrows, ncols).map((row) =>
    row.map((v) => Number(((v % 2n) + 2n) % 2n))
  );
  return new Matrix_mod2_dense(nrows, ncols, grid);
}

/**
 * Force a *degenerate* shape onto an otherwise random matrix.
 *
 * The point of the exercise: a uniformly random matrix over a small field
 * almost always has pivot columns `0, 1, ..., r-1`, which is exactly the shape
 * under which the old (transpose-echelonizing) `right_kernel_matrix` accidentally
 * agreed with SageMath.  This zeroes column 0, makes column 2 a multiple of
 * column 1 and zeroes the last row, guaranteeing pivot gaps and a rank deficit
 * no matter what the random entries were.
 */
function degenerate(rows: bigint[][], ncols: number): bigint[][] {
  const out = rows.map((r) => r.slice());
  for (const row of out) {
    if (ncols > 0) row[0] = 0n;
    if (ncols > 2) row[2] = 3n * row[1]!;
  }
  if (out.length > 0) {
    out[out.length - 1] = new Array(ncols).fill(0n);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Dispatch table
// ---------------------------------------------------------------------------

export const functions = {
  pari_lll_heuristic,
  pari_lll_proved,
  pari_lll_real_scalar,
  pari_lll_dpe,
  pari_lll_dpe_resource,
  pari_lll_fast,
  pari_ideal_matrix,
  gf_change_ring_modular: (m: bigint, p: bigint, a: bigint) => {
    const A = MatrixSpace(Zmod(m), 2, 2).__call__([
      [a, a + 1n],
      [a + 2n, a + 3n],
    ]);
    return fmtMat(change_ring(A, GF(p)));
  },
  // === matrix_operations: charpoly (the sign-error regression) ==============
  gf_charpoly: (p: bigint, n: bigint, entries: bigint[]) =>
    fmtPoly(charpoly(gfMatrix(p, Number(n), Number(n), entries))),
  /**
   * The port's second charpoly algorithm (`_charpoly_hessenberg`).
   *
   * SageMath's specialised classes reject `algorithm='hessenberg'`, so the
   * oracle is plain `A.charpoly()` -- the characteristic polynomial is unique,
   * so both of the port's algorithms must reproduce it byte for byte.
   */
  gf_charpoly_hessenberg: (p: bigint, n: bigint, entries: bigint[]) =>
    fmtPoly(charpoly(gfMatrix(p, Number(n), Number(n), entries), 'x', 'hessenberg')),
  qq_charpoly: (n: bigint, entries: bigint[]) =>
    fmtPoly(charpoly(qqMatrix(Number(n), Number(n), entries))),
  /**
   * The port has no `IntegerMatrix.charpoly`, so the integer matrix is carried
   * into QQ; the oracle is SageMath's integer (FLINT) charpoly, which is a
   * different code path from `matrix(QQ, ...).charpoly()`.
   */
  zz_charpoly: (n: bigint, entries: bigint[]) =>
    fmtPoly(charpoly(qqMatrix(Number(n), Number(n), entries))),

  // === matrix_operations: minpoly (the "minpoly of e0" regression) =========
  gf_minpoly: (p: bigint, n: bigint, entries: bigint[]) =>
    fmtPoly(minpoly(gfMatrix(p, Number(n), Number(n), entries))),
  qq_minpoly: (n: bigint, entries: bigint[]) =>
    fmtPoly(minpoly(qqMatrix(Number(n), Number(n), entries))),

  // === matrix_operations: kernels (the transpose-echelon regression) =======
  gf_rkm: (p: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    fmtMat(right_kernel_matrix(gfMatrix(p, Number(r), Number(c), entries))),
  gf_rkm_pivot: (p: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    fmtMat(right_kernel_matrix(gfMatrix(p, Number(r), Number(c), entries), { basis: 'pivot' })),
  gf_lkm: (p: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    fmtMat(left_kernel_matrix(gfMatrix(p, Number(r), Number(c), entries))),
  qq_rkm: (r: bigint, c: bigint, entries: bigint[]) =>
    fmtMat(right_kernel_matrix(qqMatrix(Number(r), Number(c), entries))),
  qq_rkm_pivot: (r: bigint, c: bigint, entries: bigint[]) =>
    fmtMat(right_kernel_matrix(qqMatrix(Number(r), Number(c), entries), { basis: 'pivot' })),
  /** Random entries forced into a pivot-gapped, rank-deficient shape. */
  gf_rkm_degenerate: (p: bigint, r: bigint, c: bigint, entries: bigint[]) => {
    const nr = Number(r);
    const nc = Number(c);
    const grid = degenerate(reshape(entries, nr, nc), nc);
    const flat = grid.flat();
    return fmtMat(right_kernel_matrix(gfMatrix(p, nr, nc, flat)));
  },
  gf_rref_degenerate: (p: bigint, r: bigint, c: bigint, entries: bigint[]) => {
    const nr = Number(r);
    const nc = Number(c);
    const grid = degenerate(reshape(entries, nr, nc), nc);
    return fmtMat(rref(gfMatrix(p, nr, nc, grid.flat())));
  },

  // === matrix_operations: misc =============================================
  gf_rank: (p: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    String(rank(gfMatrix(p, Number(r), Number(c), entries))),
  gf_right_nullity: (p: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    String(right_nullity(gfMatrix(p, Number(r), Number(c), entries))),
  gf_det: (p: bigint, n: bigint, entries: bigint[]) =>
    String(determinant(gfMatrix(p, Number(n), Number(n), entries))),
  qq_det: (n: bigint, entries: bigint[]) =>
    String(determinant(qqMatrix(Number(n), Number(n), entries))),
  gf_inverse: (p: bigint, n: bigint, entries: bigint[]) =>
    guard(() => fmtMat(inverse(gfMatrix(p, Number(n), Number(n), entries)))),
  qq_inverse: (n: bigint, entries: bigint[]) =>
    guard(() => fmtMat(inverse(qqMatrix(Number(n), Number(n), entries)))),
  gf_adjugate: (p: bigint, n: bigint, entries: bigint[]) =>
    fmtMat(adjugate(gfMatrix(p, Number(n), Number(n), entries))),
  /**
   * Guarded: SageMath rejects `m > n` with
   * `ValueError: must have m <= n, but m (=..) and n (=..)` (matrix2.pyx:1645).
   * Without the guard both sides would merely raise and `compare.ts` would score
   * the case as a vacuous pass without ever comparing the messages.
   */
  gf_permanent: (p: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    guard(() => String(permanent(gfMatrix(p, Number(r), Number(c), entries)))),
  gf_minors: (p: bigint, r: bigint, c: bigint, k: bigint, entries: bigint[]) =>
    fmtList(minors(gfMatrix(p, Number(r), Number(c), entries), Number(k))),
  /** Eigenvalues in the base field only, sorted so the order is canonical. */
  gf_eigenvalues: (p: bigint, n: bigint, entries: bigint[]) => {
    const vals = eigenvalues(gfMatrix(p, Number(n), Number(n), entries), false) as unknown[];
    const ints = vals.map((v) => BigInt(String(v)));
    ints.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
    return fmtList(ints);
  },
  gf_is_diagonalizable: (p: bigint, n: bigint, entries: bigint[]) =>
    fmtBool(is_diagonalizable(gfMatrix(p, Number(n), Number(n), entries))),
  gf_solve_right: (p: bigint, n: bigint, entriesA: bigint[], entriesB: bigint[]) =>
    guard(() => {
      const A = gfMatrix(p, Number(n), Number(n), entriesA);
      const b = gfMatrix(p, Number(n), 1, entriesB);
      return fmtMat(solve_right(A, b));
    }),

  // === matrix_decompositions ==============================================
  gf_rref: (p: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    fmtMat(rref(gfMatrix(p, Number(r), Number(c), entries))),
  qq_rref: (r: bigint, c: bigint, entries: bigint[]) =>
    fmtMat(rref(qqMatrix(Number(r), Number(c), entries))),
  gf_pivots: (p: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    fmtList(pivots(gfMatrix(p, Number(r), Number(c), entries))),
  gf_hessenberg: (p: bigint, n: bigint, entries: bigint[]) =>
    fmtMat(hessenberg_form(gfMatrix(p, Number(n), Number(n), entries))),
  qq_hessenberg: (n: bigint, entries: bigint[]) =>
    fmtMat(hessenberg_form(qqMatrix(Number(n), Number(n), entries))),
  gf_smith_form: (p: bigint, r: bigint, c: bigint, entries: bigint[]) => {
    const res = smith_form(gfMatrix(p, Number(r), Number(c), entries), true);
    return fmtMat(Array.isArray(res) ? res[0] : res);
  },
  qq_smith_form: (r: bigint, c: bigint, entries: bigint[]) => {
    const res = smith_form(qqMatrix(Number(r), Number(c), entries), true);
    return fmtMat(Array.isArray(res) ? res[0] : res);
  },
  gf_elementary_divisors: (p: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    fmtList(elementary_divisors(gfMatrix(p, Number(r), Number(c), entries))),
  gf_hermite_form: (p: bigint, r: bigint, c: bigint, entries: bigint[]) => {
    const res = hermite_form(gfMatrix(p, Number(r), Number(c), entries));
    return fmtMat(Array.isArray(res) ? res[0] : res);
  },
  gf_jordan_form: (p: bigint, n: bigint, entries: bigint[]) =>
    guard(() => fmtMat(jordan_form(gfMatrix(p, Number(n), Number(n), entries)))),
  qq_jordan_form: (n: bigint, entries: bigint[]) =>
    guard(() => fmtMat(jordan_form(qqMatrix(Number(n), Number(n), entries)))),

  // === matrix_integer ======================================================
  zz_hnf: (r: bigint, c: bigint, entries: bigint[]) =>
    fmtMat(hermite_normal_form(zzMatrix(Number(r), Number(c), entries))),
  zz_snf: (r: bigint, c: bigint, entries: bigint[]) => {
    const res = smith_form_integer(zzMatrix(Number(r), Number(c), entries), true);
    return fmtMat(Array.isArray(res) ? res[0] : res);
  },
  zz_elementary_divisors: (r: bigint, c: bigint, entries: bigint[]) =>
    fmtList(elementary_divisors_integer(zzMatrix(Number(r), Number(c), entries))),
  zz_rank: (r: bigint, c: bigint, entries: bigint[]) =>
    String(rank_integer(zzMatrix(Number(r), Number(c), entries))),
  zz_det: (n: bigint, entries: bigint[]) =>
    String(zzMatrix(Number(n), Number(n), entries).determinant()),
  zz_pivots: (r: bigint, c: bigint, entries: bigint[]) =>
    fmtList(pivots_integer(zzMatrix(Number(r), Number(c), entries))),
  /**
   * `A.right_kernel_matrix(basis='echelon')`.
   *
   * The 'computed' basis over ZZ is whatever FLINT happened to return, so it is
   * not a legitimate differential target; the *echelon* basis is the HNF of the
   * (saturated, hence unique) kernel lattice and is therefore canonical.
   */
  zz_kernel_echelon: (r: bigint, c: bigint, entries: bigint[]) => {
    const K = kernel_matrix(zzMatrix(Number(r), Number(c), entries));
    if (K.nrows === 0) {
      return fmtMat(K);
    }
    return fmtMat(hermite_normal_form(K));
  },
  zz_left_kernel_echelon: (r: bigint, c: bigint, entries: bigint[]) => {
    const K = left_kernel_matrix_integer(zzMatrix(Number(r), Number(c), entries));
    if (K.nrows === 0) {
      return fmtMat(K);
    }
    return fmtMat(hermite_normal_form(K));
  },
  zz_lll: (r: bigint, c: bigint, entries: bigint[]) =>
    fmtMat(LLL(zzMatrix(Number(r), Number(c), entries))),
  zz_is_lll_reduced: (r: bigint, c: bigint, entries: bigint[]) =>
    fmtBool(is_LLL_reduced(zzMatrix(Number(r), Number(c), entries))),
  zz_p_minimal: (n: bigint, p: bigint, sMax: bigint, entries: bigint[]) => {
    const values = p_minimal_polynomials(zzMatrix(Number(n), Number(n), entries), p, Number(sMax));
    return `[${[...values].map(([s, f]) => `[${s}, ${fmtList(f)}]`).join(', ')}]`;
  },
  zz_null_ideal: (n: bigint, b: bigint, entries: bigint[]) => {
    const generators = null_ideal(zzMatrix(Number(n), Number(n), entries), b);
    return `[${generators.map((f) => fmtList(f)).join(', ')}]`;
  },
  zz_integer_valued_polynomials: (n: bigint, entries: bigint[]) => {
    const [mu, generators] = integer_valued_polynomials_generators(
      zzMatrix(Number(n), Number(n), entries)
    );
    return `${fmtList(mu)}|[${generators.map((f) => fmtList(f)).join(', ')}]`;
  },
  zz_height: (r: bigint, c: bigint, entries: bigint[]) =>
    String(height(zzMatrix(Number(r), Number(c), entries))),
  zz_gcd: (r: bigint, c: bigint, entries: bigint[]) =>
    String(gcd_integer_matrix(zzMatrix(Number(r), Number(c), entries))),
  /**
   * SageMath's `Matrix_integer_dense.is_primitive` (matrix_integer_dense.pyx:1145)
   * is *Perron-Frobenius* primitivity: all entries nonnegative and `A^n`
   * entrywise positive for some `n > 0`.  The cases below are taken verbatim
   * from that method's doctests.
   */
  zz_is_primitive: (r: bigint, c: bigint, entries: bigint[]) =>
    fmtBool(is_primitive(zzMatrix(Number(r), Number(c), entries))),
  zz_saturation: (r: bigint, c: bigint, entries: bigint[]) => {
    const S = saturation(zzMatrix(Number(r), Number(c), entries));
    const M = Array.isArray(S) ? S[0] : S;
    return fmtMat(hermite_normal_form(M));
  },
  zz_index_in_saturation: (r: bigint, c: bigint, entries: bigint[]) =>
    String(index_in_saturation(zzMatrix(Number(r), Number(c), entries))),
  zz_frobenius: (n: bigint, entries: bigint[]) =>
    fmtMat(frobenius_form_integer(zzMatrix(Number(n), Number(n), entries)) as IntegerMatrix),

  // === matrix_modn =========================================================
  modn_echelon: (n: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    guard(() => fmtMat(modnMatrix(n, Number(r), Number(c), entries).echelon_form())),
  modn_pivots: (n: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    guard(() => fmtList(modnMatrix(n, Number(r), Number(c), entries).pivots())),
  modn_rank: (n: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    guard(() => String(modnMatrix(n, Number(r), Number(c), entries).rank())),
  modn_det: (n: bigint, k: bigint, entries: bigint[]) =>
    guard(() => String(modnMatrix(n, Number(k), Number(k), entries).determinant())),
  modn_charpoly: (n: bigint, k: bigint, entries: bigint[]) =>
    guard(() => fmtList(modnMatrix(n, Number(k), Number(k), entries).charpoly())),
  modn_minpoly: (n: bigint, k: bigint, entries: bigint[]) =>
    guard(() => fmtList(modnMatrix(n, Number(k), Number(k), entries).minpoly())),
  modn_inverse: (n: bigint, k: bigint, entries: bigint[]) =>
    guard(() => fmtMat(modnMatrix(n, Number(k), Number(k), entries).inverse())),
  modn_rkm: (n: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    guard(() => fmtMat(modnMatrix(n, Number(r), Number(c), entries).right_kernel_matrix())),
  modn_rkm_computed: (n: bigint, r: bigint, c: bigint, entries: bigint[]) =>
    guard(() =>
      fmtMat(
        modnMatrix(n, Number(r), Number(c), entries).right_kernel_matrix({ basis: 'computed' })
      )
    ),

  // === matrix_mod2 =========================================================
  mod2_echelon: (r: bigint, c: bigint, entries: bigint[]) =>
    fmtMat(mod2Matrix(Number(r), Number(c), entries).echelon_form()),
  mod2_pivots: (r: bigint, c: bigint, entries: bigint[]) =>
    fmtList(mod2Matrix(Number(r), Number(c), entries).pivots()),
  mod2_rank: (r: bigint, c: bigint, entries: bigint[]) =>
    String(mod2Matrix(Number(r), Number(c), entries).rank()),
  mod2_det: (n: bigint, entries: bigint[]) =>
    String(mod2Matrix(Number(n), Number(n), entries).determinant()),
  mod2_inverse: (n: bigint, entries: bigint[]) =>
    guard(() => fmtMat(mod2Matrix(Number(n), Number(n), entries).inverse())),
  mod2_rkm: (r: bigint, c: bigint, entries: bigint[]) =>
    fmtMat(mod2Matrix(Number(r), Number(c), entries).right_kernel_matrix()),
  mod2_charpoly: (n: bigint, entries: bigint[]) => {
    const A = mod2Matrix(Number(n), Number(n), entries);
    const flat: bigint[] = [];
    for (let i = 0; i < A.nrows; i++) {
      for (let j = 0; j < A.ncols; j++) flat.push(BigInt(A.get(i, j)));
    }
    return fmtPoly(charpoly(gfMatrix(2n, Number(n), Number(n), flat)));
  },
};

import { pari_qr } from '../pari_qr.js';
functions.pari_qr = pari_qr;

import { pari_real_linear } from '../pari_real_linear.js';
functions.pari_real_linear = pari_real_linear;

import { pari_matrix_products } from '../pari_matrix_products.js';
functions.pari_matrix_products = pari_matrix_products;

import { pari_lll_gso } from '../pari_lll_gso.js';
functions.pari_lll_gso = pari_lll_gso;

import { pari_word_linear } from '../pari_word_linear.js';
import { pari_integer_linear } from '../pari_integer_linear.js';
Object.assign(functions, { pari_word_linear, pari_integer_linear });

import { native_pari_lll_wrapper } from '../pari_lll_wrapper.js';
functions.pari_lll_wrapper = native_pari_lll_wrapper;

import { pari_lll_dependents } from '../pari_lll_dependents.js';
Object.assign(functions, { pari_lll_dependents });

import {pari_prime_matrix}from "../pari_prime_matrices.js";
Object.assign(functions,{pari_prime_matrix});

import {pari_lll_norms}from '../pari_lll_norms.js';
Object.assign(functions,{pari_lll_norms});

import {pari_lll_norms_resource,pari_lll_norms_bounded}from '../pari_lll_norms.js';
Object.assign(functions,{pari_lll_norms_resource,pari_lll_norms_bounded});

import {pari_hnfperm}from "../pari_hnfperm.js";
Object.assign(functions,{pari_hnfperm});

import {pari_lll_progress}from "../pari_lll_progress.js";
Object.assign(functions,{pari_lll_progress});

import {pari_factor_bounds}from "../pari_factor_bounds.js";
Object.assign(functions,{pari_factor_bounds});

import {pari_recombination}from "../pari_recombination.js";
Object.assign(functions,{pari_recombination});

import {pari_general_hnf}from "../pari_general_hnf.js";
Object.assign(functions,{pari_general_hnf});

import {pari_newton_sums}from "../pari_newton_sums.js";
Object.assign(functions,{pari_newton_sums});

import {pari_qx_factor}from "../pari_qx_factor.js";
Object.assign(functions,{pari_qx_factor});

import {ntl_factor_traces} from "../ntl_factor_traces.js";
functions["ntl_factor_traces"]=ntl_factor_traces;

import {ntl_multilift} from "../ntl_multilift.js";
functions["ntl_multilift"]=ntl_multilift;

import {ntl_lll} from "../ntl_lll.js";
functions["ntl_lll"]=ntl_lll;

import {ntl_factor_recovery} from "../ntl_factor_recovery.js";
functions["ntl_factor_recovery"]=ntl_factor_recovery;

import {ntl_factor_selection} from '../ntl_factor_selection.js';
functions.ntl_factor_selection = ntl_factor_selection;

import {ntl_word_quotient} from '../ntl_word_quotient.js';
functions.ntl_word_quotient = ntl_word_quotient;

import {ntl_word_matrix} from '../ntl_word_matrix.js';
functions.ntl_word_matrix = ntl_word_matrix;

import {ntl_word_composition} from '../ntl_word_composition.js';
functions.ntl_word_composition = ntl_word_composition;

import { ntl_word_ddf } from '../ntl_word_ddf.js';
functions.ntl_word_ddf = ntl_word_ddf;

import { ntl_word_minpoly } from '../ntl_word_minpoly.js';
functions.ntl_word_minpoly = ntl_word_minpoly;

import { ntl_prime_sequence } from '../ntl_prime_sequence.js';
functions.ntl_prime_sequence = ntl_prime_sequence;

import { ntl_word_projection } from '../ntl_word_projection.js';
functions.ntl_word_projection = ntl_word_projection;

import {ntl_element_minpoly} from '../ntl_element_minpoly.js';
functions.ntl_element_minpoly = ntl_element_minpoly;

import {ntl_root_prerequisites} from '../ntl_root_prerequisites.js';
functions.ntl_root_prerequisites = ntl_root_prerequisites;

import {ntl_word_factor_recovery} from '../ntl_word_factor_recovery.js';
functions.ntl_word_factor_recovery = ntl_word_factor_recovery;

import {ntl_word_context} from '../ntl_word_context.js';
functions.ntl_word_context = ntl_word_context;

import {ntl_small_prime_selection} from '../ntl_small_prime_selection.js';
functions.ntl_small_prime_selection = ntl_small_prime_selection;

import { ntl_local_info_updates } from '../ntl_local_info_updates.js';
functions.ntl_local_info_updates = ntl_local_info_updates;

import { ntl_cardinality_search } from '../ntl_cardinality_search.js';
functions.ntl_cardinality_search = ntl_cardinality_search;

import { ntl_word_linear } from '../ntl_word_linear.js';
functions.ntl_word_linear = ntl_word_linear;

import { ntl_integer_crt } from '../ntl_integer_crt.js';
functions.ntl_integer_crt = ntl_integer_crt;

import { ntl_integer_reconstruction } from '../ntl_integer_reconstruction.js';
functions.ntl_integer_reconstruction = ntl_integer_reconstruction;

import { ntl_certified_gauss } from '../ntl_certified_gauss.js';
functions.ntl_certified_gauss = ntl_certified_gauss;

import { ntl_polynomial_state } from '../ntl_polynomial_state.js';
functions.ntl_polynomial_state = ntl_polynomial_state;

import { ntl_integer_polynomial_division } from '../ntl_integer_polynomial_division.js';
functions.ntl_integer_polynomial_division = ntl_integer_polynomial_division;

import { ntl_got_them } from '../ntl_got_them.js';
functions.ntl_got_them = ntl_got_them;

import { ntl_cold_prime_selection } from '../ntl_cold_prime_selection.js';
functions['ntl_cold_prime_selection'] = ntl_cold_prime_selection;

import { ntl_cold_recombination } from '../ntl_cold_recombination.js';
functions['ntl_cold_recombination'] = ntl_cold_recombination;

import { ntl_integer_product_state } from '../ntl_integer_product_state.js';
functions['ntl_integer_product_state'] = ntl_integer_product_state;

import { ntl_word_rebuild } from '../ntl_word_rebuild.js';
functions['ntl_word_rebuild'] = ntl_word_rebuild;
functions['ntl_word_projection_rebuild'] = ntl_word_rebuild;
functions['ntl_word_multiplier_fallback'] = ntl_word_rebuild;
import { ntl_big_quotient } from '../ntl_big_quotient.js';
functions['ntl_big_quotient'] = ntl_big_quotient;

import { ntl_hensel_state } from '../ntl_hensel_state.js';
functions['ntl_hensel_state'] = ntl_hensel_state;

import { ntl_word_fft_context } from '../ntl_word_fft_context.js';
functions['ntl_word_fft_context'] = ntl_word_fft_context;

import { ntl_integer_gcd_state } from '../ntl_integer_gcd_state.js';
functions['ntl_integer_gcd_state'] = ntl_integer_gcd_state;

import { ntl_integer_factor_driver } from '../ntl_integer_factor_driver.js';
functions['ntl_integer_factor_driver'] = ntl_integer_factor_driver;
