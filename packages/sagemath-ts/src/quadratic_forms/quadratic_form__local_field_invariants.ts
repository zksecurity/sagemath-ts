/**
 * @module sage/quadratic_forms/quadratic_form__local_field_invariants
 * @description Local (`p`-adic and real) invariants of rational quadratic forms.
 *
 * Port of: `sage/quadratic_forms/quadratic_form__local_field_invariants.py`
 * Reference: `reference/sage/src/sage/quadratic_forms/quadratic_form__local_field_invariants.py`
 *
 * SageMath imports these functions into the `QuadraticForm` class body; here
 * they are plain functions taking the form as their first argument, and
 * `QuadraticForm` forwards to them.
 */

import { hilbert_symbol, kronecker_symbol, prime_divisors } from '../arith/misc.js';
import { NotImplementedError, ValueError } from '../errors.js';
import { ZZ } from '../rings/integer_ring.js';
import { Rational } from '../rings/rational.js';
import { QQ } from '../rings/rational_field.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';
import {
  QuadraticForm,
  type RationalMatrix,
  determinantQQ,
  identityQQ,
  inverseQQ,
  matrixQQ,
} from './quadratic_form.js';

/* ------------------------------------------------------------------ */
/* Rational helpers                                                    */
/* ------------------------------------------------------------------ */

/**
 * `hilbert_symbol(a, b, p)` for rational `a`, `b`.
 *
 * SageMath's `sage.arith.misc.hilbert_symbol` starts with
 * `a = QQ(a).numerator() * QQ(a).denominator()` (arith/misc.py:4985); our
 * `hilbert_symbol` only takes integers, so we do the same normalisation here.
 *
 * @see Reference: sage/arith/misc.py:4985
 */
export function hilbert_symbol_QQ(a: Rational, b: Rational, p: IntegerLike): bigint {
  const _p = toBigInt(p);
  const ai = a.numerator * a.denominator;
  const bi = b.numerator * b.denominator;
  return hilbert_symbol(ai, bi, _p);
}

/**
 * Whether the rational `x` is a square in `Q_p`.
 *
 * @see Reference: sage/rings/rational.pyx:1751 (`Rational.is_padic_square`)
 */
export function is_padic_square(x: Rational, p: IntegerLike): boolean {
  if (x.isZero()) {
    return true;
  }
  const _p = toBigInt(p);
  const [eRaw, m] = x.val_unit(_p);
  if (eRaw === 'Infinity') {
    return true;
  }
  const e = eRaw as bigint;
  if (((e % 2n) + 2n) % 2n !== 0n) {
    return false;
  }
  // `m` is a p-adic unit; SageMath evaluates `m % 8` resp. `kronecker(m, p)`
  // in QQ.  For a unit num/den, den is invertible modulo p (resp. 8) and
  // den^{-1} == den there for p = 2, so num*den has the same class as m.
  const mi = m.numerator * m.denominator;
  if (_p === 2n) {
    return ((mi % 8n) + 8n) % 8n === 1n;
  }
  return kronecker_symbol(mi, _p) === 1n;
}

/* ------------------------------------------------------------------ */
/* PARI's qfgaussred (Gauss reduction of a symmetric matrix)           */
/* ------------------------------------------------------------------ */

/**
 * PARI's `qfgaussred` for an exact (rational) symmetric matrix.
 *
 * `_rational_diagonal_form_and_transformation` calls
 * `self.__pari__().qfgaussred()`, so this is a delegation point; `parigp-ts`
 * does not export `qfgaussred` (only `qfgaussred_positive`), so the port of
 * `gaussred` lives here for now.
 *
 * The `suitable()` pivot rule of `alglin2.c:1630-1644` reduces, for exact
 * types (`t_INT`/`t_FRAC`), to "nonzero".
 *
 * @see Reference: reference/pari/src/basemath/alglin2.c:1650-1746 (`gaussred`)
 * @see Deviation: SageMath delegates to PARI; `parigp-ts` has no `qfgaussred`
 * yet, so this is a local port of the same C routine.  It should move to
 * `packages/parigp-ts/src/qfgaussred.ts` when that file exists.
 */
export function qfgaussred(a: RationalMatrix): RationalMatrix {
  const n = a.nrows;
  if (n !== a.ncols) {
    throw new ValueError('gaussred: matrix must be square');
  }
  // Work on a mutable copy (0-indexed; PARI is 1-indexed).
  const A: Rational[][] = [];
  for (let i = 0; i < n; i++) {
    A.push(a.row(i).slice());
  }
  const r = new Array<boolean>(n).fill(true); // r[k] = "row k still active"
  let t = n;

  while (t > 0) {
    let k = 0;
    while (k < n && !(r[k] && !A[k]![k]!.isZero())) {
      k++;
    }
    if (k < n) {
      const p = A[k]![k]!;
      const invp = p.inv();
      r[k] = false;
      t--;
      const ak = A[k]!.slice(); // row(a, k) is taken *before* the scaling
      for (let i = 0; i < n; i++) {
        A[k]![i] = r[i] ? A[k]![i]!.mul(invp) : Rational.zero();
      }
      for (let i = 0; i < n; i++) {
        if (!r[i]) continue;
        const c = ak[i]!;
        if (c.isZero()) continue;
        for (let j = 0; j < n; j++) {
          if (!r[j]) continue;
          A[i]![j] = A[i]![j]!.sub(c.mul(A[k]![j]!));
        }
      }
      A[k]![k] = p;
    } else {
      // All remaining diagonal coefficients are zero.
      let progressed = false;
      for (k = 0; k < n; k++) {
        if (!r[k]) continue;
        let l = k + 1;
        while (l < n && !(r[l] && !A[k]![l]!.isZero())) {
          l++;
        }
        if (l >= n) continue;

        const p = A[k]![l]!;
        const invp = p.inv();
        r[k] = false;
        r[l] = false;
        t -= 2;
        const ak = A[k]!.slice();
        const al = A[l]!.slice();
        for (let i = 0; i < n; i++) {
          if (r[i]) {
            A[k]![i] = A[k]![i]!.mul(invp);
            A[l]![i] = A[l]![i]!.mul(invp);
          } else {
            A[k]![i] = Rational.zero();
            A[l]![i] = Rational.zero();
          }
        }
        for (let i = 0; i < n; i++) {
          if (!r[i]) continue;
          const c = ak[i]!;
          const d = al[i]!;
          for (let j = 0; j < n; j++) {
            if (!r[j]) continue;
            A[i]![j] = A[i]![j]!.sub(A[l]![j]!.mul(c).add(A[k]![j]!.mul(d)));
          }
        }
        for (let i = 0; i < n; i++) {
          if (!r[i]) continue;
          const c = A[k]![i]!;
          const d = A[l]![i]!;
          A[k]![i] = c.add(d);
          A[l]![i] = c.sub(d);
        }
        A[k]![l] = Rational.one();
        A[l]![k] = Rational.one().neg();
        A[k]![k] = p.div(2n);
        A[l]![l] = A[k]![k]!.neg();
        progressed = true;
        break;
      }
      if (!progressed) break;
    }
  }

  const flat: Rational[] = [];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      flat.push(A[i]![j]!);
    }
  }
  return matrixQQ(n, n, flat);
}

/* ------------------------------------------------------------------ */
/* Rational diagonalisation                                            */
/* ------------------------------------------------------------------ */

/**
 * Return `(D, T)` with `T^t * self.matrix() * T == D.matrix()`, both over the
 * fraction field of the base ring.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:186
 */
export function _rational_diagonal_form_and_transformation(
  self: QuadraticForm
): [QuadraticForm, RationalMatrix] {
  const n = self.n;
  // K = self.base_ring().fraction_field(); for ZZ and QQ this is QQ.
  let Q = self.change_ring(QQ);

  // --- PARI path (SageMath: pariself.qfgaussred()) -----------------
  const R = qfgaussred(self.matrix());
  const Dentries: Rational[] = [];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      Dentries.push(i === j ? R.get(i, i) : Rational.zero());
    }
  }
  const D = matrixQQ(n, n, Dentries);
  // Q.parent()(D): the quadratic form whose Hessian matrix is D.
  const newQ = new QuadraticForm(QQ, D);
  const Tentries: Rational[] = [];
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      Tentries.push(i === j ? Rational.one() : R.get(i, j));
    }
  }
  const T0 = matrixQQ(n, n, Tentries);
  try {
    return [newQ, inverseQQ(T0)];
  } catch {
    // Singular case is not fully supported by PARI -- fall through.
  }

  // --- General case (SageMath's own elimination) --------------------
  let T = identityQQ(n);
  for (let i = 0; i < n; i++) {
    if (Q.get(i, i).isZero()) {
      for (let j = i + 1; j < n; j++) {
        if (!Q.get(i, j).isZero()) {
          const temp = identityQQ(n);
          if (Q.get(i, j).add(Q.get(j, j)).isZero()) {
            temp.set(j, i, Rational.one().neg());
          } else {
            temp.set(j, i, Rational.one());
          }
          Q = Q.__call__(temp) as QuadraticForm;
          T = T.mul(temp);
          break;
        }
      }
    }

    const temp = identityQQ(n);
    for (let j = i + 1; j < n; j++) {
      if (!Q.get(i, j).isZero()) {
        temp.set(i, j, Q.get(i, j).neg().div(Q.get(i, i).mul(2n)));
      }
    }
    Q = Q.__call__(temp) as QuadraticForm;
    T = T.mul(temp);
  }

  return [Q, T];
}

/**
 * A diagonal form equivalent to `self` over the fraction field.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:32
 */
export function rational_diagonal_form(self: QuadraticForm): QuadraticForm {
  return self.rational_diagonal_form();
}

/* ------------------------------------------------------------------ */
/* Signature                                                           */
/* ------------------------------------------------------------------ */

/**
 * The triple `(p, n, z)` of the numbers of positive, negative and zero
 * eigenvalues of the associated symmetric matrix.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:297
 */
export function signature_vector(self: QuadraticForm): [bigint, bigint, bigint] {
  const diag = self.rational_diagonal_form();
  let p = 0n;
  let n = 0n;
  let z = 0n;
  for (let i = 0; i < diag.n; i++) {
    const d = diag.get(i, i);
    if (d.isPositive()) {
      p += 1n;
    } else if (d.isNegative()) {
      n += 1n;
    } else {
      z += 1n;
    }
  }
  return [p, n, z];
}

/**
 * `(# positive eigenvalues) - (# negative eigenvalues)`.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:346
 */
export function signature(self: QuadraticForm): bigint {
  const [p, n] = signature_vector(self);
  return p - n;
}

/* ------------------------------------------------------------------ */
/* Hasse invariants                                                    */
/* ------------------------------------------------------------------ */

/**
 * Cassels's Hasse invariant `c_p = prod_{i<j} (a_i, a_j)_p`.
 *
 * @param p a prime, or `-1` for the infinite place
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:381
 */
export function hasse_invariant(self: QuadraticForm, p: IntegerLike): bigint {
  const Diag = self.rational_diagonal_form();
  let hasse_temp = 1n;
  const n = Diag.n;
  for (let j = 0; j < n - 1; j++) {
    for (let k = j + 1; k < n; k++) {
      hasse_temp *= hilbert_symbol_QQ(Diag.get(j, j), Diag.get(k, k), p);
    }
  }
  return hasse_temp;
}

/**
 * O'Meara's Hasse invariant `c_p = prod_{i<=j} (a_i, a_j)_p`.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:472
 */
export function hasse_invariant__OMeara(self: QuadraticForm, p: IntegerLike): bigint {
  const Diag = self.rational_diagonal_form();
  let hasse_temp = 1n;
  const n = Diag.n;
  for (let j = 0; j < n; j++) {
    for (let k = j; k < n; k++) {
      hasse_temp *= hilbert_symbol_QQ(Diag.get(j, j), Diag.get(k, k), p);
    }
  }
  return hasse_temp;
}

/* ------------------------------------------------------------------ */
/* Hyperbolicity / (an)isotropy                                        */
/* ------------------------------------------------------------------ */

/** Binomial coefficient `binomial(m, 2)`. */
function binomial2(m: bigint): bigint {
  return (m * (m - 1n)) / 2n;
}

/**
 * Whether the form is a sum of hyperbolic planes over `Q_p` (or `R` for
 * `p = -1`).
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:559
 */
export function is_hyperbolic(self: QuadraticForm, p: IntegerLike): boolean {
  const _p = toBigInt(p);
  if (self.n % 2 !== 0) {
    return false;
  }
  if (self.n === 0) {
    return true;
  }
  const m = BigInt(self.n / 2);
  if (_p === -1n) {
    return signature(self) === 0n;
  }
  // det * (-1)^m
  const detPart = self.det().mul(m % 2n === 0n ? 1n : -1n);
  if (_p === 2n) {
    // (-1)^binomial(m, 2); here -1 is hilbert_symbol(-1, -1, 2)
    const want = binomial2(m) % 2n === 0n ? 1n : -1n;
    return is_padic_square(detPart, _p) && hasse_invariant(self, _p) === want;
  }
  return is_padic_square(detPart, _p) && hasse_invariant(self, _p) === 1n;
}

/**
 * Whether the form is anisotropic over `Q_p` (or `R` for `p = -1`).
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:619
 */
export function is_anisotropic(self: QuadraticForm, p: IntegerLike): boolean {
  const _p = toBigInt(p);
  if (_p === -1n) {
    return self.is_definite();
  }
  const n = self.n;
  const D = self.det();

  if (n >= 5) {
    return false;
  }
  if (n === 4) {
    return is_padic_square(D, _p) && hasse_invariant(self, _p) === -hilbert_symbol(-1n, -1n, _p);
  }
  if (n === 3) {
    return hasse_invariant(self, _p) !== hilbert_symbol_QQ(new Rational(-1n, 1n), D.neg(), _p);
  }
  if (n === 2) {
    return !is_padic_square(D.neg(), _p);
  }
  if (n === 1) {
    return !self.get(0, 0).isZero();
  }
  throw new NotImplementedError("we have not established a convention for 0-dim'l quadratic forms");
}

/**
 * Whether the form is isotropic over `Q_p` (or `R` for `p = -1`).
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:685
 */
export function is_isotropic(self: QuadraticForm, p: IntegerLike): boolean {
  return !is_anisotropic(self, p);
}

/**
 * All anisotropic primes of the form; `-1` denotes the infinite place.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:729
 */
export function anisotropic_primes(self: QuadraticForm): bigint[] {
  const twoDet = self.det().mul(2n);
  // SageMath: prime_divisors(2 * self.det()).  For a rational argument that is
  // the support of the numerator and the denominator, which is exactly
  // `Rational.support()`; for an integer the two agree.
  const possible_primes: bigint[] = [
    ...(twoDet.isInteger() ? prime_divisors(twoDet.numerator) : twoDet.support()),
    -1n,
  ];
  return possible_primes.filter((p) => is_anisotropic(self, p));
}

/* ------------------------------------------------------------------ */
/* Definiteness                                                        */
/* ------------------------------------------------------------------ */

/**
 * The definiteness string: one of `'pos_def'`, `'neg_def'`, `'indefinite'`,
 * `'zero'`, `'degenerate'`.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:758
 */
export function compute_definiteness_string(self: QuadraticForm): string {
  if (self.base_ring() !== ZZ && self.base_ring() !== QQ) {
    throw new NotImplementedError('we can only check definiteness over ZZ, QQ, and RR for now');
  }
  const n = self.n;
  if (n === 0) {
    return 'zero';
  }
  const [sig_pos, sig_neg, sig_zer] = signature_vector(self);
  if (sig_zer > 0n) {
    return 'degenerate';
  }
  if (sig_neg === BigInt(n)) {
    return 'neg_def';
  }
  if (sig_pos === BigInt(n)) {
    return 'pos_def';
  }
  return 'indefinite';
}

/**
 * Compute and cache the definiteness of the form.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:758
 */
export function compute_definiteness(self: QuadraticForm): void {
  self.compute_definiteness();
}

/** The sign of a rational, as `-1n`, `0n` or `1n`. */
function sgn(x: Rational): bigint {
  if (x.isPositive()) return 1n;
  if (x.isNegative()) return -1n;
  return 0n;
}

/**
 * The definiteness of the form, determined from the signs of its leading
 * principal minors.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:829
 */
export function compute_definiteness_string_by_determinants(self: QuadraticForm): string {
  if (self.base_ring() !== ZZ && self.base_ring() !== QQ) {
    throw new NotImplementedError('we can only check definiteness over ZZ, QQ, and RR for now');
  }
  const n = self.n;
  const M = self.matrix();
  if (n === 0) {
    return 'zero';
  }
  if (self.det().isZero()) {
    return 'degenerate';
  }
  const first_coeff = self.get(0, 0);
  for (let r = 1; r <= n; r++) {
    const entries: Rational[] = [];
    for (let i = 0; i < r; i++) {
      for (let j = 0; j < r; j++) {
        entries.push(M.get(i, j));
      }
    }
    const new_det = determinantQQ(matrixQQ(r, r, entries));
    if (new_det.isZero()) {
      return 'indefinite';
    }
    // sgn(first_coeff)^r
    const s = sgn(first_coeff);
    const sr = s === 0n ? 0n : s === 1n ? 1n : r % 2 === 0 ? 1n : -1n;
    if (sr !== sgn(new_det)) {
      return 'indefinite';
    }
  }
  return first_coeff.isPositive() ? 'pos_def' : 'neg_def';
}

/* ------------------------------------------------------------------ */
/* Definiteness predicates                                             */
/* ------------------------------------------------------------------ */

/** @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:924 */
export function is_positive_definite(self: QuadraticForm): boolean {
  return self.is_positive_definite();
}

/** @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:957 */
export function is_negative_definite(self: QuadraticForm): boolean {
  return self.is_negative_definite();
}

/** @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:991 */
export function is_indefinite(self: QuadraticForm): boolean {
  return self.is_indefinite();
}

/** @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:1025 */
export function is_definite(self: QuadraticForm): boolean {
  return self.is_definite();
}
