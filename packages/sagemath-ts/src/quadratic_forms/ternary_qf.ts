/**
 * @module sage/quadratic_forms/ternary_qf
 * @description Ternary quadratic forms with integer coefficients.
 *
 * Port of: `sage/quadratic_forms/ternary_qf.py` together with the helper
 * routines of `sage/quadratic_forms/ternary.pyx` that it calls.
 *
 * Reference: `reference/sage/src/sage/quadratic_forms/ternary_qf.py`
 * Reference: `reference/sage/src/sage/quadratic_forms/ternary.pyx`
 *
 * A `TernaryQF` given by `[a, b, c, r, s, t]` is the form
 * `a x^2 + b y^2 + c z^2 + r y z + s x z + t x y`.
 *
 * Not ported (each throws `NotImplementedError`): the automorphism machinery
 * (`_border`, `_borders`, `_automorphisms_reduced_fast/slow`, `automorphisms`,
 * `number_of_automorphisms`) and the level/discriminant search
 * (`find_all_ternary_qf_by_level_disc`, `find_a_ternary_qf_by_level_disc`).
 */

import { gcd, inverse_mod, kronecker_symbol, sqrt_mod, squarefree_part } from '../arith/misc.js';
import { NotImplementedError, ValueError } from '../errors.js';
import { ZZ } from '../rings/integer_ring.js';
import { Rational } from '../rings/rational.js';
import { QQ } from '../rings/rational_field.js';
import { type IntegerLike, type RationalLike, toBigInt, toRational } from '../types/coercion.js';
import {
  QuadraticForm,
  type RationalMatrix,
  determinantQQ,
  isRationalMatrix,
  matrixQQ,
} from './quadratic_form.js';

/* ------------------------------------------------------------------ */
/* Python integer semantics                                            */
/* ------------------------------------------------------------------ */

/** Python's `//`: floor division. */
function fdiv(a: bigint, b: bigint): bigint {
  const q = a / b;
  if (a % b !== 0n && a < 0n !== b < 0n) {
    return q - 1n;
  }
  return q;
}

/** Python's `%`: the remainder takes the sign of the divisor. */
function pmod(a: bigint, b: bigint): bigint {
  return a - fdiv(a, b) * b;
}

function babs(a: bigint): bigint {
  return a < 0n ? -a : a;
}

/**
 * C's `%` (the remainder takes the sign of the *dividend*).
 *
 * `_find_zeros_mod_p_odd` is a Cython function whose `x0`, `y0`, `i` and `p`
 * are all `long long`, so the `%` in its body is C's remainder, not Python's.
 * That is observable: SageMath's zero list contains entries such as
 * `(4, -1, 1)` rather than `(4, 6, 1)`, and we reproduce them exactly.
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:575 (the `long long` signature)
 */
function cmod(a: bigint, b: bigint): bigint {
  return a % b;
}

/**
 * A rational reduced modulo `n` (Sage's `Rational.__mod__`): `num * den^{-1}`
 * taken in `[0, n)`.  The denominator must be invertible modulo `n`.
 */
function rational_mod(x: Rational, n: bigint): bigint {
  return pmod(x.numerator * inverse_mod(x.denominator, n), n);
}

/** Convert an exactly-integral rational to a `bigint`, or throw. */
function toZZ(x: Rational, what: string): bigint {
  if (!x.isInteger()) {
    throw new ValueError(`${what}: ${x.toString()} is not an integer`);
  }
  return x.numerator;
}

/* ------------------------------------------------------------------ */
/* ternary.pyx helpers                                                 */
/* ------------------------------------------------------------------ */

/**
 * Auxiliary reduction factor of two integers.
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:22 (`red_mfact`)
 */
export function red_mfact(a: bigint, b: bigint): bigint {
  if (a !== 0n) {
    return fdiv(-b + babs(a), 2n * a);
  }
  return 0n;
}

/** The six coefficients `[a, b, c, r, s, t]` of a ternary form. */
export type TernaryCoefficients = [bigint, bigint, bigint, bigint, bigint, bigint];

/**
 * The unique Eisenstein-reduced form equivalent to the given positive definite
 * ternary form, together with the transformation matrix.
 *
 * The coefficient updates are literally those of
 * `_reduced_ternary_form_eisenstein_with_matrix`; the matrix-free variant
 * `_reduced_ternary_form_eisenstein_without_matrix` performs exactly the same
 * updates, so one implementation serves both.
 *
 * The argument order follows `ternary.pyx`: `(a1, a2, a3, a23, a13, a12)`,
 * which is `(a, b, c, r, s, t)`.
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:47
 * @see Reference: sage/quadratic_forms/ternary.pyx:311
 */
export function _reduced_ternary_form_eisenstein_with_matrix(
  a1i: IntegerLike,
  a2i: IntegerLike,
  a3i: IntegerLike,
  a23i: IntegerLike,
  a13i: IntegerLike,
  a12i: IntegerLike
): [TernaryCoefficients, bigint[][]] {
  let a1 = toBigInt(a1i);
  let a2 = toBigInt(a2i);
  let a3 = toBigInt(a3i);
  let a23 = toBigInt(a23i);
  let a13 = toBigInt(a13i);
  let a12 = toBigInt(a12i);

  let m11 = 1n;
  let m12 = 0n;
  let m13 = 0n;
  let m21 = 0n;
  let m22 = 1n;
  let m23 = 0n;
  let m31 = 0n;
  let m32 = 0n;
  let m33 = 1n;

  let loop = true;
  while (loop) {
    // adjust
    const v = a1 + a2 + a23 + a13 + a12;
    if (v < 0n) {
      m13 = m11 + m12 + m13;
      m23 = m21 + m22 + m23;
      m33 = m31 + m32 + m33;
      a3 += v;
      a23 += a12 + 2n * a2;
      a13 += a12 + 2n * a1;
    }

    // cuadred 12
    let m = red_mfact(a1, a12);
    m12 = m * m11 + m12;
    m22 = m * m21 + m22;
    m32 = m * m31 + m32;
    let t = a1 * m;
    a12 += t;
    a2 += a12 * m;
    a12 += t;
    a23 += a13 * m;

    // cuadred 23
    m = red_mfact(a2, a23);
    m13 = m * m12 + m13;
    m23 = m * m22 + m23;
    m33 = m * m32 + m33;
    t = a2 * m;
    a23 += t;
    a3 += a23 * m;
    a23 += t;
    a13 += a12 * m;

    // cuadred 13
    m = red_mfact(a1, a13);
    m13 = m * m11 + m13;
    m23 = m * m21 + m23;
    m33 = m * m31 + m33;
    t = a1 * m;
    a13 += t;
    a3 += a13 * m;
    a13 += t;
    a23 += a12 * m;

    // order 12
    if (a1 > a2 || (a1 === a2 && babs(a23) > babs(a13))) {
      [m11, m12, m13] = [-m12, -m11, -m13];
      [m21, m22, m23] = [-m22, -m21, -m23];
      [m31, m32, m33] = [-m32, -m31, -m33];
      [a1, a2] = [a2, a1];
      [a13, a23] = [a23, a13];
    }

    // order 23
    if (a2 > a3 || (a2 === a3 && babs(a13) > babs(a12))) {
      [m11, m12, m13] = [-m11, -m13, -m12];
      [m21, m22, m23] = [-m21, -m23, -m22];
      [m31, m32, m33] = [-m31, -m33, -m32];
      [a2, a3] = [a3, a2];
      [a13, a12] = [a12, a13];
    }

    // order 12
    if (a1 > a2 || (a1 === a2 && babs(a23) > babs(a13))) {
      [m11, m12, m13] = [-m12, -m11, -m13];
      [m21, m22, m23] = [-m22, -m21, -m23];
      [m31, m32, m33] = [-m32, -m31, -m33];
      [a1, a2] = [a2, a1];
      [a13, a23] = [a23, a13];
    }

    // signs
    if (a23 * a13 * a12 > 0n) {
      if (a23 < 0n) {
        m11 = -m11;
        m21 = -m21;
        m31 = -m31;
        a23 = -a23;
      }
      if (a13 < 0n) {
        m12 = -m12;
        m22 = -m22;
        m32 = -m32;
        a13 = -a13;
      }
      if (a12 < 0n) {
        m13 = -m13;
        m23 = -m23;
        m33 = -m33;
        a12 = -a12;
      }
    } else {
      let s1 = a23 > 0n;
      let s2 = a13 > 0n;
      let s3 = a12 > 0n;
      if (((s1 ? 1 : 0) + (s2 ? 1 : 0) + (s3 ? 1 : 0)) % 2 !== 0) {
        if (a23 === 0n) {
          s1 = true;
        } else if (a13 === 0n) {
          s2 = true;
        } else if (a12 === 0n) {
          s3 = true;
        }
      }
      if (s1) {
        m11 = -m11;
        m21 = -m21;
        m31 = -m31;
        a23 = -a23;
      }
      if (s2) {
        m12 = -m12;
        m22 = -m22;
        m32 = -m32;
        a13 = -a13;
      }
      if (s3) {
        m13 = -m13;
        m23 = -m23;
        m33 = -m33;
        a12 = -a12;
      }
    }

    loop = !(
      babs(a23) <= a2 &&
      babs(a13) <= a1 &&
      babs(a12) <= a1 &&
      a1 + a2 + a23 + a13 + a12 >= 0n
    );
  }

  // adj 3
  if (a1 + a2 + a23 + a13 + a12 === 0n && 2n * a1 + 2n * a13 + a12 > 0n) {
    [m11, m12, m13] = [-m11, -m12, m11 + m12 + m13];
    [m21, m22, m23] = [-m21, -m22, m21 + m22 + m23];
    [m31, m32, m33] = [-m31, -m32, m31 + m32 + m33];
    a23 = -2n * a2 - a23 - a12;
    a13 = -2n * a1 - a13 - a12;
  }

  // adj 5.12
  if (a1 === -a12 && a13 !== 0n) {
    [m11, m12] = [-m11, -m11 - m12];
    [m21, m22] = [-m21, -m21 - m22];
    [m31, m32] = [-m31, -m31 - m32];
    a23 = -a23 - a13;
    a13 = -a13;
    a12 = -a12;
  }

  // adj 5.13
  if (a1 === -a13 && a12 !== 0n) {
    [m11, m13] = [-m11, -m11 - m13];
    [m21, m23] = [-m21, -m21 - m23];
    [m31, m33] = [-m31, -m31 - m33];
    a23 = -a23 - a12;
    a13 = -a13;
    a12 = -a12;
  }

  // adj 5.23
  if (a2 === -a23 && a12 !== 0n) {
    [m12, m13] = [-m12, -m12 - m13];
    [m22, m23] = [-m22, -m22 - m23];
    [m32, m33] = [-m32, -m32 - m33];
    a23 = -a23;
    a13 = -a13 - a12;
    a12 = -a12;
  }

  // adj 4.12
  if (a1 === a12 && a13 > 2n * a23) {
    [m11, m12, m13] = [-m11, -m11 + m12, -m13];
    [m21, m22, m23] = [-m21, -m21 + m22, -m23];
    [m31, m32, m33] = [-m31, -m31 + m32, -m33];
    a23 = -a23 + a13;
  }

  // adj 4.13
  if (a1 === a13 && a12 > 2n * a23) {
    [m11, m12, m13] = [-m11, -m12, -m11 + m13];
    [m21, m22, m23] = [-m21, -m22, -m21 + m23];
    [m31, m32, m33] = [-m31, -m32, -m31 + m33];
    a23 = -a23 + a12;
  }

  // adj 4.23
  if (a2 === a23 && a12 > 2n * a13) {
    [m11, m12, m13] = [-m11, -m12, -m12 + m13];
    [m21, m22, m23] = [-m21, -m22, -m22 + m23];
    [m31, m32, m33] = [-m31, -m32, -m32 + m33];
    a13 = -a13 + a12;
  }

  // order 12
  if (a1 === a2 && babs(a23) > babs(a13)) {
    [m11, m12, m13] = [-m12, -m11, -m13];
    [m21, m22, m23] = [-m22, -m21, -m23];
    [m31, m32, m33] = [-m32, -m31, -m33];
    [a1, a2] = [a2, a1];
    [a13, a23] = [a23, a13];
  }

  // order 23
  if (a2 === a3 && babs(a13) > babs(a12)) {
    [m11, m12, m13] = [-m11, -m13, -m12];
    [m21, m22, m23] = [-m21, -m23, -m22];
    [m31, m32, m33] = [-m31, -m33, -m32];
    [a13, a12] = [a12, a13];
  }

  // order 12
  if (a1 === a2 && babs(a23) > babs(a13)) {
    [m11, m12, m13] = [-m12, -m11, -m13];
    [m21, m22, m23] = [-m22, -m21, -m23];
    [m31, m32, m33] = [-m32, -m31, -m33];
    [a13, a23] = [a23, a13];
  }

  return [
    [a1, a2, a3, a23, a13, a12],
    [
      [m11, m12, m13],
      [m21, m22, m23],
      [m31, m32, m33],
    ],
  ];
}

/**
 * Coefficients only (SageMath's `_reduced_ternary_form_eisenstein_without_matrix`).
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:311
 */
export function _reduced_ternary_form_eisenstein_without_matrix(
  a1: IntegerLike,
  a2: IntegerLike,
  a3: IntegerLike,
  a23: IntegerLike,
  a13: IntegerLike,
  a12: IntegerLike
): TernaryCoefficients {
  return _reduced_ternary_form_eisenstein_with_matrix(a1, a2, a3, a23, a13, a12)[0];
}

/**
 * A primitive representative of the 3-tuple `v` modulo `p`.
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:470 (`primitivize`)
 */
export function primitivize(
  v0: bigint,
  v1: bigint,
  v2: bigint,
  p: bigint
): [bigint, bigint, bigint] {
  if (pmod(v2, p) !== 0n) {
    const v2_inv = inverse_mod(pmod(v2, p), p);
    return [pmod(v2_inv * v0, p), pmod(v2_inv * v1, p), 1n];
  }
  if (pmod(v1, p) !== 0n) {
    return [pmod(inverse_mod(pmod(v1, p), p) * v0, p), 1n, 0n];
  }
  return [1n, 0n, 0n];
}

/**
 * Evaluate the ternary form `(a, b, c, r, s, t)` at the 3-tuple `v`.
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:491 (`evaluate`)
 */
export function evaluate(
  a: bigint,
  b: bigint,
  c: bigint,
  r: bigint,
  s: bigint,
  t: bigint,
  v: readonly bigint[]
): bigint {
  return (
    a * v[0]! ** 2n +
    b * v[1]! ** 2n +
    c * v[2]! ** 2n +
    r * v[2]! * v[1]! +
    s * v[2]! * v[0]! +
    t * v[1]! * v[0]!
  );
}

/**
 * The zeros mod 2 of the ternary form.
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:508 (`_find_zeros_mod_p_2`)
 */
export function _find_zeros_mod_p_2(
  a: bigint,
  b: bigint,
  c: bigint,
  r: bigint,
  s: bigint,
  t: bigint
): [bigint, bigint, bigint][] {
  const zeros: [bigint, bigint, bigint][] = [];
  let v: [bigint, bigint, bigint] = [1n, 0n, 0n];
  if (pmod(evaluate(a, b, c, r, s, t, v), 2n) === 0n) {
    zeros.push(v);
  }
  for (let i = 0n; i < 2n; i++) {
    v = [i, 1n, 0n];
    if (pmod(evaluate(a, b, c, r, s, t, v), 2n) === 0n) {
      zeros.push(v);
    }
  }
  for (let i = 0n; i < 2n; i++) {
    for (let j = 0n; j < 2n; j++) {
      v = [i, j, 1n];
      if (pmod(evaluate(a, b, c, r, s, t, v), 2n) === 0n) {
        zeros.push(v);
      }
    }
  }
  return zeros;
}

/**
 * A primitive zero `(x, y, 1)` mod the odd prime `p` of the ternary form.
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:542 (`pseudorandom_primitive_zero_mod_p`)
 * @see Deviation: SageMath draws `(r1, r2)` uniformly at random; we scan
 * `r1 = 0, 1, ...` and `r2 = 0, 1, ...` in order so that the result is
 * reproducible.  As a projective set the zeros returned by
 * {@link _find_zeros_mod_p_odd} do not depend on the starting zero, so
 * `find_zeros_mod_p` returns the same set of points of `P^2(F_p)`; the order
 * and the choice of representatives (which SageMath does not always reduce
 * into `[0, p)`, see {@link cmod}) do depend on it.  Given the *same* starting
 * zero, {@link _find_zeros_mod_p_odd} reproduces SageMath's list exactly.
 */
export function primitive_zero_mod_p(
  a: bigint,
  b: bigint,
  c: bigint,
  r: bigint,
  s: bigint,
  t: bigint,
  p: bigint
): [bigint, bigint, bigint] {
  for (let r1 = 0n; r1 < p; r1++) {
    const alpha = pmod(b * r1 * r1 + t * r1 + a, p);
    if (alpha === 0n) continue;
    for (let r2 = 0n; r2 < p; r2++) {
      const beta = pmod(2n * b * r1 * r2 + t * r2 + r * r1 + s, p);
      const gamma = pmod(b * r2 * r2 + r * r2 + c, p);
      const disc = beta * beta - 4n * alpha * gamma;
      const root = sqrt_mod(pmod(disc, p), p);
      if (root === null) continue;
      const z = (-beta + root) * inverse_mod(pmod(2n * alpha, p), p);
      return [pmod(z, p), pmod(r1 * z + r2, p), 1n];
    }
  }
  throw new ValueError('no primitive zero mod p was found');
}

/**
 * All zeros mod the odd prime `p`, starting from a known zero `v`.
 *
 * The reductions use C's `%` ({@link cmod}), because `_find_zeros_mod_p_odd`
 * is declared with `long long` parameters; see the note on {@link cmod}.
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:575 (`_find_zeros_mod_p_odd`)
 */
export function _find_zeros_mod_p_odd(
  a: bigint,
  b: bigint,
  c: bigint,
  r: bigint,
  s: bigint,
  t: bigint,
  p: bigint,
  v: readonly bigint[]
): [bigint, bigint, bigint][] {
  const zeros: [bigint, bigint, bigint][] = [[v[0]!, v[1]!, v[2]!]];
  const x0 = v[0]!;
  const y0 = v[1]!;
  let more = false;
  for (let i = 0n; i < p; i++) {
    const a_i = cmod(
      a * x0 ** 2n +
        b * i ** 2n -
        2n * b * i * y0 +
        b * y0 ** 2n -
        t * i * x0 +
        t * x0 * y0 -
        2n * a * x0 +
        t * i -
        t * y0 +
        s * x0 -
        r * i +
        r * y0 +
        a +
        c -
        s,
      p
    );
    let w: [bigint, bigint, bigint];
    if (a_i === 0n) {
      w = [cmod(x0 - 1n, p), cmod(y0 - i, p), 1n];
    } else {
      const c_i = cmod(b * i ** 2n + t * i + a, p);
      const l = c_i * inverse_mod(a_i, p);
      w = primitivize(l * (x0 - 1n) + 1n, l * (y0 - i) + i, l, p);
    }
    if (w[0] === v[0] && w[1] === v[1] && w[2] === v[2]) {
      more = true;
    } else {
      zeros.push(w);
    }
  }
  if (more) {
    const a_inf = cmod(
      a * x0 ** 2n +
        b * y0 ** 2n +
        t * x0 * y0 -
        2n * b * y0 -
        t * x0 +
        s * x0 +
        r * y0 +
        b +
        c -
        r,
      p
    );
    if (a_inf === 0n) {
      zeros.push([cmod(x0, p), cmod(y0 - 1n, p), 1n]);
    } else {
      const c_inf = cmod(b, p);
      const l = c_inf * inverse_mod(a_inf, p);
      zeros.push(primitivize(l * x0, l * (y0 - 1n) + 1n, l, p));
    }
  }
  return zeros;
}

/**
 * The coefficients of a matrix `M` with `det(M) = gcd(v)` whose first column
 * is `v`.
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:933 (`extend`)
 */
export function extend(v: readonly bigint[]): bigint[] {
  const v0 = v[0]!;
  const v1 = v[1]!;
  const v2 = v[2]!;
  if (v0 === 0n && v1 === 0n) {
    if (v2 < 0n) {
      return [v0, 0n, 1n, v1, 1n, 0n, v2, 0n, 0n];
    }
    return [v0, 1n, 0n, v1, 0n, 1n, v2, 0n, 0n];
  }
  const b1 = xgcd3(v0, v1);
  const b2 = xgcd3(b1[1], b1[2]);
  const b3 = xgcd3(b1[0], v2);
  return [v0, -b1[2], -b2[1] * b3[2], v1, b1[1], -b2[2] * b3[2], v2, 0n, b3[1]];
}

/** `xgcd` returning `(g, s, t)` with `g = s*a + t*b`, as in SageMath. */
function xgcd3(a: bigint, b: bigint): [bigint, bigint, bigint] {
  // Local copy of sage.arith.misc.xgcd's sign conventions (g >= 0).
  let oldR = a;
  let r = b;
  let oldS = 1n;
  let s = 0n;
  let oldT = 0n;
  let t = 1n;
  while (r !== 0n) {
    const q = fdiv(oldR, r);
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
    [oldT, t] = [t, oldT - q * t];
  }
  if (oldR < 0n) {
    return [-oldR, -oldS, -oldT];
  }
  return [oldR, oldS, oldT];
}

/**
 * A number represented by the ternary form and coprime to the prime `n`.
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:1124 (`_basic_lemma`)
 */
export function _basic_lemma(
  a: bigint,
  b: bigint,
  c: bigint,
  r: bigint,
  s: bigint,
  t: bigint,
  n: bigint
): bigint {
  if (n === 1n) {
    return 0n;
  }
  if (pmod(a, n) !== 0n) return a;
  if (pmod(b, n) !== 0n) return b;
  if (pmod(c, n) !== 0n) return c;
  if (pmod(r, n) !== 0n) return b + c + r;
  if (pmod(s, n) !== 0n) return a + c + s;
  if (pmod(t, n) !== 0n) return a + b + t;
  throw new ValueError('not primitive form');
}

/**
 * A vector at which the ternary form takes a value coprime to `n`.
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:1089 (`_basic_lemma_vec`)
 */
export function _basic_lemma_vec(
  a: bigint,
  b: bigint,
  c: bigint,
  r: bigint,
  s: bigint,
  t: bigint,
  n: bigint
): [bigint, bigint, bigint] {
  if (n === 1n) {
    return [0n, 0n, 0n];
  }
  if (pmod(a, n) !== 0n) return [1n, 0n, 0n];
  if (pmod(b, n) !== 0n) return [0n, 1n, 0n];
  if (pmod(c, n) !== 0n) return [0n, 0n, 1n];
  if (pmod(r, n) !== 0n) return [0n, 1n, 1n];
  if (pmod(s, n) !== 0n) return [1n, 0n, 1n];
  if (pmod(t, n) !== 0n) return [1n, 1n, 0n];
  throw new ValueError('not primitive form');
}

/**
 * The reduced `p`-neighbour of the ternary form associated with the vector
 * `v`, and (optionally) the rational transformation matrix.
 *
 * @see Reference: sage/quadratic_forms/ternary.pyx:984 (`_find_p_neighbor_from_vec`)
 */
export function _find_p_neighbor_from_vec(
  a: bigint,
  b: bigint,
  c: bigint,
  r: bigint,
  s: bigint,
  t: bigint,
  p: bigint,
  v: readonly bigint[]
): [TernaryCoefficients, Rational[][]] {
  const [v0, w0, u0, v1, w1, u1, v2, w2, u2] = extend(v) as [
    bigint,
    bigint,
    bigint,
    bigint,
    bigint,
    bigint,
    bigint,
    bigint,
    bigint,
  ];

  const m00 =
    2n * a * v0 ** 2n +
    2n * b * v1 ** 2n +
    2n * c * v2 ** 2n +
    2n * r * v1 * v2 +
    2n * s * v0 * v2 +
    2n * t * v0 * v1;
  const m11 = 2n * a * w0 ** 2n + 2n * b * w1 ** 2n + 2n * t * w0 * w1;
  const m22 =
    2n * a * u0 ** 2n +
    2n * b * u1 ** 2n +
    2n * c * u2 ** 2n +
    2n * r * u1 * u2 +
    2n * s * u0 * u2 +
    2n * t * u0 * u1;
  const m01 =
    2n * a * v0 * w0 + 2n * b * v1 * w1 + r * v2 * w1 + s * v2 * w0 + t * v0 * w1 + t * v1 * w0;
  const m02 =
    2n * a * u0 * v0 +
    2n * b * u1 * v1 +
    2n * c * u2 * v2 +
    r * u1 * v2 +
    r * u2 * v1 +
    s * u0 * v2 +
    s * u2 * v0 +
    t * u0 * v1 +
    t * u1 * v0;
  const m12 =
    2n * a * u0 * w0 + 2n * b * u1 * w1 + r * u2 * w1 + s * u2 * w0 + t * u0 * w1 + t * u1 * w0;

  const R = (x: bigint, y = 1n) => new Rational(x, y);

  /** Shared tail of the two symmetric branches of `_find_p_neighbor_from_vec`. */
  const finish = (
    m0: bigint,
    m1: bigint,
    mAA: bigint, // m11 resp. m22 (the "pivot" diagonal entry)
    mAB: bigint, // m01 resp. m02 (the entry that is nonzero mod p)
    mOther: bigint, // m22 resp. m11
    mCross: bigint, // m12
    mCol: bigint, // m02 resp. m01
    e0: bigint,
    e1: bigint,
    e2: bigint, // the direction used for the first column
    f0: bigint,
    f1: bigint,
    f2: bigint // the direction used for the second column
  ): [TernaryCoefficients, Rational[][]] => {
    const p2 = new Rational(p * p, 1n);
    const b00 = R(m0 * m0 * mAA)
      .div(p2)
      .add(R(2n * m0 * mAB).div(p2))
      .add(R(m00).div(p2));
    const b11 = R(m1 * m1 * mAA + 2n * m1 * mCross + mOther);
    const b22 = R(mAA * p * p);
    const b01 = R(m0 * m1 * mAA)
      .div(p)
      .add(R(m0 * mCross).div(p))
      .add(R(mAB * m1).div(p))
      .add(R(mCol).div(p));
    const b02 = R(m0 * mAA + mAB);
    const b12 = R(m1 * mAA * p + mCross * p);

    const [q, Mr] = _reduced_ternary_form_eisenstein_with_matrix(
      toZZ(b00.div(2n), 'b00/2'),
      toZZ(b11.div(2n), 'b11/2'),
      toZZ(b22.div(2n), 'b22/2'),
      toZZ(b12, 'b12'),
      toZZ(b02, 'b02'),
      toZZ(b01, 'b01')
    );

    const T: Rational[][] = [];
    const ds = [
      [e0, f0],
      [e1, f1],
      [e2, f2],
    ];
    // v-part of the columns: (m0*d/p + v_k/p), (m1*d + f_k), p*d
    const vcomp = [v0, v1, v2];
    for (let k = 0; k < 3; k++) {
      const d = ds[k]![0]!; // e_k  (u_k resp. w_k)
      const g = ds[k]![1]!; // f_k  (w_k resp. u_k)
      const row: Rational[] = [];
      for (let j = 0; j < 3; j++) {
        const term = R(p * Mr[2]![j]! * d)
          .add(
            R(m0 * d)
              .div(p)
              .add(R(vcomp[k]!).div(p))
              .mul(Mr[0]![j]!)
          )
          .add(R(m1 * d + g).mul(Mr[1]![j]!));
        row.push(term);
      }
      T.push(row);
    }
    return [q, T];
  };

  if (pmod(m02, p) !== 0n) {
    const m0 = rational_mod(new Rational(-m00, m02).div(2n), p * p);
    const m1 = rational_mod(new Rational(-m01, m02), p);
    return finish(m0, m1, m22, m02, m11, m12, m01, u0, u1, u2, w0, w1, w2);
  }
  if (pmod(m01, p) !== 0n) {
    const m0 = rational_mod(new Rational(-m00, m01).div(2n), p * p);
    const m1 = rational_mod(new Rational(-m02, m01), p);
    return finish(m0, m1, m11, m01, m22, m12, m02, w0, w1, w2, u0, u1, u2);
  }
  throw new ValueError(
    'v is a singular point of the conic Q(v) = 0 mod p; no p-neighbor is defined'
  );
}

/* ------------------------------------------------------------------ */
/* The TernaryQF class                                                 */
/* ------------------------------------------------------------------ */

/**
 * A quadratic form in 3 variables with coefficients in `ZZ`.
 *
 * `TernaryQF([a, b, c, r, s, t])` is
 * `a x^2 + b y^2 + c z^2 + r y z + s x z + t x y`.
 *
 * @see Reference: sage/quadratic_forms/ternary_qf.py:49
 */
export class TernaryQF {
  readonly _a: bigint;
  readonly _b: bigint;
  readonly _c: bigint;
  readonly _r: bigint;
  readonly _s: bigint;
  readonly _t: bigint;

  /**
   * @param v six integers `[a, b, c, r, s, t]`
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:78 (`__init__`)
   */
  constructor(v: RationalLike[]) {
    if (v.length !== 6) {
      throw new ValueError('Ternary quadratic form must be given by a list of six coefficients');
    }
    // SageMath: `ZZ(x) for x in v` -- non-integral entries are an error.
    [this._a, this._b, this._c, this._r, this._s, this._t] = v.map((x) =>
      toZZ(toRational(x), 'TernaryQF coefficient')
    ) as TernaryCoefficients;
  }

  /**
   * The tuple `(a, b, c, r, s, t)`.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:101
   */
  coefficients(): TernaryCoefficients {
    return [this._a, this._b, this._c, this._r, this._s, this._t];
  }

  /**
   * The `n`-th coefficient, `0 <= n <= 5`.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:130
   */
  coefficient(n: number | IntegerLike): bigint {
    const i = typeof n === 'number' ? n : Number(toBigInt(n));
    const c = this.coefficients()[i];
    if (c === undefined) {
      throw new ValueError(`coefficient index ${i} out of range`);
    }
    return c;
  }

  /**
   * @see Reference: sage/quadratic_forms/ternary_qf.py:169 (`_repr_`)
   */
  _repr_(): string {
    let rep = 'Ternary quadratic form with integer coefficients:\n';
    rep += `[${this._a} ${this._b} ${this._c}]\n`;
    rep += `[${this._r} ${this._s} ${this._t}]`;
    return rep;
  }

  toString(): string {
    return this._repr_();
  }

  /**
   * Evaluate the form on a 3-vector (giving an integer) or on a matrix with 3
   * rows (giving a `TernaryQF` when it has 3 columns, else a `QuadraticForm`).
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:190 (`__call__`)
   */
  __call__(v: IntegerLike[]): bigint;
  __call__(v: RationalLike[][] | RationalMatrix): TernaryQF | QuadraticForm;
  __call__(
    v: IntegerLike[] | RationalLike[][] | RationalMatrix
  ): bigint | TernaryQF | QuadraticForm;
  __call__(
    v: IntegerLike[] | RationalLike[][] | RationalMatrix
  ): bigint | TernaryQF | QuadraticForm {
    const isMatrix = isRationalMatrix(v) || (Array.isArray(v) && Array.isArray(v[0]));
    if (isMatrix) {
      const rows: Rational[][] = isRationalMatrix(v)
        ? v.rows().map((row) => row.slice())
        : (v as RationalLike[][]).map((row) => row.map(toRational));
      if (rows.length !== 3) {
        throw new TypeError('the matrix must have 3 rows');
      }
      const ncols = rows[0]!.length;
      const H = this.matrixZZ();
      // M^t * H * M
      const prod: Rational[][] = [];
      for (let i = 0; i < ncols; i++) {
        const row: Rational[] = [];
        for (let j = 0; j < ncols; j++) {
          let sum = Rational.zero();
          for (let k = 0; k < 3; k++) {
            for (let l = 0; l < 3; l++) {
              sum = sum.add(rows[k]![i]!.mul(H[k]![l]!).mul(rows[l]![j]!));
            }
          }
          row.push(sum);
        }
        prod.push(row);
      }
      if (ncols === 3) {
        // SageMath uses `M[i,i] // 2`; the diagonal of M^t H M is even
        // whenever the result really is a ternary form over ZZ.
        const half = (x: Rational) => new Rational(x.div(2n).floor(), 1n);
        return new TernaryQF([
          half(prod[0]![0]!),
          half(prod[1]![1]!),
          half(prod[2]![2]!),
          prod[1]![2]!,
          prod[0]![2]!,
          prod[0]![1]!,
        ]);
      }
      return new QuadraticForm(ZZ, prod);
    }
    const vec = (v as IntegerLike[]).map(toBigInt);
    if (vec.length !== 3) {
      throw new TypeError('your vector needs to have length 3');
    }
    return evaluate(this._a, this._b, this._c, this._r, this._s, this._t, vec);
  }

  /**
   * The same form as a {@link QuadraticForm} over `ZZ`.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:247
   */
  quadratic_form(): QuadraticForm {
    return new QuadraticForm(ZZ, 3n, [this._a, this._t, this._s, this._b, this._r, this._c]);
  }

  /** The Hessian matrix as a plain `bigint[][]` (internal helper). */
  private matrixZZ(): bigint[][] {
    return [
      [2n * this._a, this._t, this._s],
      [this._t, 2n * this._b, this._r],
      [this._s, this._r, 2n * this._c],
    ];
  }

  /**
   * The Hessian matrix of the form.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:266
   * @see Deviation: returned over `QQ` (SageMath returns it over `ZZ`); the
   * entries are unchanged and always integral.
   */
  matrix(): RationalMatrix {
    const M = this.matrixZZ();
    const flat: Rational[] = [];
    for (const row of M) {
      for (const x of row) {
        flat.push(new Rational(x, 1n));
      }
    }
    return matrixQQ(3, 3, flat);
  }

  /**
   * The discriminant: `det(Hessian) / 2`.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:298
   */
  disc(): bigint {
    return (
      4n * this._a * this._b * this._c +
      this._r * this._s * this._t -
      this._a * this._r ** 2n -
      this._b * this._s ** 2n -
      this._c * this._t ** 2n
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/ternary_qf.py:315
   */
  is_definite(): boolean {
    const d1 = this._a;
    if (d1 === 0n) return false;
    const d2 = 4n * this._a * this._b - this._t ** 2n;
    if (d2 <= 0n) return false;
    const d3 = this.disc();
    if (d3 === 0n) return false;
    return d1 > 0n === d3 > 0n;
  }

  /**
   * @see Reference: sage/quadratic_forms/ternary_qf.py:341
   */
  is_positive_definite(): boolean {
    const d1 = this._a;
    if (d1 === 0n) return false;
    const d2 = 4n * this._a * this._b - this._t ** 2n;
    if (d2 <= 0n) return false;
    const d3 = this.disc();
    if (d3 === 0n) return false;
    return d1 > 0n && d3 > 0n;
  }

  /**
   * @see Reference: sage/quadratic_forms/ternary_qf.py:372
   */
  is_negative_definite(): boolean {
    const d1 = this._a;
    if (d1 === 0n) return false;
    const d2 = 4n * this._a * this._b - this._t ** 2n;
    if (d2 <= 0n) return false;
    const d3 = this.disc();
    if (d3 === 0n) return false;
    return d1 < 0n && d3 < 0n;
  }

  /**
   * The form with all coefficients negated.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:398 (`__neg__`)
   */
  neg(): TernaryQF {
    return new TernaryQF(this.coefficients().map((x) => -x));
  }

  /**
   * @see Reference: sage/quadratic_forms/ternary_qf.py:418
   */
  is_primitive(): boolean {
    return this.content() === 1n;
  }

  /**
   * @see Reference: sage/quadratic_forms/ternary_qf.py:440
   */
  primitive(): TernaryQF {
    const l = this.coefficients();
    const g = gcd(l);
    return new TernaryQF(l.map((x) => x / g));
  }

  /**
   * Scale the values of the form by `k`.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:465
   */
  scale_by_factor(k: RationalLike): TernaryQF | QuadraticForm {
    const _k = toRational(k);
    if (_k.mul(this.content()).isInteger()) {
      return new TernaryQF(this.coefficients().map((x) => _k.mul(x).numerator));
    }
    return new QuadraticForm(QQ, 3n, [
      _k.mul(this._a),
      _k.mul(this._t),
      _k.mul(this._s),
      _k.mul(this._b),
      _k.mul(this._r),
      _k.mul(this._c),
    ]);
  }

  /**
   * The reciprocal form: the multiple of the primitive adjoint with the same
   * content as the given form.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:506
   */
  reciprocal(): TernaryQF {
    return this.adjoint().primitive().scale_by_factor(this.content()) as TernaryQF;
  }

  /**
   * @see Reference: sage/quadratic_forms/ternary_qf.py:529
   */
  reciprocal_reduced(): TernaryQF {
    return this.reciprocal().reduced_form_eisenstein(false);
  }

  /**
   * The content of the adjoint form.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:548
   */
  divisor(): bigint {
    const A11 = 4n * this._b * this._c - this._r ** 2n;
    const A22 = 4n * this._a * this._c - this._s ** 2n;
    const A33 = 4n * this._a * this._b - this._t ** 2n;
    const A23 = this._s * this._t - 2n * this._a * this._r;
    const A13 = this._r * this._t - 2n * this._b * this._s;
    const A12 = this._r * this._s - 2n * this._c * this._t;
    return gcd([A11, A22, A33, 2n * A12, 2n * A13, 2n * A23]);
  }

  /**
   * @see Reference: sage/quadratic_forms/ternary_qf.py:567 (`__eq__`)
   */
  equals(right: unknown): boolean {
    if (!(right instanceof TernaryQF)) {
      return false;
    }
    const l = this.coefficients();
    const m = right.coefficients();
    return l.every((x, i) => x === m[i]);
  }

  /**
   * The adjoint form: its Hessian matrix is twice the classical adjoint of the
   * Hessian matrix of the given form.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:584
   */
  adjoint(): TernaryQF {
    const A11 = 4n * this._b * this._c - this._r ** 2n;
    const A22 = 4n * this._a * this._c - this._s ** 2n;
    const A33 = 4n * this._a * this._b - this._t ** 2n;
    const A23 = this._s * this._t - 2n * this._a * this._r;
    const A13 = this._r * this._t - 2n * this._b * this._s;
    const A12 = this._r * this._s - 2n * this._c * this._t;
    return new TernaryQF([A11, A22, A33, 2n * A23, 2n * A13, 2n * A12]);
  }

  /**
   * The gcd of the coefficients.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:610
   */
  content(): bigint {
    return gcd(this.coefficients());
  }

  /**
   * The content of the adjoint of the associated primitive form.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:627
   */
  omega(): bigint {
    return this.primitive().adjoint().content();
  }

  /**
   * The omega of the adjoint form.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:642
   */
  delta(): bigint {
    return this.adjoint().omega();
  }

  /**
   * `4 * disc / divisor`.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:662
   */
  level(): bigint {
    return fdiv(4n * this.disc(), this.divisor());
  }

  /**
   * Whether the form is Eisenstein reduced (Dickson's conditions).
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:676
   */
  is_eisenstein_reduced(): boolean {
    const [a, b, c, r, s, t] = this.coefficients();

    // cond 2
    if (!(r > 0n && t > 0n && s > 0n)) {
      if (!(r <= 0n && s <= 0n && t <= 0n)) {
        return false;
      }
    }

    // cond 1 & 4
    if (!(a <= b && b <= c && 0n <= a + b + r + s + t)) {
      return false;
    }

    // cond 3
    if (!(a >= babs(s) && a >= babs(t) && b >= babs(r))) {
      return false;
    }

    // cond 8
    if (a === b && babs(r) > babs(s)) return false;
    if (b === c && babs(s) > babs(t)) return false;
    if (a + b + r + s + t === 0n && 2n * a + 2n * s + t > 0n) return false;

    // cond 6 (r, s, t <= 0)
    if (r <= 0n) {
      if (a === -t && s !== 0n) return false;
      if (a === -s && t !== 0n) return false;
      if (b === -r && t !== 0n) return false;
    }

    // cond 7 (r, s, t > 0)
    if (a === t && s > 2n * r) return false;
    if (a === s && t > 2n * r) return false;
    return !(b === r && t > 2n * s);
  }

  /**
   * The unique Eisenstein-reduced form equivalent to this positive definite
   * form, with (by default) the transformation matrix.
   *
   * The form must be positive definite: as in SageMath, the reduction loop of
   * `ternary.pyx` does not terminate for indefinite or negative definite
   * forms.  Check {@link is_positive_definite} first.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:748
   */
  reduced_form_eisenstein(): [TernaryQF, bigint[][]];
  reduced_form_eisenstein(withMatrix: true): [TernaryQF, bigint[][]];
  reduced_form_eisenstein(withMatrix: false): TernaryQF;
  reduced_form_eisenstein(withMatrix = true): TernaryQF | [TernaryQF, bigint[][]] {
    const [v, M] = _reduced_ternary_form_eisenstein_with_matrix(
      this._a,
      this._b,
      this._c,
      this._r,
      this._s,
      this._t
    );
    return withMatrix ? [new TernaryQF(v), M] : new TernaryQF(v);
  }

  /**
   * A primitive zero `(x, y, 1)` of the form modulo the odd prime `p`.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:784
   * @see Deviation: deterministic (see {@link primitive_zero_mod_p}).
   */
  pseudorandom_primitive_zero_mod_p(p: IntegerLike): [bigint, bigint, bigint] {
    const _p = toBigInt(p);
    return primitive_zero_mod_p(this._a, this._b, this._c, this._r, this._s, this._t, _p);
  }

  /**
   * The zeros of the form modulo `p`, for `p` not dividing the discriminant.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:821
   */
  find_zeros_mod_p(p: IntegerLike): [bigint, bigint, bigint][] {
    const _p = toBigInt(p);
    if (_p === 2n) {
      return _find_zeros_mod_p_2(this._a, this._b, this._c, this._r, this._s, this._t);
    }
    const v = this.pseudorandom_primitive_zero_mod_p(_p);
    return _find_zeros_mod_p_odd(this._a, this._b, this._c, this._r, this._s, this._t, _p, v);
  }

  /**
   * The reduced `p`-neighbour associated with the vector `v`.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:848
   */
  find_p_neighbor_from_vec(p: IntegerLike, v: IntegerLike[]): TernaryQF;
  find_p_neighbor_from_vec(p: IntegerLike, v: IntegerLike[], mat: false): TernaryQF;
  find_p_neighbor_from_vec(
    p: IntegerLike,
    v: IntegerLike[],
    mat: true
  ): [TernaryQF, RationalMatrix];
  find_p_neighbor_from_vec(
    p: IntegerLike,
    v: IntegerLike[],
    mat = false
  ): TernaryQF | [TernaryQF, RationalMatrix] {
    const _p = toBigInt(p);
    const _v = v.map(toBigInt);
    const [q, T] = _find_p_neighbor_from_vec(
      this._a,
      this._b,
      this._c,
      this._r,
      this._s,
      this._t,
      _p,
      _v
    );
    if (!mat) {
      return new TernaryQF(q);
    }
    const M = matrixQQ(3, 3, T.flat());
    const d = determinantQQ(M);
    return [new TernaryQF(q), M.scalar_mul(d)];
  }

  /**
   * All reduced `p`-neighbours of the form.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:899
   */
  find_p_neighbors(p: IntegerLike): TernaryQF[] {
    return this.find_zeros_mod_p(p).map((v) => this.find_p_neighbor_from_vec(p, v));
  }

  /**
   * A number represented by the form and coprime to the prime `p`.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:926
   */
  basic_lemma(p: IntegerLike): bigint {
    return _basic_lemma(this._a, this._b, this._c, this._r, this._s, this._t, toBigInt(p));
  }

  /**
   * The genus character `Xi_p`.  `-1` is allowed as a "prime".
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:938
   */
  xi(p: IntegerLike): bigint {
    let _p = toBigInt(p);
    if (_p === 4n) _p = -1n;
    if (_p === 8n) _p = 2n;

    const om = this.omega();
    if (pmod(om, _p) !== 0n) {
      throw new ValueError('not a valid character');
    }
    if (_p === -1n && pmod(om, 2n ** 4n) !== 0n) {
      throw new ValueError('not a valid character');
    }
    if (_p === 2n && pmod(om, 2n ** 5n) !== 0n) {
      throw new ValueError('not a valid character');
    }
    if (_p === -1n || _p === 2n) {
      return kronecker_symbol(_p, this.basic_lemma(2n));
    }
    return kronecker_symbol(this.basic_lemma(_p), _p);
  }

  /**
   * `xi(p)` of the reciprocal form.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:975
   */
  xi_rec(p: IntegerLike): bigint {
    return this.reciprocal().xi(p);
  }

  /**
   * The reflection `A` with `A v = -v` and `A u = u` for `u` orthogonal to `v`.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:990
   */
  symmetry(v: IntegerLike[]): RationalMatrix {
    const _v = v.map(toBigInt);
    if (_v.length !== 3) {
      throw new TypeError('your vector needs to have length 3');
    }
    const Qv = this.__call__(_v) as bigint;
    if (Qv === 0n) {
      throw new ValueError('the symmetry is not defined for an isotropic vector');
    }
    const H = this.matrixZZ();
    // identity - v.column() * matrix(v) * self.matrix() / self(v)
    // (v.column() * matrix(v))[i][j] = v_i v_j
    const entries: Rational[] = [];
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        let vHj = 0n;
        for (let k = 0; k < 3; k++) {
          vHj += _v[k]! * H[k]![j]!;
        }
        const off = new Rational(_v[i]! * vHj, Qv);
        const id = i === j ? Rational.one() : Rational.zero();
        entries.push(id.sub(off));
      }
    }
    return matrixQQ(3, 3, entries);
  }

  /**
   * Two vectors whose symmetries compose to the automorphism `A` (empty for
   * the identity).
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:1026
   */
  automorphism_symmetries(A: RationalMatrix | RationalLike[][]): Rational[][] {
    const M = isRationalMatrix(A)
      ? A
      : matrixQQ(3, 3, (A as RationalLike[][]).flat().map(toRational));
    const minusOne = (X: RationalMatrix): RationalMatrix => {
      const entries: Rational[] = [];
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          entries.push(X.get(i, j).sub(i === j ? Rational.one() : Rational.zero()));
        }
      }
      return matrixQQ(3, 3, entries);
    };
    const A_minus_1 = minusOne(M);
    if (A_minus_1.is_zero()) {
      return [];
    }
    const firstNonzeroColumn = (X: RationalMatrix): Rational[] => {
      for (let j = 0; j < 3; j++) {
        const col = X.column(j);
        if (col.some((x) => !x.isZero())) {
          return col;
        }
      }
      throw new ValueError('A - 1 has no nonzero column');
    };

    const b1 = firstNonzeroColumn(A_minus_1);
    const A1 = this.symmetryQ(b1).mul(M);
    const b2 = firstNonzeroColumn(minusOne(A1));
    return [b1, b2];
  }

  /** `Q(v)` for a rational vector `v` (internal). */
  private evaluateQ(v: Rational[]): Rational {
    const H = this.matrixZZ();
    let Qv = Rational.zero();
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        Qv = Qv.add(v[i]!.mul(v[j]!).mul(H[i]![j]!));
      }
    }
    return Qv.div(2n);
  }

  /** {@link symmetry} for a rational vector (internal). */
  private symmetryQ(v: Rational[]): RationalMatrix {
    const H = this.matrixZZ();
    const Qv = this.evaluateQ(v);
    if (Qv.isZero()) {
      throw new ValueError('the symmetry is not defined for an isotropic vector');
    }
    const entries: Rational[] = [];
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        let vHj = Rational.zero();
        for (let k = 0; k < 3; k++) {
          vHj = vHj.add(v[k]!.mul(H[k]![j]!));
        }
        const off = v[i]!.mul(vHj).div(Qv);
        const id = i === j ? Rational.one() : Rational.zero();
        entries.push(id.sub(off));
      }
    }
    return matrixQQ(3, 3, entries);
  }

  /**
   * The spin norm of the automorphism `A`.
   *
   * @see Reference: sage/quadratic_forms/ternary_qf.py:1066
   */
  automorphism_spin_norm(A: RationalMatrix | RationalLike[][]): bigint {
    const syms = this.automorphism_symmetries(A);
    if (syms.length === 0) {
      return 1n;
    }
    const s = this.evaluateQ(syms[0]!).mul(this.evaluateQ(syms[1]!));
    // Rational.squarefree_part == numer().squarefree_part() * denom().squarefree_part()
    // (sage/rings/rational.pyx:1749)
    return squarefree_part(s.numerator) * squarefree_part(s.denominator);
  }

  /* --------------------------------------------------------------- */
  /* Honest stubs                                                     */
  /* --------------------------------------------------------------- */

  /**
   * @see Reference: sage/quadratic_forms/ternary_qf.py:151
   * @throws {NotImplementedError} always
   */
  polynomial(_names = 'x,y,z'): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: TernaryQF.polynomial (needs the multivariate polynomial ring plumbing)'
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/ternary_qf.py:1675
   * @throws {NotImplementedError} always
   */
  automorphisms(_slow = true): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: TernaryQF.automorphisms (requires the _border/_borders tables ' +
        'of sage/quadratic_forms/ternary_qf.py:1087-1635)'
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/ternary_qf.py:1918
   * @throws {NotImplementedError} always
   */
  number_of_automorphisms(_slow = true): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: TernaryQF.number_of_automorphisms (requires TernaryQF.automorphisms)'
    );
  }
}

/**
 * @see Reference: sage/quadratic_forms/ternary_qf.py:1957
 * @throws {NotImplementedError} always
 */
export function find_all_ternary_qf_by_level_disc(_N: IntegerLike, _d: IntegerLike): never {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: find_all_ternary_qf_by_level_disc ' +
      '(requires _find_all_ternary_qf_by_level_disc, sage/quadratic_forms/ternary.pyx:675)'
  );
}

/**
 * @see Reference: sage/quadratic_forms/ternary_qf.py:1994
 * @throws {NotImplementedError} always
 */
export function find_a_ternary_qf_by_level_disc(_N: IntegerLike, _d: IntegerLike): never {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: find_a_ternary_qf_by_level_disc ' +
      '(requires _find_a_ternary_qf_by_level_disc, sage/quadratic_forms/ternary.pyx:806)'
  );
}
