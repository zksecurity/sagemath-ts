/**
 * @module sage/quadratic_forms/binary_qf
 * @description Binary quadratic forms with integer coefficients
 *
 * This module provides a specialized class for working with a binary quadratic
 * form `a x^2 + b x y + c y^2`, stored as a triple of integers `(a, b, c)`.
 *
 * Port of: sage/quadratic_forms/binary_qf.py
 * Reference: reference/sage/src/sage/quadratic_forms/binary_qf.py
 */

import { gcd, isqrt, xgcd } from '../arith/misc.js';
import { NotImplementedError, ValueError } from '../errors.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';

/** A 2x2 integer matrix, stored row-major. */
export type Matrix2 = [[bigint, bigint], [bigint, bigint]];

const IDENTITY2: Matrix2 = [
  [1n, 0n],
  [0n, 1n],
];

/** Matrix product of two 2x2 integer matrices. */
function matmul(X: Matrix2, Y: Matrix2): Matrix2 {
  return [
    [X[0][0] * Y[0][0] + X[0][1] * Y[1][0], X[0][0] * Y[0][1] + X[0][1] * Y[1][1]],
    [X[1][0] * Y[0][0] + X[1][1] * Y[1][0], X[1][0] * Y[0][1] + X[1][1] * Y[1][1]],
  ];
}

function abs(x: bigint): bigint {
  return x < 0n ? -x : x;
}

/**
 * Floor division of integers (Python `//` semantics), as opposed to the
 * truncating division of BigInt `/`.
 */
function floorDiv(a: bigint, b: bigint): bigint {
  const q = a / b;
  if (a % b !== 0n && a < 0n !== b < 0n) {
    return q - 1n;
  }
  return q;
}

/**
 * Python/Sage `Integer.quo_rem`: floor quotient and remainder with the sign of
 * the divisor.
 */
function quoRem(a: bigint, b: bigint): [bigint, bigint] {
  const q = floorDiv(a, b);
  return [q, a - q * b];
}

/**
 * PARI's `dvmdii_round`: assuming `a > 0`, write `b = q*2a + r` with
 * `-a < r <= a`; return `[q, r]`.
 *
 * @see Reference: pari/src/basemath/Qfb.c:dvmdii_round (line 188)
 */
function dvmdiiRound(b: bigint, a: bigint): [bigint, bigint] {
  const a2 = 2n * a;
  // PARI's dvmdii truncates towards zero.
  let q = b / a2;
  let r = b - q * a2;
  if (b >= 0n) {
    if (abs(r) > a) {
      q += 1n;
      r -= a2;
    }
  } else {
    if (abs(r) >= a) {
      q -= 1n;
      r += a2;
    }
  }
  return [q, r];
}

/**
 * A binary quadratic form over the integers.
 * Represents the form `a*x^2 + b*x*y + c*y^2` for integer coefficients (a, b, c).
 */
export class BinaryQF {
  readonly a: bigint;
  readonly b: bigint;
  readonly c: bigint;
  /** Cached value of `cycle()` (Sage caches this on the form as well). */
  private _cycle_list?: BinaryQF[];

  constructor(
    a: IntegerLike | [IntegerLike, IntegerLike, IntegerLike],
    b?: IntegerLike,
    c?: IntegerLike
  ) {
    if (Array.isArray(a)) {
      if (a.length !== 3) {
        throw new TypeError(
          'binary quadratic form must be given by a quadratic homogeneous bivariate integer polynomial or its coefficients'
        );
      }
      this.a = toBigInt(a[0]!);
      this.b = toBigInt(a[1]!);
      this.c = toBigInt(a[2]!);
    } else {
      if (b === undefined || c === undefined) {
        throw new TypeError(
          'binary quadratic form must be given by a quadratic homogeneous bivariate integer polynomial or its coefficients'
        );
      }
      this.a = toBigInt(a);
      this.b = toBigInt(b);
      this.c = toBigInt(c);
    }
  }

  static principal(D: IntegerLike): BinaryQF {
    const _D = toBigInt(D);
    const D4 = ((_D % 4n) + 4n) % 4n;
    if (D4 !== 0n && D4 !== 1n) {
      throw new ValueError('discriminant must be congruent to 0 or 1 modulo 4');
    }
    return new BinaryQF(1n, D4, (D4 - _D) / 4n);
  }

  toTuple(): [bigint, bigint, bigint] {
    return [this.a, this.b, this.c];
  }

  evaluate(x: IntegerLike, y: IntegerLike): bigint {
    const _x = toBigInt(x);
    const _y = toBigInt(y);
    return (this.a * _x + this.b * _y) * _x + this.c * _y * _y;
  }

  discriminant(): bigint {
    return this.b * this.b - 4n * this.a * this.c;
  }
  content(): bigint {
    return gcd([this.a, this.b, this.c]);
  }
  is_primitive(): boolean {
    return this.content() === 1n;
  }
  is_zero(): boolean {
    return this.a === 0n && this.b === 0n && this.c === 0n;
  }
  is_positive_definite(): boolean {
    return this.discriminant() < 0n && this.a > 0n;
  }
  is_negative_definite(): boolean {
    return this.discriminant() < 0n && this.a < 0n;
  }
  is_indefinite(): boolean {
    return this.discriminant() > 0n;
  }
  is_singular(): boolean {
    return this.discriminant() === 0n;
  }
  is_nonsingular(): boolean {
    return this.discriminant() !== 0n;
  }

  is_reducible(): boolean {
    const D = this.discriminant();
    if (D < 0n) return false;
    const s = isqrt(D);
    return s * s === D;
  }

  is_reduced(): boolean {
    const D = this.discriminant();
    if (D === 0n) throw new ValueError('the quadratic form must be non-singular');
    const { a, b, c } = this;

    if (D < 0n && a > 0n) {
      return (-a < b && b <= a && a < c) || (0n <= b && b <= a && a === c);
    } else if (D < 0n && a < 0n) {
      return (a < b && b <= -a && -a < -c) || (0n <= b && b <= -a && -a === -c);
    } else {
      return (
        (b > 0n && a * c < 0n && (a - c) * (a - c) < D) ||
        (a === 0n && -b < 2n * c && 2n * c <= b) ||
        (c === 0n && -b < 2n * a && 2n * a <= b)
      );
    }
  }

  /**
   * Return a reduced form equivalent to this one.
   *
   * With `{transformation: true}` also return a matrix `M` in `SL_2(Z)` such
   * that `this.matrix_action_right(M)` is the reduced form.
   *
   * SageMath delegates definite reduction to PARI (`qfbred`/`qfbredsl2`) and
   * reduces forms of square discriminant itself (`_reduce_indef`); we port both
   * algorithms directly.
   *
   * @see Reference: sage/quadratic_forms/binary_qf.py:reduced_form (line 831)
   * @see Deviation: PARI qfb family ported in-place
   */
  reduced_form(options?: { transformation?: false }): BinaryQF;
  reduced_form(options: { transformation: true }): [BinaryQF, Matrix2];
  reduced_form(options?: { transformation?: boolean }): BinaryQF | [BinaryQF, Matrix2] {
    const transformation = options?.transformation ?? false;
    if (this.is_reduced()) {
      return transformation ? [this, IDENTITY2] : this;
    }
    const D = this.discriminant();
    if (D > 0n) {
      // Sage: algorithm 'sage' for reducible forms, PARI's qfr_red otherwise.
      // Both iterate rho until reduced (see pari/src/basemath/Qfb.c:qfr_redsl2_basecase),
      // so a single implementation covers both.
      return transformation ? this._reduce_indef(true) : this._reduce_indef(false);
    }
    if (this.is_negative_definite()) {
      // PARI does not support negative definite forms; Sage reduces
      // (-self)*M with M = diag(-1, 1) instead and conjugates the result back.
      // sage/quadratic_forms/binary_qf.py:960-966
      const negForm = new BinaryQF(-this.a, this.b, -this.c);
      if (transformation) {
        const [reduced, M] = negForm.reduced_form({ transformation: true });
        // M_diag * M * M_diag with M_diag = diag(-1, 1)
        const conjugated: Matrix2 = [
          [M[0][0], -M[0][1]],
          [-M[1][0], M[1][1]],
        ];
        return [new BinaryQF(-reduced.a, reduced.b, -reduced.c), conjugated];
      }
      const reduced = negForm.reduced_form();
      return new BinaryQF(-reduced.a, reduced.b, -reduced.c);
    }
    return transformation
      ? this._reduce_positive_definite(true)
      : this._reduce_positive_definite(false);
  }

  /**
   * Reduce a positive definite form, following PARI's `qfi_redsl2_basecase`.
   *
   * @see Reference: pari/src/basemath/Qfb.c:qfi_redsl2_basecase (line 274)
   */
  private _reduce_positive_definite(transformation: false): BinaryQF;
  private _reduce_positive_definite(transformation: true): [BinaryQF, Matrix2];
  private _reduce_positive_definite(transformation: boolean): BinaryQF | [BinaryQF, Matrix2] {
    let a = this.a;
    let b = this.b;
    let c = this.c;
    let U: Matrix2 = IDENTITY2;

    // REDBU: normalize b into (-a, a] and update c and the base change.
    // The substitution is x -> x - q*y, i.e. T = [[1, -q], [0, 1]].
    const redbu = () => {
      const [q, r] = dvmdiiRound(b, a);
      c = c - q * ((b + r) / 2n);
      b = r;
      if (transformation && q !== 0n) {
        U = matmul(U, [
          [1n, -q],
          [0n, 1n],
        ]);
      }
    };
    // Swap a and c: (a, b, c) -> (c, -b, a), i.e. T = [[0, -1], [1, 0]].
    const swapStep = () => {
      const t = a;
      a = c;
      c = t;
      b = -b;
      if (transformation) {
        U = matmul(U, [
          [0n, -1n],
          [1n, 0n],
        ]);
      }
    };

    redbu();
    let cmp = a < c ? -1 : a > c ? 1 : 0;
    while (cmp > 0) {
      swapStep();
      redbu();
      cmp = a < c ? -1 : a > c ? 1 : 0;
    }
    if (cmp === 0 && b < 0n) {
      swapStep();
    }
    const result = new BinaryQF(a, b, c);
    return transformation ? [result, U] : result;
  }

  /**
   * Reduce an indefinite, non-reduced form by repeated application of rho.
   *
   * @see Reference: sage/quadratic_forms/binary_qf.py:_reduce_indef (line 766)
   * @see Reference: pari/src/basemath/Qfb.c:qfr_redsl2_basecase (line 640)
   */
  private _reduce_indef(transformation: false): BinaryQF;
  private _reduce_indef(transformation: true): [BinaryQF, Matrix2];
  private _reduce_indef(transformation: boolean): BinaryQF | [BinaryQF, Matrix2] {
    let U: Matrix2 = IDENTITY2;
    // Sage uses the 53-bit real square root sqrt(D); floor((sqrt(D)+b)/(2|c|))
    // equals floor((isqrt(D)+b)/(2|c|)) because the numerator bound is integral,
    // so the exact integer square root gives identical results.
    const d = isqrt(this.discriminant());
    let Q: BinaryQF = this;
    while (!Q.is_reduced()) {
      const a = Q.a;
      const b = Q.b;
      const c = Q.c;
      const cabs = abs(c);
      if (cabs !== 0n) {
        // rho(f) as defined in [BUVO2007] p. 112 equation (6.12)
        const sign = c > 0n ? 1n : -1n;
        const s =
          cabs >= d ? sign * floorDiv(cabs + b, 2n * cabs) : sign * floorDiv(d + b, 2n * cabs);
        if (transformation) {
          U = matmul(U, [
            [0n, -1n],
            [1n, s],
          ]);
        }
        Q = new BinaryQF(c, -b + 2n * s * c, c * s * s - b * s + a);
      } else {
        if (b < 0n) {
          Q = new BinaryQF(a, -b, c);
          if (transformation) {
            U = matmul(U, [
              [1n, 0n],
              [0n, -1n],
            ]);
          }
        } else {
          let [q, r] = quoRem(a, b);
          if (2n * r > b) {
            [q, r] = quoRem(a, -b);
            q = -q;
          }
          if (transformation) {
            U = matmul(U, [
              [1n, 0n],
              [-q, 1n],
            ]);
          }
          Q = new BinaryQF(r, b, c);
        }
      }
    }
    return transformation ? [Q, U] : Q;
  }

  inverse(): BinaryQF {
    return new BinaryQF(this.a, -this.b, this.c);
  }

  /**
   * Compose two binary quadratic forms (Gaussian composition).
   *
   * Implements PARI's qfbcompraw algorithm from reference/pari/src/basemath/Qfb.c
   * The result is NOT reduced.
   *
   * @see Reference: pari/src/basemath/Qfb.c:qfb_comp
   */
  compose(other: BinaryQF): BinaryQF {
    const D = this.discriminant();
    if (other.discriminant() !== D) {
      throw new ValueError('forms must have the same discriminant');
    }

    // PARI's qfb_comp delegates to qfb_sqr only when both operands are the same
    // object (`if (x == y)`); Sage's `Q * Q` hits that path because __pari__()
    // caches the converted GEN on the form. Two distinct but equal forms take
    // the general path, so we compare by identity, not by value.
    if (this === other) {
      return this._square();
    }

    const { a: a1, b: b1, c: c1 } = this;
    const { a: a2, b: b2, c: c2 } = other;

    // n = (b2 - b1) / 2
    const n = (b2 - b1) / 2n;

    let v1 = a1;
    let v2 = a2;
    let c = c2;

    // d = gcd(v2, v1), with y1 such that y1*v2 + ?*v1 = d
    const [d, y1] = xgcd(v2, v1);

    let m: bigint;
    if (d === 1n) {
      m = y1 * n;
    } else {
      const s = b2 - n; // = (b1 + b2) / 2
      // d1 = gcd(s, d), with x2*s + y2*d = d1
      const [d1, x2, y2] = xgcd(s, d);

      if (d1 !== 1n) {
        v1 = v1 / d1;
        v2 = v2 / d1;
        // v1 *= gcd(c, gcd(c1, gcd(d1, n)))
        v1 = v1 * gcd(c, gcd(c1, gcd(d1, n)));
        c = c * d1;
      }
      m = y1 * y2 * n + c2 * x2;
    }

    m = -m;

    // r = m mod v1 (ensure positive remainder)
    let r = v1 !== 0n ? m % v1 : 0n;
    if (r < 0n) r += v1 < 0n ? -v1 : v1;

    const p1 = r * v2;
    const c3 = c + r * (b2 + p1);

    const newA = v1 * v2;
    const newB = b2 + 2n * p1;
    const newC = c3 / v1;

    return new BinaryQF(newA, newB, newC);
  }

  /**
   * Square this form (optimized composition with self).
   * @see Reference: pari/src/basemath/Qfb.c:qfb_sqr
   */
  private _square(): BinaryQF {
    const { a, b } = this;
    let c = this.c;

    // d1 = gcd(b, a), with x2 such that x2*b + ?*a = d1
    const [d1, x2] = xgcd(b, a);

    let m = c * x2;

    let v1: bigint;
    let v2: bigint;
    if (d1 === 1n) {
      v1 = a;
      v2 = a;
    } else {
      v1 = a / d1;
      v2 = v1 * gcd(d1, c); // = v1 iff this form is primitive
      c = c * d1;
    }

    m = -m;

    // r = m mod v2 (PARI's modii: nonnegative representative)
    let r = v2 !== 0n ? m % v2 : 0n;
    if (r < 0n) r += abs(v2);

    const p1 = r * v1;
    const c3 = c + r * (b + p1);

    const newA = v1 * v2;
    const newB = b + 2n * p1;
    const newC = c3 / v2;

    return new BinaryQF(newA, newB, newC);
  }

  /**
   * Return whether this form is equivalent to `other`.
   *
   * @see Reference: sage/quadratic_forms/binary_qf.py:is_equivalent (line 1296)
   */
  is_equivalent(other: BinaryQF, options?: { proper?: boolean }): boolean {
    const proper = options?.proper ?? true;
    if (this.discriminant() !== other.discriminant()) return false;

    if (this.is_indefinite()) {
      let selfred = this.reduced_form();
      let otherred = other.reduced_form();

      if (this.is_reducible()) {
        // Square discriminant: make sure we terminate in a form with c = 0.
        while (selfred.c !== 0n) selfred = selfred._Rho();
        while (otherred.c !== 0n) otherred = otherred._Rho();
        const b = selfred.b;
        const a = selfred.a;
        const ao = otherred.a;
        // p. 359 of Conway-Sloane [CS1999], but `2b` there is `b` here
        const isProperlyEquiv = (a - ao) % b === 0n;
        if (proper) return isProperlyEquiv;
        const g = gcd(a, b);
        return isProperlyEquiv || (gcd(ao, b) === g && (a * ao - g * g) % (b * g) === 0n);
      }

      const properCycle = otherred.cycle({ proper: true });
      const isProp = properCycle.some((f) => f.equals(selfred));
      if (proper || isProp) return isProp;
      // Note that our definition of improper equivalence differs from that of
      // Buchmann and Vollmer: their action is det(f) * q(f(x, y)), ours q(f(x, y)).
      const selfimp = new BinaryQF(selfred.c, selfred.b, selfred.a);
      return properCycle.some((f) => f.equals(selfimp));
    }

    // Definite forms.
    if (this.is_positive_definite() && !other.is_positive_definite()) return false;
    if (this.is_negative_definite() && !other.is_negative_definite()) return false;
    const Q1 = this.reduced_form();
    const Q2 = other.reduced_form();
    if (Q1.equals(Q2)) return true;
    if (!proper) {
      const Q1e = new BinaryQF(this.c, this.b, this.a).reduced_form();
      return Q1e.equals(Q2);
    }
    return false;
  }

  /**
   * Return the cycle of reduced forms to which this form belongs.
   *
   * With `{proper: true}` return the proper cycle, i.e. all reduced forms
   * properly equivalent to this one (Prop. 6.10.5 in [BUVO2007]).
   *
   * @see Reference: sage/quadratic_forms/binary_qf.py:cycle (line 1042)
   */
  cycle(options?: { proper?: boolean }): BinaryQF[] {
    const proper = options?.proper ?? false;
    if (!(this.is_indefinite() && this.is_reduced()))
      throw new ValueError(`${this} must be indefinite and reduced`);
    if (this.is_reducible())
      throw new NotImplementedError(
        'computation of cycles is only implemented for non-square discriminants'
      );
    if (proper) {
      // Prop 6.10.5 in Buchmann-Vollmer
      let C = this.cycle();
      if (C.length % 2 === 1) C = C.concat(C);
      const result = [...C];
      for (let i = 0; i < Math.floor(result.length / 2); i++) {
        result[2 * i + 1] = result[2 * i + 1]!._Tau();
      }
      return result;
    }
    if (this._cycle_list === undefined) {
      const result: BinaryQF[] = [this];
      let Q = this._RhoTau();
      while (!this.equals(Q)) {
        result.push(Q);
        Q = Q._RhoTau();
      }
      this._cycle_list = result;
    }
    return this._cycle_list;
  }

  /**
   * The `s` used by the Rho operator: `sign(c) * floor((max(|c|, sqrt(D)) + b) / (2|c|))`.
   */
  private _rho_s(): bigint {
    const d = isqrt(this.discriminant());
    const { b, c } = this;
    const cabs = abs(c);
    const sign = c > 0n ? 1n : -1n;
    return cabs >= d ? sign * floorDiv(cabs + b, 2n * cabs) : sign * floorDiv(d + b, 2n * cabs);
  }

  /**
   * Apply the Rho and Tau operators to this form.
   *
   * @see Reference: sage/quadratic_forms/binary_qf.py:_RhoTau (line 971)
   */
  private _RhoTau(): BinaryQF {
    const { a, b, c } = this;
    const s = this._rho_s();
    return new BinaryQF(-c, -b + 2n * s * c, -(a - b * s + c * s * s));
  }

  /**
   * Apply the Rho operator to this form.
   *
   * @see Reference: sage/quadratic_forms/binary_qf.py:_Rho (line 1004)
   */
  private _Rho(): BinaryQF {
    const { a, b, c } = this;
    const s = this._rho_s();
    return new BinaryQF(c, -b + 2n * s * c, a - b * s + c * s * s);
  }

  /**
   * Apply the Tau operator to this form: `(a, b, c) -> (-a, b, -c)`.
   *
   * @see Reference: sage/quadratic_forms/binary_qf.py:_Tau (line 1027)
   */
  private _Tau(): BinaryQF {
    return new BinaryQF(-this.a, this.b, -this.c);
  }

  neg(): BinaryQF {
    return new BinaryQF(-this.a, -this.b, -this.c);
  }
  add(other: BinaryQF): BinaryQF {
    return new BinaryQF(this.a + other.a, this.b + other.b, this.c + other.c);
  }
  sub(other: BinaryQF): BinaryQF {
    return new BinaryQF(this.a - other.a, this.b - other.b, this.c - other.c);
  }

  matrix_action_right(M: [[IntegerLike, IntegerLike], [IntegerLike, IntegerLike]]): BinaryQF {
    const p = toBigInt(M[0][0]);
    const q = toBigInt(M[0][1]);
    const r = toBigInt(M[1][0]);
    const s = toBigInt(M[1][1]);
    const A = this.evaluate(p, r);
    const C = this.evaluate(q, s);
    const B = this.evaluate(p + q, r + s) - A - C;
    return new BinaryQF(A, B, C);
  }

  matrix_action_left(M: [[IntegerLike, IntegerLike], [IntegerLike, IntegerLike]]): BinaryQF {
    const a = toBigInt(M[0][0]);
    const b = toBigInt(M[0][1]);
    const c = toBigInt(M[1][0]);
    const d = toBigInt(M[1][1]);
    const A = this.evaluate(a, b);
    const C = this.evaluate(c, d);
    const B = this.evaluate(a + c, b + d) - A - C;
    return new BinaryQF(A, B, C);
  }

  equals(other: BinaryQF): boolean {
    return this.a === other.a && this.b === other.b && this.c === other.c;
  }

  compare(other: BinaryQF): -1 | 0 | 1 {
    if (this.a < other.a) return -1;
    if (this.a > other.a) return 1;
    if (this.b < other.b) return -1;
    if (this.b > other.b) return 1;
    if (this.c < other.c) return -1;
    if (this.c > other.c) return 1;
    return 0;
  }

  toString(): string {
    const { a, b, c } = this;
    if (a === 0n && b === 0n && c === 0n) return '0';
    const terms: string[] = [];
    if (a !== 0n) terms.push(a === 1n ? 'x^2' : a === -1n ? '-x^2' : `${a}*x^2`);
    if (b !== 0n) {
      const pre = terms.length > 0;
      if (b === 1n) terms.push(pre ? '+ x*y' : 'x*y');
      else if (b === -1n) terms.push(pre ? '- x*y' : '-x*y');
      else if (b > 0n) terms.push(pre ? `+ ${b}*x*y` : `${b}*x*y`);
      else terms.push(pre ? `- ${-b}*x*y` : `${b}*x*y`);
    }
    if (c !== 0n) {
      const pre = terms.length > 0;
      if (c === 1n) terms.push(pre ? '+ y^2' : 'y^2');
      else if (c === -1n) terms.push(pre ? '- y^2' : '-y^2');
      else if (c > 0n) terms.push(pre ? `+ ${c}*y^2` : `${c}*y^2`);
      else terms.push(pre ? `- ${-c}*y^2` : `${c}*y^2`);
    }
    return terms.join(' ');
  }
}

function isSquare(n: bigint): boolean {
  if (n < 0n) return false;
  const s = isqrt(n);
  return s * s === n;
}

export function BinaryQF_reduced_representatives(
  D: IntegerLike,
  options?: { primitive_only?: boolean; proper?: boolean }
): BinaryQF[] {
  const _D = toBigInt(D);
  const primitive_only = options?.primitive_only ?? false;
  const proper = options?.proper ?? true;
  const D4 = ((_D % 4n) + 4n) % 4n;
  if (D4 !== 0n && D4 !== 1n) throw new ValueError(`${_D} is not a discriminant`);
  const formList: BinaryQF[] = [];

  if (_D > 0n) {
    if (isSquare(_D)) {
      // b = sqrt(D), c = 0 and -b/2 < a <= b/2
      const b = isqrt(_D);
      for (let a = floorDiv(-b, 2n) + 1n; a <= floorDiv(b, 2n); a++) {
        if (!primitive_only || gcd([a, b, 0n]) === 1n) formList.push(new BinaryQF(a, b, 0n));
      }
    } else {
      // We follow the description of Buchmann/Vollmer 6.7.1. They enumerate all
      // reduced forms; we only want representatives.
      const sqrtD = isqrt(_D);
      for (let b = 1n; b <= sqrtD; b++) {
        if ((_D - b) % 2n !== 0n) continue;
        const A = (_D - b * b) / 4n;
        // Low_a = ceil((sqrt(D) - b) / 2); D is not a square, so sqrt(D) is
        // irrational and the ceiling is floor((isqrt(D) - b)/2) + 1.
        const lowA = floorDiv(sqrtD - b, 2n) + 1n;
        const highA = isqrt(A);
        for (let a = lowA; a <= highA; a++) {
          if (a === 0n) continue;
          if (A % a !== 0n) continue;
          const c = -(A / a);
          if (!primitive_only || gcd([a, b, c]) === 1n) {
            formList.push(new BinaryQF(a, b, c));
            formList.push(new BinaryQF(-a, b, -c));
            if (abs(a) !== abs(c)) {
              formList.push(new BinaryQF(c, b, a));
              formList.push(new BinaryQF(-c, b, -a));
            }
          }
        }
      }
    }
  } else {
    const bound = isqrt(-_D / 3n) + 1n;
    for (let a = 1n; a <= bound; a++) {
      const a4 = 4n * a;
      const s = _D + a * a4;
      let w = s > 0n ? isqrt(s - 1n) + 1n : 0n;
      if (((w % 2n) + 2n) % 2n !== ((_D % 2n) + 2n) % 2n) w += 1n;
      for (let b = w; b <= a; b += 2n) {
        const t = b * b - _D;
        if (t % a4 === 0n) {
          const c = t / a4;
          if (!primitive_only || gcd([a, b, c]) === 1n) {
            if (c > a && a > b && b > 0n) formList.push(new BinaryQF(a, -b, c));
            formList.push(new BinaryQF(a, b, c));
          }
        }
      }
    }
  }

  if (!proper || _D > 0n) {
    const filtered: BinaryQF[] = [];
    for (const q of formList) {
      let dominated = false;
      for (const q1 of filtered)
        if (q.is_equivalent(q1, { proper })) {
          dominated = true;
          break;
        }
      if (!dominated) filtered.push(q);
    }
    formList.length = 0;
    formList.push(...filtered);
  }
  formList.sort((a, b) => a.compare(b));
  return formList;
}

export function class_number(D: IntegerLike): bigint {
  const _D = toBigInt(D);
  const reps = BinaryQF_reduced_representatives(_D, { primitive_only: true, proper: true });
  return BigInt(reps.length);
}

export const BinaryQuadraticForm = BinaryQF;
