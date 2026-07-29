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

import {
  type QfbForm,
  mkqfb,
  qfbcompraw,
  qfbcornacchia,
  qfbred,
  qfbredsl2,
  qfbsolve,
} from '@sagemath-ts/parigp-ts';
import { divisors, gcd, isqrt, xgcd } from '../arith/misc.js';
import {
  ArithmeticError,
  NotImplementedError,
  PariError,
  ValueError,
  ZeroDivisionError,
} from '../errors.js';
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
 * Convert a form to PARI's `t_QFB`.
 *
 * Sage sends the form to PARI through `_pari_init_` (`binary_qf.py:158-180`),
 * i.e. through PARI's `Qfb(a, b, c)` constructor. We use `mkqfb` (the unchecked
 * constructor) here and keep the domain checks in the callers, so that the
 * error raised for a negative definite or reducible form is Sage's own message
 * rather than a PARI one.
 */
function toPari(f: BinaryQF): QfbForm {
  return mkqfb(f.a, f.b, f.c, f.discriminant());
}

/**
 * Reproduce the domain validation PARI's `Qfb(a, b, c)` constructor performs.
 *
 * `reference/pari/src/basemath/Qfb.c:174-176`:
 *
 * ```c
 * if (signe(D) < 0)
 * { if (signe(a) < 0) pari_err_IMPL("negative definite t_QFB"); }
 * else if (Z_issquare(D)) pari_err_DOMAIN("Qfb","issquare(disc)","=", gen_1, q);
 * ```
 *
 * Sage reaches this on every conversion through `_pari_init_`
 * (`binary_qf.py:158-180`), so a `BinaryQF` that PARI refuses to build must
 * raise here rather than silently producing a value.
 */
function checkPariQfbDomain(f: BinaryQF): void {
  const D = f.discriminant();
  if (D < 0n) {
    if (f.a < 0n) {
      throw new PariError('Qfb: sorry, negative definite t_QFB is not yet implemented');
    }
  } else if (D >= 0n && isqrt(D) * isqrt(D) === D) {
    throw new PariError('Qfb: domain error in Qfb: issquare(disc) = 1');
  }
}

/** Convert a PARI `t_QFB` back to a {@link BinaryQF}. */
function fromPari(q: QfbForm): BinaryQF {
  return new BinaryQF(q.a, q.b, q.c);
}

/** Convert PARI's row-major 2x2 base change to our {@link Matrix2}. */
function toMatrix2(U: bigint[][]): Matrix2 {
  return [
    [U[0]![0]!, U[0]![1]!],
    [U[1]![0]!, U[1]![1]!],
  ];
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
   * Exactly as in Sage, the work is done by PARI's `qfbred`/`qfbredsl2` unless
   * the discriminant is a square, in which case Sage's own `_reduce_indef` is
   * used (PARI has no `t_QFB` of square discriminant).
   *
   * @param options.algorithm `'default'` (Sage picks), `'pari'` or `'sage'`
   *
   * @see Reference: sage/quadratic_forms/binary_qf.py:reduced_form (line 831)
   * @see Reference: pari/src/basemath/Qfb.c:qfbred (line 991), qfbredsl2 (line 889)
   */
  reduced_form(options?: {
    transformation?: false;
    algorithm?: 'default' | 'pari' | 'sage';
  }): BinaryQF;
  reduced_form(options: {
    transformation: true;
    algorithm?: 'default' | 'pari' | 'sage';
  }): [BinaryQF, Matrix2];
  reduced_form(options?: {
    transformation?: boolean;
    algorithm?: 'default' | 'pari' | 'sage';
  }): BinaryQF | [BinaryQF, Matrix2] {
    const transformation = options?.transformation ?? false;
    let algorithm = options?.algorithm ?? 'default';
    if (this.is_reduced()) {
      return transformation ? [this, IDENTITY2] : this;
    }

    // binary_qf.py:947-948
    if (algorithm === 'default') algorithm = this.is_reducible() ? 'sage' : 'pari';

    if (algorithm === 'sage') {
      // binary_qf.py:950-955
      if (this.discriminant() <= 0n) {
        throw new NotImplementedError(
          'reduction of definite binary quadratic forms is not implemented in Sage'
        );
      }
      return transformation ? this._reduce_indef(true) : this._reduce_indef(false);
    }

    if (algorithm !== 'pari') {
      throw new ValueError(
        `unknown implementation for binary quadratic form reduction: ${algorithm}`
      );
    }

    // binary_qf.py:957-974
    if (this.is_negative_definite()) {
      // PARI does not support negative definite forms; Sage reduces
      // (-self)*M with M = diag(-1, 1) instead and conjugates the result back.
      const negForm = new BinaryQF(-this.a, this.b, -this.c);
      if (transformation) {
        const [reduced, M] = negForm.reduced_form({ transformation: true, algorithm });
        // M_diag * M * M_diag with M_diag = diag(-1, 1)
        const conjugated: Matrix2 = [
          [M[0][0], -M[0][1]],
          [-M[1][0], M[1][1]],
        ];
        return [new BinaryQF(-reduced.a, reduced.b, -reduced.c), conjugated];
      }
      const reduced = negForm.reduced_form({ algorithm });
      return new BinaryQF(-reduced.a, reduced.b, -reduced.c);
    }

    if (this.is_reducible()) {
      throw new NotImplementedError('reducible forms are not supported using PARI');
    }

    if (transformation) {
      const { Q, U } = qfbredsl2(toPari(this));
      return [fromPari(Q), toMatrix2(U)];
    }
    return fromPari(qfbred(toPari(this)));
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
   * Compose two binary quadratic forms (Gaussian composition). The result is
   * NOT reduced.
   *
   * This is Sage's `__mul__` for two forms, which delegates to PARI's
   * `qfbcompraw` (`binary_qf.py:260`).
   *
   * Note on PARI's `qfb_comp` (`Qfb.c:1038-1071`): it takes the dedicated
   * squaring shortcut `qfb_sqr` only when the two operands are the *same GEN*
   * (`if (x == y)`). Sage always converts both arguments separately
   * (`self.__pari__()` and cypari2's `objtogen(right)` produce distinct GENs),
   * so even `Q * Q` goes through the general composition; we mirror that by
   * always passing two distinct `t_QFB` values.
   *
   * @see Reference: sage/quadratic_forms/binary_qf.py:__mul__ (line 224)
   * @see Reference: pari/src/basemath/Qfb.c:qfbcompraw (line 1165)
   */
  compose(other: BinaryQF): BinaryQF {
    // `binary_qf.py:224` `__mul__` converts BOTH operands through
    // `_pari_init_` -> `Qfb(a, b, c)` before calling `qfbcompraw`, so PARI's
    // constructor validation runs first and rejects negative definite operands
    // and square (including zero) discriminants.
    checkPariQfbDomain(this);
    checkPariQfbDomain(other);
    const D = this.discriminant();
    if (other.discriminant() !== D) {
      throw new ValueError('forms must have the same discriminant');
    }
    return fromPari(qfbcompraw(toPari(this), toPari(other)));
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

  /**
   * Solve `Q(x, y) = n` in integers.
   *
   * ALGORITHM: PARI's `qfbsolve` (or `qfbcornacchia`), exactly as in Sage.
   * Forms of square discriminant are not representable as a PARI `t_QFB`; for
   * those Sage has its own elementary algorithm, which we port as well.
   *
   * @param n the value to represent
   * @param options.algorithm `'general'` (default) or `'cornacchia'`
   * @param options._flag PARI's `qfbsolve` flag: 1, 2 (default) or 3
   * @param options.factorization optional known factorization of `n`
   * @returns for `_flag = 2` a pair `[x, y]` or `null`; otherwise the list of
   *          all solutions found
   *
   * @see Reference: sage/quadratic_forms/binary_qf.py:solve_integer (line 1608)
   * @see Reference: pari/src/basemath/Qfb.c:qfbsolve (line 1987)
   */
  solve_integer(
    n: IntegerLike,
    options?: {
      algorithm?: 'general' | 'cornacchia';
      _flag?: 2;
      factorization?: Array<[IntegerLike, IntegerLike]>;
    }
  ): [bigint, bigint] | null;
  solve_integer(
    n: IntegerLike,
    options: {
      algorithm?: 'general' | 'cornacchia';
      _flag: 1 | 3;
      factorization?: Array<[IntegerLike, IntegerLike]>;
    }
  ): Array<[bigint, bigint]>;
  solve_integer(
    n: IntegerLike,
    options?: {
      algorithm?: 'general' | 'cornacchia';
      _flag?: 1 | 2 | 3;
      factorization?: Array<[IntegerLike, IntegerLike]>;
    }
  ): [bigint, bigint] | null | Array<[bigint, bigint]> {
    const algorithm = options?.algorithm ?? 'general';
    const _flag = options?._flag ?? 2;
    const _n = toBigInt(n);

    // binary_qf.py:1748-1749. NB Sage's recursion drops the keyword arguments,
    // so a negative definite form always answers with the `_flag = 2` shape;
    // we reproduce that verbatim.
    if (this.is_negative_definite()) {
      return this.neg().solve_integer(-_n) as [bigint, bigint] | null;
    }

    // binary_qf.py:1751-1791: square discriminant, not supported by PARI.
    if (this.is_reducible()) {
      let M: Matrix2;
      if (this.a !== 0n) {
        // https://math.stackexchange.com/a/980075
        const w = isqrt(this.discriminant());
        // r = (-b +- w) / (2a) as an exact rational p/q with q > 0
        let p = -this.b + (w !== this.b ? w : -w);
        let q = 2n * this.a;
        const g0 = gcd(p, q);
        if (g0 !== 0n) {
          p /= g0;
          q /= g0;
        }
        if (q < 0n) {
          p = -p;
          q = -q;
        }
        const [, u, v] = xgcd(p, q);
        M = [
          [v, p],
          [-u, q],
        ];
      } else if (this.c !== 0n) {
        M = [
          [0n, 1n],
          [1n, 0n],
        ];
      } else {
        M = IDENTITY2;
      }
      const Q = this.matrix_action_right(M);
      if (Q.c !== 0n) throw new ArithmeticError('solve_integer: expected c == 0');

      if (Q.b === 0n) {
        // at this point, Q = a*x^2
        if (Q.a === 0n) {
          // Upstream (`binary_qf.py:1775-1780`) guards with
          // `if Q._a.divides(n) and (n // Q._a).is_square()`.  `ZZ(0).divides(n)`
          // is False for every n != 0, so Sage falls through to `return None`;
          // only n == 0 reaches `0 // 0`.
          if (_n !== 0n) return null;
          throw new ZeroDivisionError('Integer division by zero');
        }
        const quo = _n / Q.a;
        if (_n % Q.a === 0n && quo >= 0n && isqrt(quo) * isqrt(quo) === quo) {
          const x = isqrt(quo);
          return [M[0][0] * x, M[1][0] * x];
        }
        return null;
      }

      // at this point, Q = a*x^2 + b*x*y
      if (_n === 0n) return [M[0][1], M[1][1]];
      for (const x of divisors(_n)) {
        const yNum = _n / x - Q.a * x;
        if (yNum % Q.b === 0n) {
          const y = yNum / Q.b;
          return [M[0][0] * x + M[0][1] * y, M[1][0] * x + M[1][1] * y];
        }
      }
      return null;
    }

    if (algorithm === 'cornacchia') {
      if (!(this.a === 1n && this.b === 0n && this.c > 0n)) {
        throw new ValueError("Cornacchia's algorithm requires a=1 and b=0 and c>0");
      }
      const sol = qfbcornacchia(this.c, _n);
      return sol ? [sol[0], sol[1]] : null;
    }

    if (algorithm !== 'general') {
      throw new ValueError(`algorithm '${algorithm}' is not a valid algorithm`);
    }

    const fa = options?.factorization
      ? options.factorization.map(([p, e]) => [toBigInt(p), toBigInt(e)] as [bigint, bigint])
      : null;
    const sol = qfbsolve(toPari(this), _n, _flag, fa);
    if (_flag === 2) {
      const s = sol as bigint[];
      return s.length ? [s[0]!, s[1]!] : null;
    }
    return sol as Array<[bigint, bigint]>;
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
