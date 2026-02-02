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

/**
 * A binary quadratic form over the integers.
 * Represents the form `a*x^2 + b*x*y + c*y^2` for integer coefficients (a, b, c).
 */
export class BinaryQF {
  readonly a: bigint;
  readonly b: bigint;
  readonly c: bigint;

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

  reduced_form(options?: { transformation?: false }): BinaryQF;
  reduced_form(options: { transformation: true }): [BinaryQF, [[bigint, bigint], [bigint, bigint]]];
  reduced_form(options?: { transformation?: boolean }):
    | BinaryQF
    | [BinaryQF, [[bigint, bigint], [bigint, bigint]]] {
    const transformation = options?.transformation ?? false;
    if (this.is_reduced()) {
      return transformation
        ? [
            this,
            [
              [1n, 0n],
              [0n, 1n],
            ],
          ]
        : this;
    }
    const D = this.discriminant();
    if (D < 0n) {
      if (this.is_negative_definite()) {
        const negForm = new BinaryQF(-this.a, this.b, -this.c);
        if (transformation) {
          const [reduced, M] = negForm.reduced_form({ transformation: true });
          return [new BinaryQF(-reduced.a, reduced.b, -reduced.c), M];
        }
        const reduced = negForm.reduced_form();
        return new BinaryQF(-reduced.a, reduced.b, -reduced.c);
      }
      return transformation
        ? this._reduce_positive_definite(true)
        : this._reduce_positive_definite(false);
    }
    return transformation ? this._reduce_indefinite(true) : this._reduce_indefinite(false);
  }

  private _reduce_positive_definite(transformation: false): BinaryQF;
  private _reduce_positive_definite(
    transformation: true
  ): [BinaryQF, [[bigint, bigint], [bigint, bigint]]];
  private _reduce_positive_definite(
    transformation: boolean
  ): BinaryQF | [BinaryQF, [[bigint, bigint], [bigint, bigint]]] {
    let a = this.a;
    let b = this.b;
    let c = this.c;
    let U: [[bigint, bigint], [bigint, bigint]] = [
      [1n, 0n],
      [0n, 1n],
    ];

    while (true) {
      if (b < -a || b > a) {
        const twoA = 2n * a;
        let s = b / twoA;
        if (b >= 0n) {
          if (b % twoA > a) s += 1n;
        } else {
          if (-b % twoA > a) s -= 1n;
        }
        const newB = b - 2n * s * a;
        const newC = a * s * s - b * s + c;
        b = newB;
        c = newC;
        if (transformation)
          U = [
            [U[0][0] - s * U[0][1], U[0][1]],
            [U[1][0] - s * U[1][1], U[1][1]],
          ];
      }
      if (a > c) {
        const temp = a;
        a = c;
        c = temp;
        b = -b;
        if (transformation)
          U = [
            [U[0][1], -U[0][0]],
            [U[1][1], -U[1][0]],
          ];
      } else break;
    }
    if (a === c && b < 0n) {
      b = -b;
      if (transformation)
        U = [
          [U[0][0], -U[0][1]],
          [U[1][0], -U[1][1]],
        ];
    }
    const result = new BinaryQF(a, b, c);
    return transformation ? [result, U] : result;
  }

  private _reduce_indefinite(transformation: false): BinaryQF;
  private _reduce_indefinite(
    transformation: true
  ): [BinaryQF, [[bigint, bigint], [bigint, bigint]]];
  private _reduce_indefinite(
    transformation: boolean
  ): BinaryQF | [BinaryQF, [[bigint, bigint], [bigint, bigint]]] {
    let a = this.a;
    let b = this.b;
    let c = this.c;
    const D = this.discriminant();
    const sqrtD = isqrt(D);
    let U: [[bigint, bigint], [bigint, bigint]] = [
      [1n, 0n],
      [0n, 1n],
    ];

    const isReducedIndef = (a: bigint, b: bigint, c: bigint) =>
      (b > 0n && a * c < 0n && (a - c) * (a - c) < D) ||
      (a === 0n && -b < 2n * c && 2n * c <= b) ||
      (c === 0n && -b < 2n * a && 2n * a <= b);

    while (!isReducedIndef(a, b, c)) {
      const cAbs = c < 0n ? -c : c;
      if (c === 0n) {
        if (b < 0n) {
          b = -b;
          if (transformation)
            U = [
              [U[0][0], -U[0][1]],
              [U[1][0], -U[1][1]],
            ];
        } else {
          let q = a / b;
          if (2n * (a % b) > b) q += 1n;
          a = a - q * b;
          if (transformation)
            U = [
              [U[0][0] - q * U[0][1], U[0][1]],
              [U[1][0] - q * U[1][1], U[1][1]],
            ];
        }
      } else {
        const signC = c > 0n ? 1n : -1n;
        const s = signC * ((cAbs >= sqrtD ? cAbs + b : sqrtD + b) / (2n * cAbs));
        const newA = c;
        const newB = -b + 2n * s * c;
        const newC = c * s * s - b * s + a;
        a = newA;
        b = newB;
        c = newC;
        if (transformation)
          U = [
            [-U[0][1], U[0][0] + s * U[0][1]],
            [-U[1][1], U[1][0] + s * U[1][1]],
          ];
      }
    }
    const result = new BinaryQF(a, b, c);
    return transformation ? [result, U] : result;
  }

  inverse(): BinaryQF {
    return new BinaryQF(this.a, -this.b, this.c);
  }

  compose(other: BinaryQF): BinaryQF {
    const D = this.discriminant();
    if (other.discriminant() !== D) throw new ValueError('forms must have the same discriminant');
    const { a: a1, b: b1, c: c1 } = this;
    const { a: a2, b: b2 } = other;
    const bSum = b1 + b2;
    if (bSum % 2n !== 0n) throw new ValueError('invalid forms for composition');
    const bHalf = bSum / 2n;
    const [g, u] = xgcd(a1, a2);
    const [d, w, z] = xgcd(g, bHalf);
    const a = (a1 * a2) / (d * d);
    let b =
      b1 +
      (2n * a1 * u * w * (bHalf - b1)) / d +
      (2n * a1 * z * (c1 - (b1 * b1 - D) / (4n * a1))) / d;
    const twoA = 2n * a;
    b = ((b % twoA) + twoA) % twoA;
    if (b > a) b -= twoA;
    const c = (b * b - D) / (4n * a);
    return new BinaryQF(a, b, c);
  }

  is_equivalent(other: BinaryQF, options?: { proper?: boolean }): boolean {
    const proper = options?.proper ?? true;
    if (this.discriminant() !== other.discriminant()) return false;
    if (this.discriminant() < 0n) {
      const red1 = this.reduced_form();
      const red2 = other.reduced_form();
      if (red1.equals(red2)) return true;
      if (!proper) return new BinaryQF(red1.c, red1.b, red1.a).reduced_form().equals(red2);
      return false;
    }
    const red1 = this.reduced_form();
    const red2 = other.reduced_form();
    if (red1.equals(red2)) return true;
    const cycle = red2.cycle();
    for (const form of cycle) if (red1.equals(form)) return true;
    if (!proper) {
      const red1inv = new BinaryQF(red1.c, red1.b, red1.a);
      const r = red1inv.is_reduced() ? red1inv : red1inv.reduced_form();
      for (const form of cycle) if (r.equals(form)) return true;
    }
    return false;
  }

  cycle(): BinaryQF[] {
    if (!(this.is_indefinite() && this.is_reduced()))
      throw new ValueError('form must be indefinite and reduced');
    if (this.is_reducible())
      throw new NotImplementedError(
        'computation of cycles is only implemented for non-square discriminants'
      );
    const result: BinaryQF[] = [this];
    let Q = this._rhoTau();
    while (!this.equals(Q)) {
      result.push(Q);
      Q = Q._rhoTau();
    }
    return result;
  }

  private _rhoTau(): BinaryQF {
    const D = this.discriminant();
    const sqrtD = isqrt(D);
    const { a, b, c } = this;
    const cAbs = c < 0n ? -c : c;
    const signC = c > 0n ? 1n : -1n;
    const s = signC * ((cAbs >= sqrtD ? cAbs + b : sqrtD + b) / (2n * cAbs));
    return new BinaryQF(-c, -b + 2n * s * c, -(a - b * s + c * s * s));
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
      const sqrtD = isqrt(_D);
      for (let a = -sqrtD / 2n + 1n; a <= sqrtD / 2n; a++) {
        if (!primitive_only || gcd([a, sqrtD, 0n]) === 1n)
          formList.push(new BinaryQF(a, sqrtD, 0n));
      }
    } else {
      const sqrtD = isqrt(_D);
      for (let b = 1n; b <= sqrtD; b++) {
        if ((_D - b * b) % 2n !== 0n) continue;
        const A = (_D - b * b) / 4n;
        const lowA = (sqrtD - b + 1n) / 2n;
        const highA = isqrt(A);
        for (let a = lowA > 1n ? lowA : 1n; a <= highA; a++) {
          if (A % a !== 0n) continue;
          const c = -(A / a);
          if (!primitive_only || gcd([a, b, c]) === 1n) {
            formList.push(new BinaryQF(a, b, c));
            formList.push(new BinaryQF(-a, b, -c));
            if (a !== -c) {
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
