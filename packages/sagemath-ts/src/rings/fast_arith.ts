/**
 * Bounded arithmetic selected by sage.arith.misc's gcd/inverse factories.
 * @see Reference: sage/rings/fast_arith.pyx
 * @see Deviation: Bounded native arithmetic and Cython error sentinels
 */
import { ArithmeticError, OverflowError } from '../errors.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';

/** Cython converts through a 64-bit C long before narrowing to int. */
function nativeInteger(value: IntegerLike, int: boolean): bigint {
  const n = toBigInt(value);
  if (n < -(1n << 63n) || n >= 1n << 63n) {
    throw new OverflowError('Python int too large to convert to C long');
  }
  if (int && (n < -(1n << 31n) || n >= 1n << 31n)) {
    throw new OverflowError('value too large to convert to int');
  }
  return n;
}

function gcdNative(a: bigint, b: bigint): bigint {
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;
  while (b !== 0n) [a, b] = [b, a % b];
  return a;
}

/** c_xgcd_int/c_xgcd_longlong, retaining the original coefficient updates. */
function xgcdNative(a: bigint, b: bigint): [bigint, bigint, bigint] {
  const psign = a < 0n ? -1n : 1n;
  const qsign = b < 0n ? -1n : 1n;
  if (a === 0n) return [b * qsign, 0n, qsign];
  if (b === 0n) return [a * psign, psign, 0n];
  a *= psign;
  b *= qsign;
  let p = 1n,
    q = 0n,
    r = 0n,
    s = 1n;
  while (b !== 0n) {
    const quot = a / b;
    [a, b] = [b, a % b];
    [p, r] = [r, p - quot * r];
    [q, s] = [s, q - quot * s];
  }
  return [a, p * psign, q * qsign];
}

/** ARM native remainder by zero retains the dividend; see the module deviation. */
function nativeRemainder(a: bigint, m: bigint): bigint {
  return m === 0n ? a : a % m;
}

/** A -1 return collides with the original Cython `except -1` sentinel. */
function nativeReturn(value: bigint): bigint {
  if (value === -1n) {
    const error = new Error('error return without exception set');
    error.name = 'SystemError';
    throw error;
  }
  return value;
}

function inverseNative(a: bigint, m: bigint): bigint {
  const [g, coefficient] = xgcdNative(a, m);
  if (g !== 1n) {
    throw new ArithmeticError(`The inverse of ${a} modulo ${m} is not defined.`);
  }
  const s = nativeRemainder(coefficient, m);
  return s < 0n ? s + m : s;
}

export class arith_int {
  gcd_int(a: IntegerLike, b: IntegerLike): bigint {
    return gcdNative(nativeInteger(a, true), nativeInteger(b, true));
  }
  inverse_mod_int(a: IntegerLike, m: IntegerLike): bigint {
    const value = nativeInteger(a, true);
    const modulus = nativeInteger(m, true);
    return nativeReturn(
      value === 1n || modulus <= 1n
        ? nativeRemainder(value, modulus)
        : inverseNative(value, modulus)
    );
  }
}

export class arith_llong {
  gcd_longlong(a: IntegerLike, b: IntegerLike): bigint {
    return gcdNative(nativeInteger(a, false), nativeInteger(b, false));
  }
  inverse_mod_longlong(a: IntegerLike, m: IntegerLike): bigint {
    return nativeReturn(inverseNative(nativeInteger(a, false), nativeInteger(m, false)));
  }
}
