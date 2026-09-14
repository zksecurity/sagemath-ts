/**
 * Integer-domain totients and divisor counts from PARI basemath/arith2.c.
 * @see Deviation: Native integer totients, divisor counts and valuation routing
 */
import { Z_factor, factoru } from './ifactor.js';

import { ZV_prod } from './ZV.js';

/** Native eulerphi accepts signed integers and returns 2 at zero. */
export function eulerphi(n: bigint): bigint {
  if (n === 0n) return 2n;
  const value = n < 0n ? -n : n;
  const factors = Z_factor(value);
  if (value < 1n << 64n) {
    // eulerphiu_fact: word arithmetic, with the original power-of-two shortcut.
    let result = 1n;
    for (const [p, e] of factors) {
      if (e === 0n) continue;
      if (p === 2n) {
        if (e > 1n) result <<= e - 1n;
      } else {
        result *= p - 1n;
        if (e > 1n) result *= p ** (e - 1n);
      }
    }
    return result;
  }
  return ZV_prod(factors.map(([p, e]) => (p - 1n) * (e === 1n ? 1n : p ** (e - 1n))));
}

/** Native numdiv rejects zero and ignores an integer's sign. */
export function numdiv(n: bigint): bigint {
  if (n === 0n) {
    const error = new Error('domain error in numdiv: argument = 0');
    error.name = 'PariError';
    throw error;
  }
  const value = n < 0n ? -n : n;
  const exponents = Z_factor(value).map(([, exponent]) => exponent + 1n);
  return ZV_prod(exponents);
}

/** Native unsigned-word totient (arith2.c:eulerphiu/eulerphiu_fact).
 * The internal input contract is 0 <= n < 2^64.
 */
export function eulerphiu(n: bigint): bigint {
  if (n === 0n) return 2n;
  let result = 1n;
  for (const [p, e] of factoru(n)) {
    if (e === 0n) continue;
    if (p === 2n) {
      if (e > 1n) result <<= e - 1n;
    } else {
      result *= p - 1n;
      if (e > 1n) result *= p ** (e - 1n);
    }
  }
  return result;
}
