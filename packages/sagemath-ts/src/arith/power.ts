/** The Sage-element path of sage/arith/power.pyx and power.pxd. * @see Deviation: Number-field reduction and native array adapters
 */
import { type IntegerLike, toBigInt } from '../types/coercion.js';

type PowerElement<T> = {
  mul(other: T): T;
  inv(): T;
  parent(): { one(): T };
};

/** Exact integer powers using Sage's least-significant-set-bit schedule.
 * The internal caller supplies the port's Sage-element arithmetic protocol.
 */
export function generic_power<T extends PowerElement<T>>(a: T, exponent: IntegerLike): T {
  let n = toBigInt(exponent);
  if (n === 0n) return a.parent().one();
  if (n < 0n) {
    a = a.inv();
    n = -n;
  }
  return generic_power_pos(a, n);
}

/** Sage power.pyx generic_power_pos; the caller has handled zero and inversion. */
export function generic_power_pos<T extends { mul(other: T): T }>(a: T, n: bigint): T {
  if (n <= 0n) throw new RangeError('generic_power_pos requires a positive exponent');
  let apow = a;
  while ((n & 1n) === 0n) {
    apow = apow.mul(apow);
    n >>= 1n;
  }
  let result = apow;
  n >>= 1n;
  while (n !== 0n) {
    apow = apow.mul(apow);
    if (n & 1n) result = apow.mul(result);
    n >>= 1n;
  }
  return result;
}
