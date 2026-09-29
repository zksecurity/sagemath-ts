/** Bounded integer factorization, ported from sage/rings/factorint.pyx. */
import { type Factorization, trial_division } from '../arith/misc.js';
import { OverflowError } from '../errors.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';
import { mpz_remove } from '../types/gmp.js';

/**
 * Partial factorization using trial division, retaining an unresolved cofactor.
 * @see Deviation: Bounded factorization for elliptic hybrid orders
 */
export function factor_trial_division(
  m: IntegerLike,
  limit: IntegerLike = (1n << 63n) - 1n
): Factorization {
  const bound = toBigInt(limit);
  if (bound < -(1n << 63n) || bound >= 1n << 63n)
    throw new OverflowError('Python int too large to convert to C long');
  let n = toBigInt(m),
    p = 2n;
  const factors: Factorization = [];
  if (n <= 0n) {
    factors.push([-1n, 1n]);
    n = -n;
  }
  while (n !== 1n) {
    // mpz_get_ui in the source retains the low word of the previous divisor.
    p = trial_division(n, bound, BigInt.asUintN(64, p));
    const [e, rest] = mpz_remove(n, p);
    n = rest;
    factors.push([p, e]);
  }
  return factors;
}
