/** Fast direct square decompositions from sage/rings/sum_of_squares.pyx.
 * @see Deviation: Square decomposition input mapping and native bounds
 */
import { isqrt } from '../arith/misc.js';
import { OverflowError, ValueError } from '../errors.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';

/** Cython uint32 conversion on the original 64-bit target. */
function uint32(input: IntegerLike): bigint {
  const n = toBigInt(input);
  if (n < 0n) throw new OverflowError("can't convert negative value to uint32_t");
  if (n >= 1n << 64n) throw new OverflowError('Python int too large to convert to C unsigned long');
  if (n >= 1n << 32n) throw new OverflowError('value too large to convert to uint32_t');
  return n;
}

/**
 * Direct (factorization-free) search for `n = i^2 + j^2`, `i <= j`.
 *
 * Port of `two_squares_c` (`sage/rings/sum_of_squares.pyx:26`). Returns the
 * lexicographically smallest solution, or null when none exists. SageMath uses
 * this for every n < 2^32, so its answers -- not the Cornacchia-based ones --
 * are the reference values there.
 *
 * @see Reference: sage/rings/sum_of_squares.pyx:two_squares_c
 */
function two_squares_c(n: bigint): [bigint, bigint] | null {
  if (n === 0n) {
    return [0n, 0n];
  }

  // If n = 0 mod 4 then i and j must both be even: strip powers of 4 and
  // scale the solution back up at the end.
  let fac = 0n;
  while (n % 4n === 0n) {
    n >>= 2n;
    fac++;
  }

  // A sum of two squares is 0, 1 or 2 mod 4.
  if (n % 4n === 3n) {
    return null;
  }

  let i: bigint;
  let ii: bigint;
  let j = isqrt(n);
  let jj: bigint;

  if (n % 4n === 1n) {
    // exactly one of i, j is even
    i = 0n;
    ii = 0n;
    jj = j * j;
    while (ii <= jj) {
      const nn = n - ii;
      while (jj > nn) {
        j -= 1n;
        jj = j * j;
      }
      if (jj === nn) {
        return [i << fac, j << fac];
      }
      i += 1n;
      ii = i * i;
    }
  } else {
    // n = 2 mod 4: both i and j are odd
    i = 1n;
    ii = 1n;
    j += 1n - (j % 2n);
    jj = j * j;
    while (ii <= jj) {
      const nn = n - ii;
      while (jj > nn) {
        j -= 2n;
        jj = j * j;
      }
      if (jj === nn) {
        return [i << fac, j << fac];
      }
      i += 2n;
      ii = i * i;
    }
  }

  return null;
}

/**
 * Direct (factorization-free) search for `n = i^2 + j^2 + k^2`, `i <= j <= k`.
 *
 * Port of `three_squares_c` (`sage/rings/sum_of_squares.pyx:93`).
 *
 * @see Reference: sage/rings/sum_of_squares.pyx:three_squares_c
 */
function three_squares_c(n: bigint): [bigint, bigint, bigint] | null {
  if (n === 0n) {
    return [0n, 0n, 0n];
  }

  let fac = 0n;
  while (n % 4n === 0n) {
    n >>= 2n;
    fac++;
  }

  // Legendre: n is a sum of three squares iff it is not 4^a (8b + 7).
  if (n % 8n === 7n) {
    return null;
  }

  let i = isqrt(n);
  let res = two_squares_c(n - i * i);
  while (res === null) {
    i -= 1n;
    res = two_squares_c(n - i * i);
  }

  return [res[0] << fac, res[1] << fac, i << fac];
}

/**
 * Direct (factorization-free) search for `n = i^2 + j^2 + k^2 + l^2`.
 *
 * Port of `four_squares_pyx` (`sage/rings/sum_of_squares.pyx:274`).
 *
 * @see Reference: sage/rings/sum_of_squares.pyx:four_squares_pyx
 */
export function four_squares_pyx(input: IntegerLike): [bigint, bigint, bigint, bigint] {
  let n = uint32(input);
  if (n === 0n) {
    return [0n, 0n, 0n, 0n];
  }

  let fac = 0n;
  while (n % 4n === 0n) {
    n >>= 2n;
    fac++;
  }

  // Pick the largest square we can for j.
  let j = isqrt(n);
  let res = three_squares_c(n - j * j);
  while (res === null) {
    j -= 1n;
    res = three_squares_c(n - j * j);
  }

  return [res[0] << fac, res[1] << fac, res[2] << fac, j << fac];
}

/** Return the first ordered two-square decomposition, or raise ValueError.
 * @see Reference: sage/rings/sum_of_squares.pyx:two_squares_pyx
 */
export function two_squares_pyx(input: IntegerLike): [bigint, bigint] {
  const n = uint32(input);
  const result = two_squares_c(n);
  if (result === null) throw new ValueError(`${n} is not a sum of 2 squares`);
  return result;
}

/** Whether the uint32 input is a sum of two squares.
 * @see Reference: sage/rings/sum_of_squares.pyx:is_sum_of_two_squares_pyx
 */
export function is_sum_of_two_squares_pyx(input: IntegerLike): boolean {
  return two_squares_c(uint32(input)) !== null;
}

/** Return the original ordered three-square decomposition, or raise ValueError.
 * @see Reference: sage/rings/sum_of_squares.pyx:three_squares_pyx
 */
export function three_squares_pyx(input: IntegerLike): [bigint, bigint, bigint] {
  const n = uint32(input);
  const result = three_squares_c(n);
  if (result === null) throw new ValueError(`${n} is not a sum of 3 squares`);
  return result;
}
