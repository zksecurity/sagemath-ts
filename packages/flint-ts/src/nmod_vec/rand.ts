/** Dense-array port of nmod_vec/rand.c.
 * @see Deviation: Native FLINT factor array kernels
 */
import type { flint_rand_t } from '../flint.js';
import { n_randint, n_randtest } from '../ulong_extras/randomisation.js';
export function _nmod_vec_rand(state: flint_rand_t, len: number, p: bigint): bigint[] {
  if (!Number.isSafeInteger(len) || len < 0) throw new RangeError('length must be nonnegative');
  return Array.from({ length: len }, () => n_randint(state, p));
}
export function _nmod_vec_randtest(state: flint_rand_t, len: number, p: bigint): bigint[] {
  if (!Number.isSafeInteger(len) || len < 0) throw new RangeError('length must be nonnegative');
  if (n_randint(state, 2n)) return Array.from({ length: len }, () => n_randtest(state) % p);
  const sparseness = 1n + n_randint(state, BigInt(Math.max(2, len)));
  return Array.from({ length: len }, () =>
    n_randint(state, sparseness) ? 0n : n_randtest(state) % p
  );
}
