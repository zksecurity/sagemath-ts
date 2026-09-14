/** nmod_poly/randtest.c: normalize the native test-coefficient vector.
 * @see Deviation: Native FLINT factor array kernels
 */
import type { flint_rand_t } from '../flint.js';
import { _nmod_vec_randtest } from '../nmod_vec/rand.js';
export function nmod_poly_randtest(state: flint_rand_t, len: number, p: bigint): bigint[] {
  const result = _nmod_vec_randtest(state, len, p);
  while (result.length && result[result.length - 1] === 0n) result.pop();
  return result;
}
