import { _nmod_poly_divrem_kernel } from './gcd.js';
/** nmod_poly/divrem.c: dense arrays, nonzero divisor with an invertible leading coefficient. */
export function _nmod_poly_divrem(
  a: readonly bigint[],
  b: readonly bigint[],
  n: bigint
): [bigint[], bigint[]] {
  if (n < 2n) throw new RangeError('modulus must be at least 2');
  const normalize = (x: readonly bigint[]) => {
    const r = x.map((c) => ((c % n) + n) % n);
    while (r.length && r[r.length - 1] === 0n) r.pop();
    return r;
  };
  return _nmod_poly_divrem_kernel(normalize(a), normalize(b), n);
}
