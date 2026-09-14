/** Dense-array port of nmod_poly/add.c.
 * @see Deviation: Native Modular Polynomial Products and Fraction Fields
 */
export function _nmod_poly_add(a: readonly bigint[], b: readonly bigint[], p: bigint): bigint[] {
  if (p < 1n || p >= 1n << 64n) throw new RangeError('modulus must fit a positive unsigned word');
  const out = Array.from({ length: Math.max(a.length, b.length) }, (_, i) => {
    const r = ((a[i] ?? 0n) + (b[i] ?? 0n)) % p;
    return r < 0n ? r + p : r;
  });
  while (out.length && out[out.length - 1] === 0n) out.pop();
  return out;
}
