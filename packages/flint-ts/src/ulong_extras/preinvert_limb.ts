/** Native-division branch of FLINT ulong_extras/preinvert_limb.c.
 * BigInt division replaces udiv_qrnnd on the normalized two-limb numerator.
 */
export function n_preinvert_limb(n: bigint): bigint {
  if (n <= 0n || n >= 1n << 64n) throw new RangeError('modulus must fit a positive unsigned word');
  const normalized = n << BigInt(64 - n.toString(2).length);
  return ((1n << 128n) - 1n) / normalized - (1n << 64n);
}
