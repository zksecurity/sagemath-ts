/** FLINT ulong_extras/powmod2_ui_preinv.c binary powering.
 * @see Deviation: FLINT word square roots and arithmetic
 */
export function n_powmod2_ui_preinv(a: bigint, e: bigint, n: bigint, _ninv: bigint): bigint {
  if (a < 0n || a >= 1n << 64n || e < 0n || e >= 1n << 64n || n <= 0n || n >= 1n << 64n)
    throw new RangeError('modular power arguments must fit unsigned words with nonzero modulus');
  if (e === 0n) return n === 1n ? 0n : 1n;
  if (a === 0n) return 0n;
  a %= n;
  // BigInt residues do not need the source's word-normalization shift.
  while (!(e & 1n)) { a = a * a % n; e >>= 1n; }
  let x = a;
  while ((e >>= 1n) !== 0n) {
    a = a * a % n;
    if (e & 1n) x = x * a % n;
  }
  return x;
}
