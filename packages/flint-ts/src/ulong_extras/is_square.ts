/** FLINT ulong_extras/is_square.c residue filters and exact square check.
 * @see Deviation: FLINT word square roots and arithmetic
 */
export function n_is_square(x: bigint): boolean {
  if (x < 0n || x >= 1n << 64n) throw new RangeError('argument must fit an unsigned word');
  if (!((0x202021202030213n >> (x % 64n)) & 1n)) return false;
  if (!((0x402483012450293n >> (x % 63n)) & 1n)) return false;
  if (!((0x1218a019866014613n >> (x % 65n)) & 1n)) return false;
  if (x <= 1n) return true;
  // Keep the candidate exact instead of converting the native word to binary64.
  let root = 1n << BigInt(Math.ceil(x.toString(2).length / 2));
  let next = (root + x / root) >> 1n;
  while (next < root) {
    root = next;
    next = (root + x / root) >> 1n;
  }
  return root * root === x;
}
