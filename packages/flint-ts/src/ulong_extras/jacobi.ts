/** Binary word Jacobi algorithm from FLINT ulong_extras/jacobi.c. */
export function _n_jacobi_unsigned(x: bigint, y: bigint, r: number): number {
  if (x < 0n || x >= 1n << 64n || y <= 0n || y >= 1n << 64n || !(y & 1n))
    throw new RangeError('Jacobi arguments must be unsigned words with a positive odd denominator');
  r ^= 2;
  while (y > 1n) {
    if (x === 0n) return 0;
    let e = 0;
    while (!(x & 1n)) { x >>= 1n; e++; }
    r ^= Number((y ^ (y >> 1n)) & BigInt(2 * e));
    // The source's two-limb subtraction computes (abs(x-y), min(x,y)).
    // Only bit one of its sign accumulator affects the returned symbol.
    if (x < y) {
      r ^= Number(x & y & 2n);
      [x, y] = [y - x, x];
    } else x -= y;
  }
  return (r & 2) - 1;
}
export function n_jacobi_unsigned(x: bigint, y: bigint): number {
  return _n_jacobi_unsigned(x, y, 0);
}
export function n_jacobi(x: bigint, y: bigint): number {
  if (x < -(1n << 63n) || x >= 1n << 63n)
    throw new RangeError('Jacobi numerator must fit a signed word');
  return _n_jacobi_unsigned(x < 0n ? -x : x, y, x < 0n ? Number(y & 2n) : 0);
}
