/** NTL lip.cpp integer primitive adapter.
 * @see Deviation: NTL exact integer lattice adapters
 */
export function _ntl_gexteucl(a: bigint, b: bigint): [bigint, bigint, bigint] {
  if (b === 0n) return [a < 0n ? -1n : 1n, 0n, a < 0n ? -a : a];
  if (a === 0n) return [0n, b < 0n ? -1n : 1n, b < 0n ? -b : b];
  let A = a < 0n ? -a : a,
    B = b < 0n ? -b : b,
    rev = false;
  if (A < B) {
    [A, B] = [B, A];
    [a, b] = [b, a];
    rev = true;
  }
  let r = A,
    s = B,
    u = 1n,
    v = 0n;
  while (s) {
    const q = r / s;
    [r, s] = [s, r - q * s];
    [u, v] = [v, u - q * v];
  }
  let x = a < 0n ? -u : u,
    y = (r - a * x) / b;
  if (rev) [x, y] = [y, x];
  return [x, y, r];
}
