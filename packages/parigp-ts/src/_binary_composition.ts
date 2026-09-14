/** PARI F2x.c / RgX.c Brent–Kung blocks and native giant-step order.
 * @see Deviation: PARI packed binary-polynomial kernels
 */
import { F2x_degree, F2x_rem, F2xq_mul, F2xq_powers } from './F2x.js';
import { PariError } from './errors.js';
export function binaryComposition(
  Q: bigint,
  x: bigint,
  V: bigint[],
  T: bigint,
  direct: boolean
): bigint {
  const degree = F2x_degree(Q);
  if (degree < 0) return 0n;
  if (direct) V = F2xq_powers(x, Math.floor(Math.sqrt(degree)), T);
  const l = V.length;
  const cmul = (a: number, x: bigint) => ((Q >> BigInt(a)) & 1n ? x : 0n);
  const block = (a: number, n: number) => {
    let z = cmul(a, 1n);
    for (let i = 1; i <= n; i++) z ^= cmul(a + i, V[i]!);
    return F2x_rem(z, T);
  };
  if (degree < l) return block(0, degree);
  if (l < 2) throw new PariError('domain error in gen_RgX_bkeval_powers: #powers < 2');
  let d = degree - l,
    z = block(d + 1, l - 1);
  while (d >= l - 1) {
    d -= l - 1;
    z = block(d + 1, l - 2) ^ F2xq_mul(z, V[l - 1]!, T);
  }
  z = block(0, d) ^ F2xq_mul(z, V[d + 1]!, T);
  return F2x_rem(z, T);
}
