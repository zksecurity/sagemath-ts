/** Native integer/real branches of PARI RgV.c.
 * @see Deviation: PARI integer and real matrix products and rescaling
 */
import { type MpReal, addir, addrr, mulir, mulrr, sqrr } from './qfb.js';
import { ZM_mul } from './ZV.js';
export type RgScalar = bigint | MpReal;
const add = (x: RgScalar, y: RgScalar): RgScalar =>
  typeof x === 'bigint'
    ? typeof y === 'bigint'
      ? x + y
      : addir(x, y)
    : typeof y === 'bigint'
      ? addir(y, x)
      : addrr(x, y);
const mul = (x: RgScalar, y: RgScalar): RgScalar =>
  x === 0n || y === 0n
    ? 0n
    : typeof x === 'bigint'
      ? typeof y === 'bigint'
        ? x * y
        : mulir(x, y)
      : typeof y === 'bigint'
        ? mulir(y, x)
        : mulrr(x, y);
const sqr = (x: RgScalar): RgScalar => (typeof x === 'bigint' ? x * x : sqrr(x));

/** Native integer/real dot square with sequential summation. */
export function RgV_dotsquare(x: RgScalar[]): RgScalar {
  if (!x.length) return 0n;
  let z = sqr(x[0]!);
  for (let i = 1; i < x.length; i++) z = add(z, sqr(x[i]!));
  return z;
}
/** Native integer/real dot product and pointer-identity square dispatch. */
export function RgV_dotproduct(x: RgScalar[], y: RgScalar[]): RgScalar {
  if (x.length !== y.length) throw new RangeError('RgV_dotproduct requires equal lengths');
  if (x === y) return RgV_dotsquare(x);
  if (!x.length) return 0n;
  let z = mul(x[0]!, y[0]!);
  for (let i = 1; i < x.length; i++) z = add(z, mul(x[i]!, y[i]!));
  return z;
}
/** Native Gram matrix, preserving generic integer/real scalar arithmetic.
 * @see Deviation: PARI integer and real matrix products and rescaling
 */
export function gram_matrix(x: RgScalar[][]): RgScalar[][] {
  if (x.some((c) => c.length !== x[0]!.length))
    throw new RangeError('gram_matrix requires rectangular columns');
  const M = x.map(() => Array<RgScalar>(x.length).fill(0n));
  for (let i = 0; i < x.length; i++) {
    for (let j = 0; j < i; j++) M[j]![i] = M[i]![j] = RgV_dotproduct(x[i]!, x[j]!);
    M[i]![i] = RgV_dotsquare(x[i]!);
  }
  return M;
}
/** Native generic matrix product; integers delegate to ZM_mul.
 * @see Deviation: PARI integer and real matrix products and rescaling
 */
export function RgM_mul(x: RgScalar[][], y: RgScalar[][]): RgScalar[][] {
  if (!y.length) return [];
  if (y.some((c) => c.length !== x.length) || x.some((c) => c.length !== x[0]!.length))
    throw new RangeError('RgM_mul requires compatible rectangular columns');
  if (!x.length) return y.map(() => []);
  if ([...x, ...y].every((c) => c.every((v) => typeof v === 'bigint')))
    return ZM_mul(
      [[], ...x.map((c) => [0n, ...(c as bigint[])])],
      [[], ...y.map((c) => [0n, ...(c as bigint[])])]
    )
      .slice(1)
      .map((c) => c.slice(1));
  return y.map((c) =>
    x[0]!.map((_, i) => {
      let z = mul(x[0]![i]!, c[0]!);
      for (let j = 1; j < x.length; j++) if (x[j]![i] !== 0n) z = add(z, mul(x[j]![i]!, c[j]!));
      return z;
    })
  );
}
