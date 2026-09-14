/** Native bibli2.c vecbinomial; n is a nonnegative integer vector size.
 * @see Deviation: PARI factorization scalar and bound adapters
 */
export function vecbinomial(n: number): bigint[] {
  if (!Number.isSafeInteger(n) || n < 0)
    throw new RangeError('vecbinomial requires a nonnegative integer size');
  if (!n) return [1n];
  const C: bigint[] = new Array(n + 1);
  C[0] = 1n;
  C[1] = BigInt(n);
  const d = Math.floor((n + 1) / 2);
  let k = 2;
  for (; k <= d; k++) C[k] = (BigInt(n - k + 1) * C[k - 1]!) / BigInt(k);
  for (; k <= n; k++) C[k] = C[n - k]!;
  return C;
}
