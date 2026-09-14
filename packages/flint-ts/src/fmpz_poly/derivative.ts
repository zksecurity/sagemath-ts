/** Native fmpz_poly_derivative for dense constant-first coefficient arrays. */
export function _fmpz_poly_derivative(a: readonly bigint[]): bigint[] {
  let len = a.length;
  while (len && a[len - 1] === 0n) len--;
  return Array.from({ length: Math.max(0, len - 1) }, (_, i) => a[i + 1]! * BigInt(i + 1));
}
