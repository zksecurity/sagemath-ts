/** Dense-array port of fmpz_poly/pseudo_divrem_basecase.c.
 * Inputs have explicit nonzero leading coefficients; recursive dividends may
 * retain leading zero slots. Native GMP scalar operations map to BigInt.
 * @see Deviation: Number-field reduction and native array adapters
 */
export function _fmpz_poly_pseudo_divrem_basecase(
  a: readonly bigint[],
  b: readonly bigint[]
): [bigint[], bigint[], bigint] {
  if (!b.length || b[b.length - 1] === 0n) throw new RangeError('division by zero polynomial');
  if (a.length < b.length) return [[], a.slice(), 0n];
  const q = Array<bigint>(a.length - b.length + 1).fill(0n),
    r = a.slice();
  const lead = b[b.length - 1]!;
  let d = 0n;
  for (let ir = a.length - 1, iq = q.length - 1; iq >= 0; ir--, iq--) {
    const top = r[ir]!;
    if (top % lead === 0n) q[iq] = top / lead;
    else {
      for (let j = 0; j < q.length; j++) q[j] = q[j]! * lead;
      q[iq] = top;
      for (let j = 0; j < r.length; j++) r[j] = r[j]! * lead;
      d++;
    }
    for (let j = 0; j < b.length - 1; j++) r[iq + j] = r[iq + j]! - b[j]! * q[iq]!;
    r[ir] = 0n;
  }
  return [q, r.slice(0, b.length - 1), d];
}
