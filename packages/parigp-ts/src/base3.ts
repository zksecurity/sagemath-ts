/** PARI base3.c prime valuation using the anti-uniformizer multiplication table. */
import type { NfPrimeIdeal } from './base2.js';
import { gcd } from './ff.js';
import { Z_pvalrem } from './gen2.js';

/** Native ZC_nfvalrem without a remainder output, on a nonzero integral vector.
 * @see Deviation: Number-field ideal valuation adapters
 */
export function ZC_nfval(input: bigint[], pr: NfPrimeIdeal): bigint {
  let content = input.reduce((g, c) => gcd(g, c), 0n);
  if (!content) throw new RangeError('ZC_nfval requires a nonzero vector');
  if (pr.tau === 1n) return BigInt(Z_pvalrem(content, pr.p)[0]);
  const tau = pr.tau;
  let x = input.slice(),
    v = 0n;
  for (; ; v++) {
    const y = x.map((_, i) => x.reduce((s, c, j) => s + c * tau[j]![i]!, 0n));
    if (y.some((c) => c % pr.p !== 0n)) return v;
    x = y.map((c) => c / pr.p);
    if ((v & 15n) === 15n) {
      content = x.reduce((g, c) => gcd(g, c), 0n);
      const k = BigInt(Z_pvalrem(content, pr.p)[0]);
      x = x.map((c) => c / pr.p ** k);
      v += pr.e * k;
    }
  }
}

import { ZM_gauss } from './alglin1.js';
/** Native zkmultable_inv: solve multiplication-by-x against the first unit column.
 * Input is a nonsingular integral multiplication table, stored by columns.
 * @see Deviation: Number-field full prime decomposition adapters
 */
export function zkmultable_inv(mx: bigint[][]): [bigint[], bigint] {
  const n = mx.length;
  const result = ZM_gauss(
    [[], ...mx.map((c) => [0n, ...c])],
    [[], [0n, ...Array.from({ length: n }, (_, i) => (i === 0 ? 1n : 0n))]]
  );
  if (!result) throw new RangeError('zkmultable_inv requires a nonsingular multiplication table');
  const vector = result[0][1]!.slice(1),
    d = result[1];
  const content = vector.reduce((g, c) => gcd(g, c), d);
  return [vector.map((c) => c / content), d / content];
}
/** Native zkmultable_capZ: denominator of the inverse field element.
 * @see Deviation: Number-field full prime decomposition adapters
 */
export function zkmultable_capZ(mx: bigint[][]): bigint {
  return zkmultable_inv(mx)[1];
}
