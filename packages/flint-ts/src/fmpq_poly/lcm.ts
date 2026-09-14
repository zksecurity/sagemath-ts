import { _fmpz_poly_lcm } from '../fmpz_poly/lcm.js';
import { _fmpz_poly_resultant_kernels as k } from '../fmpz_poly/gcd.js';

/** Dense FLINT fmpq_poly/lcm.c; input denominators cancel in the monic result. */
export function _fmpq_poly_lcm(a: readonly bigint[], b: readonly bigint[]): [bigint[], bigint] {
  if (!a.length || !b.length) return [[], 1n];
  if (a.length === 1 && b.length === 1) return [[1n], 1n];
  const ca = k.content(a), cb = k.content(b);
  const r = _fmpz_poly_lcm(a.map((c) => c / ca), b.map((c) => c / cb));
  return [r, r[r.length - 1]!];
}
