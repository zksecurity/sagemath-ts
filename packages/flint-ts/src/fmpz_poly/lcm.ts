import { _fmpz_poly_gcd } from './gcd.js';
import { _fmpz_poly_mul } from './mul.js';
import { _fmpz_poly_divrem } from './divrem.js';

/** Dense FLINT fmpz_poly/lcm.c; result has positive leading coefficient. */
export function _fmpz_poly_lcm(a: readonly bigint[], b: readonly bigint[]): bigint[] {
  if (!a.length || !b.length) return [];
  const product = _fmpz_poly_mul(a, b), g = _fmpz_poly_gcd(a, b);
  if (g.length === 1) {
    const divisor = product[product.length - 1]! < 0n ? -g[0]! : g[0]!;
    return product.map((c) => c / divisor);
  }
  const q = _fmpz_poly_divrem(product, g, true)![0];
  return q[q.length - 1]! > 0n ? q : q.map((c) => -c);
}
