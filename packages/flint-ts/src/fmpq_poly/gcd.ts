import { _fmpz_poly_gcd } from '../fmpz_poly/gcd.js';

/** Monic rational GCD from dense integer numerator arrays; input denominators cancel.
 * Returns [numerator coefficients, positive common denominator], without mutation.
 * @see Deviation: Polynomial GCD Backend Coercion and Identity
 * @see Reference: fmpq_poly/gcd.c and fmpq_poly/make_monic.c
 */
export function _fmpq_poly_gcd(a: readonly bigint[], b: readonly bigint[]): [bigint[], bigint] {
  const primitive = (v: readonly bigint[]): bigint[] => {
    const r = [...v];
    while (r.length && r[r.length - 1] === 0n) r.pop();
    let content = 0n;
    for (const coefficient of r) {
      let c = coefficient < 0n ? -coefficient : coefficient;
      while (c) [content, c] = [c, content % c];
    }
    if (!content) return [];
    if (r[r.length - 1]! < 0n) content = -content;
    return r.map((c) => c / content);
  };
  const A = primitive(a),
    B = primitive(b);
  const G = !B.length ? A : !A.length ? B : _fmpz_poly_gcd(A, B);
  return [G, G.length ? G[G.length - 1]! : 1n];
}
