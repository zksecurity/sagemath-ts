/** Native fmpq_poly storage accessors; input is canonical numerator/denominator data. */
export function fmpq_poly_get_numerator(a: readonly bigint[], _den: bigint): bigint[] {
  return a.slice();
}

export function fmpq_poly_get_denominator(_a: readonly bigint[], den: bigint): bigint {
  return den;
}
