/** nmod_poly/mod_matrix_rows_evaluate.c: Horner through ten rows, then rectangular splitting. */
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
import { _nmod_poly_mul } from './mul.js';
import { _nmod_poly_preinv_remainder } from './powmod_binexp_preinv.js';
export function _nmod_poly_mod_matrix_rows_evaluate(
  rows: readonly (readonly bigint[])[],
  h: readonly bigint[],
  f: readonly bigint[],
  finv: readonly bigint[],
  p: bigint
): bigint[] {
  if (!rows.length) throw new RangeError('matrix must have at least one row');
  const add = (a: readonly bigint[], b: readonly bigint[]) =>
    k.normalized(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) + (b[i] ?? 0n)),
      p
    );
  const reduce = (a: bigint[]) => _nmod_poly_preinv_remainder(a, f, finv, p);
  const mul = (a: readonly bigint[], b: readonly bigint[]) => reduce(_nmod_poly_mul(a, b, p));
  if (rows.length <= 10) {
    let result = k.normalized(rows[rows.length - 1]!, p);
    for (let i = rows.length - 2; i >= 0; i--) result = add(mul(result, h), rows[i]!);
    return result;
  }
  const m = Math.min(Math.floor(Math.sqrt(rows.length)) + 1, 30),
    r = Math.ceil(rows.length / m);
  const powers: bigint[][] = [[1n], k.normalized(h, p)];
  for (let i = 2; i <= m; i++)
    powers.push(mul(powers[Math.floor(i / 2)]!, powers[Math.ceil(i / 2)]!));
  const evaluateBlock = (i: number) => {
    let sum = k.normalized(rows[i * m]!, p);
    for (let j = 1; j < m && i * m + j < rows.length; j++)
      sum = add(sum, _nmod_poly_mul(powers[j]!, rows[i * m + j]!, p));
    return reduce(sum);
  };
  let result = evaluateBlock(r - 1);
  for (let i = r - 2; i >= 0; i--) result = add(mul(result, powers[m]!), evaluateBlock(i));
  return result;
}
