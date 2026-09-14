import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
import { _fmpz_poly_mul } from './mul.js';
import { _fmpz_poly_evaluate_fmpz } from './evaluate_fmpz.js';
import { _fmpz_poly_taylor_shift } from './taylor_shift.js';
/** FLINT fmpz_poly/compose.c and its shared-power divide-and-conquer tree.
 * @see Deviation: Polynomial Evaluation and Composition
 */
export function _fmpz_poly_compose(a: readonly bigint[], b: readonly bigint[]): bigint[] {
  const A = k.normalized(a),
    B = k.normalized(b);
  if (!A.length) return [];
  if (A.length === 1 || !B.length) return k.normalized([A[0]!]);
  if (B.length === 1) return k.normalized([_fmpz_poly_evaluate_fmpz(A, B[0]!)]);
  const add = (a: readonly bigint[], b: readonly bigint[]) =>
    k.normalized(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) + (b[i] ?? 0n))
    );
  if (A.length <= 4) {
    let out = [A[A.length - 1]!];
    for (let i = A.length - 2; i >= 0; i--) out = add(_fmpz_poly_mul(out, B), [A[i]!]);
    return out;
  }
  if (B.length === 2) {
    const out = _fmpz_poly_taylor_shift(A, B[0]!);
    let power = 1n;
    for (let i = 1; i < out.length; i++) {
      power *= B[1]!;
      out[i] = out[i]! * power;
    }
    return k.normalized(out);
  }
  let blocks: bigint[][] = [];
  for (let i = 0; i < A.length; i += 2)
    blocks.push(
      add(
        B.map((c) => c * (A[i + 1] ?? 0n)),
        [A[i]!]
      )
    );
  let power = _fmpz_poly_mul(B, B);
  while (blocks.length > 2) {
    const next: bigint[][] = [];
    for (let i = 0; i < blocks.length; i += 2)
      next.push(
        i + 1 < blocks.length ? add(_fmpz_poly_mul(power, blocks[i + 1]!), blocks[i]!) : blocks[i]!
      );
    blocks = next;
    power = _fmpz_poly_mul(power, power);
  }
  return add(_fmpz_poly_mul(power, blocks[1]!), blocks[0]!);
}
