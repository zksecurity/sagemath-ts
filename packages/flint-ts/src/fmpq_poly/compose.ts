import { _fmpz_poly_resultant_kernels as k } from '../fmpz_poly/gcd.js';
import { _fmpz_poly_compose } from '../fmpz_poly/compose.js';
/** FLINT fmpq_poly/compose.c: rescale denominators, compose over ZZ, canonicalize.
 * @see Deviation: Polynomial Evaluation and Composition
 */
export function _fmpq_poly_compose(
  a: readonly bigint[],
  da: bigint,
  b: readonly bigint[],
  db: bigint
): [bigint[], bigint] {
  if (da <= 0n || db <= 0n) throw new RangeError('input denominators must be positive');
  let A = k.normalized(a),
    den = da;
  const B = k.normalized(b);
  const canonical = (a: bigint[], den: bigint): [bigint[], bigint] => {
    const g = k.content([...a, den]);
    return [a.map((c) => c / g), den / g];
  };
  if (!A.length) return [[], 1n];
  if (A.length === 1 || !B.length) return canonical(k.normalized([A[0]!]), den);
  if (db !== 1n) {
    let power = 1n;
    for (let i = A.length - 2; i >= 0; i--) {
      power *= db;
      A[i] = A[i]! * power;
    }
    den *= power;
    [A, den] = canonical(A, den);
  }
  return canonical(_fmpz_poly_compose(A, B), den);
}
