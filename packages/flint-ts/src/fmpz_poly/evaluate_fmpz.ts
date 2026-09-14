import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
/** FLINT evaluate_fmpz.c: Horner up to length 50, otherwise binary splitting.
 * @see Deviation: Polynomial Evaluation and Composition
 */
export function _fmpz_poly_evaluate_fmpz(a: readonly bigint[], x: bigint): bigint {
  const A = k.normalized(a),
    n = A.length;
  if (n <= 50) {
    if (!n) return 0n;
    if (n === 1 || x === 0n) return A[0]!;
    let out = A[n - 1]!;
    for (let i = n - 2; i >= 0; i--) out = A[i]! + out * x;
    return out;
  }
  const h = (n - 1).toString(2).length,
    powers = [x],
    blocks: bigint[] = [];
  for (let i = 1; i < h; i++) powers.push(powers[i - 1]! * powers[i - 1]!);
  const zeros = (n: number) => {
    let c = 0;
    while (n % 2 === 0) {
      c++;
      n /= 2;
    }
    return c;
  };
  let level = 1;
  for (let i = 0; i < n - 1; i += 2) {
    let value = A[i]! + x * A[i + 1]!;
    const c = zeros(i + 2);
    for (level = 1; level < c; level++) value = blocks[level]! + powers[level]! * value;
    blocks[level] = value;
  }
  if (n % 2) {
    let value = A[n - 1]!;
    const c = zeros(n + 1);
    for (level = 1; level < c; level++) value = blocks[level]! + powers[level]! * value;
    blocks[level] = value;
  }
  let value = blocks[level]!;
  for (; level < h; level++)
    if (Math.floor((n - 1) / 2 ** level) % 2) value = blocks[level]! + powers[level]! * value;
  return value;
}
