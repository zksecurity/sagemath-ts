import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
/** FLINT evaluate_fmpq.c, retaining the native internal unreduced denominator.
 * @see Deviation: Polynomial Evaluation and Composition
 */
export function _fmpz_poly_evaluate_fmpq(
  a: readonly bigint[],
  numerator: bigint,
  denominator: bigint
): [bigint, bigint] {
  if (denominator <= 0n) throw new RangeError('input denominator must be positive');
  const A = k.normalized(a),
    n = A.length;
  if (n < 40 || denominator.toString(2).length > 0.003 * n * n) {
    if (!n) return [0n, 1n];
    let num = A[n - 1]!,
      den = 1n;
    for (let i = n - 2; i >= 0; i--) {
      num *= numerator;
      den *= denominator;
      num += den * A[i]!;
      if (num === 0n) den = 1n;
    }
    return [num, den];
  }
  type Pair = [bigint, bigint];
  const h = (n - 1).toString(2).length,
    powers: Pair[] = [[numerator, denominator]],
    blocks: Pair[] = [];
  for (let i = 1; i < h; i++) {
    const [num, den] = powers[i - 1]!;
    powers.push([num * num, den * den]);
  }
  const combine = (low: Pair, power: Pair, high: Pair): Pair => [
    low[0] * power[1] * high[1] + power[0] * high[0] * low[1],
    low[1] * power[1] * high[1],
  ];
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
    let value: Pair = [numerator * A[i + 1]! + denominator * A[i]!, denominator];
    const c = zeros(i + 2);
    for (level = 1; level < c; level++) value = combine(blocks[level]!, powers[level]!, value);
    blocks[level] = value;
  }
  if (n % 2) {
    let value: Pair = [A[n - 1]!, 1n];
    const c = zeros(n + 1);
    for (level = 1; level < c; level++) value = combine(blocks[level]!, powers[level]!, value);
    blocks[level] = value;
  }
  let value = blocks[level]!;
  for (; level < h; level++)
    if (Math.floor((n - 1) / 2 ** level) % 2)
      value = combine(blocks[level]!, powers[level]!, value);
  return value;
}
