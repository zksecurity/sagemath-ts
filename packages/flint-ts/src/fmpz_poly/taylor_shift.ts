import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
import { _fmpz_poly_mul } from './mul.js';
/** FLINT's single-worker Horner/divide-and-conquer Taylor shift.
 * @see Deviation: Polynomial Evaluation and Composition
 */
export function _fmpz_poly_taylor_shift(a: readonly bigint[], c: bigint): bigint[] {
  const horner = (a: readonly bigint[], c: bigint): bigint[] => {
    const out = a.slice(),
      n = out.length;
    if (n <= 1 || c === 0n) return out;
    if (c !== 1n && !(n <= 4 || (c !== -1n && n <= 10))) {
      if (c === -1n)
        return horner(
          out.map((v, i) => (i % 2 ? -v : v)),
          1n
        ).map((v, i) => (i % 2 ? -v : v));
      let power = 1n;
      for (let i = 1; i < n; i++) {
        power *= c;
        out[i] = out[i]! * power;
      }
      const shifted = horner(out, 1n);
      for (let i = n - 1; i >= 1; i--) {
        shifted[i] = shifted[i]! / power;
        power /= c;
      }
      return shifted;
    }
    for (let i = n - 2; i >= 0; i--)
      for (let j = i; j < n - 1; j++) out[j] = out[j]! + c * out[j + 1]!;
    return out;
  };
  const shift = (a: readonly bigint[]): bigint[] => {
    const n = a.length;
    if (n < 50 || c === 0n) return horner(a, c);
    const cutoff = Math.min(1000, 100 + 10 * Math.floor(Math.sqrt(Math.max(k.maxBits(a) - 64, 0))));
    if (n < cutoff) return horner(a, c);
    const mid = Math.floor(n / 2),
      low = shift(a.slice(0, mid)),
      high = shift(a.slice(mid)),
      binomial = [1n];
    for (let i = 1; i <= mid; i++)
      binomial.push(
        i > mid - i ? binomial[mid - i]! : (binomial[i - 1]! * BigInt(mid + 1 - i)) / BigInt(i)
      );
    let power = c;
    for (let i = mid - 1; i >= 0; i--) {
      binomial[i] = binomial[i]! * power;
      power *= c;
    }
    const product = _fmpz_poly_mul(binomial, high);
    return Array.from({ length: n }, (_, i) => (low[i] ?? 0n) + (product[i] ?? 0n));
  };
  const A = k.normalized(a);
  return k.normalized(A.length < 64 ? horner(A, c) : shift(A));
}
