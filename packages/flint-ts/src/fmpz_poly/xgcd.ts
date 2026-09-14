import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
import { _fmpz_poly_resultant } from './resultant.js';
import { _nmod_poly_xgcd } from '../nmod_poly/xgcd.js';
import { _nmod_poly_xgcd_kernels as modK } from '../nmod_poly/gcd.js';

/** FLINT resultant-scaled integer Bezout coefficients, ordered by input length.
 * @see Deviation: Polynomial Integer and Rational Extended GCD
 */
export function _fmpz_poly_xgcd(
  a: readonly bigint[],
  b: readonly bigint[]
): [bigint, bigint[], bigint[]] {
  let A = k.normalized(a),
    B = k.normalized(b);
  const swapped = A.length < B.length;
  if (swapped) [A, B] = [B, A];
  if (!A.length || !B.length) return [0n, [], []];
  // Native constant/constant input is undefined unless the selected divisor is a unit.
  // In particular its modular reconstruction need not terminate for nonunit constants.
  if (A.length === 1 && B[0] !== 1n && B[0] !== -1n)
    throw new RangeError('constant inputs require a unit right operand');
  const resultant = _fmpz_poly_resultant(A, B);
  if (!resultant) return [0n, [], []];
  let S = Array<bigint>(B.length).fill(0n),
    T = Array<bigint>(A.length).fill(0n);
  let product = 1n,
    first = true,
    stabilised = false,
    sBits = 0,
    tBits = 0;
  for (const p of k.modularPrimes()) {
    const R = k.mod(resultant, p);
    if (!R || A[A.length - 1]! % p === 0n || B[B.length - 1]! % p === 0n) continue;
    const Ap = A.map((c) => k.mod(c, p)),
      Bp = B.map((c) => k.mod(c, p));
    if (stabilised) {
      const check = modK.add(
        modK.mul(
          Ap,
          S.map((c) => k.mod(c, p)),
          p
        ),
        modK.mul(
          Bp,
          T.map((c) => k.mod(c, p)),
          p
        ),
        p
      );
      if (check.length === 1 && check[0] === R) product *= p;
      else stabilised = false;
    }
    if (!stabilised) {
      const [, s, t] = _nmod_poly_xgcd(Ap, Bp, p);
      const sp = Array.from({ length: B.length }, (_, i) => ((s[i] ?? 0n) * R) % p);
      const tp = Array.from({ length: A.length }, (_, i) => ((t[i] ?? 0n) * R) % p);
      if (first) {
        S = sp.map((c) => (c > p / 2n ? c - p : c));
        T = tp.map((c) => (c > p / 2n ? c - p : c));
        product = p;
        stabilised = true;
        first = false;
      } else {
        const inv = k.inverse(product % p, p),
          modulus = product * p;
        const crt = (old: bigint[], values: bigint[]) =>
          old.map((c, i) => {
            const r = c + product * k.mod((values[i]! - c) * inv, p);
            return r > modulus / 2n ? r - modulus : r;
          });
        S = crt(S, sp);
        T = crt(T, tp);
        product = modulus;
        const sb = k.maxBits(S),
          tb = k.maxBits(T);
        stabilised = sb === sBits && tb === tBits;
        [sBits, tBits] = [sb, tb];
      }
    }
    if (stabilised) {
      const lenBits = k.bits(BigInt(B.length));
      const bound =
        4 +
        Math.max(
          k.bits(resultant),
          lenBits + k.maxBits(A) + k.maxBits(S),
          lenBits + k.maxBits(B) + k.maxBits(T)
        );
      if (k.bits(product) > bound) {
        const s = k.normalized(S),
          t = k.normalized(T);
        return swapped ? [resultant, t, s] : [resultant, s, t];
      }
    }
  }
  throw new Error('native word-prime stream exhausted');
}
