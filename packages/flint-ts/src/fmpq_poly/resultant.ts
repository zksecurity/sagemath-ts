import { _fmpz_poly_resultant_kernels as k, _fmpz_poly_gcd } from '../fmpz_poly/gcd.js';
import { _fmpz_poly_resultant } from '../fmpz_poly/resultant.js';
import { _nmod_poly_gcd } from '../nmod_poly/gcd.js';

/** FLINT rational resultant from integer numerator arrays and positive denominators. */
export function _fmpq_poly_resultant(
  a: readonly bigint[],
  denA: bigint,
  b: readonly bigint[],
  denB: bigint
): [bigint, bigint] {
  if (denA <= 0n || denB <= 0n) throw new RangeError('input denominators must be positive');
  let A = k.normalized(a),
    B = k.normalized(b),
    dA = denA,
    dB = denB;
  if (!A.length || !B.length) return [0n, 1n];
  // fmpq_poly_t stores canonical numerator/denominator pairs.
  const hA = k.content([...A, dA]),
    hB = k.content([...B, dB]);
  A = A.map((c) => c / hA);
  dA /= hA;
  B = B.map((c) => c / hB);
  dB /= hB;
  const swapped = A.length < B.length;
  if (swapped) {
    [A, B] = [B, A];
    [dA, dB] = [dB, dA];
  }
  const m = BigInt(A.length - 1),
    n = BigInt(B.length - 1);
  const canonical = (num: bigint, den: bigint): [bigint, bigint] => {
    const c = k.content([num, den]);
    return [swapped && m & n & 1n ? -num / c : num / c, den / c];
  };
  if (B.length === 1) return canonical(B[0]! ** m, dB ** m);
  const cA = k.content(A),
    cB = k.content(B);
  A = A.map((c) => c / cA);
  B = B.map((c) => c / cB);
  let p = 1152921504606846869n;
  while (A[A.length - 1]! % p === 0n || B[B.length - 1]! % p === 0n) p = nextPrime(p);
  const modG = _nmod_poly_gcd(
    A.map((c) => k.mod(c, p)),
    B.map((c) => k.mod(c, p)),
    p
  );
  if (modG.length > 1 && _fmpz_poly_gcd(A, B).length > 1) return [0n, 1n];
  return canonical(_fmpz_poly_resultant(A, B) * cA ** n * cB ** m, dA ** n * dB ** m);
}
// The initial 60-bit prime always uses n_nextprime's large-input mod-30 wheel.
function nextPrime(n: bigint): bigint {
  if (n >= 18446744073709551557n)
    throw new RangeError('Exception (n_nextprime). No larger single-limb prime exists.');
  const wheel = [
    1, 6, 5, 4, 3, 2, 1, 4, 3, 2, 1, 2, 1, 4, 3, 2, 1, 2, 1, 4, 3, 2, 1, 6, 5, 4, 3, 2, 1, 2,
  ];
  do {
    n += BigInt(wheel[Number(n % 30n)]!);
  } while (!k.wordIsPrime(n));
  return n;
}
