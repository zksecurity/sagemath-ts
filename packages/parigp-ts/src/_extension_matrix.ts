/** alglin1.c / Flx.c packed extension-matrix backend. */
import { ZM_mul } from './ZV.js';
import { polynomialQuotient } from './_polynomial_quotient.js';
import { extensionField } from './_extension_field.js';
import { trimPolynomial, packUnsigned, unpackUnsigned } from './_polynomial_packing.js';
/** Native FlxqM_mul schedule on columns without sentinel entries.
 * @see Deviation: PARI extension composition adapters
 */
export function wordExtensionMatrix(
  A: bigint[][][],
  B: bigint[][][],
  T: bigint[],
  p: bigint,
  providedInverse?: bigint[]
): bigint[][][] {
  T = trimPolynomial(T);
  A = A.map((c) => c.map(trimPolynomial));
  B = B.map((c) => c.map(trimPolynomial));
  let inverse = providedInverse;
  const n = A.length;
  if (!n) return [];
  if (n === 1) {
    inverse ??= polynomialQuotient(T, p, true).inverse;
    const f = extensionField(1, T, p, inverse);
    return B.map((col) => A[0]!.map((a) => f.red(f.mul(a, col[0]!, false)) as bigint[]));
  }
  const z = (p - 1n) ** 2n * BigInt(T.length - 1) * BigInt(n);
  let bits = z.toString(2).length;
  if (bits <= 32) {
    if (Math.ceil(((T.length - 1) * bits) / 64) === Math.floor(T.length / 2)) bits = 32;
  } else {
    const words = Math.ceil(bits / 64);
    if (Math.ceil(((T.length - 1) * bits) / 64) === (T.length - 1) * words) bits = words * 64;
  }
  const pack = (M: bigint[][][]) => [
    [],
    ...M.map((c) => [0n, ...c.map((v) => packUnsigned(v, bits))]),
  ];
  const product = ZM_mul(pack(A), pack(B));
  const f = extensionField(1, T, p, inverse);
  return product
    .slice(1)
    .map((c) =>
      c
        .slice(1)
        .map(
          (x) =>
            f.red(
              trimPolynomial(
                unpackUnsigned(x, bits, x === 0n ? 0 : Math.ceil(x.toString(2).length / bits)).map(
                  (v) => v % p
                )
              )
            ) as bigint[]
        )
    );
}
