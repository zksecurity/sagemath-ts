/** PARI FpX_normalize / Flx_normalize and their forced-monic scalar products.
 * @see Deviation: PARI polynomial normalization adapters
 */
import { inverseCoefficient, residue } from './_polynomial_division.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { PariError } from './errors.js';
export function polynomialNormalize(f: bigint[], p: bigint, word: boolean): bigint[] {
  if (word) {
    if (p < 1n || p >= 1n << 64n) throw new RangeError('modulus must be a positive word integer');
    if (f.some((c) => c < 0n || c >= p))
      throw new RangeError('word polynomial coefficients must be reduced');
  }
  f = trimPolynomial(f.slice());
  const lead = f.at(-1) ?? 0n;
  if ((!word && !f.length) || lead === 1n) return f;
  if (p === 0n)
    throw new PariError(`impossible inverse in Fp_inv: Mod(${lead < 0n ? -lead : lead}, 0).`);
  const m = p < 0n ? -p : p;
  let inverse: bigint;
  try {
    inverse = inverseCoefficient(lead, m, word);
  } catch (e) {
    if (p < 0n && e instanceof PariError)
      throw new PariError(e.message.replace(`, ${m}).`, `, ${p}).`));
    throw e;
  }
  if (!f.length) return [];
  return [...f.slice(0, -1).map((c) => residue(c * inverse, m)), 1n];
}
