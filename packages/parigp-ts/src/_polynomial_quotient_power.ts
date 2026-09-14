/** Native FpXQ/Flxq inverse and signed-power phases.
 * @see Deviation: PARI prime quotient cache and inverse errors
 */
import { FpX_red, FpX_Fp_mul, type FpX } from './ffinit.js';
import { _FpX_extgcd, FpX_rem, FpX_sqr } from './FpX.js';
import { _Flx_extgcd, Flx_sqr } from './Flx.js';
import { inverseCoefficient } from './_polynomial_division.js';
import { polynomialQuotient } from './_polynomial_quotient.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { displayExtensionError } from './_extension_display.js';
import { gen_pow_i } from './bb_group.js';
import { PariError } from './errors.js';

export function polynomialQuotientInverse(a: FpX, T: FpX, p: bigint, word: boolean): FpX {
  a = trimPolynomial(a);
  T = trimPolynomial(T);
  const error = () => {
    throw new PariError(
      `impossible inverse in ${word ? 'Flxq_inv' : 'FpXQ_inv'}: ${displayExtensionError(a, 0)}.`
    );
  };
  const [g, , v] = word ? _Flx_extgcd(T, a, p, false) : _FpX_extgcd(T, a, p, false);
  if (g.length !== 1) return error();
  let inverse: bigint;
  if (word) inverse = inverseCoefficient(g[0]!, p, true);
  else {
    // FpXQ_invsafe uses Fp_invsafe here; Flxq_invsafe instead calls Fl_inv.
    try {
      inverse = inverseCoefficient(g[0]!, p, false);
    } catch (e) {
      if (e instanceof PariError) return error();
      throw e;
    }
  }
  return FpX_Fp_mul(v, inverse, p);
}

export function polynomialQuotientPower(a: FpX, n: bigint, T: FpX, p: bigint): FpX {
  if (n === 0n) return [1n];
  a = trimPolynomial(a);
  T = trimPolynomial(T);
  if (n === 1n) return FpX_rem(FpX_red(a, p), T, p);
  if (n === -1n) return polynomialQuotientInverse(a, T, p, false);
  const word = p > 0n && p < 1n << 63n;
  if (word) {
    a = FpX_red(a, p);
    T = FpX_red(T, p);
  }
  const base = n < 0n ? polynomialQuotientInverse(a, T, p, word) : a;
  // Native signed powering prepares its reducer after inversion, even for n=2.
  const context = polynomialQuotient(T, p, word);
  return gen_pow_i(
    base,
    n,
    (a) => context.reduce(word ? Flx_sqr(a, p) : FpX_sqr(a, p)),
    (a, b) => context.reduce(context.multiply(a, b))
  );
}
