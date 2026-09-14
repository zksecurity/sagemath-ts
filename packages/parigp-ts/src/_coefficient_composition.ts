/** Native FpXX.c / FlxX.c / F2x.c coefficient substitution.
 * Integer coefficient tags bypass field reduction; polynomial tags use the
 * original public table evaluator unless a caller supplied a prepared cache.
 * @see Deviation: PARI coefficient-substitution adapters
 */
import { FpXQ_powers, FpX_FpXQV_eval } from './FpX.js';
import { Flxq_powers, Flx_FlxqV_eval } from './Flx.js';
import { F2x_degree, F2xq_powers } from './F2x.js';
import { F2x_F2xqV_eval } from './F2x.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { trimExtension } from './_extension_field.js';
import { brent_kung_optpow } from './RgX.js';
import { polynomialQuotient } from './_polynomial_quotient.js';
import type { ExtensionPolynomial as P } from './_extension_polynomial.js';
export function coefficientSubstitution(
  mode: 0 | 1 | 2,
  op: 0 | 1,
  p: bigint,
  T: bigint[] | bigint,
  P: P,
  x: bigint[] | bigint,
  V: bigint[][] | bigint[],
  providedInverse?: bigint[]
): P {
  P = trimExtension(P.map((c) => (Array.isArray(c) ? trimPolynomial(c) : c)));
  if (mode !== 2) T = trimPolynomial(T as bigint[]);
  if (mode === 1 && (p < 1n || p >= 1n << 64n))
    throw new RangeError('modulus must be a positive word integer');
  const ctx =
    providedInverse === undefined
      ? undefined
      : polynomialQuotient(T as bigint[], p, mode === 1, providedInverse);
  if (op === 0) {
    const n = brent_kung_optpow(
      mode === 2 ? F2x_degree(T as bigint) - 1 : (T as bigint[]).length - 2,
      P.length,
      1
    );
    V =
      mode === 2
        ? F2xq_powers(x as bigint, n, T as bigint)
        : ctx
          ? ctx.powers(x as bigint[], n)
          : mode === 0
            ? FpXQ_powers(x as bigint[], n, T as bigint[], p)
            : Flxq_powers(x as bigint[], n, T as bigint[], p);
  }
  return trimExtension(
    P.map((c) =>
      mode === 0 && typeof c === 'bigint'
        ? c
        : mode === 2
          ? F2x_F2xqV_eval(c as bigint, V as bigint[], T as bigint)
          : ctx
            ? ctx.evaluate(c as bigint[], V as bigint[][])
            : mode === 0
              ? FpX_FpXQV_eval(c as bigint[], V as bigint[][], T as bigint[], p)
              : Flx_FlxqV_eval(c as bigint[], V as bigint[][], T as bigint[], p)
    )
  );
}
