import { Flx_divrem, Flx_invBarrett, Flx_mul, Flx_rem, Flx_sqr } from './Flx.js';
import { FpX_divrem, FpX_invBarrett, FpX_mul, FpX_rem, FpX_sqr } from './FpX.js';
import { evaluatePowers, quotientPowers } from './_polynomial_composition.js';
import { divisionBarrett } from './_polynomial_division.js';
import { gen_pow_i } from './bb_group.js';
/** Reusable PARI FpX_get_red / Flx_get_red contexts, with native GMP cutoffs.
 * @see Deviation: PARI signed composition and matrix boundaries
 */
import { type FpX, FpX_red } from './ffinit.js';

export function polynomialQuotient(
  T: FpX,
  p: bigint,
  word: boolean,
  preinverse?: FpX
): PolynomialQuotient {
  const multiply = (a: FpX, b: FpX) => (word ? Flx_mul(a, b, p) : FpX_mul(a, b, p));
  const square = (a: FpX) => (word ? Flx_sqr(a, p) : FpX_sqr(a, p));
  const divide = (a: FpX, b: FpX) => (word ? Flx_divrem(a, b, p) : FpX_divrem(a, b, p));
  const cached = word ? T.length >= (p <= 3037000493n ? 92 : 31) : T.length + 2 > 38;
  const inverse =
    preinverse ?? (cached ? (word ? Flx_invBarrett(T, p) : FpX_invBarrett(T, p)) : undefined);
  const reduce = (a: FpX): FpX =>
    inverse !== undefined
      ? divisionBarrett(a, T, inverse, p, multiply, divide)[1]
      : word
        ? Flx_rem(a, T, p)
        : FpX_rem(a, T, p);
  const pp = p < 0n ? -p : p;
  const asWord = () =>
    polynomialQuotient(
      FpX_red(T, pp),
      pp,
      true,
      inverse === undefined ? undefined : FpX_red(inverse, pp)
    );
  const powers = (x: FpX, count: number): FpX[] =>
    !word && pp > 0n && pp < 1n << 64n && count > 2
      ? asWord().powers(FpX_red(x, pp), count)
      : quotientPowers(x, count, T, p, reduce, word);
  const evaluate = (x: FpX, table: FpX[]) => evaluatePowers(x, table, T, p, word, reduce);
  // These contexts serve Frobenius steps over a field, so p is at least two.
  const power = (x: FpX): FpX =>
    !word && p < 1n << 63n
      ? asWord().power(FpX_red(x, p))
      : gen_pow_i(
          x,
          p,
          (a) => reduce(square(a)),
          (a, b) => reduce(multiply(a, b))
        );
  return { T, p, word, inverse, multiply, divide, reduce, powers, evaluate, power };
}
export interface PolynomialQuotient {
  T: FpX;
  p: bigint;
  word: boolean;
  inverse: FpX | undefined;
  multiply(a: FpX, b: FpX): FpX;
  divide(a: FpX, b: FpX): [FpX, FpX];
  reduce(a: FpX): FpX;
  powers(a: FpX, count: number): FpX[];
  evaluate(a: FpX, table: FpX[]): FpX;
  power(a: FpX): FpX;
}
