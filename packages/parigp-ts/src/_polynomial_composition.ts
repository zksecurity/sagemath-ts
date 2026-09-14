import { FpX_sqr } from './FpX.js';
import { Flx_mul, Flx_sqr } from './Flx.js';
import { polynomialQuotient } from './_polynomial_quotient.js';
import { trimPolynomial } from './_polynomial_packing.js';
/** Common FpX/Flx schedules; field coefficients use exact BigInt arithmetic. */
import { FpM_mul, Flm_mul } from './FpV.js';
import { FpX_add, FpX_mul, FpX_renormalize, type FpX } from './ffinit.js';

/** bb_group.c gen_powers: unreduced initial entries and native squaring schedule. */
export function quotientPowers(
  x: FpX,
  n: number,
  T: FpX,
  p: bigint,
  reduce?: (a: FpX) => FpX,
  word = false
): FpX[] {
  if (!Number.isSafeInteger(n) || n < 0) throw new RangeError('power count must be nonnegative');
  x = trimPolynomial(x);
  T = trimPolynomial(T);
  reduce ??= polynomialQuotient(T, p, word).reduce;
  const out: FpX[] = [[1n]];
  if (!n) return out;
  out.push(x.slice());
  const useSquare = 2 * (x.length - 1) >= T.length - 1;
  for (let i = 2; i <= n; i++) {
    const square = i === 2 || (useSquare && i % 2 === 0),
      a = square ? out[i / 2]! : out[i - 1]!;
    out.push(
      reduce(
        square ? (word ? Flx_sqr(a, p) : FpX_sqr(a, p)) : word ? Flx_mul(a, x, p) : FpX_mul(a, x, p)
      )
    );
  }
  return out;
}

/** FpX.c/Flx.c blocked matrix evaluation followed by giant-step Horner folding. */
export function evaluatePowers(
  Q: FpX,
  powers: FpX[],
  T: FpX,
  p: bigint,
  word: boolean,
  reduce?: (a: FpX) => FpX
): FpX {
  Q = trimPolynomial(Q);
  T = trimPolynomial(T);
  powers = powers.map(trimPolynomial);
  if (!Q.length) return [];
  const m = T.length - 1,
    l = powers.length;
  if (m < 0) throw new RangeError('polynomial modulus must be nonzero');
  if (!l || (l === 1 && Q.length > 1)) throw new RangeError('power table is too short');
  const n = Q.length <= l ? l : l - 1,
    d = Q.length <= l ? 1 : Math.ceil(Q.length / n);
  // Flx_to_Flv does not truncate: excess coefficients overwrite native memory.
  if (word && powers.slice(0, n).some((f) => f.length > m))
    throw new RangeError('power polynomial exceeds modulus degree');
  const a = [
    [],
    ...powers.slice(0, n).map((f) => [0n, ...Array.from({ length: m }, (_, i) => f[i] ?? 0n)]),
  ];
  const b = [
    [],
    ...Array.from({ length: d }, (_, i) => [
      0n,
      ...Array.from({ length: n }, (_, j) => Q[i * n + j] ?? 0n),
    ]),
  ];
  const c = word ? Flm_mul(a, b, p) : FpM_mul(a, b, p),
    g = powers[l - 1]!;
  reduce ??= polynomialQuotient(T, p, word).reduce;
  let s = FpX_renormalize(c[d]!.slice(1));
  for (let i = d - 1; i > 0; i--)
    s = FpX_add(
      reduce(word ? Flx_mul(s, g, p) : FpX_mul(s, g, p)),
      FpX_renormalize(c[i]!.slice(1)),
      p
    );
  return s;
}
