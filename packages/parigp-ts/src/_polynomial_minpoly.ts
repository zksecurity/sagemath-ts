/** Shoup's transposed-power projection, PARI FpX.c/Flx.c minpoly. */
import { type FpX, FpX_sub, FpX_normalize } from './ffinit.js';
import { FpX_halfgcd, random_FpX } from './FpX.js';
import { Flx_halfgcd, random_Flx } from './Flx.js';
import { polynomialQuotient, type PolynomialQuotient } from './_polynomial_quotient.js';
import { residue } from './_polynomial_division.js';
import { trimPolynomial } from './_polynomial_packing.js';

/** Native cached Barrett reduction is shared by powers and transposed products. */
export function polynomialMinimalPolynomial(
  x: FpX,
  T: FpX,
  p: bigint,
  word: boolean,
  context?: PolynomialQuotient
): FpX {
  if (p < 2n) throw new RangeError('minimal polynomial requires a modulus greater than one');
  const n = T.length - 1;
  if (n < 1) throw new RangeError('minimal polynomial modulus must have positive degree');
  if (word && x.length > n) throw new RangeError('power polynomial exceeds modulus degree');
  const {
    multiply,
    divide,
    inverse: h,
    reduce,
    powers: makePowers,
    evaluate,
  } = context ?? polynomialQuotient(T, p, word);
  // ZXn_mul / Flxn_mul truncate the ordinary product after multiplication.
  const lowMultiply = (a: FpX, b: FpX, length: number): FpX =>
    trimPolynomial(multiply(a, b).slice(0, length));
  const reverse = (a: FpX, length: number): FpX =>
    trimPolynomial(Array.from({ length }, (_, i) => a[length - 1 - i] ?? 0n));
  const shift = (a: FpX, count: number): FpX =>
    !a.length ? [] : count < 0 ? a.slice(-count) : [...new Array<bigint>(count).fill(0n), ...a];
  type Transposed = [FpX, FpX, FpX];
  const transposedInit = (tau: FpX): Transposed => {
    const ft = reverse(T, n + 1),
      bt = reverse(tau, n);
    const bht =
      h !== undefined ? lowMultiply(bt, h, n - 1) : reverse(divide(shift(tau, n - 1), T)[0], n - 1);
    return [bt, bht, ft];
  };
  const transposedMultiply = ([bt, bht, ft]: Transposed, a: FpX): FpX => {
    if (!a.length) return [];
    const t2 = shift(multiply(bt, a), 1 - n);
    if (!bht.length) return t2;
    const t1 = shift(multiply(ft, a), -n);
    return FpX_sub(t2, shift(lowMultiply(t1, bht, n - 1), 1), p);
  };
  if (!word) x = reduce(x);
  const powers = makePowers(x, Math.floor(Math.sqrt(2 * n)));
  let g: FpX = [1n],
    tau: FpX = [1n];
  while (tau.length) {
    if (g.length - 1 === n) {
      tau = [1n];
      g = [1n];
    }
    let v = word ? random_Flx(n, p) : random_FpX(n, p);
    v = transposedMultiply(transposedInit(tau), v);
    const m = 2 * (n - g.length + 1),
      k = Math.floor(Math.sqrt(m));
    const tr = transposedInit(powers[k]!);
    const c = new Array<bigint>(m).fill(0n);
    for (let i = 0; i < m; i += k) {
      for (let j = 0; j < Math.min(m - i, k); j++) {
        const a = powers[j]!;
        let dot = 0n;
        for (let l = 0; l < Math.min(v.length, a.length); l++) dot += v[l]! * a[l]!;
        c[m - 1 - i - j] = residue(dot, p);
      }
      v = transposedMultiply(tr, v);
    }
    const monomial = [...new Array<bigint>(m).fill(0n), 1n];
    const M = word
      ? Flx_halfgcd(monomial, trimPolynomial(c), p)
      : FpX_halfgcd(monomial, trimPolynomial(c), p);
    const next = M[1][1];
    if (next.length < 2) continue;
    g = multiply(g, next);
    tau = reduce(multiply(tau, evaluate(next, powers)));
  }
  return FpX_normalize(g, p);
}
