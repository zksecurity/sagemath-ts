/** Shared arithmetic schedules from PARI FpX.c and Flx.c. */
import { PariError } from './errors.js';
import { pariErrorPayload } from './_error_display.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { quadratic_prec_mask } from './Zp.js';
/** @see Deviation: PARI signed scalar inverse and division residues */
export const residue = (x: bigint, p: bigint): bigint => {
  if (x === 0n) return 0n;
  if (p === 0n) throw new PariError('impossible inverse in dvmdii: 0.');
  if (p < 0n) p = -p;
  const r = x % p;
  return r < 0n ? r + p : r;
};
/** @see Deviation: PARI signed scalar inverse and division residues */
export function inverseCoefficient(a: bigint, p: bigint, word: boolean): bigint {
  const modulus = p < 0n ? -p : p;
  let r = modulus === 0n ? (a < 0n ? -a : a) : residue(a, modulus),
    s = modulus,
    u = 1n,
    v = 0n;
  while (s) {
    const q = r / s;
    [r, s] = [s, r - q * s];
    [u, v] = [v, u - q * v];
  }
  if (modulus === 0n || (r !== 1n && modulus !== 1n)) {
    const value = 'Mod(' + (word ? a : r) + ', ' + p + ')';
    throw new PariError(
      'impossible inverse in ' +
        (word ? 'Fl_inv' : 'Fp_inv') +
        ': ' +
        pariErrorPayload(value, 't_INTMOD') +
        '.'
    );
  }
  return residue(u, modulus);
}
export type PolynomialMultiply = (a: bigint[], b: bigint[]) => bigint[];

/** Back substitution, using the native last-nonzero lower divisor coefficient. */
export function divisionBasecase(
  a: bigint[],
  b: bigint[],
  p: bigint,
  word: boolean
): [bigint[], bigint[]] {
  const dy = b.length - 1,
    dx = a.length - 1,
    dz = dx - dy;
  const lead = b[dy]!,
    inv = lead === 1n ? 1n : inverseCoefficient(lead, p, word);
  let dy1 = dy - 1;
  while (dy1 >= 0 && b[dy1] === 0n) dy1--;
  const q = new Array<bigint>(dz + 1).fill(0n);
  q[dz] = !word && lead === 1n ? a[dx]! : residue(a[dx]! * inv, p);
  for (let i = dx - 1; i >= dy; i--) {
    let c = a[i]!;
    for (let j = i - dy1; j <= i && j <= dz; j++) c -= q[j]! * b[i - j]!;
    q[i - dy] = residue(c * inv, p);
  }
  const r = new Array<bigint>(dy);
  for (let i = 0; i < dy; i++) {
    let c = a[i]!;
    for (let j = Math.max(0, i - dy1); j <= i && j <= dz; j++) c -= q[j]! * b[i - j]!;
    r[i] = residue(c, p);
  }
  return [trimPolynomial(q), trimPolynomial(r)];
}

/** 1/reciprocal(T) modulo x^(degree(T)-1), native basecase/Newton schedule. */
export function inverseBarrett(
  T: bigint[],
  p: bigint,
  word: boolean,
  basecase: boolean,
  multiply: PolynomialMultiply
): bigint[] {
  const degree = T.length - 1;
  if (degree < 2) return [];
  const lead = T[degree]!;
  if (basecase) {
    const ci = lead === 1n ? 1n : inverseCoefficient(lead, p, word);
    const t = lead === 1n ? T : T.map((c) => residue(c * ci, p));
    const r = [1n];
    for (let i = 1; i < degree - 1; i++) {
      let u = t[degree - i]!;
      for (let k = 1; k < i; k++) u += t[degree - i + k]! * r[k]!;
      r[i] = residue(-u, p);
    }
    return trimPolynomial(lead === 1n ? r : r.map((c) => residue(c * ci, p)));
  }
  const q = trimPolynomial(T.slice().reverse()),
    x = new Array<bigint>(degree).fill(0n);
  x[0] = inverseCoefficient(q[0]!, p, word);
  let lx = 1;
  if (q.length > 1 && residue(q[1]!, p)) {
    x[1] = residue(-q[1]! * x[0]! * x[0]!, p);
    lx = 2;
  }
  let nold = 1,
    mask = quadratic_prec_mask(degree - 2);
  while (mask > 1n) {
    const nnew = 2 * nold - Number(mask & 1n),
      length = nnew + 1;
    mask >>= 1n;
    let z = multiply(x.slice(0, lx), trimPolynomial(q.slice(0, length))).slice(0, length);
    let i = nold;
    while (i < z.length && z[i] === 0n) i++;
    nold = nnew;
    if (i >= z.length) continue;
    z = trimPolynomial(multiply(x.slice(0, lx), trimPolynomial(z.slice(i))).slice(0, length - i));
    lx = i + z.length;
    for (let j = 0; j < z.length; j++) x[i + j] = residue(-z[j]!, p);
  }
  // Native Newton returns its logical length, including [0] when p = 1.
  return x.slice(0, lx);
}

/** One native Barrett block, with the leading term of T omitted from the product. */
function divideBlock(
  a: bigint[],
  T: bigint[],
  inverse: bigint[],
  p: bigint,
  multiply: PolynomialMultiply
): [bigint[], bigint[]] {
  const degree = T.length - 1,
    length = a.length - degree;
  const high = trimPolynomial(a.slice(degree).reverse());
  const product = multiply(high, trimPolynomial(inverse.slice(0, length)));
  const q = trimPolynomial(Array.from({ length }, (_, i) => product[length - 1 - i] ?? 0n));
  const low = multiply(q, trimPolynomial(T.slice(0, degree)));
  const r = trimPolynomial(
    Array.from({ length: degree }, (_, i) => residue((a[i] ?? 0n) - (low[i] ?? 0n), p))
  );
  return [q, r];
}

/** Descending blocks of at most 2*degree(T)-1 coefficients, as in PARI. */
export function divisionBarrett(
  a: bigint[],
  T: bigint[],
  inverse: bigint[],
  p: bigint,
  multiply: PolynomialMultiply,
  basecase: (a: bigint[], b: bigint[]) => [bigint[], bigint[]]
): [bigint[], bigint[]] {
  let r = trimPolynomial(a.map((c) => residue(c, p)));
  const degree = T.length - 1;
  if (r.length <= degree) return [[], r];
  if (degree <= 1) return basecase(r, T);
  const block = 2 * degree - 1,
    q = new Array<bigint>(r.length - degree).fill(0n);
  while (r.length > block) {
    const offset = r.length - block;
    const [zq, zr] = divideBlock(r.slice(offset), T, inverse, p, multiply);
    for (let i = 0; i < zq.length; i++) q[offset + i] = zq[i]!;
    r.length = offset;
    for (const c of zr) r.push(c);
  }
  if (r.length > degree) {
    const [zq, zr] = divideBlock(r, T, inverse, p, multiply);
    for (let i = 0; i < zq.length; i++) q[i] = zq[i]!;
    r = zr;
  }
  return [trimPolynomial(q), trimPolynomial(r)];
}
