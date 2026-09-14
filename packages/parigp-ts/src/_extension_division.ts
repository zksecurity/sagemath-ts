/** Native extension-polynomial basecase and Newton/Barrett division.
 * @see Deviation: PARI extension-polynomial division adapters
 */
import { validateExtensionInputs, type ExtensionPolynomial as P } from './_extension_polynomial.js';
import {
  extensionField as context,
  extensionZero as zero,
  extensionOne as one,
  trimExtension as trim,
  type ExtensionField as Context,
} from './_extension_field.js';
import { FpX_red } from './ffinit.js';
import { F2x_degree } from './F2x.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { quadratic_prec_mask } from './Zp.js';
import { PariError } from './errors.js';
/** GEN [reciprocal, polynomial] storage, represented without ambiguous array nesting. */
export interface ExtensionReduction<A extends P = P> {
  polynomial: A;
  inverse: A;
}
export type ExtensionModulus<A extends P = P> = A | ExtensionReduction<A>;
function base(a: P, b: P, c: Context, op: number): P | [P, P] | null {
  const dx = a.length - 1,
    dy = b.length - 1;
  if (dy < 0)
    throw new PariError(
      `impossible inverse in ${c.mode === 0 ? 'FpX_divrem' : c.mode === 1 ? 'FlxqX_divrem' : 'F2xqX_divrem'}: 0.`
    );
  if (dx < dy) {
    const r = op === 4 ? [] : c.pred(a);
    return op === 0 ? [[], r] : op === 1 ? r : op === 3 && r.length ? null : [];
  }
  const lead = b[dy]!;
  if (dy === 0) {
    if (op === 1) return [];
    const q = one(lead) ? (c.mode === 0 ? c.pred(a) : a.slice()) : c.pscale(a, c.inv(lead));
    return op === 0 ? [q, []] : q;
  }
  const ci = one(lead) ? null : c.inv(lead),
    dz = dx - dy,
    q: P = new Array(dz + 1).fill(c.z);
  let dy1 = dy - 1;
  if (c.mode === 0) while (dy1 >= 0 && zero(b[dy1]!)) dy1--;
  q[dz] = ci === null ? a[dx]! : c.mul(a[dx]!, ci);
  for (let i = dx - 1; i >= dy; i--) {
    let v = a[i]!;
    for (let j = i - dy1; j <= i && j <= dz; j++) v = c.sub(v, c.mul(q[j]!, b[i - j]!, false));
    if (ci !== null) v = c.mul(v, ci, false);
    q[i - dy] = c.red(v);
  }
  if (op === 4) return q;
  const r: P = new Array(dy).fill(c.z);
  for (let i = dy - 1; i >= 0; i--) {
    let v = a[i]!;
    for (let j = c.mode === 0 ? Math.max(0, i - dy1) : 0; j <= i && j <= dz; j++)
      v = c.sub(v, c.mul(q[j]!, b[i - j]!, false));
    r[i] = c.red(v);
    if (op === 3 && !zero(r[i]!)) return null;
  }
  return op === 0 ? [q, trim(r)] : op === 1 ? trim(r) : q;
}
function reciprocal(a: P, n: number, c: Context): P {
  return trim(Array.from({ length: n }, (_, i) => a[n - 1 - i] ?? c.z));
}
function inverseBarrett(a: P, c: Context): P {
  const degree = a.length - 1;
  if (degree < 2) return [];
  const limit = c.mode === 0 ? 40 : c.mode === 1 ? 22 : 50;
  if (a.length + 2 <= limit) {
    const lead = a[degree]!,
      ci = one(lead) ? null : c.inv(lead);
    const t = ci === null ? a : c.pscale(a, ci);
    const r: P = [c.u];
    for (let i = 1; i < degree - 1; i++) {
      if (c.mode === 0) {
        let v = t[degree - i]!;
        for (let k = 1; k < i; k++) v = c.add(v, c.mul(t[degree - i + k]!, r[k]!, false));
        r[i] = c.red(c.neg(v));
      } else {
        // Word/binary basecases reduce products, but preserve the raw initial term.
        let v = c.neg(t[degree - i]!);
        for (let k = 1; k < i; k++) v = c.sub(v, c.mul(t[degree - i + k]!, r[k]!));
        r[i] = v;
      }
    }
    return ci === null ? trim(r) : c.pscale(r, ci);
  }
  let mask = quadratic_prec_mask(degree - 2);
  const q = reciprocal(a, a.length, c),
    x: P = new Array(degree).fill(c.z);
  x[0] = c.inv(q[0]!);
  if (
    q.length > 1 &&
    (c.mode === 0 ||
      (c.mode === 1
        ? (q[1] as bigint[]).length >= trimPolynomial(c.T as bigint[]).length
        : F2x_degree(q[1] as bigint) >= F2x_degree(c.T as bigint)))
  )
    q[1] = c.red(q[1]!);
  let lx = 1;
  if (q.length > 1 && !zero(q[1]!)) {
    let v = q[1]!;
    if (!one(x[0]!)) v = c.mul(v, c.sqr(x[0]!));
    x[1] = c.neg(v);
    lx = 2;
  }
  let old = 1;
  while (mask > 1n) {
    let next = old * 2;
    if (mask & 1n) next--;
    mask >>= 1n;
    const n = next + 1;
    let z = c.pmul(x.slice(0, lx), trim(q.slice(0, n))).slice(0, n),
      i = old;
    while (i < z.length && zero(z[i]!)) i++;
    old = next;
    if (i >= z.length) continue;
    z = c.pmul(x.slice(0, lx), trim(z.slice(i)));
    if (z.length > n - i) z = trim(z.slice(0, n - i));
    lx = z.length + i;
    for (let j = 0; j < z.length; j++) x[i + j] = c.neg(z[j]!);
  }
  return x.slice(0, lx);
}
function barrettSpec(a: P, mg: P, b: P, c: Context, wantR: boolean): [P, P] {
  const lt = b.length - 1,
    ld = a.length - lt,
    lm = Math.min(ld, mg.length);
  let q = reciprocal(a.slice(lt), ld, c);
  q = c.pmul(q, trim(mg.slice(0, lm)));
  q = reciprocal(q.slice(0, ld), ld, c);
  if (!wantR) return [q, []];
  const r = c.psub(a.slice(0, lt), c.pmul(q, trim(b.slice(0, lt))).slice(0, lt));
  return [q, r];
}
function barrett(a: P, mg: P, b: P, c: Context, op: number): P | [P, P] | null {
  // Keep one backing array and move a logical length, as the native chunk loop does.
  // F2x.c retains lg(r), not l+2, in its final short remainder branch.
  let r = c.pred(a),
    l = r.length;
  const lt = b.length - 1,
    lm = 2 * lt - 1;
  if (l <= lt) return op === 0 ? [[], r] : op === 1 ? r : op === 3 && r.length ? null : [];
  if (lt <= 1) return base(c.mode === 0 ? r : a, b, c, op);
  let q: P | null = op !== 1 && l > lm ? new Array(l - lt).fill(c.z) : null;
  while (l > lm) {
    const offset = l - lm,
      [zq, zr] = barrettSpec(r.slice(offset, l), mg, b, c, true);
    if (q) for (let i = 0; i < zq.length; i++) q[offset + i] = zq[i]!;
    for (let i = 0; i < zr.length; i++) r[offset + i] = zr[i]!;
    l = offset + zr.length;
  }
  if (op === 1)
    return l > lt
      ? barrettSpec(r.slice(0, l), mg, b, c, true)[1]
      : trim(c.mode === 2 ? r : r.slice(0, l));
  if (l > lt) {
    const [zq, zr] = barrettSpec(r.slice(0, l), mg, b, c, op !== 4);
    if (q) for (let i = 0; i < zq.length; i++) q[i] = zq[i]!;
    else q = zq;
    if (op !== 4) r = zr;
  } else if (op !== 4) r = trim(r.slice(0, l));
  q = trim(q!);
  return op === 0 ? [q, r] : op === 3 && r.length ? null : q;
}
/** Internal output-pointer modes: pair, remainder, inverse, divides, quotient. */
export function extensionDivision(
  mode: 0 | 1 | 2,
  op: 0 | 1 | 2 | 3 | 4,
  p: bigint,
  T: bigint[] | bigint,
  a: P,
  source: ExtensionModulus = [],
  innerInverse?: bigint[]
): P | [P, P] | null {
  const b0 = Array.isArray(source) ? source : source.polynomial;
  validateExtensionInputs(mode, p, T, a, b0);
  a = trim(a.map((c) => (typeof c === 'bigint' ? c : trimPolynomial(c))));
  const b = trim(b0.map((c) => (typeof c === 'bigint' ? c : trimPolynomial(c))));
  const c = context(mode, T, p, innerInverse);
  if (op === 2) return inverseBarrett(a, c);
  const mg = Array.isArray(source) ? null : source.inverse;
  if (op === 1 && a.length < b.length) return c.pred(a);
  if (mode === 0 && p < 1n << 64n) {
    const convert = (v: P) => trim(v.map((c) => FpX_red(typeof c === 'bigint' ? [c] : c, p)));
    const r = extensionDivision(1, op, p, FpX_red(T as bigint[], p), convert(a), convert(b), innerInverse === undefined ? undefined : FpX_red(innerInverse,p));
    const collapse = (v: P) =>
      v.map((c) => ((c as bigint[]).length < 2 ? ((c as bigint[])[0] ?? 0n) : c));
    return r === null
      ? null
      : op === 0
        ? [collapse(r[0] as P), collapse(r[1] as P)]
        : collapse(r as P);
  }
  const limit = mode === 0 ? 30 : mode === 1 ? (op === 1 ? 48 : 46) : op === 1 ? 101 : 97;
  return mg === null && a.length - b.length + 3 < limit
    ? base(a, b, c, op)
    : barrett(a, mg ?? inverseBarrett(b, c), b, c, op);
}
/** Native get_red threshold; an existing reduction object is reused. */
export function extensionGetRed<A extends P>(
  mode: 0 | 1 | 2,
  S: ExtensionModulus<A>,
  T: bigint[] | bigint,
  p: bigint,
  innerInverse?: bigint[]
): ExtensionModulus<A> {
  const raw = Array.isArray(S) ? S : S.polynomial;
  validateExtensionInputs(mode, p, T, raw);
  if (!Array.isArray(S)) return S;
  const polynomial = trim(S.map((c) => (typeof c === 'bigint' ? c : trimPolynomial(c)))) as A;
  return polynomial.length + 2 > (mode === 0 ? 12 : mode === 1 ? 17 : 48)
    ? { polynomial, inverse: inverseBarrett(polynomial, context(mode, T, p, innerInverse)) as A }
    : polynomial;
}
