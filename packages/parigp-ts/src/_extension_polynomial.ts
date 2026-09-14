/** Native FpXX.c / FlxX.c / F2x.c / polarit3.c extension-polynomial kernels.
 * Coefficient tags preserve PARI integer versus polynomial fast paths.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
import { Flx_mul, Flx_sqr, Flx_rem } from './Flx.js';
import { FpX_mul, FpX_rem } from './FpX.js';
import { FpX_red } from './ffinit.js';
import { ZX_mul, ZX_sqr } from './ZX.js';
import { F2x_mul, F2xq_sqr, F2x_rem, F2x_degree } from './F2x.js';
import { trimPolynomial, packUnsigned, unpackUnsigned } from './_polynomial_packing.js';
import { PariError } from './errors.js';
import { extensionField, trimExtension as trim } from './_extension_field.js';
export type ExtensionCoefficient = bigint | bigint[];
export type ExtensionPolynomial = ExtensionCoefficient[];
type Coeff = ExtensionCoefficient;
type Outer = ExtensionPolynomial;
const mod = (x: bigint, p: bigint) => ((x % p) + p) % p;
function degreeT(T: bigint[]): number {
  const d = trimPolynomial(T).length - 1;
  if (d < 1) throw new RangeError('extension modulus must have positive degree');
  return d;
}
function pack(a: Outer, d: number, word: boolean): bigint[] {
  const w = 2 * d - 1,
    out: bigint[] = [];
  for (const [i, c] of a.entries()) {
    const f = typeof c === 'bigint' ? [c] : trimPolynomial(c);
    if (typeof c !== 'bigint' && f.length > d)
      throw new PariError(
        `bug in ${word ? 'zxX' : 'RgXX'}_to_Kronecker, P is not reduced mod Q, please report.`
      );
    while (out.length < i * w) out.push(0n);
    // Append without a call-argument spread: V8 has a much smaller argument limit.
    for (const coefficient of f) out.push(coefficient);
  }
  return trimPolynomial(out);
}
function unpack(z: bigint[], T: bigint[], p: bigint, word: boolean, reduce?: (a:bigint[])=>bigint[]): Outer {
  const width = 2 * degreeT(T) - 1,
    out: Outer = [];
  for (let i = 0; i < z.length; i += width)
    out.push(reduce ? reduce(z.slice(i,i+width)) : word ? Flx_rem(z.slice(i, i + width), T, p) : FpX_rem(z.slice(i, i + width), T, p));
  return trim(out);
}
function mixed(a: Outer, b: bigint[], T: bigint[], p: bigint): Outer {
  const columns: bigint[][] = [];
  for (let j = 0; j < degreeT(T); j++) {
    const c = trimPolynomial(a.map((v) => (typeof v === 'bigint' ? (j ? 0n : v) : (v[j] ?? 0n))));
    columns.push(FpX_mul(b, c, p));
  }
  const n = columns.reduce((n, c) => Math.max(n, c.length), 0);
  return trim(Array.from({ length: n }, (_, i) => trimPolynomial(columns.map((c) => c[i] ?? 0n))));
}
/** Common representation guards, without changing native coefficient tags.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function validateExtensionInputs(
  mode: 0 | 1 | 2,
  p: bigint,
  T: bigint[] | bigint,
  a: Outer,
  b: Outer = []
): void {
  if (mode === 2) {
    if ((T as bigint) < 2n) throw new RangeError('extension modulus must have positive degree');
    if ([...a, ...b].some((c) => (c as bigint) < 0n))
      throw new RangeError('polynomial bits must be nonnegative');
  } else {
    if (p < 2n) throw new RangeError('extension arithmetic requires a modulus greater than one');
    degreeT(T as bigint[]);
    if (mode === 1) {
      if (p >= 1n << 64n) throw new RangeError('modulus must be a positive word integer');
      if ([...(T as bigint[]), ...a.flat(), ...b.flat()].some((c) => c < 0n || c >= p))
        throw new RangeError('word polynomial coefficients must be reduced');
    }
  }
}
/**
 * F2xX_to_Kronecker: XOR shifted coefficient words into the native allocation.
 * Raw coefficients may overlap packing blocks; preserve that source behavior.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
function packBinary(a: bigint[], d: number): bigint {
  if (!a.length) return 0n;
  const width = 2 * d + 1,
    count = Math.ceil(((a.length - 1) * width + d + 1) / 64),
    capacity = count * 64;
  let overlap = false;
  for (let i = 0; i < a.length; i++) {
    const degree = F2x_degree(a[i]!);
    if (degree >= capacity - i * width)
      throw new RangeError('binary coefficient degree exceeds the native packing bound');
    if (degree >= width) overlap = true;
  }
  if (!overlap) return packUnsigned(a, width);
  const words = new Array<bigint>(count).fill(0n),
    mask = (1n << 64n) - 1n;
  for (let i = 0; i < a.length; i++) {
    const hex = a[i]!.toString(16),
      offset = i * width,
      shift = BigInt(offset % 64);
    let at = Math.floor(offset / 64);
    for (let end = hex.length; end > 0; end -= 16, at++) {
      const value = BigInt('0x' + hex.slice(Math.max(0, end - 16), end)) << shift;
      words[at] = words[at]! ^ (value & mask);
      const carry = value >> 64n;
      if (carry) words[at + 1] = words[at + 1]! ^ carry;
    }
  }
  return packUnsigned(words, 64);
}
/** Shared native operation dispatch: product, square, reduction, normalization. */
export function extensionPolynomial(
  mode: 0 | 1 | 2,
  op: 0 | 1 | 2 | 3,
  p: bigint,
  T: bigint[] | bigint,
  a: Outer,
  b: Outer = [],
  innerInverse?: bigint[]
): Outer {
  validateExtensionInputs(mode, p, T, a, b);

  a = trim(a.map((c) => (typeof c === 'bigint' ? c : trimPolynomial(c))));
  b = trim(b.map((c) => (typeof c === 'bigint' ? c : trimPolynomial(c))));
  if (mode === 2) {
    const modulus = T as bigint,
      aa = a as bigint[],
      bb = b as bigint[],
      d = F2x_degree(modulus);
    if (op === 3) {
      if (!aa.length) return [];
      const field = extensionField(2, T, p),
        inverse = field.inv(aa.at(-1)!);
      return [...aa.slice(0, -1).map((c) => field.mul(c, inverse)), 1n];
    }
    if (op === 2) return trim(aa.map((c) => F2x_rem(c, modulus)));
    if (op === 1) {
      if (!aa.length) return [];
      const out = new Array<bigint>(2 * aa.length - 1).fill(0n);
      for (let i = 0; i < aa.length; i++) out[2 * i] = F2xq_sqr(aa[i]!, modulus);
      return trim(out);
    }
    const width = 2 * d + 1,
      packed = F2x_mul(packBinary(bb, d), packBinary(aa, d));
    if (!packed) return [];
    const count = Math.ceil((F2x_degree(packed) + 1) / width);
    return trim(unpackUnsigned(packed, width, count).map((c) => F2x_rem(c, modulus)));
  }
  const modulus = T as bigint[],
    word = mode === 1,
    cached = innerInverse === undefined ? undefined : extensionField(mode,T,p,innerInverse);
  if (op === 3) {
    if (!a.length) return [];
    const field = cached ?? extensionField(mode, T, p),
      lead = a.at(-1)!;
    let inverse: Coeff;
    if (!word && (typeof lead === 'bigint' || lead.length < 2)) {
      const scalar = typeof lead === 'bigint' ? lead : lead[0]!;
      if (scalar === 1n) return [...a.slice(0, -1), 1n];
      inverse = field.inv(scalar);
    } else inverse = field.inv(lead);
    return [...a.slice(0, -1).map((c) => field.mul(inverse, c)), word ? [1n] : 1n];
  }
  if (op === 2 && cached) return cached.pred(a);
  if (op === 2)
    return trim(
      a.map((c) =>
        word
          ? Flx_rem(c as bigint[], modulus, p)
          : typeof c === 'bigint'
            ? mod(c, p)
            : FpX_rem(FpX_red(c, p), modulus, p)
      )
    );
  const absolute = (v: Outer) => v.every((c) => typeof c === 'bigint');
  if (!word) {
    if (op === 1 && absolute(a)) return ZX_sqr(a as bigint[]);
    if (op === 0) {
      if (absolute(b))
        return absolute(a)
          ? FpX_mul(a as bigint[], b as bigint[], p)
          : mixed(a, b as bigint[], modulus, p);
      if (absolute(a)) return mixed(b, a as bigint[], modulus, p);
    }
  }
  const d = degreeT(modulus),
    x = pack(a, d, word);
  const result =
    op === 1
      ? word
        ? Flx_sqr(x, p)
        : ZX_sqr(x)
      : word
        ? Flx_mul(pack(b, d, word), x, p)
        : ZX_mul(pack(b, d, word), x);
  return unpack(word ? result : FpX_red(result, p), modulus, p, word, cached ? a => cached.red(a) as bigint[] : undefined);
}


/** Native low products: generic packing truncates before coefficient reduction;
 * word products reduce the full product before truncating outer coefficients.
 * @see Deviation: PARI extension projection adapters
 */
export function extensionTruncatedPolynomial(
  mode: 0 | 1,
  square: boolean,
  n: number,
  T: bigint[],
  p: bigint,
  a: Outer,
  b: Outer = [],
  innerInverse?: bigint[]
): Outer {
  if (!Number.isSafeInteger(n) || n < 0)
    throw new RangeError('truncation length must be nonnegative');
  validateExtensionInputs(mode, p, T, a, square ? [] : b);
  a = trim(a.map((c) => (Array.isArray(c) ? trimPolynomial(c) : c)));
  b = square ? [] : trim(b.map((c) => (Array.isArray(c) ? trimPolynomial(c) : c)));
  if (mode === 1)
    return trim(extensionPolynomial(mode, square ? 1 : 0, p, T, a, b, innerInverse).slice(0, n));
  if (a.every((c) => typeof c === 'bigint') && (square || b.every((c) => typeof c === 'bigint'))) {
    const r = (square ? ZX_sqr(a as bigint[]) : ZX_mul(a as bigint[], b as bigint[])).slice(0, n);
    return square ? trimPolynomial(r) : FpX_red(r, p);
  }
  const d = degreeT(T),
    x = pack(a, d, false);
  const r = (square ? ZX_sqr(x) : ZX_mul(pack(b, d, false), x)).slice(0, (2 * d - 1) * n);
  const cached = innerInverse === undefined ? undefined : extensionField(mode, T, p, innerInverse);
  return unpack(FpX_red(r, p), T, p, false, cached ? (a) => cached.red(a) as bigint[] : undefined);
}
