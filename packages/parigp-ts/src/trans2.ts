import { abpq_init, abpq_sum } from './trans1.js';
import { rdivii } from './kernel/none/level1.js';
/** Factorial and its positive-integer gamma branch, PARI basemath/trans2.c. */
import {
  itor,
  rtor,
  mulrr,
  mulir,
  divrr,
  divri,
  divir,
  sqrr,
  sqrtr_abs as sqrtr,
  shiftr,
  addrr,
  subrr,
  addir,
  addrs,
  logr_abs,
  type MpReal,
} from './qfb.js';
import { rtodbl, invr } from './kernel/none/mp_indep.js';
import { mppi, mpexp, powru } from './trans1.js';
import { bernfrac, constbern } from './bern.js';
import { mpfact } from './arith1.js';
import { gen_product } from './bb_group.js';
import { PariError } from './errors.js';
type Scalar = bigint | MpReal;
const afffix = (x: MpReal, p: number) => rtor(x, Math.min(x.p, p));
function pairs(a: bigint, b: bigint, step: bigint): bigint[] {
  b -= (b - a) % step;
  const sum = a + b,
    v: bigint[] = [];
  for (let k = a; ; k += step) {
    const l = sum - k;
    if (l < k) break;
    v.push(l === k ? k : k * l);
    if (l === k) break;
  }
  return v;
}
function multiply(a: Scalar, b: Scalar): Scalar {
  if (typeof a === 'bigint') return typeof b === 'bigint' ? a * b : mulir(a, b);
  return typeof b === 'bigint' ? mulir(b, a) : mulrr(a, b);
}
function mulu_interval_step_prec(a: bigint, b: bigint, step: bigint, p: number): Scalar {
  return gen_product<Scalar>(pairs(a, b, step), (x, y) => {
    if (typeof x === 'bigint' && Math.ceil(x.toString(2).length / 64) * 64 > p) x = itor(x, p);
    if (typeof y === 'bigint' && Math.ceil(y.toString(2).length / 64) * 64 > p) y = itor(y, p);
    return multiply(x, y);
  });
}
function mpfactr_basecase(n: bigint, p: number): MpReal {
  const v: Scalar[] = [];
  for (let k = 1; ; k++) {
    const m = n >> BigInt(k - 1);
    if (m <= 2n) break;
    const l = (1n + (n >> BigInt(k))) | 1n;
    const a = mulu_interval_step_prec(l, m, 2n, p + 64);
    v.push(k === 1 ? a : typeof a === 'bigint' ? a ** BigInt(k) : powru(a, BigInt(k)));
  }
  let a = v.pop()!;
  while (v.length) a = multiply(a, v.pop()!);
  const real = typeof a === 'bigint' ? itor(a, p) : rtor(a, p);
  return shiftr(real, Number(n) - n.toString(2).replaceAll('0', '').length);
}
/** Original small-input path used by the Bernoulli cache. */
export function mpfactr_small(n: bigint, p: number): MpReal {
  return n < 410n ? itor(mpfact(n), p) : rtor(mpfactr_basecase(n, p), p);
}
function gamma_factorial(n: bigint, p = 64): MpReal<bigint> {
  const s = itor(n + 1n, p + 64),
    sd = rtodbl(s),
    u = (sd - 0.5) * Math.log(sd) - sd + Math.log(2 * Math.PI) / 2;
  const l = Math.max((p * Math.LN2 - Math.log(Math.max(u * u, 1e-6)) / 2) / 2, 0),
    lim = Math.max(Math.ceil(l / (1 + Math.log(4))), 1);
  const N = Math.max(Math.ceil(((lim - 0.5) * 4) / Math.PI - sd), 1);
  constbern(lim);
  let y = s;
  for (let i = 1; i < N; i++) y = mulir(n + 1n + BigInt(i), y);
  const nnx = addrs(s, N),
    a = invr(nnx),
    a2 = sqrr(a);
  const frac = (i: number) => {
    const [num, den] = bernfrac(2 * i);
    return [num, den * BigInt(2 * i * (2 * i - 1))] as const;
  };
  let [num, den] = frac(lim);
  // Generic rational x real: multiply numerator first, then divide denominator.
  let S: MpReal | undefined;
  for (let i = lim - 1; i >= 1; i--) {
    const product = S ? mulrr(a2, S) : divri(mulir(num, a2), den);
    [num, den] = frac(i);
    S = divri(addir(num, mulir(den, product)), den);
  }
  const correction = S ? mulrr(a, S) : divri(mulir(num, a), den);
  // Generic real - 1/2 also multiplies by the denominator before adding.
  const half = shiftr(addir(-1n, shiftr(nnx, 1)), -1);
  const B = addrr(subrr(mulrr(half, logr_abs(nnx)), nnx), correction);
  y = divrr(sqrtr(shiftr(mppi(p + 64), 1)), y);
  const result = mpexp(B);
  const value = afffix(mulrr({ ...result, e: 0 }, y), p);
  return { ...value, e: result.e + BigInt(value.e) };
}

/** PARI factorial real backend; precision is measured in bits, in 64-bit words.
 * @see Deviation: PARI factorial real representation and transcendental dependencies
 */
export function mpfactr(n: bigint, p = 64): MpReal<bigint> {
  if (!Number.isSafeInteger(p) || p < 64 || p % 64)
    throw new RangeError('PARI real precision must be a positive multiple of 64 bits');
  if (n < -(1n << 63n) || n >= 1n << 63n)
    throw new RangeError('mpfactr requires a signed-word argument');
  if (n < 410n) {
    const value = itor(mpfact(n), p);
    return { ...value, e: BigInt(value.e) };
  }
  const limit = p <= 64 ? 1930 : p <= 128 ? 2650 : p <= 192 ? 3300 : Math.trunc(p * Math.sqrt(p));
  if (n <= BigInt(limit)) {
    const value = rtor(mpfactr_basecase(n, p), p);
    return { ...value, e: BigInt(value.e) };
  }
  const result = gamma_factorial(n, p);
  if (result.e >= 1n << 61n || result.e < -(1n << 61n)) throw new PariError('overflow in expo()');
  return result;
}


/** Native atanh(u/v) by binary splitting, for unsigned words 0 < u < v. */
export function atanhuu(u: bigint, v: bigint, p: number): MpReal {
  const d = 2 * Math.log2(Number(v) / Number(u));
  const nmax = d ? Math.ceil(p / d) : -1;
  if (nmax < 0 || !Number.isFinite(nmax) || nmax >= 2 ** 63) throw new PariError('overflow in atanhuu');
  const A = abpq_init(nmax), u2 = u * u, v2 = v * v;
  A.a[0] = A.b[0] = 1n; A.p[0] = u; A.q[0] = v;
  for (let i = 1; i <= nmax; i++) {
    A.a[i] = 1n; A.b[i] = 2n * BigInt(i) + 1n; A.p[i] = u2; A.q[i] = v2;
  }
  const R = abpq_sum(0, nmax, A);
  return rdivii(R.T, R.B * R.Q, p);
}
