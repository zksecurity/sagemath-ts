/** Native Bernoulli table and zeta-based cache, PARI basemath/bern.c.
 * The supplied B2..B26 constants are available without eagerly filling unused cache entries.
 */
import {
  itor,
  rtor,
  mulrr,
  mulir,
  divrr,
  divri,
  shiftr,
  sqrr,
  addrr,
  addrs,
  type MpReal,
} from './qfb.js';
import { mppi, powru } from './trans1.js';
import { mpfactr_small as smallFactorial } from './trans2.js';
import { Z_factor, isPrime } from './ifactor.js';
import { rdivii } from './kernel/none/level1.js';
const nbits = (p: number) => Math.ceil(p / 64) * 64;
const truncate = (x: MpReal, p: number): MpReal => ({ ...x, p, m: x.m >> BigInt(x.p - p) });
const bern = [
  [1n, 6n],
  [-1n, 30n],
  [1n, 42n],
  [-1n, 30n],
  [5n, 66n],
  [-691n, 2730n],
  [7n, 6n],
  [-3617n, 510n],
  [43867n, 798n],
  [-174611n, 330n],
  [854513n, 138n],
  [-236364091n, 2730n],
  [8553103n, 6n],
];
function fracB2k(n: number): [bigint, bigint] {
  let divisors = [1n];
  for (const [p, e] of Z_factor(BigInt(n))) {
    const prior = divisors.slice();
    let power = 1n;
    for (let k = 0n; k < e; k++) {
      power *= p;
      for (const d of prior) divisors.push(d * power);
    }
  }
  divisors.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  let a = 5n,
    b = 6n;
  for (const d of divisors.slice(1)) {
    const p = 2n * d + 1n;
    if (isPrime(p)) {
      a = a * p + b;
      b *= p;
    }
  }
  return [a, b];
}
function roundr(x: MpReal): bigint {
  if (!x.s) return 0n;
  const k = x.p - 1 - x.e;
  if (k <= 0) return BigInt(x.s) * (x.m << BigInt(-k));
  const denominator = 1n << BigInt(k),
    num = BigInt(x.s) * x.m + denominator / 2n;
  return num >= 0 ? num / denominator : -((-num + denominator - 1n) / denominator);
}
function bernbitprec(n: number): number {
  return Math.ceil(((n + 4) * Math.log(n) - n * (1 + 1.83787706641) + 1.612086) / Math.LN2) + 10;
}
export function constbern(requested: number): void {
  if (bern.length >= requested) return;
  const m = bern.length + 1,
    n = Math.max(requested, bern.length === 13 ? 127 : bern.length + 128),
    N = 2 * n;
  let p = bernbitprec(N),
    prec = nbits(p),
    r = Math.ceil(N / (2 * Math.PI * Math.E)) | 1;
  const u = sqrr(shiftr(mppi(prec), 1));
  let v = shiftr(divrr(smallFactorial(BigInt(N), prec), powru(u, BigInt(n))), 1);
  const t = new Map<number, bigint>();
  for (let j = 3; j <= r; j += 2) t.set(j, (1n << BigInt(p)) / BigInt(j) ** BigInt(N));
  for (let i = n, k = N; ; i--) {
    let sum = t.get(r) ?? 0n;
    for (let j = r - 2; j >= 3; j -= 2) sum += t.get(j)!;
    const z = shiftr(itor(sum, nbits(p)), -p);
    let tail = shiftr(addrs(z, 1), -k);
    for (let h = k; h < p; h *= 2) tail = addrr(tail, shiftr(tail, -h));
    let B = addrr(v, mulrr(v, addrr(z, tail)));
    if (i % 2 === 0) B = { ...B, s: -1 };
    const [a, b] = fracB2k(i),
      rational = rdivii(a, b, 64);
    const numerator = roundr(addrr(B, rational)) * b - a;
    bern[i - 1] = [numerator, b];
    if (i === m) break;
    v = rtor(divri(mulrr(v, u), BigInt((k - 1) * k)), prec);
    for (let j = r; j >= 3; j -= 2) t.set(j, t.get(j)! * BigInt(j * j));
    k -= 2;
    if (((N - k) & 127) === 126) {
      const p2 = p,
        prec2 = prec;
      p = bernbitprec(k);
      prec = nbits(p);
      if (prec2 !== prec) {
        v = truncate(v, prec);
        r = Math.ceil(k / (2 * Math.PI * Math.E)) | 1;
        for (let j = 3; j <= r; j += 2) t.set(j, t.get(j)! >> BigInt(p2 - p));
      }
    }
  }
}
export function bernfrac(n: number): [bigint, bigint] {
  if (!Number.isSafeInteger(n) || n < 0)
    throw new RangeError('Bernoulli index must be a nonnegative safe integer');
  if (n === 0) return [1n, 1n];
  if (n === 1) return [-1n, 2n];
  if (n % 2) return [0n, 1n];
  constbern(n / 2);
  return [...bern[n / 2 - 1]!] as [bigint, bigint];
}
