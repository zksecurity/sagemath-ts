/*
fft62: a library for number-theoretic transforms

Copyright (C) 2013, David Harvey

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
import { ProbPrime, RandomBnd, type RandomStream } from './ZZ.js';

/** Mathematical fields of NTL's FFT-prime information; word preconditions are unnecessary for BigInt. */
export interface FFTPrimeInfo {
  readonly q: bigint;
  readonly qrecip: number;
  readonly RootTable: readonly [readonly bigint[], readonly bigint[]];
  readonly TwoInvTable: readonly bigint[];
}
interface FFTPrimeState {
  m: number;
  k: bigint;
  last_index: number;
  last_m: number;
  last_k: bigint;
  table: FFTPrimeInfo[];
}
const fftStates = new WeakMap<FFTPrimeContext, FFTPrimeState>();
/** Explicit counterpart of NTL's global FFT table and prime-search state.
 * @see Deviation: NTL FFT-prime cache contexts
 */
export class FFTPrimeContext {
  constructor() {
    fftStates.set(this, { m: 59, k: 0n, last_index: -1, last_m: 0, last_k: 0n, table: [] });
  }
  length(): number {
    return fftStates.get(this)!.table.length;
  }
  get(index: number): FFTPrimeInfo {
    const entry = fftStates.get(this)!.table[index];
    if (!entry) throw new RangeError('FFTPrimeContext: prime index is not initialized');
    return entry;
  }
}
function fftPower(a: bigint, e: bigint, p: bigint): bigint {
  let out = 1n;
  while (e) {
    if (e & 1n) out = (out * a) % p;
    a = (a * a) % p;
    e >>= 1n;
  }
  return out;
}
function fftInverse(a: bigint, p: bigint): bigint {
  let u = a,
    v = p,
    s = 1n,
    t = 0n;
  while (v) {
    const q = u / v;
    [u, v] = [v, u - q * v];
    [s, t] = [t, s - q * t];
  }
  if (u !== 1n) throw new Error('InvMod: inverse undefined');
  return ((s % p) + p) % p;
}
function fftSqrt(n: bigint): bigint {
  if (n < 2n) return n;
  let a = 1n << BigInt(Math.ceil(n.toString(2).length / 2));
  for (;;) {
    const b = (a + n / a) >> 1n;
    if (b >= a) return a;
    a = b;
  }
}
/** Native FFT-prime recognition; failed checks retain the previous root output.
 * @see Deviation: NTL FFT-prime cache contexts
 */
export function IsFFTPrime(n: bigint, stream: RandomStream, previousRoot = 0n): [number, bigint] {
  if (n < -(1n << 63n) || n >= 1n << 63n)
    throw new RangeError('IsFFTPrime: input must fit a signed 64-bit integer');
  if (n <= 1n || n >= 1n << 60n) return [0, previousRoot];
  for (const p of [2n, 3n, 5n, 7n]) if (n % p === 0n) return [0, previousRoot];
  let m = n - 1n,
    k = 0;
  while (!(m & 1n)) {
    m >>= 1n;
    k++;
  }
  let x = 0n;
  for (;;) {
    x = RandomBnd(n, stream, { word: true });
    if (!x) continue;
    let z = fftPower(x, m, n);
    if (z === 1n) continue;
    x = z;
    let y: bigint,
      j = 0;
    do {
      y = z;
      z = (y * y) % n;
      j++;
    } while (j !== k && z !== 1n);
    if (z !== 1n || y !== n - 1n) return [0, previousRoot];
    if (j === k) break;
  }
  let trialBound = m >> BigInt(k);
  if (trialBound > 0n) {
    if (!ProbPrime(n, stream, { NumTrials: 5, word: true })) return [0, previousRoot];
    trialBound = fftSqrt(trialBound);
    // Original special-factor trial loop, after the probabilistic prefilter.
    for (let a = 1n; a <= trialBound; a++) {
      const b = (a << BigInt(k)) + 1n;
      if (n % b === 0n) return [0, previousRoot];
    }
  }
  for (let j = 25; j < k; j++) x = (x * x) % n;
  return [1, x];
}
/** Native two-adic root limit for the 60-bit profile.
 * @see Deviation: NTL FFT-prime cache contexts
 */
export function CalcMaxRoot(p: bigint): number {
  if (p <= -(1n << 63n) || p >= 1n << 63n)
    throw new RangeError('CalcMaxRoot: p-1 must fit a signed 64-bit integer');
  if (p === 1n) throw new Error('CalcMaxRoot: p=1 has no terminating native result');
  let n = p - 1n,
    k = 0;
  while (!(n & 1n)) {
    n >>= 1n;
    k++;
  }
  return Math.min(k, 25);
}
/** Native descending candidate sequence, including same-index error rollback.
 * @see Deviation: NTL FFT-prime cache contexts
 */
export function NextFFTPrime(
  index: number,
  context: FFTPrimeContext,
  stream: RandomStream
): [bigint, bigint] {
  const s = fftStates.get(context)!;
  if (index === s.last_index) {
    s.m = s.last_m;
    s.k = s.last_k;
  } else {
    s.last_index = index;
    s.last_m = s.m;
    s.last_k = s.k;
  }
  for (;;) {
    if (s.k === 0n) {
      s.m--;
      if (s.m < 5) throw new Error('ran out of FFT primes');
      s.k = 1n << BigInt(60 - s.m - 2);
    }
    s.k--;
    const candidate = (1n << 59n) + (s.k << BigInt(s.m + 1)) + (1n << BigInt(s.m)) + 1n;
    const [ok, w] = IsFFTPrime(candidate, stream);
    if (ok) return [candidate, w];
  }
}
/** Native root and inverse-power tables, with immutable mathematical fields.
 * @see Deviation: NTL FFT-prime cache contexts
 */
export function InitFFTPrimeInfo(q: bigint, w: bigint): FFTPrimeInfo {
  if (q <= 1n || q >= 1n << 60n)
    throw new RangeError('InitFFTPrimeInfo: q must satisfy 1 < q < 2^60');
  if (w < 0n || w >= q) throw new RangeError('InitFFTPrimeInfo: root must satisfy 0 <= w < q');
  const mr = CalcMaxRoot(q),
    rt = Array<bigint>(mr + 1).fill(0n),
    rit = Array<bigint>(mr + 1).fill(0n),
    tit = Array<bigint>(mr + 1).fill(0n);
  rt[mr] = w;
  for (let j = mr - 1; j >= 0; j--) rt[j] = (rt[j + 1]! * rt[j + 1]!) % q;
  rit[mr] = fftInverse(w, q);
  for (let j = mr - 1; j >= 0; j--) rit[j] = (rit[j + 1]! * rit[j + 1]!) % q;
  const t = fftInverse(2n, q);
  tit[0] = 1n;
  for (let j = 1; j <= mr; j++) tit[j] = (tit[j - 1]! * t) % q;
  return Object.freeze({
    q,
    qrecip: 1 / Number(q),
    RootTable: Object.freeze([Object.freeze(rt), Object.freeze(rit)] as const),
    TwoInvTable: Object.freeze(tit),
  });
}
/** Native lazy cache extension. Existing entries consume no random bytes.
 * @see Deviation: NTL FFT-prime cache contexts
 */
export function UseFFTPrime(index: number, context: FFTPrimeContext, stream: RandomStream): void {
  if (index < 0) throw new Error('invalud FFT prime index');
  if (index >= 20000) throw new Error('FFT prime index too large');
  const s = fftStates.get(context)!;
  for (let i = s.table.length; i <= index; i++) {
    const [q, w] = NextFFTPrime(i, context, stream);
    s.table.push(InitFFTPrimeInfo(q, w));
  }
}
/** Native cached FFT modulus lookup.
 * @see Deviation: NTL FFT-prime cache contexts
 */
export function GetFFTPrime(index: number, context: FFTPrimeContext): bigint {
  return context.get(index).q;
}
/** Native floating reciprocal of a cached FFT modulus.
 * @see Deviation: NTL FFT-prime cache contexts
 */
export function GetFFTPrimeRecip(index: number, context: FFTPrimeContext): number {
  return context.get(index).qrecip;
}

import { FFTRoundUp } from './FFT_impl.js';
function transformInput(
  a: readonly bigint[],
  k: number,
  info: FFTPrimeInfo,
  yn: number,
  xn: number
): bigint[] {
  if (!Number.isInteger(k) || k < 0 || k > CalcMaxRoot(info.q))
    throw new RangeError('FFT transform: exponent exceeds the prime root table');
  const n = 2 ** k;
  for (const length of [yn, xn])
    if (!Number.isInteger(length) || length < 1 || length > n || FFTRoundUp(length, k) !== length)
      throw new RangeError('FFT transform: lengths must be admissible');
  if (a.length < xn)
    throw new RangeError('FFT transform: input is shorter than its declared length');
  return Array.from({ length: n }, (_, i) => (i < xn ? ((a[i]! % info.q) + info.q) % info.q : 0n));
}
// NTL's new_fft_short / new_ifft_short1 / new_ifft_short2 divide-and-conquer
// recurrences. BigInt modular butterflies replace lazy machine-word reductions.
const transformPowers = new WeakMap<FFTPrimeInfo, [bigint[][], bigint[][]]>();
function transformKernels(info: FFTPrimeInfo, k: number) {
  const q = info.q,
    mod = (x: bigint) => {
      const r = x % q;
      return r < 0n ? r + q : r;
    };
  let tables = transformPowers.get(info);
  if (!tables) {
    tables = [[[]], [[]]];
    transformPowers.set(info, tables);
  }
  for (let direction = 0; direction < 2; direction++) {
    const table = tables[direction]!;
    while (table.length <= k) {
      const level = table.length,
        a = Array<bigint>(2 ** (level - 1));
      let w = 1n;
      for (let j = 0; j < a.length; j++) {
        a[j] = w;
        w = (w * info.RootTable[direction]![level]!) % q;
      }
      table.push(a);
    }
  }
  const [forward, inverse] = tables;
  const fft = (a: bigint[], offset: number, yn: number, xn: number, level: number): void => {
    if (!level) return;
    const half = 2 ** (level - 1),
      w = forward[level]!;
    if (yn <= half) {
      if (xn > half)
        for (let j = 0; j < xn - half; j++)
          a[offset + j] = mod(a[offset + j]! + a[offset + half + j]!);
      fft(a, offset, yn, Math.min(xn, half), level - 1);
    } else {
      for (let j = 0; j < half; j++) {
        const x = j < xn ? a[offset + j]! : 0n,
          y = j + half < xn ? a[offset + half + j]! : 0n;
        a[offset + j] = mod(x + y);
        a[offset + half + j] = mod((x - y) * w[j]!);
      }
      fft(a, offset, half, Math.min(xn, half), level - 1);
      fft(a, offset + half, yn - half, Math.min(xn, half), level - 1);
    }
  };
  const ifft = (
    a: bigint[],
    offset: number,
    yn: number,
    level: number,
    fullTail: boolean
  ): void => {
    if (!level) return;
    const half = 2 ** (level - 1),
      w = forward[level]!,
      iw = inverse[level]!;
    if (yn <= half) {
      for (let j = 0; j < yn; j++) a[offset + j] = mod(2n * a[offset + j]!);
      if (fullTail)
        for (let j = yn; j < half; j++) a[offset + j] = mod(a[offset + j]! + a[offset + half + j]!);
      ifft(a, offset, yn, level - 1, fullTail);
      if (fullTail)
        for (let j = 0; j < yn; j++) a[offset + j] = mod(a[offset + j]! - a[offset + half + j]!);
    } else {
      ifft(a, offset, half, level - 1, false);
      const rest = yn - half;
      for (let j = rest; j < half; j++) {
        const x = a[offset + j]!,
          y = fullTail ? a[offset + half + j]! : 0n;
        a[offset + j] = mod(2n * x - y);
        a[offset + half + j] = mod((x - y) * w[j]!);
      }
      ifft(a, offset + half, rest, level - 1, true);
      for (let j = 0; j < rest; j++) {
        const x = a[offset + j]!,
          y = mod(a[offset + half + j]! * iw[j]!);
        a[offset + j] = mod(x + y);
        a[offset + half + j] = mod(x - y);
      }
    }
  };
  return { fft, ifft };
}
/** Native truncated forward transform in bit-reversed frequency order.
 * @see Deviation: NTL portable truncated transforms
 */
export function FFTFwd_trunc(
  a: readonly bigint[],
  k: number,
  info: FFTPrimeInfo,
  yn: number,
  xn: number
): bigint[] {
  const out = transformInput(a, k, info, yn, xn);
  transformKernels(info, k).fft(out, 0, yn, xn, k);
  return out.slice(0, yn);
}
/** Native normalized inverse truncated transform.
 * @see Deviation: NTL portable truncated transforms
 */
export function FFTRev1_trunc(
  a: readonly bigint[],
  k: number,
  info: FFTPrimeInfo,
  yn: number
): bigint[] {
  const out = transformInput(a, k, info, yn, yn);
  transformKernels(info, k).ifft(out, 0, yn, k, false);
  return out.slice(0, yn).map((x) => (x * info.TwoInvTable[k]!) % info.q);
}
const flippedPrimeInfo = new WeakMap<FFTPrimeInfo, FFTPrimeInfo>();
function flippedInfo(info: FFTPrimeInfo): FFTPrimeInfo {
  let flipped = flippedPrimeInfo.get(info);
  if (!flipped) {
    flipped = { ...info, RootTable: [info.RootTable[1], info.RootTable[0]] };
    flippedPrimeInfo.set(info, flipped);
  }
  return flipped;
}
/** Transpose of the native forward FFT.
 * @see Deviation: NTL portable transposed transforms
 */
export function FFTFwd_trans(a: readonly bigint[], k: number, info: FFTPrimeInfo): bigint[] {
  const n = 2 ** k,
    out = transformInput(a, k, info, n, n);
  // new_ifft_short1_flipped uses forward roots and omits inverse scaling.
  transformKernels(flippedInfo(info), k).ifft(out, 0, n, k, false);
  return out;
}
/** Transpose of the native normalized inverse FFT.
 * @see Deviation: NTL portable transposed transforms
 */
export function FFTRev1_trans(a: readonly bigint[], k: number, info: FFTPrimeInfo): bigint[] {
  const n = 2 ** k,
    out = transformInput(a, k, info, n, n);
  // new_fft_short_flipped uses inverted roots and scales every output by 1/N.
  transformKernels(flippedInfo(info), k).fft(out, 0, n, n, k);
  return out.map((x) => (x * info.TwoInvTable[k]!) % info.q);
}
