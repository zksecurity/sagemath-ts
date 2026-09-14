/** Batched ECM arithmetic and continuation from PARI ifactor1.c:376–1031.
 * Cell addresses retain native affii mutation and GEN pointer aliasing.
 */
import { forprime, isPrime, isqrt } from './ifactor.js';
import { Fp_pow, gcd, xgcd } from './ff.js';
const DELTA = [
  10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 6, 2, 4,
  6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2,
];
const RESIDUE = [1];
for (let i = 0; i < 47; i++) RESIDUE.push(RESIDUE[i] + DELTA[i]);
// Match the initialized native PARI prime table used by the other factor ports.
let primeTable: number[] | undefined;
function base2(n: number): boolean {
  const N = BigInt(n);
  let d = N - 1n,
    s = 0;
  while (!(d & 1n)) {
    d >>= 1n;
    s++;
  }
  let a = Fp_pow(2n, d, N);
  if (a === 1n || a === N - 1n) return true;
  for (let i = 1; i < s; i++) {
    a = (a * a) % N;
    if (a === N - 1n) return true;
    if (a === 1n) return false;
  }
  return false;
}
function get_rule(d: number, e: number): number {
  if (d <= e + Math.floor(e / 4)) {
    if ((d + e) % 3 === 0) return 0;
    if ((d - e) % 6 === 0) return 1;
  }
  if (Math.floor((d + 3) / 4) <= e) return 2;
  if (d % 2 === e % 2) return 1;
  if (d % 2 === 0) return 3;
  if (d % 3 === 0) return 4;
  if ((d + e) % 3 === 0) return 5;
  if ((d - e) % 3 === 0) return 6;
  return 7;
}
export class ECM {
  readonly nbc2: number;
  readonly X = 0;
  readonly XAUX: number;
  readonly XT: number;
  readonly XD: number;
  readonly XB: number;
  readonly XB2: number;
  readonly XH: number;
  readonly Xh: number;
  readonly Yh: number;
  private readonly addresses: number[];
  private readonly values: Array<bigint | undefined>;
  g = 1n;
  seed: bigint;
  constructor(
    readonly N: bigint,
    readonly nbc: number,
    seed: number
  ) {
    this.seed = BigInt(seed);
    const n2 = (this.nbc2 = 2 * nbc);
    this.XAUX = n2;
    this.XT = 2 * n2;
    this.XD = 3 * n2;
    this.XB = 13 * n2;
    this.XB2 = this.XB + 2048;
    this.XH = this.XB + 4096;
    // These are the original pointer offsets, including Yh within XH.
    this.Xh = this.XH + 48 * n2;
    this.Yh = this.XH + 192;
    const size = 61 * n2 + 4096 + 385;
    this.addresses = Array.from({ length: size }, (_, i) => i);
    this.values = new Array(size);
  }
  get(i: number): bigint {
    const v = this.values[this.addresses[i]];
    if (v === undefined) throw new Error('uninitialized ECM cell ' + i);
    return v;
  }
  set(i: number, v: bigint): void {
    this.values[this.addresses[i]] = v;
  }
  alias(dest: number, src: number): void {
    this.addresses[dest] = this.addresses[src];
  }
  copy(n: number, src: number, dest: number): void {
    if (src !== dest) for (let i = n; i--; ) this.set(dest + i, this.get(src + i));
  }
  mod(x: bigint): bigint {
    const r = x % this.N;
    return r < 0n ? r + this.N : r;
  }
  inverse(x: bigint): boolean {
    const [d, s] = xgcd(this.mod(x), this.N);
    this.g = d === 1n ? this.mod(s) : d;
    return d === 1n;
  }
  addPoint(
    z: bigint,
    Px: bigint,
    Py: bigint,
    Qx: bigint,
    Qy: bigint,
    Rx: number,
    Ry: number | null
  ): void {
    const slope = this.mod((Py - Qy) * z),
      x = this.mod(slope * slope - Qx - Px);
    this.set(Rx, x);
    if (Ry !== null) this.set(Ry, this.mod(slope * (Px - x) - Py));
  }
  add0(
    n: number,
    n1: number,
    X1: number,
    Y1: number,
    X2: number,
    Y2: number,
    X3: number,
    Y3: number | null
  ): number {
    const W: bigint[] = [],
      A: bigint[] = [];
    W[1] = this.get(X1) - this.get(X2);
    for (let i = 1; i < n; i++) {
      A[i] = this.get(X1 + (n1 === 4 ? i % 4 : i)) - this.get(X2 + i);
      W[i + 1] = this.mod(A[i] * W[i]);
    }
    if (!this.inverse(W[n])) {
      if (this.g !== this.N) return 2;
      this.copy(n, X2, X3);
      if (Y3 !== null) this.copy(n, Y2, Y3);
      return 1;
    }
    for (let i = n; i--; ) {
      const j = n1 === 4 ? i % 4 : i,
        z = i ? this.g * W[i] : this.g;
      this.addPoint(
        z,
        this.get(X1 + j),
        this.get(Y1 + j),
        this.get(X2 + i),
        this.get(Y2 + i),
        X3 + i,
        Y3 === null ? null : Y3 + i
      );
      if (i) this.g = this.mod(this.g * A[i]);
    }
    return 0;
  }
  add(n: number, X1: number, X2: number, X3: number): number {
    return this.add0(n, n, X1, X1 + n, X2, X2 + n, X3, X3 + n);
  }
  add2(n: number, X1: number, X2: number, X3: number, X4: number, X5: number, X6: number): number {
    const W: bigint[] = [],
      A: bigint[] = [];
    W[1] = this.get(X1) - this.get(X2);
    for (let i = 1; i < 2 * n; i++) {
      A[i] =
        i < n ? this.get(X1 + i) - this.get(X2 + i) : this.get(X4 + i - n) - this.get(X5 + i - n);
      W[i + 1] = this.mod(A[i] * W[i]);
    }
    if (!this.inverse(W[2 * n])) {
      if (this.g !== this.N) return 2;
      this.copy(2 * n, X2, X3);
      this.copy(2 * n, X5, X6);
      return 1;
    }
    for (let i = 2 * n; i--; ) {
      const j = i < n ? i : i - n,
        Pa = i < n ? X1 : X4,
        Qa = i < n ? X2 : X5,
        Ra = i < n ? X3 : X6;
      this.addPoint(
        i ? this.g * W[i] : this.g,
        this.get(Pa + j),
        this.get(Pa + n + j),
        this.get(Qa + j),
        this.get(Qa + n + j),
        Ra + j,
        Ra + n + j
      );
      if (i) this.g = this.mod(this.g * A[i]);
    }
    return 0;
  }
  double(n: number, X1: number, X2: number): number {
    const W: bigint[] = [];
    W[1] = this.get(X1 + n);
    for (let i = 1; i < n; i++) W[i + 1] = this.mod(this.get(X1 + n + i) * W[i]);
    if (!this.inverse(W[n])) {
      if (this.g !== this.N) return 2;
      this.copy(2 * n, X1, X2);
      return 1;
    }
    for (let i = n; i--; ) {
      const x = this.get(X1 + i),
        y = this.get(X1 + n + i),
        z = i ? this.g * W[i] : this.g;
      if (i) this.g = this.mod(this.g * y);
      let L = this.mod((1n + 3n * this.mod(x * x)) * z);
      if (L) L = (L % 2n ? L + this.N : L) / 2n;
      const v = this.mod(L * L - 2n * x),
        w = this.mod(L * (x - v) - y);
      this.set(X2 + i, v);
      this.set(X2 + n + i, w);
    }
    return 0;
  }
  mult(n: number, k: number, X1: number, X2: number, XAUX: number): number {
    let A = X2,
      B = XAUX;
    const T = XAUX + 2 * n;
    let res: number;
    this.copy(2 * n, X1, XAUX);
    if ((res = this.double(n, X1, X2))) return res;
    // Native PRAC heuristic: preserve its binary64 golden-ratio split.
    const r = Math.floor(k * 0.61803398875 + 0.5);
    let d = k - r,
      e = r - d;
    while (d !== e) {
      switch (get_rule(d, e)) {
        case 0: {
          if ((res = this.add(n, A, B, T))) return res;
          if ((res = this.add2(n, T, A, A, T, B, B))) return res;
          const e1 = d - e;
          d = (d + e1) / 3;
          e = (e - e1) / 3;
          break;
        }
        case 1:
          if ((res = this.add(n, A, B, B))) return res;
          if ((res = this.double(n, A, A))) return res;
          d = (d - e) / 2;
          break;
        case 3:
          if ((res = this.double(n, A, A))) return res;
          d /= 2;
          break;
        case 4:
          if ((res = this.double(n, A, T))) return res;
          if ((res = this.add(n, T, A, A))) return res;
          if ((res = this.add(n, A, B, B))) return res;
          d = d / 3 - e;
          break;
        case 2:
          if ((res = this.add(n, A, B, B))) return res;
          d -= e;
          break;
        case 5:
          if ((res = this.double(n, A, T))) return res;
          if ((res = this.add2(n, T, A, A, T, B, B))) return res;
          d = (d - 2 * e) / 3;
          break;
        case 6:
          if ((res = this.add(n, A, B, B))) return res;
          if ((res = this.double(n, A, T))) return res;
          if ((res = this.add(n, T, A, A))) return res;
          d = (d - e) / 3;
          break;
        case 7:
          if ((res = this.double(n, B, B))) return res;
          e /= 2;
          break;
      }
      if (d < e) {
        [d, e] = [e, d];
        [A, B] = [B, A];
      }
    }
    return this.add(n, XAUX, X2, X2);
  }
  round(B1: number): bigint | null {
    const B2 = 110 * B1,
      rt = Number(isqrt(BigInt(B2))),
      n = this.nbc,
      n2 = this.nbc2;
    const { X, XAUX, XT, XD, XB, XB2, XH, Xh, Yh } = this;
    for (let i = n2; i--; ) this.set(X + i, BigInt(this.seed++));
    const gse =
      B1 < 656 ? (B1 < 200 ? 5 : 6) : B1 < 10500 ? (B1 < 2625 ? 7 : 8) : B1 < 42000 ? 9 : 10;
    const gss = 2 ** gse,
      XG = XT + gse * n2,
      YG = XG + n;
    let p = 2,
      index = 0,
      rcn = 128,
      bstp = 0;
    primeTable ??= Array.from(forprime(2, 500000));
    const next = (probable = false): number => {
      const old = p;
      if (index < primeTable!.length - 1) p = primeTable![++index];
      else {
        let rc = RESIDUE.indexOf(p % 210);
        do {
          p += DELTA[rc];
          rc = (rc + 1) % 48;
        } while (!(probable ? base2(p) : isPrime(BigInt(p))));
      }
      if (rcn !== 128) {
        bstp += Math.floor(p / 210) - Math.floor(old / 210);
        rcn = RESIDUE.indexOf(p % 210);
      }
      return p;
    };
    for (let m = 1; m <= Math.floor(B2 / 2); m *= 2) {
      const fl = this.double(n, X, X);
      if (fl > 1) return this.g;
      if (fl) break;
    }
    while (p < B1 && p <= rt) {
      next();
      for (let m = 1; m <= Math.floor(B2 / p); m *= p) {
        const fl = this.mult(n, p, X, X, XAUX);
        if (fl > 1) return this.g;
        if (fl) break;
      }
    }
    while (p < B1) {
      next();
      if (this.mult(n, p, X, X, XAUX) > 1) return this.g;
    }
    if (this.double(n, X, XD) > 1) return this.g;
    if (this.double(n, XD, XD + n2) > 1) return this.g;
    if (this.add(n, XD, XD + n2, XD + 4 * n) > 1) return this.g;
    if (this.add2(n, XD, XD + 4 * n, XT + 8 * n, XD + n2, XD + 4 * n, XD + 8 * n) > 1)
      return this.g;
    next();
    if (rcn === 128) rcn = RESIDUE.indexOf(p % 210);
    if (this.mult(n, p, X, XH + rcn * n2, XAUX) > 1) return this.g;
    const p0 = p,
      index0 = index,
      rcn0 = rcn;
    for (let i = 47; i; i--) {
      const dp = DELTA[rcn];
      p += dp;
      const target = rcn === 47 ? XH : XH + (rcn + 1) * n2;
      if (this.add(n, XT + dp * n, XH + rcn * n2, target) > 1) return this.g;
      rcn = (rcn + 1) % 48;
    }
    if (this.mult(n, 3, XD + 8 * n, X, XAUX) > 1) return this.g;
    if (this.mult(n, 7, X, X, XAUX) > 1) return this.g;
    if (this.double(n, X, XAUX) > 1) return this.g;
    if (this.add(n, X, XAUX, XT) > 1) return this.g;
    if (this.add(n, X, XT, XD) > 1) return this.g;
    for (let i = 1; i <= gse; i++) if (this.double(n, XT + i * n2, XD + i * n2) > 1) return this.g;
    for (let i = n - 4; i >= 0; i -= 4) {
      for (let j = 48; j--; ) {
        const k = n2 * j + i,
          m = 4 * j;
        for (let t = 0; t < 4; t++) this.alias(Xh + m + t, XH + k + t);
        for (let t = 0; t < 4; t++) this.alias(Yh + m + t, XH + k + n + t);
      }
      for (let k = 1; k >= 0; k--) {
        const j = i + (k === 0 ? n : 0);
        let Xb = k ? XB : XB2;
        for (let t = 0; t < 4; t++) {
          this.alias(Xb + t, X + j + t);
          this.alias(Xb + 4 + t, XAUX + j + t);
          this.alias(Xb + 8 + t, XT + j + t);
        }
        Xb += 4;
        for (let m = 2; m < gse + k; m++) {
          Xb += 2 ** (m + 1);
          for (let t = 0; t < 4; t++) this.alias(Xb + t, XAUX + m * n2 + j + t);
        }
      }
      if (this.add0(12, 4, XB + 12, XB2 + 12, XB, XB2, XB + 16, XB2 + 16) > 1) return this.g;
      if (this.add0(28, 4, XB + 28, XB2 + 28, XB, XB2, XB + 32, XB2 + 32) > 1) return this.g;
      for (let m = 5; m <= gse; m++) {
        const m2 = 2 ** (m + 1);
        let j = 0;
        for (; j < m2 - 64; j += 64)
          if (
            this.add0(
              64,
              4,
              XB + m2 - 4,
              XB2 + m2 - 4,
              XB + j,
              XB2 + j,
              XB + m2 + j,
              m < gse ? XB2 + m2 + j : null
            ) > 1
          )
            return this.g;
        if (
          this.add0(
            60,
            4,
            XB + m2 - 4,
            XB2 + m2 - 4,
            XB + j,
            XB2 + j,
            XB + m2 + j,
            m < gse ? XB2 + m2 + j : null
          ) > 1
        )
          return this.g;
      }
      bstp = 0;
      p = p0;
      index = index0;
      rcn = rcn0;
      this.g = 1n;
      while (p < B2) {
        next(true);
        let k = bstp - (rcn < rcn0 ? 1 : 0);
        if (k > gss) {
          this.g = gcd(this.g, this.N);
          if (this.g !== 1n && this.g !== this.N) return this.g;
          this.g = 1n;
          while (k > gss) {
            for (let j = 0; j < 192; j += 64)
              if (this.add0(64, 4, XG + i, YG + i, Xh + j, Yh + j, Xh + j, Yh + j) > 1)
                return this.g;
            bstp -= 2 * gss;
            k = bstp - (rcn < rcn0 ? 1 : 0);
          }
        }
        if (!k) continue;
        const m = 4 * (Math.abs(k) - 1),
          j = 4 * rcn;
        for (let t = 0; t < 4; t++)
          this.g = this.mod(this.g * (this.get(XB + m + t) - this.get(Xh + j + t)));
      }
    }
    return null;
  }
}
export function ECM_loop(N: bigint, nbc: number, seed: number, B1: number): bigint | null {
  return new ECM(N, nbc, seed).round(B1);
}
