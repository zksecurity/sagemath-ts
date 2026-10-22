/**
 * Short-Weierstrass arithmetic (y² = x³ + ax + b over F_p) on plain bigints, built for
 * bulk key searches rather than for generality.
 *
 * Single scalar multiplications run in Jacobian coordinates. Searches keep many affine
 * points side by side ("lanes") and step them all with one field inversion per round
 * (Montgomery's trick), which is what makes millions of additions affordable.
 */

export interface EcPoint {
  x: bigint;
  y: bigint;
}

interface Jacobian {
  X: bigint;
  Y: bigint;
  Z: bigint;
}

export interface EcGroupParams {
  p: bigint;
  a: bigint;
  b: bigint;
  n: bigint;
  g: EcPoint;
}

/** 53 low bits of a coordinate, as a float-exact number: the key of every hash table here. */
export function ecFingerprint(x: bigint): number {
  return Number(x & 0x1fffffffffffffn);
}

export class EcGroup {
  readonly p: bigint;
  readonly a: bigint;
  readonly b: bigint;
  readonly n: bigint;
  readonly g: EcPoint;
  private readonly aIsZero: boolean;

  constructor(params: EcGroupParams) {
    this.p = params.p;
    this.a = this.mod(params.a);
    this.b = this.mod(params.b);
    this.n = params.n;
    this.g = { x: this.mod(params.g.x), y: this.mod(params.g.y) };
    this.aIsZero = this.a === 0n;
  }

  mod(v: bigint): bigint {
    const r = v % this.p;
    return r < 0n ? r + this.p : r;
  }

  inv(v: bigint): bigint {
    const p = this.p;
    let r0 = p;
    let r1 = v % p;
    if (r1 < 0n) r1 += p;
    if (r1 === 0n) throw new RangeError('inverse of zero');
    let t0 = 0n;
    let t1 = 1n;
    while (r1 !== 0n) {
      const q = r0 / r1;
      const r2 = r0 - q * r1;
      r0 = r1;
      r1 = r2;
      const t2 = t0 - q * t1;
      t0 = t1;
      t1 = t2;
    }
    return t0 < 0n ? t0 + p : t0;
  }

  isOnCurve(P: EcPoint): boolean {
    const { p, a, b } = this;
    if (P.x < 0n || P.x >= p || P.y < 0n || P.y >= p) return false;
    return (P.y * P.y - ((P.x * P.x + a) % p) * P.x - b) % p === 0n;
  }

  equals(P: EcPoint | null, Q: EcPoint | null): boolean {
    if (P === null || Q === null) return P === Q;
    return P.x === Q.x && P.y === Q.y;
  }

  neg(P: EcPoint | null): EcPoint | null {
    return P === null ? null : { x: P.x, y: P.y === 0n ? 0n : this.p - P.y };
  }

  add(P: EcPoint | null, Q: EcPoint | null): EcPoint | null {
    if (P === null) return Q;
    if (Q === null) return P;
    const p = this.p;
    let l: bigint;
    if (P.x === Q.x) {
      if (P.y !== Q.y || P.y === 0n) return null;
      l = (((3n * P.x * P.x + this.a) % p) * this.inv(2n * P.y)) % p;
    } else {
      l = (this.mod(Q.y - P.y) * this.inv(Q.x - P.x)) % p;
    }
    const x = this.mod(l * l - P.x - Q.x);
    return { x, y: this.mod(l * (P.x - x) - P.y) };
  }

  sub(P: EcPoint | null, Q: EcPoint | null): EcPoint | null {
    return this.add(P, this.neg(Q));
  }

  /** k·P for any integer k; k is not reduced, so this also tests subgroup membership. */
  mul(P: EcPoint | null, k: bigint): EcPoint | null {
    return this.toAffine(this.mulJacobian(P, k));
  }

  /** 2^e · P. */
  shl(P: EcPoint | null, e: number): EcPoint | null {
    let J = this.toJacobian(P);
    for (let i = 0; i < e; i++) J = this.dbl(J);
    return this.toAffine(J);
  }

  /** Scalar multiples of many points, normalized with a single inversion. */
  mulMany(points: readonly (EcPoint | null)[], k: bigint | readonly bigint[]): (EcPoint | null)[] {
    return this.toAffineMany(
      points.map((P, i) => this.mulJacobian(P, typeof k === 'bigint' ? k : (k[i] as bigint)))
    );
  }

  /** [P, 2^step·P, 2^(2·step)·P, …] up to `count` terms, for every point. */
  shlChain(points: readonly (EcPoint | null)[], step: number, count: number): (EcPoint | null)[][] {
    const flat: Jacobian[] = [];
    for (const P of points) {
      let J = this.toJacobian(P);
      for (let c = 0; c < count; c++) {
        if (c > 0) for (let i = 0; i < step; i++) J = this.dbl(J);
        flat.push(J);
      }
    }
    const affine = this.toAffineMany(flat);
    return points.map((_, i) => affine.slice(i * count, (i + 1) * count));
  }

  private toJacobian(P: EcPoint | null): Jacobian {
    return P === null ? { X: 1n, Y: 1n, Z: 0n } : { X: P.x, Y: P.y, Z: 1n };
  }

  private toAffine(J: Jacobian): EcPoint | null {
    if (J.Z === 0n) return null;
    const p = this.p;
    const zi = this.inv(J.Z);
    const zi2 = (zi * zi) % p;
    return { x: (J.X * zi2) % p, y: (((J.Y * zi2) % p) * zi) % p };
  }

  private toAffineMany(js: readonly Jacobian[]): (EcPoint | null)[] {
    const p = this.p;
    const prefix: bigint[] = new Array(js.length);
    let acc = 1n;
    for (let i = 0; i < js.length; i++) {
      prefix[i] = acc;
      const Z = (js[i] as Jacobian).Z;
      if (Z !== 0n) acc = (acc * Z) % p;
    }
    let ia = acc === 1n ? 1n : this.inv(acc);
    const out: (EcPoint | null)[] = new Array(js.length);
    for (let i = js.length - 1; i >= 0; i--) {
      const J = js[i] as Jacobian;
      if (J.Z === 0n) {
        out[i] = null;
        continue;
      }
      const zi = (ia * (prefix[i] as bigint)) % p;
      ia = (ia * J.Z) % p;
      const zi2 = (zi * zi) % p;
      out[i] = { x: (J.X * zi2) % p, y: (((J.Y * zi2) % p) * zi) % p };
    }
    return out;
  }

  private dbl(J: Jacobian): Jacobian {
    const { p } = this;
    const { X, Y, Z } = J;
    if (Z === 0n || Y === 0n) return { X: 1n, Y: 1n, Z: 0n };
    const XX = (X * X) % p;
    const YY = (Y * Y) % p;
    const YYYY = (YY * YY) % p;
    const S = (4n * X * YY) % p;
    let M = 3n * XX;
    if (!this.aIsZero) {
      const ZZ = (Z * Z) % p;
      M += (this.a * ((ZZ * ZZ) % p)) % p;
    }
    M %= p;
    let X3 = (M * M - 2n * S) % p;
    if (X3 < 0n) X3 += p;
    let Y3 = (M * (S - X3) - 8n * YYYY) % p;
    if (Y3 < 0n) Y3 += p;
    return { X: X3, Y: Y3, Z: (2n * Y * Z) % p };
  }

  private addMixed(J: Jacobian, Q: EcPoint): Jacobian {
    const { p } = this;
    const { X, Y, Z } = J;
    if (Z === 0n) return { X: Q.x, Y: Q.y, Z: 1n };
    const ZZ = (Z * Z) % p;
    let H = ((Q.x * ZZ) % p) - X;
    if (H < 0n) H += p;
    let r = ((((Q.y * ZZ) % p) * Z) % p) - Y;
    if (r < 0n) r += p;
    if (H === 0n) return r === 0n ? this.dbl(J) : { X: 1n, Y: 1n, Z: 0n };
    const HH = (H * H) % p;
    const HHH = (H * HH) % p;
    const V = (X * HH) % p;
    let X3 = (r * r - HHH - 2n * V) % p;
    if (X3 < 0n) X3 += p;
    let Y3 = (r * (V - X3) - Y * HHH) % p;
    if (Y3 < 0n) Y3 += p;
    return { X: X3, Y: Y3, Z: (Z * H) % p };
  }

  private mulJacobian(P: EcPoint | null, k: bigint): Jacobian {
    let Q = P;
    let e = k;
    if (e < 0n) {
      e = -e;
      Q = this.neg(Q);
    }
    if (Q === null || e === 0n) return { X: 1n, Y: 1n, Z: 0n };
    const bits = e.toString(2);
    let J: Jacobian = { X: Q.x, Y: Q.y, Z: 1n };
    for (let i = 1; i < bits.length; i++) {
      J = this.dbl(J);
      if (bits.charCodeAt(i) === 49) J = this.addMixed(J, Q);
    }
    return J;
  }
}

/**
 * Affine points held as parallel arrays. `inf[i] = 1` marks the identity, whose
 * coordinates are meaningless.
 */
export class EcLanes {
  xs: bigint[];
  ys: bigint[];
  inf: Uint8Array;

  constructor(points: readonly (EcPoint | null)[]) {
    this.xs = points.map((P) => (P === null ? 0n : P.x));
    this.ys = points.map((P) => (P === null ? 0n : P.y));
    this.inf = Uint8Array.from(points, (P) => (P === null ? 1 : 0));
  }

  get length(): number {
    return this.xs.length;
  }

  point(i: number): EcPoint | null {
    return this.inf[i] ? null : { x: this.xs[i] as bigint, y: this.ys[i] as bigint };
  }

  set(i: number, P: EcPoint | null): void {
    if (P === null) {
      this.inf[i] = 1;
    } else {
      this.xs[i] = P.x;
      this.ys[i] = P.y;
      this.inf[i] = 0;
    }
  }

  /** Keep the lanes whose index passes `keep`, in order. */
  filter(keep: (i: number) => boolean): void {
    const xs: bigint[] = [];
    const ys: bigint[] = [];
    const inf: number[] = [];
    for (let i = 0; i < this.xs.length; i++) {
      if (!keep(i)) continue;
      xs.push(this.xs[i] as bigint);
      ys.push(this.ys[i] as bigint);
      inf.push(this.inf[i] as number);
    }
    this.xs = xs;
    this.ys = ys;
    this.inf = Uint8Array.from(inf);
  }

  /**
   * Lane i ← lane i + (qx[i], qy[i]) for every lane, sharing one inversion. The added
   * points must not be the identity. Lanes that hit a doubling or the identity take the
   * generic path, which costs an inversion each but practically never happens.
   */
  add(group: EcGroup, qx: readonly bigint[], qy: readonly bigint[]): void {
    const p = group.p;
    const { xs, ys, inf } = this;
    const count = xs.length;
    const prefix = scratch(count);
    let acc = 1n;
    for (let i = 0; i < count; i++) {
      prefix[i] = acc;
      if (inf[i]) continue;
      let d = (qx[i] as bigint) - (xs[i] as bigint);
      if (d < 0n) d += p;
      if (d !== 0n) acc = (acc * d) % p;
    }
    let ia = acc === 1n ? 1n : group.inv(acc);
    for (let i = count - 1; i >= 0; i--) {
      const x2 = qx[i] as bigint;
      const y2 = qy[i] as bigint;
      if (inf[i]) {
        xs[i] = x2;
        ys[i] = y2;
        inf[i] = 0;
        continue;
      }
      const x1 = xs[i] as bigint;
      const y1 = ys[i] as bigint;
      let d = x2 - x1;
      if (d < 0n) d += p;
      if (d === 0n) {
        this.set(i, group.add({ x: x1, y: y1 }, { x: x2, y: y2 }));
        continue;
      }
      const di = (ia * (prefix[i] as bigint)) % p;
      ia = (ia * d) % p;
      let dy = y2 - y1;
      if (dy < 0n) dy += p;
      const l = (dy * di) % p;
      let x3 = (l * l - x1 - x2) % p;
      if (x3 < 0n) x3 += p;
      let y3 = (l * (x1 - x3) - y1) % p;
      if (y3 < 0n) y3 += p;
      xs[i] = x3;
      ys[i] = y3;
    }
  }
}

/**
 * x(E_t + B_k) and x(E_t − B_k) for every lane t and block entry k, sharing one
 * inversion. Results land at index t·block.length + k; −1n marks the identity. Skipping
 * y makes each neighbour cost about 3.5 multiplications instead of 6.
 */
export function ecPlusMinusX(
  group: EcGroup,
  centers: EcLanes,
  bx: readonly bigint[],
  by: readonly bigint[],
  plus: bigint[],
  minus: bigint[]
): void {
  const p = group.p;
  const m = bx.length;
  const { xs, ys, inf } = centers;
  const total = xs.length * m;
  const prefix = scratch(total);
  let acc = 1n;
  for (let t = 0, i = 0; t < xs.length; t++) {
    const x0 = xs[t] as bigint;
    const lane = inf[t] === 1;
    for (let k = 0; k < m; k++, i++) {
      prefix[i] = acc;
      if (lane) continue;
      let d = (bx[k] as bigint) - x0;
      if (d < 0n) d += p;
      if (d !== 0n) acc = (acc * d) % p;
    }
  }
  let ia = acc === 1n ? 1n : group.inv(acc);
  for (let t = xs.length - 1, i = total - 1; t >= 0; t--) {
    const x0 = xs[t] as bigint;
    const y0 = ys[t] as bigint;
    const lane = inf[t] === 1;
    for (let k = m - 1; k >= 0; k--, i--) {
      const x1 = bx[k] as bigint;
      if (lane) {
        plus[i] = x1;
        minus[i] = x1;
        continue;
      }
      let d = x1 - x0;
      if (d < 0n) d += p;
      if (d === 0n) {
        const E = { x: x0, y: y0 };
        const B = { x: x1, y: by[k] as bigint };
        plus[i] = group.add(E, B)?.x ?? -1n;
        minus[i] = group.sub(E, B)?.x ?? -1n;
        continue;
      }
      const di = (ia * (prefix[i] as bigint)) % p;
      ia = (ia * d) % p;
      const y1 = by[k] as bigint;
      const s = x0 + x1;
      const lp = ((y1 - y0) * di) % p;
      let xp = (lp * lp - s) % p;
      if (xp < 0n) xp += p;
      const lm = ((y1 + y0) * di) % p;
      let xm = (lm * lm - s) % p;
      if (xm < 0n) xm += p;
      plus[i] = xp;
      minus[i] = xm;
    }
  }
}

let scratchBuffer: bigint[] = [];

function scratch(size: number): bigint[] {
  if (scratchBuffer.length < size) scratchBuffer = new Array(size).fill(0n);
  return scratchBuffer;
}
