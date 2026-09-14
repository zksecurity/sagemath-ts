import type { zz_pXOptions } from './lzz_pX.js';
import { RandomBnd, type RandomStream } from './ZZ.js';
import {
  zz_pXModulus,
  MulMod,
  PowerMod,
  PowerXMod,
  PowerXPlusAMod,
  random as randomPolynomial,
  rem,
} from './lzz_pX.js';
import {
  BuildFromRoots,
  MinPolyMod,
  CompMod,
  Comp2Mod,
  Comp3Mod,
  zz_pXNewArgument,
  build as buildArgument,
  GCD,
  reduce,
} from './lzz_pX1.js';
/** NTL binary composition trace, preserving the native zero accumulator branches.
 * @see Deviation: NTL word polynomial composition adapters
 */
export function TraceMap(
  a: readonly bigint[],
  d: number,
  F: zz_pXModulus,
  b: readonly bigint[]
): bigint[] {
  if (d < 0) throw new Error('TraceMap: bad args');
  const k = F.arithmetic;
  let z = k.norm(b),
    y = k.norm(a),
    w: bigint[] = [],
    t: bigint[];
  while (d) {
    if (d === 1) {
      if (!w.length) w = y;
      else w = k.add(CompMod(w, z, F), y);
    } else if (d % 2 === 0) {
      [z, t] = Comp2Mod(z, y, z, F);
      y = k.add(t, y);
    } else if (!w.length) {
      w = y;
      [z, t] = Comp2Mod(z, y, z, F);
      y = k.add(t, y);
    } else {
      [z, t, w] = Comp3Mod(z, y, w, z, F);
      w = k.add(w, y);
      y = k.add(t, y);
    }
    d = Math.floor(d / 2);
  }
  return w;
}
/** NTL binary iterated composition; q=0 returns the unreduced polynomial X.
 * @see Deviation: NTL word polynomial composition adapters
 */
export function PowerCompose(h: readonly bigint[], q: number, F: zz_pXModulus): bigint[] {
  if (q < 0) throw new Error('PowerCompose: bad args');
  if (F.n < 0) throw new Error('negative length in vector::SetLength');
  let z = F.arithmetic.norm(h),
    y = [0n, 1n];
  while (q) {
    let sw = q > 1 ? 2 : 0;
    if (q % 2) {
      if (y.length === 2 && y[0] === 0n && y[1] === 1n) y = z;
      else sw |= 1;
    }
    if (sw === 1) y = CompMod(y, z, F);
    else if (sw === 2) z = CompMod(z, z, F);
    else if (sw === 3) [y, z] = Comp2Mod(y, z, z, F);
    q = Math.floor(q / 2);
  }
  return y;
}

type DegreeFactors = [bigint[], number][];
const integerSqrt = (n: number): number => {
  const a = BigInt(n);
  if (a < 2n) return n;
  let x = 1n << BigInt(Math.ceil(a.toString(2).length / 2));
  for (;;) {
    const y = (x + a / x) >> 1n;
    if (y >= x) return Number(x);
    x = y;
  }
};
/** NTL baby/giant-step distinct-degree factorization, preserving factor order.
 * The caller supplies a monic squarefree polynomial over a prime field and X^p mod f.
 * @see Deviation: NTL word polynomial distinct-degree factorization adapters
 * @see Deviation: NTL degree-dependent word contexts
 */
export function NewDDF(
  f: readonly bigint[],
  h: readonly bigint[],
  p: bigint,
  options?: zz_pXOptions
): DegreeFactors {
  const arithmetic = new zz_pXModulus(null, p, options).arithmetic;
  const ff = arithmetic.norm(f),
    hh = arithmetic.norm(h),
    n = ff.length - 1;
  if (ff[n] !== 1n) throw new Error('NewDDF: bad args');
  if (n === 0) return [];
  if (n === 1) return [[ff, 1]];
  const B = Math.floor(n / 2),
    rootn = integerSqrt(n),
    bits = p.toString(2).length;
  let k = integerSqrt(B);
  if (B >= 500) k = Math.floor(Math.min(2, Math.max(1, (0.25 * rootn) / bits)) * k);
  const babyFile: bigint[][] = [];
  let h1 = hh;
  const initial = new zz_pXModulus(ff, p, options);
  if (bits < Math.floor(rootn / 2)) {
    for (let i = 1; i < k; i++) {
      babyFile.push(h1);
      h1 = PowerMod(h1, p, initial);
    }
  } else {
    const H = new zz_pXNewArgument();
    buildArgument(H, hh, initial, 2 * rootn);
    for (let i = 1; i < k; i++) {
      babyFile.push(h1);
      h1 = CompMod(h1, H, initial);
    }
  }
  const giantArgument = new zz_pXNewArgument();
  buildArgument(giantArgument, h1, initial, 2 * rootn);
  let oldN = n;
  const giantFile = [h1];
  const fetchGiant = (gs: number, F: zz_pXModulus): bigint[] => {
    const l = giantFile.length;
    if (gs > l + 1) throw new Error('bad arg to FetchGiantStep');
    if (gs === l + 1) {
      let last = giantFile[l - 1]!;
      if (F.n < oldN) {
        last = rem(last, F);
        reduce(giantArgument, F);
        oldN = F.n;
      }
      const g = CompMod(last, giantArgument, F);
      giantFile.push(g);
      return g.slice();
    }
    const g = giantFile[gs - 1]!;
    return g.length - 1 >= F.n ? rem(g, F) : g.slice();
  };
  const fetchBaby = (): bigint[][] => [[0n, 1n], ...babyFile.map((x) => x.slice())];
  const processTable = (
    u: DegreeFactors,
    f: bigint[],
    F: zz_pXModulus,
    buf: bigint[][],
    size: number,
    start: number,
    intervalLength: number
  ): bigint[] => {
    if (size === 0) return f;
    let g = buf[size - 1]!;
    for (let i = 0; i < size - 1; i++) g = MulMod(g, buf[i]!, F);
    g = GCD(f, g, p, options);
    buf[size - 1] = g;
    if (g.length === 1) return f;
    f = arithmetic.divrem(f, g)[0];
    let d = (start - 1) * intervalLength + 1,
      i = 0,
      interval = start;
    while (i < size - 1 && 2 * d <= g.length - 1) {
      buf[i] = GCD(buf[i]!, g, p, options);
      if (buf[i]!.length > 1) {
        u.push([buf[i]!, interval]);
        g = arithmetic.divrem(g, buf[i]!)[0];
      }
      i++;
      interval++;
      d += intervalLength;
    }
    if (g.length > 1)
      u.push([g, i === size - 1 ? interval : Math.ceil((g.length - 1) / intervalLength)]);
    return f;
  };
  const u: DegreeFactors = [],
    baby = fetchBaby(),
    buf: bigint[][] = [];
  let current = ff,
    F = initial,
    g: bigint[] = [],
    size = 0,
    first = 0,
    d = 1;
  while (2 * d <= current.length - 1) {
    const oldDegree = current.length - 1,
      gs = Math.ceil(d / k),
      bs = gs * k - d;
    if (bs === k - 1) {
      size++;
      if (size === 1) first = gs;
      g = fetchGiant(gs, F);
      buf[size - 1] = arithmetic.sub(g, baby[bs]!);
    } else buf[size - 1] = MulMod(buf[size - 1]!, arithmetic.sub(g, baby[bs]!), F);
    if (size === 4 && bs === 0) {
      current = processTable(u, current, F, buf, size, first, k);
      size = 0;
    }
    d++;
    if (2 * d <= current.length - 1 && current.length - 1 < oldDegree) {
      F = new zz_pXModulus(current, p, options);
      for (let i = 1; i < k; i++) baby[i] = rem(baby[i]!, F);
    }
  }
  if (size > 0) current = processTable(u, current, F, buf, size, first, k);
  if (current.length > 1) u.push([current, 0]);
  const factors: DegreeFactors = [];
  let originalBaby: bigint[][] | undefined;
  for (const [part, gs] of u) {
    if (gs === 0 || 2 * ((gs - 1) * k + 1) > part.length - 1) {
      factors.push([part, part.length - 1]);
      continue;
    }
    originalBaby ??= fetchBaby();
    current = part;
    F = new zz_pXModulus(current, p, options);
    g = fetchGiant(gs, F);
    size = 0;
    d = (gs - 1) * k + 1;
    let bs = k - 1;
    while (bs >= 0 && 2 * d <= current.length - 1) {
      const oldDegree = current.length - 1;
      if (size === 0) first = d;
      buf[size] = arithmetic.sub(rem(originalBaby[bs]!, F), g);
      size++;
      if (size === 4) {
        current = processTable(factors, current, F, buf, size, first, 1);
        size = 0;
      }
      d++;
      bs--;
      if (bs >= 0 && 2 * d <= current.length - 1 && current.length - 1 < oldDegree) {
        F = new zz_pXModulus(current, p, options);
        g = rem(g, F);
      }
    }
    current = processTable(factors, current, F, buf, size, first, 1);
    if (current.length > 1) factors.push([current, current.length - 1]);
  }
  return factors;
}
/** NTL first squarefree Cantor–Zassenhaus stage: Frobenius and distinct-degree factors.
 * @see Deviation: NTL word polynomial distinct-degree factorization adapters
 * @see Deviation: NTL degree-dependent word contexts
 */
export function SFCanZass1(
  f: readonly bigint[],
  p: bigint,
  options?: zz_pXOptions
): [DegreeFactors, bigint[]] {
  const a = new zz_pXModulus(null, p, options).arithmetic.norm(f);
  if (a[a.length - 1] !== 1n || a.length === 1) throw new Error('SFCanZass1: bad args');
  const F = new zz_pXModulus(a, p, options),
    h = PowerXMod(p, F);
  return [NewDDF(a, h, p, options), h];
}

/** NTL recursive splitting of a monic polynomial with distinct prime-field roots. */
function RecFindRoots(
  f: bigint[],
  p: bigint,
  stream: RandomStream,
  out: bigint[],
  options?: zz_pXOptions
): void {
  const n = f.length - 1;
  if (n === 0) return;
  if (n === 1) {
    out.push((p - f[0]!) % p);
    return;
  }
  const F = new zz_pXModulus(f, p, options),
    k = F.arithmetic;
  let h: bigint[];
  do {
    const r = RandomBnd(p, stream, { word: true });
    h = GCD(k.sub(PowerXPlusAMod(r, p >> 1n, F), [1n]), f, p, options);
  } while (h.length <= 1 || h.length === f.length);
  RecFindRoots(h, p, stream, out, options);
  RecFindRoots(k.divrem(f, h)[0], p, stream, out, options);
}
/** NTL ordered roots of a monic polynomial that splits into distinct linear factors.
 * @see Deviation: NTL word factor recovery with explicit streams
 * @see Deviation: NTL degree-dependent word contexts
 */
export function FindRoots(
  f: readonly bigint[],
  p: bigint,
  stream: RandomStream,
  options?: zz_pXOptions
): bigint[] {
  const k = new zz_pXModulus(null, p, options).arithmetic,
    F = k.norm(f);
  if (!F.length) throw new Error('negative length in vector::SetLength');
  const out: bigint[] = [];
  RecFindRoots(F, p, stream, out, options);
  return out;
}
/** NTL single-root search, retaining the smaller successful split each iteration.
 * @see Deviation: NTL word factor recovery with explicit streams
 * @see Deviation: NTL degree-dependent word contexts
 */
export function FindRoot(
  f: readonly bigint[],
  p: bigint,
  stream: RandomStream,
  options?: zz_pXOptions
): bigint {
  const k = new zz_pXModulus(null, p, options).arithmetic;
  let Ff = k.norm(f);
  if (Ff.length <= 1 || Ff[Ff.length - 1] !== 1n) throw new Error('FindRoot: bad args');
  while (Ff.length > 2) {
    const F = new zz_pXModulus(Ff, p, options),
      r = RandomBnd(p, stream, { word: true });
    const h = GCD(k.sub(PowerXPlusAMod(r, p >> 1n, F), [1n]), Ff, p, options);
    if (h.length > 1 && h.length < Ff.length)
      Ff = h.length - 1 > Math.floor((Ff.length - 1) / 2) ? k.divrem(Ff, h)[0] : h;
  }
  return k.mod(-Ff[0]!);
}
/** NTL conversion of ordered roots into ordered monic linear factors.
 * @see Deviation: NTL word factor recovery with explicit streams
 * @see Deviation: NTL degree-dependent word contexts
 */
export function RootEDF(
  f: readonly bigint[],
  p: bigint,
  stream: RandomStream,
  options?: zz_pXOptions
): bigint[][] {
  return FindRoots(f, p, stream, options).map((r) => [(p - r) % p, 1n]);
}
/** NTL balanced recovery of factors from a separating element and its distinct values.
 * @see Deviation: NTL word factor recovery with explicit streams
 * @see Deviation: NTL degree-dependent word contexts
 */
export function FindFactors(
  f: readonly bigint[],
  g: readonly bigint[],
  roots: readonly bigint[],
  p: bigint,
  options?: zz_pXOptions
): bigint[][] {
  const k = new zz_pXModulus(null, p, options).arithmetic,
    values = roots.map(k.mod),
    out: bigint[][] = [];
  const rec = (f: bigint[], g: bigint[], lo: number, hi: number): void => {
    const r = hi - lo + 1;
    if (r === 0) return;
    if (r === 1) {
      out.push(f);
      return;
    }
    const mid = Math.floor((lo + hi) / 2),
      F = new zz_pXModulus(f, p, options),
      h = BuildFromRoots(values.slice(lo, mid + 1), p, options);
    const f1 = GCD(CompMod(h, g, F), f, p, options),
      f2 = k.divrem(f, f1)[0],
      g1 = k.divrem(g, f1)[1],
      g2 = k.divrem(g, f2)[1];
    rec(f1, g1, lo, mid);
    rec(f2, g2, mid + 1, hi);
  };
  rec(k.norm(f), k.norm(g), 0, roots.length - 1);
  return out;
}
/** NTL one-step equal-degree splitting through trace, minimum polynomial and roots.
 * @see Deviation: NTL word factor recovery with explicit streams
 * @see Deviation: NTL degree-dependent word contexts
 */
export function EDFSplit(
  f: readonly bigint[],
  b: readonly bigint[],
  d: number,
  p: bigint,
  stream: RandomStream,
  options?: zz_pXOptions
): bigint[][] {
  const F = new zz_pXModulus(f, p, options);
  if (d === 0) throw new Error('EDFSplit: degree must be nonzero');
  const r = Math.trunc(F.n / d),
    a = randomPolynomial(F.n, p, stream),
    g = TraceMap(a, d, F, b);
  const h = MinPolyMod(g, F, r, stream),
    roots = FindRoots(h, p, stream, options);
  return FindFactors(F.f, g, roots, p, options);
}
function RecEDF(
  f: bigint[],
  b: bigint[],
  d: number,
  p: bigint,
  stream: RandomStream,
  out: bigint[][],
  options?: zz_pXOptions
): void {
  const k = new zz_pXModulus(null, p, options).arithmetic;
  for (const v of EDFSplit(f, b, d, p, stream, options)) {
    if (v.length - 1 === d) out.push(v);
    else RecEDF(v, k.divrem(b, v)[1], d, p, stream, out, options);
  }
}
/** NTL equal-degree factorization in native output order.
 * @see Deviation: NTL word factor recovery with explicit streams
 * @see Deviation: NTL degree-dependent word contexts
 */
export function EDF(
  f: readonly bigint[],
  b: readonly bigint[],
  d: number,
  p: bigint,
  stream: RandomStream,
  options?: zz_pXOptions
): bigint[][] {
  const k = new zz_pXModulus(null, p, options).arithmetic,
    F = k.norm(f),
    B = k.norm(b);
  if (F[F.length - 1] !== 1n) throw new Error('EDF: bad args');
  if (d === 0) throw new Error('EDF: degree must be nonzero');
  const r = Math.trunc((F.length - 1) / d);
  if (r === 0) return [];
  if (r === 1) return [F];
  if (d === 1) return RootEDF(F, p, stream, options);
  const out: bigint[][] = [];
  RecEDF(F, B, d, p, stream, out, options);
  return out;
}
/** NTL completion of distinct-degree groups, preserving factor and random-draw order.
 * @see Deviation: NTL word factor recovery with explicit streams
 * @see Deviation: NTL degree-dependent word contexts
 */
export function SFCanZass2(
  groups: readonly (readonly [readonly bigint[], number])[],
  h: readonly bigint[],
  p: bigint,
  stream: RandomStream,
  options?: zz_pXOptions
): bigint[][] {
  const k = new zz_pXModulus(null, p, options).arithmetic,
    H = k.norm(h),
    out: bigint[][] = [];
  for (const [f, d] of groups) {
    const g = k.norm(f);
    if (d === 0) throw new Error('SFCanZass2: factor degree must be nonzero');
    const r = Math.trunc((g.length - 1) / d);
    if (r === 1) out.push(g);
    else if (d === 1) out.push(...RootEDF(g, p, stream, options));
    else out.push(...EDF(g, k.divrem(H, g)[1], d, p, stream, options));
  }
  return out;
}
/** NTL square-free monic polynomial factorization over a prime field.
 * @see Deviation: NTL word factor recovery with explicit streams
 * @see Deviation: NTL degree-dependent word contexts
 */
export function SFCanZass(
  f: readonly bigint[],
  p: bigint,
  stream: RandomStream,
  options?: zz_pXOptions
): bigint[][] {
  const k = new zz_pXModulus(null, p, options).arithmetic,
    poly = k.norm(f);
  if (poly[poly.length - 1] !== 1n) throw new Error('SFCanZass: bad args');
  if (poly.length === 1) return [];
  if (poly.length === 2) return [poly];
  const F = new zz_pXModulus(poly, p, options),
    h = PowerXMod(p, F),
    groups = NewDDF(poly, h, p, options);
  return SFCanZass2(groups, h, p, stream, options);
}
