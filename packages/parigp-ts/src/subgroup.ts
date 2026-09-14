import { absZ_factor_limit_strict } from './ifactor1.js';
import { logint0 } from './ispower.js';
/** PARI subgroup.c: Birkhoff enumeration, exact/scalar/unbounded modes and native order. */
import { PariError } from './errors.js';
import { Z_factor } from './ifactor.js';
import { gcd } from './ff.js';
import { ZM_hnfmodid } from './hnf_snf.js';
const val = (a: bigint, p: bigint) => {
  let v = 0;
  while (a !== 0n && a % p === 0n) {
    a /= p;
    v++;
  }
  return v;
};
const eye = (n: number): bigint[][] =>
  Array.from({ length: n }, (_, j) => Array.from({ length: n }, (_, i) => (i === j ? 1n : 0n)));
function boundedSubgroups(CYC: bigint[], bound: bigint | null, exact: boolean): bigint[][][] {
  const N = CYC.length;
  let cyc = CYC.slice();
  if (bound !== null) {
    if (exact) cyc = cyc.map((c) => gcd(c, bound));
    else if (bound < 1n << 63n) {
      const [f] = absZ_factor_limit_strict(cyc[0] ?? 1n, bound + 1n),
        smooth = f.reduce((v, [p, e]) => v * p ** e, 1n);
      cyc = cyc.map((c) => gcd(c, smooth));
    }
  }
  while (cyc.at(-1) === 1n) cyc.pop();
  const n = cyc.length;
  if (!n) return !exact || bound === 1n ? [eye(N)] : [];
  const primes = Z_factor(cyc[0]!).map(([p]) => p);
  let p = primes[0]!,
    L: number[] = [],
    rank = -1;
  for (const q of primes) {
    const l = cyc.map((c) => val(c, q));
    const r = l.filter((v) => v > 0).length;
    if (r > rank) {
      rank = r;
      p = q;
      L = l;
    }
  }
  const powers = Array.from({ length: L[0]! + 1 }, (_, i) => p ** BigInt(i));
  let subq: bigint[][][] | null = null,
    expoI = 1n,
    indexq: bigint[] = [];
  if (primes.length === 1) {
    if (exact && bound! / p ** BigInt(val(bound!, p)) !== 1n) return [];
  } else {
    const I = cyc.map((c, i) => c / powers[L[i]!]!);
    while (I.at(-1) === 1n) I.pop();
    expoI = I[0]!;
    subq = boundedSubgroups(I, exact ? bound! / p ** BigInt(val(bound!, p)) : bound, exact).map(
      (h) => {
        const H = eye(n);
        for (let j = 0; j < h.length; j++)
          for (let i = 0; i < n; i++) H[j]![i] = i < h.length ? h[j]![i]! : 0n;
        return H;
      }
    );
    indexq = subq.map((h) => h.reduce((v, c, i) => v * c[i]!, 1n));
    subq = subq.map((H) => H.map((c) => c.map((x) => x * powers[L[0]!]!)));
  }
  const out: bigint[][][] = [],
    wG = L.reduce((a, b) => a + b, 0),
    target = bound === null ? 0 : wG - (exact ? val(bound, p) : logint0(bound, p));
  const record = (H: bigint[][], weight: number) => {
    const matrices: bigint[][][] =
      subq === null
        ? [H]
        : subq
            .filter((_, i) => bound === null || indexq[i]! <= bound / p ** BigInt(wG - weight))
            .map((I) => [...H.map((c) => c.map((x) => x * expoI)), ...I]);
    for (const matrix of matrices) {
      const h = ZM_hnfmodid(matrix, cyc),
        result = eye(N);
      for (let j = 0; j < n; j++)
        for (let i = 0; i < N; i++) result[j]![i] = i < n ? packedInteger(h[j]![i]!) : 0n;
      out.push(result);
    }
  };
  const types = new Array(n).fill(0);
  types[0] = -1;
  for (;;) {
    types[0]++;
    if (types[0] > L[0]!) {
      let j = 1;
      while (j < n && !(types[j] < L[j]! && types[j] < types[j - 1])) j++;
      if (j === n) break;
      types[j]++;
      for (let k = 0; k < j; k++) types[k] = types[j];
    }
    const weight = types.reduce((a, b) => a + b, 0);
    if (exact ? weight !== target : weight < target) continue;
    const t = types.findIndex((v) => v === 0) < 0 ? n : types.findIndex((v) => v === 0);
    if (!t) {
      record([new Array(n).fill(0n)], weight);
      continue;
    }
    if (n === 1) {
      record([[powers[L[0]! - types[0]]!]], weight);
      continue;
    }
    const c = new Array(n).fill(0),
      available = new Array(n).fill(true),
      maxc = types.slice(0, t).map((m) => L.filter((l) => l >= m).length);
    const dogroup = () => {
      let next = t;
      for (let i = 0; i < n; i++) if (available[i]) c[next++] = i;
      const entries: { i: number; r: number; max: bigint }[] = [];
      for (let i = 0; i < t; i++)
        for (let r = i + 1; r < n; r++) {
          const exponent =
            c[r] < c[i]
              ? types[i] - types[r] - 1
              : L[c[r]]! < types[i]
                ? L[c[r]]! - types[r]
                : types[i] - types[r];
          entries.push({ i, r, max: powers[exponent]! });
        }
      const a = entries.map(() => 1n);
      a[a.length - 1] = 0n;
      for (;;) {
        let j = a.length - 1;
        a[j] = a[j]! + 1n;
        if (a[j]! > entries[j]!.max) {
          j--;
          while (j >= 0 && a[j] === entries[j]!.max) j--;
          if (j < 0) break;
          a[j] = a[j]! + 1n;
          for (let k = j + 1; k < a.length; k++) a[k] = 1n;
        }
        const H = Array.from({ length: t }, () => new Array<bigint>(n).fill(0n));
        for (let i = 0; i < t; i++) H[i]![c[i]] = powers[L[c[i]]! - types[i]]!;
        for (let k = 0; k < entries.length; k++) {
          const { i, r } = entries[k]!,
            d = L[c[r]]! - types[i];
          H[i]![c[r]] = a[k]! * (c[r] < c[i] ? powers[d + 1]! : d > 0 ? powers[d]! : 1n);
        }
        record(H, weight);
      }
    };
    const loop = (r: number) => {
      if (r === t) {
        dogroup();
        return;
      }
      let j = r && types[r - 1] === types[r] ? c[r - 1] + 1 : 0;
      for (; j < maxc[r]!; j++)
        if (available[j]) {
          c[r] = j;
          available[j] = false;
          loop(r + 1);
          available[j] = true;
        }
    };
    loop(0);
  }
  return out;
}

/** addcell/packtoi discard initial zero limbs in the packed cell roundtrip. */
function packedInteger(value: bigint): bigint {
  while (value !== 0n && (value & ((1n << 64n) - 1n)) === 0n) value >>= 64n;
  return value;
}
/** PARI subgrouplist: [index] is exact, a scalar bounds the index, omission lists all.
 * @see Deviation: PARI unit-subgroup dependency boundaries
 */
export function subgrouplist(cyc: bigint[], bound?: bigint | readonly bigint[]): bigint[][][] {
  if (cyc.some((c, i) => c <= 0n || (i + 1 < cyc.length && c % cyc[i + 1] !== 0n)))
    throw new PariError('incorrect type in subgrouplist [not a finite group] (t_VEC).');
  if (Array.isArray(bound) && bound.length !== 1)
    throw new PariError('incorrect type in subgroup (t_VEC).');
  const exact = Array.isArray(bound),
    B = exact ? (bound[0] as bigint) : bound === undefined ? null : (bound as bigint);
  if (B !== null && B <= 0n) throw new PariError('domain error in subgroup: index bound <= 0');
  return boundedSubgroups(cyc, B, exact);
}
