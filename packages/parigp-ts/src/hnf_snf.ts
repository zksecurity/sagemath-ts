import { gcd, xgcd } from './ff.js';
import { PariError } from './errors.js';
/** Column echelon reduction from PARI hnf_snf.c:1018–1154. */
import { inverseCoefficient, residue } from './_polynomial_division.js';
import { diviiround } from './gen3.js';

/** Center the residues of a square column HNF (native ZM_hnfcenter).
 * The typed adapter copies its input; native PARI may reuse its allocation.
 * @see Deviation: Number-field ideal coercion and centered integral bases
 */
export function ZM_hnfcenter(input: bigint[][]): bigint[][] {
  const n = input.length;
  if (input.some((c, j) => c.length !== n || c[j]! <= 0n))
    throw new RangeError('ZM_hnfcenter requires a square HNF with positive diagonal');
  const M = input.map((c) => c.slice());
  for (let j = n - 2; j >= 0; j--) {
    const a = M[j]![j]!;
    for (let k = j + 1; k < n; k++) {
      const q = diviiround(M[k]![j]!, a);
      if (q === 0n) continue;
      for (let i = 0; i <= j; i++) M[k]![i] -= q * M[j]![i]!;
    }
  }
  return M;
}

/** Matrices use an array of columns, as in PARI. */
function echelon(
  input: bigint[][],
  early: boolean,
  p: bigint,
  pm: bigint,
  word: boolean
): bigint[][] | null {
  const x = input.map((c) => c.slice());
  if (!x.length) return x;
  const rows = x[0]!.length;
  let m = 0,
    power = pm;
  while (power !== 0n && power % p === 0n) {
    m++;
    power /= p;
  }
  let lower = Math.max(rows - x.length, 0),
    active = x.length;
  for (let i = rows - 1; i >= lower; i--) {
    let min = Infinity,
      pivot = -1,
      unit = 0n;
    for (let k = 0; k < active; k++) {
      let u = x[k]![i]!;
      if (u === 0n) continue;
      let v = 0;
      while (u % p === 0n) {
        v++;
        u /= p;
      }
      if (v >= m) x[k]![i] = 0n;
      else if (v < min) {
        min = v;
        pivot = k;
        unit = u;
        if (v === 0) break;
      }
    }
    if (pivot < 0) {
      if (early) return null;
      x[active - 1]![i] = 0n;
      lower = Math.max(lower - 1, 0);
      continue;
    }
    [x[pivot], x[active - 1]] = [x[active - 1]!, x[pivot]!];
    const q = min ? p ** BigInt(m - min) : pm;
    unit = residue(unit, q);
    const col = x[active - 1]!;
    if (unit !== 1n) {
      const inv = inverseCoefficient(unit, q, word);
      for (let r = 0; r < i; r++) if (col[r]) col[r] = residue(col[r]! * inv, pm);
    }
    const pv = p ** BigInt(min);
    col[i] = pv;
    for (let j = active - 2; j >= 0; j--) {
      const a = word ? x[j]![i]! : (x[j]![i] = residue(x[j]![i]!, pm));
      if (a === 0n) continue;
      const t = word ? residue(-a / pv, q) : -a / pv;
      for (let r = 0; r < rows; r++) {
        const c = x[j]![r]! + t * col[r]!;
        x[j]![r] = word ? residue(c, pm) : c;
      }
    }
    active--;
  }
  return x.length > rows ? x.slice(x.length - rows) : x;
}
export function ZpM_echelon(
  x: bigint[][],
  early: boolean,
  p: bigint,
  pm: bigint
): bigint[][] | null {
  return echelon(x, early, p, pm, false);
}
export function zlm_echelon(
  x: bigint[][],
  early: boolean,
  p: bigint,
  pm: bigint
): bigint[][] | null {
  return echelon(x, early, p, pm, true);
}

/** PARI hnf_snf.c:1180-1378, diagonal-modulus HNF, column-major storage.
 * Coefficients reduce eagerly after column elimination; native storage-size
 * thresholds only postpone these same reductions.
 */
export function ZM_hnfmodid(input: bigint[][], moduli: bigint[]): bigint[][] {
  const n = moduli.length;
  if (input.length && input.some((c) => c.length !== n))
    throw new PariError('inconsistent dimensions in ZM_hnfmod.');
  if (!n) return [];
  let x: bigint[][] = [[], ...input.map((c) => [0n, ...c])],
    D = [0n, ...moduli];
  const diagonal = (i: number) => Array.from({ length: n + 1 }, (_, r) => (r === i ? D[i]! : 0n));
  if (!input.length) x = [[], ...moduli.map((_, i) => diagonal(i + 1))];
  const elem = (a: bigint, b: bigint, j: number, k: number) => {
    if (!b) {
      [x[j], x[k]] = [x[k]!, x[j]!];
      return;
    }
    const [d, u, v] = xgcd(a, b),
      A = x[j]!,
      B = x[k]!;
    if (!u) {
      x[j] = A.map((c, r) => c - (a / b) * B[r]!);
      return;
    }
    if (!v) {
      x[k] = A;
      x[j] = B.map((c, r) => c - (b / a) * A[r]!);
      return;
    }
    x[k] = A.map((c, r) => u * c + v * B[r]!);
    x[j] = A.map((c, r) => (b / d) * c - (a / d) * B[r]!);
  };
  const reduce = (c: bigint[], rows: number) => {
    for (let r = 1; r <= rows; r++) c[r] = c[r]! % D[r]!;
  };
  const invgen = (a: bigint, d: bigint): [bigint, bigint] => {
    a = residue(a, d);
    if (!a) return [d, 0n];
    let [g, v] = xgcd(a, d);
    if (g === 1n) return [g, v];
    let e = d / g,
      d0 = g;
    for (;;) {
      const t = gcd(d0, e);
      if (t === 1n) break;
      d0 /= t;
    }
    if (d0 === 1n) return [g, v];
    if (g !== d0) e = (e / gcd(e, g / d0)) * (g / d0);
    const inv = xgcd(e, d0)[1];
    return [g, residue(v + e * residue((1n - v) * inv, d0), e * d0)];
  };
  const same = moduli.every((d) => d === moduli[0]);
  let co = x.length,
    li = n + 1,
    ldef = Math.max(li - co, 0);
  for (let def = co - 1, i = n; i > ldef; i--, def--) {
    const d = D[i]!;
    let addN = same;
    for (let j = 1; j < def; j++) {
      let a = (x[j]![i] = x[j]![i]! % d);
      if (!a) continue;
      const k = j + 1;
      const b = (x[k]![i] = x[k]![i]! % d);
      if (!b) {
        [x[j], x[k]] = [x[k]!, x[j]!];
        continue;
      }
      if (addN) {
        addN = false;
        if (a !== 1n) {
          const [g, u] = invgen(a, d);
          a = g;
          for (let r = 1; r < i; r++) x[j]![r] = x[j]![r]! * u;
          reduce(x[j]!, i - 1);
          x[j]![i] = a;
        }
      }
      elem(a, b, j, k);
      reduce(x[j]!, i - 1);
      reduce(x[k]!, i - 1);
    }
    if (!x[def]![i]) {
      x.splice(def + 1, 0, diagonal(i));
      ldef = Math.max(ldef - 1, 0);
      co++;
      def++;
    }
  }
  if (co < li)
    x = [
      [],
      ...Array.from({ length: li - co }, (_, k) => diagonal(k + 1)),
      ...x.slice(1),
      new Array(n + 1).fill(0n),
    ];
  else x = [[], ...x.slice(co - li + 1), new Array(n + 1).fill(0n)];
  x[1]![1] = gcd(x[1]![1]!, D[1]!);
  const D2 = D.slice();
  D2[1] = x[1]![1]!;
  for (let i = 2; i < n; i++) {
    let c = D2[i - 1]! * x[i]![i]!;
    if (c < 0n) c = -c;
    if (c >= D[i]!) break;
    D2[i] = c;
  }
  for (let i = n; i > 0; i--) {
    x[li]![i] = D[i]!;
    for (let j = i; j > 0; j--) {
      const a = x[li]![j]!;
      if (!a) continue;
      elem(a, x[j]![j]!, li, j);
      reduce(x[li]!, j - 1);
      reduce(x[j]!, j - 1);
    }
  }
  D = D2;
  x.pop();
  for (let i = n; i > 0; i--) {
    let diag = x[i]![i]!;
    if (diag < 0n) {
      x[i] = x[i]!.map((v) => -v);
      diag = -diag;
    }
    if (i !== n) reduce(x[i]!, i - 1);
    for (let j = i + 1; j < li; j++) {
      const b = (x[j]![i] = x[j]![i]! % D[i]!),
        r = residue(b, diag),
        q = (b - r) / diag;
      if (!q) continue;
      for (let k = 1; k < i; k++) x[j]![k] = x[j]![k]! - q * x[i]![k]!;
      x[j]![i] = r;
    }
  }
  return x.slice(1).map((c) => c.slice(1));
}

/** Native integral triangular solve; both arguments use zero-indexed columns. */
export function hnf_divscale(A: bigint[][], B: bigint[][], t: bigint): bigint[][] {
  const n = A.length;
  if (A.some((c) => c.length !== n) || B.length !== n || B.some((c) => c.length !== n))
    throw new RangeError('hnf_divscale requires square matrices of the same dimension');
  return B.map((b) => {
    const u = Array<bigint>(n).fill(0n);
    for (let i = n - 1; i >= 0; i--) {
      let m = b[i]! * t;
      for (let j = i + 1; j < n; j++) m -= A[j]![i]! * u[j]!;
      const d = A[i]![i]!;
      if (!d || m % d !== 0n) throw new RangeError('hnf_divscale requires an integral solution');
      u[i] = m / d;
    }
    return u;
  });
}

import {
  residue as primeResidue,
  inverseCoefficient as primeInverse,
} from './_polynomial_division.js';
/** Native FpM_echelon/FpM_hnfend for prime modulus HNF, in zero-indexed columns.
 * @see Deviation: PARI prime-decomposition matrix adapters
 */
export function ZM_hnfmodprime(input: bigint[][], p: bigint): bigint[][] {
  if (!input.length) return [];
  const rows = input[0]!.length,
    cols = input.length,
    x = input.map((c) => c.map((v) => primeResidue(v, p))),
    P: number[] = [];
  let def = cols - 1,
    ldef = Math.max(0, rows - cols);
  for (let i = rows - 1; i >= ldef; i--) {
    let k = def;
    for (; k >= 0; k--) if (x[k]![i] !== 0n) break;
    if (k < 0) {
      ldef = Math.max(0, ldef - 1);
      continue;
    }
    P.push(i);
    if (k !== def) [x[k], x[def]] = [x[def]!, x[k]!];
    const inv = primeInverse(x[def]![i]!, p, false);
    for (let j = 0; j < i; j++) x[def]![j] = primeResidue(x[def]![j]! * inv, p);
    x[def]![i] = 1n;
    for (let j = def - 1; j >= 0; j--) {
      const u = x[j]![i]!;
      if (!u) continue;
      for (let k = 0; k < rows; k++) x[j]![k]! -= u * x[def]![k]!;
      for (let k = 0; k < i; k++) x[j]![k] = primeResidue(x[j]![k]!, p);
    }
    def--;
  }
  P.sort((a, b) => a - b);
  const y = Array.from({ length: rows }, (_, j) =>
    Array.from({ length: rows }, (_, i) => (i === j ? p : 0n))
  );
  for (let j = 0; j < P.length; j++) y[P[j]!] = x[cols - P.length + j]!;
  for (let i = rows - 1; i >= 0; i--) {
    if (y[i]![i] === 1n || y[i]![i] === -1n) {
      for (let j = i + 1; j < rows; j++) {
        const b = y[j]![i]!;
        if (!b) continue;
        for (let k = 0; k < rows; k++) y[j]![k]! -= b * y[i]![k]!;
        for (let k = 0; k < i; k++)
          if (y[j]![k]! >= 1n << 64n || y[j]![k]! <= -(1n << 64n)) y[j]![k]! %= p;
      }
    } else for (let j = i + 1; j < rows; j++) y[j]![i] = primeResidue(y[j]![i]!, p);
  }
  return y;
}

/** Native column Hermite reduction with row permutation.
 * @see Deviation: PARI permuted Hermite and knapsack adapters
 */
export function ZM_hnfperm(
  input: bigint[][],
  withU = false,
  withPerm = false
): [bigint[][], bigint[][] | null, number[] | null] {
  const n = input.length,
    m = input[0]?.length ?? 0;
  if (input.some((c) => c.length !== m))
    throw new RangeError('ZM_hnfperm requires rectangular columns');
  if (!n) return [[], withU ? [] : null, withPerm ? [] : null];
  const A = input.map((c) => c.slice());
  const U: bigint[][] | null = withU
    ? Array.from({ length: n }, (_, j) => Array.from({ length: n }, (_, i) => (i === j ? 1n : 0n)))
    : null;
  const used = new Array<boolean>(m).fill(false),
    pivot = new Array<number>(n).fill(-1),
    perm: number[] = [];
  const abs = (x: bigint) => (x < 0n ? -x : x);
  const floor = (a: bigint, b: bigint) => a / b - (a < 0n && a % b ? 1n : 0n);
  const {add, neg, elem} = hnfOperations(A, U);
  for (let k = 0; k < n; k++) {
    for (let j = 0; j < k; j++) {
      const t = pivot[j]!;
      if (t < 0 || !A[k]![t]) continue;
      elem(A[k]![t]!, A[j]![t]!, k, j);
      if (A[j]![t]! < 0n) neg(j);
      const d = A[j]![t]!;
      for (let j1 = 0; j1 < j; j1++) {
        if (pivot[j1]! < 0) continue;
        const q = floor(A[j1]![t]!, d);
        if (q) add(j1, j, -q);
      }
    }
    let t = m - 1;
    while (t >= 0 && (used[t] || !A[k]![t])) t--;
    if (t < 0) continue;
    let p = A[k]![t]!;
    for (let i = t - 1; i >= 0; i--) {
      const q = A[k]![i]!;
      if (q && abs(p) > abs(q)) {
        p = q;
        t = i;
      }
    }
    perm.push(t);
    pivot[k] = t;
    used[t] = true;
    if (p < 0n) {
      neg(k);
      p = A[k]![t]!;
    }
    for (let j = 0; j < k; j++) {
      if (pivot[j]! < 0) continue;
      const q = floor(A[j]![t]!, p);
      if (q) add(j, k, -q);
    }
  }
  const rank = perm.length;
  for (let i = 0; i < m; i++) if (!used[i]) perm.push(i);
  perm.reverse();
  const H: bigint[][] = new Array(rank),
    resultU: bigint[][] | null = U ? new Array(n) : null;
  for (let j = 0, k = rank - 1, t = 0; j < n; j++) {
    if (pivot[j]! >= 0) {
      if (resultU) resultU[k + n - rank] = U![j]!;
      H[k--] = perm.map((i) => A[j]![i]!);
    } else if (resultU) resultU[t++] = U![j]!;
  }
  return [H, resultU, withPerm ? perm.map((i) => i + 1) : null];
}
/** Native hnfperm wrapper: H, full transformation U and one-based row permutation. */
export function hnfperm(input: bigint[][]): [bigint[][], bigint[][], number[]] {
  return ZM_hnfperm(input, true, true) as [bigint[][], bigint[][], number[]];
}
/** Native knapsack lattice recognition; null means the HNF fails its row test.
 * @see Deviation: PARI permuted Hermite and knapsack adapters
 */
export function ZM_hnf_knapsack(input: bigint[][]): bigint[][] | null {
  const [H, , perm] = ZM_hnfperm(input, false, true);
  if (!H.length) throw new PariError('bug in PARI/GP (Segmentation Fault), please report');
  const m = H[0]!.length;
  for (let i = 0; i < m; i++) {
    let found = false;
    for (const c of H) {
      const t = c[i]!;
      if (t) {
        if (found || (t !== 1n && t !== -1n)) return null;
        found = true;
      }
    }
  }
  const inverse = new Array<number>(m);
  for (let i = 0; i < m; i++) inverse[perm![i]! - 1] = i;
  return H.map((c) => inverse.map((i) => c[i]!));
}

/** Shared native column operations, including the exact Bézout special cases. */
function hnfOperations(A: bigint[][], U: bigint[][] | null) {
  const add = (j: number, k: number, q: bigint) => {
    A[j] = A[j]!.map((v, i) => v + q * A[k]![i]!);
    if (U) U[j] = U[j]!.map((v, i) => v + q * U[k]![i]!);
  };
  const swap = (j: number, k: number) => {
    [A[j], A[k]] = [A[k]!, A[j]!];
    if (U) [U[j], U[k]] = [U[k]!, U[j]!];
  };
  const neg = (j: number) => {
    A[j] = A[j]!.map(v => -v);
    if (U) U[j] = U[j]!.map(v => -v);
  };
  const elem = (a: bigint, b: bigint, j: number, k: number) => {
    if (!b) { swap(j,k); return; }
    const [d,u,v] = xgcd(a,b);
    if (!u) { add(j,k,-a/b); return; }
    if (!v) { add(k,j,-b/a); swap(j,k); return; }
    for (const M of U ? [A,U] : [A]) {
      const J=M[j]!,K=M[k]!;
      M[k]=J.map((c,i)=>u*c+v*K[i]!);
      M[j]=J.map((c,i)=>(b/d)*c-(a/d)*K[i]!);
    }
  };
  const reduce = (i: number, j0: number) => {
    if (A[j0]![i]! < 0n) neg(j0);
    const d=A[j0]![i]!;
    for(let j=j0+1;j<A.length;j++) {
      const x=A[j]![i]!,q=x/d-(x<0n&&x%d?1n:0n);
      if(q)add(j,j0,-q);
    }
  };
  return {add,swap,neg,elem,reduce};
}

/** Native hnf_i, retaining dependent columns when remove is false.
 * @see Deviation: PARI general Hermite form adapters
 */
function hnf_i(input: bigint[][], remove = true): bigint[][] {
  const n=input.length,m=input[0]?.length??0;
  if(input.some(c=>c.length!==m))throw new RangeError('hnf_i requires rectangular columns');
  if(!n)return [];
  const A=input.map(c=>c.slice()),{elem,neg,reduce}=hnfOperations(A,null);
  let def=n,lower=Math.max(m-n,0);
  for(let row=m-1;row>=lower;row--) {
    for(let j=def-2;j>=0;j--) {
      const a=A[j]![row]!;
      if(!a)continue;
      const k=j===0?def-1:j-1;
      elem(a,A[k]![row]!,j,k);
    }
    const s=A[def-1]![row]!;
    if(s) {
      if(s<0n)neg(def-1);
      reduce(row,def-1);def--;
    } else if(lower)lower--;
  }
  return remove?A.slice(def):A;
}

/** Native general Hermite form with optional transformation output.
 * remove=0 retains zero columns; 1 removes them only from H; 2 from H and U.
 * @see Deviation: PARI general Hermite form adapters
 */
export function ZM_hnfall_i(input: bigint[][], withU=false, remove:0|1|2=1):[bigint[][],bigint[][]|null] {
  const n=input.length,m=input[0]?.length??0;
  if(input.some(c=>c.length!==m))throw new RangeError('ZM_hnfall_i requires rectangular columns');
  if(!n)return [[],withU?[]:null];
  const A=input.map(c=>c.slice()),U:bigint[][]|null=withU?Array.from({length:n},(_,j)=>Array.from({length:n},(_,i)=>i===j?1n:0n)):null;
  const {elem,neg,reduce,swap}=hnfOperations(A,U),
    c=new Array<number>(m).fill(-1),h=new Array<number>(n).fill(m-1);
  let r=n;
  for(let row=m-1;row>=0;row--) {
    let j=0;
    for(;j<r;j++) {
      for(let i=h[j]!;i>row;i--) {
        const a=A[j]![i]!,k=c[i]!;
        if(a)elem(a,A[k]![i]!,j,k);
        reduce(i,k);
      }
      if(A[j]![row])break;
      h[j]=row-1;
    }
    if(j===r)continue;
    r--;
    if(j<r) {
      swap(j,r);h[j]=h[r]!;h[r]=row;c[row]=r;
    }
    if(A[r]![row]!<0n)neg(r);
    reduce(row,r);
  }
  for(let j=0;j<r;j++)for(let i=h[j]!;i>=0;i--) {
    const a=A[j]![i]!,k=c[i]!;
    if(a)elem(a,A[k]![i]!,j,k);
    reduce(i,k);
  }
  return [remove?A.slice(r):A,U&&remove===2?U.slice(r):U];
}
/** Native ownership wrapper; JavaScript arrays already own the returned storage.
 * @see Deviation: PARI general Hermite form adapters
 */
export function ZM_hnfall(input:bigint[][],withU=false,remove:0|1|2=1):[bigint[][],bigint[][]|null] {
  return ZM_hnfall_i(input,withU,remove);
}
/** Native full transformation wrapper, with zero columns removed only from H.
 * @see Deviation: PARI general Hermite form adapters
 */
export function hnfall(input:bigint[][]):[bigint[][],bigint[][]] {
  return ZM_hnfall(input,true,1) as [bigint[][],bigint[][]];
}
/** Native dispatch: hnf_i for up to seven columns, general hnfall thereafter.
 * @see Deviation: PARI general Hermite form adapters
 */
export function ZM_hnf(input:bigint[][]):bigint[][] {
  return input.length>7?ZM_hnfall(input,false,1)[0]:hnf_i(input,true);
}
