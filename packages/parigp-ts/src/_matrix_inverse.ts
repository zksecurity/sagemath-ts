/** PARI Flv.c CUP/adjoint and alglin1.c modular integer inverse dependencies.
 * Matrices are row-oriented internally; products delegate to the native matrix port.
 * @see Deviation: Number-field ideal backend adapters
 */
import { Flm_mul } from './FpV.js';
import { fromColumns, toColumns, type Rows } from './_matrix_mul.js';
import { inverseCoefficient, residue } from './_polynomial_division.js';

const zero = (m: number, n: number): Rows =>
  Array.from({ length: m }, () => Array<bigint>(n).fill(0n));
const id = (n: number): Rows =>
  Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1n : 0n)));
const slice = (a: Rows, r: number, m: number, c: number, n: number): Rows =>
  a.slice(r, r + m).map((row) => row.slice(c, c + n));
const sub = (a: Rows, b: Rows, p: bigint): Rows =>
  a.map((row, i) => row.map((c, j) => residue(c - b[i]![j]!, p)));
const scale = (a: Rows, c: bigint, p: bigint): Rows =>
  a.map((row) => row.map((v) => residue(v * c, p)));
const mul = (a: Rows, b: Rows, m: number, n: number, k: number, p: bigint): Rows =>
  !m || !n || !k ? zero(m, k) : fromColumns(Flm_mul(toColumns(a, n), toColumns(b, k), p));
const inverse = (a: bigint, p: bigint): bigint => inverseCoefficient(a, p, true);

/** Flm_rsolve_upper_pre and Flm_lsolve_upper_pre, with their two-row basecases. */
function upper(U: Rows, B: Rows, p: bigint, left = false): Rows {
  const n = U.length,
    k = left ? B.length : (B[0]?.length ?? 0);
  if (!n) return B.map((row) => row.slice());
  if (n === 1) return scale(B, inverse(U[0]![0]!, p), p);
  if (n === 2) {
    const a = U[0]![0]!,
      b = U[0]![1]!,
      d = U[1]![1]!,
      di = inverse(residue(a * d, p), p),
      ai = residue(d * di, p),
      bi = residue(a * di, p);
    if (left)
      return B.map((row) => {
        const x = residue(row[0]! * ai, p);
        return [x, residue((row[1]! - x * b) * bi, p)];
      });
    const x2 = scale([B[1]!], bi, p),
      x1 = scale(sub([B[0]!], scale(x2, b, p), p), ai, p);
    return [...x1, ...x2];
  }
  const n1 = Math.ceil(n / 2),
    n2 = n - n1,
    U11 = slice(U, 0, n1, 0, n1),
    U12 = slice(U, 0, n1, n1, n2),
    U22 = slice(U, n1, n2, n1, n2);
  if (left) {
    const X1 = upper(U11, slice(B, 0, k, 0, n1), p, true),
      B2 = sub(slice(B, 0, k, n1, n2), mul(X1, U12, k, n1, n2, p), p),
      X2 = upper(U22, B2, p, true);
    return X1.map((row, i) => row.concat(X2[i]!));
  }
  const X2 = upper(U22, slice(B, n1, n2, 0, k), p),
    B1 = sub(slice(B, 0, n1, 0, k), mul(U12, X2, n1, n2, k, p), p),
    X1 = upper(U11, B1, p);
  return [...X1, ...X2];
}

/** Flm_rsolve_lower_unit_pre / Flm_lsolve_lower_unit_pre. */
function lower(L: Rows, B: Rows, p: bigint, left = false): Rows {
  const n = L.length,
    k = left ? B.length : (B[0]?.length ?? 0);
  if (n <= 1) return B.map((row) => row.slice());
  if (n === 2) {
    const c = L[1]![0]!;
    return left
      ? B.map((row) => [residue(row[0]! - row[1]! * c, p), row[1]!])
      : [B[0]!.slice(), B[1]!.map((v, j) => residue(v - c * B[0]![j]!, p))];
  }
  const n1 = Math.ceil(n / 2),
    n2 = n - n1,
    L11 = slice(L, 0, n1, 0, n1),
    L21 = slice(L, n1, n2, 0, n1),
    L22 = slice(L, n1, n2, n1, n2);
  if (left) {
    const X2 = lower(L22, slice(B, 0, k, n1, n2), p, true),
      B1 = sub(slice(B, 0, k, 0, n1), mul(X2, L21, k, n2, n1, p), p),
      X1 = lower(L11, B1, p, true);
    return X1.map((row, i) => row.concat(X2[i]!));
  }
  const X1 = lower(L11, slice(B, 0, n1, 0, k), p),
    B2 = sub(slice(B, n1, n2, 0, k), mul(L21, X1, n2, n1, k, p), p),
    X2 = lower(L22, B2, p);
  return [...X1, ...X2];
}

type CUP = { R: number[]; C: Rows; U: Rows; P: number[] };
/** Flm_CUP_pre, cutoff eight and the original rank-profile/column permutation. */
function cup(A: Rows, m: number, n: number, p: bigint): CUP {
  if (m < 8 || n < 8) {
    const a = A.map((row) => row.slice()),
      R: number[] = [],
      P = Array.from({ length: n }, (_, i) => i);
    let pr = -1;
    for (let j = 0; j < n; j++) {
      let pc = -1;
      for (pr++; pr < m; pr++) {
        for (let k = j; k < n; k++) if (pc < 0 && a[pr]![k] !== 0n) pc = k;
        if (pc >= 0) break;
      }
      if (pc < 0) break;
      R.push(pr);
      if (pc !== j) {
        for (const row of a) [row[j], row[pc]] = [row[pc]!, row[j]!];
        [P[j], P[pc]] = [P[pc]!, P[j]!];
      }
      const u = inverse(a[pr]![j]!, p);
      for (let i = pr + 1; i < m; i++) {
        const v = residue(a[i]![j]! * u, p);
        a[i]![j] = v;
        for (let k = j + 1; k < n; k++) a[i]![k] = residue(a[i]![k]! - a[pr]![k]! * v, p);
      }
    }
    return { R, C: a.map((row) => row.slice(0, R.length)), U: R.map((i) => a[i]!.slice()), P };
  }
  const m1 = Math.ceil(Math.min(m, n) / 2),
    m2 = m - m1,
    first = cup(A.slice(0, m1), m1, n, p),
    r1 = first.R.length,
    A2 = A.slice(m1);
  if (!r1) {
    const second = cup(A2, m2, n, p);
    return {
      ...second,
      R: second.R.map((i) => i + m1),
      C: [...zero(m1, second.R.length), ...second.C],
    };
  }
  const U11 = slice(first.U, 0, r1, 0, r1),
    U12 = slice(first.U, 0, r1, r1, n - r1),
    T21 = A2.map((row) => first.P.slice(0, r1).map((j) => row[j]!)),
    T22 = A2.map((row) => first.P.slice(r1).map((j) => row[j]!));
  const C21 = upper(U11, T21, p, true),
    B2 = sub(T22, mul(C21, U12, m2, r1, n - r1, p), p),
    second = cup(B2, m2, n - r1, p),
    r2 = second.R.length;
  return {
    R: first.R.concat(second.R.map((i) => i + m1)),
    C: [
      ...first.C.map((row) => row.concat(Array<bigint>(r2).fill(0n))),
      ...C21.map((row, i) => row.concat(second.C[i]!)),
    ],
    U: [
      ...U11.map((row, i) => row.concat(second.P.map((j) => U12[i]![j]!))),
      ...second.U.map((row) => Array<bigint>(r1).fill(0n).concat(row)),
    ],
    P: first.P.slice(0, r1).concat(second.P.map((j) => first.P[r1 + j]!)),
  };
}
function sign(P: number[]): bigint {
  let s = 1n;
  for (let i = 0; i < P.length; i++) for (let j = 0; j < i; j++) if (P[j]! > P[i]!) s = -s;
  return s;
}
function unpermute(a: Rows, P: number[]): Rows {
  const out: Rows = [];
  for (let i = 0; i < P.length; i++) out[P[i]!] = a[i]!;
  return out;
}

/** Flm_adjoint: full-rank solve or native rank-one outer product. */
export function wordMatrixAdjoint(A: Rows, p: bigint): Rows {
  const m = A.length,
    n = A[0]?.length ?? 0;
  if (p <= 1n || p >= 1n << 64n || A.some((row) => row.length !== n) || m < n)
    throw new RangeError(
      'adjoint requires a rectangular matrix with rows >= columns and a word prime'
    );
  if (!n) return [];
  const { R, C, U, P } = cup(
      A.map((row) => row.map((v) => residue(v, p))),
      m,
      n,
      p
    ),
    r = R.length;
  if (r === n) {
    const Y = lower(
        R.map((i) => C[i]!),
        R.map((i) => id(m)[i]!),
        p
      ),
      X = unpermute(upper(U, Y, p), P);
    let d = sign(P);
    for (let i = 0; i < n; i++) d = residue(d * U[i]![i]!, p);
    return scale(X, d, p);
  }
  if (r < n - 1) return zero(n, m);
  let q = n - 1;
  for (let i = 0; i < n - 1; i++)
    if (R[i] !== i) {
      q = i;
      break;
    }
  const c = lower(slice(C, 0, q, 0, q), [C[q]!.slice(0, q)], p, true)[0] ?? [],
    d = Array<bigint>(m).fill(0n);
  for (let i = 0; i < q; i++) d[i] = c[i]!;
  d[q] = p - 1n;
  const U1 = slice(U, 0, n - 1, 0, n - 1),
    v = upper(
      U1,
      U.map((row) => [row[n - 1]!]),
      p
    ).map((row) => row[0]!);
  v.push(p - 1n);
  let D = sign(P) * ((q + 1 + n) % 2 === 0 ? 1n : -1n);
  for (let i = 0; i < n - 1; i++) D = residue(D * U1[i]![i]!, p);
  return unpermute(
    v.map((c) => d.map((x) => residue(c * x * D, p))),
    P
  );
}

import { gcd } from './ff.js';
import { nextprime } from './ifactor.js';
import { sqrtremi } from './kernel/gmp/mp.js';
import { halfgcdii } from './kernel/none/halfgcd.js';
import { quadratic_prec_mask } from './Zp.js';
import { ZM_mul } from './ZV.js';
import { PariError } from './errors.js';
const absolute = (a: bigint) => (a < 0n ? -a : a);

/** Equal numerator/denominator bounds: Fp_ratlift_hgcd. */
export function rationalLift(x: bigint, m: bigint, bound: bigint): [bigint, bigint] | null {
  x = residue(x, m);
  if (x === 0n) return [0n, 1n];
  if (bound === 0n) return null;
  if (x <= bound) return [x, 1n];
  const [M, V] = halfgcdii(x, m);
  let [u, a] = V,
    y = M[0][0],
    b = M[1][0];
  while (absolute(b) <= bound) {
    if (absolute(a) <= bound) return b < 0n ? [-a, -b] : [a, b];
    const q = u / a;
    [u, a] = [a, u % a];
    [y, b] = [b, y - b * q];
  }
  return null;
}
function primitive(A: Rows): Rows {
  let g = 0n;
  for (const row of A) for (const x of row) g = gcd(g, x);
  return g > 1n ? A.map((row) => row.map((x) => x / g)) : A;
}
function integerProduct(A: Rows, B: Rows): Rows {
  return fromColumns(ZM_mul(toColumns(A, B.length), toColumns(B, B[0]?.length ?? 0)));
}
function reconstruct(A: Rows, H: Rows, mod: bigint): Rows | null {
  let D = 0n;
  for (let i = 0; i < A.length; i++) D += H[0]![i]! * A[i]![0]!;
  const g = gcd(D, mod);
  if (g !== 1n) mod /= g;
  if (mod < 1n) return null;
  const inv = inverseCoefficient(residue(D, mod), mod, false),
    bound = sqrtremi(mod >> 1n)[0];
  if (!rationalLift(residue(inv * A[0]![0]!, mod), mod, bound)) return null;
  const fractions: [bigint, bigint][][] = [];
  let denominator = 1n;
  for (const row of H) {
    const out: [bigint, bigint][] = [];
    for (const x of row) {
      const r = rationalLift(residue(x * inv, mod), mod, bound);
      if (!r) return null;
      const g = gcd(r[0], r[1]);
      r[0] /= g;
      r[1] /= g;
      denominator = (denominator / gcd(denominator, r[1])) * r[1];
      out.push(r);
    }
    fractions.push(out);
  }
  return primitive(fractions.map((row) => row.map(([n, d]) => n * (denominator / d))));
}

/** ZM_inv_i: small closed forms, modular adjoints, centered CRT and rational lifting.
 * The Hadamard iteration bound is evaluated exactly instead of through t_REAL.
 */
export function integerMatrixInverse(A: Rows): [Rows, bigint] | null {
  const m = A.length,
    n = A[0]?.length ?? 0;
  if (A.some((row) => row.length !== n)) throw new RangeError('inconsistent matrix dimensions');
  if (!n) return [[], 1n];
  if (m < n) return null;
  if (n === 1 && m === 1) {
    const x = A[0]![0]!;
    return x === 0n ? null : [[[x > 0n ? 1n : -1n]], absolute(x)];
  }
  if (n === 2 && m === 2) {
    let c = 0n;
    for (const row of A) for (const x of row) c = gcd(c, x);
    if (c === 0n) return null;
    const [[a, b], [cc, d]] = A.map((row) => row.map((x) => x / c)) as [
        [bigint, bigint],
        [bigint, bigint],
      ],
      D = a * d - b * cc;
    if (!D) return null;
    const s = D < 0n ? -1n : 1n;
    return [
      [
        [s * d, -s * b],
        [-s * cc, s * a],
      ],
      absolute(D) * c,
    ];
  }
  let squaredBound = 1n;
  for (let j = 0; j < n; j++) {
    let norm = 0n;
    for (let i = 0; i < m; i++) norm += A[i]![j]! ** 2n;
    squaredBound *= norm;
  }
  const [root, remainder] = sqrtremi(squaredBound),
    bound = root + (remainder ? 1n : 0n),
    bits = bound === 0n ? 0 : bound.toString(2).length - 1;
  let p = 1n << 63n,
    mod = 1n,
    H = zero(n, m),
    k1 = 1,
    k2 = 0;
  const combine = (count: number) => {
    for (let t = 0; t < count; t++) {
      p = nextprime(p + 1n);
      if (p >= 1n << 64n) throw new PariError('exhausted word primes in ZM_inv');
      const Hp = wordMatrixAdjoint(A, p),
        u = inverseCoefficient(mod % p, p, true),
        nextMod = mod * p;
      H = H.map((row, i) =>
        row.map((x, j) => {
          const r = x + mod * residue((Hp[i]![j]! - x) * u, p);
          return r > nextMod / 2n ? r - nextMod : r;
        })
      );
      mod = nextMod;
    }
  };
  combine(1);
  const primes = Math.floor((bits + 1) / (p.toString(2).length - 1)) + 1;
  let mask = primes <= 1 ? 3n : quadratic_prec_mask(primes);
  for (;;) {
    if (k2 > 0) {
      combine(k2);
      k1 += k2;
    }
    if (mask === 1n) break;
    k2 = mask & 1n ? k1 - 1 : k1;
    mask >>= 1n;
    let candidate = reconstruct(A, H, mod);
    if (candidate) {
      let R = integerProduct(candidate, A),
        d = R[0]![0]!;
      if (d < 0n) {
        d = -d;
        candidate = candidate.map((row) => row.map((x) => -x));
        R = R.map((row) => row.map((x) => -x));
      }
      if (R.every((row, i) => row.every((x, j) => x === (i === j ? d : 0n)))) return [candidate, d];
    }
  }
  let D = 0n;
  for (let i = 0; i < m; i++) D += H[0]![i]! * A[i]![0]!;
  if (D === 0n) throw new PariError('impossible inverse in ZM_inv');
  let g = absolute(D);
  for (const row of H) for (const x of row) g = gcd(g, x);
  if (D < 0n) g = -g;
  return [H.map((row) => row.map((x) => x / g)), D / g];
}


/** Native Flm_pivots, returning one-based pivot rows and the nullity.
 * @see Deviation: PARI word matrix pivot and solve adapters
 */
export function wordMatrixPivots(
  A: Rows,
  m: number,
  n: number,
  p: bigint
): [number[] | null, number] {
  if (A.length !== m || A.some((row) => row.length !== n))
    throw new RangeError('inconsistent word matrix dimensions');
  if (!n) return [null, 0];
  const a = A.map((row) => row.map((x) => residue(x, p))),
    d = Array<number>(n).fill(0);
  if (m >= 8 && n >= 8) {
    const { R, P } = cup(a, m, n, p);
    for (let i = 0; i < R.length; i++) d[P[i]!] = R[i]! + 1;
    return [d, n - R.length];
  }
  const used = Array<boolean>(m).fill(false);
  let r = 0;
  for (let k = 0; k < n; k++) {
    let j = 0;
    while (j < m && (used[j] || a[j]![k] === 0n)) j++;
    if (j === m) {
      r++;
      continue;
    }
    const pivot = p - inverse(a[j]![k]!, p);
    used[j] = true;
    d[k] = j + 1;
    for (let i = k + 1; i < n; i++) a[j]![i] = residue(pivot * a[j]![i]!, p);
    for (let t = 0; t < m; t++)
      if (!used[t]) {
        const c = a[t]![k]!;
        if (c) {
          a[t]![k] = 0n;
          for (let i = k + 1; i < n; i++) a[t]![i] = residue(a[t]![i]! + c * a[j]![i]!, p);
        }
      }
    for (let i = k; i < n; i++) a[j]![i] = 0n;
  }
  return [d, r];
}

/** Native Flm_gauss, with the cutoff-eight CUP and small elimination branches.
 * @see Deviation: PARI word matrix pivot and solve adapters
 */
export function wordMatrixSolve(
  A: Rows,
  B: Rows,
  m: number,
  n: number,
  k: number,
  p: bigint
): Rows | null {
  if (
    A.length !== m ||
    B.length !== m ||
    A.some((row) => row.length !== n) ||
    B.some((row) => row.length !== k)
  )
    throw new RangeError('inconsistent word solve dimensions');
  if (!n) return [];
  const a = A.map((row) => row.map((x) => residue(x, p))),
    b = B.map((row) => row.map((x) => residue(x, p)));
  if (n >= 8) {
    if (m < n) return null;
    const { R, C, U, P } = cup(a, m, n, p);
    if (R.length < n) return null;
    const Y = lower(
      R.map((i) => C[i]!),
      R.map((i) => b[i]!),
      p
    );
    return unpermute(upper(U, Y, p), P);
  }
  for (let i = 0; i < n; i++) {
    let j = i;
    while (j < m && a[j]![i] === 0n) j++;
    if (j === m) return null;
    a[j]![i] = inverse(a[j]![i]!, p);
    if (j !== i) {
      for (let col = i; col < n; col++) [a[i]![col], a[j]![col]] = [a[j]![col]!, a[i]![col]!];
      [b[i], b[j]] = [b[j]!, b[i]!];
    }
    if (i === n - 1) break;
    const invpiv = p - a[i]![i]!;
    for (let row = i + 1; row < m; row++) {
      if (a[row]![i] === 0n) continue;
      const factor = residue(a[row]![i]! * invpiv, p);
      for (let col = i + 1; col < n; col++)
        a[row]![col] = residue(
          a[row]![col]! + (factor === 1n ? a[i]![col]! : factor * a[i]![col]!),
          p
        );
      for (let col = 0; col < k; col++)
        b[row]![col] = residue(
          b[row]![col]! + (factor === 1n ? b[i]![col]! : factor * b[i]![col]!),
          p
        );
    }
  }
  const u = zero(n, k);
  for (let col = 0; col < k; col++)
    for (let i = n - 1; i >= 0; i--) {
      let s = b[i]![col]!;
      for (let j = i + 1; j < n; j++) s -= a[i]![j]! * u[j]![col]!;
      u[i]![col] = residue(s * a[i]![i]!, p);
    }
  return u;
}

/** Native column-wise rational lifting with denominator reuse (nffactor.c). */
function liftLinearSolution(H: Rows, modulus: bigint, n: number, k: number): [Rows, bigint] | null {
  const bound = sqrtremi(modulus >> 1n)[0];
  const fractions: [bigint, bigint][][] = Array.from({ length: n }, () => []);
  let denominator = 1n;
  const normalize = (a: bigint, b: bigint): [bigint, bigint] => {
    const g = gcd(a, b);
    return [a / g, b / g];
  };
  for (let j = 0; j < k; j++) {
    let previous = 0n;
    for (let i = 0; i < n; i++) {
      let t = H[i]![j]!,
        r: [bigint, bigint] | null = null;
      if (t < 0n) {
        if (-t < bound) r = [t, 1n];
        else t += modulus;
      }
      if (!r && previous) {
        let a = residue(t * previous, modulus);
        if (a > modulus / 2n) a -= modulus;
        if (absolute(a) < bound) r = normalize(a, previous);
      }
      if (!r) {
        r = rationalLift(t, modulus, bound);
        if (!r || gcd(r[0], r[1]) !== 1n) return null;
      }
      fractions[i]![j] = r;
      if (r[1] > previous && r[1] > 1n) previous = r[1];
      denominator = (denominator / gcd(denominator, r[1])) * r[1];
    }
  }
  return [fractions.map((row) => row.map(([a, b]) => a * (denominator / b))), denominator];
}
type LinearResidues = { value: Rows; modulus: bigint };
function combineLinearResidues(a: LinearResidues, b: LinearResidues): LinearResidues {
  if (a.modulus === 1n) return b;
  if (b.modulus === 1n) return a;
  const p = b.modulus,
    u = inverseCoefficient(residue(a.modulus, p), p, false),
    modulus = a.modulus * p;
  return {
    modulus,
    value: a.value.map((row, i) =>
      row.map((v, j) => {
        const x = v + a.modulus * residue((b.value[i]![j]! - v) * u, p);
        return x > modulus / 2n ? x - modulus : x;
      })
    ),
  };
}
function combineLinearTree(values: LinearResidues[], lo: number, hi: number): LinearResidues {
  if (hi - lo === 1) return values[lo]!;
  const mid = Math.floor((lo + hi) / 2);
  return combineLinearResidues(
    combineLinearTree(values, lo, mid),
    combineLinearTree(values, mid, hi)
  );
}
/** Native ZM_gauss modular solve, CRT doubling, rational lifting and certification.
 * Returns numerator rows and their positive common denominator, or null on rank failure.
 * @see Deviation: PARI modular integer rank and solve adapters
 */
export function integerMatrixSolve(
  A: Rows,
  B: Rows,
  m: number,
  n: number,
  k: number
): [Rows, bigint] | null {
  if (
    A.length !== m ||
    B.length !== m ||
    A.some((row) => row.length !== n) ||
    B.some((row) => row.length !== k)
  )
    throw new RangeError('inconsistent integer solve dimensions');
  if (!n) {
    if (k) throw new PariError('inconsistent dimensions in gauss');
    return [[], 1n];
  }
  if (m < n) throw new RangeError('integer matrix solve requires rows >= columns');
  if (!k) return [zero(n, 0), 1n];
  const [d, nullity] = integerMatrixPivots(A, m, n);
  if (nullity) return null;
  const selected = d!
    .filter((v) => v)
    .map((v) => v - 1)
    .sort((a, b) => a - b);
  const a = selected.map((i) => A[i]!),
    b = selected.map((i) => B[i]!);
  let prime = 1n << 63n,
    H: LinearResidues = { value: zero(n, k), modulus: 1n };
  for (let step = 1; ; step *= 2) {
    const batch: LinearResidues[] = [];
    for (let j = 0; j < Math.ceil(step / 2); j++) {
      prime = nextprime(prime + 1n);
      if (prime >= 1n << 64n) throw new PariError('overflow in ZM_gauss [ran out of primes]');
      const hp = wordMatrixSolve(a, b, n, n, k, prime);
      batch.push(
        hp === null
          ? { value: zero(n, k), modulus: 1n }
          : {
              value: hp.map((row) => row.map((v) => (v > prime / 2n ? v - prime : v))),
              modulus: prime,
            }
      );
    }
    H = combineLinearResidues(H, combineLinearTree(batch, 0, batch.length));
    const candidate = liftLinearSolution(H.value, H.modulus, n, k);
    if (!candidate) continue;
    const [X, denominator] = candidate,
      product = integerProduct(a, X);
    if (product.every((row, i) => row.every((v, j) => v === b[i]![j]! * denominator)))
      return candidate;
  }
}
/** Native ZM_pivots: word-prime rank profiles followed by exact kernel certification.
 * @see Deviation: PARI modular integer rank and solve adapters
 */
export function integerMatrixPivots(A: Rows, m: number, n: number): [number[] | null, number] {
  if (A.length !== m || A.some((row) => row.length !== n))
    throw new RangeError('inconsistent integer matrix dimensions');
  if (!n) return [null, 0];
  let zeroColumns = 0;
  for (let j = 0; j < n; j++) if (A.every((row) => row[j] === 0n)) zeroColumns++;
  if (zeroColumns === n) return [Array<number>(n).fill(0), n];
  const minimum = Math.max(zeroColumns, n - m),
    nn = Math.min(m, n),
    imax = nn < 16 ? 1 : nn < 64 ? 2 : 3;
  let best = n,
    pivots: number[] | null = null,
    beenThere = false,
    prime = 1n << 31n;
  const complete = (selected: number[], length: number): number[] => {
    const used = new Set(selected),
      first = [...used].sort((a, b) => a - b),
      rest: number[] = [];
    for (let i = length - 1; i >= 0; i--) if (!used.has(i)) rest.push(i);
    return first.concat(rest);
  };
  for (;;) {
    for (let i = 0; ; i++) {
      prime = nextprime(prime + 1n);
      if (prime >= 1n << 64n) throw new PariError('overflow in ZM_pivots [ran out of primes]');
      const [d, r] = wordMatrixPivots(A, m, n, prime);
      if (r === minimum) return [d, r];
      if (r < best) {
        best = r;
        pivots = d;
        if (beenThere) break;
      }
      if (!beenThere && i >= imax) break;
    }
    beenThere = true;
    const rows = complete(
      (pivots ?? []).filter((v) => v).map((v) => v - 1),
      m
    );
    const cols = complete(
      (pivots ?? []).flatMap((v, j) => (v ? [j] : [])),
      n
    );
    let M = rows.map((i) => cols.map((j) => A[i]![j]!));
    if (n > m) M = cols.map((_, j) => rows.map((_, i) => M[i]![j]!));
    const rank = n - best,
      IM = M.map((row) => row.slice(0, rank)),
      KM = M.map((row) => row.slice(rank));
    const solution = integerMatrixSolve(
      IM.slice(0, rank),
      KM.slice(0, rank),
      rank,
      rank,
      nn - rank
    )!;
    const [X, denominator] = solution,
      rhs = KM.slice(rank),
      product = integerProduct(IM.slice(rank), X);
    if (
      product.length === rhs.length &&
      product.every(
        (row, i) =>
          row.length === rhs[i]!.length && row.every((v, j) => v === rhs[i]![j]! * denominator)
      )
    )
      return [pivots, best];
  }
}
