/** Native PARI lll.c floating-point reduction stages, with column-oriented bases. */
import {
  itor,
  setexpo,
  type MpReal,
  addrr,
  subrr,
  mulrr,
  sqrr,
  mulir,
  rtor,
  divrr,
  shiftr,
  subir,
} from './qfb.js';
import { roundr_safe } from './gen3.js';
import { mantissa2nr } from './kernel/gmp/mp.js';
import { cmprr, abscmprr } from './kernel/none/cmp.js';
import { PariError } from './errors.js';
import { fma, frexp, ldexp } from './_binary64.js';
import { dbltor, rtodbl } from './kernel/none/mp_indep.js';

type Matrix = bigint[][];
// AArch64 FCVTZS saturates overflow and maps NaN to zero. Native C leaves these
// out-of-range casts undefined; preserve the audited ABI in partial failure states.
const signedWord = (x: number): bigint =>
  Number.isNaN(x)
    ? 0n
    : x >= 2 ** 63
      ? (1n << 63n) - 1n
      : x <= -(2 ** 63)
        ? -(1n << 63n)
        : BigInt(Math.trunc(x));
const rint = (x: number): number => {
  const a = Math.abs(x);
  if (a >= 2 ** 52) return x;
  const r = 2 ** 52 + a - 2 ** 52;
  return x < 0 ? -r : r;
};
function rotate<T>(v: T[], last: number, first: number): void {
  const x = v[last]!;
  for (let i = last; i > first; i--) v[i] = v[i - 1]!;
  v[first] = x;
}
function updateAlpha(alpha: number[], kappa: number, kappa2: number, max: number): void {
  for (let i = kappa; i < kappa2; i++) if (kappa <= alpha[i]!) alpha[i] = kappa;
  for (let i = kappa2; i > kappa; i--) alpha[i] = alpha[i - 1]!;
  for (let i = kappa2 + 1; i <= max; i++) if (kappa < alpha[i]!) alpha[i] = kappa;
  alpha[kappa] = kappa;
}
function rotateGram<T>(G: T[][], last: number, first: number, max: number): void {
  const tmp: T[] = [];
  for (let i = 1; i <= last; i++) tmp[i] = G[last]![i]!;
  for (let i = last + 1; i <= max; i++) tmp[i] = G[i]![last]!;
  for (let i = last; i > first; i--) {
    for (let j = 1; j < first; j++) G[i]![j] = G[i - 1]![j]!;
    G[i]![first] = tmp[i - 1]!;
    for (let j = first + 1; j <= i; j++) G[i]![j] = G[i - 1]![j - 1]!;
    for (let j = last + 1; j <= max; j++) G[j]![i] = G[j]![i - 1]!;
  }
  for (let i = 1; i < first; i++) G[first]![i] = tmp[i]!;
  G[first]![first] = tmp[last]!;
  for (let i = last + 1; i <= max; i++) G[i]![first] = tmp[i]!;
}

/** fplll_fast and Babai_fast. Return the native status, transformed basis and U.
 * The native routine mutates B/U; this adapter starts with copies and returns them,
 * including partially reduced matrices when the native stage reports -1.
 * This is an uncertified internal stage, not a complete LLL implementation. Native
 * keep-first reduction can fail to terminate on dependent bases.
 * @see Deviation: PARI fast LLL floating-point profile
 * @see reference/pari/src/basemath/lll.c:698-1127
 */
export function fplll_fast(
  input: Matrix,
  delta = 0.99,
  eta = 0.51,
  keepfirst = false,
  trackTransform = true
): [number, Matrix, Matrix | null] {
  const d = input.length,
    n = input[0]?.length ?? 0;
  if (!d || !n || input.some((c) => c.length !== n))
    throw new RangeError('fplll_fast requires a nonempty rectangular basis');
  const B: Matrix = [[], ...input.map((c) => [0n, ...c])],
    U: Matrix | null = trackTransform
      ? [
          [],
          ...Array.from({ length: d }, (_, j) => [
            0n,
            ...Array.from({ length: d }, (_, i) => (i === j ? 1n : 0n)),
          ]),
        ]
      : null;
  const mat = () => Array.from({ length: d + 1 }, () => Array<number>(d + 1).fill(0));
  const G = mat(),
    mu = mat(),
    r = mat(),
    s = Array<number>(d + 1).fill(0),
    alpha = Array<number>(d + 1).fill(0),
    expoB = Array<number>(d + 1).fill(0),
    appB = Array.from({ length: d + 1 }, () => Array<number>(n + 1).fill(0));
  const setLine = (j: number) => {
    let max = 0;
    const exponents: number[] = [];
    for (let i = 1; i <= n; i++) {
      const x = itor(B[j]![i]!, 64);
      appB[j]![i] = rtodbl(setexpo(x, 0));
      exponents[i] = x.e;
      if (x.e > max) max = x.e;
    }
    for (let i = 1; i <= n; i++) appB[j]![i] = ldexp(appB[j]![i]!, exponents[i]! - max);
    expoB[j] = max;
  };
  // AArch64 Clang 21 -O3 contracts scalar tails but computes blocks of eight
  // products separately. dbldotsquare contracts every term. These choices alter
  // native partial failure states; do not replace with an arbitrary dot product.
  const dot = (a: number, b: number, square = false) => {
    let sum = appB[a]![1]! * appB[b]![1]!;
    const separate = square ? 1 : 1 + Math.floor((n - 1) / 8) * 8;
    for (let i = 2; i <= n; i++)
      sum =
        i <= separate ? sum + appB[a]![i]! * appB[b]![i]! : fma(appB[a]![i]!, appB[b]![i]!, sum);
    return sum;
  };
  const setG = (k: number, a: number, b: number) => {
    for (let i = a; i <= b; i++) G[k]![i] = dot(k, i);
  };
  const setG2 = (k: number, a: number, b: number) => {
    for (let i = a; i <= b; i++) G[i]![k] = dot(k, i);
  };
  const result = (status: number): [number, Matrix, Matrix | null] => [
    status,
    B.slice(1).map((c) => c.slice(1)),
    U ? U.slice(1).map((c) => c.slice(1)) : null,
  ];
  for (let i = 1; i <= d; i++) setLine(i);
  let kappamax = 1,
    maxG = d,
    i = 1;
  do {
    G[i]![i] = dot(i, i, true);
  } while (G[i]![i]! <= 0 && ++i <= d);
  let zeros = i - 1,
    kappa = i;
  if (zeros < d) r[zeros + 1]![zeros + 1] = G[zeros + 1]![zeros + 1]!;
  for (i = zeros + 1; i <= d; i++) alpha[i] = 1;

  const babai = (stage: number, a: number): boolean => {
    let aa = a > zeros ? a : zeros + 1,
      emaxmu = -2,
      emax2mu = -2,
      didSomething = false;
    for (;;) {
      let goOn = false;
      const emax3mu = emax2mu;
      emax2mu = emaxmu;
      emaxmu = -2;
      for (let j = aa; j < stage; j++) {
        let g = G[stage]![j]!;
        // The native GSO loop has the same eight-product contraction boundary.
        const separate = zeros + Math.floor((j - zeros - 1) / 8) * 8;
        for (let k = zeros + 1; k < j; k++)
          g = k <= separate ? g - mu[j]![k]! * r[stage]![k]! : fma(-mu[j]![k]!, r[stage]![k]!, g);
        r[stage]![j] = g;
        mu[stage]![j] = g / r[j]![j]!;
        emaxmu = Math.max(emaxmu, expoB[stage]! - expoB[j]!);
      }
      if (emax3mu !== -2 && emax3mu <= emax2mu + 5) return true;
      for (let j = stage - 1; j > zeros; j--)
        if (Math.abs(ldexp(mu[stage]![j]!, expoB[stage]! - expoB[j]!)) > eta) {
          goOn = true;
          break;
        }
      if (goOn)
        for (let j = stage - 1; j > zeros; j--) {
          const e = expoB[j]! - expoB[stage]!,
            tmp = ldexp(mu[stage]![j]!, -e),
            absolute = Math.abs(tmp);
          if (absolute <= 0.5) continue;
          didSomething = true;
          let coefficient: bigint,
            coefficientScale = 0n,
            scaled: number,
            scaleExponent = e;
          if (absolute <= 1.5) {
            coefficient = mu[stage]![j]! > 0 ? 1n : -1n;
            for (let k = zeros + 1; k < j; k++) {
              const term = ldexp(mu[j]![k]!, e);
              mu[stage]![k] = coefficient > 0n ? mu[stage]![k]! - term : mu[stage]![k]! + term;
            }
          } else {
            if (absolute < 2 ** 53) {
              scaled = rint(tmp);
              coefficient = BigInt(scaled);
            } else {
              const [fraction, exponent] = frexp(mu[stage]![j]!);
              let xx = signedWord(ldexp(fraction, 53)),
                E = exponent - e - 53;
              if (E <= 0) {
                xx = BigInt.asIntN(64, xx << BigInt(-E & 63));
                coefficient = xx;
                scaled = Number(xx);
              } else {
                coefficient = xx;
                coefficientScale = BigInt(E);
                scaled = Number(xx);
                scaleExponent = E + e;
              }
            }
            for (let k = zeros + 1; k < j; k++)
              mu[stage]![k] -= ldexp(scaled! * mu[j]![k]!, scaleExponent);
          }
          for (let row = 1; row <= n; row++)
            B[stage]![row] -= (coefficient * B[j]![row]!) << coefficientScale;
          if (U)
            for (let row = 1; row <= d; row++)
              U[stage]![row] -= (coefficient * U[j]![row]!) << coefficientScale;
        }
      if (!goOn) break;
      setLine(stage);
      setG(stage, zeros + 1, stage - 1);
      aa = zeros + 1;
    }
    if (didSomething) setG2(stage, stage, maxG);
    s[zeros + 1] = G[stage]![stage]!;
    for (let k = zeros + 1; k <= stage - 2; k++)
      s[k + 1] = fma(-mu[stage]![k]!, r[stage]![k]!, s[k]!);
    return false;
  };

  while (++kappa <= d) {
    if (kappa > kappamax) {
      maxG = kappamax = kappa;
      setG(kappa, zeros + 1, kappa);
    }
    if (babai(kappa, alpha[kappa]!)) return result(-1);
    let tmp = ldexp(r[kappa - 1]![kappa - 1]! * delta, 2 * (expoB[kappa - 1]! - expoB[kappa]!));
    if ((keepfirst && kappa === 2) || tmp <= s[kappa - 1]!) {
      alpha[kappa] = kappa;
      tmp = mu[kappa]![kappa - 1]! * r[kappa]![kappa - 1]!;
      r[kappa]![kappa] = s[kappa - 1]! - tmp;
      continue;
    }
    const kappa2 = kappa;
    do {
      kappa--;
      if (kappa < zeros + 2 + (keepfirst ? 1 : 0)) break;
      tmp = ldexp(r[kappa - 1]![kappa - 1]! * delta, 2 * (expoB[kappa - 1]! - expoB[kappa2]!));
    } while (s[kappa - 1]! <= tmp);
    updateAlpha(alpha, kappa, kappa2, kappamax);
    rotate(mu, kappa2, kappa);
    rotate(r, kappa2, kappa);
    r[kappa]![kappa] = s[kappa]!;
    rotate(B, kappa2, kappa);
    rotate(appB, kappa2, kappa);
    if (U) rotate(U, kappa2, kappa);
    rotate(expoB, kappa2, kappa);
    rotateGram(G, kappa2, kappa, maxG);
    if (kappa === zeros + 1 && G[kappa]![kappa]! <= 0) {
      zeros++;
      kappa++;
      G[kappa]![kappa] = dot(kappa, kappa, true);
      r[kappa]![kappa] = G[kappa]![kappa]!;
    }
  }
  return result(zeros);
}

/** A normalized double significand with PARI's signed-word exponent. */
type Dpe = { d: number; e: bigint };
const DPE_ZERO_EXPONENT = -((1n << 63n) - 1n);
function dpe_normalize(x: Dpe): Dpe {
  if (x.d === 0) return { d: x.d, e: DPE_ZERO_EXPONENT };
  const [d, e] = frexp(x.d);
  return { d, e: BigInt.asIntN(64, x.e + BigInt(e)) };
}
function affidpe(x: bigint): Dpe {
  const r = itor(x, 64);
  return { d: rtodbl(setexpo(r, -1)), e: BigInt(r.e + 1) };
}
function affdbldpe(x: number): Dpe {
  return dpe_normalize({ d: x, e: 0n });
}
function dpe_mulz(x: Dpe, y: Dpe): Dpe {
  return dpe_normalize({ d: x.d * y.d, e: BigInt.asIntN(64, x.e + y.e) });
}
function dpe_divz(x: Dpe, y: Dpe): Dpe {
  return dpe_normalize({ d: x.d / y.d, e: BigInt.asIntN(64, x.e - y.e) });
}
function dpe_addz(y: Dpe, z: Dpe): Dpe {
  if (y.e > z.e + 53n) return y;
  if (z.e > y.e + 53n) return z;
  const e = y.e - z.e;
  return dpe_normalize(
    e >= 0n
      ? { d: y.d + ldexp(z.d, -Number(e)), e: y.e }
      : { d: z.d + ldexp(y.d, Number(e)), e: z.e }
  );
}
function dpe_subz(y: Dpe, z: Dpe): Dpe {
  if (y.e > z.e + 53n) return y;
  if (z.e > y.e + 53n) return { d: -z.d, e: z.e };
  const e = y.e - z.e;
  return dpe_normalize(
    e >= 0n
      ? { d: y.d - ldexp(z.d, -Number(e)), e: y.e }
      : { d: ldexp(y.d, Number(e)) - z.d, e: z.e }
  );
}
function dpe_muluz(x: Dpe, u: bigint): Dpe {
  return dpe_normalize({ d: x.d * Number(u), e: x.e });
}
function dpe_submulz(x: Dpe, y: Dpe, z: Dpe): Dpe {
  return dpe_subz(x, dpe_mulz(y, z));
}
function dpe_cmp(x: Dpe, y: Dpe): number {
  const sx = x.d < 0 ? -1 : +(x.d > 0),
    sy = y.d < 0 ? -1 : +(y.d > 0);
  if (sx !== sy) return sx - sy;
  if (x.e > y.e) return sx > 0 ? 1 : -1;
  if (y.e > x.e) return sx > 0 ? -1 : 1;
  return x.d < y.d ? -1 : +(x.d > y.d);
}
function dpe_abscmp(x: Dpe, y: Dpe): number {
  if (x.e > y.e) return 1;
  if (y.e > x.e) return -1;
  return Math.abs(x.d) < Math.abs(y.d) ? -1 : +(Math.abs(x.d) > Math.abs(y.d));
}
function dpetor(x: Dpe): MpReal {
  const r = dbltor(x.d);
  return r.s ? setexpo(r, Number(x.e - 1n)) : r;
}

/** PARI's DPE reduction stage, with optional exact Gram matrix and returned norms.
 * Copies the supplied column matrices and returns [status,G,B,U,norms]. A null G
 * requests incremental construction; B may be null when a Gram matrix is supplied.
 * @see Deviation: PARI DPE LLL stage and resource boundary
 * @see reference/pari/src/basemath/lll.c:1430-2022
 */
export function fplll_dpe(
  input: Matrix | null,
  gram: Matrix | null = null,
  delta = 0.99,
  eta = 0.51,
  keepfirst = false,
  trackTransform = true,
  wantNorms = false
): [number, Matrix | null, Matrix | null, Matrix | null, MpReal[] | null] {
  const d = gram?.length ?? input?.length ?? 0,
    n = input?.[0]?.length ?? 0;
  if (
    !d ||
    (input && (!n || input.length !== d || input.some((c) => c.length !== n))) ||
    (gram && gram.some((c) => c.length !== d))
  )
    throw new RangeError('fplll_dpe requires a nonempty rectangular basis or square Gram matrix');
  const incgram = !gram;
  const B: Matrix | null = input ? [[], ...input.map((c) => [0n, ...c])] : null;
  const U: Matrix | null = trackTransform
    ? [
        [],
        ...Array.from({ length: d }, (_, j) => [
          0n,
          ...Array.from({ length: d }, (_, i) => (i === j ? 1n : 0n)),
        ]),
      ]
    : null;
  const G: Matrix = gram
    ? [[], ...gram.map((c) => [0n, ...c])]
    : Array.from({ length: d + 1 }, () => Array<bigint>(d + 1).fill(0n));
  const mat = () =>
    Array.from({ length: d + 1 }, () => Array.from({ length: d + 1 }, () => ({ d: 0, e: 0n })));
  const mu = mat(),
    r = mat(),
    s: Dpe[] = Array.from({ length: d + 1 }, () => ({ d: 0, e: 0n })),
    alpha = Array<number>(d + 1).fill(0),
    D = affdbldpe(delta),
    E = affdbldpe(eta);
  const strip = (A: Matrix | null): Matrix | null => (A ? A.slice(1).map((c) => c.slice(1)) : null);
  const result = (
    status: number
  ): [number, Matrix | null, Matrix | null, Matrix | null, MpReal[] | null] => [
    status,
    status < 0 && incgram ? null : strip(G),
    strip(B),
    strip(U),
    wantNorms && status >= 0
      ? Array.from({ length: d }, (_, i) => dpetor(r[i + 1]![i + 1]!))
      : null,
  ];
  const dot = (a: number, b: number): bigint => {
    let sum = 0n;
    for (let i = 1; i <= n; i++) sum += B![a]![i]! * B![b]![i]!;
    return sum;
  };
  let maxG = incgram ? 2 : d,
    kappamax = 1,
    i = 1;
  do {
    if (incgram) G[i]![i] = dot(i, i);
    r[i]![i] = affidpe(G[i]![i]!);
  } while (G[i]![i] === 0n && ++i <= d);
  let zeros = i - 1,
    kappa = i;
  for (i = zeros + 1; i <= d; i++) alpha[i] = 1;
  const babai = (stage: number, a: number): boolean => {
    let aa = a > zeros ? a : zeros + 1,
      emaxmu = -2n,
      emax2mu = -2n;
    for (;;) {
      let goOn = false;
      const emax3mu = emax2mu;
      emax2mu = emaxmu;
      emaxmu = -2n;
      for (let j = aa; j < stage; j++) {
        let g = affidpe(G[stage]![j]!);
        for (let k = zeros + 1; k < j; k++) g = dpe_submulz(g, mu[j]![k]!, r[stage]![k]!);
        r[stage]![j] = g;
        const x = dpe_divz(g, r[j]![j]!);
        mu[stage]![j] = x;
        if (x.e > emaxmu) emaxmu = x.e;
      }
      if (emax3mu !== -2n && emax3mu <= emax2mu + 5n) return true;
      for (let j = stage - 1; j > zeros; j--)
        if (dpe_abscmp(mu[stage]![j]!, E) > 0) {
          goOn = true;
          break;
        }
      if (goOn)
        for (let j = stage - 1; j > zeros; j--) {
          const tmp = mu[stage]![j]!;
          if (tmp.e < 0n) continue;
          const positive = tmp.d > 0,
            small = tmp.e <= 0n || (tmp.e === 1n && Math.abs(tmp.d) <= 0.75);
          let coefficient: bigint,
            scale = 0n;
          if (small) {
            coefficient = positive ? 1n : -1n;
            for (let k = zeros + 1; k < j; k++)
              mu[stage]![k] = positive
                ? dpe_subz(mu[stage]![k]!, mu[j]![k]!)
                : dpe_addz(mu[stage]![k]!, mu[j]![k]!);
          } else {
            const shift = tmp.e < 63n ? 0n : tmp.e - 63n;
            scale = shift;
            // A nonfinite GSO coefficient can carry an exponent near LONG_MAX.
            // Native temporary reservation then overflows its PARI stack, even for u=0.
            if (shift > BigInt(Number.MAX_SAFE_INTEGER))
              throw new RangeError(
                'fplll_dpe: Gram-Schmidt coefficient requires an unrepresentable shift'
              );
            const value = rint(ldexp(positive ? tmp.d : -tmp.d, Number(tmp.e - shift)));
            const u =
              Number.isNaN(value) || value <= 0
                ? 0n
                : value >= 2 ** 64
                  ? (1n << 64n) - 1n
                  : BigInt(value);
            coefficient = positive ? u : -u;
            for (let k = zeros + 1; k < j; k++) {
              const x = dpe_muluz(mu[j]![k]!, u);
              const scaled = { d: x.d, e: BigInt.asIntN(64, x.e + shift) };
              mu[stage]![k] = positive
                ? dpe_subz(mu[stage]![k]!, scaled)
                : dpe_addz(mu[stage]![k]!, scaled);
            }
          }
          // Native mului/sqru multiply by the word mantissa before shifting.
          // Expanding X = u*2^e first turns these into huge dense products.
          const product = (x: bigint): bigint => (coefficient * x) << scale;
          if (B) for (let row = 1; row <= n; row++) B[stage]![row] -= product(B[j]![row]!);
          if (U) for (let row = 1; row <= d; row++) U[stage]![row] -= product(U[j]![row]!);
          G[stage]![stage] +=
            ((coefficient * coefficient * G[j]![j]!) << (2n * scale)) -
            ((coefficient * G[stage]![j]!) << (scale + 1n));
          for (let row = 1; row <= j; row++) G[stage]![row] -= product(G[j]![row]!);
          for (let row = j + 1; row < stage; row++) G[stage]![row] -= product(G[row]![j]!);
          for (let row = stage + 1; row <= maxG; row++) G[row]![stage] -= product(G[row]![j]!);
        }
      if (!goOn) break;
      aa = zeros + 1;
    }
    s[zeros + 1] = affidpe(G[stage]![stage]!);
    for (let k = zeros + 1; k <= stage - 2; k++)
      s[k + 1] = dpe_submulz(s[k]!, mu[stage]![k]!, r[stage]![k]!);
    return false;
  };
  while (++kappa <= d) {
    if (kappa > kappamax) {
      kappamax = kappa;
      if (incgram) {
        for (i = zeros + 1; i <= kappa; i++) G[kappa]![i] = dot(kappa, i);
        maxG = kappamax;
      }
    }
    if (babai(kappa, alpha[kappa]!)) return result(-1);
    if (
      (keepfirst && kappa === 2) ||
      dpe_cmp(dpe_mulz(r[kappa - 1]![kappa - 1]!, D), s[kappa - 1]!) <= 0
    ) {
      alpha[kappa] = kappa;
      r[kappa]![kappa] = dpe_submulz(s[kappa - 1]!, mu[kappa]![kappa - 1]!, r[kappa]![kappa - 1]!);
      continue;
    }
    const kappa2 = kappa;
    do {
      kappa--;
      if (kappa < zeros + 2 + (keepfirst ? 1 : 0)) break;
    } while (dpe_cmp(dpe_mulz(r[kappa - 1]![kappa - 1]!, D), s[kappa - 1]!) >= 0);
    updateAlpha(alpha, kappa, kappa2, kappamax);
    rotate(mu, kappa2, kappa);
    rotate(r, kappa2, kappa);
    r[kappa]![kappa] = s[kappa]!;
    if (U) rotate(U, kappa2, kappa);
    if (B) rotate(B, kappa2, kappa);
    rotateGram(G, kappa2, kappa, maxG);
    if (kappa === zeros + 1 && G[kappa]![kappa] === 0n) {
      zeros++;
      kappa++;
      r[kappa]![kappa] = affidpe(G[kappa]![kappa]!);
    }
  }
  return result(zeros);
}

/** Native <= 1/2 predicate, including finite-accuracy zeros. */
function absrsmall(x: MpReal): boolean {
  if (!x.s || x.e < -1) return true;
  return x.e === -1 && x.m === 1n << BigInt(x.p - 1);
}
/** Preserve the native leading-word AND zero-tail test, even below 3/2. */
function absrsmall2(x: MpReal): boolean {
  if (x.e < 0) return true;
  if (x.e > 0) return false;
  const shift = BigInt(x.p - 64);
  return x.m >> shift <= 3n << 62n && (x.m & ((1n << shift) - 1n)) === 0n;
}
function truncexpo(x: MpReal, bits: number): [bigint, number] {
  const e = x.e + 1 - bits;
  return e >= 0 ? [mantissa2nr(x, 0), e] : [roundr_safe(x), 0];
}
/** Raw GEN boundary around the typed divrr dependency, preserving cypari2's payload. */
function lll_divrr(x: MpReal, y: MpReal): MpReal {
  if (!y.s) {
    // es.c real0tostr/ex10. LLL's assigned zero has e <= -working precision.
    const decimal = Math.trunc(y.e * 0.30102999566398119521 - (y.e < 0 ? 1 : 0)) + 1;
    throw new PariError(`impossible inverse in divrr: 0.E${decimal}`);
  }
  return divrr(x, y);
}

/** Native heuristic LLL stage with separate GSO and approximate Gram precision.
 * Precision arguments count bits and must be positive multiples of 64.
 * @see Deviation: PARI real LLL stage adapters
 * @see reference/pari/src/basemath/lll.c:1129-1425
 */
export function fplll_heuristic(
  input: Matrix,
  delta = 0.99,
  eta = 0.51,
  keepfirst = false,
  trackTransform = true,
  precision = 64,
  gramPrecision = precision
): [number, Matrix, Matrix | null] {
  const d = input.length,
    n = input[0]?.length ?? 0;
  if (!d || !n || input.some((c) => c.length !== n))
    throw new RangeError('fplll_heuristic requires a nonempty rectangular basis');
  if (![precision, gramPrecision].every((p) => Number.isSafeInteger(p) && p >= 64 && p % 64 === 0))
    throw new RangeError('LLL real precision must be a positive multiple of 64 bits');
  const B: Matrix = [[], ...input.map((c) => [0n, ...c])];
  const U: Matrix | null = trackTransform
    ? [
        [],
        ...Array.from({ length: d }, (_, j) => [
          0n,
          ...Array.from({ length: d }, (_, i) => (i === j ? 1n : 0n)),
        ]),
      ]
    : null;
  const assign = (x: MpReal, p = precision): MpReal => ({ ...rtor(x, p) });
  const mat = (p: number) =>
    Array.from({ length: d + 1 }, () => Array.from({ length: d + 1 }, () => itor(0n, p)));
  const G = mat(gramPrecision),
    mu = mat(precision),
    r = mat(precision),
    s = Array.from({ length: d + 1 }, () => itor(0n, precision)),
    alpha = Array<number>(d + 1).fill(0),
    appB = [
      [],
      ...input.map((c) => [itor(0n, gramPrecision), ...c.map((x) => itor(x, gramPrecision))]),
    ],
    D = dbltor(delta),
    E = dbltor(eta);
  const dot = (a: number, b: number): MpReal => {
    let sum = a === b ? sqrr(appB[a]![1]!) : mulrr(appB[a]![1]!, appB[b]![1]!);
    for (let row = 2; row <= n; row++)
      sum = addrr(sum, a === b ? sqrr(appB[a]![row]!) : mulrr(appB[a]![row]!, appB[b]![row]!));
    return sum;
  };
  const setG = (stage: number, a: number, b: number) => {
    for (let j = a; j <= b; j++) G[stage]![j] = assign(dot(stage, j), gramPrecision);
  };
  const setG2 = (stage: number, a: number, b: number) => {
    for (let j = a; j <= b; j++) G[j]![stage] = assign(dot(stage, j), gramPrecision);
  };
  const result = (status: number): [number, Matrix, Matrix | null] => [
    status,
    B.slice(1).map((c) => c.slice(1)),
    U ? U.slice(1).map((c) => c.slice(1)) : null,
  ];
  let kappamax = 1,
    maxG = d,
    i = 1;
  do {
    G[i]![i] = assign(dot(i, i), gramPrecision);
  } while (!G[i]![i]!.s && ++i <= d);
  let zeros = i - 1,
    kappa = i;
  if (zeros < d) r[zeros + 1]![zeros + 1] = assign(G[zeros + 1]![zeros + 1]!);
  for (i = zeros + 1; i <= d; i++) alpha[i] = 1;
  const babai = (stage: number, a: number): boolean => {
    let aa = a > zeros ? a : zeros + 1,
      emaxmu = -2,
      emax2mu = -2,
      didSomething = false;
    for (;;) {
      let goOn = false;
      const emax3mu = emax2mu;
      emax2mu = emaxmu;
      emaxmu = -2;
      for (let j = aa; j < stage; j++) {
        let g = G[stage]![j]!;
        for (let k = zeros + 1; k < j; k++) g = subrr(g, mulrr(mu[j]![k]!, r[stage]![k]!));
        r[stage]![j] = assign(g);
        mu[stage]![j] = assign(lll_divrr(r[stage]![j]!, r[j]![j]!));
        emaxmu = Math.max(emaxmu, mu[stage]![j]!.e);
      }
      if (emax3mu !== -2 && emax3mu <= emax2mu + 5) return true;
      for (let j = stage - 1; j > zeros; j--)
        if (abscmprr(mu[stage]![j]!, E) > 0) {
          goOn = true;
          break;
        }
      if (goOn)
        for (let j = stage - 1; j > zeros; j--) {
          const tmp = mu[stage]![j]!;
          if (absrsmall(tmp)) continue;
          didSomething = true;
          let coefficient: bigint,
            scale = 0n;
          if (absrsmall2(tmp)) {
            coefficient = tmp.s > 0 ? 1n : -1n;
            for (let k = zeros + 1; k < j; k++)
              mu[stage]![k] = assign(
                tmp.s > 0 ? subrr(mu[stage]![k]!, mu[j]![k]!) : addrr(mu[stage]![k]!, mu[j]![k]!)
              );
          } else if (tmp.e < 64) {
            const rounded = roundr_safe(tmp),
              u = BigInt.asUintN(64, rounded < 0n ? -rounded : rounded);
            coefficient = tmp.s > 0 ? u : -u;
            for (let k = zeros + 1; k < j; k++) {
              const product = mulir(u, mu[j]![k]!);
              mu[stage]![k] = assign(
                tmp.s > 0 ? subrr(mu[stage]![k]!, product) : addrr(mu[stage]![k]!, product)
              );
            }
          } else {
            const [X, e] = truncexpo(tmp, precision);
            coefficient = X;
            scale = BigInt(e);
            for (let k = zeros + 1; k < j; k++) {
              let product = mulir(X, mu[j]![k]!);
              if (e) product = shiftr(product, e);
              mu[stage]![k] = assign(subrr(mu[stage]![k]!, product));
            }
          }
          for (let row = 1; row <= n; row++) B[stage]![row] -= (coefficient * B[j]![row]!) << scale;
          if (U)
            for (let row = 1; row <= d; row++)
              U[stage]![row] -= (coefficient * U[j]![row]!) << scale;
        }
      if (!goOn) break;
      for (let row = 1; row <= n; row++) appB[stage]![row] = itor(B[stage]![row]!, gramPrecision);
      setG(stage, zeros + 1, stage - 1);
      aa = zeros + 1;
    }
    if (didSomething) setG2(stage, stage, maxG);
    s[zeros + 1] = assign(G[stage]![stage]!);
    for (let k = zeros + 1; k <= stage - 2; k++)
      s[k + 1] = assign(subrr(s[k]!, mulrr(mu[stage]![k]!, r[stage]![k]!)));
    return false;
  };
  while (++kappa <= d) {
    if (kappa > kappamax) {
      maxG = kappamax = kappa;
      setG(kappa, zeros + 1, kappa);
    }
    if (babai(kappa, alpha[kappa]!)) return result(-1);
    if (
      (keepfirst && kappa === 2) ||
      cmprr(mulrr(r[kappa - 1]![kappa - 1]!, D), s[kappa - 1]!) <= 0
    ) {
      alpha[kappa] = kappa;
      r[kappa]![kappa] = assign(
        subrr(s[kappa - 1]!, mulrr(mu[kappa]![kappa - 1]!, r[kappa]![kappa - 1]!))
      );
      continue;
    }
    const kappa2 = kappa;
    do {
      kappa--;
      if (kappa < zeros + 2 + (keepfirst ? 1 : 0)) break;
    } while (cmprr(s[kappa - 1]!, mulrr(r[kappa - 1]![kappa - 1]!, D)) <= 0);
    updateAlpha(alpha, kappa, kappa2, kappamax);
    rotate(mu, kappa2, kappa);
    rotate(r, kappa2, kappa);
    r[kappa]![kappa] = assign(s[kappa]!);
    rotate(B, kappa2, kappa);
    rotate(appB, kappa2, kappa);
    if (U) rotate(U, kappa2, kappa);
    rotateGram(G, kappa2, kappa, maxG);
    if (kappa === zeros + 1 && !G[kappa]![kappa]!.s) {
      zeros++;
      kappa++;
      G[kappa]![kappa] = assign(dot(kappa, kappa), gramPrecision);
      r[kappa]![kappa] = assign(G[kappa]![kappa]!);
    }
  }
  return result(zeros);
}

/** PARI arbitrary-precision proved LLL stage, with optional Gram/B/U/norm state.
 * Precision counts bits and must be a positive multiple of 64.
 * @see Deviation: PARI real LLL stage adapters
 * @see reference/pari/src/basemath/lll.c:2025-2314
 */
export function fplll(
  input: Matrix | null,
  gram: Matrix | null = null,
  delta = 0.99,
  eta = 0.51,
  keepfirst = false,
  trackTransform = true,
  wantNorms = false,
  precision = 64
): [number, Matrix | null, Matrix | null, Matrix | null, MpReal[] | null] {
  const d = gram?.length ?? input?.length ?? 0,
    n = input?.[0]?.length ?? 0;
  if (
    !d ||
    (input && (!n || input.length !== d || input.some((c) => c.length !== n))) ||
    (gram && gram.some((c) => c.length !== d))
  )
    throw new RangeError('fplll requires a nonempty rectangular basis or square Gram matrix');
  if (!Number.isSafeInteger(precision) || precision < 64 || precision % 64 !== 0)
    throw new RangeError('LLL real precision must be a positive multiple of 64 bits');
  const incgram = !gram;
  const B: Matrix | null = input ? [[], ...input.map((c) => [0n, ...c])] : null;
  const U: Matrix | null = trackTransform
    ? [
        [],
        ...Array.from({ length: d }, (_, j) => [
          0n,
          ...Array.from({ length: d }, (_, i) => (i === j ? 1n : 0n)),
        ]),
      ]
    : null;
  const G: Matrix = gram
    ? [[], ...gram.map((c) => [0n, ...c])]
    : Array.from({ length: d + 1 }, () => Array<bigint>(d + 1).fill(0n));
  const assign = (x: MpReal): MpReal => ({ ...rtor(x, precision) });
  const mat = () =>
    Array.from({ length: d + 1 }, () => Array.from({ length: d + 1 }, () => itor(0n, precision)));
  const mu = mat(),
    r = mat(),
    s = Array.from({ length: d + 1 }, () => itor(0n, precision)),
    alpha = Array<number>(d + 1).fill(0),
    D = dbltor(delta),
    E = dbltor(eta);
  const strip = (A: Matrix | null): Matrix | null => (A ? A.slice(1).map((c) => c.slice(1)) : null);
  const result = (
    status: number
  ): [number, Matrix | null, Matrix | null, Matrix | null, MpReal[] | null] => [
    status,
    status < 0 && incgram ? null : strip(G),
    strip(B),
    strip(U),
    wantNorms && status >= 0 ? Array.from({ length: d }, (_, i) => r[i + 1]![i + 1]!) : null,
  ];
  const dot = (a: number, b: number): bigint => {
    let sum = 0n;
    for (let i = 1; i <= n; i++) sum += B![a]![i]! * B![b]![i]!;
    return sum;
  };
  let maxG = incgram ? 2 : d,
    kappamax = 1,
    i = 1;
  do {
    if (incgram) G[i]![i] = dot(i, i);
    r[i]![i] = itor(G[i]![i]!, precision);
  } while (G[i]![i] === 0n && ++i <= d);
  let zeros = i - 1,
    kappa = i;
  for (i = zeros + 1; i <= d; i++) alpha[i] = 1;
  const babai = (stage: number, a: number): boolean => {
    let aa = a > zeros ? a : zeros + 1,
      emaxmu = -2,
      emax2mu = -2;
    for (;;) {
      let goOn = false;
      const emax3mu = emax2mu;
      emax2mu = emaxmu;
      emaxmu = -2;
      for (let j = aa; j < stage; j++) {
        let g: bigint | MpReal = G[stage]![j]!;
        for (let k = zeros + 1; k < j; k++) {
          const product = mulrr(mu[j]![k]!, r[stage]![k]!);
          g = typeof g === 'bigint' ? subir(g, product) : subrr(g, product);
        }
        r[stage]![j] = typeof g === 'bigint' ? itor(g, precision) : assign(g);
        mu[stage]![j] = assign(lll_divrr(r[stage]![j]!, r[j]![j]!));
        emaxmu = Math.max(emaxmu, mu[stage]![j]!.e);
      }
      if (emax3mu !== -2 && emax3mu <= emax2mu + 5) return true;
      for (let j = stage - 1; j > zeros; j--)
        if (abscmprr(mu[stage]![j]!, E) > 0) {
          goOn = true;
          break;
        }
      if (goOn)
        for (let j = stage - 1; j > zeros; j--) {
          const tmp = mu[stage]![j]!;
          if (absrsmall(tmp)) continue;
          let coefficient: bigint,
            scale = 0n;
          if (absrsmall2(tmp)) {
            coefficient = tmp.s > 0 ? 1n : -1n;
            for (let k = zeros + 1; k < j; k++)
              mu[stage]![k] = assign(
                tmp.s > 0 ? subrr(mu[stage]![k]!, mu[j]![k]!) : addrr(mu[stage]![k]!, mu[j]![k]!)
              );
          } else if (tmp.e < 64) {
            const rounded = roundr_safe(tmp),
              u = BigInt.asUintN(64, rounded < 0n ? -rounded : rounded);
            coefficient = tmp.s > 0 ? u : -u;
            for (let k = zeros + 1; k < j; k++) {
              const product = mulir(u, mu[j]![k]!);
              mu[stage]![k] = assign(
                tmp.s > 0 ? subrr(mu[stage]![k]!, product) : addrr(mu[stage]![k]!, product)
              );
            }
          } else {
            const [X, e] = truncexpo(tmp, precision);
            coefficient = X;
            scale = BigInt(e);
            for (let k = zeros + 1; k < j; k++) {
              let product = mulir(X, mu[j]![k]!);
              if (e) product = shiftr(product, e);
              mu[stage]![k] = assign(subrr(mu[stage]![k]!, product));
            }
          }
          if (B)
            for (let row = 1; row <= n; row++)
              B[stage]![row] -= (coefficient * B[j]![row]!) << scale;
          if (U)
            for (let row = 1; row <= d; row++)
              U[stage]![row] -= (coefficient * U[j]![row]!) << scale;
          G[stage]![stage] +=
            ((coefficient * coefficient * G[j]![j]!) << (2n * scale)) -
            ((coefficient * G[stage]![j]!) << (scale + 1n));
          for (let row = 1; row <= j; row++) G[stage]![row] -= (coefficient * G[j]![row]!) << scale;
          for (let row = j + 1; row < stage; row++)
            G[stage]![row] -= (coefficient * G[row]![j]!) << scale;
          for (let row = stage + 1; row <= maxG; row++)
            G[row]![stage] -= (coefficient * G[row]![j]!) << scale;
        }
      if (!goOn) break;
      aa = zeros + 1;
    }
    s[zeros + 1] = itor(G[stage]![stage]!, precision);
    for (let k = zeros + 1; k <= stage - 2; k++)
      s[k + 1] = assign(subrr(s[k]!, mulrr(mu[stage]![k]!, r[stage]![k]!)));
    return false;
  };
  while (++kappa <= d) {
    if (kappa > kappamax) {
      kappamax = kappa;
      if (incgram) {
        for (i = zeros + 1; i <= kappa; i++) G[kappa]![i] = dot(kappa, i);
        maxG = kappamax;
      }
    }
    if (babai(kappa, alpha[kappa]!)) return result(-1);
    if (
      (keepfirst && kappa === 2) ||
      cmprr(mulrr(r[kappa - 1]![kappa - 1]!, D), s[kappa - 1]!) <= 0
    ) {
      alpha[kappa] = kappa;
      r[kappa]![kappa] = assign(
        subrr(s[kappa - 1]!, mulrr(mu[kappa]![kappa - 1]!, r[kappa]![kappa - 1]!))
      );
      continue;
    }
    const kappa2 = kappa;
    do {
      kappa--;
      if (kappa < zeros + 2 + (keepfirst ? 1 : 0)) break;
    } while (cmprr(s[kappa - 1]!, mulrr(r[kappa - 1]![kappa - 1]!, D)) <= 0);
    updateAlpha(alpha, kappa, kappa2, kappamax);
    rotate(mu, kappa2, kappa);
    rotate(r, kappa2, kappa);
    r[kappa]![kappa] = assign(s[kappa]!);
    if (U) rotate(U, kappa2, kappa);
    if (B) rotate(B, kappa2, kappa);
    rotateGram(G, kappa2, kappa, maxG);
    if (kappa === zeros + 1 && G[kappa]![kappa] === 0n) {
      zeros++;
      kappa++;
      r[kappa]![kappa] = itor(G[kappa]![kappa]!, precision);
    }
  }
  return result(zeros);
}

export {
  drop,
  potential,
  spread,
  condition_bound,
  GS_extraprec,
  gramschmidt_upper,
  gramschmidt_dynprec,
  RgM_Cholesky_dynprec,
} from './_lll_gso.js';

export {
  flat,
  ZM_flatter,
  ZM_flatter_rank,
  ZM_flattergram,
  ZM_lll,
  ZM_lll_norms,
  lllfp,
  LLL_KER,
  LLL_IM,
  LLL_ALL,
  LLL_GRAM,
  LLL_KEEP_FIRST,
  LLL_INPLACE,
  LLL_COMPATIBLE,
  LLL_UPPER,
  LLL_NOCERTIFY,
  LLL_NOFLATTER,
  type LllResult,
} from './_lll_wrapper.js';
