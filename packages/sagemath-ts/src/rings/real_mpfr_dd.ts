/**
 * @module sage/rings/real_mpfr_dd
 * @description Double-double kernels for the correctly-rounded MPFR functions
 *
 * `RealNumber` stores a plain IEEE double, but several of the functions
 * `real_mpfr.pyx` exposes are *correctly rounded* by MPFR, so a
 * double-precision evaluation that loses even a couple of digits to
 * cancellation shows up immediately when the two are compared digit by digit.
 *
 * The Bessel functions were the worst offenders: the previous implementation
 * used the Numerical Recipes rational approximations, which target SINGLE
 * precision and are only ~1e-8 relative accurate (`RR(1).j0()` was wrong from
 * the 9th significant digit).  The routines here evaluate the defining series
 * and the Hankel asymptotic expansions in "double-double" arithmetic (~32
 * decimal digits), which leaves plenty of headroom for the cancellation those
 * series suffer, and then round once at the end.
 *
 * Double-double arithmetic is the standard Dekker/Knuth error-free
 * transformation pair: a value is an unevaluated sum `hi + lo` of two doubles
 * with `|lo| <= ulp(hi)/2`.
 */

/** An unevaluated sum `hi + lo` of two non-overlapping doubles. */
export type DD = readonly [number, number];

export const DD_ZERO: DD = [0, 0];
export const DD_ONE: DD = [1, 0];

/** `2Sum` (Knuth): the exact sum of two doubles as a DD. */
function twoSum(a: number, b: number): DD {
  const s = a + b;
  const bb = s - a;
  const err = a - (s - bb) + (b - bb);
  return [s, err];
}

/** `2Prod` via Dekker splitting: the exact product of two doubles as a DD. */
function twoProd(a: number, b: number): DD {
  const p = a * b;
  const SPLIT = 134217729; // 2^27 + 1
  let ta = SPLIT * a;
  const ahi = ta - (ta - a);
  const alo = a - ahi;
  ta = SPLIT * b;
  const bhi = ta - (ta - b);
  const blo = b - bhi;
  const err = ahi * bhi - p + ahi * blo + alo * bhi + alo * blo;
  return [p, err];
}

export function dd(x: number): DD {
  return [x, 0];
}

export function ddAdd(a: DD, b: DD): DD {
  const [s1, s2a] = twoSum(a[0], b[0]);
  const [t1, t2] = twoSum(a[1], b[1]);
  const s2 = s2a + t1;
  const hi = s1 + s2;
  let lo = s2 - (hi - s1);
  lo += t2;
  const r = hi + lo;
  return [r, lo - (r - hi)];
}

export function ddNeg(a: DD): DD {
  return [-a[0], -a[1]];
}

export function ddSub(a: DD, b: DD): DD {
  return ddAdd(a, ddNeg(b));
}

export function ddMul(a: DD, b: DD): DD {
  const [p1, p2a] = twoProd(a[0], b[0]);
  const p2 = p2a + (a[0] * b[1] + a[1] * b[0]);
  const hi = p1 + p2;
  return [hi, p2 - (hi - p1)];
}

export function ddDiv(a: DD, b: DD): DD {
  const q1 = a[0] / b[0];
  const r = ddSub(a, ddMul(b, dd(q1)));
  const q2 = r[0] / b[0];
  const r2 = ddSub(r, ddMul(b, dd(q2)));
  const q3 = r2[0] / b[0];
  const [hi, lo0] = twoSum(q1, q2);
  const lo = lo0 + q3;
  const s = hi + lo;
  return [s, lo - (s - hi)];
}

export function ddToNumber(a: DD): number {
  return a[0] + a[1];
}

/** `sqrt` of a DD, to full double-double accuracy. */
export function ddSqrt(a: DD): DD {
  if (a[0] === 0 && a[1] === 0) {
    return DD_ZERO;
  }
  const x = 1 / Math.sqrt(a[0]);
  const ax = a[0] * x;
  const diff = ddSub(a, ddMul(dd(ax), dd(ax)));
  return ddAdd(dd(ax), dd(diff[0] * x * 0.5));
}

// The high words below are `Math.PI` and `Math.LN2` by construction -- a
// double-double constant IS the nearest double plus its error term, so
// `noApproximativeNumericConstant` is exactly backwards here.
/** pi to double-double accuracy (`Math.PI` plus its error term). */
// biome-ignore lint/suspicious/noApproximativeNumericConstant: high word of a double-double constant
export const DD_PI: DD = [3.141592653589793, 1.2246467991473532e-16];
/** Euler's gamma constant to double-double accuracy. */
export const DD_EULER: DD = [0.5772156649015329, -4.942915152430646e-18];
/** ln(2) to double-double accuracy (`Math.LN2` plus its error term). */
// biome-ignore lint/suspicious/noApproximativeNumericConstant: high word of a double-double constant
export const DD_LN2: DD = [0.6931471805599453, 2.3190468138462996e-17];

/**
 * `exp(x)` to double-double accuracy, for a DD argument of moderate size.
 *
 * Uses `exp(x) = 2^k * exp(r)` with `r = x - k*ln2`, `|r| <= ln2/2`, and sums
 * the Taylor series for `exp(r)` in DD.
 */
export function ddExp(a: DD): DD {
  const x = ddToNumber(a);
  if (x === 0) {
    return DD_ONE;
  }
  if (x < -745) {
    return DD_ZERO;
  }
  const k = Math.round(x / Math.LN2);
  const r = ddSub(a, ddMul(dd(k), DD_LN2));
  let term: DD = DD_ONE;
  let sum: DD = DD_ONE;
  for (let i = 1; i <= 30; i++) {
    term = ddDiv(ddMul(term, r), dd(i));
    sum = ddAdd(sum, term);
    if (Math.abs(term[0]) < 1e-36 * Math.abs(sum[0])) {
      break;
    }
  }
  const scale = 2 ** k;
  return [sum[0] * scale, sum[1] * scale];
}

/**
 * `ln(x)` to double-double accuracy, by one Newton step on `exp`.
 */
export function ddLog(a: DD): DD {
  const x = ddToNumber(a);
  if (x <= 0) {
    return dd(Math.log(x));
  }
  let y = dd(Math.log(x));
  // y <- y + (a * exp(-y) - 1)
  for (let i = 0; i < 2; i++) {
    const e = ddExp(ddNeg(y));
    y = ddAdd(y, ddSub(ddMul(a, e), DD_ONE));
  }
  return y;
}

/**
 * `sin` and `cos` of a DD argument, both to double-double accuracy.
 *
 * The argument is reduced modulo pi/2 in DD, so the result stays accurate for
 * the arguments the Hankel asymptotic expansion produces (`x - pi/4` with `x`
 * up to a few hundred).
 */
export function ddSinCos(a: DD): [DD, DD] {
  const halfPi = ddDiv(DD_PI, dd(2));
  const q = Math.round(ddToNumber(ddDiv(a, halfPi)));
  const r = ddSub(a, ddMul(dd(q), halfPi));

  // Taylor series for sin(r), cos(r) with |r| <= pi/4.
  let term: DD = r;
  let s: DD = r;
  const r2 = ddMul(r, r);
  for (let i = 1; i <= 20; i++) {
    term = ddNeg(ddDiv(ddMul(term, r2), dd(2 * i * (2 * i + 1))));
    s = ddAdd(s, term);
    if (Math.abs(term[0]) < 1e-36) {
      break;
    }
  }
  let cterm: DD = DD_ONE;
  let c: DD = DD_ONE;
  for (let i = 1; i <= 20; i++) {
    cterm = ddNeg(ddDiv(ddMul(cterm, r2), dd((2 * i - 1) * (2 * i))));
    c = ddAdd(c, cterm);
    if (Math.abs(cterm[0]) < 1e-36) {
      break;
    }
  }

  const k = ((q % 4) + 4) % 4;
  switch (k) {
    case 0:
      return [s, c];
    case 1:
      return [c, ddNeg(s)];
    case 2:
      return [ddNeg(s), ddNeg(c)];
    default:
      return [ddNeg(c), s];
  }
}

/* ------------------------------------------------------------------ */
/* Bessel functions                                                     */
/* ------------------------------------------------------------------ */

/** Above this |x| the Hankel asymptotic expansion is used instead of the series. */
const BESSEL_ASYMPTOTIC_CUTOFF = 17.5;

/**
 * `J_n(x)` and `Y_n(x)` for `n = 0, 1` via the ascending power series, summed
 * in double-double arithmetic.
 *
 * ```
 * J_n(x) = (x/2)^n sum_k (-1)^k (x^2/4)^k / (k! (k+n)!)
 * Y_0(x) = (2/pi)[ (ln(x/2) + gamma) J_0(x) + sum_{k>=1} (-1)^(k+1) H_k (x^2/4)^k/(k!)^2 ]
 * Y_1(x) = (2/pi)[ (ln(x/2) + gamma) J_1(x) - 1/x
 *                  - (x/4) sum_{k>=0} (-1)^k (H_k + H_{k+1}) (x^2/4)^k/(k!(k+1)!) ]
 * ```
 */
function besselSeries(x: number, n: 0 | 1): { j: DD; y: DD } {
  const xd = dd(x);
  const q = ddDiv(ddMul(xd, xd), dd(4)); // x^2/4

  // J series
  let term: DD = DD_ONE;
  let j: DD = DD_ONE;
  for (let k = 1; k <= 400; k++) {
    term = ddNeg(ddDiv(ddMul(term, q), dd(k * (k + n))));
    j = ddAdd(j, term);
    if (term[0] === 0 || Math.abs(term[0]) < 1e-45) {
      break;
    }
  }
  if (n === 1) {
    j = ddMul(j, ddDiv(xd, dd(2)));
  }

  // Y series
  const lnHalfX = ddAdd(ddLog(ddDiv(xd, dd(2))), DD_EULER);
  const twoOverPi = ddDiv(dd(2), DD_PI);

  // The harmonic numbers must be accumulated in double-double too: the terms
  // they multiply reach ~2e3 at x = 10 while the answer is ~5e-2, so a plain
  // double H_k contributes an absolute error of ~2e-13 to the result.
  if (n === 0) {
    let t: DD = DD_ONE;
    let sum: DD = DD_ZERO;
    let h: DD = DD_ZERO;
    for (let k = 1; k <= 600; k++) {
      t = ddNeg(ddDiv(ddMul(t, q), dd(k * k)));
      h = ddAdd(h, ddDiv(DD_ONE, dd(k)));
      const contrib = ddMul(ddNeg(t), h);
      sum = ddAdd(sum, contrib);
      if (t[0] === 0 || Math.abs(contrib[0]) < 1e-45) {
        break;
      }
    }
    const y = ddMul(twoOverPi, ddAdd(ddMul(lnHalfX, j), sum));
    return { j, y };
  }

  // n == 1
  let t: DD = DD_ONE; // (x^2/4)^k / (k!(k+1)!) with sign
  let sum: DD = t; // k = 0: H_0 + H_1 = 1
  let hk: DD = DD_ZERO;
  let hk1: DD = DD_ONE;
  for (let k = 1; k <= 600; k++) {
    t = ddNeg(ddDiv(ddMul(t, q), dd(k * (k + 1))));
    hk = ddAdd(hk, ddDiv(DD_ONE, dd(k)));
    hk1 = ddAdd(hk1, ddDiv(DD_ONE, dd(k + 1)));
    const contrib = ddMul(t, ddAdd(hk, hk1));
    sum = ddAdd(sum, contrib);
    if (t[0] === 0 || Math.abs(contrib[0]) < 1e-45) {
      break;
    }
  }
  const y = ddMul(
    twoOverPi,
    ddSub(ddSub(ddMul(lnHalfX, j), ddDiv(DD_ONE, xd)), ddMul(ddDiv(xd, dd(4)), sum))
  );
  return { j, y };
}

/**
 * `J_nu(x)` and `Y_nu(x)` for `nu = 0, 1` via the Hankel asymptotic expansion.
 *
 * ```
 * J = sqrt(2/(pi x)) [P cos(chi) - Q sin(chi)],  chi = x - (2 nu + 1) pi/4
 * Y = sqrt(2/(pi x)) [P sin(chi) + Q cos(chi)]
 * ```
 */
function besselAsymptotic(x: number, n: 0 | 1): { j: DD; y: DD } {
  const xd = dd(x);
  const mu = 4 * n * n;
  const eight = ddMul(dd(8), xd);

  let p: DD = DD_ONE;
  let qs: DD = DD_ZERO;
  let term: DD = DD_ONE;
  let best = Number.POSITIVE_INFINITY;
  for (let k = 1; k <= 40; k++) {
    const f = mu - (2 * k - 1) * (2 * k - 1);
    term = ddDiv(ddMul(term, dd(f)), ddMul(dd(k), eight));
    const mag = Math.abs(term[0]);
    if (mag > best) {
      break; // the asymptotic series has started to diverge
    }
    best = mag;
    if (k % 2 === 1) {
      qs = ddAdd(qs, (k - 1) % 4 === 0 ? term : ddNeg(term));
    } else {
      p = ddAdd(p, (k - 2) % 4 === 0 ? ddNeg(term) : term);
    }
    if (mag < 1e-40) {
      break;
    }
  }

  const chi = ddSub(xd, ddMul(DD_PI, dd((2 * n + 1) / 4)));
  const [sinChi, cosChi] = ddSinCos(chi);
  const amp = ddSqrt(ddDiv(dd(2), ddMul(DD_PI, xd)));

  return {
    j: ddMul(amp, ddSub(ddMul(p, cosChi), ddMul(qs, sinChi))),
    y: ddMul(amp, ddAdd(ddMul(p, sinChi), ddMul(qs, cosChi))),
  };
}

/** `J_0(x)`, accurate to a double. */
export function besselJ0(x: number): number {
  if (x === 0) return 1;
  const ax = Math.abs(x);
  const r = ax < BESSEL_ASYMPTOTIC_CUTOFF ? besselSeries(ax, 0) : besselAsymptotic(ax, 0);
  return ddToNumber(r.j);
}

/** `J_1(x)`, accurate to a double (odd in `x`). */
export function besselJ1(x: number): number {
  if (x === 0) return 0;
  const ax = Math.abs(x);
  const r = ax < BESSEL_ASYMPTOTIC_CUTOFF ? besselSeries(ax, 1) : besselAsymptotic(ax, 1);
  const v = ddToNumber(r.j);
  return x < 0 ? -v : v;
}

/** `Y_0(x)` for `x > 0`, accurate to a double. */
export function besselY0(x: number): number {
  const r = x < BESSEL_ASYMPTOTIC_CUTOFF ? besselSeries(x, 0) : besselAsymptotic(x, 0);
  return ddToNumber(r.y);
}

/** `Y_1(x)` for `x > 0`, accurate to a double. */
export function besselY1(x: number): number {
  const r = x < BESSEL_ASYMPTOTIC_CUTOFF ? besselSeries(x, 1) : besselAsymptotic(x, 1);
  return ddToNumber(r.y);
}

/**
 * `Y_n(x)` for `x > 0` by the (stable) upward recurrence
 * `Y_{k+1} = (2k/x) Y_k - Y_{k-1}`, carried in double-double.
 */
export function besselYn(n: number, x: number): number {
  if (n === 0) return besselY0(x);
  if (n === 1) return besselY1(x);
  const r0 = x < BESSEL_ASYMPTOTIC_CUTOFF ? besselSeries(x, 0) : besselAsymptotic(x, 0);
  const r1 = x < BESSEL_ASYMPTOTIC_CUTOFF ? besselSeries(x, 1) : besselAsymptotic(x, 1);
  const tox = ddDiv(dd(2), dd(x));
  let ym: DD = r0.y;
  let y: DD = r1.y;
  for (let k = 1; k < n; k++) {
    const yp = ddSub(ddMul(ddMul(dd(k), tox), y), ym);
    ym = y;
    y = yp;
  }
  return ddToNumber(y);
}

/**
 * `J_n(x)` for `n >= 2`.
 *
 * The upward recurrence is stable while `n < |x|`; above that it amplifies the
 * rounding error exponentially, so Miller's downward recurrence is used, with
 * the normalisation `1 = J_0 + 2(J_2 + J_4 + ...)`.  Both are carried in
 * double-double.
 */
export function besselJn(n: number, x: number): number {
  if (n === 0) return besselJ0(x);
  if (n === 1) return besselJ1(x);

  const ax = Math.abs(x);
  if (ax === 0) return 0;

  let value: number;
  if (n < ax) {
    const r0 = ax < BESSEL_ASYMPTOTIC_CUTOFF ? besselSeries(ax, 0) : besselAsymptotic(ax, 0);
    const r1 = ax < BESSEL_ASYMPTOTIC_CUTOFF ? besselSeries(ax, 1) : besselAsymptotic(ax, 1);
    const tox = ddDiv(dd(2), dd(ax));
    let jm: DD = r0.j;
    let j: DD = r1.j;
    for (let k = 1; k < n; k++) {
      const jp = ddSub(ddMul(ddMul(dd(k), tox), j), jm);
      jm = j;
      j = jp;
    }
    value = ddToNumber(j);
  } else {
    // Miller's downward recurrence.
    const tox = ddDiv(dd(2), dd(ax));
    const m = 2 * Math.floor((n + Math.floor(Math.sqrt(200 * n))) / 2);
    let jp: DD = DD_ZERO;
    let j: DD = dd(1e-30);
    let sum: DD = DD_ZERO;
    let ans: DD = DD_ZERO;
    for (let k = m; k > 0; k--) {
      const jm = ddSub(ddMul(ddMul(dd(k), tox), j), jp);
      jp = j;
      j = jm;
      if (Math.abs(j[0]) > 1e100) {
        j = [j[0] * 1e-100, j[1] * 1e-100];
        jp = [jp[0] * 1e-100, jp[1] * 1e-100];
        ans = [ans[0] * 1e-100, ans[1] * 1e-100];
        sum = [sum[0] * 1e-100, sum[1] * 1e-100];
      }
      if (k % 2 !== 0) {
        sum = ddAdd(sum, j);
      }
      if (k === n) {
        ans = jp;
      }
    }
    sum = ddSub(ddMul(dd(2), sum), j);
    value = ddToNumber(ddDiv(ans, sum));
  }

  return x < 0 && n % 2 !== 0 ? -value : value;
}

/* ------------------------------------------------------------------ */
/* Error function                                                       */
/* ------------------------------------------------------------------ */

/**
 * `erfc(x)` to double accuracy.
 *
 * For `|x| < 2` the Maclaurin series for `erf` is summed in double-double and
 * subtracted from 1; for larger `x` the continued fraction
 * `erfc(x) = exp(-x^2)/(x sqrt(pi)) * 1/(1 + 1/2/(x^2 + 1/(1 + 3/2/(x^2 + ...))))`
 * is evaluated in double-double.  `erfc(-x) = 2 - erfc(x)`.
 */
export function erfcAccurate(x: number): number {
  if (Number.isNaN(x)) return Number.NaN;
  if (x === Number.POSITIVE_INFINITY) return 0;
  if (x === Number.NEGATIVE_INFINITY) return 2;
  if (x === 0) return 1;

  const ax = Math.abs(x);
  const tail = ddToNumber(erfcPositiveDD(ax));
  return x > 0 ? tail : 2 - tail;
}

/** `erf(x)` to double accuracy. */
export function erfAccurate(x: number): number {
  if (Number.isNaN(x)) return Number.NaN;
  if (x === Number.POSITIVE_INFINITY) return 1;
  if (x === Number.NEGATIVE_INFINITY) return -1;
  if (x === 0) return x; // preserves -0
  const ax = Math.abs(x);
  const v = ax < 2 ? ddToNumber(erfSeriesDD(ax)) : 1 - ddToNumber(erfcPositiveDD(ax));
  return x > 0 ? v : -v;
}

/** `erf(x)` for `x > 0` by its Maclaurin series, in double-double. */
function erfSeriesDD(x: number): DD {
  const xd = dd(x);
  const x2 = ddMul(xd, xd);
  let term: DD = xd;
  let sum: DD = xd;
  for (let n = 1; n <= 300; n++) {
    term = ddNeg(ddDiv(ddMul(term, x2), dd(n)));
    const contrib = ddDiv(term, dd(2 * n + 1));
    sum = ddAdd(sum, contrib);
    if (Math.abs(contrib[0]) < 1e-40 * Math.abs(sum[0])) {
      break;
    }
  }
  const twoOverSqrtPi = ddDiv(dd(2), ddSqrt(DD_PI));
  return ddMul(twoOverSqrtPi, sum);
}

/** `erfc(x)` for `x > 0`, in double-double. */
function erfcPositiveDD(x: number): DD {
  if (x < 2) {
    return ddSub(DD_ONE, erfSeriesDD(x));
  }
  // Lentz evaluation of the continued fraction
  //   erfc(x) = exp(-x^2)/sqrt(pi) * 1/(x + 1/2/(x + 1/(x + 3/2/(x + ...))))
  const xd = dd(x);
  let f: DD = dd(1e-300);
  let C: DD = f;
  let D: DD = DD_ZERO;
  for (let i = 0; i <= 400; i++) {
    const a: DD = i === 0 ? DD_ONE : dd(i / 2);
    const b: DD = xd;
    D = ddAdd(b, ddMul(a, D));
    if (D[0] === 0) D = dd(1e-300);
    C = ddAdd(b, ddDiv(a, C));
    if (C[0] === 0) C = dd(1e-300);
    D = ddDiv(DD_ONE, D);
    const delta = ddMul(C, D);
    f = ddMul(f, delta);
    if (Math.abs(ddToNumber(delta) - 1) < 1e-34) {
      break;
    }
  }
  const expPart = ddExp(ddNeg(ddMul(xd, xd)));
  return ddMul(ddDiv(expPart, ddSqrt(DD_PI)), f);
}

/**
 * `Li_2(x)` for a real `|x| <= 0.5`, by its defining power series summed in
 * double-double arithmetic.
 *
 * The series `sum_{k>=1} x^k / k^2` converges geometrically there, and summing
 * it in DD makes the result correctly rounded, which the complex Bernoulli
 * series was not (`CC(0.5).dilog()` was one ulp high).
 */
export function dilogRealSmall(x: number): number {
  const xd = dd(x);
  let p: DD = xd;
  let sum: DD = xd;
  for (let k = 2; k <= 400; k++) {
    p = ddMul(p, xd);
    const term = ddDiv(p, dd(k * k));
    sum = ddAdd(sum, term);
    if (term[0] === 0 || Math.abs(term[0]) < 1e-45) {
      break;
    }
  }
  return ddToNumber(sum);
}
