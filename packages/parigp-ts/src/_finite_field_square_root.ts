/** FF_issquareall -> FF_ispower(x,2), with the native root output requested.
 * References: FF.c:913; Flx.c:3546; FpX.c:2746; F2x.c:1700.
 * @see Deviation: PARI finite-field scalar square-root adapter
 */
import type { PariFfelt } from './types.js';
import { random_Flx, Flx_sqr } from './Flx.js';
import { random_FpX, FpX_sqr } from './FpX.js';
import { FpX_add, FpX_Fp_mul } from './ffinit.js';
import { Fp_sqrt, Fl_sqrt } from './ff.js';
import { polynomialQuotient } from './_polynomial_quotient.js';
import { polynomialQuotientInverse } from './_polynomial_quotient_power.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { gen_pow_i, gen_pow_fold } from './bb_group.js';
import { brent_kung_optpow } from './RgX.js';
import { F2x_sqr, F2x_rem, F2x_mul } from './F2x.js';

export function finiteFieldSquareRoot(x: PariFfelt): PariFfelt | null {
  const a = trimPolynomial(typeof x.value === 'bigint' ? [x.value] : [...x.value]);
  const T = [...(x.definingPoly ?? [0n, 1n])],
    p = x.p,
    d = T.length - 1;
  const result = (value: bigint[]): PariFfelt => ({
    ...x,
    value,
    definingPoly: x.definingPoly?.slice(),
  });
  if (!a.length) return result(a);
  // Native Flxq_sqrtn_spec_pre uses gequal1 on a t_VECSMALL, which is false
  // even for polynomial one. FpXQ uses a t_POL; binary/linear Shanks uses its
  // field-specific equal1 hook. Keep the word extension's random draw.
  if (a.length === 1 && a[0] === 1n && (p === 2n || p >= 1n << 64n || d === 1)) return result(a);
  if (p === 2n) {
    // gen_Shanks_sqrtn: gcd(2,2^d-1)=1, so exponent 2^{-1} mod (2^d-1).
    const pack = (v: bigint[]) => v.reduce((n, c, i) => n | (c << BigInt(i)), 0n);
    const t = pack(T),
      v = pack(a);
    const r = gen_pow_i(
      v,
      1n << BigInt(d - 1),
      (b) => F2x_rem(F2x_sqr(b), t),
      (b, c) => F2x_rem(F2x_mul(b, c), t)
    );
    return result(
      Array.from({ length: r === 0n ? 0 : r.toString(2).length }, (_, i) => (r >> BigInt(i)) & 1n)
    );
  }
  const word = p < 1n << 64n,
    ctx = polynomialQuotient(T, p, word);
  const mul = (b: bigint[], c: bigint[]) => ctx.reduce(ctx.multiply(b, c));
  const sqr = (b: bigint[]) => ctx.reduce(word ? Flx_sqr(b, p) : FpX_sqr(b, p));
  const pow = (b: bigint[], n: bigint): bigint[] =>
    n === 0n ? [1n] : gen_pow_i(n < 0n ? polynomialQuotientInverse(b, T, p, word) : b, n, sqr, mul);
  const sample = () => (word ? random_Flx(d, p) : random_FpX(d, p));
  const one = (b: bigint[]) => b.length === 1 && b[0] === 1n;
  if (d === 1) {
    // F=1 leaves the n=2 factor for gen_Shanks_sqrtn. Preserve its random
    // Sylow-generator selection and signed Bezout exponent, even for squares.
    let r = p - 1n,
      e = 0;
    while (!(r & 1n)) {
      r >>= 1n;
      e++;
    }
    let y: bigint[];
    for (;;) {
      let c: bigint[];
      do c = sample();
      while (!c.length);
      y = pow(c, r);
      let m = y;
      if (one(m)) continue;
      let i = 1;
      for (; i < e; i++) {
        m = sqr(m);
        if (one(m)) break;
      }
      if (i === e) break;
    }
    let v = pow(a, (1n - r) / 2n),
      w = mul(sqr(v), pow(a, -1n));
    while (!one(w)) {
      let k = 0,
        b = w;
      do {
        b = sqr(b);
        k++;
      } while (!one(b));
      if (k === e) return null;
      const correction = pow(pow(y!, -1n), 1n << BigInt(e - k - 1));
      e = k;
      v = mul(correction, v);
      y = sqr(correction);
      w = mul(y, w);
    }
    return result(v);
  }
  // For n=2 and odd p, ord_2(p)=1 divides every extension degree >1.
  // Both native sqrtn dispatchers therefore choose sqrtl_spec, with s=1.
  const xi = pow([0n, 1n], p);
  const sumautsum = (b: bigint[], i: number): bigint[] => {
    const zeta = ctx.evaluate(b, ctx.powers(xi, Math.floor(Math.sqrt(b.length - 1))));
    const xiPowers = ctx.powers(
      xi,
      brent_kung_optpow(d - 1, 2 * (BigInt(i).toString(2).replaceAll('0', '').length - 1), 1)
    );
    type Triple = [bigint[], bigint[], bigint[]];
    const square = ([xx, zz, dd]: Triple): Triple => {
      const powers = ctx.powers(xx, brent_kung_optpow(d - 1, 3, 1));
      return [
        ctx.evaluate(xx, powers),
        mul(zz, ctx.evaluate(zz, powers)),
        FpX_add(dd, mul(zz, ctx.evaluate(dd, powers)), p),
      ];
    };
    const multiplySquare = (v: Triple): Triple => {
      const s = square(v),
        zz = mul(zeta, ctx.evaluate(s[1], xiPowers));
      return [ctx.evaluate(s[0], xiPowers), zz, FpX_add(s[2], zz, p)];
    };
    const v = gen_pow_fold<Triple>([xi, zeta, zeta], BigInt(i), square, multiplySquare);
    return mul(b, FpX_add([1n], v[2], p));
  };
  let c: bigint[], b: bigint[], z: bigint[];
  do {
    do c = sample();
    while (!c.length);
    z = mul(a, pow(c, 2n));
    const alpha = pow(z, (p - 1n) / 2n);
    b = FpX_add([1n], d === 2 ? alpha : sumautsum(alpha, d - 2), p);
  } while (!b.length);
  const scalar = mul(z, pow(b, 2n));
  if (scalar.length > 1) return null;
  // The inner Fp/Fl_sqrtn exponent is exactly two, selecting the scalar sqrt.
  const beta = word ? Fl_sqrt(scalar[0] ?? 0n, p) : Fp_sqrt(scalar[0] ?? 0n, p);
  if (beta === null) return null;
  return result(FpX_Fp_mul(polynomialQuotientInverse(mul(b, c), T, p, word), beta, p));
}
