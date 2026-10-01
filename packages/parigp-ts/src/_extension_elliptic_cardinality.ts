/** Constant-j twist counts, PARI FpE.c:2064–2132 / FlxqE.c:1378–1440.
 * Copyright (C) The PARI group; adapted under GPL-2.0-or-later.
 */
import { extensionField } from './_extension_field.js';
import { gen_pow_i } from './bb_group.js';
import { polynomialQuotientPower } from './_polynomial_quotient_power.js';
import { Flxq_issquare } from './Flx.js';
import { FpXQ_issquare } from './FpX.js';
import { Fp_ffellcard } from './FpE.js';
import { Fp_div, Fp_pow } from './ff.js';
import { isqrt } from './ifactor.js';

/** Valid reduced field inputs, p>3, q=p^n and constant j. */
export function constantJCard(
  mode: 0 | 1, a4: bigint[], a6: bigint[], j: bigint,
  T: bigint[], q: bigint, p: bigint, n: number
): bigint {
  const F = extensionField(mode, T, p), q1 = q + 1n;
  const pow = (a: bigint[], e: bigint): bigint[] => mode === 0
    ? polynomialQuotientPower(a, e, T, p)
    : e === 0n ? [1n] : gen_pow_i(a, e, x => F.sqr(x) as bigint[], (x,y) => F.mul(x,y) as bigint[]);
  const square = (a: bigint[]) => mode === 1 ? Flxq_issquare(a,T,p) : FpXQ_issquare(a,T,p);
  const red = (a: bigint) => ((a % p) + p) % p;
  if (j === 0n) {
    if (q % 6n !== 1n) return q1;
    const N = Fp_ffellcard(0n,1n,q,n,p), t = q1 - N;
    const W = pow(a6, (q / 2n) / 3n);
    if (W.length > 1) {
      const cube = pow(W, 3n);
      return cube.length === 1 && cube[0] === 1n ? q1 + t/2n : q1 - t/2n;
    }
    const w = W[0]!;
    if (w === 1n) return N;
    if (w === p-1n) return q1+t;
    const u = t/2n, v = isqrt((q-u*u)/3n), a = u+v, b = 2n*v;
    if (Fp_pow(w,3n,p) === 1n)
      return red(a+w*b) === 0n ? q1-(2n*b-a) : q1+a+b;
    return red(a-w*b) === 0n ? q1-(a-2n*b) : q1-a-b;
  }
  if (j === 1728n % p) {
    if (q % 4n === 3n) return q1;
    const W = pow(a4,q/4n);
    if (W.length > 1) return q1;
    const w = W[0]!, N = Fp_ffellcard(1n,0n,q,n,p);
    if (w === 1n) return N;
    const t = q1-N;
    if (w === p-1n) return q1+t;
    const u = t/2n, v = isqrt(q-u*u);
    return red(u+w*v) === 0n ? q1-2n*v : q1+2n*v;
  }
  const g = Fp_div(j,red(1728n-j),p);
  const numerator = F.mul(a6,3n), denominator = F.mul(a4,2n);
  // The native word routine uses a product, while FpXQ uses a quotient.
  // Keep both routes even though their quadratic characters agree.
  let l: bigint[], N: bigint;
  if (mode === 1) {
    N = Fp_ffellcard(3n*g%p,2n*g%p,q,n,p);
    l = F.mul(numerator,denominator) as bigint[];
  } else {
    l = F.mul(numerator,F.inv(denominator)) as bigint[];
    N = Fp_ffellcard(3n*g%p,2n*g%p,q,n,p);
  }
  return square(l) ? N : 2n*q1-N;
}
