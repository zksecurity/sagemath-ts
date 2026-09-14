import { Flx_gcd } from './Flx.js';
import { FpX_gcd } from './FpX.js';
import { brent_kung_optpow } from './RgX.js';
import type { PolynomialQuotient } from './_polynomial_quotient.js';
/** PARI FpX_factor.c Shoup distinct-degree splitting, using baby/giant steps. */
import { FpX_normalize, FpX_sub, type FpX as P } from './ffinit.js';

export function polynomialDistinctDegrees(context: PolynomialQuotient, XP: P): P[] {
  const { T, p, word, multiply, divide, reduce, powers, evaluate, power } = context;
  const n = T.length - 1;
  if (n === 0) return [];
  if (n === 1) return [T];
  const gcd = (a: P, b: P) => (word ? Flx_gcd(a, b, p) : FpX_gcd(a, b, p));
  const B = Math.floor(n / 2);
  const l = Math.floor(Math.sqrt(B));
  const m = Math.ceil(B / l);
  const b: P[] = [[0n, 1n], XP];
  const bo = brent_kung_optpow(n, l - 1, 1);
  const ro = l <= 1 ? 0 : Math.floor((bo - 1) / (l - 1)) + Math.floor((n - 1) / bo);
  if (word && p.toString(2).length - 1 <= ro) {
    for (let i = 2; i <= l; i++) b.push(power(b[i - 1]!));
  } else {
    const v = powers(b[1]!, bo);
    for (let i = 2; i <= l; i++) b.push(evaluate(b[i - 1]!, v));
  }
  const v = powers(b[l]!, brent_kung_optpow(n, m - 1, 1));
  const g: P[] = [v[1]!];
  for (let i = 1; i < m; i++) g.push(evaluate(g[i - 1]!, v));
  const H: P[] = g.map((a) => {
    let e = FpX_sub(a, b[0]!, p);
    for (let i = 1; i < l; i++) e = reduce(multiply(e, FpX_sub(a, b[i]!, p)));
    return e;
  });
  let Tr = T;
  const F: P[] = H.map((a) => {
    let u = gcd(Tr, a);
    if (u.length !== 1) {
      u = FpX_normalize(u, p);
      Tr = divide(Tr, u)[0];
    }
    return u;
  });
  const out: P[] = Array.from({ length: n }, () => [1n]);
  for (let j = 1; j <= m; j++) {
    let e = F[j - 1]!;
    for (let i = l - 1; i >= 0; i--) {
      let u = gcd(e, FpX_sub(g[j - 1]!, b[i]!, p));
      if (u.length !== 1) {
        if (!word) u = FpX_normalize(u, p);
        out[l * j - i - 1] = u;
        e = divide(e, u)[0];
      }
      if (e.length === 1) break;
    }
  }
  if (Tr.length !== 1) out[Tr.length - 2] = Tr;
  return out;
}
