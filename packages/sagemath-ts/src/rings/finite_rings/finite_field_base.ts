/** Algorithms from Sage's abstract FiniteField base.
 * @see Deviation: Polynomial squarefree backend dispatch
 */
import {Polynomial, type RingElement} from '../polynomial/polynomial_element.js';
/** @internal finite_field_base.pyx:565; explicit characteristic/order replace the field context. */
export function _squarefree_decomposition_univariate_polynomial<C extends RingElement>(f: Polynomial<C>, p: bigint, q: bigint): Array<[Polynomial<C>, number]> {
  if (f.degree() === 0) return [[f, 1]];
  const factors: Array<[Polynomial<C>, number]> = [], R = f.parent;
  const unit = f.leading_coefficient();
  const pthRoot = (a: Polynomial<C>) => new Polynomial(Array.from({length: Math.floor(a.degree() / Number(p)) + 1}, (_, i) => {
    const c = a.getCoeff(i * Number(p)) as C & {pow(n: bigint): C};
    return q === p ? c : c.pow(q / p);
  }), R);
  let t0 = f._monic(), e = 1, degree = t0.degree(), derivative = t0.derivative();
  while (derivative.isZero()) { t0 = pthRoot(t0); degree = Math.floor(degree / Number(p)); derivative = t0.derivative(); e *= Number(p); }
  let t = t0.gcd(derivative), v = t0.quo_rem(t)[0], k = 0;
  while (degree > 0) {
    k++;
    if (BigInt(k) % p === 0n) { t = t.quo_rem(v)[0]; k++; }
    const w = v.gcd(t);
    if (w.degree() < v.degree()) {
      factors.push([v.quo_rem(w)[0], e * k]); v = w; t = t.quo_rem(v)[0];
      if (v.degree() === 0) {
        if (t.degree() === 0) break;
        t0 = pthRoot(t); degree = Math.floor(degree / Number(p)); derivative = t0.derivative(); e *= Number(p);
        while (derivative.isZero()) { t0 = pthRoot(t0); degree = Math.floor(degree / Number(p)); derivative = t0.derivative(); e *= Number(p); }
        t = t0.gcd(derivative); v = t0.quo_rem(t)[0]; k = 0;
      }
    } else t = t.quo_rem(v)[0];
  }
  if (!unit.eq(R.base_ring.one())) factors.unshift([new Polynomial([unit], R), 1]);
  return factors;
}
