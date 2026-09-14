/** Algorithms supplied by Sage's Fields.ParentMethods category.
 * @see Deviation: Polynomial squarefree backend dispatch
 */
import {NotImplementedError} from '../errors.js';
import {Polynomial, type RingElement} from '../rings/polynomial/polynomial_element.js';
/** @internal categories/fields.py:665; the receiver is f's coefficient field. */
export function _squarefree_decomposition_univariate_polynomial<C extends RingElement>(f: Polynomial<C>): Array<[Polynomial<C>, number]> {
  const base = f.parent.base_ring;
  if (f.degree() === 0) return [[f, 1]];
  const characteristic = (base as unknown as {characteristic: bigint | (() => bigint)}).characteristic;
  const p = typeof characteristic === 'function' ? characteristic.call(base) : characteristic;
  if (p !== 0n) throw new NotImplementedError('square-free decomposition not implemented for this polynomial');
  const chain = [f];
  let current = f;
  while (current.degree() > 0) { current = current.gcd(current.derivative()); chain.push(current); }
  const g = chain.slice(0, -1).map((v, i) => v.quo_rem(chain[i + 1]!)[0]);
  const a = g.slice(0, -1).map((v, i) => v.quo_rem(g[i + 1]!)[0]);
  a.push(g.at(-1)!);
  let unit = chain.at(-1)!;
  const factors: Array<[Polynomial<C>, number]> = [];
  for (let i = 0; i < a.length; i++) {
    if (a[i]!.degree() > 0) factors.push([a[i]!, i + 1]);
    else unit = unit.mul(a[i]!.pow(i + 1) as Polynomial<C>);
  }
  if (!unit.eq(1n)) factors.unshift([unit, 1]);
  return factors;
}
