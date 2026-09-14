/** Mixed integer/rational coercion and generic-group dispatch. */
import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import * as g from '../../../packages/sagemath-ts/src/groups/generic.js';
export function integer_rational(
  fn: bigint,
  n: bigint,
  m: bigint,
  den: bigint,
  kind: bigint
): string {
  const a = new Integer(n),
    b = kind === 0n ? new Rational(m, den) : kind === 1n ? m : new Integer(m);
  const methods = [
    () => a.add(b),
    () => a.sub(b),
    () => a.mul(b),
    () => a.eq(b),
    () => a.lt(b),
    () => a.le(b),
    () => a.gt(b),
    () => a.ge(b),
  ];
  const r = methods[Number(fn)]!();
  return JSON.stringify([
    String(r),
    typeof r === 'boolean' ? 'bool' : r instanceof Rational ? 'Rational' : 'Integer',
  ]);
}
export function integer_group(
  kind: bigint,
  fn: bigint,
  mode: bigint,
  value: bigint,
  n: bigint
): string {
  const a = kind === 0n ? new Integer(value) : new Rational(value),
    operation = mode === 0n ? '+' : '*';
  const target = operation === '+' ? a.mul(n) : a.pow(n);
  let r: unknown;
  if (fn === 0n) r = g.multiple(a, n, operation);
  else if (fn === 1n) r = g.bsgs(a, target, [0n, 64n], operation);
  else if (fn === 2n) r = g.has_order(a, n, operation);
  else if (fn === 3n) r = g.order_from_multiple(a, n, undefined, operation);
  else if (fn === 4n) r = g.discrete_log(target, a, 4n, operation);
  else if (fn === 5n) r = g.order_from_bounds(a, [1n, 16n], undefined, operation);
  else {
    const ops = g.parseGroupOps(operation, undefined, undefined, undefined, a);
    return JSON.stringify([
      String(ops.identity),
      String(ops.inverse(a)),
      String(ops.power(a, n)),
      ops.isIdentity(a),
    ]);
  }
  return JSON.stringify(typeof r === 'boolean' ? r : String(r));
}
