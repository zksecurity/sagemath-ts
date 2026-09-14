import { ValueError } from '../../../packages/sagemath-ts/src/errors.js';
import * as g from '../../../packages/sagemath-ts/src/groups/generic.js';
export function group_parent(
  fn: bigint,
  mode: bigint,
  modulus: bigint,
  value: bigint,
  n: bigint,
  lb: bigint,
  ub: bigint,
  d: bigint,
  variant: bigint = 0n,
  impl = g
): string {
  let active = false;
  const trace: unknown[][] = [],
    additive = mode % 2n === 0n;
  const norm = (v: bigint) => ((v % modulus) + modulus) % modulus;
  const pow = (v: bigint, k: bigint) => {
    if (k < 0n) {
      v = pow(v, modulus - 2n);
      k = -k;
    }
    let r = 1n;
    while (k > 0n) {
      if (k & 1n) r = norm(r * v);
      v = norm(v * v);
      k >>= 1n;
    }
    return r;
  };
  const parent = {
    zero: () => {
      trace.push(['zero']);
      if (active && variant === 1n) throw new ValueError('identity sentinel');
      return new Box(0n);
    },
    one: () => {
      trace.push(['one']);
      if (active && variant === 1n) throw new ValueError('identity sentinel');
      return new Box(1n);
    },
  };
  class Box {
    v: bigint;
    parent = parent;
    constructor(v: bigint) {
      this.v = norm(v);
    }
    add(b: Box) {
      trace.push(['add', Number(this.v), Number(b.v)]);
      return new Box(this.v + b.v);
    }
    mul(b: Box | bigint) {
      trace.push([
        typeof b === 'bigint' ? 'scale' : 'mul',
        Number(this.v),
        Number(typeof b === 'bigint' ? b : b.v),
      ]);
      return new Box(this.v * (typeof b === 'bigint' ? b : b.v));
    }
    pow(k: bigint) {
      trace.push(['pow', Number(this.v), Number(k)]);
      return new Box(pow(this.v, k));
    }
    neg() {
      trace.push(['neg', Number(this.v)]);
      return new Box(-this.v);
    }
    inv() {
      trace.push(['inv', Number(this.v)]);
      return new Box(pow(this.v, modulus - 2n));
    }
    eq(b: Box) {
      trace.push(['eq', Number(this.v), Number(b.v)]);
      if (active && variant === 2n) throw new TypeError('equality sentinel');
      return this.v === b.v;
    }
    is_zero() {
      trace.push(['bool', Number(this.v)]);
      return this.v === 0n;
    }
    is_one() {
      trace.push(['is_one', Number(this.v)]);
      return this.v === 1n;
    }
    toString() {
      return mode >= 4n ? 'collision' : String(this.v);
    }
  }
  const a = new Box(value),
    target = new Box(additive ? value * n : pow(value, n));
  const operation = mode % 4n >= 2n ? 'other' : additive ? '+' : '*';
  const identity = mode % 4n >= 2n ? (additive ? parent.zero() : parent.one()) : undefined;
  const inverse = mode % 4n >= 2n ? (x: Box) => (additive ? x.neg() : x.inv()) : undefined;
  const op = mode % 4n >= 2n ? (x: Box, y: Box) => (additive ? x.add(y) : x.mul(y)) : undefined;
  let r: unknown,
    error = null;
  trace.length = 0;
  active = true;
  try {
    if (fn === 0n) r = String(impl.multiple(a, n, operation, identity, inverse, op));
    else if (fn === 1n)
      r = String(impl.bsgs(a, target, [lb, ub], operation, identity, inverse, op));
    else if (fn === 2n) r = impl.has_order(a, n, operation);
    else if (fn === 3n)
      r = String(impl.order_from_bounds(a, [lb, ub], d, operation, identity, inverse, op));
    else if (fn === 4n)
      r = String(impl.order_from_multiple(a, n, undefined, operation, identity, inverse, op));
    else if (fn === 6n)
      r = String(impl.order_from_bounds(a, undefined, d, operation, identity, inverse, op));
    else
      r = String(
        impl.discrete_log(
          target,
          a,
          additive ? modulus : modulus - 1n,
          operation,
          identity,
          inverse,
          op
        )
      );
  } catch (e) {
    r = null;
    error = e.name + ': ' + e.message;
  }
  return JSON.stringify([r, error, trace]);
}
