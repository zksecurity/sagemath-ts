import * as g from '../../../packages/sagemath-ts/src/groups/generic.js';
import {
  set_random_seed,
  current_randstate,
} from '../../../packages/sagemath-ts/src/misc/randstate.js';
import { ValueError } from '../../../packages/sagemath-ts/src/errors.js';
export function pollard(
  algorithm: bigint,
  seed: bigint,
  mode: bigint,
  modulus: bigint,
  value: bigint,
  n: bigint,
  order: bigint,
  lb: bigint,
  ub: bigint,
  mutable: bigint,
  collision: bigint,
  failParent: bigint,
  impl = g
): string {
  set_random_seed(seed);
  const trace: unknown[][] = [];
  let active = false;
  const additive = mode % 2n === 0n;
  const norm = (v: bigint) => ((v % modulus) + modulus) % modulus;
  const pow = (v: bigint, k: bigint): bigint => {
    if (k < 0n) {
      v = pow(v, modulus - 2n);
      k = -k;
    }
    let r = 1n;
    while (k) {
      if (k & 1n) r = norm(r * v);
      v = norm(v * v);
      k >>= 1n;
    }
    return r;
  };
  const check = (x: Box) => {
    if (mutable === 2n && !x.frozen) throw new TypeError('mutable element');
  };
  const parent = {
    zero: () => {
      if (active) {
        trace.push(['zero']);
        if (failParent) throw new ValueError('identity sentinel');
      }
      return new Box(0n);
    },
    one: () => {
      if (active) {
        trace.push(['one']);
        if (failParent) throw new ValueError('identity sentinel');
      }
      return new Box(1n);
    },
  };
  class Box {
    v: bigint;
    frozen = false;
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
      return this.v === b.v;
    }
    toString() {
      check(this);
      return collision === 1n
        ? 'collision'
        : String(collision === 2n ? this.v % 4n : collision === 3n ? this.v % 17n : this.v);
    }
    set_immutable() {
      trace.push(['freeze', Number(this.v)]);
      this.frozen = true;
      if (mutable === 3n) throw new ValueError('freeze sentinel');
    }
  }
  if (!mutable) delete (Box.prototype as { set_immutable?: () => void }).set_immutable;
  const base = new Box(value),
    target = new Box(additive ? value * n : pow(value, n));
  const operation = mode >= 2n ? 'other' : additive ? '+' : '*',
    identity = mode >= 2n ? (additive ? parent.zero() : parent.one()) : undefined;
  const inverse = mode >= 2n ? (x: Box) => (additive ? x.neg() : x.inv()) : undefined,
    op = mode >= 2n ? (x: Box, y: Box) => (additive ? x.add(y) : x.mul(y)) : undefined;
  trace.length = 0;
  active = true;
  let r = null,
    error = null;
  try {
    r = String(
      algorithm === 0n
        ? impl.discrete_log_lambda(
            target,
            base,
            [lb, ub],
            operation,
            identity,
            inverse,
            op,
            (x) => {
              check(x);
              return x.v * x.v + 1n;
            }
          )
        : impl.discrete_log_rho(target, base, order, operation, identity, inverse, op)
    );
  } catch (e) {
    error = e.name + ': ' + e.message;
  }
  return JSON.stringify([
    r,
    error,
    trace,
    base.frozen,
    target.frozen,
    String(current_randstate().python_random().getrandbits(64)),
  ]);
}

import {
  matrix_gf2_from_entries,
  zero_matrix_gf2,
  identity_matrix_gf2,
} from '../../../packages/sagemath-ts/src/matrix/matrix_mod2.js';
export function pollard_matrix(
  seed: bigint,
  mode: bigint,
  n: bigint,
  lb: bigint,
  ub: bigint,
  impl = g
): string {
  set_random_seed(seed);
  const base = matrix_gf2_from_entries([
      [0, 1],
      [1, 1],
    ]),
    identity = mode === 0n ? zero_matrix_gf2(2, 2) : identity_matrix_gf2(2);
  const op = (x: typeof base, y: typeof base) => (mode === 0n ? x.add(y) : x.mul(y)),
    inverse = (x: typeof base) => (mode === 0n ? x.neg() : x.inverse());
  let target = identity;
  for (let i = 0n; i < n; i++) target = op(target, base);
  let r = null,
    error = null;
  try {
    r = String(
      impl.discrete_log_lambda(target, base, [lb, ub], 'other', identity, inverse, op, (x) => {
        let v = 0n;
        for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) v = 2n * v + BigInt(x.get(i, j));
        return v;
      })
    );
  } catch (e) {
    error = e.name + ': ' + e.message;
  }
  return JSON.stringify([
    r,
    error,
    base.is_immutable(),
    target.is_immutable(),
    String(current_randstate().python_random().getrandbits(64)),
  ]);
}
