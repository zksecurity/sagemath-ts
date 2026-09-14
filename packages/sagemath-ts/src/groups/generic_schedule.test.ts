import { expect, test } from 'bun:test';
import * as g from './generic.js';
import { bsgs, discrete_log, order_from_bounds, order_from_multiple } from './generic.js';
import { RDF } from '../rings/real_double.js';

// These operation traces were captured from the bundled Sage generic.py.
test('BSGS inverts first and reuses the baby-step state', () => {
  expect(JSON.parse(group_schedule(1n, 0n, 36n, 1n, 35n, 2n, 72n, 1n))).toEqual([
    '35',
    null,
    [
      ['neg', 35],
      ['add', 1, 1],
      ['add', 1, 2],
      ['add', 3, 1],
      ['add', 4, 1],
      ['add', 5, 1],
      ['add', 6, 1],
      ['add', 7, 1],
      ['add', 8, 1],
      ['add', 9, 1],
      ['add', 10, 1],
      ['add', 11, 1],
      ['neg', 12],
      ['add', 3, 24],
      ['add', 27, 0],
      ['add', 27, 27],
      ['add', 27, 18],
    ],
  ]);
});

test('has_order computes both recursive powers before testing', () => {
  expect(JSON.parse(group_schedule(2n, 0n, 36n, 1n, 30n, 1n, 72n, 1n))).toEqual([
    false,
    null,
    [
      ['add', 1, 1],
      ['add', 2, 1],
      ['add', 1, 1],
      ['add', 2, 2],
      ['add', 4, 4],
      ['add', 2, 8],
      ['add', 3, 3],
      ['add', 6, 6],
      ['add', 12, 3],
      ['add', 3, 3],
      ['add', 15, 15],
    ],
  ]);
});

test('order_from_bounds uses binary scaling and native order reduction', () => {
  expect(JSON.parse(group_schedule(3n, 0n, 36n, 1n, 1n, 1n, 72n, 3n))).toEqual([
    '36',
    null,
    [
      ['add', 1, 1],
      ['add', 2, 1],
      ['neg', 0],
      ['add', 0, 3],
      ['add', 3, 3],
      ['add', 3, 6],
      ['add', 3, 9],
      ['add', 3, 12],
      ['add', 3, 15],
      ['add', 3, 18],
      ['add', 3, 21],
      ['add', 3, 24],
      ['add', 3, 27],
      ['add', 3, 30],
      ['add', 3, 33],
      ['scale', 1, 9],
      ['scale', 9, 2],
      ['scale', 1, 4],
      ['scale', 4, 3],
    ],
  ]);
});

test('order_from_multiple retains native scalar powers', () => {
  expect(JSON.parse(group_schedule(4n, 0n, 36n, 1n, 36n, 1n, 72n, 1n))).toEqual([
    '36',
    null,
    [
      ['scale', 1, 36],
      ['scale', 1, 9],
      ['scale', 9, 2],
      ['scale', 1, 4],
      ['scale', 4, 3],
    ],
  ]);
});

test('discrete_log uses zero and negative native powers', () => {
  expect(JSON.parse(group_schedule(5n, 0n, 36n, 1n, 3n, 0n, 72n, 1n))).toEqual([
    '3',
    null,
    [
      ['scale', 1, 0],
      ['scale', 1, 18],
      ['scale', 18, 0],
      ['scale', 1, 0],
      ['add', 3, 0],
      ['scale', 3, 18],
      ['neg', 18],
      ['add', 18, 0],
      ['add', 18, 18],
      ['scale', 1, -1],
      ['add', 3, 35],
      ['scale', 2, 9],
      ['neg', 18],
      ['add', 18, 0],
      ['add', 18, 18],
      ['scale', 1, 12],
      ['scale', 12, 0],
      ['scale', 1, 0],
      ['add', 3, 0],
      ['scale', 3, 12],
      ['neg', 0],
      ['add', 0, 0],
      ['scale', 1, 0],
      ['add', 3, 0],
      ['scale', 3, 4],
      ['neg', 12],
      ['add', 24, 0],
      ['add', 12, 24],
      ['scale', 1, 3],
    ],
  ]);
});

test('BSGS resolves equal printed keys with element equality', () => {
  expect(JSON.parse(group_schedule(1n, 4n, 36n, 1n, 35n, 0n, 72n, 1n))).toEqual([
    '35',
    null,
    [
      ['neg', 35],
      ['add', 1, 0],
      ['add', 1, 1],
      ['add', 2, 1],
      ['add', 3, 1],
      ['add', 4, 1],
      ['add', 5, 1],
      ['add', 6, 1],
      ['add', 7, 1],
      ['add', 8, 1],
      ['add', 9, 1],
      ['neg', 10],
      ['add', 1, 26],
      ['add', 27, 0],
      ['add', 27, 27],
      ['add', 27, 18],
    ],
  ]);
});

test('BSGS does not equate distinct NaN intermediates', () => {
  expect(() => bsgs(RDF.__call__(Infinity), RDF.__call__(Infinity), [0n, 64n], '+')).toThrow(
    'log of +infinity to the base +infinity does not exist in (0, 64)'
  );
});

test('BSGS zero powers do not create a false logarithm through NaN', () => {
  expect(() => bsgs(RDF.zero(), RDF.zero(), [0n, 64n], '*')).toThrow(
    'log of 0.0 to the base 0.0 does not exist in (0, 64)'
  );
});

test('order bounds reject a NaN base through native BSGS failure', () => {
  expect(() => order_from_bounds(RDF.__call__(NaN), [1n, 64n], 2n, '+')).toThrow(
    'log of 0.0 to the base NaN does not exist in (1, 32)'
  );
});

test('generic f-string failures use RDF float formatting', () => {
  expect(() => discrete_log(RDF.__call__(Infinity), RDF.__call__(Infinity), 12n, '+')).toThrow(
    'no discrete log of inf found to base inf'
  );
  expect(() => order_from_multiple(RDF.__call__(NaN), 12n, undefined, '*')).toThrow(
    'The order of P(=nan) does not divide 12'
  );
});

// Finite modular group with explicit operation tracing for the Sage transcripts.
function group_schedule(
  fn: bigint,
  mode: bigint,
  modulus: bigint,
  value: bigint,
  n: bigint,
  lb: bigint,
  ub: bigint,
  d: bigint,
  impl = g
): string {
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
  const parent = { zero: () => new Box(0n), one: () => new Box(1n) };
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
      return this.v === b.v;
    }
    is_zero() {
      return this.v === 0n;
    }
    is_one() {
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
  try {
    if (fn === 0n) r = String(impl.multiple(a, n, operation, identity, inverse, op));
    else if (fn === 1n)
      r = String(impl.bsgs(a, target, [lb, ub], operation, identity, inverse, op));
    else if (fn === 2n) r = impl.has_order(a, n, operation);
    else if (fn === 3n)
      r = String(impl.order_from_bounds(a, [lb, ub], d, operation, identity, inverse, op));
    else if (fn === 4n)
      r = String(impl.order_from_multiple(a, n, undefined, operation, identity, inverse, op));
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
    if (!(e instanceof Error)) throw e;
    r = null;
    error = e.name + ': ' + e.message;
  }
  return JSON.stringify([r, error, trace]);
}
