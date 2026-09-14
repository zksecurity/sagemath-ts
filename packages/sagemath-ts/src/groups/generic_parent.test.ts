import { expect, test } from 'bun:test';
import { ValueError } from '../errors.js';
import { discrete_log, order_from_bounds, order_from_multiple } from './generic.js';

function trackedGroup(failParent = false, failEquality = false) {
  const trace: string[] = [];
  const parent = {
    zero: () => {
      trace.push('zero');
      if (failParent) throw new ValueError('identity sentinel');
      return new Element(0n);
    },
  };
  class Element {
    parent = parent;
    constructor(readonly value: bigint) {}
    add(other: Element) {
      return new Element((this.value + other.value) % 5n);
    }
    neg() {
      return new Element((5n - this.value) % 5n);
    }
    mul(n: bigint) {
      trace.push(`scale ${n}`);
      return new Element((((this.value * n) % 5n) + 5n) % 5n);
    }
    eq(other: Element) {
      trace.push(`eq ${this.value} ${other.value}`);
      if (failEquality) throw new TypeError('equality sentinel');
      return this.value === other.value;
    }
    is_zero(): boolean {
      throw new Error('order reduction must use equality with the parent identity');
    }
    toString() {
      return String(this.value);
    }
  }
  return { trace, element: (n: bigint) => new Element(n) };
}

test('order_from_multiple obtains the identity before its equality shortcut', () => {
  const { trace, element } = trackedGroup();
  expect(order_from_multiple(element(0n), 5n, undefined, '+')).toBe(1n);
  expect(trace).toEqual(['zero', 'eq 0 0']);
});

test('order_from_multiple uses equality throughout native order reduction', () => {
  const { trace, element } = trackedGroup();
  expect(order_from_multiple(element(1n), 10n, undefined, '+')).toBe(5n);
  expect(trace[0]).toBe('zero');
  expect(trace.filter((x) => x === 'zero')).toHaveLength(1);
});

test('order_from_bounds shares its parsed parent identity with both nested algorithms', () => {
  const { trace, element } = trackedGroup();
  expect(order_from_bounds(element(1n), [1n, 10n], 2n, '+')).toBe(5n);
  expect(trace.filter((x) => x === 'zero')).toHaveLength(1);
});

test('automatic order bounds retain the original native power function', () => {
  const { trace, element } = trackedGroup();
  expect(order_from_bounds(element(1n), undefined, 3n, '+')).toBe(5n);
  expect(trace.filter((x) => x === 'zero')).toHaveLength(1);
  expect(trace.some((x) => x.startsWith('scale '))).toBe(true);
});

test('discrete_log resolves its parent before powers and shares it with BSGS', () => {
  const { trace, element } = trackedGroup();
  expect(discrete_log(element(3n), element(1n), 5n, '+')).toBe(3n);
  expect(trace[0]).toBe('zero');
  expect(trace.filter((x) => x === 'zero')).toHaveLength(1);
});

test('parent errors propagate before order or logarithm arithmetic', () => {
  for (const run of [
    (x: ReturnType<ReturnType<typeof trackedGroup>['element']>) =>
      order_from_multiple(x, 5n, undefined, '+'),
    (x: ReturnType<ReturnType<typeof trackedGroup>['element']>) =>
      order_from_bounds(x, undefined, 1n, '+'),
    (x: ReturnType<ReturnType<typeof trackedGroup>['element']>) => discrete_log(x, x, 5n, '+'),
  ]) {
    const { trace, element } = trackedGroup(true);
    expect(() => run(element(0n))).toThrow('identity sentinel');
    expect(trace).toEqual(['zero']);
  }
});

test('identity equality errors are not bypassed by a zero predicate', () => {
  const { trace, element } = trackedGroup(false, true);
  expect(() => order_from_multiple(element(0n), 5n, undefined, '+')).toThrow('equality sentinel');
  expect(trace).toEqual(['zero', 'eq 0 0']);
});
