import { expect, test } from 'bun:test';
import { generic_power } from './power.js';
import { Integer } from '../rings/integer_ring.js';

test('generic_power follows Sage left multiplication and omits identity products', () => {
  const trace: unknown[][] = [];
  const ring = {
    one: () => {
      trace.push(['one']);
      return new Box(1n);
    },
  };
  class Box {
    constructor(readonly value: bigint) {}
    parent() {
      return ring;
    }
    mul(b: Box) {
      trace.push([this.value, b.value]);
      return new Box(this.value * b.value);
    }
    inv(): Box {
      throw new Error('unused');
    }
  }
  const a = new Box(2n);
  expect(generic_power(a, 1n)).toBe(a);
  expect(trace).toEqual([]);
  expect(generic_power(a, new Integer(13n)).value).toBe(8192n);
  expect(trace).toEqual([
    [2n, 2n],
    [4n, 4n],
    [16n, 2n],
    [16n, 16n],
    [256n, 32n],
  ]);
  trace.length = 0;
  expect(generic_power(a, 0n).value).toBe(1n);
  expect(trace).toEqual([['one']]);
});
