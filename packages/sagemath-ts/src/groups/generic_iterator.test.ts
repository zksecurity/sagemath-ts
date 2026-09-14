import { describe, expect, test } from 'bun:test';
import { ValueError } from '../errors.js';
import { QuadraticField } from '../rings/number_field/number_field.js';
import { multiples, parseGroupOps } from './generic.js';

class Box {
  constructor(public value: number) {}
  eq(other: Box): boolean {
    return this.value === other.value;
  }
}
const add = (x: Box, y: Box) => new Box((x.value + y.value) % 11);

describe('Native multiple-iterator state and copies', () => {
  test('validates during construction', () => {
    expect(() => multiples(new Box(2), -1n, new Box(5), false, 'other', add)).toThrow(
      'n cannot be negative in multiples'
    );
    expect(() => multiples(new Box(2), 1n, undefined, false, 'other', add)).toThrow(
      'P0 must be supplied when operation is neither addition nor multiplication'
    );
  });

  test('copies arguments and advances before returning a value', () => {
    const p = new Box(2),
      p0 = new Box(5);
    let calls = 0;
    const iterator = multiples(p, 2n, p0, false, 'other', (x, y) => {
      calls++;
      return add(x, y);
    });
    p.value = 9;
    p0.value = 8;
    const first = iterator.next();
    expect(first.done).toBe(false);
    expect((first.value as Box).value).toBe(5);
    expect(calls).toBe(1);
    expect((iterator.next().value as Box).value).toBe(7);
    expect(calls).toBe(2);
  });

  test('a callback error advances the index but leaves the iterator usable', () => {
    let calls = 0;
    const iterator = multiples(new Box(2), 3n, new Box(5), true, 'other', (x, y) => {
      if (++calls === 2) throw new ValueError('step blocked');
      return add(x, y);
    });
    const first = iterator.next().value as [bigint, Box];
    expect([first[0], first[1].value]).toEqual([0n, 5]);
    expect(() => iterator.next()).toThrow('step blocked');
    const third = iterator.next().value as [bigint, Box];
    expect([third[0], third[1].value]).toEqual([2n, 7]);
    expect(iterator.next().done).toBe(true);
    expect(calls).toBe(3);
  });

  test('initial-value copy failures happen after the step is copied', () => {
    for (const fail of [false, true]) {
      let copies = 0;
      class CopyBox extends Box {
        __copy__(): Box | null {
          copies++;
          if (this.value === 5) {
            if (fail) throw new ValueError('copy blocked');
            return null;
          }
          return new Box(this.value);
        }
      }
      expect(() => multiples(new CopyBox(2), 0n, new CopyBox(5), false, 'other', add)).toThrow(
        fail ? 'copy blocked' : 'P and Q must not be None'
      );
      expect(copies).toBe(2);
    }
  });

  test('uses a parent method for the standard identity', () => {
    let zeroCalls = 0;
    class AddBox extends Box {
      parent() {
        return {
          zero() {
            zeroCalls++;
            return new AddBox(0);
          },
        };
      }
      add(other: Box): AddBox {
        return new AddBox((this.value + other.value) % 11);
      }
    }
    const iterator = multiples(new AddBox(2), 1n);
    expect(zeroCalls).toBe(1);
    expect((iterator.next().value as Box).value).toBe(0);
  });

  test('parses number-field elements whose parent is a method', () => {
    const K = QuadraticField.create(2n, 'a'),
      a = K.gen();
    const ops = parseGroupOps('*', undefined, undefined, undefined, a);
    expect(ops.identity.eq(K.one())).toBe(true);
    expect(ops.inverse(a).mul(a).eq(K.one())).toBe(true);
    expect(ops.power(a, -3n).list().map(String)).toEqual(['0', '1/4']);
  });

  test('retains the JavaScript return and throw closure protocol', () => {
    const iterator = multiples(new Box(2), 3n, new Box(5), false, 'other', add);
    expect(iterator.return().done).toBe(true);
    expect(iterator.next().done).toBe(true);
    const other = multiples(new Box(2), 3n, new Box(5), false, 'other', add);
    const error = new Error('cancelled');
    expect(() => other.throw(error)).toThrow(error);
    expect(other.next().done).toBe(true);
  });

  test('retains JavaScript for-of closure and iterator identity', () => {
    const iterator = multiples(new Box(2), 3n, new Box(5), false, 'other', add);
    expect(iterator[Symbol.iterator]()).toBe(iterator);
    for (const value of iterator) {
      expect(value.value).toBe(5);
      break;
    }
    expect(iterator.next().done).toBe(true);
  });
});
