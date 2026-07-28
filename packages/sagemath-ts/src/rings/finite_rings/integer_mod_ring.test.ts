/**
 * Tests for sage/rings/finite_rings/integer_mod_ring
 *
 * The module had no test file; the generator methods below are ports of
 * `IntegerModRing_generic.multiplicative_group_is_cyclic` / `unit_gens` /
 * `multiplicative_generator` (`integer_mod_ring.py:810-895`, `:1442`).
 */

import { describe, expect, test } from 'bun:test';
import { Zmod, Integers, IntegerModRing } from './integer_mod_ring.js';

describe('IntegerModRing', () => {
  test('exposes the generator methods', () => {
    expect(Zmod(7n).multiplicative_generator().value).toBe(3n);
    expect(Zmod(9n).multiplicative_generator().value).toBe(2n);
    expect(Zmod(75n).multiplicative_group_is_cyclic()).toBe(false);
    expect(Zmod(75n).unit_gens().map((g) => g.value)).toEqual([26n, 52n]);
    expect(() => Zmod(8n).multiplicative_generator()).toThrow(
      'multiplicative group of this ring is not cyclic'
    );
  });

  test('basic ring operations', () => {
    const R = Zmod(12n);
    expect(R.__call__(7n).add(R.__call__(8n)).value).toBe(3n);
    expect(R.is_field()).toBe(false);
    expect(Zmod(7n).is_field()).toBe(true);
    expect(R.cardinality()).toBe(12n);
    expect(R.units().map((u) => u.value)).toEqual([1n, 5n, 7n, 11n]);
    expect(R.toString()).toBe('Ring of integers modulo 12');
  });
});

describe('IntegerModRing construction', () => {
  test('aliases and validation', () => {
    expect(Integers(5n).modulus).toBe(5n);
    expect(new IntegerModRing(5).modulus).toBe(5n);
    expect(() => Zmod(0n)).toThrow('modulus must be positive');
    expect(() => Zmod(-3n)).toThrow('modulus must be positive');
  });

  test('iteration and coercion', () => {
    const R = Zmod(5n);
    expect([...R].map((x) => x.value)).toEqual([0n, 1n, 2n, 3n, 4n]);
    expect(R.__call__(-1n).value).toBe(4n);
    expect(R.__call__(true).value).toBe(1n);
    expect(R.zero().value).toBe(0n);
    expect(R.one().value).toBe(1n);
    expect(R.gen().value).toBe(1n);
  });
});
