/**
 * Tests for sage/rings/finite_rings/integer_mod (and integer_mod_ring)
 *
 * Neither module had a test file, which is how the truncated CRT in `log()`
 * (H15), the tautological local-solvability guard (L33) and the residue scan
 * in the primitive-root search (M25) survived.
 *
 * The expected values below are SageMath's, cited from the doctests in
 * `sage/rings/finite_rings/integer_mod.pyx` and
 * `sage/rings/finite_rings/integer_mod_ring.py`.
 */

import { describe, expect, test } from 'bun:test';
import { ValueError } from '../../errors.js';
import {
  Mod,
  multiplicative_generator,
  multiplicative_group_is_cyclic,
  unit_gens,
} from './integer_mod.js';

describe('IntegerMod.log', () => {
  test('running CRT covers every prime power (H15)', () => {
    // 105 = 3 * 5 * 7 has three prime factors; only the first two used to be
    // combined, so every exponent from 4 up came back wrong.
    const b = Mod(2n, 105n);
    expect(b.multiplicative_order()).toBe(12n);
    for (let e = 0n; e < 12n; e++) {
      expect(b.pow(e).log(b)).toBe(e);
    }
    // The exponent is reduced modulo the order of the base.
    expect(b.pow(12n).log(b)).toBe(0n);
    expect(b.pow(17n).log(b)).toBe(5n);
  });

  test('four prime factors', () => {
    const n = 3n * 5n * 7n * 11n;
    const b = Mod(2n, n);
    const ord = b.multiplicative_order();
    for (let e = 0n; e < ord; e++) {
      expect(b.pow(e).log(b)).toBe(e);
    }
  });

  test('SageMath doctests', () => {
    // sage: Mod(5, 9).log(Mod(2, 9))  ->  5
    expect(Mod(5n, 9n).log(Mod(2n, 9n))).toBe(5n);
    // sage: GF(7)(5).log()  ->  5   (default base = multiplicative generator 3)
    expect(Mod(5n, 7n).log()).toBe(5n);
    // sage: Mod(3, 7).log(Mod(2, 7))
    // ValueError: no logarithm of 3 found to base 2 modulo 7
    expect(() => Mod(3n, 7n).log(Mod(2n, 7n))).toThrow(
      'no logarithm of 3 found to base 2 modulo 7'
    );
    // sage: Mod(16, 100).log(Mod(4, 100))
    // ValueError: logarithm of 16 is not defined since it is not a unit modulo 100
    expect(() => Mod(16n, 100n).log(Mod(4n, 100n))).toThrow(
      'logarithm of 16 is not defined since it is not a unit modulo 100'
    );
    // sage: Mod(1111, 1234567).log(1111^3)
    // ValueError: ... (no solution modulo 9721)
    expect(() => Mod(1111n, 1234567n).log(Mod(1111n ** 3n, 1234567n))).toThrow(
      'no logarithm of 1111 found to base 961261 modulo 1234567 (no solution modulo 9721)'
    );
    // sage: Mod(230, 323).log(173)
    // ValueError: ... (incompatible local solutions)
    expect(() => Mod(230n, 323n).log(Mod(173n, 323n))).toThrow(
      'no logarithm of 230 found to base 173 modulo 323 (incompatible local solutions)'
    );
  });

  test('local solvability guard is not a tautology (L33)', () => {
    // Sage's precondition is `ord(self) | ord(base)`: modulo 9721 the base
    // 1111^3 has order 3240 while 1111 itself has order 9720, so 1111 is not a
    // power of it.  The old guard `ord % gcd(ord, ordb) != 0` was always false
    // and never fired; the failure only surfaced later, from discrete_log.
    const a = Mod(1111n, 9721n);
    const b = Mod(1111n ** 3n, 9721n);
    expect(a.multiplicative_order()).toBe(9720n);
    expect(b.multiplicative_order()).toBe(3240n);
    expect(3240n % 9720n).not.toBe(0n);
    expect(() => a.log(b)).toThrow('no logarithm of 1111 found to base 8603 modulo 9721');
  });

  test('non-unit base is rejected', () => {
    expect(() => Mod(3n, 10n).log(Mod(2n, 10n))).toThrow(
      'logarithm with base 2 is not defined since it is not a unit modulo 10'
    );
  });

  test('check option verifies the claimed order', () => {
    // sage: t.log(t, 57, check=True)
    // ValueError: base does not have the provided order
    const t = Mod(3n, 127n); // a primitive element of GF(127), order 126
    expect(t.multiplicative_order()).toBe(126n);
    expect(() => t.log(t, 57n, { check: true })).toThrow('base does not have the provided order');
    expect(t.log(t, 126n, { check: true })).toBe(1n);
  });
});

describe('multiplicative generators (M25)', () => {
  test('cyclicity test matches SageMath', () => {
    // integer_mod_ring.py:810-846 doctests
    expect(multiplicative_group_is_cyclic(7n)).toBe(true);
    expect(multiplicative_group_is_cyclic(9n)).toBe(true);
    expect(multiplicative_group_is_cyclic(8n)).toBe(false);
    expect(multiplicative_group_is_cyclic(4n)).toBe(true);
    expect(multiplicative_group_is_cyclic(75n)).toBe(false);
    expect(multiplicative_group_is_cyclic(162n)).toBe(true);
    expect(multiplicative_group_is_cyclic(1n << 20n)).toBe(false);
  });

  test('generators match SageMath', () => {
    // sage: Integers(7).multiplicative_generator()  ->  3
    expect(multiplicative_generator(7n)).toBe(3n);
    // sage: Integers(9).multiplicative_generator()  ->  2
    expect(multiplicative_generator(9n)).toBe(2n);
    // sage: Integers(4).multiplicative_generator()  ->  3
    expect(multiplicative_generator(4n)).toBe(3n);
    // sage: Integers(8).multiplicative_generator()
    // ValueError: multiplicative group of this ring is not cyclic
    expect(() => multiplicative_generator(8n)).toThrow(
      'multiplicative group of this ring is not cyclic'
    );
    expect(() => multiplicative_generator(75n)).toThrow(
      'multiplicative group of this ring is not cyclic'
    );
  });

  test('unit_gens matches SageMath', () => {
    // sage: Integers(25*3).unit_gens()  ->  (26, 52)
    expect(unit_gens(75n).map(([g]) => g)).toEqual([26n, 52n]);
    // sage: Integers(162).unit_gens()  ->  (83,)
    expect(unit_gens(162n).map(([g]) => g)).toEqual([83n]);
    // sage: _unit_gens_primepowercase(2, 3)  ->  [(7, 2), (5, 2)]
    expect(unit_gens(8n)).toEqual([
      [7n, 2n],
      [5n, 2n],
    ]);
    // generators really generate: every unit is a product of powers
    const gens = unit_gens(75n);
    const seen = new Set<bigint>();
    for (let i = 0n; i < gens[0]![1]!; i++) {
      for (let j = 0n; j < gens[1]![1]!; j++) {
        seen.add(Mod(gens[0]![0], 75n).pow(i).mul(Mod(gens[1]![0], 75n).pow(j)).value);
      }
    }
    expect(seen.size).toBe(40); // phi(75)
  });

  test('log without a base reports non-cyclic groups (M25)', () => {
    // sage: Mod(3, 16).log()
    // ValueError: multiplicative group of this ring is not cyclic
    expect(() => Mod(3n, 16n).log()).toThrow('multiplicative group of this ring is not cyclic');
    expect(() => Mod(3n, 16n).log()).toThrow(ValueError);
  });

  test('log without a base is fast for large non-cyclic moduli', () => {
    // The old residue scan computed a full multiplicative order for every one
    // of the 2^19 units before giving up (6.4 s).
    const start = Date.now();
    expect(() => Mod(3n, 1n << 20n).log()).toThrow(
      'multiplicative group of this ring is not cyclic'
    );
    expect(Date.now() - start).toBeLessThan(1000);
  });

  test('log without a base works for cyclic prime-power moduli', () => {
    const b = Mod(multiplicative_generator(2187n), 2187n); // 3^7
    for (const e of [0n, 1n, 5n, 100n, 1457n]) {
      expect(b.pow(b.pow(e).log()).value).toBe(b.pow(e).value);
    }
  });
});
