/**
 * Tests for sage/crypto/lattice
 */

import { describe, expect, test } from 'bun:test';
import { gen_lattice } from './lattice.js';

describe('gen_lattice', () => {
  test('generates modular lattice with default parameters', () => {
    const B = gen_lattice({ seed: 42 });

    // Default: n=4, m=8, q=11
    expect(B.length).toBe(8);
    expect(B[0]!.length).toBe(8);

    // First n rows should have q on diagonal
    for (let i = 0; i < 4; i++) {
      expect(B[i]![i]).toBe(11n);
      for (let j = 0; j < 8; j++) {
        if (j !== i) {
          expect(B[i]![j]).toBe(0n);
        }
      }
    }

    // Last (m-n) rows should have identity in the right part
    for (let i = 4; i < 8; i++) {
      expect(B[i]![i]).toBe(1n);
    }
  });

  test('generates random lattice with n=1', () => {
    const B = gen_lattice({ type: 'random', n: 1n, m: 10n, q: 14641n, seed: 42 });

    expect(B.length).toBe(10);
    expect(B[0]![0]).toBe(14641n); // q on diagonal
  });

  test('throws for random lattice with n != 1', () => {
    expect(() => gen_lattice({ type: 'random', n: 2n })).toThrow('random bases require n = 1');
  });

  test('generates dual lattice', () => {
    const B = gen_lattice({ m: 10n, seed: 42, dual: true });

    expect(B.length).toBe(10);
    expect(B[0]!.length).toBe(10);

    // For dual, the structure is different - q should appear in last columns
    // The matrix is reversed vertically
  });

  test('custom q parameter', () => {
    const q = 17n;
    const B = gen_lattice({ q, seed: 42 });

    // First row should have q on diagonal
    expect(B[0]![0]).toBe(q);
  });

  test('reproducible with seed', () => {
    const B1 = gen_lattice({ m: 6n, seed: 123 });
    const B2 = gen_lattice({ m: 6n, seed: 123 });

    for (let i = 0; i < B1.length; i++) {
      for (let j = 0; j < B1[i]!.length; j++) {
        expect(B1[i]![j]).toBe(B2[i]![j]);
      }
    }
  });

  test('different seeds give different lattices', () => {
    const B1 = gen_lattice({ m: 6n, seed: 123 });
    const B2 = gen_lattice({ m: 6n, seed: 456 });

    // At least some elements in the random part should differ
    let foundDiff = false;
    for (let i = 4; i < B1.length && !foundDiff; i++) {
      for (let j = 0; j < 4 && !foundDiff; j++) {
        if (B1[i]![j] !== B2[i]![j]) {
          foundDiff = true;
        }
      }
    }
    expect(foundDiff).toBe(true);
  });

  test('throws for ntl=true and lattice=true', () => {
    expect(() => gen_lattice({ ntl: true, lattice: true })).toThrow(
      'Cannot specify ntl=True and lattice=True at the same time'
    );
  });

  test('cyclotomic lattice with n=4', () => {
    // euler_phi(5) = 4, so n=4 should work
    const B = gen_lattice({ type: 'cyclotomic', n: 4n, m: 8n, seed: 42 });

    expect(B.length).toBe(8);
    expect(B[0]![0]).toBe(11n); // q on diagonal for first 4 rows
  });

  test('elements are within expected range', () => {
    const q = 11n;
    const B = gen_lattice({ q, seed: 42 });

    // Check random elements are in [-(q-1)/2, (q-1)/2]
    const halfQ = (q - 1n) / 2n;
    for (let i = 4; i < B.length; i++) {
      for (let j = 0; j < 4; j++) {
        expect(B[i]![j]!).toBeGreaterThanOrEqual(-halfQ);
        expect(B[i]![j]!).toBeLessThanOrEqual(halfQ);
      }
    }
  });

  test('lattice determinant matches expected value', () => {
    // For primal lattice: det(L) = q^n
    // We can verify this by checking the structure
    const n = 4n;
    const q = 11n;
    const B = gen_lattice({ n, q, m: 8n, seed: 42 });

    // First n rows have q on diagonal
    // The matrix is lower triangular, so det = product of diagonal elements
    // = q^n * 1^(m-n) = q^n = 11^4 = 14641
    let diagonalProduct = 1n;
    for (let i = 0; i < B.length; i++) {
      diagonalProduct *= B[i]![i]!;
    }
    expect(diagonalProduct).toBe(14641n); // 11^4
  });
});
