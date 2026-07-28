/**
 * Tests for sage/crypto/lattice
 */

import { describe, expect, test } from 'bun:test';
import { cyclotomic_polynomial } from '../rings/finite_rings/roots_of_unity.js';
import { gen_lattice } from './lattice.js';

/**
 * Exact (fraction-free Bareiss) determinant of an integer matrix.
 */
function det(M: bigint[][]): bigint {
  const A = M.map((r) => [...r]);
  const n = A.length;
  let sign = 1n;
  let prev = 1n;
  for (let k = 0; k < n - 1; k++) {
    if (A[k]![k] === 0n) {
      let s = -1;
      for (let i = k + 1; i < n; i++) {
        if (A[i]![k] !== 0n) {
          s = i;
          break;
        }
      }
      if (s === -1) return 0n;
      [A[k], A[s]] = [A[s]!, A[k]!];
      sign = -sign;
    }
    for (let i = k + 1; i < n; i++) {
      for (let j = k + 1; j < n; j++) {
        A[i]![j] = (A[i]![j]! * A[k]![k]! - A[i]![k]! * A[k]![j]!) / prev;
      }
    }
    prev = A[k]![k]!;
  }
  return sign * A[n - 1]![n - 1]!;
}

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
    // Sage builds block_matrix([[1, -A'^T], [0, q]]) with A' of shape
    // (m-n) x n, then reverses the rows.  So the identity block is I_n and the
    // q block is q*I_{m-n}: for m=10, n=4, q=11 exactly six rows contain 11
    // (Sage's doctest shows six) and det(B) = q^(m-n) = 11^6.
    const m = 10;
    const n = 4;
    const q = 11n;
    const B = gen_lattice({ m: 10n, seed: 42, dual: true });

    expect(B.length).toBe(m);
    expect(B[0]!.length).toBe(m);

    const rowsWithQ = B.filter((row) => row.some((v) => v === q || v === -q)).length;
    expect(rowsWithQ).toBe(m - n);

    const d = det(B);
    expect(d < 0n ? -d : d).toBe(q ** BigInt(m - n));
  });

  test('dual basis is the exact q-dual of the primal basis', () => {
    // Sage's documented relation:
    //   B_dual.hermite_form() == transpose(q*B_primal.inverse()).hermite_form()
    // With B_primal = [[q, 0], [A', 1]] we have q*B_primal^{-T} = [[1, -A'^T],
    // [0, q]], which is exactly the (row-reversed) dual basis.
    const m = 10;
    const n = 4;
    const q = 11n;
    const P = gen_lattice({ m: 10n, seed: 42 });
    const D = gen_lattice({ m: 10n, seed: 42, dual: true });

    const Aprime = P.slice(n).map((row) => row.slice(0, n)); // (m-n) x n
    const E = [...D].reverse();

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < m; j++) {
        const want = j < n ? (i === j ? 1n : 0n) : -Aprime[j - n]![i]!;
        expect(E[i]![j]).toBe(want);
      }
    }
    for (let i = 0; i < m - n; i++) {
      for (let j = 0; j < m; j++) {
        const want = j < n ? 0n : i === j - n ? q : 0n;
        expect(E[n + i]![j]).toBe(want);
      }
    }
  });

  test('primal determinant is q^n', () => {
    const d = det(gen_lattice({ m: 10n, seed: 42 }));
    expect(d < 0n ? -d : d).toBe(11n ** 4n);
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

  test('cyclotomic blocks are multiplication matrices modulo Phi_k, not x^n+1', () => {
    // n = 6 -> k = 9 and Phi_9 = x^6 + x^3 + 1, which is NOT x^6 + 1.  Row j of
    // each stacked block must be the coefficient vector of x^j * a mod Phi_9.
    const n = 6;
    const q = 11n;
    const B = gen_lattice({ type: 'cyclotomic', n: 6n, m: 12n, q, seed: 7 });
    const phi9 = (cyclotomic_polynomial(9) as number[]).map(BigInt);
    expect(phi9).toEqual([1n, 0n, 0n, 1n, 0n, 0n, 1n]);

    const mod = (a: bigint) => ((a % q) + q) % q;
    const mulX = (c: bigint[]): bigint[] => {
      const lead = c[n - 1]!;
      const out = new Array<bigint>(n);
      out[0] = mod(-lead * phi9[0]!);
      for (let i = 1; i < n; i++) out[i] = mod(c[i - 1]! - lead * phi9[i]!);
      return out;
    };

    let cur = B[n]!.slice(0, n).map(mod); // the random element a
    for (let row = 0; row < n; row++) {
      expect(B[n + row]!.slice(0, n).map(mod)).toEqual(cur);
      cur = mulX(cur);
    }
    // The pure negacyclic shift (x^6 = -1) would give row 1 = (-a5, a0, a1, a2,
    // a3, a4).  Modulo Phi_9 the x^3 coefficient is a2 - a5 instead of a2, so
    // the two differ whenever a5 != 0.
    const a = B[n]!.slice(0, n).map(mod);
    expect(a[5]).not.toBe(0n);
    expect(B[n + 1]!.slice(0, n).map(mod)[3]).toBe(mod(a[2]! - a[5]!));
    expect(B[n + 1]!.slice(0, n).map(mod)[3]).not.toBe(a[2]!);
  });

  test('ideal lattice with quotient x^4 - 1 is a circulant NTRU basis', () => {
    // sage: sage.crypto.gen_lattice(type='ideal', seed=42, quotient=x^4 - 1)
    const n = 4;
    const B = gen_lattice({ type: 'ideal', seed: 42, quotient: [-1n, 0n, 0n, 0n, 1n] });
    expect(B.length).toBe(8);
    for (let i = 0; i < n; i++) {
      expect(B[i]![i]).toBe(11n);
    }
    // x^4 = 1 makes each block a circulant matrix.
    const a = B[n]!.slice(0, n);
    for (let row = 0; row < n; row++) {
      for (let col = 0; col < n; col++) {
        expect(B[n + row]![col]).toBe(a[(col - row + n) % n]!);
      }
    }
  });

  test('ideal lattice rejects a missing or wrong-degree quotient', () => {
    expect(() => gen_lattice({ type: 'ideal', seed: 1234 })).toThrow(
      'ideal bases require a quotient polynomial'
    );
    // sage: gen_lattice(type='ideal', seed=1234, quotient=x^23-1)
    // ValueError: ideal basis requires n = quotient.degree()
    const x23m1 = [-1n, ...Array<bigint>(22).fill(0n), 1n];
    expect(() => gen_lattice({ type: 'ideal', seed: 1234, quotient: x23m1 })).toThrow(
      'ideal basis requires n = quotient.degree()'
    );
  });

  test('ntl and lattice output flags', () => {
    // Sage returns B._ntl_(), printed in the bracketed NTL matrix format.
    const s = gen_lattice({ m: 10n, q: 11n, seed: 42, ntl: true });
    expect(typeof s).toBe('string');
    expect(s.startsWith('[\n[')).toBe(true);
    expect(s.endsWith(']\n]')).toBe(true);
    expect(s.split('\n').length).toBe(12); // '[' + 10 rows + ']'

    // Sage returns IntegerLattice(B).
    const L = gen_lattice({ m: 10n, q: 11n, seed: 42, lattice: true });
    expect(L.rank()).toBe(10);
    expect(L.degree()).toBe(10);
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
