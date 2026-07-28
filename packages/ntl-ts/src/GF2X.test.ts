/**
 * Tests for ntl-ts GF2X.
 *
 * Oracles:
 *  - NTL: ntl/src/GF2XFactoring.cpp (BuildIrred / BuildSparseIrred / IterIrredTest)
 *  - SageMath doctests: sage/rings/polynomial/polynomial_gf2x.pyx:298-341
 *      GF2X_BuildIrred_list(2) == [1, 1, 1]
 *      GF2X_BuildIrred_list(3) == [1, 1, 0, 1]
 *      GF2X_BuildIrred_list(4) == [1, 1, 0, 0, 1]
 *      GF(2)['x'](GF2X_BuildIrred_list(33)) == x^33 + x^6 + x^3 + x + 1
 *      GF2X_BuildSparseIrred_list(n) == GF2X_BuildIrred_list(n) for n in 1..32
 *      GF(2)['x'](GF2X_BuildSparseIrred_list(33)) == x^33 + x^10 + 1
 */

import { describe, expect, it } from 'bun:test';
import { GF2X } from './GF2X.js';
import { GF2 } from './GF2.js';

/** Sage's GF2X_BuildIrred_list: coefficients of BuildIrred(n), constant first. */
function BuildIrred_list(n: number): number[] {
  const f = GF2X.BuildIrred(n);
  return Array.from({ length: n + 1 }, (_, i) => f.coeff(i).rep());
}

/** Sage's GF2X_BuildSparseIrred_list. */
function BuildSparseIrred_list(n: number): number[] {
  const f = GF2X.BuildSparseIrred(n);
  return Array.from({ length: n + 1 }, (_, i) => f.coeff(i).rep());
}

describe('GF2 basics', () => {
  it('implements NTL arithmetic', () => {
    expect(new GF2(1).add(new GF2(1)).rep()).toBe(0);
    expect(new GF2(1).mul(new GF2(1)).rep()).toBe(1);
    expect(new GF2(3).rep()).toBe(1); // GF2(INIT_VAL, a) keeps a & 1
    expect(new GF2(1).div(new GF2(1)).rep()).toBe(1);
    expect(() => new GF2(1).div(new GF2(0))).toThrow('GF2: division by zero');
    expect(new GF2(0).power(0).rep()).toBe(1); // NTL GF2.cpp:10
    expect(() => new GF2(0).power(-1)).toThrow('GF2: division by zero');
    expect(GF2.conv(4n).rep()).toBe(0);
    expect(GF2.conv(7n).rep()).toBe(1);
  });
});

describe('GF2X arithmetic', () => {
  it('adds, multiplies, squares and divides', () => {
    const a = new GF2X([1, 1, 1]); // x^2 + x + 1
    const b = new GF2X([1, 1]); // x + 1
    expect(a.add(b).rep()).toBe(0b100n); // x^2
    expect(a.mul(b).rep()).toBe(0b1001n); // x^3 + 1
    expect(a.sqr().rep()).toBe(0b10101n); // x^4 + x^2 + 1
    expect(a.sqr().equals(a.mul(a))).toBe(true);
    const [q, r] = a.mul(b).DivRem(b);
    expect(q.equals(a)).toBe(true);
    expect(r.IsZero()).toBe(true);
    expect(a.deg()).toBe(2);
    expect(GF2X.zero().deg()).toBe(-1);
    expect(() => a.DivRem(GF2X.zero())).toThrow('GF2X: division by zero');
  });

  it('computes GCD and XGCD with s*a + t*b = d', () => {
    const a = new GF2X([1, 0, 1, 1, 0, 1]); // x^5 + x^3 + x^2 + 1
    const b = new GF2X([1, 1, 0, 1]); // x^3 + x + 1
    const d = GF2X.GCD(a, b);
    const [d2, s, t] = GF2X.XGCD(a, b);
    expect(d2.equals(d)).toBe(true);
    expect(s.mul(a).add(t.mul(b)).equals(d)).toBe(true);
    // XGCD(a, 0) == (a, 1, 0)  (NTL BaseXGCD)
    const [d3, s3, t3] = GF2X.XGCD(a, GF2X.zero());
    expect(d3.equals(a)).toBe(true);
    expect(s3.IsOne()).toBe(true);
    expect(t3.IsZero()).toBe(true);
  });

  it('does modular arithmetic against a GF(2^8) modulus', () => {
    // AES modulus x^8 + x^4 + x^3 + x + 1
    const f = new GF2X(0b100011011n);
    expect(f.isIrreducible()).toBe(true);
    const a = new GF2X(0b10000011n);
    const inv = GF2X.InvMod(a, f);
    expect(GF2X.MulMod(a, inv, f).IsOne()).toBe(true);
    // a^(2^8 - 1) = 1 in GF(2^8)^*
    expect(GF2X.PowerMod(a, 255n, f).IsOne()).toBe(true);
    expect(GF2X.PowerMod(a, -1n, f).equals(inv)).toBe(true);
    expect(GF2X.SqrMod(a, f).equals(GF2X.MulMod(a, a, f))).toBe(true);
    expect(() => GF2X.InvMod(f, f)).toThrow('InvMod: bad args');
    // x^2 + 1 = (x+1)^2 is not invertible modulo x+1... use a reducible modulus
    const g = new GF2X([1, 1]).mul(new GF2X([1, 1, 1])); // (x+1)(x^2+x+1)
    expect(() => GF2X.InvMod(new GF2X([1, 1]), g)).toThrow('InvMod: inverse undefined');
  });

  it('implements diff, reverse, shifts, trunc, hex', () => {
    const a = new GF2X([1, 1, 0, 1, 1]); // x^4 + x^3 + x + 1
    expect(a.diff().rep()).toBe(0b101n); // x^2 + 1
    expect(a.reverse().rep()).toBe(0b11011n);
    expect(a.LeftShift(2).rep()).toBe(a.rep() << 2n);
    expect(a.RightShift(2).rep()).toBe(a.rep() >> 2n);
    expect(a.trunc(3).rep()).toBe(0b011n);
    // NTL hex format packs coefficients four at a time, low nibble first
    expect(new GF2X([1, 1, 0, 0, 1]).toHex()).toBe('0x31'); // x^4 + x + 1
    expect(new GF2X([1, 1, 1]).toHex()).toBe('0x7');
    expect(GF2X.fromHex('0x31').rep()).toBe(0b10011n);
    expect(new GF2X([1, 1, 0, 1]).toString()).toBe('[1 1 0 1]');
  });
});

describe('GF2X.isIrreducible (NTL IterIrredTest)', () => {
  it('matches brute-force irreducibility for every polynomial of degree <= 10', () => {
    // Brute force: f (deg d >= 1) is irreducible iff no divisor of degree 1..d/2
    const isIrredBrute = (rep: bigint, d: number): boolean => {
      const f = new GF2X(rep);
      for (let e = 1; e <= d / 2; e++) {
        for (let g = 1n << BigInt(e); g < 1n << BigInt(e + 1); g++) {
          if (f.rem(new GF2X(g)).IsZero()) return false;
        }
      }
      return true;
    };
    let irredCount = 0;
    for (let d = 1; d <= 10; d++) {
      for (let low = 0n; low < 1n << BigInt(d); low++) {
        const rep = (1n << BigInt(d)) | low;
        const expected = isIrredBrute(rep, d);
        expect(new GF2X(rep).isIrreducible()).toBe(expected);
        if (expected) irredCount++;
      }
    }
    // Number of monic irreducibles of degree d over GF(2), d = 1..10:
    // 2, 1, 2, 3, 6, 9, 18, 30, 56, 99
    expect(irredCount).toBe(2 + 1 + 2 + 3 + 6 + 9 + 18 + 30 + 56 + 99);
  });

  it('rejects constants and zero', () => {
    expect(GF2X.zero().isIrreducible()).toBe(false);
    expect(GF2X.one().isIrreducible()).toBe(false);
    expect(GF2X.X().isIrreducible()).toBe(true);
  });
});

describe('GF2X.BuildIrred (H120)', () => {
  it("matches Sage's GF2X_BuildIrred_list doctests", () => {
    // sage/rings/polynomial/polynomial_gf2x.pyx:306-313
    expect(BuildIrred_list(2)).toEqual([1, 1, 1]);
    expect(BuildIrred_list(3)).toEqual([1, 1, 0, 1]);
    expect(BuildIrred_list(4)).toEqual([1, 1, 0, 0, 1]);
    // GF(2)['x'](GF2X_BuildIrred_list(33)) == x^33 + x^6 + x^3 + x + 1
    expect(GF2X.BuildIrred(33).rep()).toBe(
      (1n << 33n) | (1n << 6n) | (1n << 3n) | (1n << 1n) | 1n,
    );
  });

  it('returns X for n = 1 (NTL GF2XFactoring.cpp:479)', () => {
    expect(GF2X.BuildIrred(1).rep()).toBe(2n);
    expect(BuildIrred_list(1)).toEqual([0, 1]);
  });

  it('is the lexicographically smallest irreducible of its degree', () => {
    for (let n = 1; n <= 12; n++) {
      const f = GF2X.BuildIrred(n);
      expect(f.deg()).toBe(n);
      expect(f.isIrreducible()).toBe(true);
      // NTL enumerates x^n + ConvertBits(2i+1); nothing smaller can be irreducible
      const low = f.rep() ^ (1n << BigInt(n));
      for (let i = 1n; i < low; i += 2n) {
        expect(new GF2X((1n << BigInt(n)) | i).isIrreducible()).toBe(false);
      }
    }
  });

  it('rejects non-positive degrees', () => {
    expect(() => GF2X.BuildIrred(0)).toThrow('BuildIrred: n must be positive');
    expect(() => GF2X.BuildIrred(-1)).toThrow('BuildIrred: n must be positive');
  });
});

describe('GF2X.BuildSparseIrred (H120)', () => {
  it("matches Sage's GF2X_BuildSparseIrred_list doctests", () => {
    // all([GF2X_BuildSparseIrred_list(n) == GF2X_BuildIrred_list(n) for n in 1..32])
    for (let n = 1; n < 33; n++) {
      expect(BuildSparseIrred_list(n)).toEqual(BuildIrred_list(n));
    }
    // GF(2)['x'](GF2X_BuildSparseIrred_list(33)) == x^33 + x^10 + 1
    expect(GF2X.BuildSparseIrred(33).rep()).toBe((1n << 33n) | (1n << 10n) | 1n);
  });

  it('produces an irreducible polynomial of the requested degree for n in [1, 200]', () => {
    for (let n = 1; n <= 200; n++) {
      const f = GF2X.BuildSparseIrred(n);
      expect(f.deg()).toBe(n);
      expect(f.isIrreducible()).toBe(true);
    }
  });

  it('produces minimal-weight polynomials (trinomial when one exists)', () => {
    const weight = (f: GF2X): number => {
      let w = 0;
      for (let i = 0; i <= f.deg(); i++) w += f.coeff(i).rep();
      return w;
    };
    // The table is NTL's; spot-check that a trinomial is used exactly when one
    // is irreducible, by an independent search.
    for (let n = 2; n <= 60; n++) {
      let trinom = 0;
      for (let k = 1; k < n; k++) {
        const g = new GF2X((1n << BigInt(n)) | (1n << BigInt(k)) | 1n);
        if (g.isIrreducible()) {
          trinom = k;
          break;
        }
      }
      const w = weight(GF2X.BuildSparseIrred(n));
      if (trinom) expect(w).toBe(3);
      else expect(w).toBe(5);
    }
  });

  it('rejects non-positive degrees', () => {
    expect(() => GF2X.BuildSparseIrred(0)).toThrow('SparseIrred: n <= 0');
  });
});
