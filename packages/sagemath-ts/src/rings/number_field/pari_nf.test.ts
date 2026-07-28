/**
 * Tests for the number-field kernel routines (the port of PARI's nfbasis,
 * nfdisc, idealprimedec support and nfgaloisconj).
 *
 * Reference values are PARI/GP outputs quoted in the SageMath doctests.
 */

import { describe, expect, it } from 'bun:test';
import { Rational } from '../rational.js';
import {
  type ZPoly,
  fpFactor,
  fpRoots,
  hnf,
  hnfLower,
  integralDefiningPolynomial,
  nfbasis,
  nfdisc,
  nfgaloisconj,
  rationalReconstruct,
  zpDiscriminant,
  zpFactorSquarefree,
  zpIsIrreducibleOverQ,
  zpResultant,
} from './pari_nf.js';

describe('integer polynomial arithmetic', () => {
  it('computes discriminants', () => {
    // PARI: poldisc(x^3 + x^2 - 2*x + 8) = -2012
    expect(zpDiscriminant([8n, -2n, 1n, 1n])).toBe(-2012n);
    expect(zpDiscriminant([-2n, 0n, 1n])).toBe(8n);
    expect(zpDiscriminant([1n, 0n, 1n])).toBe(-4n);
    expect(zpDiscriminant([-5n, 0n, 1n])).toBe(20n);
    // PARI: poldisc(x^5 - 1) = 125 * ... ; use the cyclotomic Phi_5
    expect(zpDiscriminant([1n, 1n, 1n, 1n, 1n])).toBe(125n);
    // PARI: poldisc(x^3 - 2) = -108
    expect(zpDiscriminant([-2n, 0n, 0n, 1n])).toBe(-108n);
  });

  it('computes resultants', () => {
    // PARI: polresultant(x^2 - 2, x^2 - 3) = 1
    expect(zpResultant([-2n, 0n, 1n], [-3n, 0n, 1n])).toBe(1n);
    // PARI: polresultant(x^2 - 1, x - 1) = 0
    expect(zpResultant([-1n, 0n, 1n], [-1n, 1n])).toBe(0n);
    // res(x^3 + x + 1, x^2 + 1) = f(i) * f(-i) with f = x^3 + x + 1:
    // f(i) = -i + i + 1 = 1 and f(-i) = i - i + 1 = 1, so the resultant is 1.
    expect(zpResultant([1n, 1n, 0n, 1n], [1n, 0n, 1n])).toBe(1n);
    // res(x^2 + 1, x^2 + 4) = (i^2 + 4)(( -i)^2 + 4) = 3 * 3 = 9
    expect(zpResultant([1n, 0n, 1n], [4n, 0n, 1n])).toBe(9n);
  });
});

describe('factorisation over F_p', () => {
  it('splits x^2 + 1 mod 5', () => {
    const f = fpFactor([1n, 0n, 1n], 5n);
    expect(f.length).toBe(2);
    expect(f.map(([g]) => g)).toEqual([
      [2n, 1n],
      [3n, 1n],
    ]);
  });

  it('finds the repeated factor of x^2 - 2 mod 2', () => {
    const f = fpFactor([-2n, 0n, 1n], 2n);
    expect(f).toEqual([[[0n, 1n], 2]]);
  });

  it('splits x^3 - 2 mod 5 into a linear and a quadratic factor', () => {
    const f = fpFactor([-2n, 0n, 0n, 1n], 5n);
    expect(f.map(([g]) => g.length - 1).sort()).toEqual([1, 2]);
  });

  it('recognises x^3 - 2 as irreducible mod 7', () => {
    // 2 is not a cube mod 7, so x^3 - 2 has no root and is irreducible
    expect(fpRoots([-2n, 0n, 0n, 1n], 7n)).toEqual([]);
    expect(fpFactor([-2n, 0n, 0n, 1n], 7n).length).toBe(1);
  });

  it('splits x^3 - 2 completely mod 31', () => {
    const roots = fpRoots([-2n, 0n, 0n, 1n], 31n);
    expect(roots.length).toBe(3);
    for (const r of roots) {
      expect((r * r * r - 2n) % 31n).toBe(0n);
    }
  });

  it('reconstructs the polynomial from its factors', () => {
    const p = 13n;
    const f: ZPoly = [6n, 5n, 4n, 3n, 2n, 1n];
    let prod: bigint[] = [1n];
    for (const [g, e] of fpFactor(f, p)) {
      for (let i = 0; i < e; i++) {
        const r: bigint[] = new Array(prod.length + g.length - 1).fill(0n);
        for (let a = 0; a < prod.length; a++) {
          for (let b = 0; b < g.length; b++) {
            r[a + b] = (r[a + b]! + prod[a]! * g[b]!) % p;
          }
        }
        prod = r;
      }
    }
    // f is monic already (leading coefficient 1)
    expect(prod.map((c) => ((c % p) + p) % p)).toEqual(f.map((c) => ((c % p) + p) % p));
  });
});

describe('factorisation over Z', () => {
  it('factors x^4 + 4 as (x^2 - 2x + 2)(x^2 + 2x + 2)', () => {
    const factors = zpFactorSquarefree([4n, 0n, 0n, 0n, 1n]);
    expect(factors.length).toBe(2);
    const sorted = factors.map((f) => f.join(',')).sort();
    expect(sorted).toEqual(['2,-2,1', '2,2,1']);
  });

  it('decides irreducibility over Q', () => {
    expect(zpIsIrreducibleOverQ([-2n, 0n, 1n])).toBe(true);
    expect(zpIsIrreducibleOverQ([-1n, 0n, 1n])).toBe(false);
    expect(zpIsIrreducibleOverQ([-4n, 0n, 1n])).toBe(false);
    expect(zpIsIrreducibleOverQ([1n, 0n, 2n, 0n, 1n])).toBe(false); // (x^2+1)^2
    expect(zpIsIrreducibleOverQ([2n, 0n, 3n, 0n, 1n])).toBe(false); // (x^2+1)(x^2+2)
    expect(zpIsIrreducibleOverQ([1n, 0n, 0n, 0n, 1n])).toBe(true); // x^4 + 1
    expect(zpIsIrreducibleOverQ([-2n, 0n, 0n, 1n])).toBe(true); // x^3 - 2
    expect(zpIsIrreducibleOverQ([1n, 1n, 1n, 1n, 1n, 1n, 1n])).toBe(true); // Phi_7
  });
});

describe('nfbasis / nfdisc', () => {
  it('matches PARI on the Sage doctest x^3 + x^2 - 2x + 8', () => {
    const r = nfbasis([8n, -2n, 1n, 1n]);
    expect(r.disc).toBe(-503n);
    expect(r.index).toBe(2n);
    expect(r.den).toBe(2n);
    // PARI: nfbasis(x^3 + x^2 - 2*x + 8) = [1, x, 1/2*x^2 + 1/2*x]
    expect(r.basis).toEqual([
      [2n, 0n, 0n],
      [0n, 2n, 0n],
      [0n, 1n, 1n],
    ]);
  });

  it('matches PARI on quadratic fields', () => {
    expect(nfdisc([-5n, 0n, 1n])).toBe(5n);
    expect(nfdisc([-2n, 0n, 1n])).toBe(8n);
    expect(nfdisc([1n, 0n, 1n])).toBe(-4n);
    expect(nfdisc([3n, 0n, 1n])).toBe(-3n);
    expect(nfdisc([-8n, 0n, 1n])).toBe(8n);
    expect(nfdisc([30n, 0n, 1n])).toBe(-120n);
    expect(nfdisc([9n, 0n, 1n])).toBe(-4n);
    expect(nfdisc([27n, 0n, 1n])).toBe(-3n);
    expect(nfbasis([-5n, 0n, 1n]).basis).toEqual([
      [2n, 0n],
      [1n, 1n],
    ]);
  });

  it('matches PARI on cyclotomic fields', () => {
    // disc(Q(zeta_p)) = (-1)^((p-1)/2) p^(p-2)
    expect(nfdisc([1n, 1n, 1n, 1n, 1n])).toBe(125n); // p = 5
    expect(nfdisc([1n, 1n, 1n, 1n, 1n, 1n, 1n])).toBe(-16807n); // p = 7
    expect(nfdisc([1n, 0n, 0n, 0n, 1n])).toBe(256n); // Q(zeta_8)
  });

  it('matches PARI on cubic fields', () => {
    expect(nfdisc([-2n, 0n, 0n, 1n])).toBe(-108n);
    expect(nfdisc([1n, -3n, 0n, 1n])).toBe(81n);
    // x^3 - x - 1 has squarefree discriminant -23
    expect(nfdisc([-1n, -1n, 0n, 1n])).toBe(-23n);
  });

  it('satisfies disc(poly) = index^2 * disc(field)', () => {
    const polys: ZPoly[] = [
      [8n, -2n, 1n, 1n],
      [-5n, 0n, 1n],
      [-8n, 0n, 1n],
      [30n, 0n, 1n],
      [-2n, 0n, 0n, 1n],
      [1n, 1n, 1n, 1n, 1n],
      [-12n, 0n, 0n, 1n],
    ];
    for (const g of polys) {
      const r = nfbasis(g);
      expect(zpDiscriminant(g)).toBe(r.index * r.index * r.disc);
    }
  });
});

describe('integralDefiningPolynomial', () => {
  it('clears denominators of a monic rational polynomial', () => {
    const { g, scale } = integralDefiningPolynomial([
      new Rational(-1n, 2n),
      Rational.zero(),
      Rational.one(),
    ]);
    expect(scale).toBe(2n);
    expect(g).toEqual([-2n, 0n, 1n]);
  });

  it('leaves integral polynomials alone', () => {
    const { g, scale } = integralDefiningPolynomial([
      new Rational(-2n),
      Rational.zero(),
      Rational.one(),
    ]);
    expect(scale).toBe(1n);
    expect(g).toEqual([-2n, 0n, 1n]);
  });

  it('handles a cubic with several denominators', () => {
    // x^3 + (1/3)x + 1/9 -> y = 3x gives y^3 + 3y + 3
    const { g, scale } = integralDefiningPolynomial([
      new Rational(1n, 9n),
      new Rational(1n, 3n),
      Rational.zero(),
      Rational.one(),
    ]);
    expect(scale).toBe(3n);
    expect(g).toEqual([3n, 3n, 0n, 1n]);
  });
});

describe('nfgaloisconj', () => {
  it('finds the three conjugates of the cyclic cubic x^3 - 3x + 1', () => {
    const conj = nfgaloisconj([1n, -3n, 0n, 1n]);
    expect(conj.length).toBe(3);
    // PARI: nfgaloisconj(x^3-3*x+1) = [-x^2-x+2, x^2-2, x]
    const asStrings = conj.map((c) => c.map((r) => r.toString()).join(',')).sort();
    expect(asStrings).toEqual(['0,1,0', '2,-1,-1', '-2,0,1'].sort());
  });

  it('returns only the identity for the non-Galois x^3 - 2', () => {
    expect(nfgaloisconj([-2n, 0n, 0n, 1n]).length).toBe(1);
  });

  it('handles quadratic and quartic fields', () => {
    expect(nfgaloisconj([-2n, 0n, 1n]).length).toBe(2);
    // x^4 + 1 is the 8th cyclotomic polynomial: Galois with group C2 x C2
    expect(nfgaloisconj([1n, 0n, 0n, 0n, 1n]).length).toBe(4);
    // x^4 - 2 is not Galois; only a |-> +/- a lie in the field
    expect(nfgaloisconj([-2n, 0n, 0n, 0n, 1n]).length).toBe(2);
  });
});

describe('lattice helpers', () => {
  it('computes an upper-triangular HNF', () => {
    const H = hnf(
      [
        [2n, 0n],
        [1n, 1n],
      ],
      2
    );
    expect(H).toEqual([
      [1n, 1n],
      [0n, 2n],
    ]);
  });

  it('computes a lower-triangular HNF', () => {
    const H = hnfLower(
      [
        [2n, 0n],
        [1n, 1n],
      ],
      2
    );
    expect(H).toEqual([
      [2n, 0n],
      [1n, 1n],
    ]);
  });

  it('reconstructs rationals from residues', () => {
    const m = 1000003n;
    // 3/7 mod m, with 1/7 built directly from a multiple of m
    let inv7 = 0n;
    for (let k = 1n; k < 7n; k++) {
      if ((k * m + 1n) % 7n === 0n) {
        inv7 = (k * m + 1n) / 7n;
        break;
      }
    }
    const residue = (3n * inv7) % m;
    const q = rationalReconstruct(residue, m, 1000n, 1000n);
    expect(q).not.toBeNull();
    expect(q!.toString()).toBe('3/7');
  });
});
