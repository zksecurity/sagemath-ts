/**
 * Tests for the number-field kernel routines (the port of PARI's nfbasis,
 * nfdisc, idealprimedec support and nfgaloisconj).
 *
 * Reference values are PARI/GP outputs quoted in the SageMath doctests.
 */

import { describe, expect, it } from 'bun:test';
import { isqrt } from '../../arith/misc.js';
import { Rational } from '../rational.js';
import {
  type NumberField,
  NumberFieldConstructor,
  type NumberFieldElement,
} from './number_field.js';
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
  numberofconjugates,
  primedec,
  quadunit,
  quadunitnorm,
  ratInverse,
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

  // Audit H24: the port used to cap nfgaloisconj at degree 8 because it
  // enumerated all n! permutations of the p-adic roots.  The LLL
  // reconstruction removes the cap entirely.
  describe('beyond the old degree-8 cap', () => {
    /** Phi_n, by dividing x^n - 1 by the lower cyclotomics. */
    const cyclotomic = (n: number): ZPoly => {
      const divExact = (a: ZPoly, b: ZPoly): ZPoly => {
        const r = [...a];
        const q = new Array<bigint>(a.length - b.length + 1).fill(0n);
        for (let d = r.length - 1; d >= b.length - 1; d--) {
          const c = r[d]! / b[b.length - 1]!;
          q[d - b.length + 1] = c;
          for (let i = 0; i < b.length; i++) {
            r[d - b.length + 1 + i] = r[d - b.length + 1 + i]! - c * b[i]!;
          }
        }
        return q;
      };
      let num: ZPoly = new Array<bigint>(n + 1).fill(0n);
      num[0] = -1n;
      num[n] = 1n;
      for (let d = 1; d < n; d++) if (n % d === 0) num = divExact(num, cyclotomic(d));
      return num;
    };

    /** beta is a root of g in Q[x]/(g)? (independent of nfgaloisconj's own check) */
    const isRoot = (c: Rational[], g: ZPoly): boolean => {
      const n = g.length - 1;
      const mul = (a: Rational[], b: Rational[]): Rational[] => {
        const out: Rational[] = new Array(2 * n - 1).fill(Rational.zero());
        for (let i = 0; i < a.length; i++) {
          for (let j = 0; j < b.length; j++) out[i + j] = out[i + j]!.add(a[i]!.mul(b[j]!));
        }
        for (let d = out.length - 1; d >= n; d--) {
          const t = out[d]!;
          if (t.sign === 0n) continue;
          out[d] = Rational.zero();
          for (let i = 0; i < n; i++) {
            out[d - n + i] = out[d - n + i]!.sub(t.mul(new Rational(g[i]!)));
          }
        }
        return out.slice(0, n);
      };
      const acc: Rational[] = new Array(n).fill(Rational.zero());
      let pw: Rational[] = new Array(n).fill(Rational.zero());
      pw[0] = Rational.one();
      for (let k = 0; k <= n; k++) {
        const gk = new Rational(g[k] ?? 0n);
        for (let i = 0; i < n; i++) acc[i] = acc[i]!.add(gk.mul(pw[i]!));
        if (k < n) pw = mul(pw, c);
      }
      return acc.every((r) => r.sign === 0n);
    };

    // Q(zeta_n) is Galois of degree phi(n), so nfgaloisconj must find phi(n)
    // conjugates.  Degrees 10..20 were all rejected before the fix.
    for (const [n, phi] of [
      [11, 10],
      [13, 12],
      [15, 8],
      [16, 8],
      [17, 16],
      [21, 12],
    ] as Array<[number, number]>) {
      it(`finds all ${phi} conjugates of Q(zeta_${n})`, () => {
        const G = cyclotomic(n);
        expect(G.length - 1).toBe(phi);
        const conj = nfgaloisconj(G);
        expect(conj.length).toBe(phi);
        // every one is genuinely a root of G in the field, and they are distinct
        const keys = new Set<string>();
        for (const c of conj) {
          expect(isRoot(c, G)).toBe(true);
          keys.add(c.map((r) => r.toString()).join(','));
        }
        expect(keys.size).toBe(phi);
        // the identity comes first, as in PARI
        expect(conj[0]!.map((r) => r.toString()).join(',')).toBe(
          ['0', '1', ...Array(phi - 2).fill('0')].join(',')
        );
      });
    }

    it('still reports the right (small) automorphism group of non-Galois fields', () => {
      // x^n - 2 has as many automorphisms as it has real roots in Q(2^(1/n)):
      // 2 for even n (+/- the real root), 1 for odd n.
      const xn2 = (n: number): ZPoly => {
        const g = new Array<bigint>(n + 1).fill(0n);
        g[0] = -2n;
        g[n] = 1n;
        return g;
      };
      expect(nfgaloisconj(xn2(9)).length).toBe(1);
      expect(nfgaloisconj(xn2(10)).length).toBe(2);
      expect(nfgaloisconj(xn2(12)).length).toBe(2);
      // x^5 - x - 1 has Galois group S5, so Aut(K) is trivial
      expect(nfgaloisconj([-1n, -1n, 0n, 0n, 0n, 1n]).length).toBe(1);
    });
  });
});

describe('numberofconjugates (PARI galconj.c:3113)', () => {
  it('is an upper bound reached exactly on these fields', () => {
    // #Aut(K) divides the answer; here the bound is tight.
    expect(numberofconjugates([1n, -3n, 0n, 1n])).toBe(3n); // cyclic cubic
    expect(numberofconjugates([-2n, 0n, 0n, 1n])).toBe(1n); // x^3 - 2
    expect(numberofconjugates([-2n, 0n, 0n, 0n, 1n])).toBe(2n); // x^4 - 2
    expect(numberofconjugates([1n, 0n, 0n, 0n, 1n])).toBe(4n); // Q(zeta_8)
    expect(numberofconjugates([-1n, -1n, 0n, 0n, 0n, 1n])).toBe(1n); // S5
  });

  it('always dominates the true number of conjugates', () => {
    const polys: ZPoly[] = [
      [1n, -3n, 0n, 1n],
      [-2n, 0n, 0n, 1n],
      [-2n, 0n, 0n, 0n, 1n],
      [1n, 0n, 0n, 0n, 1n],
      [-1n, -1n, 0n, 0n, 0n, 1n],
      [8n, -2n, 1n, 1n],
      [-8n, -2n, -1n, 1n],
      [1n, 1n, 1n, 1n, 1n, 1n, 1n],
    ];
    for (const g of polys) {
      const m = BigInt(nfgaloisconj(g).length);
      const c = numberofconjugates(g);
      expect(c % m).toBe(0n);
    }
  });
});

describe('quadunit (PARI quad.c:281)', () => {
  it('reproduces PARI test-suite values', () => {
    // reference/pari/src/test/in/number:81 -> 32/number:238  quadunit(17) = 3 + 2*w
    expect(quadunit(17n)).toEqual([3n, 2n]);
    // classical fundamental units, epsilon = u + v*w_D
    expect(quadunit(5n)).toEqual([0n, 1n]); // (1+sqrt5)/2
    expect(quadunit(8n)).toEqual([1n, 1n]); // 1 + sqrt2
    expect(quadunit(12n)).toEqual([2n, 1n]); // 2 + sqrt3
    expect(quadunit(376n)).toEqual([2143295n, 221064n]); // 2143295 + 221064*sqrt94
  });

  it('reproduces the big values of reference/pari/src/test/32/quadclassunit', () => {
    // quadunit(74881)
    expect(quadunit(74881n)).toEqual([
      BigInt(
        '13313836214635923679350727498327540498727524435249418680727535284881984230516155252521913978327933871319513567734674995952944882291564414920804355520376637396859279722757945014417671528542599695402642789674901376031360881547644988472350241908693541649220052571440234192695534' +
          '5567'
      ),
      BigInt(
        '97664644524691035450385949078519979040360731478061598807512211227050101784793847267130640656779998175421277270552234496536985990846172799874910264740903918762169611136867754945467502719195263170711683406552475763065043739487933741699661116184342816008228619899979763439296' +
          '8896'
      ),
    ]);
  });

  it('always produces a unit: N(u + v w_D) = +/- 1', () => {
    // reference/pari/src/test/in/quadclassunit:11-16 runs exactly this loop
    for (let i = 1n; i <= 2000n; i++) {
      const D = 2n * i - (i % 2n);
      const r = isqrt(D);
      if (r * r === D) continue;
      const nrm = quadunitnorm(D);
      expect(nrm === 1n || nrm === -1n).toBe(true);
    }
  });

  it('rejects the arguments PARI rejects', () => {
    expect(() => quadunit(-3n)).toThrow('disc <= 0');
    expect(() => quadunit(7n)).toThrow('disc % 4 > 1');
    expect(() => quadunit(9n)).toThrow('issquare(disc) = 1');
  });
});

describe('primedec: Buchmann-Lenstra (PARI base2.c:2248)', () => {
  /** Multiplication table of O_K in the integral basis. */
  const mulTable = (K: NumberField): bigint[][][] => {
    const n = K.degree();
    const basis = K._pari_integral_basis();
    const Winv = ratInverse(basis.map((b) => b.list()));
    const coords = (x: NumberFieldElement): bigint[] => {
      const l = x.list();
      const row: bigint[] = [];
      for (let k = 0; k < n; k++) {
        let acc = Rational.zero();
        for (let t = 0; t < n; t++) acc = acc.add(l[t]!.mul(Winv[t]![k]!));
        expect(acc.denominator).toBe(1n);
        row.push(acc.numerator);
      }
      return row;
    };
    const T: bigint[][][] = [];
    for (let i = 0; i < n; i++) {
      const r: bigint[][] = [];
      for (let j = 0; j < n; j++) r.push(coords(basis[i]!.mul(basis[j]!)));
      T.push(r);
    }
    return T;
  };

  const latticeMul = (A: bigint[][], B: bigint[][], mt: bigint[][][]): bigint[][] => {
    const n = A.length;
    const rows: bigint[][] = [];
    for (const a of A) {
      for (const b of B) {
        const out = new Array<bigint>(n).fill(0n);
        for (let i = 0; i < n; i++) {
          for (let j = 0; j < n; j++) {
            const c = a[i]! * b[j]!;
            if (c === 0n) continue;
            const r = mt[i]![j]!;
            for (let k = 0; k < n; k++) out[k] = out[k]! + c * r[k]!;
          }
        }
        rows.push(out);
      }
    }
    return hnf(rows, n);
  };

  const primeLattice = (gens: bigint[][], p: bigint, n: number): bigint[][] => {
    const rows = gens.map((g) => [...g]);
    for (let i = 0; i < n; i++) {
      const r = new Array<bigint>(n).fill(0n);
      r[i] = p;
      rows.push(r);
    }
    return hnf(rows, n);
  };

  const key = (M: bigint[][]): string => M.map((r) => r.join(',')).join(';');

  const CUBICS: ZPoly[] = [
    [-8n, -2n, -1n, 1n], // Dedekind's x^3 - x^2 - 2x - 8: 2 is inessential
    [8n, -2n, 1n, 1n],
    [-2n, 0n, 0n, 1n],
    [-10n, -9n, -6n, 1n], // 2 inessential
    [-2n, -1n, -6n, 1n], // 2 inessential
    [3n, 3n, 0n, 1n],
    [-1n, -1n, 0n, 1n],
  ];

  it('sum e_i f_i = n, and prod P_i^e_i = p O_K exactly', () => {
    let roundFour = 0;
    for (const poly of CUBICS) {
      const K = NumberFieldConstructor(poly, 'a');
      const n = K.degree();
      const mt = mulTable(K);
      for (const p of [2n, 3n, 5n, 7n, 11n, 13n]) {
        const dec = primedec(mt, p);
        let s = 0n;
        for (const d of dec) s += d.e * d.f;
        expect(s).toBe(BigInt(n));
        if (nfbasis(K.pari_polynomial()).index % p === 0n) roundFour++;
        let prod: bigint[][] | null = null;
        for (const d of dec) {
          const P = primeLattice(d.gens, p, n);
          let acc = P;
          for (let t = 1n; t < d.e; t++) acc = latticeMul(acc, P, mt);
          prod = prod === null ? acc : latticeMul(prod, acc, mt);
        }
        const pO = hnf(
          Array.from({ length: n }, (_, i) =>
            Array.from({ length: n }, (_, j) => (i === j ? p : 0n))
          ),
          n
        );
        expect(key(prod!)).toBe(key(pO));
      }
    }
    // the round-4 branch really was exercised
    expect(roundFour).toBeGreaterThan(0);
  });

  it('agrees with Dedekind-Kummer wherever that theorem applies', () => {
    for (const poly of CUBICS) {
      const K = NumberFieldConstructor(poly, 'a');
      const mt = mulTable(K);
      for (const p of [2n, 3n, 5n, 7n, 11n, 13n]) {
        const dec = primedec(mt, p);
        const fromField = K.decomposition(p);
        const norm = (a: Array<[bigint, bigint]>): string =>
          a
            .map((x) => x.join(':'))
            .sort()
            .join(',');
        expect(norm(dec.map((d) => [p ** d.f, d.e] as [bigint, bigint]))).toBe(
          norm(fromField.map(([P, e]) => [P.norm().numerator, e] as [bigint, bigint]))
        );
      }
    }
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
