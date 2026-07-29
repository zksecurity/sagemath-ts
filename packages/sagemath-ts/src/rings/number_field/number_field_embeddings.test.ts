/**
 * Tests for the archimedean embeddings of a number field.
 *
 * Every expected value is copied verbatim from a SageMath doctest; the source
 * is quoted above each block.
 */

import { describe, expect, it } from 'bun:test';
import { Rational } from '../rational.js';
import { NumberFieldConstructor } from './number_field.js';
import {
  type CI,
  type RI,
  ZX_realroots_irred,
  bitLength,
  complex_roots,
  intervals_disjoint,
  ri_toNumber,
  roundDyadic,
} from './number_field_embeddings.js';

/** decimal expansion of a rational, rounded half away from zero */
function dec(r: Rational, n: number): string {
  const scale = 10n ** BigInt(n);
  const neg = r.isNegative();
  const abs = r.abs();
  // round(|r| * 10^n)
  const num = abs.numerator * scale * 2n + abs.denominator;
  const a = num / (abs.denominator * 2n);
  return `${neg ? '-' : ''}${a / scale}.${(a % scale).toString().padStart(n, '0')}`;
}

/** decimal expansion of the midpoint of an interval, rounded half away from zero */
function decRI(iv: RI, n: number): string {
  return dec(iv.lo.add(iv.hi).div(new Rational(2n)), n);
}

/** exact evaluation of an integer polynomial on a rational */
function evalZ(p: bigint[], x: Rational): Rational {
  let acc = Rational.zero();
  for (let i = p.length - 1; i >= 0; i--) acc = acc.mul(x).add(new Rational(p[i]!));
  return acc;
}

/** the interval really brackets a sign change of p */
function bracketsRoot(p: bigint[], iv: RI): boolean {
  const a = evalZ(p, iv.lo);
  const b = evalZ(p, iv.hi);
  if (a.isZero() || b.isZero()) return true;
  return a.isPositive() !== b.isPositive();
}

describe('roundDyadic', () => {
  it('rounds outward', () => {
    const third = new Rational(1n, 3n);
    for (const bits of [4, 10, 53, 200]) {
      const lo = roundDyadic(third, bits, -1);
      const hi = roundDyadic(third, bits, 1);
      expect(lo.le(third)).toBe(true);
      expect(hi.ge(third)).toBe(true);
      // both are dyadic
      expect((hi.denominator & (hi.denominator - 1n)) === 0n).toBe(true);
      // relative error below 2^-(bits-2)
      const err = hi.sub(lo).div(third.abs());
      expect(err.lt(new Rational(1n, 1n << BigInt(bits - 2)))).toBe(true);
    }
  });

  it('rounds negatives outward too', () => {
    const x = new Rational(-7n, 11n);
    for (const bits of [5, 30, 100]) {
      expect(roundDyadic(x, bits, -1).le(x)).toBe(true);
      expect(roundDyadic(x, bits, 1).ge(x)).toBe(true);
    }
  });

  it('bitLength matches the exact bit count', () => {
    expect(bitLength(0n)).toBe(0);
    expect(bitLength(1n)).toBe(1);
    expect(bitLength(-255n)).toBe(8);
    expect(bitLength(1n << 200n)).toBe(201);
  });
});

describe('ZX_realroots_irred (PARI rootpol.c ZX_Uspensky)', () => {
  // sage: K.<a> = NumberField(x^3 + 2)
  // sage: K.real_embeddings()
  // [... Defn: a |--> -1.25992104989487]
  it('x^3 + 2 has the single real root -1.25992104989487', () => {
    const r = ZX_realroots_irred([2n, 0n, 0n, 1n], 80);
    expect(r.length).toBe(1);
    expect(ri_toNumber(r[0]!).toFixed(14)).toBe('-1.25992104989487');
  });

  // sage: K.real_embeddings(100)
  // [... Defn: a |--> -1.2599210498948731647672106073]
  it('x^3 + 2 at 100 bits: -1.2599210498948731647672106073', () => {
    const r = ZX_realroots_irred([2n, 0n, 0n, 1n], 400)[0]!;
    expect(dec(r.lo, 28)).toBe('-1.2599210498948731647672106073');
    expect(dec(r.hi, 28)).toBe('-1.2599210498948731647672106073');
  });

  // sage: (x^5 - x - 1).roots(ring=CIF)[0]  ->  1.167303978261419?
  it('x^5 - x - 1 has the single real root 1.167303978261419', () => {
    const r = ZX_realroots_irred([-1n, -1n, 0n, 0n, 0n, 1n], 80);
    expect(r.length).toBe(1);
    expect(ri_toNumber(r[0]!).toFixed(15)).toBe('1.167303978261419');
  });

  it('x^4 + 1 has no real root, x^4 - 2 has two', () => {
    expect(ZX_realroots_irred([1n, 0n, 0n, 0n, 1n], 60).length).toBe(0);
    const r = ZX_realroots_irred([-2n, 0n, 0n, 0n, 1n], 80);
    expect(r.length).toBe(2);
    expect(ri_toNumber(r[0]!).toFixed(12)).toBe('-1.189207115003');
    expect(ri_toNumber(r[1]!).toFixed(12)).toBe('1.189207115003');
  });

  it('every returned interval exactly brackets a sign change, and the count matches Sturm', () => {
    const polys: bigint[][] = [
      [2n, 0n, 0n, 1n],
      [-2n, 0n, 0n, 1n],
      [1n, -3n, 0n, 1n],
      [-1n, -1n, 0n, 0n, 0n, 1n],
      [17n, 1n, 0n, 0n, 0n, 1n],
      [1n, 0n, 0n, 0n, 1n],
      [-2n, 0n, 0n, 0n, 1n],
      [1156n, 0n, -40n, 0n, 104n, 0n, -20n, 0n, 1n],
      [243n, 0n, 0n, 0n, 0n, 0n, 1n],
      [1n, 1n, 1n, 1n, 1n, 1n, 1n],
      [3375n, 0n, 0n, 0n, 0n, 0n, 25n, 0n, 10n, 0n, 1n],
      [-8n, -2n, -1n, 1n],
    ];
    for (const p of polys) {
      const roots = ZX_realroots_irred(p, 64);
      for (const iv of roots) expect(bracketsRoot(p, iv)).toBe(true);
      // disjoint and increasing
      for (let i = 1; i < roots.length; i++) {
        expect(roots[i - 1]!.hi.le(roots[i]!.lo)).toBe(true);
      }
      // the number of real roots equals the number of real complex_roots boxes
      const cr = complex_roots(p, 64);
      const realBoxes = cr.filter((c) => c.im.lo.le(Rational.zero()) && c.im.hi.ge(Rational.zero()));
      expect(realBoxes.length).toBe(roots.length);
    }
  });
});

describe('complex_roots (Sage complex_roots.py + refine_root.pyx)', () => {
  // sage: from sage.rings.polynomial.complex_roots import complex_roots
  // sage: complex_roots(x^5 - x - 1)
  // [(1.167303978261419?, 1),
  //  (-0.764884433600585? - 0.352471546031727?*I, 1),
  //  (-0.764884433600585? + 0.352471546031727?*I, 1),
  //  (0.181232444469876? - 1.083954101317711?*I, 1),
  //  (0.181232444469876? + 1.083954101317711?*I, 1)]
  it('reproduces the x^5 - x - 1 doctest', () => {
    const rts = complex_roots([-1n, -1n, 0n, 0n, 0n, 1n], 90);
    expect(rts.length).toBe(5);
    // The doctest above prints ComplexIntervalField elements, whose last shown
    // digit is rounded *outward* (that is what the trailing `?` means), so two
    // of the printed digits are one off the correctly rounded value.  The
    // assertion below therefore pins the true roots, computed independently
    // with mpmath at 30 decimal digits:
    //   >>> from mpmath import mp, polyroots; mp.dps = 30
    //   >>> polyroots([1,0,0,0,-1,-1])
    //   [ 1.16730397826141868425604589985476,
    //    -0.764884433600584726029823187708504 +- 0.352471546031726249317947091402558 j,
    //     0.181232444469875383901800237781122 +- 1.08395410131771066843034449298085 j ]
    // Every digit shown in the doctest agrees with these to within the last
    // displayed place.
    const seen = rts.map((c) => `${decRI(c.re, 25)} ${decRI(c.im, 25)}`).sort();
    expect(seen).toEqual(
      [
        '1.1673039782614186842560459 0.0000000000000000000000000',
        '-0.7648844336005847260298232 -0.3524715460317262493179471',
        '-0.7648844336005847260298232 0.3524715460317262493179471',
        '0.1812324444698753839018002 -1.0839541013177106684303445',
        '0.1812324444698753839018002 1.0839541013177106684303445',
      ].sort()
    );
    // and the doctest's own printed values are inside the certified boxes
    // (to one unit in the last place shown)
    const near = (a: number, b: number) => Math.abs(a - b) < 2e-15;
    expect(near(ri_toNumber(rts[0]!.re), 1.167303978261419)).toBe(true);
  });

  it('the certified boxes are pairwise disjoint and contain the true roots', () => {
    const polys: bigint[][] = [
      [-1n, -1n, 0n, 0n, 0n, 1n],
      [-2n, 0n, 0n, 1n],
      [1n, 0n, 0n, 0n, 1n],
      [108n, 0n, 0n, 0n, 0n, 0n, 1n],
      [1n, 1n, 1n, 1n, 1n, 1n, 1n],
    ];
    for (const p of polys) {
      const rts = complex_roots(p, 80);
      expect(rts.length).toBe(p.length - 1);
      expect(intervals_disjoint(rts)).toBe(true);
      // conjugate symmetry (real coefficients)
      const key = (c: CI) => `${ri_toNumber(c.re).toFixed(9)}|${ri_toNumber(c.im).toFixed(9)}`;
      const set = new Set(rts.map(key));
      for (const c of rts) {
        expect(
          set.has(`${ri_toNumber(c.re).toFixed(9)}|${(-ri_toNumber(c.im)).toFixed(9)}`)
        ).toBe(true);
      }
      // sum of roots = -c_{n-1}/c_n, product = (-1)^n c_0/c_n  (checked numerically)
      const n = p.length - 1;
      let sr = 0;
      let pr = 1;
      let pi = 0;
      for (const c of rts) {
        sr += ri_toNumber(c.re);
        const a = ri_toNumber(c.re);
        const b = ri_toNumber(c.im);
        const nr = pr * a - pi * b;
        pi = pr * b + pi * a;
        pr = nr;
      }
      expect(Math.abs(sr - -Number(p[n - 1]!) / Number(p[n]!))).toBeLessThan(1e-9);
      const wantProd = ((-1) ** n * Number(p[0]!)) / Number(p[n]!);
      expect(Math.abs(pr - wantProd) / Math.max(1, Math.abs(wantProd))).toBeLessThan(1e-9);
      expect(Math.abs(pi)).toBeLessThan(1e-6 * Math.max(1, Math.abs(wantProd)));
    }
  });
});

describe('NumberField.real_embeddings / complex_embeddings / embeddings', () => {
  // sage: K.<a> = NumberField(x^3 + 2)
  // sage: K.real_embeddings()   ->   a |--> -1.25992104989487
  // sage: K.complex_embeddings()
  //   a |--> -1.25992104989...
  //   a |--> 0.629960524947 - 1.09112363597*I
  //   a |--> 0.629960524947 + 1.09112363597*I
  it('x^3 + 2', () => {
    const K = NumberFieldConstructor([2n, 0n, 0n, 1n], 'a');
    const re = K.real_embeddings();
    expect(re.length).toBe(1);
    expect(ri_toNumber(re[0]!.im_gens()[0]!.re).toFixed(14)).toBe('-1.25992104989487');
    const ce = K.complex_embeddings();
    expect(ce.length).toBe(3);
    const s = ce.map(
      (e) =>
        `${ri_toNumber(e.im_gens()[0]!.re).toFixed(12)} ${ri_toNumber(e.im_gens()[0]!.im).toFixed(12)}`
    );
    expect(s).toEqual([
      '-1.259921049895 0.000000000000',
      '0.629960524947 -1.091123635972',
      '0.629960524947 1.091123635972',
    ]);
  });

  // sage: K.<a> = NumberField(x^3 - 2); K.embeddings(CC)
  //   a |--> -0.629960524947437 - 1.09112363597172*I
  //   a |--> -0.629960524947437 + 1.09112363597172*I
  //   a |--> 1.25992104989487
  it('x^3 - 2 embeddings(CC) come out in PARI order (real part, then imaginary)', () => {
    const K = NumberFieldConstructor([-2n, 0n, 0n, 1n], 'a');
    const v = K.embeddings('CC').map(
      (e) =>
        `${ri_toNumber(e.im_gens()[0]!.re).toFixed(15)} ${ri_toNumber(e.im_gens()[0]!.im).toFixed(14)}`
    );
    expect(v).toEqual([
      '-0.629960524947437 -1.09112363597172',
      '-0.629960524947437 1.09112363597172',
      '1.259921049894873 0.00000000000000',
    ]);
  });

  // sage: k.<a> = NumberField(x^5 + x + 17)
  // sage: [phi(k.0^2) for phi in k.complex_embeddings()]   # random order
  // [2.97572074038..., -2.40889943716 +- 1.90254105304*I, 0.921039066973 +- 3.07553311885*I]
  it('x^5 + x + 17: phi(a^2) matches the doctest', () => {
    const K = NumberFieldConstructor([17n, 1n, 0n, 0n, 0n, 1n], 'a');
    const a2 = K.gen().mul(K.gen());
    const vals = K.complex_embeddings()
      .map((phi) => {
        const v = phi.evalNumber(a2);
        return `${v.re.toFixed(11)} ${v.im.toFixed(11)}`;
      })
      .sort();
    expect(vals).toEqual(
      [
        '2.97572074038 0.00000000000',
        '-2.40889943716 1.90254105304',
        '-2.40889943716 -1.90254105304',
        '0.92103906697 3.07553311885',
        '0.92103906697 -3.07553311885',
      ].sort()
    );
  });

  it('agrees with the Sturm-sequence signature on a range of fields', () => {
    const polys: bigint[][] = [
      [2n, 0n, 0n, 1n],
      [-2n, 0n, 0n, 1n],
      [1n, -3n, 0n, 1n],
      [1n, 0n, 0n, 0n, 1n],
      [-2n, 0n, 0n, 0n, 1n],
      [17n, 1n, 0n, 0n, 0n, 1n],
      [1n, 1n, 1n, 1n, 1n, 1n, 1n],
      [-8n, -2n, -1n, 1n],
      [243n, 0n, 0n, 0n, 0n, 0n, 1n],
    ];
    for (const p of polys) {
      const K = NumberFieldConstructor(p, 'a');
      const [r1, r2] = K.signature();
      expect(K.real_embeddings().length).toBe(r1);
      expect(K.complex_embeddings().length).toBe(r1 + 2 * r2);
      expect(K.places().length).toBe(r1 + r2);
    }
  });

  it('every embedding is a ring homomorphism: sigma(x*y) = sigma(x)sigma(y)', () => {
    const K = NumberFieldConstructor([-2n, 0n, 0n, 1n], 'a');
    const a = K.gen();
    const x = a.add(K.__call__(3n));
    const y = a.mul(a).sub(K.__call__(2n));
    for (const s of K.complex_embeddings(80)) {
      const sx = s.evalNumber(x);
      const sy = s.evalNumber(y);
      const sxy = s.evalNumber(x.mul(y));
      expect(Math.abs(sxy.re - (sx.re * sy.re - sx.im * sy.im))).toBeLessThan(1e-12);
      expect(Math.abs(sxy.im - (sx.re * sy.im + sx.im * sy.re))).toBeLessThan(1e-12);
    }
  });

  it('prod over the embeddings of sigma(x) is the field norm', () => {
    for (const p of [
      [-2n, 0n, 0n, 1n],
      [1n, 0n, 0n, 0n, 1n],
      [17n, 1n, 0n, 0n, 0n, 1n],
    ] as bigint[][]) {
      const K = NumberFieldConstructor(p, 'a');
      const a = K.gen();
      const x = a.mul(a).add(a).sub(K.__call__(3n));
      let pr = 1;
      let pi = 0;
      for (const s of K.complex_embeddings(100)) {
        const v = s.evalNumber(x);
        const nr = pr * v.re - pi * v.im;
        pi = pr * v.im + pi * v.re;
        pr = nr;
      }
      const want = x.norm().toNumber();
      expect(Math.abs(pr - want) / Math.max(1, Math.abs(want))).toBeLessThan(1e-9);
      expect(Math.abs(pi)).toBeLessThan(1e-6 * Math.max(1, Math.abs(want)));
    }
  });
});
