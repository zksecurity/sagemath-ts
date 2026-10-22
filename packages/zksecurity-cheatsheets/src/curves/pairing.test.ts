import { describe, expect, test } from 'bun:test';
import { isqrt } from '../../../sagemath-ts/src/arith/misc.js';
import {
  ate_pairing,
  tate_pairing,
  weil_pairing,
} from '../../../sagemath-ts/src/schemes/elliptic_curves/index.js';
import { namedCurve } from './named.js';
import { type PairingCurveName, pairingCurve } from './pairing.js';
import { bls12_381, bn254 } from './params.js';

const seeds: Record<PairingCurveName, { trace: bigint }> = {
  // t = 6x² + 1 and t = x + 1 for the curve seeds x.
  bn254: { trace: 6n * 4965661367192848881n ** 2n + 1n },
  bls12_381: { trace: -0xd201000000010000n + 1n },
};

for (const name of ['bn254', 'bls12_381'] as const) {
  describe(name, () => {
    const pc = pairingCurve(name);
    const { g1, g2 } = pc;
    const p = g1.field.characteristic;
    const r = g1.order;

    test('G2 generator has order r', () => {
      expect(pc.g1).toBe(namedCurve(name));
      expect(g2.order).toBe(r);
      expect(g2.generator.is_zero()).toBe(false);
      expect(g2.generator.mul(r).is_zero()).toBe(true);
      expect(g2.generator.order()).toBe(r);
      expect(g2.curve.cardinality()).toBe(g2.cofactor * r);
    });

    test('trace and cofactors agree with the curve seed', () => {
      expect(pc.trace).toBe(seeds[name].trace);
      // The sextic twists of E over Fp2 have order p² + 1 - t' with t' in this list.
      const t2 = pc.trace ** 2n - 2n * p;
      const f = isqrt((4n * p * p - t2 * t2) / 3n);
      expect(3n * f * f).toBe(4n * p * p - t2 * t2);
      const traces = [t2, (t2 + 3n * f) / 2n, (t2 - 3n * f) / 2n].flatMap((t) => [t, -t]);
      expect(traces.map((t) => p * p + 1n - t)).toContain(g2.cofactor * r);
      let x = 1n;
      while (!g2.curve.is_x_coord(pc.fp2(x, 1n))) x++;
      const R = g2.curve.lift_x(pc.fp2(x, 1n));
      expect(R.mul(r).is_zero()).toBe(false);
      expect(R.mul(g2.cofactor * r).is_zero()).toBe(true);
      expect(R.mul(g2.cofactor).mul(r).is_zero()).toBe(true);
    });

    test('Fp2 embeds in Fp12', () => {
      const u = pc.Fp2.gen();
      const a = pc.fp2(123n, 456n);
      const b = pc.fp2(p - 7n, 11n);
      const emb = pc.fp2ToFp12;
      expect(emb(u).pow(2n).eq(pc.Fp12.__call__(-1n))).toBe(true);
      expect(emb(a.mul(b)).eq(emb(a).mul(emb(b)))).toBe(true);
      expect(emb(a.add(b)).eq(emb(a).add(emb(b)))).toBe(true);
      expect(emb(pc.xi).eq(pc.Fp12.gen().pow(6n))).toBe(true);
    });

    test('untwist and liftG1 are homomorphisms into E(Fp12)', () => {
      const Q1 = g2.generator.mul(3n);
      const Q2 = g2.generator.mul(5n);
      expect(pc.untwist(Q1.add(Q2)).eq(pc.untwist(Q1).add(pc.untwist(Q2)))).toBe(true);
      expect(pc.untwist(g2.curve.zero()).is_zero()).toBe(true);
      const P = g1.generator.mul(9n);
      expect(pc.liftG1(P.add(g1.generator)).eq(pc.liftG1(P).add(pc.liftG1(g1.generator)))).toBe(
        true
      );
      expect(pc.liftG1(g1.curve.zero()).is_zero()).toBe(true);
    });

    test('the standard coordinates parse only in the u² + 1 tower order', () => {
      const params = name === 'bn254' ? bn254 : bls12_381;
      const [x0, x1] = params.g2.x;
      const [y0, y1] = params.g2.y;
      expect(() => pc.g2Point([x1, x0], [y1, y0])).toThrow();
    });

    test('non-degenerate and bilinear on planted scalars', () => {
      const P = g1.generator;
      const Q = g2.generator;
      const a = 0x1f2e3d4c5b6a7988n;
      const b = 0x0123456789abcdefn;
      const aP = pc.liftG1(P.mul(a));
      const bQ = pc.untwist(Q.mul(b));
      const P12 = pc.liftG1(P);
      const Q12 = pc.untwist(Q);

      const e = pc.pairing(P, Q);
      expect(e.eq(ate_pairing(P12, Q12, r, 12, pc.trace))).toBe(true);
      expect(e.isOne()).toBe(false);
      expect(e.pow(r).isOne()).toBe(true);
      expect(ate_pairing(aP, bQ, r, 12, pc.trace).eq(e.pow(a * b))).toBe(true);

      const t = tate_pairing(P12, Q12, r, 12);
      expect(t.isOne()).toBe(false);
      expect(tate_pairing(aP, bQ, r, 12).eq(t.pow(a * b))).toBe(true);

      const w = weil_pairing(P12, Q12, r);
      expect(w.isOne()).toBe(false);
      expect(weil_pairing(aP, bQ, r).eq(w.pow(a * b))).toBe(true);
    });

    test('foreign points are rejected', () => {
      const other = pairingCurve(name === 'bn254' ? 'bls12_381' : 'bn254');
      expect(() => pc.untwist(other.g2.generator)).toThrow('build it with g2Point');
      expect(() => pc.liftG1(other.g1.generator)).toThrow('is not a point of');
    });
  });
}
