import { describe, expect, test } from 'bun:test';
import { inverse_mod } from '../../../sagemath-ts/src/arith/misc.js';
import type { EllipticCurvePoint } from '../../../sagemath-ts/src/schemes/elliptic_curves/index.js';
import {
  type EdwardsPoint,
  type MontgomeryPoint,
  type TwistedEdwardsCurve,
  createMontgomeryCurve,
  createTwistedEdwardsCurve,
} from './models.js';
import { namedCurve } from './named.js';
import { babyJubJub, curve25519, ed25519, jubjub } from './params.js';

const modp = (p: bigint) => (a: bigint) => ((a % p) + p) % p;

function edwardsAdd(c: TwistedEdwardsCurve, P: EdwardsPoint, Q: EdwardsPoint): EdwardsPoint {
  const p = c.field.characteristic;
  const m = modp(p);
  const t = m(c.d * P.x * Q.x * P.y * Q.y);
  return {
    x: m((P.x * Q.y + P.y * Q.x) * inverse_mod(m(1n + t), p)),
    y: m((P.y * Q.y - c.a * P.x * Q.x) * inverse_mod(m(1n - t), p)),
  };
}

/** A point of order exactly `cofactor` (all named curves here have cyclic torsion). */
function fullTorsionPoint(E: {
  curve: { random_point(): EllipticCurvePoint };
  order: bigint;
  cofactor: bigint;
}) {
  for (;;) {
    const T = E.curve.random_point().mul(E.order);
    if (!T.mul(E.cofactor / 2n).isZero()) return T;
  }
}

describe('named Edwards curves', () => {
  const cases = [
    ['ed25519', ed25519],
    ['babyJubJub', babyJubJub],
    ['jubjub', jubjub],
  ] as const;

  for (const [name, params] of cases) {
    test(`${name}: generator and group law`, () => {
      const c = namedCurve(name);
      const G = params.generator;
      expect(c.fromWeierstrass(c.generator)).toEqual(G);
      expect(c.toWeierstrass(G).eq(c.generator)).toBe(true);
      const P = c.fromWeierstrass(c.generator.mul(0x1234567n));
      const Q = c.fromWeierstrass(c.generator.mul(0xabcdef0123n));
      expect(c.isOnCurve(P) && c.isOnCurve(Q)).toBe(true);
      expect(
        c.toWeierstrass(edwardsAdd(c, P, Q)).eq(c.toWeierstrass(P).add(c.toWeierstrass(Q)))
      ).toBe(true);
      expect(c.toWeierstrass(edwardsAdd(c, P, P)).eq(c.toWeierstrass(P).double())).toBe(true);
    });

    test(`${name}: the torsion points round-trip`, () => {
      const c = namedCurve(name);
      const p = c.field.characteristic;
      const T = fullTorsionPoint(c);
      let W = T;
      const seen: EdwardsPoint[] = [];
      for (let k = 1n; k <= c.cofactor; k++) {
        const E = c.fromWeierstrass(W);
        expect(c.isOnCurve(E)).toBe(true);
        expect(c.toWeierstrass(E).eq(W)).toBe(true);
        seen.push(E);
        W = W.add(T);
      }
      expect(seen).toContainEqual({ x: 0n, y: 1n });
      expect(seen).toContainEqual({ x: 0n, y: p - 1n });
      expect(c.toWeierstrass({ x: 0n, y: 1n }).isZero()).toBe(true);
      const two = c.toWeierstrass({ x: 0n, y: p - 1n });
      expect(two.isZero()).toBe(false);
      expect(two.double().isZero()).toBe(true);
      expect(c.toMontgomery({ x: 0n, y: p - 1n })).toEqual({ u: 0n, v: 0n });
      expect(c.fromMontgomery(null)).toEqual({ x: 0n, y: 1n });
    });
  }

  test('ed25519 shares its Weierstrass curve with curve25519 (RFC 7748 §4.1)', () => {
    const ed = namedCurve('ed25519');
    const x = namedCurve('curve25519');
    const p = ed.field.characteristic;
    expect(ed.montgomery).toBe(x);
    expect(ed.curve).toBe(x.curve);
    expect(ed.toMontgomery(ed25519.generator)).toEqual(curve25519.generator);
    // v = sqrt(-486664)·u/x
    const { u, v } = curve25519.generator;
    const s = modp(p)(v * ed25519.generator.x * inverse_mod(u, p));
    expect(modp(p)(s * s + 486664n)).toBe(0n);
  });

  test('babyJubJub maps to the EIP-2494 Montgomery form', () => {
    const c = namedCurve('babyJubJub');
    expect([c.montgomery.A, c.montgomery.B]).toEqual([168698n, 1n]);
  });

  test('invalid input is rejected', () => {
    const c = namedCurve('ed25519');
    expect(() => c.toWeierstrass({ x: 1n, y: 1n })).toThrow('is not on ed25519');
    expect(() => c.fromMontgomery({ u: 9n, v: 1n })).toThrow('is not on curve25519');
    expect(() => c.fromWeierstrass(namedCurve('secp256k1').generator)).toThrow(
      "not a point of curve25519's curve"
    );
    expect(() => createTwistedEdwardsCurve(babyJubJub, namedCurve('curve25519'))).toThrow(
      'does not correspond'
    );
  });
});

describe('X25519 u-coordinates', () => {
  const x = namedCurve('curve25519');
  const p = x.field.characteristic;

  test('the low-order u values', () => {
    const lowOrder = [
      0n,
      1n,
      325606250916557431795983626356110631294008115727848805560023387167927233504n,
      39382357235489614581723060781553021112529911719440698176882885853963445705823n,
    ];
    for (const u of lowOrder) {
      const P = x.liftU(u)!;
      expect(P.mul(8n).isZero()).toBe(true);
      expect(x.fromWeierstrass(P)!.u).toBe(u);
    }
    expect(x.liftU(0n)!.double().isZero()).toBe(true);
    // -1 is on the twist since 486660 is not a square; p and p + 1 reduce to 0 and 1.
    expect(x.liftU(p - 1n)).toBeNull();
    expect(x.liftU(p)!.eq(x.liftU(0n)!)).toBe(true);
    expect(x.liftU(p + 1n)!.eq(x.liftU(1n)!)).toBe(true);
  });

  test('liftU recovers ±P and the identity maps to null', () => {
    const P = x.generator.mul(0x5eedn);
    const L = x.liftU(x.fromWeierstrass(P)!.u)!;
    expect(L.eq(P) || L.eq(P.neg())).toBe(true);
    expect(x.fromWeierstrass(x.curve.zero())).toBeNull();
    expect(x.toWeierstrass(null).isZero()).toBe(true);
    expect(x.liftU(9n)!.eq(x.generator) || x.liftU(9n)!.eq(x.generator.neg())).toBe(true);
  });

  test('Montgomery and Weierstrass group laws agree', () => {
    const m = modp(p);
    const P = x.fromWeierstrass(x.generator.mul(77n))!;
    const Q = x.fromWeierstrass(x.generator.mul(1000003n))!;
    const l = m((Q.v - P.v) * inverse_mod(m(Q.u - P.u), p));
    const u = m(x.B * l * l - x.A - P.u - Q.u);
    const sum: MontgomeryPoint = { u, v: m(l * (P.u - u) - P.v) };
    expect(x.toWeierstrass(sum).eq(x.generator.mul(1000080n))).toBe(true);
  });
});

describe('toy curves, exhaustively', () => {
  const p = 103n;
  const m = modp(p);
  const isSquare = (a: bigint) => m(a) === 0n || modPowSmall(m(a), (p - 1n) / 2n) === 1n;
  function modPowSmall(b: bigint, e: bigint): bigint {
    let r = 1n;
    for (let i = 0n; i < e; i++) r = m(r * b);
    return r;
  }

  function weierstrassPoints(E: TwistedEdwardsCurve): EllipticCurvePoint[] {
    const points = [E.curve.zero()];
    for (let x = 0n; x < p; x++) {
      for (const P of E.curve.lift_x(x, true)) points.push(P);
    }
    return points;
  }

  function edwardsPoints(a: bigint, d: bigint): EdwardsPoint[] {
    const points: EdwardsPoint[] = [];
    for (let x = 0n; x < p; x++) {
      for (let y = 0n; y < p; y++) {
        if (m(a * x * x + y * y - 1n - d * x * x * y * y) === 0n) points.push({ x, y });
      }
    }
    return points;
  }

  // [label, a square, d square, d/a square]; the last two decide which points lie at infinity.
  const shapes = [
    ['complete', true, false],
    ['d and d/a square', true, true],
    ['d square, d/a not', false, true],
    ['d/a square, d not', false, false],
  ] as const;

  for (const [label, aSquare, dSquare] of shapes) {
    test(label, () => {
      let a = 2n;
      while (isSquare(a) !== aSquare) a++;
      let d = a + 1n;
      while (isSquare(d) !== dSquare) d++;
      const edwards = edwardsPoints(a, d);
      const G = edwards.find((P) => P.x !== 0n)!;
      // The maps do not depend on the orders.
      const c = createTwistedEdwardsCurve({
        name: 'toy',
        field: { p },
        a,
        d,
        order: 1n,
        cofactor: 1n,
        subgroupOrder: 1n,
        generator: G,
      });

      const images = new Set<string>();
      for (const P of edwards) {
        const W = c.toWeierstrass(P);
        images.add(String(W));
        expect(c.fromWeierstrass(W)).toEqual(P);
      }
      expect(images.size).toBe(edwards.length);

      let atInfinity = 0;
      const all = weierstrassPoints(c);
      for (const W of all) {
        let E: EdwardsPoint;
        try {
          E = c.fromWeierstrass(W);
        } catch (e) {
          expect((e as Error).message).toContain('lies at infinity on toy');
          atInfinity++;
          continue;
        }
        expect(c.toWeierstrass(E).eq(W)).toBe(true);
      }
      const expected = (dSquare ? 2 : 0) + (isSquare(m(d * inverse_mod(a, p))) ? 2 : 0);
      expect(atInfinity).toBe(expected);
      expect(all.length).toBe(edwards.length + expected);

      for (const P of edwards.slice(0, 12)) {
        for (const Q of edwards.slice(0, 12)) {
          const t = m(d * P.x * Q.x * P.y * Q.y);
          if (t === 1n || t === p - 1n) continue;
          expect(
            c.toWeierstrass(edwardsAdd(c, P, Q)).eq(c.toWeierstrass(P).add(c.toWeierstrass(Q)))
          ).toBe(true);
        }
      }
    });
  }

  test('a Montgomery curve', () => {
    const A = 5n;
    const B = 7n;
    const rhs = (u: bigint) => m(u * (u * u + A * u + 1n) * inverse_mod(B, p));
    const affine: MontgomeryPoint[] = [];
    for (let u = 0n; u < p; u++) {
      for (let v = 0n; v < p; v++) if (m(v * v) === rhs(u)) affine.push({ u, v });
    }
    const c = createMontgomeryCurve({
      name: 'toy',
      field: { p },
      A,
      B,
      order: 1n,
      cofactor: 1n,
      generator: affine.find((P) => P.v !== 0n)!,
    });
    for (const P of affine) expect(c.fromWeierstrass(c.toWeierstrass(P))).toEqual(P);
    for (let u = 0n; u < p; u++) expect(c.liftU(u) === null).toBe(!isSquare(rhs(u)));
    expect(c.curve.cardinality()).toBe(BigInt(affine.length + 1));
  });
});
