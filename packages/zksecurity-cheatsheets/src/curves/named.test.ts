import { describe, expect, test } from 'bun:test';
import { is_prime } from '../../../sagemath-ts/src/arith/misc.js';
import { createPrimeFieldCurve as createFromBarrel, curves } from '../curves.js';
import { namedCurve, namedCurveNames } from './named.js';
import { bls12_377_g1, bls12_381, bn254, grumpkin, p256 } from './params.js';
import { createPrimeFieldCurve, hasseCofactor } from './weierstrass.js';

describe('every named curve', () => {
  for (const name of namedCurveNames) {
    test(name, () => {
      const c = namedCurve(name);
      const p = c.field.characteristic;
      expect(c.generator.isZero()).toBe(false);
      expect(c.curve.is_on_curve(c.generator.x!, c.generator.y!)).toBe(true);
      expect(c.generator.mul(c.order).isZero()).toBe(true);
      // With `order` prime and above 4√p, the Hasse interval holds one multiple of it.
      expect(is_prime(c.order)).toBe(true);
      expect(hasseCofactor(p, c.order)).toBe(c.cofactor);
      const R = c.curve.random_point();
      expect(R.mul(c.cofactor * c.order).isZero()).toBe(true);
      expect(namedCurve(name)).toBe(c);
      // Seeded caches: neither call runs SEA.
      expect(c.generator.order()).toBe(c.order);
      expect(c.curve.cardinality()).toBe(c.cofactor * c.order);
    });
  }

  test('models', () => {
    expect(namedCurveNames.map((n) => [n, namedCurve(n).model])).toEqual([
      ['secp256k1', 'short-weierstrass'],
      ['p256', 'short-weierstrass'],
      ['p384', 'short-weierstrass'],
      ['p521', 'short-weierstrass'],
      ['bn254', 'short-weierstrass'],
      ['bls12_381', 'short-weierstrass'],
      ['bls12_377', 'short-weierstrass'],
      ['grumpkin', 'short-weierstrass'],
      ['pallas', 'short-weierstrass'],
      ['vesta', 'short-weierstrass'],
      ['curve25519', 'montgomery'],
      ['ed25519', 'twisted-edwards'],
      ['babyJubJub', 'twisted-edwards'],
      ['jubjub', 'twisted-edwards'],
    ]);
  });

  test('unknown names are rejected', () => {
    expect(() => namedCurve('p224' as 'p256')).toThrow('unknown curve p224');
  });
});

describe('parameters', () => {
  test('BLS cofactors follow from the curve seed', () => {
    const x381 = -0xd201000000010000n;
    expect(bls12_381.field.p + 1n - (x381 + 1n)).toBe(bls12_381.g1.cofactor! * bls12_381.r);
    const x377 = 0x8508c00000000001n;
    expect(bls12_377_g1.order).toBe(x377 ** 4n - x377 ** 2n + 1n);
    expect(bls12_377_g1.field.p + 1n - (x377 + 1n)).toBe(
      bls12_377_g1.cofactor! * bls12_377_g1.order
    );
  });

  test('the cycles swap base and scalar fields', () => {
    expect(grumpkin.field.p).toBe(bn254.r);
    expect(grumpkin.order).toBe(bn254.field.p);
    expect(curves.pallas.order).toBe(curves.vesta.field.p);
    expect(curves.vesta.order).toBe(curves.pallas.field.p);
  });

  test('Jubjub d = -10240/10241 over the BLS12-381 scalar field', () => {
    const { p } = curves.jubjub.field;
    expect((curves.jubjub.d * 10241n + 10240n) % p).toBe(0n);
    expect(p).toBe(bls12_381.r);
  });
});

describe('createPrimeFieldCurve', () => {
  test('keeps its field, curve and generator', () => {
    const c = createFromBarrel(p256);
    expect(c.field.characteristic).toBe(p256.field.p);
    expect(c.generator.mul(p256.order).isZero()).toBe(true);
    expect(c.order).toBe(p256.order);
    expect(c.cofactor).toBe(1n);
  });

  test('derives a missing cofactor only when the order fixes it', () => {
    const { cofactor: _, ...noCofactor } = curves.bls12_381_g1;
    expect(createPrimeFieldCurve(noCofactor).cofactor).toBe(bls12_381.g1.cofactor!);
    // y² = x³ + 2x + 3 over GF(97) has 100 points; (3, 6) has order 5 < 4√97.
    const toy = {
      name: 'toy',
      field: { p: 97n },
      a: 2n,
      b: 3n,
      order: 5n,
      generator: { x: 3n, y: 6n },
    };
    expect(() => createPrimeFieldCurve(toy)).toThrow('does not determine the cofactor');
    expect(createPrimeFieldCurve({ ...toy, cofactor: 20n }).cofactor).toBe(20n);
  });

  test('seeds the order caches only with consistent values', () => {
    const toy = { name: 'toy', field: { p: 97n }, a: 2n, b: 3n, generator: { x: 3n, y: 6n } };
    const right = createPrimeFieldCurve({ ...toy, order: 5n, cofactor: 20n });
    expect(right.generator.order()).toBe(5n);
    const wrong = createPrimeFieldCurve({ ...toy, order: 7n, cofactor: 20n });
    expect(wrong.generator.order()).toBe(5n);
    expect(wrong.curve.cardinality()).toBe(100n);
  });

  test('hasseCofactor', () => {
    expect(hasseCofactor(97n, 5n)).toBeUndefined();
    expect(hasseCofactor(97n, 50n)).toBe(2n);
    // 60 > 4√97, but no multiple of it lies in [79, 117].
    expect(hasseCofactor(97n, 60n)).toBeUndefined();
  });
});
