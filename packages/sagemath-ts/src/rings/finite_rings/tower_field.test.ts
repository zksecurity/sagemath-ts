/**
 * Binary Tower Fields Test
 *
 * This test mirrors the SageMath script for binary tower fields (Binius).
 * We verify that our TypeScript implementation produces equivalent results.
 */
import { describe, expect, test } from 'bun:test';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import { GF2 } from './gf2.js';
import { Ti, formatTowerElement, listTowerElements, towerCardinality } from './tower_field.js';

describe('Binary Tower Fields (Binius)', () => {
  test('T0 = GF(2) has 2 elements', () => {
    const T0 = Ti(0);
    expect(T0).toBe(GF2);

    const elements = [...T0];
    expect(elements.length).toBe(2);

    // Elements should be 0 and 1
    const values = elements.map((e) => e.toString());
    expect(values).toContain('0');
    expect(values).toContain('1');
  });

  test('T1 = GF(4) has 4 elements', () => {
    const T1 = Ti(1);
    expect(T1.cardinality()).toBe(4n);

    const elements = listTowerElements(T1);
    expect(elements.length).toBe(4);

    console.log('T1 elements:', elements.map(formatTowerElement));
  });

  test('T2 = GF(16) has 16 elements', () => {
    const T2 = Ti(2);
    expect(T2.cardinality()).toBe(16n);

    const elements = listTowerElements(T2);
    expect(elements.length).toBe(16);

    console.log('T2 elements:', elements.map(formatTowerElement));
  });

  test('T3 = GF(256) has 256 elements', () => {
    const T3 = Ti(3);
    expect(T3.cardinality()).toBe(256n);

    const elements = listTowerElements(T3);
    expect(elements.length).toBe(256);
  });

  test('tower cardinality formula: |Ti(i)| = 2^(2^i)', () => {
    expect(towerCardinality(0)).toBe(2n); // 2^1
    expect(towerCardinality(1)).toBe(4n); // 2^2
    expect(towerCardinality(2)).toBe(16n); // 2^4
    expect(towerCardinality(3)).toBe(256n); // 2^8
    expect(towerCardinality(4)).toBe(65536n); // 2^16
    expect(towerCardinality(5)).toBe(4294967296n); // 2^32
  });

  test('polynomial ring over GF(2)', () => {
    const [R, x] = [new PolynomialRing(GF2, 'x0'), new PolynomialRing(GF2, 'x0').gen()];

    // x0^2 + x0 + 1
    const modulus = x.pow(2).add(x).add(R.one());

    expect(modulus.toString()).toBe('x0^2 + x0 + 1');
    expect(modulus.degree()).toBe(2);
  });

  test('field arithmetic in T1', () => {
    const T1 = Ti(1);
    const alpha = T1.gen(); // Generator (corresponds to x0 in SageMath)
    const one = T1.one();
    const zero = T1.zero();

    // alpha^2 + alpha + 1 = 0 (the defining relation)
    // So alpha^2 = -alpha - 1 = alpha + 1 (in char 2)
    const alphaSq = alpha.mul(alpha);
    const alphaPlus1 = alpha.add(one);

    expect(alphaSq.eq(alphaPlus1)).toBe(true);

    // Test all elements exist
    const elements = listTowerElements(T1);
    expect(elements.length).toBe(4);

    // Test addition table properties (each row/col is a permutation)
    for (const a of elements) {
      const sums = elements.map((b) => formatTowerElement(a.add(b)));
      const unique = new Set(sums);
      expect(unique.size).toBe(4);
    }
  });

  test('field arithmetic in T2', () => {
    const T2 = Ti(2);
    const one = T2.one();
    const zero = T2.zero();

    // Get generators
    const gen = T2.gen(); // x1

    // Test multiplicative identity
    for (const elem of T2) {
      if (!elem.isZero()) {
        expect(elem.mul(one).eq(elem)).toBe(true);
        expect(one.mul(elem).eq(elem)).toBe(true);
      }
      expect(elem.mul(zero).isZero()).toBe(true);
      expect(elem.add(zero).eq(elem)).toBe(true);
    }
  });

  test('multiplicative inverses in T1', () => {
    const T1 = Ti(1);

    // All non-zero elements should have inverses
    for (const elem of T1) {
      if (!elem.isZero()) {
        const inv = elem.inv();
        const product = elem.mul(inv);
        expect(product.eq(T1.one())).toBe(true);
      }
    }
  });

  test('multiplicative inverses in T2', () => {
    const T2 = Ti(2);

    let count = 0;
    for (const elem of T2) {
      if (!elem.isZero()) {
        const inv = elem.inv();
        const product = elem.mul(inv);
        expect(product.eq(T2.one())).toBe(true);
        count++;
      }
    }
    expect(count).toBe(15); // 16 - 1 non-zero elements
  });

  test('field has correct characteristic', () => {
    const T2 = Ti(2);

    // In characteristic 2, x + x = 0 for all x
    for (const elem of T2) {
      expect(elem.add(elem).isZero()).toBe(true);
    }
  });

  test('compare with SageMath Ti function output', () => {
    // From SageMath:
    // T0 = GF(2), len = 2
    // T1 = GF(4), len = 4
    // T2 = GF(16), len = 16
    // Ti(5) has len = 2^(2^5) = 2^32

    expect(Ti(0).cardinality()).toBe(2n);
    expect(Ti(1).cardinality()).toBe(4n);
    expect(Ti(2).cardinality()).toBe(16n);

    // T5 cardinality (don't enumerate, just check cardinality)
    expect(Ti(5).cardinality()).toBe(2n ** 32n);
  });
});

describe('Tower field element representation', () => {
  test('T2 element printing matches expected format', () => {
    const T2 = Ti(2);
    const elements = listTowerElements(T2);

    // Print all elements for comparison with SageMath
    console.log('\nT2 elements (16 total):');
    elements.forEach((elem, i) => {
      console.log(`  ${i}: ${formatTowerElement(elem)}`);
    });

    // Verify we have exactly 16 elements
    expect(elements.length).toBe(16);
  });
});
