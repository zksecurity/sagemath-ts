/**
 * Tests for complex multiplication of elliptic curves.
 *
 * Values are taken from the doctests of
 * sage/schemes/elliptic_curves/cm.py.
 */

import { describe, expect, it } from 'bun:test';
import { NotImplementedError, ValueError } from '../../errors.js';
import {
  cm_j_invariants,
  cm_j_invariants_and_orders,
  cm_orders,
  discriminants_with_bounded_class_number,
  is_cm_j_invariant,
  largest_disc_with_class_number,
  largest_fundamental_disc_with_class_number,
} from './cm.js';

const pairs = (l: Array<[bigint, bigint]>): string =>
  l.map(([d, f]) => `(${d},${f})`).join(', ');

describe('cm_j_invariants', () => {
  // cm.py:cm_j_invariants doctest for QQ (sorted ascending)
  it('matches SageMath over QQ', () => {
    expect(cm_j_invariants('QQ')).toEqual([
      -262537412640768000n,
      -147197952000n,
      -884736000n,
      -12288000n,
      -884736n,
      -32768n,
      -3375n,
      0n,
      1728n,
      8000n,
      54000n,
      287496n,
      16581375n,
    ]);
  });

  it('throws for general number fields', () => {
    expect(() => cm_j_invariants({ toString: () => 'Number Field' })).toThrow(NotImplementedError);
  });
});

describe('cm_j_invariants_and_orders', () => {
  // cm.py:cm_j_invariants_and_orders doctest for QQ: 3-tuples (D, f, j)
  it('returns (D, f, j) triples in SageMath order', () => {
    expect(cm_j_invariants_and_orders('QQ')).toEqual([
      [-3n, 3n, -12288000n],
      [-3n, 2n, 54000n],
      [-3n, 1n, 0n],
      [-4n, 2n, 287496n],
      [-4n, 1n, 1728n],
      [-7n, 2n, 16581375n],
      [-7n, 1n, -3375n],
      [-8n, 1n, 8000n],
      [-11n, 1n, -32768n],
      [-19n, 1n, -884736n],
      [-43n, 1n, -884736000n],
      [-67n, 1n, -147197952000n],
      [-163n, 1n, -262537412640768000n],
    ]);
  });
});

describe('discriminants_with_bounded_class_number', () => {
  // cm.py doctest: v = discriminants_with_bounded_class_number(3)
  it('matches SageMath for hmax = 3', () => {
    const v = discriminants_with_bounded_class_number(3);
    expect([...v.keys()].sort((a, b) => Number(a - b))).toEqual([1n, 2n, 3n]);
    expect(pairs(v.get(1n)!)).toBe(
      '(-3,1), (-3,2), (-3,3), (-4,1), (-4,2), (-7,1), (-7,2), (-8,1), (-11,1), (-19,1), (-43,1), (-67,1), (-163,1)'
    );
    expect(pairs(v.get(2n)!)).toBe(
      '(-3,4), (-3,5), (-3,7), (-4,3), (-4,4), (-4,5), (-7,4), (-8,2), (-8,3), (-11,3), (-15,1), (-15,2), (-20,1), (-24,1), (-35,1), (-40,1), (-51,1), (-52,1), (-88,1), (-91,1), (-115,1), (-123,1), (-148,1), (-187,1), (-232,1), (-235,1), (-267,1), (-403,1), (-427,1)'
    );
    expect(pairs(v.get(3n)!)).toBe(
      '(-3,6), (-3,9), (-11,2), (-19,2), (-23,1), (-23,2), (-31,1), (-31,2), (-43,2), (-59,1), (-67,2), (-83,1), (-107,1), (-139,1), (-163,2), (-211,1), (-283,1), (-307,1), (-331,1), (-379,1), (-499,1), (-547,1), (-643,1), (-883,1), (-907,1)'
    );
  });

  // cm.py doctest: sorted(len(v[h]) for h in v) with hmax = 8
  it('matches SageMath cardinalities for hmax = 8', () => {
    const v = discriminants_with_bounded_class_number(8);
    const lens = [...v.values()].map((l) => l.length).sort((a, b) => a - b);
    expect(lens).toEqual([13, 25, 29, 29, 38, 84, 101, 208]);
  });

  // H101: this used to take 6.26 s because of an O(|D|) character sum per
  // discriminant with no caching.
  it('computes hmax = 8 quickly', () => {
    const t0 = Date.now();
    discriminants_with_bounded_class_number(8);
    expect(Date.now() - t0).toBeLessThan(2000);
  });

  // cm.py doctest: discriminants_with_bounded_class_number(5, 50)
  it('matches SageMath for hmax = 5 with a discriminant bound', () => {
    const v = discriminants_with_bounded_class_number(5, 50);
    expect(pairs(v.get(1n)!)).toBe(
      '(-3,1), (-3,2), (-3,3), (-4,1), (-4,2), (-7,1), (-7,2), (-8,1), (-11,1), (-19,1), (-43,1)'
    );
    expect(pairs(v.get(2n)!)).toBe(
      '(-3,4), (-4,3), (-8,2), (-15,1), (-20,1), (-24,1), (-35,1), (-40,1)'
    );
    expect(pairs(v.get(3n)!)).toBe('(-11,2), (-23,1), (-31,1)');
    expect(pairs(v.get(4n)!)).toBe('(-39,1)');
    expect(pairs(v.get(5n)!)).toBe('(-47,1)');
  });

  it('requires a bound when hmax > 100', () => {
    expect(() => discriminants_with_bounded_class_number(101)).toThrow(ValueError);
  });
});

describe('cm_orders', () => {
  it('returns the empty list for h = 0', () => {
    expect(cm_orders(0)).toEqual([]);
  });

  // cm.py:cm_orders doctest for h = 1
  it('matches SageMath for h = 1', () => {
    expect(pairs(cm_orders(1))).toBe(
      '(-3,1), (-3,2), (-3,3), (-4,1), (-4,2), (-7,1), (-7,2), (-8,1), (-11,1), (-19,1), (-43,1), (-67,1), (-163,1)'
    );
  });

  // cm.py:cm_orders doctest: len(cm_orders(2)) == 29 and all have class number 2
  it('matches SageMath for h = 2', () => {
    const v = cm_orders(2);
    expect(v.length).toBe(29);
    expect(pairs(v)).toBe(
      '(-3,4), (-3,5), (-3,7), (-4,3), (-4,4), (-4,5), (-7,4), (-8,2), (-8,3), (-11,3), (-15,1), (-15,2), (-20,1), (-24,1), (-35,1), (-40,1), (-51,1), (-52,1), (-88,1), (-91,1), (-115,1), (-123,1), (-148,1), (-187,1), (-232,1), (-235,1), (-267,1), (-403,1), (-427,1)'
    );
  });

  // H101: hDf_dict is now populated, so a repeat query is served from cache.
  it('serves repeat queries from the cache', () => {
    cm_orders(3);
    const t0 = Date.now();
    const again = cm_orders(3);
    expect(Date.now() - t0).toBeLessThan(50);
    expect(again.length).toBe(25);
  });
});

describe('largest discriminants with given class number', () => {
  // cm.py doctests: both return the pair (|D|, number of such discriminants)
  it('matches SageMath', () => {
    expect(largest_fundamental_disc_with_class_number(1)).toEqual([163n, 9n]);
    expect(largest_fundamental_disc_with_class_number(2)).toEqual([427n, 18n]);
    expect(largest_fundamental_disc_with_class_number(10)).toEqual([13843n, 87n]);
    expect(largest_fundamental_disc_with_class_number(100)).toEqual([1856563n, 1736n]);
    expect(largest_fundamental_disc_with_class_number(6)).toEqual([3763n, 51n]);
    expect(largest_disc_with_class_number(1)).toEqual([163n, 13n]);
    expect(largest_disc_with_class_number(2)).toEqual([427n, 29n]);
    expect(largest_disc_with_class_number(10)).toEqual([13843n, 123n]);
    expect(largest_disc_with_class_number(100)).toEqual([1856563n, 2311n]);
    expect(largest_disc_with_class_number(6)).toEqual([4075n, 101n]);
  });

  it('returns (0, 0) for non-positive h', () => {
    expect(largest_fundamental_disc_with_class_number(0)).toEqual([0n, 0n]);
    expect(largest_disc_with_class_number(0)).toEqual([0n, 0n]);
  });

  it('throws for h > 100', () => {
    expect(() => largest_fundamental_disc_with_class_number(101)).toThrow(NotImplementedError);
    expect(() => largest_disc_with_class_number(101)).toThrow(NotImplementedError);
  });
});

describe('is_cm_j_invariant', () => {
  // cm.py:is_cm_j_invariant doctests return (True, (D, f)) / (False, None)
  it('recognises the rational CM j-invariants', () => {
    expect(is_cm_j_invariant(0n)).toEqual([true, [-3n, 1n]]);
    expect(is_cm_j_invariant(8000n)).toEqual([true, [-8n, 1n]]);
    expect(is_cm_j_invariant(1728n)).toEqual([true, [-4n, 1n]]);
    expect(is_cm_j_invariant(54000n)).toEqual([true, [-3n, 2n]]);
    expect(is_cm_j_invariant(287496n)).toEqual([true, [-4n, 2n]]);
    expect(is_cm_j_invariant(-12288000n)).toEqual([true, [-3n, 3n]]);
  });

  // cm.py doctest: all(is_cm_j_invariant(j) == (True, (d,f)) for d,f,j in
  //                   cm_j_invariants_and_orders(QQ))
  it('agrees with cm_j_invariants_and_orders over QQ', () => {
    for (const [d, f, j] of cm_j_invariants_and_orders('QQ')) {
      expect(is_cm_j_invariant(j)).toEqual([true, [d, f]]);
    }
  });

  it('returns (false, null) for non-CM j-invariants', () => {
    expect(is_cm_j_invariant(389n)).toEqual([false, null]);
    expect(is_cm_j_invariant(1n)).toEqual([false, null]);
  });
});
