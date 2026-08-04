/**
 * Tests for free module elements (vectors).
 *
 * The expected values are the ones SageMath produces; the corresponding
 * doctests are quoted above each test.
 *
 * @see Reference: sage/modules/free_module_element.pyx
 */

import { describe, expect, it } from 'bun:test';
import { Rational } from '../rings/rational.js';
import { FreeModule } from './free_module.js';
import { vector, zeroVector } from './free_module_element.js';

const ZZ = {
  zero: () => 0n,
  one: () => 1n,
  is_field: () => false,
  toString: () => 'Integer Ring',
};

const r = (n: number, d: number = 1) => new Rational(BigInt(n), BigInt(d));

describe('FreeModuleElement indexing', () => {
  it('wraps negative indices Python-style', () => {
    // sage: v = vector([1,2,3]); v[-1] == 3; v[-3] == 1
    const v = vector([1n, 2n, 3n]);
    expect(v.getItem(-1)).toBe(3n);
    expect(v.getItem(-2)).toBe(2n);
    expect(v.getItem(-3)).toBe(1n);
  });

  it('raises IndexError out of range', () => {
    // sage: vector([1,2,3])[3]
    // IndexError: vector index out of range
    const v = vector([1n, 2n, 3n]);
    expect(() => v.getItem(3)).toThrow('vector index out of range');
    expect(() => v.getItem(-4)).toThrow('vector index out of range');
    try {
      v.getItem(3);
      throw new Error('expected an IndexError');
    } catch (e) {
      expect((e as Error).name).toBe('IndexError');
    }
  });

  it('assigns through negative indices', () => {
    // sage: v = vector([1,2,3]); v[-1] = 9; v == (1,2,9)
    const v = vector([1n, 2n, 3n]);
    v.setItem(-1, 9n);
    expect(v.list()).toEqual([1n, 2n, 9n]);
    expect(() => v.setItem(3, 1n)).toThrow('vector index out of range');
  });

  it('refuses to modify an immutable vector', () => {
    const v = vector([1n, 2n, 3n]);
    v.setImmutable();
    expect(() => v.setItem(0, 5n)).toThrow('vector is immutable');
  });
});

describe('FreeModuleElement.norm', () => {
  it('computes the Euclidean norm exactly', () => {
    // sage: v = vector(ZZ, [3,4]); v.norm() == 5 with parent Rational Field
    expect(vector([3n, 4n]).norm()).toBe(5n);
    expect(vector([1n, 2n, 2n]).norm()).toBe(3n);
  });

  it('computes the taxicab norm', () => {
    // sage: vector([1,-2,3]).norm(1) == 6
    expect(vector([1n, -2n, 3n]).norm(1)).toBe(6n);
  });

  it('computes the infinity norm', () => {
    // sage: v = vector([1,2,-3]); v.norm(Infinity) == 3
    expect(vector([1n, 2n, -3n]).norm(Number.POSITIVE_INFINITY)).toBe(3n);
  });

  it('falls back to double precision for irrational norms', () => {
    // sage: vector([1,2,-3]).norm() == sqrt(14); vector([1,2,-3]).norm(5) == 276^(1/5)
    expect(Number(vector([1n, 2n, -3n]).norm())).toBeCloseTo(Math.sqrt(14), 12);
    expect(Number(vector([1n, 2n, -3n]).norm(5))).toBeCloseTo(276 ** 0.2, 12);
  });

  it('handles rational entries exactly', () => {
    // sage: vector(QQ, [3/5, 4/5]).norm() == 1
    expect(String(vector([r(3, 5), r(4, 5)]).norm())).toBe('1');
  });

  it('rejects p < 1', () => {
    // sage: vector(QQ,[1,2]).norm(0.99)
    // ValueError: 0.990000000000000 is not greater than or equal to 1
    expect(() => vector([1n, 2n]).norm(0.99)).toThrow('is not greater than or equal to 1');
  });
});

describe('FreeModuleElement.normalized', () => {
  it('divides by the exact norm', () => {
    // sage: v = vector([3,4]); v.normalized() == (3/5, 4/5)
    const u = vector([3n, 4n]).normalized();
    expect(u.list().map(String)).toEqual(['3/5', '4/5']);
  });

  it('supports other p', () => {
    // sage: v = vector(QQ, [4,1,3,2]); sum(v.normalized(1)) == 1
    const u = vector([4n, 1n, 3n, 2n]).normalized(1);
    expect(u.list().map(String)).toEqual(['2/5', '1/10', '3/10', '1/5']);
    const total = (u.list() as Rational[]).reduce((a, b) => a.add(b), Rational.zero());
    expect(String(total)).toBe('1');
  });

  it('throws for the zero vector', () => {
    expect(() => vector([0n, 0n, 0n]).normalized()).toThrow();
  });
});

describe('FreeModuleElement.crossProduct', () => {
  it('computes the three dimensional cross product', () => {
    // sage: v = vector([1,2,3]); w = vector([0,5,-9]); v.cross_product(w) == (-33,9,5)
    const v = vector([1n, 2n, 3n]);
    const w = vector([0n, 5n, -9n]);
    const u = v.crossProduct(w);
    expect(u.list()).toEqual([-33n, 9n, 5n]);
    expect(u.dotProduct(v)).toBe(0n);
    expect(u.dotProduct(w)).toBe(0n);
    expect(v.crossProduct(v).list()).toEqual([0n, 0n, 0n]);
  });

  it('computes the seven dimensional cross product', () => {
    // sage: u = vector(QQ, [1,-1/3,57,-9,56/4,-4,1])
    //       v = vector(QQ, [37,55,-99/57,9,-12,11/3,4/98])
    //       u.cross_product(v) == (1394815/2793, -2808401/2793, 39492/49,
    //                              -48737/399, -9151880/2793, 62513/2793, -326603/171)
    const u = vector([r(1), r(-1, 3), r(57), r(-9), r(56, 4), r(-4), r(1)]);
    const v = vector([r(37), r(55), r(-99, 57), r(9), r(-12), r(11, 3), r(4, 98)]);
    expect(u.crossProduct(v).list().map(String)).toEqual([
      '1394815/2793',
      '-2808401/2793',
      '39492/49',
      '-48737/399',
      '-9151880/2793',
      '62513/2793',
      '-326603/171',
    ]);
  });

  it('is anticommutative in seven dimensions', () => {
    // sage: u.cross_product(v) + v.cross_product(u) == 0
    const u = vector([1n, -1n, 57n, -9n, 14n, -4n, 1n]);
    const v = vector([37n, 55n, -99n, 9n, -12n, 11n, 4n]);
    expect(u.crossProduct(v).add(v.crossProduct(u)).list()).toEqual([0n, 0n, 0n, 0n, 0n, 0n, 0n]);
  });

  it('is orthogonal to both arguments in seven dimensions', () => {
    const u = vector([2n, 6n, -7n, -9n, -7n, 12n, 9n]);
    const v = vector([22n, -7n, -9n, 12n, 15n, 15n, 11n]);
    const c = u.crossProduct(v);
    expect(c.dotProduct(u)).toBe(0n);
    expect(c.dotProduct(v)).toBe(0n);
  });

  it('raises a TypeError for other degrees', () => {
    // sage: vector(range(7)).cross_product(vector(range(3)))
    // TypeError: Cross product only defined for vectors of length three or
    //            seven, not (7 and 3)
    const u = vector([0n, 1n, 2n, 3n, 4n, 5n, 6n]);
    const v = vector([0n, 1n, 2n]);
    expect(() => u.crossProduct(v)).toThrow(
      'Cross product only defined for vectors of length three or seven, not (7 and 3)'
    );
    expect(() => vector([1n, 2n]).crossProduct(vector([1n, 2n]))).toThrow(TypeError);
  });
});

describe('FreeModuleElement.innerProduct', () => {
  it('is the dot product without an inner product matrix', () => {
    // sage: v = vector(QQ,[1,2,3]); w = vector(QQ,[-1,2,-3]); v.inner_product(w) == -6
    const v = vector([1n, 2n, 3n]);
    const w = vector([-1n, 2n, -3n]);
    expect(v.innerProduct(w)).toBe(-6n);
    expect(v.innerProduct(w)).toBe(v.dotProduct(w));
  });

  it('applies the inner product matrix of the parent', () => {
    // sage: ipm = matrix(ZZ,[[2,0,-1],[0,2,0],[-1,0,6]])
    //       M = FreeModule(ZZ, 3, inner_product_matrix=ipm)
    //       v = M([1,0,0]); v.inner_product(v) == 2
    const ipm = [
      [2n, 0n, -1n],
      [0n, 2n, 0n],
      [-1n, 0n, 6n],
    ];
    const M = FreeModule(ZZ, 3, { innerProductMatrix: ipm });
    const v = M.createElement([1n, 0n, 0n]);
    expect(v.innerProduct(v)).toBe(2n);

    // sage: w = M([1,3,-1]); v = M([2,-4,5]); w.row()*ipm*v.column() == w.inner_product(v)
    const w = M.createElement([1n, 3n, -1n]);
    const u = M.createElement([2n, -4n, 5n]);
    expect(w.innerProduct(u)).toBe(-53n);
  });

  it('passes the inner product matrix through to submodules', () => {
    // sage: K = M.span_of_basis([[0,-1/2,-1/2],[0,1/2,-1/2],[2,0,0]])
    //       (K.0).inner_product(K.0) == 2
    const ipm = [
      [2n, 0n, -1n],
      [0n, 2n, 0n],
      [-1n, 0n, 6n],
    ];
    const M = FreeModule(ZZ, 3, { innerProductMatrix: ipm });
    const K = (
      M as unknown as { spanOfBasis: (b: unknown[]) => { gen: (i: number) => never } }
    ).spanOfBasis([
      M.createElement([r(0), r(-1, 2), r(-1, 2)]),
      M.createElement([r(0), r(1, 2), r(-1, 2)]),
      M.createElement([r(2), r(0), r(0)]),
    ]);
    const b0 = K.gen(0) as unknown as { innerProduct: (v: unknown) => unknown };
    expect(String(b0.innerProduct(K.gen(0)))).toBe('2');
  });
});

describe('FreeModuleElement basics', () => {
  it('adds, subtracts, negates and scales', () => {
    const v = vector([1n, 2n, 3n]);
    const w = vector([4n, 5n, 6n]);
    expect(v.add(w).list()).toEqual([5n, 7n, 9n]);
    expect(w.sub(v).list()).toEqual([3n, 3n, 3n]);
    expect(v.neg().list()).toEqual([-1n, -2n, -3n]);
    expect(v.mul(2n).list()).toEqual([2n, 4n, 6n]);
    expect(v.pairwiseProduct(w).list()).toEqual([4n, 10n, 18n]);
  });

  it('reports support, weight and zeroness', () => {
    const v = vector([0n, 2n, 0n, 3n]);
    expect(v.support()).toEqual([1, 3]);
    expect(v.hammingWeight()).toBe(2);
    expect(v.isZero()).toBe(false);
    expect(zeroVector(ZZ, 4).isZero()).toBe(true);
  });
});
