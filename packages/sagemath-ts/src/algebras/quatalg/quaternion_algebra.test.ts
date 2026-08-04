/**
 * Tests for quaternion algebras, their orders and their fractional ideals.
 *
 * Every expected value below is either
 *
 *  - copied verbatim from a doctest of
 *    `reference/sage/src/sage/algebras/quatalg/quaternion_algebra.py` /
 *    `quaternion_algebra_element.pyx` / `quaternion_algebra_cython.pyx`
 *    (marked `sage:`), or
 *  - produced by running the installed SageMath 10.3 (marked `sage 10.3:`);
 *    where SageMath 10.3 predates the vendored reference the two are noted.
 *
 * @see Reference: sage/algebras/quatalg/quaternion_algebra.py
 */

import { describe, expect, it } from 'bun:test';
import { is_squarefree } from '../../arith/misc.js';
import { QQ } from '../../rings/rational_field.js';
import {
  QuaternionAlgebra,
  type QuaternionFractionalIdeal_rational,
  basis_for_quaternion_lattice,
  hilbert_symbol_QQ,
  maxord_solve_aux_eq,
  normalize_basis_at_p,
} from './index.js';
import {
  integral_matrix_and_denom_from_rational_quaternions,
  rational_matrix_from_rational_quaternions,
  rational_quaternions_from_integral_matrix_and_denom,
} from './quaternion_algebra_cython.js';

const strs = (xs: { toString(): string }[]) => xs.map((x) => x.toString());
const lat = (I: QuaternionFractionalIdeal_rational) =>
  I.free_module()
    .basis_matrix()
    .rows()
    .map((r) => r.map(String));

describe('quaternion algebra elements', () => {
  it('multiplies according to i^2 = a, j^2 = b, ji = -ij', () => {
    // sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2)
    const Q = QuaternionAlgebra(-5n, -2n);
    const [i, j, k] = Q.gens();
    expect(i.mul(i).toString()).toBe('-5');
    expect(j.mul(j).toString()).toBe('-2');
    expect(i.mul(j).toString()).toBe('k');
    expect(j.mul(i).toString()).toBe('-k');
    expect(k.mul(k).toString()).toBe('-10');
    expect(Q.__call__([1n, 2n, 3n, 4n]).toString()).toBe('1 + 2*i + 3*j + 4*k');
    expect(Q.__call__('-3/5').toString()).toBe('-3/5');
    expect(Q.basis().map(String)).toEqual(['1', 'i', 'j', 'k']);
  });

  it('prints like SageMath', () => {
    // sage: Q.<i,j,k> = QuaternionAlgebra(-17,-19); str(i+j+k-3/4)
    const Q = QuaternionAlgebra(-17n, -19n);
    const [i, j, k] = Q.gens();
    expect(i.add(j).add(k).sub(Q.__call__('3/4')).toString()).toBe('-3/4 + i + j + k');
    expect(Q.zero().toString()).toBe('0');
    // sage: Q.<i,j,k> = QuaternionAlgebra(QQ,-5,-2); list(1/2 + 2/3*i - 3/4*j + 5/7*k)
    const Q2 = QuaternionAlgebra(-5n, -2n);
    const theta = Q2.__call__(['1/2', '2/3', '-3/4', '5/7']);
    expect(strs(theta.list())).toEqual(['1/2', '2/3', '-3/4', '5/7']);
  });

  it('computes conjugate, reduced norm/trace and inverse', () => {
    // sage: A.<i,j,k> = QuaternionAlgebra(QQ,-5,-2); (3*i - j + 2).conjugate()
    const A = QuaternionAlgebra(-5n, -2n);
    const [i, j] = A.gens();
    const a = i.scalar_mul(3n).sub(j).add(A.__call__(2n));
    expect(a.conjugate().toString()).toBe('2 - 3*i + j');

    // sage 10.3: a = 1/2 + 2/3*i - 3/4*j + 5/7*k in (-5,-2)
    const t = A.__call__(['1/2', '2/3', '-3/4', '5/7']);
    expect(t.reduced_norm().toString()).toBe('30691/3528');
    expect(t.reduced_trace().toString()).toBe('1');
    expect(t.inverse().toString()).toBe('1764/30691 - 2352/30691*i + 2646/30691*j - 2520/30691*k');
    expect(t.mul(t.inverse()).toString()).toBe('1');
    expect(t.denominator()).toBe(84n);
    expect(t.denominator_and_integer_coefficient_tuple()).toEqual([84n, 42n, 56n, -63n, 60n]);

    // sage: A.<i,j,k>=QuaternionAlgebra(-1,-2); (2 + 3*i + 4/3*j - 5*k)....
    const B = QuaternionAlgebra(-1n, -2n);
    const [bi, bj, bk] = B.gens();
    const q = B.__call__(2n)
      .add(bi.scalar_mul(3n))
      .add(bj.scalar_mul('4/3'))
      .sub(bk.scalar_mul(5n));
    expect(q.denominator_and_integer_coefficient_tuple()).toEqual([3n, 6n, 9n, 4n, -15n]);
    expect(q.integer_coefficient_tuple()).toEqual([6n, 9n, 4n, -15n]);
    expect(B.__call__(1n).is_constant()).toBe(true);
    expect(B.one().add(bi).is_constant()).toBe(false);
  });

  it('computes the reduced characteristic polynomial', () => {
    // sage: A.<i,j,k>=QuaternionAlgebra(-1,-2)
    const A = QuaternionAlgebra(-1n, -2n);
    const [i, j, k] = A.gens();
    expect(strs(i.reduced_characteristic_polynomial().coeffs as never)).toEqual(['1', '0', '1']);
    expect(strs(j.reduced_characteristic_polynomial().coeffs as never)).toEqual(['2', '0', '1']);
    expect(strs(i.add(j).reduced_characteristic_polynomial().coeffs as never)).toEqual([
      '3',
      '0',
      '1',
    ]);
    // sage: (2+j+k).reduced_trace() == 4, char poly T^2 - 4*T + 8
    const x = A.__call__(2n).add(j).add(k);
    expect(x.reduced_trace().toString()).toBe('4');
    expect(strs(x.reduced_characteristic_polynomial('T').coeffs as never)).toEqual([
      '8',
      '-4',
      '1',
    ]);
  });

  it('computes multiplication matrices and the pairing', () => {
    // sage: Q.<i,j,k> = QuaternionAlgebra(-3,-19); a = 2/3 -1/2*i + 3/5*j - 4/3*k
    const Q = QuaternionAlgebra(-3n, -19n);
    const a = Q.__call__(['2/3', '-1/2', '3/5', '-4/3']);
    expect(
      a
        .matrix()
        .rows()
        .map((r) => r.map(String))
    ).toEqual([
      ['2/3', '-1/2', '3/5', '-4/3'],
      ['3/2', '2/3', '4', '3/5'],
      ['-57/5', '-76/3', '2/3', '1/2'],
      ['76', '-57/5', '-3/2', '2/3'],
    ]);
    expect(
      a
        .matrix('left')
        .rows()
        .map((r) => r.map(String))
    ).toEqual([
      ['2/3', '-1/2', '3/5', '-4/3'],
      ['3/2', '2/3', '-4', '-3/5'],
      ['-57/5', '76/3', '2/3', '-1/2'],
      ['76', '57/5', '3/2', '2/3'],
    ]);
    // sage: A.<i,j,k>=QuaternionAlgebra(-1,-2); (1+i+j-2*k).pair(2/3+5*i-3*j+k)
    const A = QuaternionAlgebra(-1n, -2n);
    const x = A.__call__([1n, 1n, 1n, -2n]);
    const y = A.__call__(['2/3', 5n, -3n, 1n]);
    expect(x.pair(y).toString()).toBe('-26/3');
    expect(y.pair(x).toString()).toBe('-26/3');
    expect(x.conjugate().mul(y).reduced_trace().toString()).toBe('-26/3');
  });

  it('is multiplicative on the reduced norm (property)', () => {
    const Q = QuaternionAlgebra(-11n, -3n);
    const B = Q.basis();
    let seed = 12345n;
    const next = (): bigint => {
      seed = (seed * 6364136223846793005n + 1442695040888963407n) % (1n << 62n);
      return (seed % 21n) - 10n;
    };
    for (let t = 0; t < 60; t++) {
      let x = Q.zero();
      let y = Q.zero();
      for (let m = 0; m < 4; m++) {
        x = x.add(B[m]!.scalar_mul(next()));
        y = y.add(B[m]!.scalar_mul(next()));
      }
      expect(x.mul(y).reduced_norm().eq(x.reduced_norm().mul(y.reduced_norm()))).toBe(true);
      // conjugation is an anti-automorphism
      expect(x.mul(y).conjugate().eq(y.conjugate().mul(x.conjugate()))).toBe(true);
      // x satisfies its own reduced characteristic polynomial
      const p = x.mul(x).sub(x.scalar_mul(x.reduced_trace())).add(Q.__call__(x.reduced_norm()));
      expect(p.is_zero()).toBe(true);
      if (!x.is_zero()) {
        expect(x.mul(x.inverse()).eq(Q.one())).toBe(true);
      }
    }
  });
});

describe('quaternion algebras', () => {
  it('constructs from invariants and from a discriminant', () => {
    // sage: QuaternionAlgebra(-2,-3)
    expect(QuaternionAlgebra(-2n, -3n).toString()).toBe(
      'Quaternion Algebra (-2, -3) with base ring Rational Field'
    );
    // sage: Q.<i,j,k> = QuaternionAlgebra(15); Q.invariants()
    expect(strs(QuaternionAlgebra(15n).invariants())).toEqual(['-3', '5']);
    // sage: QuaternionAlgebra(QQ, -7, -21)
    expect(QuaternionAlgebra(QQ as never, -7n, -21n).toString()).toBe(
      'Quaternion Algebra (-7, -21) with base ring Rational Field'
    );
    // the factory is cached, as SageMath's UniqueFactory is
    expect(QuaternionAlgebra(-1n, -1n)).toBe(QuaternionAlgebra(-1n, -1n));
    // sage: QuaternionAlgebra(0, 1)
    expect(() => QuaternionAlgebra(0n, 1n)).toThrow(
      'defining elements of quaternion algebra (0, 1) are not invertible in Rational Field'
    );
  });

  it('constructs from ramification data over QQ', () => {
    // sage: QuaternionAlgebra(QQ, [2,3], []) == QuaternionAlgebra(6)
    expect(QuaternionAlgebra(QQ as never, [2n, 3n] as never, [] as never).toString()).toBe(
      QuaternionAlgebra(6n).toString()
    );
    expect(QuaternionAlgebra(QQ as never, [2n, 3n, 5n] as never, ['1/2'] as never).toString()).toBe(
      QuaternionAlgebra(30n).toString()
    );
    expect(() => QuaternionAlgebra(QQ as never, [2n, 3n] as never, ['1/2'] as never)).toThrow(
      'quaternion algebra over the rationals must have an even number of ramified places'
    );
  });

  it('computes discriminants and ramified places', () => {
    // sage: QuaternionAlgebra(210,-22).ramified_places()
    expect(QuaternionAlgebra(210n, -22n).ramified_primes()).toEqual([2n, 3n, 5n, 7n]);
    expect(QuaternionAlgebra(210n, -22n).discriminant()).toBe(210n);
    // sage: QuaternionAlgebra(-58, -69).ramified_primes()
    expect(QuaternionAlgebra(-58n, -69n).ramified_primes()).toEqual([3n, 23n, 29n]);
    // sage: QuaternionAlgebra(19).discriminant(); QuaternionAlgebra(-1,-1).discriminant()
    expect(QuaternionAlgebra(19n).discriminant()).toBe(19n);
    expect(QuaternionAlgebra(-1n, -1n).discriminant()).toBe(2n);
    // sage: QuaternionAlgebra(-1,-1).ramified_places()  ->  ([2], [<the real place>])
    expect(QuaternionAlgebra(-1n, -1n).ramified_places()).toEqual([[2n], ['infinity']]);
    expect(QuaternionAlgebra(1n, 1n).ramified_places()).toEqual([[], []]);
    // sage 10.3: rational invariants
    expect(QuaternionAlgebra('-1/4', -3n).discriminant()).toBe(3n);
    expect(QuaternionAlgebra('2/3', '-5/7').ramified_primes()).toEqual([2n, 7n]);
    expect(QuaternionAlgebra(-1n, '1/2').discriminant()).toBe(1n);
  });

  it('decides definiteness, division algebras and matrix rings', () => {
    // sage: QuaternionAlgebra(QQ,-5,-2).is_definite(); QuaternionAlgebra(1).is_definite()
    expect(QuaternionAlgebra(-5n, -2n).is_definite()).toBe(true);
    expect(QuaternionAlgebra(1n).is_definite()).toBe(false);
    // sage: QuaternionAlgebra(QQ,-5,-2).is_division_algebra() / .is_matrix_ring()
    expect(QuaternionAlgebra(-5n, -2n).is_division_algebra()).toBe(true);
    expect(QuaternionAlgebra(2n, 9n).is_division_algebra()).toBe(false);
    expect(QuaternionAlgebra(-5n, -2n).is_matrix_ring()).toBe(false);
    expect(QuaternionAlgebra(2n, 9n).is_matrix_ring()).toBe(true);
    expect(QuaternionAlgebra(-5n, -2n).is_totally_definite()).toBe(true);
    // sage: B = QuaternionAlgebra(-46,-87); A = QuaternionAlgebra(-58,-69); B.is_isomorphic(A)
    expect(QuaternionAlgebra(-46n, -87n).is_isomorphic(QuaternionAlgebra(-58n, -69n))).toBe(true);
    expect(QuaternionAlgebra(-1n, -1n).is_isomorphic(QuaternionAlgebra(-1n, -11n))).toBe(false);
  });

  it('computes the inner product matrix', () => {
    // sage: Q.<i,j,k> = QuaternionAlgebra(-5,-19); Q.inner_product_matrix()
    expect(
      QuaternionAlgebra(-5n, -19n)
        .inner_product_matrix()
        .rows()
        .map((r) => r.map(String))
    ).toEqual([
      ['2', '0', '0', '0'],
      ['0', '10', '0', '0'],
      ['0', '0', '38', '0'],
      ['0', '0', '0', '190'],
    ]);
  });

  it('agrees with SageMath on the ramified primes of rational invariants', () => {
    // Cross-check of hilbert_symbol over QQ: values from SageMath 10.3.
    const cases: [string, string, string][] = [
      ['47/7', '13/4', '329'],
      ['-1', '45/4', '1'],
      ['-6', '-38/3', '114'],
      ['-31/5', '13/3', '1'],
      ['-11/2', '-2', '2'],
      ['-19/12', '16/3', '38'],
    ];
    for (const [a, b, d] of cases) {
      expect(QuaternionAlgebra(a, b).discriminant().toString()).toBe(d);
    }
    // hilbert_symbol_QQ itself
    expect(hilbert_symbol_QQ(QQ.__call__('-1/4'), QQ.__call__(-1n), 2n)).toBe(-1n);
    expect(hilbert_symbol_QQ(QQ.__call__('-1/4'), QQ.__call__(-1n), 3n)).toBe(1n);
  });
});

describe('maximal orders', () => {
  it('reproduces the SageMath doctests (Pizer shortcuts)', () => {
    const expected: [bigint, bigint, string][] = [
      // sage: QuaternionAlgebra(-1,-7).maximal_order()
      [-1n, -7n, '1/2 + 1/2*j, 1/2*i + 1/2*k, j, k'],
      [-1n, -1n, '1/2 + 1/2*i + 1/2*j + 1/2*k, i, j, k'],
      [-1n, -11n, '1/2 + 1/2*j, 1/2*i + 1/2*k, j, k'],
      [-1n, -3n, '1/2 + 1/2*j, 1/2*i + 1/2*k, j, k'],
      [-3n, -1n, '1/2 + 1/2*i, 1/2*j - 1/2*k, i, -k'],
      [-2n, -5n, '1/2 + 1/2*j + 1/2*k, 1/4*i + 1/2*j + 1/4*k, j, k'],
      [-5n, -2n, '1/2 + 1/2*i - 1/2*k, 1/2*i + 1/4*j - 1/4*k, i, -k'],
      [-17n, -3n, '1/2 + 1/2*j, 1/2*i + 1/2*k, -1/3*j - 1/3*k, k'],
      [-3n, -17n, '1/2 + 1/2*i, 1/2*j - 1/2*k, -1/3*i + 1/3*k, -k'],
      [-153n, -3n, '1, 1/3*i, 1/6*i + 1/2*j, 1/2 + 1/3*j + 1/18*k'],
      [-2n, -389n, '1/2 + 1/2*j + 1/2*k, 1/4*i + 1/2*j + 1/4*k, j, k'],
    ];
    for (const [a, b, basis] of expected) {
      expect(strs(QuaternionAlgebra(a, b).maximal_order().basis()).join(', ')).toBe(basis);
    }
  });

  it('reproduces the SageMath doctests (Voight, take_shortcuts=false)', () => {
    // sage: QuaternionAlgebra(-3,-89).maximal_order(take_shortcuts=False)
    expect(
      strs(QuaternionAlgebra(-3n, -89n).maximal_order({ take_shortcuts: false }).basis()).join(', ')
    ).toBe('1, 1/2 + 1/2*i, j, 1/2 + 1/6*i + 1/2*j + 1/6*k');
    // sage: QuaternionAlgebra(1,1).maximal_order(take_shortcuts=False)
    expect(
      strs(QuaternionAlgebra(1n, 1n).maximal_order({ take_shortcuts: false }).basis()).join(', ')
    ).toBe('1, 1/2 + 1/2*i, j, 1/2*j + 1/2*k');
    // sage: QuaternionAlgebra(-22,210).maximal_order(take_shortcuts=False)
    expect(
      strs(QuaternionAlgebra(-22n, 210n).maximal_order({ take_shortcuts: false }).basis()).join(
        ', '
      )
    ).toBe('1, i, 1/2*i + 1/2*j, 1/2 + 17/22*i + 1/44*k');
  });

  it('is maximal for every squarefree discriminant below 350', () => {
    // sage: all(Q(d).maximal_order(take_shortcuts=False).is_maximal()
    //           for d in range(1, 350) if is_squarefree(d))
    for (let d = 1n; d < 350n; d++) {
      if (!is_squarefree(d)) continue;
      const A = QuaternionAlgebra(d);
      expect(A.maximal_order({ take_shortcuts: false }).is_maximal()).toBe(true);
      expect(A.maximal_order().is_maximal()).toBe(true);
      expect(A.maximal_order().discriminant().toString()).toBe(A.discriminant().toString());
    }
  });

  it('handles the invariants of Sage issues 37417 and 37217', () => {
    // sage: all(QuaternionAlgebra(a, b).maximal_order().is_maximal() for a, b in invars)
    const invars: [bigint, bigint][] = [
      [-4n, -28n],
      [-292n, -732n],
      [-48n, -564n],
      [-436n, -768n],
      [-752n, -708n],
      [885n, 545n],
      [411n, -710n],
      [-411n, 593n],
      [805n, -591n],
      [-921n, 353n],
      [409n, 96n],
      [394n, 873n],
      [353n, -722n],
      [730n, 830n],
      [-466n, -427n],
      [-213n, -630n],
      [-511n, 608n],
      [493n, 880n],
      [105n, -709n],
      [-213n, 530n],
      [97n, 745n],
    ];
    for (const [a, b] of invars) {
      expect(QuaternionAlgebra(a, b).maximal_order().is_maximal()).toBe(true);
    }
    // sage 10.3 values (the ones it can compute):
    expect(strs(QuaternionAlgebra(885n, 545n).maximal_order().basis()).join(', ')).toBe(
      '1, 1/2 + 1/2*i, j, 1/2*j + 1/10*k'
    );
    expect(strs(QuaternionAlgebra(409n, 96n).maximal_order().basis()).join(', ')).toBe(
      '1, 1/2 + 1/2*i, 1/4*j, 131/409*i + 5/48*j + 1/19632*k'
    );
  });

  it('extends a given order to a maximal one', () => {
    // sage: O = A.quaternion_order(basis=order_basis); R = A.maximal_order(order_basis=...)
    const A = QuaternionAlgebra(-292n, -732n);
    const [i, j, k] = A.gens();
    const alpha = A.__call__([1n, 2n, 3n, 4n]);
    const conj = [k, i, j].map((b) => alpha.mul(b).mul(alpha.inverse()));
    const order_basis = [...conj, A.one()];
    const O = A.quaternion_order(order_basis);
    const R = A.maximal_order({ order_basis });
    expect(O.le(R)).toBe(true);
    expect(R.is_maximal()).toBe(true);
  });

  it('computes orders with a given level', () => {
    // sage: A.<i,j,k> = QuaternionAlgebra(5); A.order_with_level(2*5*17)
    const A = QuaternionAlgebra(5n);
    const O = A.order_with_level(2n * 5n * 17n);
    expect(strs(O.basis()).join(', ')).toBe('1/2 + 1/2*j + 7/2*k, 1/2*i + 19/2*k, j + 7*k, 17*k');
    // sage: L = O.free_module(); N = A.maximal_order().free_module(); L.index_in(N) == level/5
    expect(O.free_module().index_in(A.maximal_order().free_module()).toString()).toBe('34');
  });
});

describe('quaternion orders', () => {
  it('validates the defining basis exactly as SageMath does', () => {
    // sage: A.<i,j,k> = QuaternionAlgebra(-1,-1)
    const A = QuaternionAlgebra(-1n, -1n);
    const [i, j, k] = A.gens();
    expect(() => A.quaternion_order([A.one(), i, j, i.sub(j)])).toThrow('basis must have rank 4');
    expect(() => A.quaternion_order([A.__call__(2n), i, j, k])).toThrow('lattice must contain 1');
    expect(() =>
      A.quaternion_order([A.one(), i.scalar_mul('1/2'), j.scalar_mul('1/2'), k.scalar_mul('1/2')])
    ).toThrow('given lattice must be a ring');
    // sage: Q.quaternion_order([i,j,k], check=False)
    expect(strs(A.quaternion_order([i, j, k], false).basis()).join(', ')).toBe('i, j, k');
  });

  it('computes discriminants, unit ideals and membership', () => {
    // sage: QuaternionAlgebra(-11,-1).maximal_order().discriminant()
    const Q = QuaternionAlgebra(-11n, -1n);
    const R = Q.maximal_order();
    expect(R.discriminant().toString()).toBe('11');
    expect(R.is_maximal()).toBe(true);
    // sage: R.unit_ideal()
    expect(R.unit_ideal().toString()).toBe('Fractional ideal (1/2 + 1/2*i, 1/2*j - 1/2*k, i, -k)');
    // sage: O = Q.quaternion_order([1,i,j,k]); O(1/2) raises; 1/5 in O; j/2 in O
    const [i, j] = Q.gens();
    const O = Q.quaternion_order([Q.one(), i, j, i.mul(j)]);
    expect(() => O.__call__('1/2')).toThrow('does not lie in');
    expect(O.unit_ideal().contains('1/5')).toBe(false);
    expect(O.unit_ideal().contains(j.scalar_mul('1/2'))).toBe(false);
    expect(O.unit_ideal().contains(Q.one().add(i))).toBe(true);
    // sage: Q.quaternion_order([1,-i,k,j+i*7]) == Q.quaternion_order([1,i,j,k])
    const QQ2 = QuaternionAlgebra(-1n, -19n);
    const [a, b, c] = QQ2.gens();
    expect(
      QQ2.quaternion_order([QQ2.one(), a.neg(), c, b.add(a.scalar_mul(7n))]).eq(
        QQ2.quaternion_order([QQ2.one(), a, b, c])
      )
    ).toBe(true);
  });

  it('compares orders by containment', () => {
    // sage: B = QuaternionAlgebra(-1, -11); O = B.quaternion_order([1,i,j,k])
    //       R = B.quaternion_order([1,i,(i+j)/2,(1+k)/2]); O <= R
    const B = QuaternionAlgebra(-1n, -11n);
    const [i, j, k] = B.gens();
    const O = B.quaternion_order([B.one(), i, j, k]);
    const R = B.quaternion_order([
      B.one(),
      i,
      i.add(j).scalar_mul('1/2'),
      B.one().add(k).scalar_mul('1/2'),
    ]);
    expect(O.le(R)).toBe(true);
    expect(O.ge(R)).toBe(false);
    expect(O.eq(R)).toBe(false);
    expect(O.lt(R)).toBe(true);
    expect(O.le(O)).toBe(true);
    expect(R.ge(R)).toBe(true);
  });

  it('computes left ideal bases and right orders', () => {
    // sage: A.<i,j,k> = QuaternionAlgebra(-17, -3)
    const A = QuaternionAlgebra(-17n, -3n);
    const [i, j, k] = A.gens();
    expect(
      strs(
        A.maximal_order()._left_ideal_basis([i.add(j), i.sub(j), k.scalar_mul(2n), A.__call__(3n)])
      ).join(', ')
    ).toBe('1, i, 1/2 + 1/2*j, 1/2 + 1/2*i + 1/6*j + 1/6*k');
    expect(
      strs(
        A.maximal_order()._left_ideal_basis([
          i.add(j).scalar_mul(3n),
          i.sub(j).scalar_mul(3n),
          k.scalar_mul(6n),
          A.__call__(3n),
        ])
      ).join(', ')
    ).toBe('3, 3*i, 3/2 + 3/2*j, 3/2 + 3/2*i + 1/2*j + 1/2*k');

    // sage: A.<i,j,k> = QuaternionAlgebra(17); O = A.maximal_order()
    const A17 = QuaternionAlgebra(17n);
    const O17 = A17.maximal_order();
    const basis = O17._left_ideal_basis([A17.one()]);
    expect(strs(basis).join(', ')).toBe('1, 1/2 + 1/2*i, j, 1/3*i + 1/2*j + 1/6*k');
    expect(O17._right_order_from_ideal_basis(basis).toString()).toBe(
      'Order of Quaternion Algebra (-3, -17) with base ring Rational Field with basis ' +
        '(1/2 + 1/6*i + 1/3*k, 1/3*i + 2/3*k, 1/2*j + 1/2*k, k)'
    );
    const [x, y] = A17.gens();
    const basis2 = O17._left_ideal_basis([x.mul(y).sub(y)]);
    expect(strs(basis2).join(', ')).toBe('34, 17 + 17*i, 2*j, 17 + 17/3*i + j + 1/3*k');
    expect(O17._right_order_from_ideal_basis(basis2).toString()).toBe(
      'Order of Quaternion Algebra (-3, -17) with base ring Rational Field with basis ' +
        '(1/2 + 1/6*i + 1/3*k, 1/3*i + 2/3*k, 1/2*j + 1/2*k, k)'
    );
  });

  it('intersects orders', () => {
    // sage: R = QuaternionAlgebra(-11,-1).maximal_order(); R.intersection(R)
    const R = QuaternionAlgebra(-11n, -1n).maximal_order();
    expect(strs(R.intersection(R).basis()).join(', ')).toBe('1/2 + 1/2*i, i, 1/2*j + 1/2*k, k');
  });

  it('computes quadratic forms of orders', () => {
    // sage 10.3: R = QuaternionAlgebra(-11,-1).maximal_order()
    const R = QuaternionAlgebra(-11n, -1n).maximal_order();
    const q = R.quadratic_form();
    expect([0, 1, 2, 3].map((r) => [0, 1, 2, 3].map((c) => q.get(r, c).value))).toEqual([
      [6n, 0n, 11n, 0n],
      [0n, 6n, 0n, 11n],
      [11n, 0n, 22n, 0n],
      [0n, 11n, 0n, 22n],
    ]);
    const [m, basis] = R.ternary_quadratic_form(true) as [
      { rows(): { toString(): string }[][] },
      { toString(): string }[],
    ];
    expect(m.rows().map((r) => r.map(String))).toEqual([
      ['22', '0', '0'],
      ['0', '24', '44'],
      ['0', '44', '88'],
    ]);
    expect(strs(basis)).toEqual(['i', 'j + k', '2*k']);
  });

  it('computes an isomorphism between conjugate maximal orders', () => {
    // sage: Quat.<i,j,k> = QuaternionAlgebra(-1, -19)
    //       O0 = Quat.quaternion_order([1, i, (i+j)/2, (1+k)/2])
    //       O1 = Quat.quaternion_order([1, 667*i, 1/2+j/2+9*i, (222075/2*i+333*j+k/2)/667])
    const Quat = QuaternionAlgebra(-1n, -19n);
    const [i, j, k] = Quat.gens();
    const O0 = Quat.quaternion_order([
      Quat.one(),
      i,
      i.add(j).scalar_mul('1/2'),
      Quat.one().add(k).scalar_mul('1/2'),
    ]);
    const O1 = Quat.quaternion_order([
      Quat.one(),
      i.scalar_mul(667n),
      Quat.__call__('1/2').add(j.scalar_mul('1/2')).add(i.scalar_mul(9n)),
      i.scalar_mul('222075/2').add(j.scalar_mul(333n)).add(k.scalar_mul('1/2')).scalar_mul('1/667'),
    ]);
    const gamma = O0.isomorphism_to(O1, { conjugator: true }) as ReturnType<typeof Quat.__call__>;
    // gamma^-1 * O0 * gamma == O1
    const conj = Quat.quaternion_order(O0.basis().map((b) => gamma.inverse().mul(b).mul(gamma)));
    expect(conj.eq(O1)).toBe(true);
    // sage 10.3: O0.isomorphism_to(O1, conjugator=True) == -36*i - j + k (up to sign)
    expect(gamma.toString() === '36*i + j - k' || gamma.toString() === '-36*i - j + k').toBe(true);
    const iso = O0.isomorphism_to(O1) as { im_gens: { toString(): string }[] };
    expect(strs(iso.im_gens)).toEqual([
      '629/667*i + 36/667*j - 36/667*k',
      '684/667*i - 648/667*j - 19/667*k',
      '-684/667*i - 19/667*j - 648/667*k',
    ]);
    // sage: O1.isomorphism_to(O2) for a non-isomorphic order raises
    const O2 = Quat.quaternion_order([Quat.one(), i, j, k]);
    expect(() => O2.isomorphism_to(O0)).toThrow('only implemented for maximal orders');
  });
});

describe('quaternion fractional ideals', () => {
  it('constructs left and right ideals', () => {
    // sage: Q.<i,j,k> = QuaternionAlgebra(-11,-1); R = Q.maximal_order()
    const Q = QuaternionAlgebra(-11n, -1n);
    const [i, j, k] = Q.gens();
    const R = Q.maximal_order();
    expect(
      R.left_ideal(
        R.basis().map((a) => a.scalar_mul(2n)),
        true,
        { is_basis: true }
      ).toString()
    ).toBe('Fractional ideal (1 + i, 2*i, j + k, 2*k)');
    expect(
      R.left_ideal(
        R.basis().map((a) => a.mul(i.add(j))),
        true,
        { is_basis: true }
      ).toString()
    ).toBe('Fractional ideal (1/2 + 1/2*i + 1/2*j + 13/2*k, i + j, 6*j + 6*k, 12*k)');
    // sage: R.left_ideal([i+j])
    expect(R.left_ideal([i.add(j)]).toString()).toBe(
      'Fractional ideal (12, 6 + 6*i, i + j, 13/2 + 1/2*i + 1/2*j + 1/2*k)'
    );
    // sage: R.right_ideal([i+j])
    expect(R.right_ideal([i.add(j)]).toString()).toBe(
      'Fractional ideal (12, 6 + 6*i, i + j, 11/2 + 1/2*i + 1/2*j + 1/2*k)'
    );
    // sage: R.left_ideal([2, 1+j]) == R*2 + R*(1+j)
    expect(
      R.left_ideal([Q.__call__(2n), Q.one().add(j)]).eq(R.mul(2n).add(R.mul(Q.one().add(j))))
    ).toBe(true);
    expect(
      R.right_ideal([Q.__call__(2n), Q.one().add(j)]).eq(R.rmul(2n).add(R.rmul(Q.one().add(j))))
    ).toBe(true);
  });

  it('computes ideal basis matrices and free modules', () => {
    // sage: QuaternionAlgebra(-11,-1).maximal_order().unit_ideal().basis_matrix()
    expect(
      QuaternionAlgebra(-11n, -1n)
        .maximal_order()
        .unit_ideal()
        .basis_matrix()
        .rows()
        .map((r) => r.map(String))
    ).toEqual([
      ['1/2', '1/2', '0', '0'],
      ['0', '1', '0', '0'],
      ['0', '0', '1/2', '1/2'],
      ['0', '0', '0', '1'],
    ]);
    // sage: R.<i,j,k> = QuaternionAlgebra(-1, -13); I = R.ideal([2+i, 3*i, 5*j, j+k])
    const R = QuaternionAlgebra(-1n, -13n);
    const [i, j, k] = R.gens();
    const I = R.ideal([R.__call__(2n).add(i), i.scalar_mul(3n), j.scalar_mul(5n), j.add(k)]);
    expect(I.toString()).toBe('Fractional ideal (2 + i, 3*i, j + k, 5*k)');
    expect(
      I.free_module()
        .basis_matrix()
        .rows()
        .map((r) => r.map(String))
    ).toEqual([
      ['2', '1', '0', '0'],
      ['0', '3', '0', '0'],
      ['0', '0', '1', '1'],
      ['0', '0', '0', '5'],
    ]);
  });

  it('computes left and right orders of an ideal', () => {
    // sage: R.<i,j,k> = QuaternionAlgebra(-1,-11)
    //       I = R.ideal([2 + 2*j + 140*k, 2*i + 4*j + 150*k, 8*j + 104*k, 152*k])
    const R = QuaternionAlgebra(-1n, -11n);
    const I = R.ideal([
      R.__call__([2n, 0n, 2n, 140n]),
      R.__call__([0n, 2n, 4n, 150n]),
      R.__call__([0n, 0n, 8n, 104n]),
      R.__call__([0n, 0n, 0n, 152n]),
    ]);
    expect(I._compute_order('left').toString()).toBe(
      'Order of Quaternion Algebra (-1, -11) with base ring Rational Field with basis ' +
        '(1/2 + 1/2*j + 35*k, 1/4*i + 1/2*j + 75/4*k, j + 32*k, 38*k)'
    );
    expect(I._compute_order('right').toString()).toBe(
      'Order of Quaternion Algebra (-1, -11) with base ring Rational Field with basis ' +
        '(1/2 + 1/2*j + 16*k, 1/2*i + 11/2*k, j + 13*k, 19*k)'
    );
    expect(I.left_order().discriminant().toString()).toBe('209');
    expect(I.right_order().discriminant().toString()).toBe('209');
    expect(I.is_integral()).toBe(true);
  });

  it('computes norms, conjugates, products, sums and intersections', () => {
    // sage: A.<i,j,k> = QuaternionAlgebra(-1,-1)
    //       R = A.ideal([i,j,k,1/2 + 1/2*i + 1/2*j + 1/2*k]); R.norm()
    const A = QuaternionAlgebra(-1n, -1n);
    const [i, j, k] = A.gens();
    const Rideal = A.ideal([i, j, k, A.__call__(['1/2', '1/2', '1/2', '1/2'])]);
    expect(Rideal.norm().toString()).toBe('1');
    // sage: [J.norm() for J in R.cyclic_right_subideals(3)]
    expect(strs(Rideal.cyclic_right_subideals(3n).map((J) => J.norm()))).toEqual([
      '3',
      '3',
      '3',
      '3',
    ]);

    // sage 10.3: Q.<i,j,k> = QuaternionAlgebra(-1,-11); O = Q.maximal_order()
    const Q = QuaternionAlgebra(-1n, -11n);
    const [qi, qj] = Q.gens();
    const O = Q.maximal_order();
    const I = O.left_ideal([qi.add(qj)]);
    expect(I.norm().toString()).toBe('12');
    expect(I.conjugate().norm().toString()).toBe('12');
    expect(I.mul(I).norm().toString()).toBe('12');
    const I2 = O.left_ideal([Q.__call__(2n).add(qj)]);
    // I + J is the unit ideal here (sage 10.3)
    expect(I.add(I2).eq(O.unit_ideal())).toBe(true);
    expect(I.multiply_by_conjugate(I2).eq(O.unit_ideal())).toBe(true);
    // the intersection is a common sublattice
    const inter = I.intersection(I2);
    expect(inter.free_module().is_submodule(I.free_module())).toBe(true);
    expect(inter.free_module().is_submodule(I2.free_module())).toBe(true);
    // sage 10.3 lattice of the intersection
    expect(lat(inter)).toEqual([
      ['3/2', '3/2', '39/2', '69/2'],
      ['0', '3', '15', '24'],
      ['0', '0', '30', '30'],
      ['0', '0', '0', '60'],
    ]);
  });

  it('computes theta series, reduced bases and minimal elements', () => {
    // sage 10.3: I = QuaternionAlgebra(-1,-11).maximal_order().left_ideal([i+j])
    const Q = QuaternionAlgebra(-1n, -11n);
    const [i, j] = Q.gens();
    const I = Q.maximal_order().left_ideal([i.add(j)]);
    expect(I.theta_series_vector(12)).toEqual([
      1n,
      4n,
      4n,
      8n,
      20n,
      16n,
      32n,
      16n,
      36n,
      28n,
      40n,
      4n,
    ]);
    expect(I.theta_series(6).toString()).toContain('1 + 4*q + 4*q^2 + 8*q^3 + 20*q^4 + 16*q^5');
    expect(I.minimal_element().toString()).toBe('i + j');
    expect(I.minimal_element().reduced_norm().toString()).toBe('12');
    expect(I.is_principal(false)).toBe(true);
    const [ok, alpha] = I.is_principal(true);
    expect(ok).toBe(true);
    expect(alpha!.reduced_norm().toString()).toBe('12');
    // a reduced basis consists of four elements of the ideal with increasing norm
    const red = I.reduced_basis();
    expect(red.length).toBe(4);
    for (const x of red) expect(I.contains(x)).toBe(true);
    for (let t = 0; t + 1 < red.length; t++) {
      expect(red[t]!.reduced_norm().le(red[t + 1]!.reduced_norm())).toBe(true);
    }

    // sage: Quat.<i,j,k> = QuaternionAlgebra(-3,-101); O = Quat.maximal_order()
    const Quat = QuaternionAlgebra(-3n, -101n);
    const O = Quat.maximal_order();
    expect(O.toString()).toBe(
      'Order of Quaternion Algebra (-3, -101) with base ring Rational Field with basis ' +
        '(1/2 + 1/2*i, 1/2*j - 1/2*k, -1/3*i + 1/3*k, -k)'
    );
    expect(O.mul(5n).minimal_element().toString()).toBe('5');
    // sage: alpha = 1/2 + 1/6*i + j + 55/3*k; I = O*141 + O*alpha; I.norm()
    const alphaQ = Quat.__call__(['1/2', '1/6', 1n, '55/3']);
    const Ia = O.mul(141n).add(O.mul(alphaQ));
    expect(Ia.norm().toString()).toBe('141');
    // sage: I.minimal_element()
    expect(Ia.minimal_element().toString()).toBe('13/2 - 7/6*i + j + 2/3*k');
    expect(Ia.minimal_element().reduced_norm().toString()).toBe('282');
  });

  it('decides integrality, primitivity and equivalence', () => {
    // sage: A.<i,j,k> = QuaternionAlgebra(QQ, -1,-11)
    //       I = A.ideal([1/2 + 1/2*i + 1/2*j + 3/2*k, i + k, j + k, 2*k])
    const A = QuaternionAlgebra(-1n, -11n);
    const [i, j, k] = A.gens();
    const I = A.ideal([
      A.__call__(['1/2', '1/2', '1/2', '3/2']),
      i.add(k),
      j.add(k),
      k.scalar_mul(2n),
    ]);
    expect(I.is_primitive()).toBe(true);
    const [Jp, g] = I.primitive_decomposition();
    expect(g).toBe(1n);
    expect(Jp.eq(I)).toBe(true);
    // sage: (2*I).is_primitive() is False
    const twoI = I.scale(2n);
    expect(twoI.is_primitive()).toBe(false);
    const [J2, g2] = twoI.primitive_decomposition();
    expect(g2).toBe(2n);
    expect(J2.scale(2n).eq(twoI)).toBe(true);

    // sage: J = A.ideal([7/2 + 7/2*i + 49/2*j + 91/2*k, 7*i + 21*k, 35*j + 35*k, 70*k])
    const J = A.ideal([
      A.__call__(['7/2', '7/2', '49/2', '91/2']),
      i.scalar_mul(7n).add(k.scalar_mul(21n)),
      j.scalar_mul(35n).add(k.scalar_mul(35n)),
      k.scalar_mul(70n),
    ]);
    const [Jeq, gj] = J.primitive_decomposition();
    expect(gj).toBe(7n);
    expect(Jeq.scale(7n).eq(J)).toBe(true);
    expect(Jeq.is_primitive()).toBe(true);
    // sage: Jequiv
    expect(Jeq.toString()).toBe(
      'Fractional ideal (10, 5 + 5*i, 3 + j, 13/2 + 7/2*i + 1/2*j + 1/2*k)'
    );

    // sage: I = R.ideal([1/2 + 2*j + 140*k, ...]).is_integral() is False
    const R = QuaternionAlgebra(-1n, -11n);
    const nonint = R.ideal([
      R.__call__(['1/2', 0n, 2n, 140n]),
      R.__call__([0n, 2n, 4n, 150n]),
      R.__call__([0n, 0n, 8n, 104n]),
      R.__call__([0n, 0n, 0n, 152n]),
    ]);
    expect(nonint.is_integral()).toBe(false);

    // right/left equivalence of an ideal with a scaled copy
    const O = A.maximal_order();
    const I0 = O.left_ideal([i.add(j)]);
    expect(I0.is_right_equivalent(I0)).toBe(true);
    const [eq, cert] = I0.is_right_equivalent(I0, 10, true) as [boolean, ReturnType<typeof A.one>];
    expect(eq).toBe(true);
    expect(cert.mul(cert.conjugate()).is_constant()).toBe(true);
    expect(I0.conjugate().is_left_equivalent(I0.conjugate())).toBe(true);
  });

  it('computes cyclic right subideals', () => {
    // sage 10.3: Q.<i,j,k> = QuaternionAlgebra(-2,-5); I = Q.ideal([Q(1),i,j,k])
    const Q = QuaternionAlgebra(-2n, -5n);
    const [i, j, k] = Q.gens();
    const I = Q.ideal([Q.one(), i, j, k]);
    const subs = I.cyclic_right_subideals(3n);
    // sage: I.cyclic_right_subideals(3)
    expect(subs.map(String)).toEqual([
      'Fractional ideal (3, 3*i, 2 + j, i + k)',
      'Fractional ideal (3, 3*i, 1 + j, 2*i + k)',
      'Fractional ideal (3, 2 + i, 3*j, 2*j + k)',
      'Fractional ideal (3, 1 + i, 3*j, j + k)',
    ]);
    // each subideal has index p^2 and is contained in I
    for (const J of subs) {
      expect(J.free_module().is_submodule(I.free_module())).toBe(true);
      expect(J.free_module().index_in(I.free_module()).toString()).toBe('9');
    }
    // sage: I.cyclic_right_subideals(3)[0].cyclic_right_subideals(3) raises
    expect(() => subs[0]!.cyclic_right_subideals(3n)).toThrow('general algorithm not implemented');
  });

  it('computes pushforwards and pullbacks', () => {
    // sage: B = QuaternionAlgebra(419)
    //       I1 = B.ideal([1/2 + 3/2*j + 2*k, 1/2*i + j + 3/2*k, 3*j, 3*k])
    //       I2 = B.ideal([1/2 + 9/2*j, 1/2*i + 9/2*k, 5*j, 5*k])
    const B = QuaternionAlgebra(419n);
    const [i, j, k] = B.gens();
    const I1 = B.ideal([
      B.__call__(['1/2', 0n, '3/2', 2n]),
      B.__call__([0n, '1/2', 1n, '3/2']),
      j.scalar_mul(3n),
      k.scalar_mul(3n),
    ]);
    const I2 = B.ideal([
      B.__call__(['1/2', 0n, '9/2', 0n]),
      B.__call__([0n, '1/2', 0n, '9/2']),
      j.scalar_mul(5n),
      k.scalar_mul(5n),
    ]);
    expect(I1.left_order().eq(I2.left_order())).toBe(true);
    // sage: I1.pushforward(I2, side='left')
    const I3 = I1.pushforward(I2, 'left');
    expect(I3.toString()).toBe(
      'Fractional ideal (3, 15*i, 3/2 + 10*i + 1/2*j, 1 + 9/10*i + 1/10*k)'
    );
    // sage: I3.left_order() == I2.right_order(); I3.pullback(I2, side='left') == I1
    expect(I3.left_order().eq(I2.right_order())).toBe(true);
    expect(I3.pullback(I2, 'left').eq(I1)).toBe(true);
    // sage: O0.unit_ideal().pushforward(O0.unit_ideal(), "left")
    const U = B.maximal_order().unit_ideal();
    expect(U.pushforward(U, 'left').toString()).toBe(
      'Fractional ideal (1, i, 1/2 + 1/2*j, 1/2*i + 1/2*k)'
    );
    // sage: O0.unit_ideal().pullback(O0.unit_ideal(), "left")
    expect(U.pullback(U, 'left').toString()).toBe(
      'Fractional ideal (1/2 + 1/2*j, 1/2*i + 1/2*k, j, k)'
    );
    expect(() => U.pullback(U)).toThrow(
      'self and J have same left and right orders, side of pullback must be specified'
    );
    // sage: O0.unit_ideal().pushforward(O0.unit_ideal()) raises
    expect(() => U.pushforward(U)).toThrow(
      'self and J have same left and right orders, side of pushforward must be specified'
    );
    expect(() => I1.pushforward(I2, 'right')).toThrow('self and J must have the same right orders');
    expect(() => I1.pushforward(I1, 'left')).toThrow('self and J must have coprime norms');
    // sage: I2b.pullback(I1)
    const I2b = B.ideal([
      B.__call__(['1/2', 0n, '15/2', 2n]),
      B.__call__([0n, '1/6', '43/3', '5/2']),
      j.scalar_mul(15n),
      k.scalar_mul(5n),
    ]);
    // sage: I2.pullback(I1)
    expect(I2b.pullback(I1).toString()).toBe(
      'Fractional ideal (1/2 + 5/2*j + 2*k, 1/2*i + 3*j + 5/2*k, 5*j, 5*k)'
    );
    // sage: I2.pullback(I1, side='right')
    expect(() => I2b.pullback(I1, 'right')).toThrow(
      'right order of self should be left order of J'
    );
    // sage: I2.conjugate().pullback(I1.conjugate(), side='right')
    expect(I2b.conjugate().pullback(I1.conjugate(), 'right').toString()).toBe(
      'Fractional ideal (1/2 + 5/2*j + 3*k, 1/2*i + 3*j + 5/2*k, 5*j, 5*k)'
    );
    // sage: I1.pullback(I1.conjugate(), side='left')
    expect(() => I1.pullback(I1.conjugate(), 'left')).toThrow('self and J must have coprime norms');
  });

  it('scales ideals on either side', () => {
    // sage: I.scale(1) == I  (issue 32245)
    const Q = QuaternionAlgebra(-1n, -19n);
    const [i, j, k] = Q.gens();
    const I = Q.ideal([Q.one(), i, j, k]);
    expect(I.scale(1n).eq(I)).toBe(true);
    expect(() => I.scale(0n)).toThrow('the scaling factor must be nonzero');
    // sage: (5+i-j)*I == I.scale(5+i-j, left=True)
    const a = Q.__call__(5n).add(i).sub(j);
    const left = I.scale(a, true);
    const right = I.scale(a, false);
    for (const b of I.basis()) {
      expect(left.contains(a.mul(b))).toBe(true);
      expect(right.contains(b.mul(a))).toBe(true);
    }
    // Ol * J subset J and J * Or subset J (issue 32726)
    for (const J of [left, right]) {
      const Ol = J.left_order();
      const Or = J.right_order();
      expect(Ol.unit_ideal().mul(J).free_module().is_submodule(J.free_module())).toBe(true);
      expect(J.mul(Or.unit_ideal()).free_module().is_submodule(J.free_module())).toBe(true);
    }
  });

  it('satisfies the norm identities of ideals (property)', () => {
    for (const d of [2n, 3n, 5n, 7n, 11n, 13n, 101n]) {
      const A = QuaternionAlgebra(d);
      const O = A.maximal_order();
      const B = O.basis();
      for (const cs of [
        [1n, 2n, 3n, 4n],
        [0n, 1n, -1n, 2n],
        [3n, -1n, 0n, 5n],
      ]) {
        let alpha = A.zero();
        for (let t = 0; t < 4; t++) alpha = alpha.add(B[t]!.scalar_mul(cs[t]!));
        if (alpha.is_zero()) continue;
        const N = alpha.reduced_norm();
        const I = O.mul(alpha).add(O.mul(N));
        // I * Ibar = N(I) * (left order of I)
        const prod = I.multiply_by_conjugate(I);
        expect(prod.eq(I.left_order().unit_ideal().scale(I.norm(), true))).toBe(true);
        // the left and right orders both have the discriminant of a maximal order
        expect(I.left_order().discriminant().toString()).toBe(A.discriminant().toString());
        expect(I.right_order().discriminant().toString()).toBe(A.discriminant().toString());
        // the norm is multiplicative along compatible orders
        expect(prod.norm().eq(I.norm().mul(I.norm()))).toBe(true);
        // theta series counts elements of small norm; the constant term is 1
        expect(I.theta_series_vector(4)[0]).toBe(1n);
      }
    }
  });
});

describe('mod p splitting data', () => {
  it('reproduces the SageMath doctests', () => {
    // sage: Q = QuaternionAlgebra(-15, -19); Q.modp_splitting_data(7)
    const Q = QuaternionAlgebra(-15n, -19n);
    const [I, J, K] = Q.modp_splitting_data(7n);
    expect(I).toEqual([
      [0n, 6n],
      [1n, 0n],
    ]);
    expect(J).toEqual([
      [6n, 1n],
      [1n, 1n],
    ]);
    expect(K).toEqual([
      [6n, 6n],
      [6n, 1n],
    ]);
    // sage: I,J,K = Q.modp_splitting_data(23); I; J
    const [I23, J23] = Q.modp_splitting_data(23n);
    expect(I23).toEqual([
      [0n, 8n],
      [1n, 0n],
    ]);
    expect(J23).toEqual([
      [19n, 2n],
      [17n, 4n],
    ]);
    // sage: Q.modp_splitting_data(5) / (2)
    expect(() => Q.modp_splitting_data(5n)).toThrow(
      'algorithm for computing local splittings not implemented in general ' +
        '(currently require the first invariant to be coprime to p)'
    );
    expect(() => Q.modp_splitting_data(2n)).toThrow('p must be odd');
    expect(() => QuaternionAlgebra(-1n, -11n).modp_splitting_data(11n)).toThrow(
      'p (=11) must be an unramified prime'
    );
    expect(() => Q.modp_splitting_data(4n)).toThrow('p (=4) must be prime');
  });

  it('produces a ring homomorphism to 2x2 matrices mod p', () => {
    // sage: Q.<i,j,k> = QuaternionAlgebra(-1, -7); f = Q.modp_splitting_map(13)
    //       a = 2+i-j+3*k; b = 7+2*i-4*j+k; f(a*b) == f(a)*f(b)
    const Q = QuaternionAlgebra(-1n, -7n);
    const [i, j, k] = Q.gens();
    const f = Q.modp_splitting_map(13n);
    const a = Q.__call__(2n).add(i).sub(j).add(k.scalar_mul(3n));
    const b = Q.__call__(7n).add(i.scalar_mul(2n)).sub(j.scalar_mul(4n)).add(k);
    const mul = (X: bigint[][], Y: bigint[][]) => [
      [
        (X[0]![0]! * Y[0]![0]! + X[0]![1]! * Y[1]![0]!) % 13n,
        (X[0]![0]! * Y[0]![1]! + X[0]![1]! * Y[1]![1]!) % 13n,
      ],
      [
        (X[1]![0]! * Y[0]![0]! + X[1]![1]! * Y[1]![0]!) % 13n,
        (X[1]![0]! * Y[0]![1]! + X[1]![1]! * Y[1]![1]!) % 13n,
      ],
    ];
    expect(f(a.mul(b))).toEqual(mul(f(a), f(b)));
    // sage 10.3 values
    expect(f(a).flat()).toEqual([12n, 4n, 6n, 5n]);
    expect(f(b).flat()).toEqual([2n, 6n, 10n, 12n]);
  });

  it('splits the algebra at many unramified primes (property)', () => {
    for (const [a, b] of [
      [-1n, -1n],
      [-1n, -7n],
      [-2n, -5n],
      [-11n, -1n],
    ] as [bigint, bigint][]) {
      const A = QuaternionAlgebra(a, b);
      for (const p of [3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n]) {
        if (A.discriminant() % p === 0n) continue;
        let I: bigint[][];
        let J: bigint[][];
        let K: bigint[][];
        try {
          [I, J, K] = A.modp_splitting_data(p);
        } catch {
          continue; // i^2 = 0 mod p is not implemented upstream either
        }
        const mul = (X: bigint[][], Y: bigint[][]) => [
          [
            (((X[0]![0]! * Y[0]![0]! + X[0]![1]! * Y[1]![0]!) % p) + p) % p,
            (((X[0]![0]! * Y[0]![1]! + X[0]![1]! * Y[1]![1]!) % p) + p) % p,
          ],
          [
            (((X[1]![0]! * Y[0]![0]! + X[1]![1]! * Y[1]![0]!) % p) + p) % p,
            (((X[1]![0]! * Y[0]![1]! + X[1]![1]! * Y[1]![1]!) % p) + p) % p,
          ],
        ];
        const scalar = (c: bigint) => [
          [((c % p) + p) % p, 0n],
          [0n, ((c % p) + p) % p],
        ];
        expect(mul(I, I)).toEqual(scalar(a));
        expect(mul(J, J)).toEqual(scalar(b));
        expect(mul(I, J)).toEqual(K);
        expect(mul(J, I)).toEqual(K.map((r) => r.map((e) => (p - e) % p)));
      }
    }
  });
});

describe('lattice utilities', () => {
  it('computes a basis for a quaternion lattice', () => {
    // sage: A.<i,j,k> = QuaternionAlgebra(-1,-7)
    //       basis_for_quaternion_lattice([i+j, i-j, 2*k, A(1/3)])
    const A = QuaternionAlgebra(-1n, -7n);
    const [i, j, k] = A.gens();
    expect(
      strs(basis_for_quaternion_lattice([i.add(j), i.sub(j), k.scalar_mul(2n), A.__call__('1/3')]))
    ).toEqual(['1/3', '2*i', 'i + j', '2*k']);
    expect(strs(basis_for_quaternion_lattice([A.__call__(1n), i, j, k]))).toEqual([
      '1',
      'i',
      'j',
      'k',
    ]);
    // sage 10.3 (reverse=False, the old default)
    expect(
      strs(
        basis_for_quaternion_lattice(
          [i.add(j), i.sub(j), k.scalar_mul(2n), A.__call__('1/3')],
          false
        )
      )
    ).toEqual(['1/3', 'i + j', '2*j', '2*k']);
    expect(basis_for_quaternion_lattice([])).toEqual([]);
  });

  it('normalizes a basis at p', () => {
    // sage: A.<i,j,k> = QuaternionAlgebra(-1, -1); normalize_basis_at_p([A(1), i, j, k], 2)
    const A = QuaternionAlgebra(-1n, -1n);
    expect(normalize_basis_at_p(A.basis().slice(), 2n).map(([v, n]) => `(${v}, ${n})`)).toEqual([
      '(1, 0)',
      '(i, 0)',
      '(j, 0)',
      '(k, 0)',
    ]);
    // sage: A.<i,j,k> = QuaternionAlgebra(210); normalize_basis_at_p([A(1), i, j, k], 2)
    const A210 = QuaternionAlgebra(210n);
    expect(normalize_basis_at_p(A210.basis().slice(), 2n).map(([v, n]) => `(${v}, ${n})`)).toEqual([
      '(1, 0)',
      '(i, 1)',
      '(j, 1)',
      '(k, 2)',
    ]);
    // sage: A.<i,j,k> = QuaternionAlgebra(286)
    //       normalize_basis_at_p([A(1), k, 1/2*j + 1/2*k, 1/2 + 1/2*i + 1/2*k], 5)
    const A286 = QuaternionAlgebra(286n);
    const [i3, j3, k3] = A286.gens();
    expect(
      normalize_basis_at_p(
        [
          A286.one(),
          k3,
          j3.add(k3).scalar_mul('1/2'),
          A286.one().add(i3).add(k3).scalar_mul('1/2'),
        ],
        5n
      ).map(([v, n]) => `(${v}, ${n})`)
    ).toEqual(['(1, 0)', '(1/2*j + 1/2*k, 0)', '(-5/6*j + 1/6*k, 1)', '(1/2*i, 1)']);
    // sage: A.<i,j,k> = QuaternionAlgebra(-1,-7) (issue 37217)
    const A17 = QuaternionAlgebra(-1n, -7n);
    const [i4, j4, k4] = A17.gens();
    expect(
      normalize_basis_at_p(
        [A17.one(), k4, j4, A17.one().add(i4).add(j4).add(k4).scalar_mul('1/2')],
        2n
      ).map(([v, n]) => `(${v}, ${n})`)
    ).toEqual([
      '(1, 0)',
      '(1/2 + 1/2*i + 1/2*j + 1/2*k, 0)',
      '(-34/105*i - 463/735*j + 71/105*k, 1)',
      '(1/7*i - 8/49*j + 1/7*k, 1)',
    ]);
  });

  it('solves the auxiliary equation at 2', () => {
    // sage: for a in [1,3]: for b in [1,2,3]: maxord_solve_aux_eq(a, b, 2)
    const expected: Record<string, [bigint, bigint, bigint]> = {
      '1,1': [1n, 1n, 1n],
      '1,2': [1n, 0n, 0n],
      '1,3': [1n, 0n, 0n],
      '3,1': [1n, 1n, 1n],
      '3,2': [1n, 0n, 1n],
      '3,3': [1n, 1n, 1n],
    };
    for (const a of [1n, 3n]) {
      for (const b of [1n, 2n, 3n]) {
        const [y, z, w] = maxord_solve_aux_eq(a, b, 2n);
        expect([y, z, w]).toEqual(expected[`${a},${b}`]);
        // sage: assert mod(1 - a*y^2 - b*z^2 + a*b*w^2, 4) == 0
        expect((((1n - a * y * y - b * z * z + a * b * w * w) % 4n) + 4n) % 4n).toBe(0n);
        expect(((y % 4n) + 4n) % 4n === 1n || ((y % 4n) + 4n) % 4n === 3n).toBe(true);
      }
    }
    expect(() => maxord_solve_aux_eq(1n, 1n, 3n)).toThrow(
      'algorithm only implemented over ZZ at the moment'
    );
  });

  it('converts between quaternions and matrices (quaternion_algebra_cython)', () => {
    // sage: A.<i,j,k>=QuaternionAlgebra(-4,-5)
    //       integral_matrix_and_denom_from_rational_quaternions([i/2,1/3+j+k])
    const A = QuaternionAlgebra(-4n, -5n);
    const [i, j, k] = A.gens();
    const v = [i.scalar_mul('1/2'), A.__call__('1/3').add(j).add(k)];
    const [M, d] = integral_matrix_and_denom_from_rational_quaternions(v);
    expect(d).toBe(6n);
    expect([0, 1].map((r) => [0, 1, 2, 3].map((c) => M.get(r, c).value))).toEqual([
      [0n, 3n, 0n, 0n],
      [2n, 0n, 6n, 6n],
    ]);
    const [Mr, dr] = integral_matrix_and_denom_from_rational_quaternions(v, true);
    expect(dr).toBe(6n);
    expect([0, 1].map((r) => [0, 1, 2, 3].map((c) => Mr.get(r, c).value))).toEqual([
      [6n, 6n, 0n, 2n],
      [0n, 0n, 3n, 0n],
    ]);
    // sage: rational_matrix_from_rational_quaternions([i/2,1/3+j+k])
    expect(
      rational_matrix_from_rational_quaternions(v)
        .rows()
        .map((r) => r.map(String))
    ).toEqual([
      ['0', '1/2', '0', '0'],
      ['1/3', '0', '1', '1'],
    ]);
    expect(
      rational_matrix_from_rational_quaternions(v, true)
        .rows()
        .map((r) => r.map(String))
    ).toEqual([
      ['1', '1', '0', '1/3'],
      ['0', '0', '1/2', '0'],
    ]);
    // sage: A.<i,j,k>=QuaternionAlgebra(-1,-2)
    //       f(A, matrix([[1,2,3,4],[-1,2,-4,3]]), 3)
    const A2 = QuaternionAlgebra(-1n, -2n);
    const [Z] = integral_matrix_and_denom_from_rational_quaternions([
      A2.__call__([1n, 2n, 3n, 4n]),
      A2.__call__([-1n, 2n, -4n, 3n]),
    ]);
    expect(strs(rational_quaternions_from_integral_matrix_and_denom(A2, Z, 3n))).toEqual([
      '1/3 + 2/3*i + j + 4/3*k',
      '-1/3 + 2/3*i - 4/3*j + k',
    ]);
  });
});
