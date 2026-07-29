/**
 * Tests for Galois group functionality
 */

import { describe, expect, it } from 'bun:test';
import { NotImplementedError } from '../../errors.js';
import {
  GaloisGroup,
  type GaloisGroupElement,
  GaloisSubgroup,
  Permutation,
  galois_group,
} from './galois_group.js';
import { Rational } from '../rational.js';
import {
  CyclotomicField,
  NumberField,
  NumberFieldConstructor,
  QuadraticField,
  RationalPolynomial,
} from './number_field.js';

describe('GaloisGroup', () => {
  describe('construction', () => {
    it('should create Galois group of quadratic field', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      expect(G.degree()).toBe(2);
      expect(G.order()).toBe(2n);
    });

    it('should create Galois group of imaginary quadratic field', () => {
      const K = QuadraticField.create(-1n);
      const G = galois_group(K);
      expect(G.degree()).toBe(2);
      expect(G.order()).toBe(2n);
    });

    it('should create Galois group of cyclotomic field Q(zeta_3)', () => {
      const K = CyclotomicField.create(3n);
      const G = galois_group(K);
      expect(G.degree()).toBe(2);
      expect(G.order()).toBe(2n);
    });

    it('should create Galois group of cyclotomic field Q(zeta_5)', () => {
      const K = CyclotomicField.create(5n);
      const G = galois_group(K);
      expect(G.degree()).toBe(4);
      expect(G.order()).toBe(4n);
    });
  });

  describe('elements', () => {
    it('should list elements of quadratic field Galois group', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      const elements = G.list();
      expect(elements.length).toBe(2);
    });

    it('should have identity element', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      const id = G.identity();
      expect(id.is_identity()).toBe(true);
      expect(id.order()).toBe(1n);
    });

    it('should have non-identity element with order 2 for quadratic', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      const elements = G.list();
      const nonId = elements.find((e) => !e.is_identity());
      expect(nonId).toBeDefined();
      expect(nonId!.order()).toBe(2n);
    });
  });

  describe('generators', () => {
    it('should return generators for quadratic field', () => {
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const gens = G.gens();
      expect(gens.length).toBeGreaterThanOrEqual(1);
    });

    it('should return number of generators', () => {
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      expect(G.ngens()).toBeGreaterThanOrEqual(1);
    });

    it('should return specific generator by index', () => {
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const gen0 = G.gen(0);
      expect(gen0).toBeDefined();
    });
  });

  describe('group properties', () => {
    it('should identify abelian group', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      // All groups of order 2 are abelian
      expect(G.is_abelian()).toBe(true);
    });

    it('should identify cyclic group', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      // Order 2 groups are cyclic
      expect(G.is_cyclic()).toBe(true);
    });

    it('should return transitive for Galois group', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      expect(G.is_transitive()).toBe(true);
    });
  });

  describe('random element', () => {
    it('should return a valid element', () => {
      const K = QuadraticField.create(3n);
      const G = galois_group(K);
      const elem = G.random_element();
      expect(elem.parent()).toBe(G);
    });
  });

  describe('iterator', () => {
    it('should iterate over all elements', () => {
      const K = QuadraticField.create(7n);
      const G = galois_group(K);
      const elements: GaloisGroupElement[] = [];
      for (const g of G) {
        elements.push(g);
      }
      expect(elements.length).toBe(2);
    });
  });

  describe('conjugacy classes', () => {
    it('should compute conjugacy classes', () => {
      const K = QuadraticField.create(11n);
      const G = galois_group(K);
      const classes = G.conjugacy_classes();
      // For abelian groups, each element is its own conjugacy class
      expect(classes.length).toBe(2);
    });
  });
});

describe('GaloisGroupElement', () => {
  describe('arithmetic', () => {
    it('should multiply elements', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      const elements = G.list();
      const g = elements[0]!;
      const h = elements[1]!;
      const product = g.mul(h);
      expect(product.parent()).toBe(G);
    });

    it('should compute inverse', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      const elements = G.list();
      for (const g of elements) {
        const gInv = g.inv();
        const product = g.mul(gInv);
        expect(product.is_identity()).toBe(true);
      }
    });

    it('should compute powers', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      const elements = G.list();
      const nonId = elements.find((e) => !e.is_identity())!;

      // g^0 = identity
      expect(nonId.pow(0n).is_identity()).toBe(true);

      // g^1 = g
      expect(nonId.pow(1n).eq(nonId)).toBe(true);

      // g^2 = identity for order 2 element
      expect(nonId.pow(2n).is_identity()).toBe(true);

      // g^(-1) = g for order 2 element
      expect(nonId.pow(-1n).eq(nonId)).toBe(true);
    });
  });

  describe('equality', () => {
    it('should check equality correctly', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      const id = G.identity();
      const id2 = G.identity();
      expect(id.eq(id2)).toBe(true);

      const elements = G.list();
      if (elements.length === 2) {
        expect(elements[0]!.eq(elements[1]!)).toBe(false);
      }
    });
  });

  describe('order', () => {
    it('should compute element order', () => {
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const id = G.identity();
      expect(id.order()).toBe(1n);

      const nonId = G.list().find((e) => !e.is_identity())!;
      expect(nonId.order()).toBe(2n);
    });
  });

  describe('permutation', () => {
    it('should return as permutation', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      const id = G.identity();
      const perm = id.as_permutation();
      expect(perm).toEqual([0, 1]);
    });

    it('should return non-identity permutation', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      const nonId = G.list().find((e) => !e.is_identity())!;
      const perm = nonId.as_permutation();
      expect(perm).toEqual([1, 0]);
    });
  });

  describe('string representation', () => {
    it('should produce cycle notation for identity', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      const id = G.identity();
      expect(id.toString()).toBe('()');
    });

    it('should produce cycle notation for transposition', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      const nonId = G.list().find((e) => !e.is_identity())!;
      expect(nonId.toString()).toBe('(1,2)');
    });
  });
});

describe('Cyclotomic Field Galois Groups', () => {
  it('should have correct order for Q(zeta_4)', () => {
    const K = CyclotomicField.create(4n);
    const G = galois_group(K);
    // Gal(Q(zeta_4)/Q) = (Z/4Z)* = {1, 3} has order 2
    expect(G.order()).toBe(2n);
  });

  it('should have correct order for Q(zeta_7)', () => {
    const K = CyclotomicField.create(7n);
    const G = galois_group(K);
    // Gal(Q(zeta_7)/Q) = (Z/7Z)* has order phi(7) = 6
    expect(G.order()).toBe(6n);
  });

  it('should be abelian for cyclotomic fields', () => {
    const K = CyclotomicField.create(5n);
    const G = galois_group(K);
    expect(G.is_abelian()).toBe(true);
  });

  it('should be cyclic for prime cyclotomic fields', () => {
    const K = CyclotomicField.create(5n);
    const G = galois_group(K);
    // (Z/pZ)* is cyclic for prime p
    expect(G.is_cyclic()).toBe(true);
  });
});

describe('Cubic Field Galois Groups', () => {
  it('should refuse to invent a group for the non-Galois x^3 - 2', () => {
    // Sage returns S_3 (order 6) via the Galois closure; the port has no
    // galoisinit, so it must throw rather than report the cyclic group of
    // order 3 that it used to fabricate.
    const poly = RationalPolynomial.fromBigInts([-2n, 0n, 0n, 1n]);
    const K = new NumberField(poly, 'a');
    const G = galois_group(K);
    expect(G.degree()).toBe(3);
    expect(() => G.order()).toThrow(NotImplementedError);
    expect(() => G.list()).toThrow(NotImplementedError);
  });

  it('should create Galois group for cyclic cubic', () => {
    // x^3 - 3x + 1 is a cyclic cubic (totally real, Galois)
    const poly = RationalPolynomial.fromBigInts([1n, -3n, 0n, 1n]);
    const K = new NumberField(poly, 'a');
    const G = galois_group(K);
    expect(G.degree()).toBe(3);
    expect(G.order()).toBe(3n);
    expect(G.is_cyclic()).toBe(true);
    expect(G.is_abelian()).toBe(true);
    // the nontrivial elements must act as ring homomorphisms
    const a = K.gen();
    for (const g of G.list()) {
      expect(g.__call__(a.mul(a)).eq(g.__call__(a).mul(g.__call__(a)))).toBe(true);
      expect(g.__call__(a.add(K.one())).eq(g.__call__(a).add(K.one()))).toBe(true);
    }
  });
});

describe('Galois groups of cyclotomic fields are (Z/nZ)* (audit H25)', () => {
  it('is C2 x C2 for Q(zeta_8), not cyclic of order 4', () => {
    const K = CyclotomicField.create(8n);
    const G = galois_group(K);
    expect(G.order()).toBe(4n);
    expect(G.is_abelian()).toBe(true);
    expect(G.is_cyclic()).toBe(false);
    for (const g of G.list()) {
      expect(g.order() === 1n || g.order() === 2n).toBe(true);
    }
  });

  it('sends zeta to a primitive n-th root of unity', () => {
    for (const n of [5n, 7n, 8n, 12n]) {
      const K = CyclotomicField.create(n);
      const G = galois_group(K);
      const zeta = K.gen();
      for (const g of G.list()) {
        const img = g.__call__(zeta);
        expect(img.pow(n).is_one()).toBe(true);
        for (const p of [2n, 3n, 5n]) {
          if (n % p === 0n) {
            expect(img.pow(n / p).is_one()).toBe(false);
          }
        }
      }
    }
  });

  it('is cyclic exactly when (Z/nZ)* is', () => {
    expect(galois_group(CyclotomicField.create(5n)).is_cyclic()).toBe(true);
    expect(galois_group(CyclotomicField.create(7n)).is_cyclic()).toBe(true);
    expect(galois_group(CyclotomicField.create(9n)).is_cyclic()).toBe(true);
    expect(galois_group(CyclotomicField.create(8n)).is_cyclic()).toBe(false);
    expect(galois_group(CyclotomicField.create(12n)).is_cyclic()).toBe(false);
  });
});

describe('Frobenius at p = 2 (audit H16 follow-up)', () => {
  it('uses the real prime decomposition, not a Legendre symbol', () => {
    // disc(Q(sqrt 5)) = 5 = 5 mod 8, so 2 is inert and Frob_2 is nontrivial.
    const K5 = new NumberField(RationalPolynomial.fromBigInts([-5n, 0n, 1n]), 'a');
    expect(K5.galois_group().frobenius(2n).is_identity()).toBe(false);
    // disc(Q(sqrt -7)) = -7 = 1 mod 8, so 2 splits and Frob_2 is trivial.
    const K7 = new NumberField(RationalPolynomial.fromBigInts([7n, 0n, 1n]), 'a');
    expect(K7.galois_group().frobenius(2n).is_identity()).toBe(true);
    // 2 ramifies in Q(i).  SageMath's message names the *ideal*:
    //   sage: G.artin_symbol(K.primes_above(2)[0])
    //   ValueError: Fractional ideal (...) is ramified
    // (galois_group.py:767 artin_symbol); the old expectation 'Prime 2 is
    // ramified' was this port's own wording and has been corrected.
    expect(() => QuadraticField.create(-1n).galois_group().frobenius(2n)).toThrow('is ramified');
  });
});

describe('Integration tests', () => {
  it('should work end-to-end for Q(sqrt(5))', () => {
    // Create field
    const K = QuadraticField.create(5n);
    expect(K.degree()).toBe(2);

    // Create Galois group
    const G = galois_group(K);
    expect(G.order()).toBe(2n);

    // List elements
    const elements = G.list();
    expect(elements.length).toBe(2);

    // Check group structure
    expect(G.is_abelian()).toBe(true);
    expect(G.is_cyclic()).toBe(true);

    // Test element operations
    const id = G.identity();
    const conj = elements.find((e) => !e.is_identity())!;

    expect(id.mul(conj).eq(conj)).toBe(true);
    expect(conj.mul(conj).is_identity()).toBe(true);
    expect(conj.inv().eq(conj)).toBe(true);
  });

  it('should support group operations for Q(zeta_5)', () => {
    const K = CyclotomicField.create(5n);
    const G = galois_group(K);

    // Order should be phi(5) = 4
    expect(G.order()).toBe(4n);

    // List elements
    const elements = G.list();
    expect(elements.length).toBe(4);

    // Check closure under multiplication
    for (const g of elements) {
      for (const h of elements) {
        const product = g.mul(h);
        const found = elements.some((e) => e.eq(product));
        expect(found).toBe(true);
      }
    }

    // Check inverses
    for (const g of elements) {
      const gInv = g.inv();
      const product = g.mul(gInv);
      expect(product.is_identity()).toBe(true);
    }
  });
});

describe('Galois Group Arithmetic Features', () => {
  describe('automorphisms', () => {
    it('should return automorphisms for quadratic field', () => {
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const auts = G.automorphisms();
      expect(auts.length).toBe(2);
    });

    it('should return automorphisms with __call__ method', () => {
      const K = QuadraticField.create(2n);
      const G = galois_group(K);
      const auts = G.automorphisms();
      expect(auts[0]).toHaveProperty('__call__');
    });
  });

  describe('subgroups', () => {
    it('should compute subgroups for quadratic field', () => {
      const K = QuadraticField.create(3n);
      const G = galois_group(K);
      const subs = G.subgroups();
      // Order 2 group has 2 subgroups: trivial and whole group
      expect(subs.length).toBe(2);
    });
  });

  describe('fixed_field', () => {
    it('should return whole field for trivial subgroup', () => {
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const fixed = G.fixed_field([G.identity()]);
      expect(fixed).toBe(K);
    });

    it('should return Q for whole group', () => {
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const fixed = G.fixed_field(G.list());
      expect(fixed).toBe('Q');
    });
  });

  describe('decomposition_group', () => {
    it('should compute decomposition group for split prime', () => {
      // In Q(sqrt(5)), 11 splits since (5/11) = 1
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const D = G.decomposition_group(11n);
      // Split prime: D_P = {1}
      expect(D.order()).toBe(1n);
    });

    it('should compute decomposition group for inert prime', () => {
      // In Q(sqrt(5)), 3 is inert since (5/3) = -1
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const D = G.decomposition_group(3n);
      // Inert prime: D_P = G
      expect(D.order()).toBe(2n);
    });

    it('should compute decomposition group for ramified prime', () => {
      // In Q(sqrt(5)), 5 is ramified
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const D = G.decomposition_group(5n);
      // Ramified prime: D_P = G
      expect(D.order()).toBe(2n);
    });
  });

  describe('inertia_group', () => {
    it('should return trivial group for unramified prime', () => {
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const I = G.inertia_group(11n);
      // Unramified: I_P = {1}
      expect(I.order()).toBe(1n);
    });

    it('should return whole group for ramified prime', () => {
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const I = G.inertia_group(5n);
      // Ramified: I_P = G
      expect(I.order()).toBe(2n);
    });
  });

  describe('frobenius', () => {
    it('should return identity for split prime', () => {
      // In Q(sqrt(5)), 11 splits
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const frob = G.frobenius(11n);
      expect(frob.is_identity()).toBe(true);
    });

    it('should return non-identity for inert prime', () => {
      // In Q(sqrt(5)), 3 is inert
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const frob = G.frobenius(3n);
      expect(frob.is_identity()).toBe(false);
    });

    it('should throw for ramified prime', () => {
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      expect(() => G.frobenius(5n)).toThrow();
    });
  });

  describe('artin_symbol', () => {
    it('should equal frobenius for unramified primes', () => {
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      const artin = G.artin_symbol(11n);
      const frob = G.frobenius(11n);
      expect(artin.eq(frob)).toBe(true);
    });

    it('should throw for ramified prime', () => {
      const K = QuadraticField.create(5n);
      const G = galois_group(K);
      expect(() => G.artin_symbol(5n)).toThrow();
    });
  });
});

describe('GaloisGroupElement fixed_field', () => {
  it('should return Q for non-identity in quadratic field', () => {
    const K = QuadraticField.create(7n);
    const G = galois_group(K);
    const elements = G.list();
    const nonId = elements.find((e) => !e.is_identity())!;
    expect(nonId.fixed_field()).toBe('Q');
  });

  it('should return whole field for identity', () => {
    const K = QuadraticField.create(7n);
    const G = galois_group(K);
    const id = G.identity();
    expect(id.fixed_field()).toBe(K);
  });
});

/* ==================================================================== */
/* Degree > 2: everything below goes through PARI galoisinit /          */
/* galoisfixedfield, and every expected value is a SageMath doctest     */
/* copied verbatim from reference/sage/src/sage/rings/number_field/     */
/* galois_group.py.                                                     */
/* ==================================================================== */

describe('GaloisGroup for degree > 2 (galoisinit)', () => {
  // sage: L.<a> = NumberField(x^4 + 1)
  // sage: G = L.galois_group()
  // sage: H = G.decomposition_group(L.primes_above(3)[0])
  // sage: H.fixed_field()
  // (Number Field in a0 with defining polynomial x^2 + 2 with a0 = a^3 + a, ...)
  it('x^4 + 1: the decomposition group at 3 fixes Q(sqrt(-2)), generated by a^3 + a', () => {
    const L = NumberFieldConstructor([1n, 0n, 0n, 0n, 1n], 'a');
    const G = L.galois_group();
    expect(G.is_galois()).toBe(true);
    expect(G.order()).toBe(4n);
    const P3 = L.primes_above(3n)[0]!;
    const H = G.decomposition_group(P3);
    expect(H.order()).toBe(2n);
    const d = G.fixed_field_data(H.list());
    expect(d.polynomial).toEqual([2n, 0n, 1n]); // x^2 + 2
    expect(d.gen!.toString()).toBe('a^3 + a');
    // the generator really is a root of the returned polynomial
    const g = d.gen!;
    expect(g.mul(g).add(L.__call__(2n)).is_zero()).toBe(true);
  });

  // sage: G = NumberField(x^5 - 5*x^2 - 3, 'a').galois_group()
  // sage: sigma, tau = G.gens(); H = G.subgroup([tau])
  // sage: H.fixed_field(polred=False)
  // (Number Field in a0 with defining polynomial x^2 + 84375
  //   with a0 = 5*ac^5 + 25*ac^3, ...)
  // sage: G.splitting_field()
  // Number Field in ac with defining polynomial x^10 + 10*x^8 + 25*x^6 + 3375
  it('the splitting field of x^5 - 5x^2 - 3: <order 5> fixes x^2 + 84375 = (5ac^5+25ac^3)^2', () => {
    const L = NumberFieldConstructor([3375n, 0n, 0n, 0n, 0n, 0n, 25n, 0n, 10n, 0n, 1n], 'ac');
    const G = L.galois_group();
    expect(G.order()).toBe(10n);
    const order5 = G.list().filter((e) => e.order() === 5n);
    expect(order5.length).toBe(4);
    for (const s of order5) {
      const d = s.fixed_field_data();
      expect(d.polynomial).toEqual([84375n, 0n, 1n]);
      expect(d.gen!.toString()).toBe('5*ac^5 + 25*ac^3');
    }
  }, 60000);

  // sage: K.<a> = NumberField(x^4 - 2*x^2 + 2, 'b').galois_closure()
  //   -> x^8 - 20*x^6 + 104*x^4 - 40*x^2 + 1156
  // sage: P = K.ideal([17, a^2]); G.decomposition_group(P)
  // Subgroup generated by [(1,8)(2,7)(3,6)(4,5)] of (Galois group 8T4 ...)
  it('x^8 - 20x^6 + 104x^4 - 40x^2 + 1156: D(17, a^2) has order 2', () => {
    const K = NumberFieldConstructor(
      [1156n, 0n, -40n, 0n, 104n, 0n, -20n, 0n, 1n],
      'a'
    );
    const G = K.galois_group();
    expect(G.order()).toBe(8n);
    const a = K.gen();
    const P = K.ideal(K.__call__(17n), a.mul(a));
    expect(P.is_prime()).toBe(true);
    expect(P.residue_class_degree()).toBe(2n);
    const D = G.decomposition_group(P);
    expect(D.order()).toBe(2n);
    // the nontrivial element is an involution moving every root
    const s = D.list().find((e) => !e.is_identity())!;
    expect(s.order()).toBe(2n);
  }, 60000);

  // sage: K.<b> = NumberField(x^8 - 20*x^6 + 104*x^4 - 40*x^2 + 1156)
  // sage: P = K.primes_above(2)[0]; G.ramification_breaks(P)  ->  {1, 3, 5}
  // sage: min(G.ramification_group(P, i).order() / G.ramification_group(P, i+1).order()
  // ....:     for i in G.ramification_breaks(P))  ->  2
  it('ramification breaks above 2 are {1, 3, 5} with every jump of index 2', () => {
    const K = NumberFieldConstructor(
      [1156n, 0n, -40n, 0n, 104n, 0n, -20n, 0n, 1n],
      'a'
    );
    const G = K.galois_group();
    const P2 = K.primes_above(2n)[0]!;
    expect(G.ramification_breaks(P2)).toEqual([1, 3, 5]);
    const ratios = G.ramification_breaks(P2).map(
      (i) =>
        Number(G.ramification_group(P2, i).order()) /
        Number(G.ramification_group(P2, i + 1).order())
    );
    expect(Math.min(...ratios)).toBe(2);
    // the whole filtration
    expect([-1, 0, 1, 2, 3, 4, 5, 6].map((v) => Number(G.ramification_group(P2, v).order()))).toEqual(
      [8, 8, 8, 4, 4, 2, 2, 1]
    );
  }, 120000);

  // sage: K.<b> = NumberField(x^3 - 3, 'a').galois_closure()  ->  x^6 + 243
  // sage: P = K.primes_above(3)[0]
  // sage: G.ramification_group(P, 3)  ->  Subgroup generated by [(1,2,4)(3,5,6)]
  // sage: G.ramification_group(P, 5)  ->  Subgroup generated by [()]
  // sage: s = hom(K, K, 1/18*b^4 - 1/2*b); G(s).ramification_degree(P)  ->  4
  it('x^6 + 243: G_3 has order 3, G_5 is trivial, and ramification_degree is 4', () => {
    const K = NumberFieldConstructor([243n, 0n, 0n, 0n, 0n, 0n, 1n], 'b');
    const G = K.galois_group();
    expect(G.order()).toBe(6n);
    const P = K.primes_above(3n)[0]!;
    expect(G.ramification_group(P, 3).order()).toBe(3n);
    expect(G.ramification_group(P, 5).order()).toBe(1n);
    for (const s of G.list()) {
      if (s.is_identity()) continue;
      expect(s.ramification_degree(P)).toBe(s.order() === 3n ? 4n : 1n);
    }
  });

  // sage: K.<b> = NumberField(x^4 - 2*x^2 + 2, 'a').galois_closure()
  // sage: sorted([G.artin_symbol(P) for P in K.primes_above(7)])  # random
  // [(1,4)(2,3)(5,8)(6,7), (1,4)(2,3)(5,8)(6,7),
  //  (1,5)(2,6)(3,7)(4,8), (1,5)(2,6)(3,7)(4,8)]
  // sage: G.artin_symbol(K.primes_above(2)[0])
  // ValueError: Fractional ideal (...) is ramified
  it('Artin symbols above 7 are two involutions, each twice; above 2 it is ramified', () => {
    const K = NumberFieldConstructor(
      [1156n, 0n, -40n, 0n, 104n, 0n, -20n, 0n, 1n],
      'a'
    );
    const G = K.galois_group();
    const above7 = K.primes_above(7n);
    expect(above7.length).toBe(4);
    const symbols = above7.map((P) => G.artin_symbol(P));
    for (const s of symbols) expect(s.order()).toBe(2n);
    const distinct = new Set(symbols.map((s) => s.toString()));
    expect(distinct.size).toBe(2);
    // the Frobenius generates the decomposition group, of order f = 2
    for (let i = 0; i < above7.length; i++) {
      const P = above7[i]!;
      expect(P.residue_class_degree()).toBe(2n);
      expect(G.decomposition_group(P).order()).toBe(2n);
    }
    expect(() => G.artin_symbol(K.primes_above(2n)[0]!)).toThrow('is ramified');
  }, 120000);

  // sage: QuadraticField(-7,'c').galois_group().artin_symbol(13)  ->  (1,2)
  it('QuadraticField(-7).artin_symbol(13) is the nontrivial element', () => {
    const K = QuadraticField.create(-7n, 'c');
    expect(K.galois_group().artin_symbol(13n).toString()).toBe('(1,2)');
  });

  // sage: K.<b> = NumberField(x^2 - 3, 'a'); G = K.galois_group()
  // sage: G.inertia_group(K.primes_above(2)[0])  ->  [(1,2)]
  // sage: G.inertia_group(K.primes_above(5)[0])  ->  [()]
  it('x^2 - 3: inertia at 2 is the full group, at 5 it is trivial', () => {
    const K = NumberFieldConstructor([-3n, 0n, 1n], 'b');
    const G = K.galois_group();
    expect(G.inertia_group(K.primes_above(2n)[0]!).order()).toBe(2n);
    expect(G.inertia_group(K.primes_above(5n)[0]!).order()).toBe(1n);
  });

  it('|D(P)| = e*f and |I(P)| = e for every prime of every test field', () => {
    const fields: bigint[][] = [
      [1n, 0n, 0n, 0n, 1n], // x^4 + 1
      [243n, 0n, 0n, 0n, 0n, 0n, 1n], // x^6 + 243
      [1n, 1n, 1n, 1n, 1n, 1n, 1n], // Phi_7
      [1n, -3n, 0n, 1n], // x^3 - 3x + 1
    ];
    for (const f of fields) {
      const K = NumberFieldConstructor(f, 'a');
      const G = K.galois_group();
      for (const p of [2n, 3n, 5n, 7n]) {
        const dec = K.decomposition(p);
        for (const [P, e] of dec) {
          const fdeg = P.residue_class_degree();
          expect(G.decomposition_group(P).order()).toBe(e * fdeg);
          expect(G.inertia_group(P).order()).toBe(e);
          // G_v is a decreasing chain of subgroups of D
          let prev = Number(G.decomposition_group(P).order());
          for (let v = 0; v <= 3; v++) {
            const o = Number(G.ramification_group(P, v).order());
            expect(o).toBeLessThanOrEqual(prev);
            expect(prev % o).toBe(0);
            prev = o;
          }
        }
      }
    }
  }, 120000);

  it('the fixed field of every subgroup has degree [G:H] and its generator is a root', () => {
    const fields: bigint[][] = [
      [1n, 0n, 0n, 0n, 1n], // x^4 + 1
      [243n, 0n, 0n, 0n, 0n, 0n, 1n], // x^6 + 243
      [1n, -3n, 0n, 1n], // x^3 - 3x + 1
      [1n, 1n, 1n, 1n, 1n, 1n, 1n], // Phi_7
    ];
    for (const f of fields) {
      const K = NumberFieldConstructor(f, 'a');
      const G = K.galois_group();
      const n = Number(G.order());
      for (const H of G.subgroups()) {
        const d = G.fixed_field_data(H);
        if (d.polynomial === null) continue; // Q or K
        expect(d.polynomial.length - 1).toBe(n / H.length);
        // P(gen) = 0 in K
        let acc = K.zero();
        let pw = K.one();
        for (const c of d.polynomial) {
          acc = acc.add(pw.scalarMul(new Rational(c)));
          pw = pw.mul(d.gen!);
        }
        expect(acc.is_zero()).toBe(true);
        // and gen is fixed by exactly the subgroup H
        for (const s of G.list()) {
          const fixed = s.__call__(d.gen!).eq(d.gen!);
          const inH = H.some((h) => h.eq(s));
          expect(fixed).toBe(inH);
        }
      }
    }
  }, 120000);

  it('galois_group of a non-Galois field names what is missing', () => {
    const K = NumberFieldConstructor([-2n, 0n, 0n, 1n], 'a'); // x^3 - 2
    expect(() => K.galois_group().order()).toThrow(NotImplementedError);
    expect(() => K.galois_group().order()).toThrow('nffactor');
  });
});

describe('complex_conjugation', () => {
  // sage: L.<z> = CyclotomicField(7)
  // sage: G = L.galois_group()
  // sage: conj = G.complex_conjugation(); conj
  // (1,4)(2,5)(3,6)
  // sage: conj(z)
  // -z^5 - z^4 - z^3 - z^2 - z - 1
  it('Q(zeta_7): conj(z) = -z^5 - z^4 - z^3 - z^2 - z - 1 at every place', () => {
    const L = NumberFieldConstructor([1n, 1n, 1n, 1n, 1n, 1n, 1n], 'z');
    const G = L.galois_group();
    const places = L.complex_embeddings();
    expect(places.length).toBe(6);
    for (const pl of places) {
      const c = G.complex_conjugation(pl);
      expect(c.__call__(L.gen()).toString()).toBe('-z^5 - z^4 - z^3 - z^2 - z - 1');
      expect(c.order()).toBe(2n);
    }
    // Sage raises when no default embedding is specified
    expect(() => G.complex_conjugation()).toThrow('No default complex embedding specified');
  }, 60000);

  // sage: L = NumberField(x^6 + 40*x^3 + 1372, 'a')
  // sage: G = L.galois_group()
  // sage: [G.complex_conjugation(x) for x in L.places()]
  // [(1,3)(2,6)(4,5), (1,5)(2,4)(3,6), (1,2)(3,4)(5,6)]
  it('x^6 + 40x^3 + 1372 is not CM: the three places give three different conjugations', () => {
    const M = NumberFieldConstructor([1372n, 0n, 0n, 40n, 0n, 0n, 1n], 'a');
    const G = M.galois_group();
    expect(M.signature()).toEqual([0, 3]);
    const conjs = M.places().map((pl) => G.complex_conjugation(pl));
    expect(conjs.length).toBe(3);
    for (const c of conjs) expect(c.order()).toBe(2n);
    expect(new Set(conjs.map((c) => c.toString())).size).toBe(3);
  }, 60000);

  it('refuses a totally real field, as SageMath does', () => {
    const K = NumberFieldConstructor([1n, -3n, 0n, 1n], 'a'); // x^3 - 3x + 1, totally real
    const G = K.galois_group();
    expect(K.is_totally_real()).toBe(true);
    const pl = K.complex_embeddings()[0]!;
    expect(() => G.complex_conjugation(pl)).toThrow('No complex conjugation');
  });
});
