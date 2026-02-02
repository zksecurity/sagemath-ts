/**
 * Tests for Galois group functionality
 */

import { describe, expect, it } from 'bun:test';
import {
  GaloisGroup,
  type GaloisGroupElement,
  GaloisSubgroup,
  Permutation,
  galois_group,
} from './galois_group.js';
import {
  CyclotomicField,
  NumberField,
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
  it('should create Galois group for cubic field', () => {
    // x^3 - 2 (not Galois over Q, Galois closure has degree 6)
    const poly = RationalPolynomial.fromBigInts([-2n, 0n, 0n, 1n]);
    const K = new NumberField(poly, 'a');
    const G = galois_group(K);
    expect(G.degree()).toBe(3);
  });

  it('should create Galois group for cyclic cubic', () => {
    // x^3 - 3x + 1 is a cyclic cubic (totally real, Galois)
    const poly = RationalPolynomial.fromBigInts([1n, -3n, 0n, 1n]);
    const K = new NumberField(poly, 'a');
    const G = galois_group(K);
    expect(G.degree()).toBe(3);
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
