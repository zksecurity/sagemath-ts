/**
 * Tests for Weierstrass morphisms between elliptic curves
 */
import { describe, expect, it } from 'vitest';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import { EllipticCurve } from './constructor.js';
import type { FiniteFieldElement } from '../../rings/finite_rings/finite_field_prime.js';
import { baseWI, WeierstrassIsomorphism, _isomorphisms, identity_morphism, negation_morphism } from './weierstrass_morphism.js';

describe('baseWI', () => {
  describe('constructor', () => {
    it('should create baseWI with identity parameters', () => {
      const F = GF(7n);
      const w = new baseWI(F.__call__(1), F.__call__(0), F.__call__(0), F.__call__(0));
      expect(w.is_identity()).toBe(true);
    });

    it('should create baseWI with custom parameters', () => {
      const F = GF(7n);
      const w = new baseWI(F.__call__(2), F.__call__(3), F.__call__(4), F.__call__(5));
      const [u, r, s, t] = w.tuple();
      expect(u.eq(F.__call__(2))).toBe(true);
      expect(r.eq(F.__call__(3))).toBe(true);
      expect(s.eq(F.__call__(4))).toBe(true);
      expect(t.eq(F.__call__(5))).toBe(true);
    });

    it('should throw error for u=0', () => {
      const F = GF(7n);
      expect(() => { new baseWI(F.__call__(0), F.__call__(0), F.__call__(0), F.__call__(0)); }).toThrow('u!=0 required for baseWI');
    });
  });

  describe('composition (mul)', () => {
    it('should compose with identity', () => {
      const F = GF(7n);
      const w = new baseWI(F.__call__(2), F.__call__(3), F.__call__(4), F.__call__(5));
      const identity = new baseWI(F.__call__(1), F.__call__(0), F.__call__(0), F.__call__(0));
      const composed1 = w.mul(identity);
      const composed2 = identity.mul(w);
      expect(composed1.u.eq(w.u)).toBe(true);
      expect(composed2.u.eq(w.u)).toBe(true);
    });
  });

  describe('inverse (invert)', () => {
    it('should compute inverse correctly', () => {
      const F = GF(97n);
      const w = new baseWI(F.__call__(2), F.__call__(3), F.__call__(4), F.__call__(5));
      const winv = w.invert();
      const composed = w.mul(winv);
      expect(composed.is_identity()).toBe(true);
    });
  });
});

describe('_isomorphisms', () => {
  it('should return no isomorphisms for different j-invariants', () => {
    const F = GF(23n);
    const E1 = EllipticCurve<FiniteFieldElement>(F, [1n, 0n]);
    const E2 = EllipticCurve<FiniteFieldElement>(F, [0n, 1n]);
    const isos = [..._isomorphisms(E1, E2)];
    expect(isos.length).toBe(0);
  });

  it('should find automorphisms of a curve', () => {
    const F = GF(23n);
    const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
    const isos = [..._isomorphisms(E, E)];
    expect(isos.length).toBeGreaterThan(0);
    const hasIdentity = isos.some(([u, r, s, t]) => u.eq(F.one()) && r.isZero() && s.isZero() && t.isZero());
    expect(hasIdentity).toBe(true);
  });
});

describe('WeierstrassIsomorphism', () => {
  describe('constructor', () => {
    it('should create isomorphism from curve and (u,r,s,t)', () => {
      const F = GF(97n);
      const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
      const iso = new WeierstrassIsomorphism<FiniteFieldElement>(E, [F.__call__(2) as FiniteFieldElement, F.__call__(0) as FiniteFieldElement, F.__call__(0) as FiniteFieldElement, F.__call__(0) as FiniteFieldElement]);
      expect(iso.domain()).toBe(E);
      expect(iso.degree()).toBe(1n);
    });

    it('should throw for non-isomorphic curves', () => {
      const F = GF(97n);
      const E1 = EllipticCurve<FiniteFieldElement>(F, [1n, 0n]);
      const E2 = EllipticCurve<FiniteFieldElement>(F, [0n, 1n]);
      expect(() => { new WeierstrassIsomorphism<FiniteFieldElement>(E1, null, E2); }).toThrow('elliptic curves not isomorphic');
    });
  });

  describe('invert', () => {
    it('should return inverse isomorphism', () => {
      const F = GF(97n);
      const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
      const iso = new WeierstrassIsomorphism<FiniteFieldElement>(E, [F.__call__(2) as FiniteFieldElement, F.__call__(3) as FiniteFieldElement, F.__call__(4) as FiniteFieldElement, F.__call__(5) as FiniteFieldElement]);
      const isoInv = iso.invert();
      expect(isoInv.domain()).toBe(iso.codomain());
      expect(isoInv.codomain()).toBe(iso.domain());
      const composed = WeierstrassIsomorphism._composition_impl(isoInv, iso);
      expect(composed).not.toBeNull();
      expect(composed!.is_identity()).toBe(true);
    });
  });

  describe('properties', () => {
    it('should have degree 1', () => {
      const F = GF(97n);
      const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
      const iso = new WeierstrassIsomorphism<FiniteFieldElement>(E, [F.__call__(2) as FiniteFieldElement, F.__call__(0) as FiniteFieldElement, F.__call__(0) as FiniteFieldElement, F.__call__(0) as FiniteFieldElement]);
      expect(iso.degree()).toBe(1n);
      expect(iso.inseparable_degree()).toBe(1n);
    });

    it('should return scaling factor', () => {
      const F = GF(97n);
      const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
      const iso = new WeierstrassIsomorphism<FiniteFieldElement>(E, [F.__call__(2) as FiniteFieldElement, F.__call__(3) as FiniteFieldElement, F.__call__(4) as FiniteFieldElement, F.__call__(5) as FiniteFieldElement]);
      expect(iso.scaling_factor().eq(F.__call__(2))).toBe(true);
    });
  });

  describe('order', () => {
    it('should return 1 for identity', () => {
      const F = GF(97n);
      const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
      const id = identity_morphism(E);
      expect(id.order()).toBe(1n);
    });

    it('should return 2 for negation', () => {
      const F = GF(97n);
      const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
      const neg = negation_morphism(E);
      expect(neg.order()).toBe(2n);
    });
  });

  describe('point transformation', () => {
    it('should map infinity to infinity', () => {
      const F = GF(97n);
      const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
      const iso = new WeierstrassIsomorphism<FiniteFieldElement>(E, [F.__call__(2) as FiniteFieldElement, F.__call__(3) as FiniteFieldElement, F.__call__(4) as FiniteFieldElement, F.__call__(5) as FiniteFieldElement]);
      const O = E.zero();
      const imageO = iso._call_(O);
      expect(imageO.is_zero()).toBe(true);
    });
  });
});

describe('identity_morphism', () => {
  it('should create identity isomorphism', () => {
    const F = GF(97n);
    const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
    const id = identity_morphism(E);
    expect(id.is_identity()).toBe(true);
    expect(id.domain()).toBe(E);
    // Codomain may be a new curve object with same invariants
    expect(id.codomain().j_invariant().eq(E.j_invariant())).toBe(true);
  });
});

describe('negation_morphism', () => {
  it('should create negation isomorphism', () => {
    const F = GF(97n);
    const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
    const neg = negation_morphism(E);
    expect(neg.is_identity()).toBe(false);
    expect(neg.domain()).toBe(E);
    // Codomain may be a new curve object with same invariants
    expect(neg.codomain().j_invariant().eq(E.j_invariant())).toBe(true);
  });
});

describe('automorphisms', () => {
  it('should find 2 automorphisms for generic j', () => {
    const F = GF(23n);
    const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
    const auts = [..._isomorphisms(E, E)];
    expect(auts.length).toBe(2);
  });
});
