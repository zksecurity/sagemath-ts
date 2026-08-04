/**
 * Tests for Weierstrass morphisms between elliptic curves
 */
import { describe, expect, it } from 'bun:test';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import type { FiniteFieldElement } from '../../rings/finite_rings/finite_field_prime.js';
import { EllipticCurve } from './constructor.js';
import {
  WeierstrassIsomorphism,
  _isomorphisms,
  baseWI,
  identity_morphism,
  negation_morphism,
} from './weierstrass_morphism.js';

describe('baseWI', () => {
  describe('constructor', () => {
    it('should create baseWI with identity parameters', () => {
      const F = GF(7n);
      const w = new baseWI(F.__call__(1n), F.__call__(0n), F.__call__(0n), F.__call__(0n));
      expect(w.is_identity()).toBe(true);
    });

    it('should create baseWI with custom parameters', () => {
      const F = GF(7n);
      const w = new baseWI(F.__call__(2n), F.__call__(3n), F.__call__(4n), F.__call__(5n));
      const [u, r, s, t] = w.tuple();
      expect(u.eq(F.__call__(2n))).toBe(true);
      expect(r.eq(F.__call__(3n))).toBe(true);
      expect(s.eq(F.__call__(4n))).toBe(true);
      expect(t.eq(F.__call__(5n))).toBe(true);
    });

    it('should throw error for u=0', () => {
      const F = GF(7n);
      expect(() => {
        new baseWI(F.__call__(0n), F.__call__(0n), F.__call__(0n), F.__call__(0n));
      }).toThrow('u!=0 required for baseWI');
    });
  });

  describe('composition (mul)', () => {
    it('should compose with identity', () => {
      const F = GF(7n);
      const w = new baseWI(F.__call__(2n), F.__call__(3n), F.__call__(4n), F.__call__(5n));
      const identity = new baseWI(F.__call__(1n), F.__call__(0n), F.__call__(0n), F.__call__(0n));
      const composed1 = w.mul(identity);
      const composed2 = identity.mul(w);
      expect(composed1.u.eq(w.u)).toBe(true);
      expect(composed2.u.eq(w.u)).toBe(true);
    });
  });

  describe('inverse (invert)', () => {
    it('should compute inverse correctly', () => {
      const F = GF(97n);
      const w = new baseWI(F.__call__(2n), F.__call__(3n), F.__call__(4n), F.__call__(5n));
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
    const hasIdentity = isos.some(
      ([u, r, s, t]) => u.eq(F.one()) && r.isZero() && s.isZero() && t.isZero()
    );
    expect(hasIdentity).toBe(true);
  });
});

describe('WeierstrassIsomorphism', () => {
  describe('constructor', () => {
    it('should create isomorphism from curve and (u,r,s,t)', () => {
      const F = GF(97n);
      const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
      const iso = new WeierstrassIsomorphism<FiniteFieldElement>(E, [
        F.__call__(2n) as FiniteFieldElement,
        F.__call__(0n) as FiniteFieldElement,
        F.__call__(0n) as FiniteFieldElement,
        F.__call__(0n) as FiniteFieldElement,
      ]);
      expect(iso.domain()).toBe(E);
      expect(iso.degree()).toBe(1n);
    });

    it('should throw for non-isomorphic curves', () => {
      const F = GF(97n);
      const E1 = EllipticCurve<FiniteFieldElement>(F, [1n, 0n]);
      const E2 = EllipticCurve<FiniteFieldElement>(F, [0n, 1n]);
      expect(() => {
        new WeierstrassIsomorphism<FiniteFieldElement>(E1, null, E2);
      }).toThrow('elliptic curves not isomorphic');
    });
  });

  describe('invert', () => {
    it('should return inverse isomorphism', () => {
      const F = GF(97n);
      const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
      const iso = new WeierstrassIsomorphism<FiniteFieldElement>(E, [
        F.__call__(2n) as FiniteFieldElement,
        F.__call__(3n) as FiniteFieldElement,
        F.__call__(4n) as FiniteFieldElement,
        F.__call__(5n) as FiniteFieldElement,
      ]);
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
      const iso = new WeierstrassIsomorphism<FiniteFieldElement>(E, [
        F.__call__(2n) as FiniteFieldElement,
        F.__call__(0n) as FiniteFieldElement,
        F.__call__(0n) as FiniteFieldElement,
        F.__call__(0n) as FiniteFieldElement,
      ]);
      expect(iso.degree()).toBe(1n);
      expect(iso.inseparable_degree()).toBe(1n);
    });

    it('should return scaling factor', () => {
      const F = GF(97n);
      const E = EllipticCurve<FiniteFieldElement>(F, [1n, 1n]);
      const iso = new WeierstrassIsomorphism<FiniteFieldElement>(E, [
        F.__call__(2n) as FiniteFieldElement,
        F.__call__(3n) as FiniteFieldElement,
        F.__call__(4n) as FiniteFieldElement,
        F.__call__(5n) as FiniteFieldElement,
      ]);
      expect(iso.scaling_factor().eq(F.__call__(2n))).toBe(true);
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
      const iso = new WeierstrassIsomorphism<FiniteFieldElement>(E, [
        F.__call__(2n) as FiniteFieldElement,
        F.__call__(3n) as FiniteFieldElement,
        F.__call__(4n) as FiniteFieldElement,
        F.__call__(5n) as FiniteFieldElement,
      ]);
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

// ---------------------------------------------------------------------------
// M109: WeierstrassIsomorphism(E, None, F) must pick the same representative
// as SageMath, which enumerates the roots of x^m - u^m in the order produced by
// Polynomial.roots() (i.e. sorted by the constant coefficient -r).
// ---------------------------------------------------------------------------

import { NotImplementedError, ValueError } from '../../errors.js';

describe('WeierstrassIsomorphism.order (M109)', () => {
  // weierstrass_morphism.py:WeierstrassIsomorphism.order doctest
  it('matches the GF(97) doctest', () => {
    const Fp = GF(97n);
    const E = EllipticCurve<FiniteFieldElement>(Fp, [1n, 28n]);
    const ws = new WeierstrassIsomorphism<FiniteFieldElement>(E, null, E);
    expect(ws.order()).toBe(2n);
  });

  it('enumerates automorphisms in SageMath order', () => {
    const Fp = GF(97n);
    const E = EllipticCurve<FiniteFieldElement>(Fp, [1n, 28n]);
    const auts = [..._isomorphisms(E, E)].map((t) => t.map(String).join(','));
    // (x^2 - 1).roots() over GF(97) is [96, 1]: x + 1 sorts before x + 96
    expect(auts).toEqual(['96,0,0,0', '1,0,0,0']);
  });

  it('raises ValueError when the domain differs from the codomain', () => {
    const Fp = GF(97n);
    const E = EllipticCurve<FiniteFieldElement>(Fp, [1n, 28n]);
    const E1 = EllipticCurve<FiniteFieldElement>(Fp, [1n, 69n]);
    const urst = E.isomorphism_to(E1);
    const ws = new WeierstrassIsomorphism<FiniteFieldElement>(E, urst, E1);
    expect(() => ws.order()).toThrow(ValueError);
  });

  it('matches the j = 0 doctests', () => {
    const Fp = GF(97n);
    const E = EllipticCurve<FiniteFieldElement>(Fp, [0n, 1n]); // j = 0
    const ws = new WeierstrassIsomorphism<FiniteFieldElement>(
      E,
      [Fp.__call__(36n), Fp.__call__(0n), Fp.__call__(0n), Fp.__call__(0n)],
      null
    );
    expect(ws.order()).toBe(6n);
    const ws2 = WeierstrassIsomorphism._composition_impl(ws, ws)!;
    expect(ws2.order()).toBe(3n);
  });

  it('raises NotImplementedError for orders outside {1,2,3,4,6}', () => {
    // Sage raises NotImplementedError (not ValueError) in this branch.  Genuine
    // automorphism groups only have order 1, 2, 3, 4 or 6, so the branch is
    // unreachable for real data; we force it by making is_identity() always
    // false on a valid automorphism.
    const Fp = GF(97n);
    const E = EllipticCurve<FiniteFieldElement>(Fp, [1n, 28n]);
    const ws = new WeierstrassIsomorphism<FiniteFieldElement>(E, null, E);
    const proto = WeierstrassIsomorphism.prototype as unknown as {
      is_identity: () => boolean;
    };
    const original = proto.is_identity;
    proto.is_identity = () => false;
    const thrown = (() => {
      let err: unknown = null;
      // eslint-disable-next-line no-empty
      try {
        ws.order();
      } catch (e) {
        err = e;
      }
      proto.is_identity = original;
      return err;
    })();
    expect(thrown).toBeInstanceOf(NotImplementedError);
  });
});
