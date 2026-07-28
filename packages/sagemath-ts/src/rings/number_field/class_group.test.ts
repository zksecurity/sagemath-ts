/**
 * Tests for ClassGroup functionality
 */

import { describe, expect, it } from 'bun:test';
import {
  ClassGroup,
  ClassGroupElement,
  NarrowClassGroup,
  RayClassGroup,
  quadraticClassGroupInvariants,
  quadraticClassNumber,
} from './class_group.js';
import { NumberField, QuadraticField, RationalPolynomial } from './number_field.js';

describe('ClassGroup', () => {
  describe('construction', () => {
    it('should create trivial class group', () => {
      const K = QuadraticField.create(-1n); // Q(i) has class number 1
      const Cl = new ClassGroup(K, [], []);
      expect(Cl.order()).toBe(1n);
      expect(Cl.is_trivial()).toBe(true);
    });

    it('should create class group with invariants', () => {
      const K = QuadraticField.create(-23n); // Q(sqrt(-23)) has class number 3
      const Cl = new ClassGroup(K, [3n], []);
      expect(Cl.order()).toBe(3n);
      expect(Cl.is_cyclic()).toBe(true);
    });

    it('should create class group with multiple invariants', () => {
      // Simulating a class group with structure Z/2 x Z/4
      const K = QuadraticField.create(-5n);
      const Cl = new ClassGroup(K, [2n, 4n], []);
      expect(Cl.order()).toBe(8n);
      expect(Cl.is_cyclic()).toBe(false);
    });
  });

  describe('class_number', () => {
    it('should return 1 for trivial class group', () => {
      const K = QuadraticField.create(-1n);
      const Cl = new ClassGroup(K, [], []);
      expect(Cl.class_number()).toBe(1n);
    });

    it('should return correct class number', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      expect(Cl.class_number()).toBe(3n);
    });
  });

  describe('class_group_generators', () => {
    it('should return empty for trivial group', () => {
      const K = QuadraticField.create(-1n);
      const Cl = new ClassGroup(K, [], []);
      expect(Cl.class_group_generators().length).toBe(0);
    });
  });

  describe('exponent', () => {
    it('should return 1 for trivial group', () => {
      const K = QuadraticField.create(-1n);
      const Cl = new ClassGroup(K, [], []);
      expect(Cl.exponent()).toBe(1n);
    });

    it('should return largest invariant', () => {
      const K = QuadraticField.create(-5n);
      const Cl = new ClassGroup(K, [2n, 6n], []);
      expect(Cl.exponent()).toBe(6n);
    });
  });

  describe('invariants', () => {
    it('should return copy of invariants', () => {
      const K = QuadraticField.create(-5n);
      const Cl = new ClassGroup(K, [2n, 3n], []);
      const invs = Cl.invariants();
      expect(invs).toEqual([2n, 3n]);
    });
  });

  describe('identity', () => {
    it('should return identity element', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const id = Cl.identity();
      expect(id.is_trivial()).toBe(true);
      expect(id.is_principal()).toBe(true);
    });
  });

  describe('random_element', () => {
    it('should return valid element', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const elem = Cl.random_element();
      expect(elem.parent()).toBe(Cl);
    });

    it('should return identity for trivial group', () => {
      const K = QuadraticField.create(-1n);
      const Cl = new ClassGroup(K, [], []);
      const elem = Cl.random_element();
      expect(elem.is_trivial()).toBe(true);
    });
  });

  describe('is_cyclic', () => {
    it('should return true for cyclic groups', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      expect(Cl.is_cyclic()).toBe(true);
    });

    it('should return false for non-cyclic groups', () => {
      const K = QuadraticField.create(-5n);
      const Cl = new ClassGroup(K, [2n, 2n], []);
      expect(Cl.is_cyclic()).toBe(false);
    });
  });

  describe('list', () => {
    it('should list all elements of trivial group', () => {
      const K = QuadraticField.create(-1n);
      const Cl = new ClassGroup(K, [], []);
      const elements = Cl.list();
      expect(elements.length).toBe(1);
    });

    it('should list all elements of cyclic group', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const elements = Cl.list();
      expect(elements.length).toBe(3);
    });
  });

  describe('iterator', () => {
    it('should iterate over all elements', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const elements: ClassGroupElement[] = [];
      for (const elem of Cl) {
        elements.push(elem);
      }
      expect(elements.length).toBe(3);
    });
  });
});

describe('ClassGroupElement', () => {
  describe('construction', () => {
    it('should create identity element', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const elem = new ClassGroupElement(Cl, [0n]);
      expect(elem.is_trivial()).toBe(true);
    });

    it('should create non-trivial element', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const elem = new ClassGroupElement(Cl, [1n]);
      expect(elem.is_trivial()).toBe(false);
    });

    it('should normalize exponents', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const elem = new ClassGroupElement(Cl, [4n]); // 4 mod 3 = 1
      expect(elem.exponents()).toEqual([1n]);
    });
  });

  describe('order', () => {
    it('should return 1 for identity', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const id = Cl.identity();
      expect(id.order()).toBe(1n);
    });

    it('should compute order of non-trivial element', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const elem = new ClassGroupElement(Cl, [1n]);
      expect(elem.order()).toBe(3n);
    });
  });

  describe('multiplication', () => {
    it('should multiply elements', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const a = new ClassGroupElement(Cl, [1n]);
      const b = new ClassGroupElement(Cl, [1n]);
      const product = a.mul(b);
      expect(product.exponents()).toEqual([2n]);
    });

    it('should satisfy group axioms', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const a = new ClassGroupElement(Cl, [1n]);
      const id = Cl.identity();

      // a * id = a
      expect(a.mul(id).eq(a)).toBe(true);

      // id * a = a
      expect(id.mul(a).eq(a)).toBe(true);
    });
  });

  describe('inverse', () => {
    it('should compute inverse', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const a = new ClassGroupElement(Cl, [1n]);
      const aInv = a.inv();

      // a * a^(-1) = identity
      expect(a.mul(aInv).is_trivial()).toBe(true);
    });
  });

  describe('power', () => {
    it('should compute powers', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const a = new ClassGroupElement(Cl, [1n]);

      // a^0 = identity
      expect(a.pow(0n).is_trivial()).toBe(true);

      // a^1 = a
      expect(a.pow(1n).eq(a)).toBe(true);

      // a^3 = identity (since order is 3)
      expect(a.pow(3n).is_trivial()).toBe(true);
    });
  });

  describe('equality', () => {
    it('should check equality', () => {
      const K = QuadraticField.create(-23n);
      const Cl = new ClassGroup(K, [3n], []);
      const a = new ClassGroupElement(Cl, [1n]);
      const b = new ClassGroupElement(Cl, [1n]);
      const c = new ClassGroupElement(Cl, [2n]);

      expect(a.eq(b)).toBe(true);
      expect(a.eq(c)).toBe(false);
    });
  });
});

describe('NarrowClassGroup', () => {
  it('should extend ClassGroup', () => {
    const K = QuadraticField.create(-5n);
    const NCl = new NarrowClassGroup(K, [2n], []);
    expect(NCl.order()).toBe(2n);
  });
});

describe('RayClassGroup', () => {
  it('should store modulus', () => {
    const K = QuadraticField.create(-5n);
    const modulus = { prime: 2n };
    const RCl = new RayClassGroup(K, modulus, [4n], []);
    expect(RCl.modulus()).toEqual(modulus);
    expect(RCl.order()).toBe(4n);
  });
});

describe('Integration tests', () => {
  it('should work with number field class_group method', () => {
    // Q(i) has class number 1
    const K = QuadraticField.create(-1n);
    const Cl = K.class_group();
    expect(Cl.order()).toBe(1n);
    expect(Cl.is_trivial()).toBe(true);
  });

  it('should work with number field class_number method', () => {
    // Q(sqrt(-23)) has class number 3
    const K = QuadraticField.create(-23n);
    expect(K.class_number()).toBe(3n);
  });
});

// ---------------------------------------------------------------------------
// Quadratic class groups computed from binary quadratic forms (audit H17/H18)
// ---------------------------------------------------------------------------

describe('quadraticClassNumber', () => {
  it('matches PARI qfbclassno on imaginary quadratic discriminants', () => {
    const known: Array<[bigint, bigint]> = [
      [-3n, 1n],
      [-4n, 1n],
      [-7n, 1n],
      [-8n, 1n],
      [-11n, 1n],
      [-15n, 2n],
      [-20n, 2n],
      [-23n, 3n],
      [-24n, 2n],
      [-31n, 3n],
      [-39n, 4n],
      [-47n, 5n],
      [-71n, 7n],
      [-95n, 8n],
      [-119n, 10n],
      [-120n, 4n],
      [-163n, 1n],
      [-20072n, 76n],
    ];
    for (const [D, h] of known) {
      expect(quadraticClassNumber(D)).toBe(h);
    }
  });

  it('matches PARI quadclassunit on real quadratic discriminants', () => {
    const known: Array<[bigint, bigint]> = [
      [5n, 1n],
      [8n, 1n],
      [12n, 1n],
      [40n, 2n],
      [60n, 2n],
      [65n, 2n],
      [85n, 2n],
      [104n, 2n],
      [105n, 2n],
      [120n, 2n],
      [229n, 3n],
      [257n, 3n],
      [473n, 3n],
    ];
    for (const [D, h] of known) {
      expect(quadraticClassNumber(D)).toBe(h);
    }
  });

  it('also handles non-fundamental (order) discriminants', () => {
    // h(-44) = 3, h(-36) = 2, h(-108) = 3 (class numbers of the orders)
    expect(quadraticClassNumber(-44n)).toBe(3n);
    expect(quadraticClassNumber(-36n)).toBe(2n);
    expect(quadraticClassNumber(-108n)).toBe(3n);
  });
});

describe('quadraticClassGroupInvariants', () => {
  it('returns the elementary divisors in Sage/PARI order', () => {
    expect(quadraticClassGroupInvariants(-120n)).toEqual([2n, 2n]);
    expect(quadraticClassGroupInvariants(-84n)).toEqual([2n, 2n]);
    expect(quadraticClassGroupInvariants(-20072n)).toEqual([38n, 2n]);
    expect(quadraticClassGroupInvariants(-23n)).toEqual([3n]);
    expect(quadraticClassGroupInvariants(-47n)).toEqual([5n]);
    expect(quadraticClassGroupInvariants(-39n)).toEqual([4n]);
    expect(quadraticClassGroupInvariants(-4n)).toEqual([]);
    expect(quadraticClassGroupInvariants(40n)).toEqual([2n]);
    expect(quadraticClassGroupInvariants(8n)).toEqual([]);
  });

  it('has product equal to the class number', () => {
    for (const D of [-120n, -84n, -39n, -47n, -20072n, 40n, 229n]) {
      const invs = quadraticClassGroupInvariants(D);
      expect(invs.reduce((a, b) => a * b, 1n)).toBe(quadraticClassNumber(D));
      // each divides the previous
      for (let i = 1; i < invs.length; i++) {
        expect(invs[i - 1]! % invs[i]!).toBe(0n);
      }
    }
  });
});
