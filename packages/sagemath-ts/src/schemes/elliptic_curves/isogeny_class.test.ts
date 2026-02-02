/**
 * Tests for isogeny classes of elliptic curves
 *
 * Reference: sage/schemes/elliptic_curves/isogeny_class.py
 *
 * NOTE: Many tests are limited because the full implementation requires
 * isogeny computation and database access which are not yet available.
 */

import { describe, expect, it, test } from 'bun:test';
import { NotImplementedError, ValueError } from '../../errors.js';
import type { EllipticCurveGeneric } from './ell_generic.js';
import {
  IsogenyClass,
  IsogenyClassNumberField,
  IsogenyClassRational,
  isogeny_degrees_cm,
  possible_isogeny_degrees,
} from './isogeny_class.js';
import type { FieldElement, FieldRing } from './types.js';

class MockFieldRing implements FieldRing {
  readonly characteristic = 0n;

  zero(): MockFieldElement {
    return new MockFieldElement(0n, this);
  }

  one(): MockFieldElement {
    return new MockFieldElement(1n, this);
  }

  __call__(value: bigint | number | FieldElement): MockFieldElement {
    if (value instanceof MockFieldElement) {
      return new MockFieldElement(value.value, this);
    }
    if (typeof value === 'bigint' || typeof value === 'number') {
      return new MockFieldElement(BigInt(value), this);
    }
    return new MockFieldElement(BigInt(value.toString()), this);
  }

  toString(): string {
    return 'MockField';
  }
}

class MockFieldElement implements FieldElement {
  readonly parent: FieldRing;

  constructor(
    readonly value: bigint,
    parent: FieldRing
  ) {
    this.parent = parent;
  }

  private toBigInt(other: FieldElement | number | bigint): bigint {
    if (typeof other === 'bigint' || typeof other === 'number') {
      return BigInt(other);
    }
    if (other instanceof MockFieldElement) {
      return other.value;
    }
    return BigInt(other.toString());
  }

  add(other: FieldElement | number | bigint): MockFieldElement {
    return new MockFieldElement(this.value + this.toBigInt(other), this.parent);
  }

  sub(other: FieldElement | number | bigint): MockFieldElement {
    return new MockFieldElement(this.value - this.toBigInt(other), this.parent);
  }

  mul(other: FieldElement | number | bigint): MockFieldElement {
    return new MockFieldElement(this.value * this.toBigInt(other), this.parent);
  }

  div(other: FieldElement | number | bigint): MockFieldElement {
    return new MockFieldElement(this.value / this.toBigInt(other), this.parent);
  }

  neg(): MockFieldElement {
    return new MockFieldElement(-this.value, this.parent);
  }

  inv(): MockFieldElement {
    return new MockFieldElement(1n, this.parent);
  }

  pow(_n: bigint | number): MockFieldElement {
    return new MockFieldElement(1n, this.parent);
  }

  isZero(): boolean {
    return this.value === 0n;
  }

  eq(other: FieldElement): boolean {
    if (other instanceof MockFieldElement) {
      return this.value === other.value;
    }
    return this.value === BigInt(other.toString());
  }

  toString(): string {
    return this.value.toString();
  }
}

type MockCurve = {
  _id: string;
  a_invariants: () => [
    MockFieldElement,
    MockFieldElement,
    MockFieldElement,
    MockFieldElement,
    MockFieldElement,
  ];
  is_isomorphic: (other: { _id?: string }) => boolean;
  j_invariant: () => MockFieldElement;
  toString: () => string;
  has_cm?: () => boolean;
  cm_discriminant?: () => bigint;
  has_rational_cm?: () => boolean;
  base_field?: () => { degree: () => number };
};

const mockField = new MockFieldRing();

const asCurve = (curve: MockCurve): EllipticCurveGeneric<FieldElement> =>
  curve as unknown as EllipticCurveGeneric<FieldElement>;

const createMockCurveRaw = (id: string, jInvariant: bigint = 0n): MockCurve => ({
  _id: id,
  a_invariants: () => [
    mockField.__call__(0n),
    mockField.__call__(0n),
    mockField.__call__(1n),
    mockField.__call__(-1n),
    mockField.__call__(0n),
  ],
  is_isomorphic: (other: { _id?: string }) => other._id === id,
  j_invariant: () => mockField.__call__(jInvariant),
  toString: () => `Elliptic Curve ${id}`,
});

// Mock elliptic curve for testing
const createMockCurve = (id: string, jInvariant: bigint = 0n): EllipticCurveGeneric<FieldElement> =>
  asCurve(createMockCurveRaw(id, jInvariant));

describe('IsogenyClass', () => {
  describe('constructor', () => {
    it('should create an isogeny class from a curve', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      expect(iso.E).toBe(E);
    });

    it('should store the label if provided', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E, '11a');
      expect(iso.toString()).toContain('11a');
    });

    it('should initialize with at least the original curve', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      expect(iso.length()).toBeGreaterThanOrEqual(1);
    });
  });

  describe('length', () => {
    it('should return the number of curves', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      expect(typeof iso.length()).toBe('number');
      expect(iso.length()).toBeGreaterThanOrEqual(1);
    });
  });

  describe('iterator', () => {
    it('should iterate over curves', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);

      const curves: Array<EllipticCurveGeneric<FieldElement>> = [];
      for (const curve of iso) {
        curves.push(curve);
      }

      expect(curves.length).toBe(iso.length());
    });

    it('should support spread operator', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const curves = [...iso];
      expect(curves.length).toBe(iso.length());
    });
  });

  describe('get', () => {
    it('should return the i-th curve', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const first = iso.get(0);
      expect(first).toBeDefined();
    });

    it('should throw for out of range index', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      expect(() => iso.get(-1)).toThrow(ValueError);
      expect(() => iso.get(100)).toThrow(ValueError);
    });
  });

  describe('index', () => {
    it('should return the index of a curve in the class', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const idx = iso.index(E);
      expect(idx).toBe(0);
    });

    it('should throw for curve not in class', () => {
      const E = createMockCurve('11a1');
      const E2 = createMockCurve('37a1');
      const iso = new IsogenyClass(E);
      expect(() => iso.index(E2)).toThrow(ValueError);
    });
  });

  describe('matrix', () => {
    it('should return a matrix of isogeny degrees', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const mat = iso.matrix();

      expect(Array.isArray(mat)).toBe(true);
      expect(mat.length).toBe(iso.length());
      expect(mat[0]!.length).toBe(iso.length());
    });

    it('should have 1s on the diagonal when filled', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const mat = iso.matrix(true);

      for (let i = 0; i < mat.length; i++) {
        expect(mat[i]![i]).toBe(1n);
      }
    });

    it('should have 0s on the diagonal when unfilled', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const mat = iso.matrix(false);

      for (let i = 0; i < mat.length; i++) {
        expect(mat[i]![i]).toBe(0n);
      }
    });
  });

  describe('qf_matrix', () => {
    it('should throw for non-CM curves', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      expect(() => iso.qf_matrix()).toThrow(ValueError);
    });
  });

  describe('isogenies', () => {
    it('should return a 2D array', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const maps = iso.isogenies();

      expect(Array.isArray(maps)).toBe(true);
      expect(maps.length).toBe(iso.length());
    });

    it('should throw for fill=true (not implemented)', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      expect(() => iso.isogenies(true)).toThrow(NotImplementedError);
    });
  });

  describe('graph', () => {
    it('should return a graph representation', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const graph = iso.graph();

      // Graph should have the expected structure
      expect(graph.vertices).toBeInstanceOf(Map);
      expect(graph.edges).toBeInstanceOf(Array);
      expect(graph.positions).toBeInstanceOf(Map);
      expect(graph.numVertices).toBe(iso.length());
      expect(typeof graph.numEdges).toBe('number');
    });

    it('should have vertices labeled 1 to n', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const graph = iso.graph();

      // Vertices should be 1-indexed (not 0-indexed)
      expect(graph.vertices.has(1)).toBe(true);
      expect(graph.vertices.has(0)).toBe(false);
    });

    it('should have positions for all vertices', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const graph = iso.graph();

      for (let i = 1; i <= graph.numVertices; i++) {
        expect(graph.positions.has(i)).toBe(true);
        const pos = graph.positions.get(i);
        expect(Array.isArray(pos)).toBe(true);
        expect(pos!.length).toBe(2);
      }
    });
  });

  describe('reorder', () => {
    it('should accept "lmfdb" ordering', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const reordered = iso.reorder('lmfdb');
      expect(reordered.length()).toBe(iso.length());
    });

    it('should accept array of indices', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const reordered = iso.reorder([0]);
      expect(reordered.length()).toBe(iso.length());
    });

    it('should throw for incorrect length', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      expect(() => iso.reorder([0, 1, 2, 3])).toThrow(ValueError);
    });
  });

  describe('copy', () => {
    it('should create an independent copy', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const copy = iso.copy();

      expect(copy).not.toBe(iso);
      expect(copy.length()).toBe(iso.length());
      expect(copy.E).toBe(iso.E);
    });

    it('should not share the curves array', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      const copy = iso.copy();

      // The curves array should be a different reference
      expect(copy.curves).not.toBe(iso.curves);
    });
  });

  describe('contains', () => {
    it('should return true for curves in the class', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      expect(iso.contains(E)).toBe(true);
    });

    it('should return false for curves not in the class', () => {
      const E = createMockCurve('11a1');
      const E2 = createMockCurve('37a1');
      const iso = new IsogenyClass(E);
      expect(iso.contains(E2)).toBe(false);
    });
  });

  describe('toString', () => {
    it('should include the label if present', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E, '11a');
      expect(iso.toString()).toContain('11a');
    });

    it('should include the curve if no label', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E);
      expect(iso.toString()).toContain('Isogeny class of');
    });
  });
});

describe('IsogenyClassNumberField', () => {
  it('should accept algorithm option', () => {
    const E = createMockCurve('11a1');
    const iso = new IsogenyClassNumberField(E, { algorithm: 'Billerey' });
    expect(iso.length()).toBeGreaterThanOrEqual(1);
  });

  it('should accept minimal_models option', () => {
    const E = createMockCurve('11a1');
    const iso = new IsogenyClassNumberField(E, { minimal_models: false });
    expect(iso.length()).toBeGreaterThanOrEqual(1);
  });

  describe('copy', () => {
    it('should return an IsogenyClassNumberField', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClassNumberField(E);
      const copy = iso.copy();
      expect(copy).toBeInstanceOf(IsogenyClassNumberField);
    });
  });
});

describe('IsogenyClassRational', () => {
  it('should accept algorithm parameter', () => {
    const E = createMockCurve('11a1');
    const iso = new IsogenyClassRational(E, 'sage');
    expect(iso.length()).toBeGreaterThanOrEqual(1);
  });

  it('should accept label parameter', () => {
    const E = createMockCurve('11a1');
    const iso = new IsogenyClassRational(E, 'sage', '11a');
    expect(iso.toString()).toContain('11a');
  });

  describe('copy', () => {
    it('should return an IsogenyClassRational', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClassRational(E);
      const copy = iso.copy();
      expect(copy).toBeInstanceOf(IsogenyClassRational);
    });
  });
});

describe('isogeny_degrees_cm', () => {
  it('should throw ValueError for non-CM curves', () => {
    const E = createMockCurve('27a1');
    // Mock curves don't have CM, so this should throw ValueError
    expect(() => isogeny_degrees_cm(E)).toThrow(ValueError);
  });

  it('should return primes for CM curves', () => {
    // Create a mock curve with CM
    const E = asCurve({
      ...createMockCurveRaw('CM'),
      has_cm: () => true,
      cm_discriminant: () => -4n,
      has_rational_cm: () => true,
      base_field: () => ({ degree: () => 1 }),
    });
    const result = isogeny_degrees_cm(E);
    expect(Array.isArray(result)).toBe(true);
    // Should include 2 at minimum
    expect(result.includes(2n)).toBe(true);
  });
});

describe('possible_isogeny_degrees', () => {
  it('should return Mazur primes for non-CM curves over Q', () => {
    const E = createMockCurve('11a1');
    const result = possible_isogeny_degrees(E);

    expect(Array.isArray(result)).toBe(true);
    // Should contain the Mazur primes (for curves over Q without CM)
    expect(result).toContain(2n);
    expect(result).toContain(3n);
    expect(result).toContain(5n);
    expect(result).toContain(7n);
  });

  it('should accept algorithm option', () => {
    const E = createMockCurve('11a1');

    // All algorithms should return results for curves over Q
    const billerey = possible_isogeny_degrees(E, { algorithm: 'Billerey' });
    const larson = possible_isogeny_degrees(E, { algorithm: 'Larson' });
    const heuristic = possible_isogeny_degrees(E, { algorithm: 'heuristic' });

    expect(Array.isArray(billerey)).toBe(true);
    expect(Array.isArray(larson)).toBe(true);
    expect(Array.isArray(heuristic)).toBe(true);
  });

  it('should respect max_l option', () => {
    const E = createMockCurve('11a1');
    const result = possible_isogeny_degrees(E, { max_l: 10 });

    // Should only contain primes <= 10
    for (const p of result) {
      expect(p <= 10n).toBe(true);
    }
  });

  it('should delegate to isogeny_degrees_cm for CM curves', () => {
    const E = asCurve({
      ...createMockCurveRaw('CM'),
      has_cm: () => true,
      cm_discriminant: () => -3n,
      has_rational_cm: () => true,
      base_field: () => ({ degree: () => 1 }),
    });

    const result = possible_isogeny_degrees(E);
    expect(Array.isArray(result)).toBe(true);
  });
});
