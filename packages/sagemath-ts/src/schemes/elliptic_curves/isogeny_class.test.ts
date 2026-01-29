/**
 * Tests for isogeny classes of elliptic curves
 *
 * Reference: sage/schemes/elliptic_curves/isogeny_class.py
 *
 * NOTE: Many tests are limited because the full implementation requires
 * isogeny computation and database access which are not yet available.
 */

import { describe, expect, it, test } from 'bun:test';
import {
  IsogenyClass,
  IsogenyClassNumberField,
  IsogenyClassRational,
  isogeny_degrees_cm,
  possible_isogeny_degrees,
} from './isogeny_class.js';
import { ValueError, NotImplementedError } from '../../errors.js';

// Mock field element
const mockFieldElement = (value: bigint) => ({
  isZero: () => value === 0n,
  eq: (other: any) => value === other.value,
  toString: () => value.toString(),
  value,
  add: (other: any) => mockFieldElement(value + (typeof other === 'bigint' ? other : other.value)),
  sub: (other: any) => mockFieldElement(value - (typeof other === 'bigint' ? other : other.value)),
  mul: (other: any) => mockFieldElement(value * (typeof other === 'bigint' ? other : other.value)),
  div: (other: any) => mockFieldElement(value / (typeof other === 'bigint' ? other : other.value)),
  neg: () => mockFieldElement(-value),
  inv: () => mockFieldElement(1n), // Simplified
  pow: (_n: bigint) => mockFieldElement(1n), // Simplified
  parent: {} as any,
});

// Mock elliptic curve for testing
const createMockCurve = (id: string, jInvariant: bigint = 0n) => ({
  _id: id,
  a_invariants: () => [
    mockFieldElement(0n),
    mockFieldElement(0n),
    mockFieldElement(1n),
    mockFieldElement(-1n),
    mockFieldElement(0n),
  ],
  is_isomorphic: (other: any) => other._id === id,
  j_invariant: () => mockFieldElement(jInvariant),
  toString: () => `Elliptic Curve ${id}`,
});

describe('IsogenyClass', () => {
  describe('constructor', () => {
    it('should create an isogeny class from a curve', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      expect(iso.E).toBe(E);
    });

    it('should store the label if provided', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any, '11a');
      expect(iso.toString()).toContain('11a');
    });

    it('should initialize with at least the original curve', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      expect(iso.length()).toBeGreaterThanOrEqual(1);
    });
  });

  describe('length', () => {
    it('should return the number of curves', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      expect(typeof iso.length()).toBe('number');
      expect(iso.length()).toBeGreaterThanOrEqual(1);
    });
  });

  describe('iterator', () => {
    it('should iterate over curves', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);

      const curves: any[] = [];
      for (const curve of iso) {
        curves.push(curve);
      }

      expect(curves.length).toBe(iso.length());
    });

    it('should support spread operator', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      const curves = [...iso];
      expect(curves.length).toBe(iso.length());
    });
  });

  describe('get', () => {
    it('should return the i-th curve', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      const first = iso.get(0);
      expect(first).toBeDefined();
    });

    it('should throw for out of range index', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      expect(() => iso.get(-1)).toThrow(ValueError);
      expect(() => iso.get(100)).toThrow(ValueError);
    });
  });

  describe('index', () => {
    it('should return the index of a curve in the class', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      const idx = iso.index(E as any);
      expect(idx).toBe(0);
    });

    it('should throw for curve not in class', () => {
      const E = createMockCurve('11a1');
      const E2 = createMockCurve('37a1');
      const iso = new IsogenyClass(E as any);
      expect(() => iso.index(E2 as any)).toThrow(ValueError);
    });
  });

  describe('matrix', () => {
    it('should return a matrix of isogeny degrees', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      const mat = iso.matrix();

      expect(Array.isArray(mat)).toBe(true);
      expect(mat.length).toBe(iso.length());
      expect(mat[0]!.length).toBe(iso.length());
    });

    it('should have 1s on the diagonal when filled', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      const mat = iso.matrix(true);

      for (let i = 0; i < mat.length; i++) {
        expect(mat[i]![i]).toBe(1n);
      }
    });

    it('should have 0s on the diagonal when unfilled', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      const mat = iso.matrix(false);

      for (let i = 0; i < mat.length; i++) {
        expect(mat[i]![i]).toBe(0n);
      }
    });
  });

  describe('qf_matrix', () => {
    it('should throw for non-CM curves', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      expect(() => iso.qf_matrix()).toThrow(ValueError);
    });
  });

  describe('isogenies', () => {
    it('should return a 2D array', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      const maps = iso.isogenies();

      expect(Array.isArray(maps)).toBe(true);
      expect(maps.length).toBe(iso.length());
    });

    it('should throw for fill=true (not implemented)', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      expect(() => iso.isogenies(true)).toThrow(NotImplementedError);
    });
  });

  describe('graph', () => {
    it('should return a graph representation', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
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
      const iso = new IsogenyClass(E as any);
      const graph = iso.graph();

      // Vertices should be 1-indexed (not 0-indexed)
      expect(graph.vertices.has(1)).toBe(true);
      expect(graph.vertices.has(0)).toBe(false);
    });

    it('should have positions for all vertices', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
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
      const iso = new IsogenyClass(E as any);
      const reordered = iso.reorder('lmfdb');
      expect(reordered.length()).toBe(iso.length());
    });

    it('should accept array of indices', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      const reordered = iso.reorder([0]);
      expect(reordered.length()).toBe(iso.length());
    });

    it('should throw for incorrect length', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      expect(() => iso.reorder([0, 1, 2, 3])).toThrow(ValueError);
    });
  });

  describe('copy', () => {
    it('should create an independent copy', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      const copy = iso.copy();

      expect(copy).not.toBe(iso);
      expect(copy.length()).toBe(iso.length());
      expect(copy.E).toBe(iso.E);
    });

    it('should not share the curves array', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      const copy = iso.copy();

      // The curves array should be a different reference
      expect(copy.curves).not.toBe(iso.curves);
    });
  });

  describe('contains', () => {
    it('should return true for curves in the class', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      expect(iso.contains(E as any)).toBe(true);
    });

    it('should return false for curves not in the class', () => {
      const E = createMockCurve('11a1');
      const E2 = createMockCurve('37a1');
      const iso = new IsogenyClass(E as any);
      expect(iso.contains(E2 as any)).toBe(false);
    });
  });

  describe('toString', () => {
    it('should include the label if present', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any, '11a');
      expect(iso.toString()).toContain('11a');
    });

    it('should include the curve if no label', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClass(E as any);
      expect(iso.toString()).toContain('Isogeny class of');
    });
  });
});

describe('IsogenyClassNumberField', () => {
  it('should accept algorithm option', () => {
    const E = createMockCurve('11a1');
    const iso = new IsogenyClassNumberField(E as any, { algorithm: 'Billerey' });
    expect(iso.length()).toBeGreaterThanOrEqual(1);
  });

  it('should accept minimal_models option', () => {
    const E = createMockCurve('11a1');
    const iso = new IsogenyClassNumberField(E as any, { minimal_models: false });
    expect(iso.length()).toBeGreaterThanOrEqual(1);
  });

  describe('copy', () => {
    it('should return an IsogenyClassNumberField', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClassNumberField(E as any);
      const copy = iso.copy();
      expect(copy).toBeInstanceOf(IsogenyClassNumberField);
    });
  });
});

describe('IsogenyClassRational', () => {
  it('should accept algorithm parameter', () => {
    const E = createMockCurve('11a1');
    const iso = new IsogenyClassRational(E as any, 'sage');
    expect(iso.length()).toBeGreaterThanOrEqual(1);
  });

  it('should accept label parameter', () => {
    const E = createMockCurve('11a1');
    const iso = new IsogenyClassRational(E as any, 'sage', '11a');
    expect(iso.toString()).toContain('11a');
  });

  describe('copy', () => {
    it('should return an IsogenyClassRational', () => {
      const E = createMockCurve('11a1');
      const iso = new IsogenyClassRational(E as any);
      const copy = iso.copy();
      expect(copy).toBeInstanceOf(IsogenyClassRational);
    });
  });
});

describe('isogeny_degrees_cm', () => {
  it('should throw ValueError for non-CM curves', () => {
    const E = createMockCurve('27a1');
    // Mock curves don't have CM, so this should throw ValueError
    expect(() => isogeny_degrees_cm(E as any)).toThrow(ValueError);
  });

  it('should return primes for CM curves', () => {
    // Create a mock curve with CM
    const E = {
      ...createMockCurve('CM'),
      has_cm: () => true,
      cm_discriminant: () => -4n,
      has_rational_cm: () => true,
      base_field: () => ({ degree: () => 1 }),
    };
    const result = isogeny_degrees_cm(E as any);
    expect(Array.isArray(result)).toBe(true);
    // Should include 2 at minimum
    expect(result.includes(2n)).toBe(true);
  });
});

describe('possible_isogeny_degrees', () => {
  it('should return Mazur primes for non-CM curves over Q', () => {
    const E = createMockCurve('11a1');
    const result = possible_isogeny_degrees(E as any);

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
    const billerey = possible_isogeny_degrees(E as any, { algorithm: 'Billerey' });
    const larson = possible_isogeny_degrees(E as any, { algorithm: 'Larson' });
    const heuristic = possible_isogeny_degrees(E as any, { algorithm: 'heuristic' });

    expect(Array.isArray(billerey)).toBe(true);
    expect(Array.isArray(larson)).toBe(true);
    expect(Array.isArray(heuristic)).toBe(true);
  });

  it('should respect max_l option', () => {
    const E = createMockCurve('11a1');
    const result = possible_isogeny_degrees(E as any, { max_l: 10 });

    // Should only contain primes <= 10
    for (const p of result) {
      expect(p <= 10n).toBe(true);
    }
  });

  it('should delegate to isogeny_degrees_cm for CM curves', () => {
    const E = {
      ...createMockCurve('CM'),
      has_cm: () => true,
      cm_discriminant: () => -3n,
      has_rational_cm: () => true,
      base_field: () => ({ degree: () => 1 }),
    };

    const result = possible_isogeny_degrees(E as any);
    expect(Array.isArray(result)).toBe(true);
  });
});
