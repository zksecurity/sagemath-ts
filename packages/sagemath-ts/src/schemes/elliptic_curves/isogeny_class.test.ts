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
  Frobenius_filter,
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

/** Duck-typed curve stand-in for the functions that only read a few methods. */
const asAnyCurve = (curve: Record<string, unknown>): EllipticCurveGeneric<FieldElement> =>
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
    // Curve 32a1, y^2 = x^3 + 4x: it really does have CM by the order of
    // discriminant -4, and it really does admit a rational 2-isogeny.
    //
    // This test used to use `createMockCurveRaw('CM')`, whose a-invariants are
    // those of 37a1 -- a curve with no CM and no isogenies at all -- while
    // claiming `cm_discriminant() == -4`.  With Frobenius_filter now applied
    // (as SageMath does, isogeny_class.py:1326) that inconsistent input
    // correctly loses its 2: 37a1's 2-division polynomial 4*x^3 - 4*x + 1 is
    // irreducible, so 2 is discarded by gal_reps_number_field.py:556-559.
    const E = asAnyCurve({
      ainvs: () => [0n, 0n, 0n, 4n, 0n],
      a_invariants: () => [0n, 0n, 0n, 4n, 0n],
      has_cm: () => true,
      cm_discriminant: () => -4n,
      has_rational_cm: () => true,
      base_field: () => ({ degree: () => 1 }),
      is_isomorphic: () => false,
      j_invariant: () => mockField.__call__(1728n),
      toString: () => 'Elliptic Curve 32a1',
    });
    const result = isogeny_degrees_cm(E);
    expect(Array.isArray(result)).toBe(true);
    // Should include 2 at minimum
    expect(result.includes(2n)).toBe(true);
  });

  // isogeny_class.py:1193-1203 -- SageMath's own doctest, printed verbatim:
  //   CM case, discriminant = -23
  //   initial primes: {2}
  //   upward primes: {}
  //   downward ramified primes: {}
  //   downward split primes: {2, 3}
  //   downward inert primes: {5}
  //   primes generating the class group: [2]
  //   Set of primes before filtering: {2, 3, 5}
  //   List of primes after filtering: [2, 3]
  //   [2, 3]
  it("reproduces SageMath's candidate set for d = -23", () => {
    const lines: string[] = [];
    const log = console.log;
    console.log = (...args: unknown[]) => {
      lines.push(args.join(' '));
    };
    let result: bigint[];
    try {
      result = isogeny_degrees_cm(
        asAnyCurve({
          has_cm: () => true,
          cm_discriminant: () => -23n,
          has_rational_cm: () => true,
          base_field: () => ({ degree: () => 6 }),
        }),
        true
      );
    } finally {
      console.log = log;
    }
    expect(lines).toEqual([
      'CM case, discriminant = -23',
      'initial primes: {2}',
      'upward primes: {}',
      'downward ramified primes: {}',
      'downward split primes: {2, 3}',
      'downward inert primes: {5}',
      'primes generating the class group: [2]',
      'Set of primes before filtering: {2, 3, 5}',
      'List of primes: [2, 3, 5] (Frobenius_filter not applicable)',
    ]);
    // SageMath's post-filter answer is [2, 3]; the filter needs the primes of
    // the degree-6 number field, which is not ported, so we return the
    // (sufficient, non-minimal) unfiltered set.  See the deviation note on
    // Frobenius_filter.
    expect(result).toEqual([2n, 3n, 5n]);
  });

  // The class-group step (isogeny_class.py:1309-1317).  For a class group with
  // a single non-principal class (h = 2) or a single class up to inverses,
  // PARI's quadclassunit generator is forced, so these values are SageMath's.
  // small_prime_value is sage/quadratic_forms/binary_qf.py:1572.
  it('finds a prime represented by each class group generator', () => {
    const gensLine = (d: bigint): string => {
      const lines: string[] = [];
      const log = console.log;
      console.log = (...args: unknown[]) => {
        lines.push(args.join(' '));
      };
      try {
        isogeny_degrees_cm(
          asAnyCurve({
            has_cm: () => true,
            cm_discriminant: () => d,
            has_rational_cm: () => true,
            base_field: () => ({ degree: () => 12 }),
          }),
          true
        );
      } finally {
        console.log = log;
      }
      return lines.find((l) => l.startsWith('primes generating')) ?? '';
    };
    // class number 1: no generators at all
    for (const d of [-3n, -4n, -7n, -8n, -11n, -19n, -43n, -67n, -163n]) {
      expect(`${d}: ${gensLine(d)}`).toBe(`${d}: primes generating the class group: []`);
    }
    // h = 2, one non-principal class: [2,1,2] for -15, [2,2,3] for -20,
    // [3,2,3] for -32, [3,1,3] for -35, [3,3,5] for -51
    for (const [d, l] of [
      [-15n, 2n],
      [-20n, 2n],
      [-32n, 3n],
      [-35n, 3n],
      [-51n, 3n],
    ] as const) {
      expect(`${d}: ${gensLine(d)}`).toBe(`${d}: primes generating the class group: [${l}]`);
    }
    // isogeny_class.py:1200 -- SageMath prints exactly this for d = -23
    expect(gensLine(-23n)).toBe('primes generating the class group: [2]');
  });
});

describe('Frobenius_filter', () => {
  const curve = (ainvs: bigint[]) => asAnyCurve({ ainvs: () => ainvs });
  const primes40 = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n];

  // gal_reps_number_field.py:520-522
  //   sage: E = EllipticCurve('11a1')
  //   sage: Frobenius_filter(E, primes(40))
  //   [5]
  it("matches SageMath's doctest for 11a1", () => {
    expect(Frobenius_filter(curve([0n, -1n, 1n, -10n, -20n]), primes40)).toEqual([5n]);
  });

  // gal_reps_number_field.py:528-530
  //   sage: E = EllipticCurve_from_j(2268945/128)
  //   sage: Frobenius_filter(E, [7, 11])
  //   [7]
  // j determines E up to quadratic twist, and disc(x^2 - a_p*x + p) = a_p^2 - 4p
  // is twist invariant, so any model with this j gives the same answer.  Here
  // y^2 = x^3 - 3*j*(j-1728)*x - 2*j*(j-1728)^2 scaled by u = 16 to clear the
  // denominator 128.
  it("matches SageMath's doctest for j = 2268945/128", () => {
    const jn = 2268945n;
    const jd = 128n;
    const kn = 2268945n - 1728n * 128n; // j - 1728 = kn/jd
    const u = 16n;
    const a4 = (-3n * jn * kn * u ** 4n) / (jd * jd);
    const a6 = (-2n * jn * kn * kn * u ** 6n) / (jd * jd * jd);
    expect(Frobenius_filter(curve([0n, 0n, 0n, a4, a6]), [7n, 11n])).toEqual([7n]);
  });

  // The filter is never allowed to discard a prime for which the mod-l
  // representation really is reducible; these isogeny degrees are Cremona's.
  it('keeps every prime that really is an isogeny degree', () => {
    const cases: Array<[string, bigint[], bigint[]]> = [
      ['11a1', [0n, -1n, 1n, -10n, -20n], [5n]],
      ['14a1', [1n, 0n, 1n, 4n, -6n], [2n, 3n]],
      ['15a1', [1n, 1n, 1n, -10n, -10n], [2n]],
      ['17a1', [1n, -1n, 1n, -1n, -14n], [2n]],
      ['27a1', [0n, 0n, 1n, 0n, -7n], [3n]],
      ['37b1', [0n, 1n, 1n, -23n, -50n], [3n]],
      ['26b1', [1n, -1n, 1n, -3n, 3n], [7n]],
      ['50a1', [1n, 0n, 1n, -1n, -2n], [3n, 5n]],
    ];
    for (const [, ainvs, degrees] of cases) {
      const got = Frobenius_filter(curve(ainvs), primes40);
      // exactly the known reducible primes survive
      expect(got).toEqual(degrees);
    }
  });

  // 37a1 has a trivial isogeny class: every mod-l representation is surjective.
  it('discards everything for a curve with no isogenies', () => {
    expect(Frobenius_filter(curve([0n, 0n, 1n, -1n, 0n]), primes40)).toEqual([]);
  });

  it('clears denominators (global integral model)', () => {
    // y^2 = x^3 - x/4 - 1/4 is the model of a curve with rational coefficients;
    // scaling by u = 2 gives the integral [0,0,0,-4,-16], the same curve.
    expect(
      Frobenius_filter(
        asAnyCurve({ ainvs: () => [0n, 0n, 0n, { num: -1n, den: 4n }, { num: -1n, den: 4n }] }),
        primes40
      )
    ).toEqual(Frobenius_filter(curve([0n, 0n, 0n, -4n, -16n]), primes40));
  });

  // The include_2 branch (gal_reps_number_field.py:556-559) needs an exact
  // decision on whether 4*x^3 + b2*x^2 + 2*b4*x + b6 has a rational root.
  // Compared here against a brute-force search over the integers.
  it('decides 2-division-polynomial reducibility exactly', () => {
    const abs = (x: bigint) => (x < 0n ? -x : x);
    const bruteHasIntRoot = (a2: bigint, a4: bigint, a6: bigint): boolean => {
      if (a6 === 0n) return true;
      for (let r = -abs(a6); r <= abs(a6); r++) {
        if (r * r * r + a2 * r * r + a4 * r + a6 === 0n) return true;
      }
      return false;
    };
    let seed = 12345n;
    const rnd = (m: bigint) => {
      seed = (seed * 6364136223846793005n + 1442695040888963407n) & ((1n << 63n) - 1n);
      return (seed >> 20n) % m;
    };
    let checked = 0;
    for (let i = 0; i < 300; i++) {
      const a2 = rnd(401n) - 200n;
      const a4 = rnd(401n) - 200n;
      const a6 = rnd(401n) - 200n;
      const b2 = 4n * a2;
      const b4 = 2n * a4;
      const b6 = 4n * a6;
      const b8 = 4n * a2 * a6 - a4 * a4;
      const disc = -b2 * b2 * b8 - 8n * b4 * b4 * b4 - 27n * b6 * b6 + 9n * b2 * b4 * b6;
      if (disc === 0n) continue;
      checked++;
      // With l = [2] the filter returns [2] exactly when the 2-division
      // polynomial 4*(x^3 + a2*x^2 + a4*x + a6) is reducible over Q.
      const got = Frobenius_filter(curve([0n, a2, 0n, a4, a6]), [2n]);
      expect(`${a2},${a4},${a6}: ${got.length === 1}`).toBe(
        `${a2},${a4},${a6}: ${bruteHasIntRoot(a2, a4, a6)}`
      );
    }
    expect(checked).toBeGreaterThan(250);
  });

  it('refuses base fields other than QQ', () => {
    expect(() =>
      Frobenius_filter(asAnyCurve({ ainvs: () => [0n, 0n, 0n, 'a', 0n] }), [3n])
    ).toThrow(NotImplementedError);
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

// ---------------------------------------------------------------------------
// M104: IsogenyClass.matrix() must delegate to fill_isogeny_matrix /
// unfill_isogeny_matrix, which take the *minimal* path degree.  The private
// helper it used before took the first path it found, so Sage's 6x6 doctest
// matrix produced 36 where the minimum is 9.
// ---------------------------------------------------------------------------

import { fill_isogeny_matrix, unfill_isogeny_matrix } from './ell_curve_isogeny.js';

describe('isogeny matrix helpers used by IsogenyClass.matrix (M104)', () => {
  const M: bigint[][] = [
    [0n, 2n, 3n, 3n, 0n, 0n],
    [2n, 0n, 0n, 0n, 3n, 3n],
    [3n, 0n, 0n, 0n, 2n, 0n],
    [3n, 0n, 0n, 0n, 0n, 2n],
    [0n, 3n, 2n, 0n, 0n, 0n],
    [0n, 3n, 0n, 2n, 0n, 0n],
  ];

  // ell_curve_isogeny.py:fill_isogeny_matrix doctest
  it('takes minimal path degrees', () => {
    const filled = fill_isogeny_matrix(M);
    expect(filled).toEqual([
      [1n, 2n, 3n, 3n, 6n, 6n],
      [2n, 1n, 6n, 6n, 3n, 3n],
      [3n, 6n, 1n, 9n, 2n, 18n],
      [3n, 6n, 9n, 1n, 18n, 2n],
      [6n, 3n, 2n, 18n, 1n, 9n],
      [6n, 3n, 18n, 2n, 9n, 1n],
    ]);
    // the old private helper returned 36 here
    expect(filled[2]![3]).toBe(9n);
  });

  it('round-trips through unfill', () => {
    expect(unfill_isogeny_matrix(fill_isogeny_matrix(M))).toEqual(M);
  });
});
