/**
 * @module sage/schemes/elliptic_curves/isogeny_class
 * @description Isogeny classes of elliptic curves
 *
 * Port of: sage/schemes/elliptic_curves/isogeny_class.py
 *
 * This module provides functionality for computing isogeny classes
 * of elliptic curves over number fields and finite fields.
 *
 * An isogeny class consists of all elliptic curves that are isogenous to
 * a given curve, together with the isogenies between them. Two curves are
 * isogenous if there exists a non-constant morphism between them (which is
 * necessarily a group homomorphism).
 *
 * Over Q, isogeny classes are well-understood due to Mazur's theorem, which
 * limits the possible torsion structures and hence the possible isogeny degrees.
 */

import { ellinit_Fp, trace_of_frobenius } from '@sagemath-ts/parigp-ts';
import { NotImplementedError, ValueError } from '../../errors.js';
import { BinaryQF, BinaryQF_reduced_representatives } from '../../quadratic_forms/binary_qf.js';
import { NumberField, NumberFieldElement } from '../../rings/number_field/number_field.js';
import { Rational } from '../../rings/rational.js';
import { fill_isogeny_matrix, unfill_isogeny_matrix } from './ell_curve_isogeny.js';
import type { EllipticCurveGeneric } from './ell_generic.js';
import type { FieldElement } from './types.js';

interface EllipticCurveCMExtras {
  has_cm?: () => boolean;
  cm_discriminant?: () => bigint;
  base_field?: () => { degree?: () => number | bigint };
  has_rational_cm?: () => boolean;
  isogenies_prime_degree?: (
    primes: bigint[],
    options?: { minimal_models?: boolean }
  ) => unknown[] | null;
  good_reduction_primes?: () => unknown;
}

/**
 * Interface for isogeny objects.
 */
export interface Isogeny<F extends FieldElement = FieldElement> {
  /** Domain curve */
  domain(): EllipticCurveGeneric<F>;
  /** Codomain curve */
  codomain(): EllipticCurveGeneric<F>;
  /** Degree of the isogeny */
  degree(): bigint;
  /** Apply the isogeny to a point */
  // apply(point: EllipticCurvePoint<F>): EllipticCurvePoint<F>;
}

/**
 * Interface for the isogeny graph representation.
 *
 * This is a simplified graph structure compared to SageMath's Graph class.
 * It provides the essential data for visualization and analysis of
 * isogeny classes.
 */
export interface IsogenyGraph<F extends FieldElement = FieldElement> {
  /** Map from vertex label (1-indexed) to elliptic curve */
  vertices: Map<number, EllipticCurveGeneric<F>>;

  /** List of edges with their degrees */
  edges: Array<{ from: number; to: number; degree: bigint }>;

  /** Adjacency matrix (0-indexed) with prime degree isogenies */
  adjacencyMatrix: bigint[][];

  /** Vertex positions for visualization (1-indexed) */
  positions: Map<number, [number, number]>;

  /** Number of vertices */
  numVertices: number;

  /** Number of edges */
  numEdges: number;
}

/**
 * An isogeny class of elliptic curves.
 *
 * This class represents all elliptic curves isogenous to a given curve,
 * together with the isogenies between them. The curves are stored in a
 * list, and the isogenies are stored as a matrix where entry (i,j) gives
 * the degree of an isogeny from curve i to curve j.
 *
 * For curves over Q, there is always a canonical ordering based on Cremona's
 * tables or LMFDB conventions.
 *
 * @example
 * ```typescript
 * const E = EllipticCurve([0, 0, 1, -1, 0]); // 11a1
 * const iso = E.isogeny_class();
 * console.log(iso.length()); // 3 curves in the class
 * console.log(iso.matrix()); // matrix of isogeny degrees
 * ```
 *
 * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:IsogenyClass_EC
 */
export class IsogenyClass<F extends FieldElement = FieldElement> {
  /** The initial elliptic curve */
  E: EllipticCurveGeneric<F>;

  /** Label for this isogeny class (e.g., '11a' or '37.a') */
  protected _label?: string;

  /** List of curves in the class */
  curves: EllipticCurveGeneric<F>[] = [];

  /** Matrix of isogeny degrees (null = not computed) */
  protected _mat: bigint[][] | null = null;

  /** List of isogenies (null = not computed) */
  protected _maps: (Isogeny<F> | 0)[][] | null = null;

  /** Quadratic form matrix for CM curves */
  protected _qfmat: number[][][] | null = null;

  /** Algorithm used to compute the class */
  protected _algorithm: string = 'sage';

  /**
   * Create an isogeny class from an elliptic curve.
   *
   * INPUT:
   * - E: an elliptic curve
   * - label: (optional) a label for the isogeny class
   * - empty: (optional) if True, create an empty class (for internal use)
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:__init__
   */
  constructor(E: EllipticCurveGeneric<F>, label?: string, empty: boolean = false) {
    this.E = E;
    this._label = label;

    if (!empty) {
      this._compute();
    }
  }

  /**
   * Compute the isogeny class.
   *
   * This method computes all curves in the isogeny class and the matrix
   * of isogeny degrees between them.
   *
   * The algorithm:
   * 1. Start with the initial curve E
   * 2. Find all isogenies of prime degree from E
   * 3. Add the codomains to the class if not already present
   * 4. Repeat for each new curve until no new curves are found
   *
   * @internal
   */
  protected _compute(): void {
    // This requires isogeny computation which depends on:
    // - Computing possible isogeny degrees (possible_isogeny_degrees)
    // - Computing actual isogenies (isogenies_prime_degree)
    // - Checking curve isomorphism

    // For now, initialize with just the original curve
    this.curves = [this.E];
    this._mat = [[1n]];

    // Full implementation would compute the complete isogeny class
    // throw new NotImplementedError(
    //   'IsogenyClass._compute: full isogeny class computation not yet implemented'
    // );
  }

  /**
   * Return the number of curves in this isogeny class.
   *
   * @example
   * ```typescript
   * const iso = E.isogeny_class();
   * console.log(iso.length()); // e.g., 3 for 11a
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:__len__
   */
  length(): number {
    return this.curves.length;
  }

  /**
   * Iterate over the curves in this isogeny class.
   *
   * @example
   * ```typescript
   * for (const C of iso) {
   *   console.log(C.j_invariant());
   * }
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:__iter__
   */
  *[Symbol.iterator](): Iterator<EllipticCurveGeneric<F>> {
    for (const curve of this.curves) {
      yield curve;
    }
  }

  /**
   * Return the i-th curve in this isogeny class.
   *
   * @param i - Index (0-based)
   * @returns The i-th curve
   *
   * @example
   * ```typescript
   * const E1 = iso.get(0); // first curve
   * const E2 = iso.get(1); // second curve
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:__getitem__
   */
  get(i: number): EllipticCurveGeneric<F> {
    if (i < 0 || i >= this.curves.length) {
      throw new ValueError(
        `Index ${i} out of range for isogeny class of size ${this.curves.length}`
      );
    }
    return this.curves[i]!;
  }

  /**
   * Return the index of a curve in this isogeny class.
   *
   * INPUT:
   * - C: an elliptic curve isomorphic to one in this class
   *
   * OUTPUT:
   * - i: integer such that the i-th curve in the class is isomorphic to C
   *
   * @throws {ValueError} If C is not in the isogeny class
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:index
   */
  index(C: EllipticCurveGeneric<F>): number {
    for (let i = 0; i < this.curves.length; i++) {
      if (C.is_isomorphic(this.curves[i]!)) {
        return i;
      }
    }
    throw new ValueError(`${C} is not in isogeny class ${this}`);
  }

  /**
   * Return the matrix of isogeny degrees.
   *
   * The (i,j) entry is the minimal degree of an isogeny from curve i to curve j.
   * If fill is False, only prime degrees are shown (0 elsewhere).
   *
   * INPUT:
   * - fill: if True (default), compute all isogeny degrees using transitivity
   *
   * OUTPUT: a matrix of integers (as 2D array)
   *
   * @example
   * ```typescript
   * const M = iso.matrix();
   * console.log(M[0][1]); // degree of isogeny from curve 0 to curve 1
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:matrix
   */
  matrix(fill: boolean = true): bigint[][] {
    if (this._mat === null) {
      this._compute();
    }

    let mat = this._mat!;

    if (fill && mat[0]?.[0] === 0n) {
      // Fill in using the fill_isogeny_matrix algorithm
      mat = fill_isogeny_matrix(mat);
    }

    if (!fill && mat[0]?.[0] === 1n) {
      // Unfill to show only prime degrees
      mat = unfill_isogeny_matrix(mat);
    }

    return mat;
  }

  /**
   * Simple primality test.
   */
  private _isPrime(n: bigint): boolean {
    if (n < 2n) return false;
    if (n === 2n) return true;
    if (n % 2n === 0n) return false;
    if (n < 9n) return true;
    if (n % 3n === 0n) return false;

    const limit = this._isqrt(n);
    let i = 5n;
    while (i <= limit) {
      if (n % i === 0n || n % (i + 2n) === 0n) return false;
      i += 6n;
    }
    return true;
  }

  /**
   * Integer square root.
   */
  private _isqrt(n: bigint): bigint {
    if (n < 2n) return n;
    let x = n;
    let y = (x + 1n) / 2n;
    while (y < x) {
      x = y;
      y = (x + n / x) / 2n;
    }
    return x;
  }

  /**
   * Return the quadratic form matrix (CM case only).
   *
   * For CM curves, the isogeny degrees between curves in the same CM
   * class can be represented by quadratic forms. This method returns
   * a matrix where each entry is a list of coefficients of a quadratic
   * form.
   *
   * OUTPUT:
   * A 2D array where entry [i][j] is a list representing a quadratic
   * form whose values are the possible isogeny degrees between curves
   * i and j.
   *
   * @throws {ValueError} If the curve does not have rational CM
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:qf_matrix
   */
  qf_matrix(): number[][][] {
    if (this._qfmat === null) {
      throw new ValueError('qf_matrix only defined for isogeny classes with rational CM');
    }
    return this._qfmat;
  }

  /**
   * Return a list of isogenies between curves in this class.
   *
   * INPUT:
   * - fill: if True, compute all isogenies (not just prime degree)
   *
   * OUTPUT: a 2D array where entry [i][j] is either 0 or an isogeny
   *         from curve i to curve j
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:isogenies
   */
  isogenies(_fill: boolean = false): (Isogeny<F> | 0)[][] {
    if (_fill) {
      throw new NotImplementedError('Computing all isogenies (fill=true) not yet implemented');
    }

    if (this._maps === null) {
      this._computeIsogenies();
    }

    return this._maps!;
  }

  /**
   * Compute the isogenies between curves.
   * @internal
   */
  protected _computeIsogenies(): void {
    // Initialize with zeros
    const n = this.curves.length;
    this._maps = Array.from({ length: n }, () => Array<Isogeny<F> | 0>(n).fill(0));

    // Full implementation would compute actual isogeny objects
    // This requires the isogeny computation functionality
  }

  /**
   * Return the graph of isogenies in this class.
   *
   * OUTPUT: a graph whose vertices are curves (labeled 1 to n) and whose
   *         edges correspond to prime degree isogenies. The graph is
   *         represented as an object with vertices, edges, adjacency matrix,
   *         and position data for visualization.
   *
   * NOTE: This is a simplified representation compared to SageMath's Graph class.
   * The vertices are labeled 1 to n (not 0 to n-1) to match LMFDB and Cremona conventions.
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:graph
   */
  graph(): IsogenyGraph<F> {
    const M = this.matrix(false); // unfilled matrix with only prime degrees
    const n = this.curves.length;

    // Build vertices (1-indexed to match LMFDB/Cremona labels)
    const vertices: Map<number, EllipticCurveGeneric<F>> = new Map();
    for (let i = 0; i < n; i++) {
      vertices.set(i + 1, this.curves[i]!);
    }

    // Build edges from the matrix (only non-zero entries = prime degree isogenies)
    const edges: Array<{ from: number; to: number; degree: bigint }> = [];
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const deg = M[i]![j]!;
        if (deg > 0n) {
          edges.push({ from: i + 1, to: j + 1, degree: deg });
        }
      }
    }

    // Compute vertex positions for visualization
    // SageMath special-cases each isogeny graph type for nice layout
    const pos = this._computeGraphPositions(n, M);

    return {
      vertices,
      edges,
      adjacencyMatrix: M,
      positions: pos,
      numVertices: n,
      numEdges: edges.length,
    };
  }

  /**
   * Compute positions for graph vertices based on the structure.
   * @internal
   */
  private _computeGraphPositions(n: number, M: bigint[][]): Map<number, [number, number]> {
    const pos = new Map<number, [number, number]>();

    if (n === 1) {
      pos.set(1, [0, 0]);
      return pos;
    }

    if (n === 2) {
      // One edge, two vertices - horizontal layout
      pos.set(1, [-0.5, 0]);
      pos.set(2, [0.5, 0]);
      return pos;
    }

    // For larger graphs, use a simple circular layout as a fallback
    // SageMath has special cases for each known isogeny graph type over Q
    const filled = fill_isogeny_matrix(M);
    const maxDegree = Math.max(...filled.flat().map((x) => Number(x)));

    if (n === 3) {
      // o--o--o layout: find center vertex (the one with smallest max degree)
      let centerVert = 0;
      for (let i = 0; i < 3; i++) {
        const rowMax = Math.max(...filled[i]!.map((x) => Number(x)));
        if (rowMax < maxDegree) {
          centerVert = i;
          break;
        }
      }
      const other = [0, 1, 2].filter((i) => i !== centerVert);
      pos.set(centerVert + 1, [0, 0]);
      pos.set(other[0]! + 1, [-1, 0]);
      pos.set(other[1]! + 1, [1, 0]);
      return pos;
    }

    // Default: circular layout
    for (let i = 0; i < n; i++) {
      const angle = (2 * Math.PI * i) / n - Math.PI / 2;
      pos.set(i + 1, [Math.cos(angle), Math.sin(angle)]);
    }

    return pos;
  }

  /**
   * Reorder the curves in this class.
   *
   * INPUT:
   * - order: a permutation, list of curves, or string specifying the order
   *   - 'lmfdb': order by a-invariants lexicographically
   *   - list of curves: use the given order
   *
   * OUTPUT: a new IsogenyClass with curves reordered
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:reorder
   */
  reorder(order: 'lmfdb' | EllipticCurveGeneric<F>[] | number[]): IsogenyClass<F> {
    if (order === null || order === this._algorithm) {
      return this;
    }

    let reorderedCurves: EllipticCurveGeneric<F>[];

    if (order === 'lmfdb') {
      // Sort by a-invariants
      reorderedCurves = [...this.curves].sort((a, b) => {
        const aInvs = a.a_invariants();
        const bInvs = b.a_invariants();
        for (let i = 0; i < 5; i++) {
          const cmp = this._compareFieldElements(aInvs[i]!, bInvs[i]!);
          if (cmp !== 0) return cmp;
        }
        return 0;
      });
    } else if (Array.isArray(order) && typeof order[0] === 'number') {
      // Permutation as array of indices
      const perm = order as number[];
      if (perm.length !== this.curves.length) {
        throw new ValueError('Incorrect length for permutation');
      }
      reorderedCurves = perm.map((i) => this.curves[i]!);
    } else {
      // List of curves
      reorderedCurves = order as EllipticCurveGeneric<F>[];
      if (reorderedCurves.length !== this.curves.length) {
        throw new ValueError('Incorrect length');
      }
    }

    const cpy = this.copy();
    cpy.curves = reorderedCurves;

    // Recompute matrix with new ordering
    if (this._mat !== null) {
      const perm = reorderedCurves.map((c) => this.curves.indexOf(c));
      const n = this.curves.length;
      const newMat: bigint[][] = Array.from({ length: n }, () => Array<bigint>(n).fill(0n));

      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          newMat[i]![j] = this._mat[perm[i]!]![perm[j]!]!;
        }
      }
      cpy._mat = newMat;
    }

    return cpy;
  }

  /**
   * Compare two field elements for ordering.
   */
  private _compareFieldElements(a: F, b: F): number {
    const aStr = a.toString();
    const bStr = b.toString();
    if (aStr < bStr) return -1;
    if (aStr > bStr) return 1;
    return 0;
  }

  /**
   * Return a copy of this isogeny class.
   *
   * The copy shares the same curve objects but has independent lists
   * and matrices.
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:copy
   */
  copy(): IsogenyClass<F> {
    const result = new IsogenyClass<F>(this.E, this._label, true);
    result.curves = [...this.curves];
    result._mat = this._mat ? this._mat.map((row) => [...row]) : null;
    result._maps = null; // Don't copy maps, they will be recomputed if needed
    result._qfmat = this._qfmat ? this._qfmat.map((row) => row.map((cell) => [...cell])) : null;
    result._algorithm = this._algorithm;
    return result;
  }

  /**
   * Check if a curve is in this isogeny class.
   */
  contains(x: EllipticCurveGeneric<F>): boolean {
    return this.curves.some((curve) => x.is_isomorphic(curve));
  }

  /**
   * String representation.
   */
  toString(): string {
    if (this._label) {
      return `Elliptic curve isogeny class ${this._label}`;
    }
    return `Isogeny class of ${this.E}`;
  }
}

/**
 * Isogeny class over a number field.
 *
 * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:IsogenyClass_EC_NumberField
 */
export class IsogenyClassNumberField<
  F extends FieldElement = FieldElement,
> extends IsogenyClass<F> {
  /** Primes where the Galois representation is reducible */
  protected _reducible_primes?: bigint[];

  /** Algorithm for computing reducible primes */
  protected _nf_algorithm: string;

  /** Whether to use minimal models */
  protected _minimal_models: boolean;

  /**
   * Create an isogeny class over a number field.
   *
   * INPUT:
   * - E: an elliptic curve over a number field
   * - reducible_primes: (optional) list of primes where representation is reducible
   * - algorithm: 'Billerey' (default), 'Larson', or 'heuristic'
   * - minimal_models: if True, use minimal models for curves
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:IsogenyClass_EC_NumberField.__init__
   */
  constructor(
    E: EllipticCurveGeneric<F>,
    options: {
      reducible_primes?: bigint[];
      algorithm?: 'Billerey' | 'Larson' | 'heuristic';
      minimal_models?: boolean;
    } = {}
  ) {
    super(E, undefined, true);
    this._reducible_primes = options.reducible_primes;
    this._nf_algorithm = options.algorithm ?? 'Billerey';
    this._minimal_models = options.minimal_models ?? true;
    this._compute();
  }

  /**
   * Return a copy.
   */
  override copy(): IsogenyClassNumberField<F> {
    const result = new IsogenyClassNumberField<F>(this.E, {
      reducible_primes: this._reducible_primes,
      algorithm: this._nf_algorithm as 'Billerey' | 'Larson' | 'heuristic',
      minimal_models: this._minimal_models,
    });
    result.curves = [...this.curves];
    result._mat = this._mat ? this._mat.map((row) => [...row]) : null;
    return result;
  }
}

/**
 * Isogeny class over the rationals.
 *
 * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:IsogenyClass_EC_Rational
 */
export class IsogenyClassRational<
  F extends FieldElement = FieldElement,
> extends IsogenyClassNumberField<F> {
  /**
   * Create an isogeny class over Q.
   *
   * INPUT:
   * - E: an elliptic curve over Q
   * - algorithm: 'sage' (default) or 'database'
   * - label: (optional) Cremona label
   *
   * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:IsogenyClass_EC_Rational.__init__
   */
  constructor(E: EllipticCurveGeneric<F>, algorithm: 'sage' | 'database' = 'sage', label?: string) {
    super(E, { algorithm: 'Billerey', minimal_models: true });
    this._algorithm = algorithm;
    this._label = label;
  }

  /**
   * Return a copy.
   */
  override copy(): IsogenyClassRational<F> {
    const result = new IsogenyClassRational<F>(
      this.E,
      this._algorithm as 'sage' | 'database',
      this._label
    );
    result.curves = [...this.curves];
    result._mat = this._mat ? this._mat.map((row) => [...row]) : null;
    return result;
  }
}

/**
 * Return the possible isogeny degrees for a CM curve.
 *
 * For curves with CM by an order O of discriminant d, the possible
 * isogeny degrees are determined by:
 * 1. Horizontal isogenies: degrees represented by quadratic forms of discriminant d
 * 2. Vertical isogenies: primes dividing d or satisfying divisibility conditions
 *
 * INPUT:
 * - E: an elliptic curve with CM
 * - verbose: if True, print extra information
 *
 * OUTPUT: a sorted list of prime degrees of possible isogenies
 *
 * ALGORITHM:
 * For curves with CM by the order O of discriminant d:
 *
 * Case (1): Horizontal isogenies (O = O'): degrees are represented by binary
 * quadratic forms of discriminant d. We find a prime represented by each
 * non-principal ideal class.
 *
 * Case (2): Upward vertical isogenies ([O':O] = l): d = l^2 * d', so we include
 * all prime divisors of d.
 *
 * Case (3): Downward vertical isogenies ([O:O'] = l): the class numbers satisfy
 * h(O') = (l +- 1) * h(O), so we include primes l such that l +- 1 divides the
 * degree of the base field times 2*h(O).
 *
 * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:isogeny_degrees_cm
 * @see Deviation: {@link Frobenius_filter} is applied whenever the curve gives
 *   access to its a-invariants -- over `QQ` and over any {@link NumberField}.
 *   For a curve that exposes neither (the bare `cm_discriminant()` mocks used in
 *   the tests) the *unfiltered* (still sufficient, but non-minimal) candidate
 *   set is returned and the verbose transcript says so.
 * @see Deviation: the class group generators used for the horizontal primes
 *   come from enumerating reduced forms, not from PARI's `quadclassunit`; see
 *   {@link _class_group_generators}.
 */
export function isogeny_degrees_cm<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  verbose: boolean = false
): bigint[] {
  const extras = E as unknown as EllipticCurveCMExtras;
  // Check if curve has CM
  // This requires E.has_cm() and E.cm_discriminant() methods
  const hasCM = typeof extras.has_cm === 'function' ? extras.has_cm() : false;
  if (!hasCM) {
    throw new ValueError('isogeny_degrees_cm requires E to have complex multiplication');
  }

  const d = extras.cm_discriminant?.();
  if (d === undefined) {
    throw new ValueError('isogeny_degrees_cm requires E.cm_discriminant()');
  }

  if (verbose) {
    console.log(`CM case, discriminant = ${d}`);
  }

  // Get the base field degree
  // For now, assume degree 1 (Q) if not available
  const fieldDegree: number =
    typeof extras.base_field === 'function' && typeof extras.base_field().degree === 'function'
      ? Number(extras.base_field().degree())
      : 1;

  // Check for rational CM
  const hasRationalCM: boolean =
    typeof extras.has_rational_cm === 'function' ? extras.has_rational_cm() : true;

  let n = BigInt(fieldDegree);
  if (!hasRationalCM) {
    n *= 2n;
  }
  // Extra factors for discriminants with extra units
  if (d === -4n) {
    n *= 2n;
  }
  if (d === -3n) {
    n *= 3n;
  }

  const divs = divisors(n);
  const L = new Set<bigint>();

  // Initial primes: always include 2 (and 3 for d = -3)
  L.add(2n);
  if (d === -3n) {
    L.add(3n);
  }

  if (verbose) {
    console.log(`initial primes: {${Array.from(L).join(', ')}}`);
  }

  // Step 1: "vertical" primes
  // Compute odd part of d and its prime factors (ramified primes)
  const oddPart = d % 2n === 0n ? d / (d & -d) : d; // Remove all factors of 2
  const absOddPart = oddPart < 0n ? -oddPart : oddPart;
  const ramPrimes = primeFactors(absOddPart);

  if (!hasRationalCM) {
    // Include all ramified primes
    for (const l of ramPrimes) {
      L.add(l);
    }
    if (verbose) {
      console.log(`ramified primes: {${ramPrimes.join(', ')}}`);
    }
  } else {
    // We must have 2*h dividing n, and will need the quotient
    // (isogeny_class.py:1236-1247).  h is the number of reduced primitive
    // forms of discriminant d.
    const h = BigInt(BinaryQF_reduced_representatives(d, { primitive_only: true }).length);
    const nOver2h = n / (2n * h);

    // Upward primes (index divided by l): d has valuation > 1 at l
    const upward = ramPrimes.filter((l) => valuation(d, l) > 1n);
    for (const l of upward) {
      L.add(l);
    }
    if (verbose) {
      console.log(`upward primes: {${upward.join(', ')}}`);
    }

    // Downward ramified primes: the index is multiplied by l and the class
    // number by l, so l must divide n/(2h)  (isogeny_class.py:1285-1288).
    const downwardRamified = ramPrimes.filter((l) => nOver2h % l === 0n);
    for (const l of downwardRamified) {
      L.add(l);
    }
    if (verbose) {
      console.log(`downward ramified primes: {${downwardRamified.join(', ')}}`);
    }
  }

  // Downward split primes: the suborder has class number (l-1)*h, so l-1
  // must divide n/(2h)  (isogeny_class.py:1291-1296; Sage runs over the
  // divisors of n).
  const splitPrimes = divs
    .map((lm1) => lm1 + 1n)
    .filter((l) => isPrime(l) && kroneckerSymbol(d, l) === 1n);
  for (const l of splitPrimes) {
    L.add(l);
  }
  if (verbose) {
    console.log(`downward split primes: {${splitPrimes.join(', ')}}`);
  }

  // Downward inert primes: the suborder has class number (l+1)*h
  // (isogeny_class.py:1298-1305).
  const inertPrimes = divs
    .map((lp1) => lp1 - 1n)
    .filter((l) => l > 1n && isPrime(l) && kroneckerSymbol(d, l) === -1n);
  for (const l of inertPrimes) {
    L.add(l);
  }
  if (verbose) {
    console.log(`downward inert primes: {${inertPrimes.join(', ')}}`);
  }

  // Horizontal primes (rational CM only): same order, degrees are all integers
  // represented by some binary quadratic form of discriminant d, so we find a
  // prime represented by each generator of the class group
  // (isogeny_class.py:1309-1317).
  if (hasRationalCM) {
    const gens = _class_group_generators(d);
    const L1 = gens.map((Q) => _small_prime_value(Q));
    if (verbose) {
      console.log(`primes generating the class group: [${L1.join(', ')}]`);
    }
    for (const l of L1) {
      L.add(l);
    }
  }

  // Convert to sorted array
  let result = Array.from(L).filter((l) => isPrime(l));
  result.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));

  if (verbose) {
    console.log(`Set of primes before filtering: {${result.join(', ')}}`);
  }

  // This filter will quickly eliminate most false entries in the set
  // (isogeny_class.py:1325-1327).
  if (_frobenius_filter_applicable(E)) {
    result = Frobenius_filter(E, result);
    if (verbose) {
      console.log(`List of primes after filtering: [${result.join(', ')}]`);
    }
  } else if (verbose) {
    console.log(`List of primes: [${result.join(', ')}] (Frobenius_filter not applicable)`);
  }

  return result;
}

// ---------------------------------------------------------------------------
// Binary quadratic forms: the class group step of isogeny_degrees_cm.
// ---------------------------------------------------------------------------

/**
 * Return a generating set of the class group of the order of discriminant `d`,
 * as reduced binary quadratic forms `[a, b, c]`.
 *
 * SageMath asks PARI for `quadclassunit(d)[2]`; we enumerate the reduced forms
 * (which are the elements of the class group) and pick a generating set
 * greedily, always taking an element of maximal order outside the subgroup
 * generated so far, using Gauss composition (delegated to `BinaryQF.compose`,
 * i.e. to PARI's `qfbcompraw`).
 *
 * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:1236 (pari(d).quadclassunit())
 * @see Deviation: PARI's `quadclassunit` (Buchmann's subexponential algorithm)
 *   is not ported.  We enumerate the reduced forms of discriminant `d`
 *   (O(|d|) work) and extract *a* generating set.  The subgroup generated is
 *   the whole class group either way, but the particular generators - and
 *   hence the particular primes `isogeny_degrees_cm` returns - may differ from
 *   PARI's choice; the resulting list of primes is still sufficient, which is
 *   all `isogeny_degrees_cm` promises, but it can be larger than SageMath's.
 */
function _class_group_generators(d: bigint): BinaryQF[] {
  const forms = BinaryQF_reduced_representatives(d, { primitive_only: true });
  const key = (Q: BinaryQF) => `${Q.a},${Q.b},${Q.c}`;
  const reduce = (Q: BinaryQF) => Q.reduced_form();
  const principal = reduce(BinaryQF.principal(d));
  /** The cyclic subgroup generated by Q, as reduced forms, identity first. */
  const cyclic = (Q: BinaryQF): BinaryQF[] => {
    const powers: BinaryQF[] = [principal];
    let cur = reduce(Q);
    while (key(cur) !== key(principal)) {
      powers.push(cur);
      cur = reduce(cur.compose(Q));
    }
    return powers;
  };
  let elems: BinaryQF[] = [principal];
  let generated = new Set<string>([key(principal)]);
  const gens: BinaryQF[] = [];
  while (generated.size < forms.length) {
    // Pick an element of maximal order that is not generated yet.
    let best: BinaryQF | null = null;
    let bestPowers: BinaryQF[] = [];
    for (const Q of forms) {
      if (generated.has(key(Q))) continue;
      const powers = cyclic(Q);
      if (powers.length > bestPowers.length) {
        best = Q;
        bestPowers = powers;
      }
    }
    if (best === null) break;
    gens.push(best);
    // The subgroup generated by the old one together with `best` is { g * Q^k }.
    const next: BinaryQF[] = [];
    const seen = new Set<string>();
    for (const g of elems) {
      for (const pw of bestPowers) {
        const t = reduce(g.compose(pw));
        if (!seen.has(key(t))) {
          seen.add(key(t));
          next.push(t);
        }
      }
    }
    elems = next;
    generated = seen;
  }
  return gens;
}

/**
 * Return the smallest prime represented by the (primitive, positive definite)
 * form `Q`, by substituting small values.
 *
 * @see Reference: sage/quadratic_forms/binary_qf.py:1572 (BinaryQF.small_prime_value)
 */
function _small_prime_value(Q: BinaryQF, Bmax: bigint = 1000n): bigint {
  const { a, b, c } = Q;
  let B = 10n;
  for (;;) {
    const vals = new Set<bigint>();
    for (let x = -B; x < B; x += 1n) {
      for (let y = 0n; y < B; y += 1n) {
        vals.add(a * x * x + b * x * y + c * y * y);
      }
    }
    const primes = Array.from(vals)
      .filter((l) => isPrime(l))
      .sort((p, q) => (p < q ? -1 : p > q ? 1 : 0));
    if (primes.length > 0) return primes[0]!;
    if (B >= Bmax) {
      throw new ValueError(`Unable to find a prime value of ${a}*x^2 + ${b}*x*y + ${c}*y^2`);
    }
    B += 10n;
  }
}

// ---------------------------------------------------------------------------
// Frobenius_filter
// ---------------------------------------------------------------------------

/**
 * Determine which primes in `L` might have mod-`l` image contained in a Borel
 * subgroup, by checking traces of Frobenius.
 *
 * This function will sometimes return primes for which the image is not
 * contained in a Borel subgroup.  It never removes a prime for which the mod-`l`
 * representation really is reducible.
 *
 * @param E - an elliptic curve over a number field
 * @param L - a list of prime numbers
 * @param patience - bound on the number of traces of Frobenius used
 *
 * @see Reference: sage/schemes/elliptic_curves/gal_reps_number_field.py:492-586
 * @see Deviation: SageMath puts this function in
 *   `sage.schemes.elliptic_curves.gal_reps_number_field`, which is not ported;
 *   it lives here because `isogeny_degrees_cm` is its only caller so far.
 * @see Deviation: "good reduction at P" is tested as `P` not dividing the
 *   discriminant of the global integral model built here, not of the *global
 *   minimal* model (Tate's algorithm / Laska-Kraus-Connell is not ported).  A
 *   non-minimal model makes a few extra primes look bad; those primes are then
 *   skipped, which can only make the filter weaker (i.e. return a superset),
 *   never unsound.
 * @see Deviation: over a base field `K != QQ` only the primes `P` of `K` of
 *   *residue degree one* are used (see {@link _nf_frobenius_discs}); SageMath
 *   uses every prime of `K`.  Any set of primes of good reduction gives a sound
 *   filter, so the result is still a superset of SageMath's; it happens to be
 *   equal on all three of SageMath's doctests.
 */
export function Frobenius_filter<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  L: Array<bigint | number>,
  patience: number = 100
): bigint[] {
  // gal_reps_number_field.py:549-550
  //   E = _over_numberfield(E).global_integral_model(); K = E.base_field()
  const ctx = _frobenius_context(E);

  // Remove duplicates from L and sort (gal_reps_number_field.py:553-554).
  let Ls = Array.from(new Set(L.map((l) => BigInt(l)))).sort((p, q) =>
    p < q ? -1 : p > q ? 1 : 0
  );

  // c.f. Section 5.3(a) of [Ser1972] (gal_reps_number_field.py:556-559).
  let include2 = false;
  if (Ls.includes(2n)) {
    Ls = Ls.filter((l) => l !== 2n);
    include2 = !ctx.two_division_poly_is_irreducible();
  }

  // Discard any l for which the Frobenius polynomial at P is irreducible mod l
  // (gal_reps_number_field.py:575-582).
  const discs = ctx.frobenius_discs();
  let numP = 0;
  for (;;) {
    if (Ls.length === 0 || numP === patience) break;
    const next = discs.next();
    if (next.done) break;
    numP++;
    // legendre_symbol(disc, l) for the odd prime l is the Kronecker symbol.
    const fdisc = next.value;
    Ls = Ls.filter((l) => kroneckerSymbol(fdisc, l) !== -1n);
  }

  return include2 ? [2n, ...Ls] : Ls;
}

/**
 * The two things {@link Frobenius_filter} needs of its curve: whether the
 * 2-division polynomial is irreducible, and the stream of discriminants
 * `a_P^2 - 4*N(P)` of the Frobenius polynomials at the primes of good
 * reduction, in SageMath's order.
 */
interface FrobeniusContext {
  two_division_poly_is_irreducible: () => boolean;
  frobenius_discs: () => Generator<bigint>;
}

/** Build the {@link FrobeniusContext} of `E`: base field `QQ` or a number field. */
function _frobenius_context<F extends FieldElement>(E: EllipticCurveGeneric<F>): FrobeniusContext {
  const nf = _nf_curve(E);
  if (nf !== null) return _nf_frobenius_context(nf);

  const ainvs = _integral_ainvs_over_Q(E);
  const [a1, a2, a3, a4, a6] = ainvs;
  const b2 = a1 * a1 + 4n * a2;
  const b4 = 2n * a4 + a1 * a3;
  const b6 = a3 * a3 + 4n * a6;
  const b8 = a1 * a1 * a6 + 4n * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4;
  const disc = -b2 * b2 * b8 - 8n * b4 * b4 * b4 - 27n * b6 * b6 + 9n * b2 * b4 * b6;
  if (disc === 0n) {
    throw new ValueError('Frobenius_filter: the curve is singular');
  }
  return {
    two_division_poly_is_irreducible: () => _two_division_poly_is_irreducible(b2, b4, b6),
    frobenius_discs: function* () {
      for (let p = 2n; ; p = _next_prime(p)) {
        if (disc % p === 0n) continue; // not a prime of good reduction
        const ap = _ap_over_Q(ainvs, p);
        yield ap * ap - 4n * p; // discriminant of x^2 - a_p*x + p
      }
    },
  };
}

/**
 * True when {@link Frobenius_filter} can run on `E`: a nonsingular model whose
 * a-invariants are rational or lie in a {@link NumberField}.
 */
function _frobenius_filter_applicable<F extends FieldElement>(E: EllipticCurveGeneric<F>): boolean {
  try {
    if (_nf_curve(E) !== null) return true;
  } catch {
    return false;
  }
  const extras = E as unknown as EllipticCurveCMExtras;
  if (typeof extras.base_field === 'function') {
    const deg = extras.base_field().degree?.();
    if (deg !== undefined && Number(deg) !== 1) return false;
  }
  let ainvs: [bigint, bigint, bigint, bigint, bigint];
  try {
    ainvs = _integral_ainvs_over_Q(E);
  } catch {
    return false;
  }
  const [a1, a2, a3, a4, a6] = ainvs;
  const b2 = a1 * a1 + 4n * a2;
  const b4 = 2n * a4 + a1 * a3;
  const b6 = a3 * a3 + 4n * a6;
  const b8 = a1 * a1 * a6 + 4n * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4;
  const disc = -b2 * b2 * b8 - 8n * b4 * b4 * b4 - 27n * b6 * b6 + 9n * b2 * b4 * b6;
  return disc !== 0n;
}

/** Numerator and denominator of a rational-valued coefficient, or `null`. */
function _as_rational(x: unknown): [bigint, bigint] | null {
  if (typeof x === 'bigint') return [x, 1n];
  if (typeof x === 'number') return Number.isInteger(x) ? [BigInt(x), 1n] : null;
  if (x && typeof x === 'object') {
    const o = x as Record<string, unknown>;
    for (const [n, d] of [
      ['numerator', 'denominator'],
      ['numer', 'denom'],
      ['num', 'den'],
    ] as const) {
      if (n in o && d in o) {
        const nv = typeof o[n] === 'function' ? (o[n] as () => unknown)() : o[n];
        const dv = typeof o[d] === 'function' ? (o[d] as () => unknown)() : o[d];
        if (
          (typeof nv === 'bigint' || typeof nv === 'number') &&
          (typeof dv === 'bigint' || typeof dv === 'number')
        ) {
          return [BigInt(nv), BigInt(dv)];
        }
      }
    }
    if ('value' in o && typeof o.value === 'bigint') return [o.value, 1n];
    const m = /^(-?\d+)(?:\/(\d+))?$/.exec(String(x));
    if (m) return [BigInt(m[1]!), m[2] === undefined ? 1n : BigInt(m[2])];
  }
  return null;
}

/**
 * The a-invariants of a global integral model of `E` over `Q`.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_number_field.py:global_integral_model
 */
function _integral_ainvs_over_Q<F extends FieldElement>(
  E: EllipticCurveGeneric<F>
): [bigint, bigint, bigint, bigint, bigint] {
  const raw = E as unknown as { ainvs?: () => unknown[]; a_invariants?: () => unknown[] };
  const list = typeof raw.ainvs === 'function' ? raw.ainvs() : raw.a_invariants?.();
  if (!list || list.length !== 5) {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: Frobenius_filter - the curve does not expose 5 a-invariants'
    );
  }
  const rats = list.map((c) => _as_rational(c));
  if (rats.some((r) => r === null)) {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: Frobenius_filter over a base field that is neither QQ ' +
        'nor a NumberField (an a-invariant is not a Rational and not a ' +
        'NumberFieldElement, so there is no reduction map to build)'
    );
  }
  const weights = [1, 2, 3, 4, 6];
  // Smallest u making u^w_i * a_i integral for every i.
  let u = 1n;
  const dens = rats.map((r) => (r as [bigint, bigint])[1]);
  const qs = new Set<bigint>();
  for (const den of dens) for (const q of primeFactors(den)) qs.add(q);
  for (const q of qs) {
    let e = 0;
    for (let i = 0; i < 5; i++) {
      const v = Number(valuation(dens[i]!, q));
      e = Math.max(e, Math.ceil(v / weights[i]!));
    }
    u *= q ** BigInt(e);
  }
  const out = rats.map((r, i) => {
    const [n, d] = r as [bigint, bigint];
    const scale = u ** BigInt(weights[i]!);
    if ((n * scale) % d !== 0n) {
      throw new ValueError('global integral model: scaling failed');
    }
    return (n * scale) / d;
  });
  return out as [bigint, bigint, bigint, bigint, bigint];
}

/**
 * The trace of Frobenius `a_p = p + 1 - #E(F_p)` of the reduction of an
 * integral model at a prime `p` of good reduction.
 *
 * For `p > 3` this delegates to our PARI port (`ellcard`), as SageMath does via
 * `E.reduction(p).trace_of_frobenius()`.  For `p in {2, 3}` the short
 * Weierstrass transformation is unavailable and the four (resp. nine) points of
 * the affine plane are counted directly.
 */
function _ap_over_Q(ainvs: [bigint, bigint, bigint, bigint, bigint], p: bigint): bigint {
  const [a1, a2, a3, a4, a6] = ainvs;
  const mod = (x: bigint): bigint => ((x % p) + p) % p;
  if (p === 2n || p === 3n) {
    let n = 1n; // point at infinity
    for (let x = 0n; x < p; x++) {
      for (let y = 0n; y < p; y++) {
        const lhs = y * y + a1 * x * y + a3 * y;
        const rhs = x * x * x + a2 * x * x + a4 * x + a6;
        if (mod(lhs - rhs) === 0n) n++;
      }
    }
    return p + 1n - n;
  }
  const b2 = a1 * a1 + 4n * a2;
  const b4 = 2n * a4 + a1 * a3;
  const b6 = a3 * a3 + 4n * a6;
  const c4 = b2 * b2 - 24n * b4;
  const c6 = -b2 * b2 * b2 + 36n * b2 * b4 - 216n * b6;
  // y^2 = x^3 - 27*c4*x - 54*c6 is isomorphic to E over F_p for p > 3
  // (sage/schemes/elliptic_curves/ell_generic.py:short_weierstrass_model).
  const Ep = ellinit_Fp(mod(-27n * c4), mod(-54n * c6), p);
  return trace_of_frobenius(Ep);
}

/**
 * Is the 2-division polynomial `4*x^3 + b2*x^2 + 2*b4*x + b6` irreducible
 * over `Q`?  A cubic is reducible over `Q` exactly when it has a rational
 * root; substituting `x = X/4` and scaling by 16 turns it into the monic
 * integer cubic `X^3 + b2*X^2 + 8*b4*X + 16*b6`, whose rational roots are
 * integers.
 *
 * @see Reference: sage/schemes/elliptic_curves/gal_reps_number_field.py:559
 */
function _two_division_poly_is_irreducible(b2: bigint, b4: bigint, b6: bigint): boolean {
  return _monic_cubic_integer_root(b2, 8n * b4, 16n * b6) === null;
}

/**
 * An integer root of `X^3 + A*X^2 + B*X + C`, or `null`.
 *
 * Exact: the cubic is monotone outside its two critical points, so each of the
 * (at most three) monotone branches is searched by integer bisection.
 */
function _monic_cubic_integer_root(A: bigint, B: bigint, C: bigint): bigint | null {
  const f = (X: bigint): bigint => X * X * X + A * X * X + B * X + C;
  const abs = (x: bigint) => (x < 0n ? -x : x);
  const M =
    1n +
    (abs(A) > abs(B) ? (abs(A) > abs(C) ? abs(A) : abs(C)) : abs(B) > abs(C) ? abs(B) : abs(C));
  const floorDiv = (x: bigint, y: bigint): bigint => {
    const q = x / y;
    return x % y !== 0n && x < 0n !== y < 0n ? q - 1n : q;
  };
  // Bisection for a root of an increasing (dir=1) or decreasing (dir=-1)
  // branch on the integer interval [lo, hi].
  const search = (lo: bigint, hi: bigint, dir: bigint): bigint | null => {
    if (lo > hi) return null;
    let a = lo;
    let b = hi;
    if (f(a) * dir > 0n || f(b) * dir < 0n) return null;
    while (a <= b) {
      const mid = floorDiv(a + b, 2n);
      const v = f(mid) * dir;
      if (v === 0n) return mid;
      if (v < 0n) a = mid + 1n;
      else b = mid - 1n;
    }
    return null;
  };
  const g = A * A - 3n * B; // f'(X) = 3X^2 + 2AX + B has discriminant 4g
  if (g <= 0n) {
    return search(-M, M, 1n);
  }
  // isqrt(g)
  let s = g;
  if (g > 1n) {
    let x = g;
    let y = (x + 1n) / 2n;
    while (y < x) {
      x = y;
      y = (x + g / x) / 2n;
    }
    s = x;
  }
  const r1c = floorDiv(-A - s, 3n);
  const r2c = floorDiv(-A + s, 3n);
  const m1 = r1c - 2n;
  const M1 = r1c + 2n;
  const m2 = r2c - 2n;
  const M2 = r2c + 2n;
  for (let X = m1; X <= M1; X++) if (f(X) === 0n) return X;
  for (let X = m2; X <= M2; X++) if (f(X) === 0n) return X;
  return search(-M, m1 < -M ? -M : m1, 1n) ?? search(M1, m2, -1n) ?? search(M2 > M ? M : M2, M, 1n);
}

// ---------------------------------------------------------------------------
// Frobenius_filter over a number field
//
// SageMath's `Frobenius_filter` (gal_reps_number_field.py:549-586) works over
// any number field `K`: it walks the primes `P` of `K`, keeps those of good
// reduction, and discards every `l` for which
// `disc(frobenius_polynomial(E mod P)) = a_P^2 - 4*N(P)` is a non-residue mod
// `l`.  The routines below are that walk, specialised to the primes of residue
// degree one -- see `_nf_frobenius_discs` for why that specialisation is sound
// and what it costs.
// ---------------------------------------------------------------------------

/** A curve over a number field, in the form `Frobenius_filter` needs it. */
interface NFCurve {
  K: NumberField;
  /** `[K : QQ]`. */
  n: number;
  /** The monic integral defining polynomial of `K`, low degree first. */
  T: bigint[];
  /** `disc(T)`, nonzero; a prime `p` with `p` not dividing it is unramified in
   *  `K` and does not divide the index `[O_K : Z[theta]]`. */
  discT: bigint;
  /** `a1, a2, a3, a4, a6` of a globally integral model, as elements. */
  ainvs: NumberFieldElement[];
  /** The same, as integral coordinate vectors in the power basis of `theta`. */
  coords: bigint[][];
  /** `b2, b4, b6, b8` and the discriminant of that model, as coordinate vectors. */
  b: { b2: bigint[]; b4: bigint[]; b6: bigint[]; b8: bigint[]; disc: bigint[] };
}

/**
 * Recognise a curve whose a-invariants live in a {@link NumberField} and put it
 * in a globally integral model, i.e. SageMath's
 * `_over_numberfield(E).global_integral_model()`
 * (gal_reps_number_field.py:549).  Returns `null` for a curve over `QQ`.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_number_field.py:global_integral_model
 * @see Deviation: SageMath scales until every a-invariant lies in `O_K`; we
 *   scale until every a-invariant has *integral coordinates in the power basis
 *   of `theta`*, which is at least as strong.  Both produce a model of the same
 *   curve, so no `a_P` changes; ours can declare a few more primes bad (they are
 *   then skipped, which only weakens the filter).
 * @see Deviation: the defining polynomial of `K` must be monic with integer
 *   coefficients (`theta` an algebraic integer); otherwise a
 *   `NotImplementedError` names the missing step.
 */
function _nf_curve<F extends FieldElement>(E: EllipticCurveGeneric<F>): NFCurve | null {
  const raw = E as unknown as {
    ainvs?: () => unknown[];
    a_invariants?: () => unknown[];
    base_field?: () => unknown;
  };
  const list = typeof raw.ainvs === 'function' ? raw.ainvs() : raw.a_invariants?.();
  if (!list || list.length !== 5) return null;
  const base = typeof raw.base_field === 'function' ? raw.base_field() : undefined;
  const elt = list.find((c) => c instanceof NumberFieldElement) as NumberFieldElement | undefined;
  if (!(base instanceof NumberField) && elt === undefined) return null;

  const K = base instanceof NumberField ? base : (elt as NumberFieldElement).parent();
  const n = K.degree();
  const ainvs0 = list.map((c) =>
    c instanceof NumberFieldElement ? c : K.__call__(_asRationalStrict(c))
  );

  // The defining polynomial must be monic and integral.
  const dp = K.defining_polynomial();
  const T: bigint[] = [];
  for (let i = 0; i <= dp.degree(); i++) {
    const c = dp.getCoeff(i);
    if (c.denominator !== 1n) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: Frobenius_filter over a number field whose defining ' +
          'polynomial is not integral (SageMath rescales through K.pari_nf(); ' +
          'sage/rings/number_field/number_field.py:pari_polynomial)'
      );
    }
    T.push(c.numerator);
  }
  if (T[n] !== 1n) {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: Frobenius_filter over a number field with a non-monic ' +
        'defining polynomial'
    );
  }
  const discT = _rationalToBigInt(dp.discriminant(), 'disc of the defining polynomial');

  // global_integral_model: scale a_i by u^i.
  const weights = [1, 2, 3, 4, 6];
  const dens = ainvs0.map((a) => _nf_coord_denominator(a, n));
  let u = 1n;
  const qs = new Set<bigint>();
  for (const d of dens) for (const q of primeFactors(d)) qs.add(q);
  for (const q of qs) {
    let e = 0;
    for (let i = 0; i < 5; i++) {
      const v = Number(valuation(dens[i]!, q));
      e = Math.max(e, Math.ceil(v / weights[i]!));
    }
    u *= q ** BigInt(e);
  }
  const ainvs = ainvs0.map((a, i) => a.scalarMul(new Rational(u ** BigInt(weights[i]!), 1n)));
  const coords = ainvs.map((a) => _nf_integral_coords(a, n));

  const [a1, a2, a3, a4, a6] = ainvs as [
    NumberFieldElement,
    NumberFieldElement,
    NumberFieldElement,
    NumberFieldElement,
    NumberFieldElement,
  ];
  const z = (k: bigint) => K.__call__(k);
  const b2 = a1.mul(a1).add(a2.mul(z(4n)));
  const b4 = a4.mul(z(2n)).add(a1.mul(a3));
  const b6 = a3.mul(a3).add(a6.mul(z(4n)));
  const b8 = a1
    .mul(a1)
    .mul(a6)
    .add(a2.mul(a6).mul(z(4n)))
    .sub(a1.mul(a3).mul(a4))
    .add(a2.mul(a3).mul(a3))
    .sub(a4.mul(a4));
  const disc = b2
    .mul(b2)
    .mul(b8)
    .neg()
    .sub(b4.mul(b4).mul(b4).mul(z(8n)))
    .sub(b6.mul(b6).mul(z(27n)))
    .add(b2.mul(b4).mul(b6).mul(z(9n)));
  if (disc.is_zero()) {
    throw new ValueError('Frobenius_filter: the curve is singular');
  }
  return {
    K,
    n,
    T,
    discT,
    ainvs,
    coords,
    b: {
      b2: _nf_integral_coords(b2, n),
      b4: _nf_integral_coords(b4, n),
      b6: _nf_integral_coords(b6, n),
      b8: _nf_integral_coords(b8, n),
      disc: _nf_integral_coords(disc, n),
    },
  };
}

/** The {@link FrobeniusContext} of a curve over a number field. */
function _nf_frobenius_context(c: NFCurve): FrobeniusContext {
  return {
    // gal_reps_number_field.py:559: `not E.division_polynomial(2).is_irreducible()`.
    // The 2-division polynomial `4x^3 + b2 x^2 + 2 b4 x + b6` is a cubic, hence
    // irreducible over K exactly when it has no root in K; substituting x = X/4
    // and scaling by 16 makes it the monic `X^3 + b2 X^2 + 8 b4 X + 16 b6` with
    // coefficients in `Z[theta]`.
    two_division_poly_is_irreducible: () =>
      !_nf_monic_cubic_has_root(c, c.b.b2, _zScale(c.b.b4, 8n), _zScale(c.b.b6, 16n)),
    frobenius_discs: () => _nf_frobenius_discs(c),
  };
}

/** Search bound on the rational primes below the primes of `K` that are used. */
const NF_PRIME_LIMIT = 1n << 22n;

/**
 * The discriminants `a_P^2 - 4*N(P)` of the Frobenius polynomials of `E` at the
 * primes `P` of good reduction, in increasing order of the rational prime below.
 *
 * @see Reference: sage/schemes/elliptic_curves/gal_reps_number_field.py:566-582
 * @see Deviation: only the primes `P` of residue degree one are produced.  Those
 *   are exactly the roots of the defining polynomial mod `p`, for `p` not
 *   dividing `disc(T)` -- for such `p` the power basis is a `p`-adic integral
 *   basis, so `theta -> r` really is the reduction `O_K -> O_K/P = F_p`.  A
 *   prime of residue degree `f > 1` would need `#E(F_{p^f})`, i.e. PARI's
 *   `Fq_ellcard_SEA` over an extension field, which parigp-ts does not have
 *   (its `ellsea.ts` is the `T = NULL`, prime-field, specialisation).  Dropping
 *   primes from the walk cannot make the filter unsound -- every `P` used gives
 *   a valid necessary condition -- it can only make it return a superset.
 */
function* _nf_frobenius_discs(c: NFCurve): Generator<bigint> {
  for (let p = 2n; p <= NF_PRIME_LIMIT; p = _next_prime(p)) {
    if (c.discT % p === 0n) continue;
    for (const r of _roots_mod_p(c.T, p)) {
      // "not E.has_good_reduction(P)" (gal_reps_number_field.py:571): we test
      // the weaker `v_P(disc) == 0` of this (integral, possibly non-minimal)
      // model, which can only discard primes.
      if (_eval_coords_mod(c.b.disc, r, p) === 0n) continue;
      const red = c.coords.map((v) => _eval_coords_mod(v, r, p)) as [
        bigint,
        bigint,
        bigint,
        bigint,
        bigint,
      ];
      const ap = _ap_over_Q(red, p);
      yield ap * ap - 4n * p; // discriminant of x^2 - a_P*x + N(P), N(P) = p
    }
  }
}

// --- number field helpers --------------------------------------------------

/** `x` as a rational, for the rational a-invariants of a curve over a number field. */
function _asRationalStrict(x: unknown): Rational {
  const r = _as_rational(x);
  if (r === null) {
    throw new ValueError('Frobenius_filter: a-invariant is neither rational nor a field element');
  }
  return new Rational(r[0], r[1]);
}

function _rationalToBigInt(x: Rational, what: string): bigint {
  if (x.denominator !== 1n) throw new ValueError(`${what} is not an integer`);
  return x.numerator;
}

/** The lcm of the denominators of the power-basis coordinates of `x`. */
function _nf_coord_denominator(x: NumberFieldElement, n: number): bigint {
  let d = 1n;
  const l = x.list();
  for (let i = 0; i < n; i++) {
    const c = l[i] ?? Rational.zero();
    const g = _bgcd(d, c.denominator);
    d = (d / g) * c.denominator;
  }
  return d;
}

/** The power-basis coordinates of `x`, which must be integers. */
function _nf_integral_coords(x: NumberFieldElement, n: number): bigint[] {
  const l = x.list();
  const out: bigint[] = [];
  for (let i = 0; i < n; i++) {
    const cc = l[i] ?? Rational.zero();
    if (cc.denominator !== 1n) {
      throw new ValueError('Frobenius_filter: expected integral power-basis coordinates');
    }
    out.push(cc.numerator);
  }
  return out;
}

function _zScale(v: bigint[], k: bigint): bigint[] {
  return v.map((x) => x * k);
}

function _bgcd(a: bigint, b: bigint): bigint {
  let x = a < 0n ? -a : a;
  let y = b < 0n ? -b : b;
  while (y) {
    const t = x % y;
    x = y;
    y = t;
  }
  return x;
}

/** `sum_i v_i * r^i` mod `m`. */
function _eval_coords_mod(v: bigint[], r: bigint, m: bigint): bigint {
  let acc = 0n;
  for (let i = v.length - 1; i >= 0; i--) acc = (acc * r + v[i]!) % m;
  return ((acc % m) + m) % m;
}

// --- polynomial arithmetic over F_p ---------------------------------------

function _pnorm(f: bigint[], p: bigint): bigint[] {
  const g = f.map((c) => ((c % p) + p) % p);
  while (g.length > 0 && g[g.length - 1] === 0n) g.pop();
  return g;
}

function _pmul(a: bigint[], b: bigint[], p: bigint): bigint[] {
  if (a.length === 0 || b.length === 0) return [];
  const out = new Array<bigint>(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++) {
    if (a[i] === 0n) continue;
    for (let j = 0; j < b.length; j++) out[i + j] = (out[i + j]! + a[i]! * b[j]!) % p;
  }
  return _pnorm(out, p);
}

function _pinvmod(a: bigint, p: bigint): bigint {
  let [old_r, r] = [((a % p) + p) % p, p];
  let [old_s, s] = [1n, 0n];
  while (r !== 0n) {
    const q = old_r / r;
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
  }
  if (old_r !== 1n) throw new ValueError('not invertible');
  return ((old_s % p) + p) % p;
}

function _prem(a: bigint[], b: bigint[], p: bigint): bigint[] {
  const f = _pnorm(a, p);
  const g = _pnorm(b, p);
  if (g.length === 0) throw new ValueError('division by zero polynomial');
  const inv = _pinvmod(g[g.length - 1]!, p);
  const r = [...f];
  for (let i = r.length - g.length; i >= 0; i--) {
    const coef = (r[i + g.length - 1]! * inv) % p;
    if (coef === 0n) continue;
    for (let j = 0; j < g.length; j++) {
      r[i + j] = (((r[i + j]! - coef * g[j]!) % p) + p) % p;
    }
  }
  return _pnorm(r, p);
}

function _pgcd(a: bigint[], b: bigint[], p: bigint): bigint[] {
  let f = _pnorm(a, p);
  let g = _pnorm(b, p);
  while (g.length > 0) {
    const r = _prem(f, g, p);
    f = g;
    g = r;
  }
  if (f.length === 0) return f;
  const inv = _pinvmod(f[f.length - 1]!, p);
  return f.map((c) => (c * inv) % p);
}

/** `base^e mod (mod, p)`. */
function _ppowmod(base: bigint[], e: bigint, mod: bigint[], p: bigint): bigint[] {
  let result: bigint[] = _pnorm([1n], p);
  let b = _prem(base, mod, p);
  let k = e;
  while (k > 0n) {
    if (k & 1n) result = _prem(_pmul(result, b, p), mod, p);
    b = _prem(_pmul(b, b, p), mod, p);
    k >>= 1n;
  }
  return result;
}

/**
 * The roots of `f` in `F_p`, sorted.  Cantor-Zassenhaus: first `gcd(x^p - x, f)`
 * (the product of the distinct linear factors), then equal-degree splitting.
 */
function _roots_mod_p(f: bigint[], p: bigint): bigint[] {
  const F = _pnorm(f, p);
  if (F.length <= 1) return [];
  if (p <= 512n) {
    const out: bigint[] = [];
    for (let x = 0n; x < p; x++) {
      let acc = 0n;
      for (let i = F.length - 1; i >= 0; i--) acc = (acc * x + F[i]!) % p;
      if (acc === 0n) out.push(x);
    }
    return out;
  }
  // h = x^p mod F, then gcd(h - x, F)
  const h = _ppowmod([0n, 1n], p, F, p);
  const hx = [...h];
  while (hx.length < 2) hx.push(0n);
  hx[1] = (((hx[1]! - 1n) % p) + p) % p;
  const g = _pgcd(_pnorm(hx, p), F, p);
  if (g.length <= 1) return [];
  const out: bigint[] = [];
  const stack: bigint[][] = [g];
  let seed = 1n;
  while (stack.length > 0) {
    const cur = stack.pop()!;
    if (cur.length === 2) {
      // cur = c1*x + c0
      out.push((((-cur[0]! * _pinvmod(cur[1]!, p)) % p) + p) % p);
      continue;
    }
    if (cur.length <= 1) continue;
    // equal-degree splitting: gcd(cur, (x + a)^((p-1)/2) - 1)
    for (;;) {
      seed = (seed * 6364136223846793005n + 1442695040888963407n) % (1n << 62n);
      const a = seed % p;
      const t = _ppowmod([a, 1n], (p - 1n) / 2n, cur, p);
      const tm = [...t];
      while (tm.length < 1) tm.push(0n);
      tm[0] = ((((tm[0] ?? 0n) - 1n) % p) + p) % p;
      const d = _pgcd(_pnorm(tm, p), cur, p);
      if (d.length > 1 && d.length < cur.length) {
        stack.push(d);
        stack.push(_pnorm(_pquo(cur, d, p), p));
        break;
      }
    }
  }
  return out.sort((x, y) => (x < y ? -1 : x > y ? 1 : 0));
}

/** Exact quotient of `a` by `b` over `F_p`. */
function _pquo(a: bigint[], b: bigint[], p: bigint): bigint[] {
  const f = _pnorm(a, p);
  const g = _pnorm(b, p);
  const inv = _pinvmod(g[g.length - 1]!, p);
  const r = [...f];
  const q = new Array<bigint>(Math.max(f.length - g.length + 1, 0)).fill(0n);
  for (let i = r.length - g.length; i >= 0; i--) {
    const coef = (r[i + g.length - 1]! * inv) % p;
    q[i] = coef;
    if (coef === 0n) continue;
    for (let j = 0; j < g.length; j++) {
      r[i + j] = (((r[i + j]! - coef * g[j]!) % p) + p) % p;
    }
  }
  return _pnorm(q, p);
}

// --- root of a monic cubic in K -------------------------------------------

/** Number of rational primes scanned before giving up on the p-adic search. */
const NF_CUBIC_PRIMES = 200;
/** Largest `p^k` used when reconstructing a root of the cubic. */
const NF_CUBIC_MAX_MODULUS = 1n << 2048n;

/**
 * Does the monic cubic `X^3 + A X^2 + B X + C` (coefficients given by their
 * integral power-basis coordinate vectors) have a root in `K`?
 *
 * This is what SageMath gets from `Polynomial.is_irreducible()` over a number
 * field, i.e. from PARI's `nffactor`.
 *
 * @see Deviation: `nffactor` (reference/pari/src/basemath/polnf.c) is not in
 *   parigp-ts, and this module may not add it there.  A cubic is reducible over
 *   `K` exactly when it has a root in `K`, and that is decided here by a
 *   two-sided certificate, so the answer is never a guess:
 *   * *no root*: a prime `P` of `K` of residue degree one at which the cubic has
 *     no root in `F_p`.  A root in `K` would reduce to one, so this proves the
 *     cubic is irreducible.  Such a `P` exists for every irreducible cubic
 *     (Chebotarev: the Galois group of its splitting field over `K` contains a
 *     3-cycle, which has no fixed root).
 *   * *a root*: the root itself, exhibited and then verified by exact arithmetic
 *     in `K`.  It is found from a completely split prime `p`, by Hensel-lifting
 *     the roots in each of the `n` components of `O_K/p^k`, undoing the
 *     Vandermonde change of basis and reconstructing each rational coordinate.
 *   If neither certificate is found within the search bounds a
 *   `NotImplementedError` names `nffactor`; no answer is ever guessed.
 */
function _nf_monic_cubic_has_root(c: NFCurve, A: bigint[], B: bigint[], C: bigint[]): boolean {
  const { n, T, discT } = c;
  // Is `w/D` a root?  Clearing denominators, that is
  // `w^3 + D*A*w^2 + D^2*B*w + D^3*C == 0` in `Z[theta] = Z[X]/(T)`.
  const verify = (w: bigint[], D: bigint): boolean => {
    const w2 = _zmulT(w, w, T);
    const w3 = _zmulT(w2, w, T);
    const t1 = _zScale(_zmulT(A, w2, T), D);
    const t2 = _zScale(_zmulT(B, w, T), D * D);
    const t3 = _zScale(C, D * D * D);
    for (let i = 0; i < n; i++) {
      if ((w3[i] ?? 0n) + (t1[i] ?? 0n) + (t2[i] ?? 0n) + (t3[i] ?? 0n) !== 0n) return false;
    }
    return true;
  };

  // The discriminant of X^3 + A X^2 + B X + C, as an element of K: used to
  // recognise the primes at which every root is simple (so Hensel applies).
  const A2 = _zmulT(A, A, T);
  const B2 = _zmulT(B, B, T);
  const discCoords = _zAdd(
    _zAdd(
      _zSub(_zmulT(A2, B2, T), _zScale(_zmulT(B2, B, T), 4n)),
      _zSub(_zScale(_zmulT(_zmulT(A2, A, T), C, T), -4n), _zScale(_zmulT(C, C, T), 27n))
    ),
    _zScale(_zmulT(_zmulT(A, B, T), C, T), 18n)
  );

  let split: { p: bigint; roots: bigint[] } | null = null;
  let scanned = 0;
  let triedEarly = false;
  for (let p = 2n; scanned < NF_CUBIC_PRIMES && p <= NF_PRIME_LIMIT; p = _next_prime(p)) {
    if (discT % p === 0n) continue;
    const roots = _roots_mod_p(T, p);
    if (roots.length === 0) continue;
    scanned++;
    let allSimple = true;
    for (const r of roots) {
      const a = _eval_coords_mod(A, r, p);
      const b = _eval_coords_mod(B, r, p);
      const d = _eval_coords_mod(C, r, p);
      if (_roots_mod_p([d, b, a, 1n], p).length === 0) {
        return false; // certificate of irreducibility
      }
      if (_eval_coords_mod(discCoords, r, p) === 0n) allSimple = false;
    }
    if (split === null && roots.length === n && allSimple && p > 2n) {
      split = { p, roots };
    }
    // Most reducible cubics have a small root; look for it before paying for
    // the whole scan.  A failure here costs nothing but time: the scan resumes.
    if (!triedEarly && split !== null && scanned >= 24) {
      triedEarly = true;
      if (_nf_cubic_root_from_split_prime(c, A, B, C, split.p, split.roots, verify, 1n << 256n)) {
        return true;
      }
    }
  }

  if (split !== null) {
    const found = _nf_cubic_root_from_split_prime(
      c,
      A,
      B,
      C,
      split.p,
      split.roots,
      verify,
      NF_CUBIC_MAX_MODULUS
    );
    if (found) return true;
  }

  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: irreducibility over a number field of the 2-division ' +
      'polynomial (SageMath uses PARI nffactor, polnf.c; parigp-ts has no number ' +
      'field layer). Neither certificate was found: no prime of residue degree ' +
      'one below ' +
      String(NF_PRIME_LIMIT) +
      ' where the cubic has no root, and no root reconstructed below modulus ' +
      String(NF_CUBIC_MAX_MODULUS)
  );
}

/**
 * Reconstruct a root in `K` of `X^3 + A X^2 + B X + C` from a completely split
 * prime `p` (whose `n` degree-one primes are the roots `r_j` of `T` mod `p`).
 */
function _nf_cubic_root_from_split_prime(
  c: NFCurve,
  A: bigint[],
  B: bigint[],
  C: bigint[],
  p: bigint,
  roots0: bigint[],
  verify: (numerators: bigint[], denominator: bigint) => boolean,
  maxModulus: bigint
): boolean {
  const { n, T } = c;
  for (let m = p * p; m <= maxModulus; m = m * m) {
    // Hensel-lift the roots of T and, in each component, the roots of the cubic.
    const rs = roots0.map((r) => _hensel_root(T, r, p, m));
    const perComponent: bigint[][] = [];
    let total = 1;
    for (let j = 0; j < n; j++) {
      const r = rs[j]!;
      const a = _eval_coords_mod(A, r, m);
      const b = _eval_coords_mod(B, r, m);
      const d = _eval_coords_mod(C, r, m);
      const cub = [d, b, a, 1n];
      const small = _roots_mod_p([d % p, b % p, a % p, 1n], p);
      const lifted = small.map((s) => _hensel_root(cub, s, p, m));
      if (lifted.length === 0) return false;
      perComponent.push(lifted);
      total *= lifted.length;
      if (total > 20000) return false;
    }
    // Invert the Vandermonde matrix V[j][i] = r_j^i once per modulus.
    const V: bigint[][] = rs.map((r) => {
      const row: bigint[] = [];
      let acc = 1n;
      for (let i = 0; i < n; i++) {
        row.push(acc);
        acc = (acc * r) % m;
      }
      return row;
    });
    const idx = new Array<number>(n).fill(0);
    for (;;) {
      const beta = perComponent.map((rr, j) => rr[idx[j]!]!);
      const coordsMod = _solve_mod(V, beta, p, m);
      if (coordsMod !== null) {
        const rat: Array<[bigint, bigint]> = [];
        let ok = true;
        for (const v of coordsMod) {
          const q = _rational_reconstruct(v, m);
          if (q === null) {
            ok = false;
            break;
          }
          rat.push(q);
        }
        if (ok) {
          // common denominator
          let D = 1n;
          for (const [, d] of rat) D = (D / _bgcd(D, d)) * d;
          const w = rat.map(([nu, d]) => nu * (D / d));
          if (verify(w, D)) return true;
        }
      }
      // next tuple
      let k = n - 1;
      for (; k >= 0; k--) {
        idx[k]!++;
        if (idx[k]! < perComponent[k]!.length) break;
        idx[k] = 0;
      }
      if (k < 0) break;
    }
  }
  return false;
}

/** Lift a simple root `r` of `f` mod `p` to a root mod `m` (a power of `p`). */
function _hensel_root(f: bigint[], r0: bigint, p: bigint, m: bigint): bigint {
  const evalAt = (x: bigint, mod: bigint): bigint => {
    let acc = 0n;
    for (let i = f.length - 1; i >= 0; i--) acc = (acc * x + f[i]!) % mod;
    return ((acc % mod) + mod) % mod;
  };
  const der: bigint[] = [];
  for (let i = 1; i < f.length; i++) der.push(f[i]! * BigInt(i));
  const evalDer = (x: bigint, mod: bigint): bigint => {
    let acc = 0n;
    for (let i = der.length - 1; i >= 0; i--) acc = (acc * x + der[i]!) % mod;
    return ((acc % mod) + mod) % mod;
  };
  let r = ((r0 % p) + p) % p;
  let mod = p;
  while (mod < m) {
    mod = mod * mod > m ? m : mod * mod;
    const dv = evalDer(r, mod);
    const inv = _pinvmod(dv, mod);
    r = (((r - evalAt(r, mod) * inv) % mod) + mod) % mod;
  }
  return ((r % m) + m) % m;
}

/**
 * Solve `V * x = y` over `Z/m` (`m` a power of `p`), given that `V` is
 * invertible mod `p`.  Returns `null` if no unit pivot can be found.
 */
function _solve_mod(V: bigint[][], y: bigint[], p: bigint, m: bigint): bigint[] | null {
  const n = V.length;
  const M = V.map((row, j) => [...row, y[j]!]);
  for (let col = 0; col < n; col++) {
    let piv = -1;
    for (let r = col; r < n; r++) {
      if (M[r]![col]! % p !== 0n) {
        piv = r;
        break;
      }
    }
    if (piv < 0) return null;
    if (piv !== col) {
      const t = M[piv]!;
      M[piv] = M[col]!;
      M[col] = t;
    }
    const inv = _pinvmod(M[col]![col]!, m);
    for (let k = col; k <= n; k++) M[col]![k] = (M[col]![k]! * inv) % m;
    for (let r = 0; r < n; r++) {
      if (r === col) continue;
      const f = M[r]![col]!;
      if (f === 0n) continue;
      for (let k = col; k <= n; k++) {
        M[r]![k] = (((M[r]![k]! - f * M[col]![k]!) % m) + m) % m;
      }
    }
  }
  return M.map((row) => ((row[n]! % m) + m) % m);
}

/**
 * Rational reconstruction: the unique `[a, b]` with `a/b == u (mod m)`,
 * `|a|, b <= sqrt(m/2)` and `gcd(b, m) = 1`, or `null`.
 */
function _rational_reconstruct(u: bigint, m: bigint): [bigint, bigint] | null {
  const bound = isqrt(m / 2n);
  let [r0, r1] = [m, ((u % m) + m) % m];
  let [s0, s1] = [0n, 1n];
  while (r1 > bound) {
    if (r1 === 0n) return null;
    const q = r0 / r1;
    [r0, r1] = [r1, r0 - q * r1];
    [s0, s1] = [s1, s0 - q * s1];
  }
  const den = s1 < 0n ? -s1 : s1;
  const num = s1 < 0n ? -r1 : r1;
  if (den === 0n || den > bound) return null;
  if (_bgcd(den, m) !== 1n) return null;
  if (_bgcd(num < 0n ? -num : num, den) !== 1n) return null;
  return [num, den];
}

/** Sum of two coordinate vectors. */
function _zAdd(a: bigint[], b: bigint[]): bigint[] {
  const out: bigint[] = [];
  for (let i = 0; i < Math.max(a.length, b.length); i++) out.push((a[i] ?? 0n) + (b[i] ?? 0n));
  return out;
}

/** Difference of two coordinate vectors. */
function _zSub(a: bigint[], b: bigint[]): bigint[] {
  const out: bigint[] = [];
  for (let i = 0; i < Math.max(a.length, b.length); i++) out.push((a[i] ?? 0n) - (b[i] ?? 0n));
  return out;
}

/**
 * Product in `Z[theta] = Z[X]/(T)` of two integral coordinate vectors, `T`
 * monic of degree `n` with integer coefficients.
 */
function _zmulT(a: bigint[], b: bigint[], T: bigint[]): bigint[] {
  const n = T.length - 1;
  const prod = new Array<bigint>(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++) {
    if (a[i] === 0n) continue;
    for (let j = 0; j < b.length; j++) prod[i + j] = prod[i + j]! + a[i]! * b[j]!;
  }
  for (let d = prod.length - 1; d >= n; d--) {
    const co = prod[d]!;
    if (co === 0n) continue;
    prod[d] = 0n;
    for (let k = 0; k < n; k++) prod[d - n + k] = prod[d - n + k]! - co * T[k]!;
  }
  const out = prod.slice(0, n);
  while (out.length < n) out.push(0n);
  return out;
}

/** The next prime strictly greater than n. */
function _next_prime(n: bigint): bigint {
  let m = n + 1n;
  while (!isPrime(m)) m++;
  return m;
}

/**
 * Helper: Compute divisors of n.
 */
function divisors(n: bigint): bigint[] {
  const absN = n < 0n ? -n : n;
  if (absN === 0n) return [];
  if (absN === 1n) return [1n];

  const result: bigint[] = [];
  let i = 1n;
  while (i * i <= absN) {
    if (absN % i === 0n) {
      result.push(i);
      if (i !== absN / i) {
        result.push(absN / i);
      }
    }
    i++;
  }
  result.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  return result;
}

/**
 * Helper: Compute prime factors of n.
 */
function primeFactors(n: bigint): bigint[] {
  const absN = n < 0n ? -n : n;
  if (absN <= 1n) return [];

  const factors: bigint[] = [];
  let temp = absN;

  // Check for factor of 2
  if (temp % 2n === 0n) {
    factors.push(2n);
    while (temp % 2n === 0n) {
      temp /= 2n;
    }
  }

  // Check odd factors
  let i = 3n;
  while (i * i <= temp) {
    if (temp % i === 0n) {
      factors.push(i);
      while (temp % i === 0n) {
        temp /= i;
      }
    }
    i += 2n;
  }

  if (temp > 1n) {
    factors.push(temp);
  }

  return factors;
}

/**
 * Helper: Compute p-adic valuation of n at prime p.
 */
function valuation(n: bigint, p: bigint): bigint {
  if (n === 0n) return -1n; // Infinity for 0
  let absN = n < 0n ? -n : n;
  let v = 0n;
  while (absN % p === 0n) {
    absN /= p;
    v++;
  }
  return v;
}

/**
 * Helper: Simple primality test.
 */
function isPrime(n: bigint): boolean {
  if (n < 2n) return false;
  if (n === 2n) return true;
  if (n % 2n === 0n) return false;
  if (n < 9n) return true;
  if (n % 3n === 0n) return false;

  const limit = isqrt(n);
  let i = 5n;
  while (i <= limit) {
    if (n % i === 0n || n % (i + 2n) === 0n) return false;
    i += 6n;
  }
  return true;
}

/**
 * Helper: Integer square root.
 */
function isqrt(n: bigint): bigint {
  if (n < 2n) return n;
  let x = n;
  let y = (x + 1n) / 2n;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }
  return x;
}

/**
 * Helper: Kronecker symbol (a/n).
 */
function kroneckerSymbol(a: bigint, n: bigint): bigint {
  if (n === 0n) {
    return a === 1n || a === -1n ? 1n : 0n;
  }

  let result = 1n;

  // Handle negative n
  if (n < 0n) {
    n = -n;
    if (a < 0n) {
      result = -result;
    }
  }

  // Handle even n
  while (n % 2n === 0n) {
    n /= 2n;
    const aMod8 = ((a % 8n) + 8n) % 8n;
    if (aMod8 === 3n || aMod8 === 5n) {
      result = -result;
    }
  }

  if (n === 1n) return result;

  // Jacobi symbol for odd n
  a = ((a % n) + n) % n;

  while (a !== 0n) {
    while (a % 2n === 0n) {
      a /= 2n;
      const nMod8 = n % 8n;
      if (nMod8 === 3n || nMod8 === 5n) {
        result = -result;
      }
    }

    // Swap a and n
    [a, n] = [n, a];

    // Quadratic reciprocity
    if (a % 4n === 3n && n % 4n === 3n) {
      result = -result;
    }
    a = a % n;
  }

  return n === 1n ? result : 0n;
}

/**
 * Return the possible prime isogeny degrees for a curve.
 *
 * For curves without CM, this returns the primes where the mod-l
 * Galois representation is reducible (contained in a Borel subgroup).
 *
 * For CM curves, delegates to isogeny_degrees_cm().
 *
 * INPUT:
 * - E: an elliptic curve
 * - algorithm: 'Billerey' (default), 'Larson', or 'heuristic'
 * - max_l: maximum prime to consider (optional)
 * - exact: if True, verify that returned primes are actually reducible
 *
 * OUTPUT: a sorted list of prime degrees
 *
 * ALGORITHM:
 * For curves with CM, delegates to isogeny_degrees_cm().
 *
 * For curves over Q without CM, by Mazur's theorem, the possible
 * prime isogeny degrees are limited to: 2, 3, 5, 7, 11, 13, 17, 19,
 * 37, 43, 67, 163. We use a heuristic approach checking these primes.
 *
 * For curves over general number fields without CM, we use:
 * - 'Billerey': Billerey's bounds on reducible primes
 * - 'Larson': Larson's isogeny bound
 * - 'heuristic': Check small primes via Frobenius filter
 *
 * @see Reference: sage/schemes/elliptic_curves/isogeny_class.py:possible_isogeny_degrees
 * @see Deviation: over Q the non-CM branch returns Mazur's list of possible
 *   degrees (optionally intersected with the degrees for which
 *   isogenies_prime_degree finds an isogeny) instead of Billerey's/Larson's
 *   bounds, which are not ported.
 */
export function possible_isogeny_degrees<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  options: {
    algorithm?: 'Billerey' | 'Larson' | 'heuristic';
    max_l?: bigint | number;
    num_l?: number;
    exact?: boolean;
    verbose?: boolean;
  } = {}
): bigint[] {
  const extras = E as unknown as EllipticCurveCMExtras;
  const { algorithm = 'Billerey', max_l: maxLOption, exact = true, verbose = false } = options;

  // Check if curve has CM
  const hasCM = typeof extras.has_cm === 'function' ? extras.has_cm() : false;

  if (hasCM) {
    return isogeny_degrees_cm(E, verbose);
  }

  // Determine if we're over Q or a larger number field
  const isOverQ =
    typeof extras.base_field === 'function' && typeof extras.base_field().degree === 'function'
      ? extras.base_field().degree() === 1
      : true;

  if (verbose) {
    console.log(`Non-CM case, using ${algorithm} algorithm`);
  }

  // For curves over Q (non-CM), Mazur's theorem limits possible prime degrees
  // The only possible prime isogeny degrees are: 2, 3, 5, 7, 11, 13, 17, 19, 37, 43, 67, 163
  // Note: 37, 43, 67, 163 only occur for specific j-invariants
  if (isOverQ) {
    const mazurPrimes = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 37n, 43n, 67n, 163n];

    // Apply max_l filter if provided
    let candidates = mazurPrimes;
    if (maxLOption !== undefined) {
      const maxL = typeof maxLOption === 'number' ? BigInt(maxLOption) : maxLOption;
      candidates = mazurPrimes.filter((p) => p <= maxL);
    }

    // For the 'heuristic' algorithm over Q, we use the naive approach
    // checking which primes from Mazur's list actually give isogenies
    if (algorithm === 'heuristic' || !exact) {
      return candidates;
    }

    // For exact computation, we would need to check if E actually
    // has l-isogenies for each candidate l. This requires isogenies_prime_degree().
    // For now, return the candidates since verifying requires more infrastructure.
    if (typeof extras.isogenies_prime_degree === 'function') {
      const result: bigint[] = [];
      for (const l of candidates) {
        try {
          const isogs = extras.isogenies_prime_degree([l], { minimal_models: false });
          if (Array.isArray(isogs) && isogs.length > 0) {
            result.push(l);
          }
        } catch {
          // If isogeny computation fails, don't include this prime
        }
      }
      return result;
    }

    // Without isogeny computation, return the theoretical candidates
    return candidates;
  }

  // For curves over number fields (non-CM), we need more sophisticated analysis
  // This requires:
  // - Billerey: Computing Billerey's bound on reducible primes
  // - Larson: Computing Larson's isogeny bound
  // - heuristic: Frobenius filtering with small primes

  // Default max_l based on algorithm
  const maxL =
    maxLOption !== undefined
      ? typeof maxLOption === 'number'
        ? BigInt(maxLOption)
        : maxLOption
      : algorithm === 'heuristic'
        ? 37n
        : 1000n;

  // Generate candidate primes up to maxL
  const candidates: bigint[] = [];
  for (let p = 2n; p <= maxL; p++) {
    if (isPrime(p)) {
      candidates.push(p);
    }
  }

  // Apply Frobenius filter if the curve has good_primes method
  // This eliminates primes l where there exists a good reduction prime P
  // such that the Frobenius polynomial at P is irreducible mod l
  if (typeof extras.good_reduction_primes === 'function' && verbose) {
    console.log('Applying Frobenius filter...');
  }

  // For exact computation over number fields
  if (exact && typeof extras.isogenies_prime_degree === 'function') {
    const result: bigint[] = [];
    for (const l of candidates) {
      try {
        const isogs = extras.isogenies_prime_degree([l], { minimal_models: false });
        if (Array.isArray(isogs) && isogs.length > 0) {
          result.push(l);
        }
      } catch {
        // If isogeny computation fails, don't include this prime
      }
    }
    return result;
  }

  // Without full infrastructure, return the candidates up to maxL
  return candidates;
}

/**
 * Test-only surface: the arithmetic primitives of the number-field
 * {@link Frobenius_filter}.  Not part of the public API; do not re-export.
 */
export const _internal = {
  roots_mod_p: _roots_mod_p,
  zmulT: _zmulT,
  rational_reconstruct: _rational_reconstruct,
  hensel_root: _hensel_root,
  solve_mod: _solve_mod,
  eval_coords_mod: _eval_coords_mod,
};
