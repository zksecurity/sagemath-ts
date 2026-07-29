/**
 * @module sage/modules/free_module
 * @description Free modules over commutative rings
 * @see Reference: sage/modules/free_module.py
 *
 * Sage supports computation with free modules over an arbitrary commutative ring.
 * Nontrivial functionality is available over ZZ, fields, and some principal
 * ideal domains (e.g. QQ[x] and rings of integers of number fields).
 */

import { ArithmeticError, NotImplementedError, ValueError, ZeroDivisionError } from '../errors.js';
import {
  IntegerMatrixFromEntries,
  hermite_normal_form,
  saturation as matrixSaturation,
  smith_form_integer,
} from '../matrix/matrix_integer.js';
import { Rational } from '../rings/rational.js';
import {
  type FreeModuleElement,
  FreeModuleElementDense,
  FreeModuleElementSparse,
  type FreeModuleParent,
  type RingLike,
} from './free_module_element.js';

/**
 * Options for creating a free module.
 */
export interface FreeModuleOptions {
  sparse?: boolean;
  innerProductMatrix?: unknown;
  withBasis?: 'standard' | null;
}

/**
 * Create a free module over the given base ring.
 *
 * @param baseRing - A commutative ring
 * @param rank - The rank (dimension) of the module
 * @param options - Additional options
 * @returns A free module of the given rank over the base ring
 *
 * @example
 * ```typescript
 * // Create ZZ^3
 * const M = FreeModule(ZZ, 3);
 *
 * // Create a sparse module
 * const S = FreeModule(ZZ, 100, { sparse: true });
 *
 * // Create with inner product matrix
 * const Q = FreeModule(ZZ, 2, { innerProductMatrix: [[1, 0], [0, -1]] });
 * ```
 *
 * @see Reference: sage/modules/free_module.py:FreeModule
 */
export function FreeModule(
  baseRing: RingLike,
  rank: number | bigint,
  options?: FreeModuleOptions
): FreeModuleGeneric {
  const rankNum = typeof rank === 'bigint' ? Number(rank) : rank;
  if (rankNum < 0) {
    throw new ValueError(`rank (=${rankNum}) must be nonnegative`);
  }

  const sparse = options?.sparse ?? false;
  const innerProductMatrix = options?.innerProductMatrix;

  // Determine what type of module to create based on base ring
  // For now, we use the generic ambient module
  if (isField(baseRing)) {
    return new FreeModuleAmbientField(baseRing, rankNum, sparse, innerProductMatrix);
  }
  if (isPID(baseRing)) {
    return new FreeModuleAmbientPID(baseRing, rankNum, sparse, innerProductMatrix);
  }

  return new FreeModuleAmbient(baseRing, rankNum, sparse, innerProductMatrix);
}

/**
 * Create a vector space over a field.
 *
 * @param field - A field
 * @param dimension - The dimension of the vector space
 * @param options - Additional options
 * @returns A vector space of the given dimension
 *
 * @see Reference: sage/modules/free_module.py:VectorSpace
 */
export function VectorSpace(
  field: RingLike,
  dimension: number | bigint,
  options?: FreeModuleOptions
): FreeModuleField {
  if (!isField(field)) {
    throw new TypeError(`Argument K (= ${field}) must be a field.`);
  }

  const dimensionNum = typeof dimension === 'bigint' ? Number(dimension) : dimension;
  if (dimensionNum < 0) {
    throw new ValueError(`dimension (=${dimensionNum}) must be nonnegative`);
  }

  const sparse = options?.sparse ?? false;
  const innerProductMatrix = options?.innerProductMatrix;

  return new FreeModuleAmbientField(field, dimensionNum, sparse, innerProductMatrix);
}

/**
 * Return the span of the vectors in gens using scalars from baseRing.
 *
 * @param gens - List of vectors or lists of ring elements
 * @param baseRing - Optional base ring (default: inferred from gens)
 * @param options - Additional options
 * @returns The R-span of the generators
 *
 * @example
 * ```typescript
 * // Create span over QQ
 * const V = span([[1, 2, 5], [2, 2, 2]], QQ);
 *
 * // Create span over ZZ
 * const M = span([[1, 2, 3], [4, 5, 6]], ZZ);
 * ```
 *
 * @see Reference: sage/modules/free_module.py:span
 */
export function span(
  gens: unknown[][],
  baseRing?: RingLike,
  options?: { check?: boolean; alreadyEchelonized?: boolean }
): FreeModuleGeneric {
  if (gens.length === 0) {
    if (!baseRing) {
      throw new ValueError('base_ring must be specified for empty generators');
    }
    return FreeModule(baseRing, 0);
  }

  // Infer degree from first generator
  const degree = gens[0]!.length;

  // Infer base ring if not provided
  const ring = baseRing ?? inferRing(gens[0]!);

  // Create the ambient module
  const ambient = FreeModule(ring, degree);

  // Create the span
  const vectors = gens.map((g) => ambient.createElement(g));

  return ambient.span(vectors, ring, options);
}

// ============================================================================
// Helper functions
// ============================================================================

/**
 * Check if a ring is a field.
 */
function isField(ring: RingLike): boolean {
  if ('is_field' in ring && typeof ring.is_field === 'function') {
    return ring.is_field();
  }
  return false;
}

/**
 * Check if a ring is a PID.
 */
function isPID(ring: RingLike): boolean {
  // Check for common PIDs
  if ('is_principal_ideal_domain' in ring) {
    return (ring as { is_principal_ideal_domain: () => boolean }).is_principal_ideal_domain();
  }
  // ZZ is a PID
  if (ring.toString?.() === 'Integer Ring') {
    return true;
  }
  // A univariate polynomial ring over a field is a Euclidean domain, hence a
  // PID (SageMath: PolynomialRing_field lies in PrincipalIdealDomains()).
  return isEuclideanRing(ring);
}

/**
 * Infer a ring from a list of elements.
 */
function inferRing(elements: unknown[]): RingLike {
  if (elements.length === 0) {
    // Default to integers
    return {
      zero: () => 0n,
      one: () => 1n,
      is_field: () => false,
    };
  }

  const first = elements[0];

  // Check for bigint
  if (typeof first === 'bigint') {
    return {
      zero: () => 0n,
      one: () => 1n,
      is_field: () => false,
      toString: () => 'Integer Ring',
    };
  }

  // Check for number (assume rationals)
  if (typeof first === 'number') {
    return {
      zero: () => 0,
      one: () => 1,
      is_field: () => true,
      toString: () => 'Rational Field',
    };
  }

  // Check for objects with parent method
  if (typeof first === 'object' && first !== null && 'parent' in first) {
    const parent = (first as { parent: () => RingLike }).parent();
    return parent;
  }

  // Default to integers
  return {
    zero: () => 0n,
    one: () => 1n,
    is_field: () => false,
  };
}
// ============================================================================
// Exact arithmetic helpers
//
// SageMath performs all of the linear algebra below over the base ring (or its
// fraction field) with exact arithmetic: Hermite normal form over ZZ, reduced
// row echelon form over a field.  These helpers provide the same exactness for
// the loosely typed `RingLike` rings used by this port.
// ============================================================================

/**
 * Exact arithmetic in the fraction field of a base ring.
 *
 * `lift` maps an entry of the module into the fraction field, `lower` maps a
 * fraction field element back to the representation used for entries.
 */
interface FractionFieldArithmetic {
  /** Whether the base ring is its own fraction field. */
  readonly isField: boolean;
  /** Whether the base ring is ZZ (entries are bigints, echelon form is HNF). */
  readonly isIntegral: boolean;
  /**
   * Whether the base ring is a Euclidean domain that is neither ZZ nor a field
   * (e.g. `QQ[x]`); the echelon form is then `_echelon_form_PID`.
   */
  readonly isEuclidean: boolean;
  /** Whether exact arithmetic in the fraction field is available at all. */
  readonly exact: boolean;
  zero(): unknown;
  one(): unknown;
  add(a: unknown, b: unknown): unknown;
  sub(a: unknown, b: unknown): unknown;
  mul(a: unknown, b: unknown): unknown;
  div(a: unknown, b: unknown): unknown;
  neg(a: unknown): unknown;
  isZero(a: unknown): boolean;
  eq(a: unknown, b: unknown): boolean;
  lift(x: unknown): unknown;
  lower(x: unknown): unknown;
  /** Whether a fraction field element belongs to the base ring. */
  inBaseRing(x: unknown): boolean;
  /**
   * Denominator of a fraction field element, as an element of the base ring.
   * This is a bigint over ZZ/QQ and a ring element over a Euclidean domain.
   */
  denominator(x: unknown): unknown;
  /** LCM of two denominators, in the base ring. */
  denominatorLcm(a: unknown, b: unknown): unknown;
  /** The denominator `1`, in the base ring. */
  denominatorOne(): unknown;
}

/**
 * Convert an arbitrary entry to an exact rational.
 */
function toRational(x: unknown): Rational {
  if (x instanceof Rational) {
    return x;
  }
  if (typeof x === 'bigint') {
    return new Rational(x);
  }
  if (typeof x === 'number') {
    return Rational.from(x);
  }
  if (typeof x === 'object' && x !== null) {
    const value = (x as { value?: unknown }).value;
    if (typeof value === 'bigint') {
      return new Rational(value);
    }
    const num = (x as { numerator?: unknown }).numerator;
    const den = (x as { denominator?: unknown }).denominator;
    if (typeof num === 'bigint' && typeof den === 'bigint') {
      return new Rational(num, den);
    }
  }
  return Rational.from(String(x));
}

/**
 * Exact arithmetic over QQ, used for base rings whose elements are bigints
 * (ZZ), JavaScript numbers, or {@link Rational}s.
 */
class RationalArithmetic implements FractionFieldArithmetic {
  readonly isField: boolean;
  readonly isIntegral: boolean;
  readonly isEuclidean = false;
  readonly exact = true;
  private readonly mode: 'bigint' | 'number' | 'rational';

  constructor(mode: 'bigint' | 'number' | 'rational', isField: boolean) {
    this.mode = mode;
    this.isField = isField;
    this.isIntegral = mode === 'bigint';
  }

  zero(): unknown {
    return Rational.zero();
  }
  one(): unknown {
    return Rational.one();
  }
  add(a: unknown, b: unknown): unknown {
    return (a as Rational).add(b as Rational);
  }
  sub(a: unknown, b: unknown): unknown {
    return (a as Rational).sub(b as Rational);
  }
  mul(a: unknown, b: unknown): unknown {
    return (a as Rational).mul(b as Rational);
  }
  div(a: unknown, b: unknown): unknown {
    return (a as Rational).div(b as Rational);
  }
  neg(a: unknown): unknown {
    return (a as Rational).neg();
  }
  isZero(a: unknown): boolean {
    return (a as Rational).isZero();
  }
  eq(a: unknown, b: unknown): boolean {
    return (a as Rational).eq(b as Rational);
  }
  lift(x: unknown): unknown {
    return toRational(x);
  }
  lower(x: unknown): unknown {
    const r = x as Rational;
    if (this.mode === 'number') {
      return r.toNumber();
    }
    if (this.mode === 'rational') {
      return r;
    }
    return r.isInteger() ? r.numerator : r;
  }
  inBaseRing(x: unknown): boolean {
    if (this.mode === 'bigint') {
      return (x as Rational).isInteger();
    }
    return true;
  }
  denominator(x: unknown): unknown {
    return (x as Rational).denominator;
  }
  denominatorLcm(a: unknown, b: unknown): unknown {
    return bigintLcm(a as bigint, b as bigint);
  }
  denominatorOne(): unknown {
    return 1n;
  }
}

/**
 * Exact arithmetic using the ring elements themselves; used when the base ring
 * is a field whose elements provide `div` (or `inv`).
 */
class RingElementArithmetic implements FractionFieldArithmetic {
  readonly isField = true;
  readonly isIntegral = false;
  readonly isEuclidean = false;
  readonly exact = true;
  private readonly ring: RingLike;

  constructor(ring: RingLike) {
    this.ring = ring;
  }

  zero(): unknown {
    return this.ring.zero();
  }
  one(): unknown {
    return this.ring.one();
  }
  add(a: unknown, b: unknown): unknown {
    return (a as { add: (x: unknown) => unknown }).add(b);
  }
  sub(a: unknown, b: unknown): unknown {
    return (a as { sub: (x: unknown) => unknown }).sub(b);
  }
  mul(a: unknown, b: unknown): unknown {
    return (a as { mul: (x: unknown) => unknown }).mul(b);
  }
  div(a: unknown, b: unknown): unknown {
    const x = a as { div?: (y: unknown) => unknown; mul: (y: unknown) => unknown };
    if (typeof x.div === 'function') {
      return x.div(b);
    }
    const y = b as { inv?: () => unknown; inverse?: () => unknown };
    if (typeof y.inv === 'function') {
      return x.mul(y.inv());
    }
    if (typeof y.inverse === 'function') {
      return x.mul(y.inverse());
    }
    throw new NotImplementedError('base ring elements do not support division');
  }
  neg(a: unknown): unknown {
    return (a as { neg: () => unknown }).neg();
  }
  isZero(a: unknown): boolean {
    return (a as { isZero: () => boolean }).isZero();
  }
  eq(a: unknown, b: unknown): boolean {
    return (a as { eq: (x: unknown) => boolean }).eq(b);
  }
  lift(x: unknown): unknown {
    if (typeof x === 'object' && x !== null && typeof (x as { add?: unknown }).add === 'function') {
      return x;
    }
    if (this.ring.__call__) {
      return this.ring.__call__(x);
    }
    return x;
  }
  lower(x: unknown): unknown {
    return x;
  }
  inBaseRing(_x: unknown): boolean {
    return true;
  }
  denominator(_x: unknown): unknown {
    return 1n;
  }
  denominatorLcm(_a: unknown, _b: unknown): unknown {
    return 1n;
  }
  denominatorOne(): unknown {
    return 1n;
  }
}

/**
 * Arithmetic for rings for which no fraction field is available.  Echelon
 * forms and linear solving are not implemented over such rings, exactly as in
 * SageMath, where the corresponding module class stores its generators
 * verbatim (`Submodule_free_ambient`).
 */
class InexactArithmetic implements FractionFieldArithmetic {
  readonly isField = false;
  readonly isIntegral = false;
  readonly isEuclidean = false;
  readonly exact = false;
  private readonly ring: RingLike;

  constructor(ring: RingLike) {
    this.ring = ring;
  }

  private fail(): never {
    throw new NotImplementedError('exact linear algebra is not implemented over this base ring');
  }

  zero(): unknown {
    return this.ring.zero();
  }
  one(): unknown {
    return this.ring.one();
  }
  add(_a: unknown, _b: unknown): unknown {
    this.fail();
  }
  sub(_a: unknown, _b: unknown): unknown {
    this.fail();
  }
  mul(_a: unknown, _b: unknown): unknown {
    this.fail();
  }
  div(_a: unknown, _b: unknown): unknown {
    this.fail();
  }
  neg(_a: unknown): unknown {
    this.fail();
  }
  isZero(_a: unknown): boolean {
    this.fail();
  }
  eq(_a: unknown, _b: unknown): boolean {
    this.fail();
  }
  lift(_x: unknown): unknown {
    this.fail();
  }
  lower(x: unknown): unknown {
    return x;
  }
  inBaseRing(_x: unknown): boolean {
    return true;
  }
  denominator(_x: unknown): unknown {
    return 1n;
  }
  denominatorLcm(_a: unknown, _b: unknown): unknown {
    return 1n;
  }
  denominatorOne(): unknown {
    return 1n;
  }
}

// ============================================================================
// Euclidean base rings (e.g. QQ[x]) and their fraction field
// ============================================================================

/**
 * The operations a base ring element must provide for the ring to be treated
 * as a Euclidean domain (SageMath's `EuclideanDomains` category).
 *
 * `sage.rings.polynomial.polynomial_element.Polynomial` provides all of them.
 */
interface EuclideanElement {
  add(other: EuclideanElement): EuclideanElement;
  sub(other: EuclideanElement): EuclideanElement;
  mul(other: EuclideanElement): EuclideanElement;
  neg(): EuclideanElement;
  eq(other: EuclideanElement): boolean;
  isZero(): boolean;
  degree(): number;
  monic(): EuclideanElement;
  quo_rem(other: EuclideanElement): [EuclideanElement, EuclideanElement];
  gcd(other: EuclideanElement): EuclideanElement;
  xgcd(other: EuclideanElement): [EuclideanElement, EuclideanElement, EuclideanElement];
  toString(): string;
}

/**
 * Return whether `ring` is a univariate polynomial ring over a field, i.e. a
 * Euclidean domain (hence a PID) whose elements support `quo_rem`, `gcd` and
 * `xgcd`.
 *
 * This is exactly SageMath's criterion: `PolynomialRing(K, 'x')` lies in
 * `EuclideanDomains()` (and therefore in `PrincipalIdealDomains()`) precisely
 * when `K` is a field.
 *
 * @see Reference: sage/rings/polynomial/polynomial_ring.py:PolynomialRing_field
 */
function isEuclideanRing(ring: RingLike): boolean {
  const r = ring as {
    variable_name?: unknown;
    base_ring?: { is_field?: () => boolean };
    zero?: () => unknown;
  };
  if (typeof r.variable_name !== 'string') {
    return false;
  }
  if (!r.base_ring || typeof r.base_ring.is_field !== 'function' || !r.base_ring.is_field()) {
    return false;
  }
  if (isField(ring)) {
    return false;
  }
  let zero: unknown;
  try {
    zero = ring.zero();
  } catch {
    return false;
  }
  const z = zero as Record<string, unknown>;
  for (const m of [
    'add',
    'sub',
    'mul',
    'neg',
    'eq',
    'isZero',
    'degree',
    'monic',
    'quo_rem',
    'gcd',
    'xgcd',
  ]) {
    if (typeof z[m] !== 'function') {
      return false;
    }
  }
  return true;
}

/**
 * An element `num/den` of the fraction field of a Euclidean domain `R`
 * (for `R = QQ[x]` this is the field of rational functions `QQ(x)`).
 *
 * Instances are always normalized: `gcd(num, den) = 1` and `den` is monic,
 * so `0` is `0/1` and equality is entrywise.
 *
 * @see Reference: sage/rings/fraction_field_element.pyx:FractionFieldElement
 */
export class FractionFieldElement {
  readonly num: EuclideanElement;
  readonly den: EuclideanElement;

  private constructor(num: EuclideanElement, den: EuclideanElement) {
    this.num = num;
    this.den = den;
  }

  /** Create `num/den`, reduced to lowest terms with a monic denominator. */
  static make(num: EuclideanElement, den: EuclideanElement): FractionFieldElement {
    if (den.isZero()) {
      throw new ZeroDivisionError('fraction has denominator zero');
    }
    if (num.isZero()) {
      const [one] = den.quo_rem(den);
      return new FractionFieldElement(num, one);
    }
    let n = num;
    let d = den;
    const g = n.gcd(d);
    if (!(g.degree() === 0 && g.eq(g.monic()))) {
      n = exactQuotient(n, g);
      d = exactQuotient(d, g);
    }
    // Make the denominator monic: d = u * monic(d) with u a unit of R.
    const dm = d.monic();
    if (!d.eq(dm)) {
      const [u] = d.quo_rem(dm);
      n = exactQuotient(n, u);
      d = dm;
    }
    return new FractionFieldElement(n, d);
  }

  /** The element `x/1`, where `one` is the identity of the base ring. */
  static integral(x: EuclideanElement, one: EuclideanElement): FractionFieldElement {
    return new FractionFieldElement(x, one);
  }

  isIntegral(): boolean {
    return this.den.degree() === 0;
  }

  add(other: FractionFieldElement): FractionFieldElement {
    return FractionFieldElement.make(
      this.num.mul(other.den).add(other.num.mul(this.den)),
      this.den.mul(other.den)
    );
  }
  sub(other: FractionFieldElement): FractionFieldElement {
    return FractionFieldElement.make(
      this.num.mul(other.den).sub(other.num.mul(this.den)),
      this.den.mul(other.den)
    );
  }
  mul(other: FractionFieldElement): FractionFieldElement {
    return FractionFieldElement.make(this.num.mul(other.num), this.den.mul(other.den));
  }
  div(other: FractionFieldElement): FractionFieldElement {
    if (other.num.isZero()) {
      throw new ZeroDivisionError('division by zero');
    }
    return FractionFieldElement.make(this.num.mul(other.den), this.den.mul(other.num));
  }
  neg(): FractionFieldElement {
    return new FractionFieldElement(this.num.neg(), this.den);
  }
  inv(): FractionFieldElement {
    if (this.num.isZero()) {
      throw new ZeroDivisionError('division by zero');
    }
    return FractionFieldElement.make(this.den, this.num);
  }
  isZero(): boolean {
    return this.num.isZero();
  }
  eq(other: FractionFieldElement): boolean {
    return this.num.eq(other.num) && this.den.eq(other.den);
  }
  toString(): string {
    if (this.isIntegral() && this.den.eq(this.den.monic())) {
      return this.num.toString();
    }
    return `(${this.num.toString()})/(${this.den.toString()})`;
  }
}

/**
 * Divide exactly in a Euclidean domain, raising if the division is inexact.
 */
function exactQuotient(a: EuclideanElement, b: EuclideanElement): EuclideanElement {
  const [q, r] = a.quo_rem(b);
  if (!r.isZero()) {
    throw new ArithmeticError(`${a.toString()} is not divisible by ${b.toString()}`);
  }
  return q;
}

/**
 * Exact arithmetic in the fraction field of a Euclidean domain such as
 * `QQ[x]`, whose fraction field is `QQ(x)`.
 */
class EuclideanArithmetic implements FractionFieldArithmetic {
  readonly isField = false;
  readonly isIntegral = false;
  readonly isEuclidean = true;
  readonly exact = true;
  readonly ring: RingLike;

  constructor(ring: RingLike) {
    this.ring = ring;
  }

  /** The zero of the base ring, as a Euclidean element. */
  ringZero(): EuclideanElement {
    return this.ring.zero() as EuclideanElement;
  }
  /** The one of the base ring, as a Euclidean element. */
  ringOne(): EuclideanElement {
    return this.ring.one() as EuclideanElement;
  }

  zero(): unknown {
    return FractionFieldElement.integral(this.ringZero(), this.ringOne());
  }
  one(): unknown {
    const one = this.ringOne();
    return FractionFieldElement.integral(one, one);
  }
  add(a: unknown, b: unknown): unknown {
    return (a as FractionFieldElement).add(b as FractionFieldElement);
  }
  sub(a: unknown, b: unknown): unknown {
    return (a as FractionFieldElement).sub(b as FractionFieldElement);
  }
  mul(a: unknown, b: unknown): unknown {
    return (a as FractionFieldElement).mul(b as FractionFieldElement);
  }
  div(a: unknown, b: unknown): unknown {
    return (a as FractionFieldElement).div(b as FractionFieldElement);
  }
  neg(a: unknown): unknown {
    return (a as FractionFieldElement).neg();
  }
  isZero(a: unknown): boolean {
    return (a as FractionFieldElement).isZero();
  }
  eq(a: unknown, b: unknown): boolean {
    return (a as FractionFieldElement).eq(b as FractionFieldElement);
  }
  lift(x: unknown): unknown {
    if (x instanceof FractionFieldElement) {
      return x;
    }
    if (
      typeof x === 'object' &&
      x !== null &&
      typeof (x as { quo_rem?: unknown }).quo_rem === 'function'
    ) {
      return FractionFieldElement.integral(x as EuclideanElement, this.ringOne());
    }
    if (this.ring.__call__) {
      return FractionFieldElement.integral(
        this.ring.__call__(x) as EuclideanElement,
        this.ringOne()
      );
    }
    throw new NotImplementedError(`cannot lift ${String(x)} into the fraction field`);
  }
  lower(x: unknown): unknown {
    const f = x as FractionFieldElement;
    return f.isIntegral() ? this.baseRingElement(f) : f;
  }
  /** `f` with denominator a unit, as an element of the base ring. */
  private baseRingElement(f: FractionFieldElement): EuclideanElement {
    if (f.den.eq(f.den.monic())) {
      return f.num;
    }
    return exactQuotient(f.num, f.den);
  }
  inBaseRing(x: unknown): boolean {
    return (x as FractionFieldElement).isIntegral();
  }
  denominator(x: unknown): unknown {
    return (x as FractionFieldElement).den;
  }
  denominatorLcm(a: unknown, b: unknown): unknown {
    const x = a as EuclideanElement;
    const y = b as EuclideanElement;
    if (x.isZero() || y.isZero()) {
      return this.ringZero();
    }
    return exactQuotient(x.mul(y), x.gcd(y)).monic();
  }
  denominatorOne(): unknown {
    return this.ringOne();
  }

  /** `f * d` as an element of the base ring; `d` must clear the denominator. */
  scaleToRing(f: FractionFieldElement, d: EuclideanElement): EuclideanElement {
    return exactQuotient(f.num.mul(d), f.den);
  }
  /** `x / d` as an element of the fraction field. */
  divideByRing(x: EuclideanElement, d: EuclideanElement): FractionFieldElement {
    return FractionFieldElement.make(x, d);
  }
}

// ----------------------------------------------------------------------------
// Echelon form over a PID that is neither ZZ nor a field.
//
// SageMath computes it with Matrix._echelon_form_PID (matrix2.pyx:17305),
// which recursively applies _generic_clear_column (matrix2.pyx:20613).  The
// reduction of the entries *above* each pivot at the end of
// _echelon_form_PID is guarded by `except AttributeError` on
// `Ideal.small_residue`, which univariate polynomial ring ideals do not have,
// so for K[x] it never runs; we omit it for the same reason.
// ----------------------------------------------------------------------------

/**
 * Invert `a` modulo `m` in a Euclidean domain.
 *
 * @see Reference: sage/rings/polynomial/polynomial_element.pyx:1524 (inverse_mod)
 */
function euclideanInverseMod(
  a: EuclideanElement,
  m: EuclideanElement,
  ar: EuclideanArithmetic
): EuclideanElement {
  const [g, s] = a.xgcd(m);
  const one = ar.ringOne();
  if (g.eq(one)) {
    // s is already reduced modulo m by the Euclidean algorithm.
    return m.isZero() ? s : s.quo_rem(m)[1];
  }
  if (g.degree() === 0 && !g.isZero()) {
    // g is a unit: multiply through by its inverse.
    const inv = exactQuotient(one, g);
    const t = inv.mul(s);
    return m.isZero() ? t : t.quo_rem(m)[1];
  }
  throw new ArithmeticError('Impossible inverse modulo');
}

/** The n x n identity matrix over the base ring. */
function euclideanIdentity(n: number, ar: EuclideanArithmetic): EuclideanElement[][] {
  const zero = ar.ringZero();
  const one = ar.ringOne();
  const I: EuclideanElement[][] = [];
  for (let i = 0; i < n; i++) {
    const row: EuclideanElement[] = [];
    for (let j = 0; j < n; j++) {
      row.push(i === j ? one : zero);
    }
    I.push(row);
  }
  return I;
}

/** Matrix product over the base ring. */
function euclideanMatMul(
  A: EuclideanElement[][],
  B: EuclideanElement[][],
  ar: EuclideanArithmetic
): EuclideanElement[][] {
  const n = A.length;
  const k = B.length;
  const m = k === 0 ? 0 : B[0]!.length;
  const zero = ar.ringZero();
  const C: EuclideanElement[][] = [];
  for (let i = 0; i < n; i++) {
    const row: EuclideanElement[] = [];
    for (let j = 0; j < m; j++) {
      let acc = zero;
      for (let t = 0; t < k; t++) {
        const x = A[i]![t]!;
        if (x.isZero()) continue;
        acc = acc.add(x.mul(B[t]![j]!));
      }
      row.push(acc);
    }
    C.push(row);
  }
  return C;
}

/**
 * Reduce the first column of `m` to canonical form by left multiplication with
 * an invertible matrix over the base ring.
 *
 * @returns `[left, a]` with `left * m == a`
 * @see Reference: sage/matrix/matrix2.pyx:20613 (_generic_clear_column)
 */
function genericClearColumn(
  m: EuclideanElement[][],
  ar: EuclideanArithmetic
): { left: EuclideanElement[][]; a: EuclideanElement[][] } {
  const nrows = m.length;
  const ncols = nrows === 0 ? 0 : m[0]!.length;
  if (nrows <= 1 || ncols <= 0) {
    return { left: euclideanIdentity(nrows, ar), a: m.map((row) => [...row]) };
  }

  let a = m.map((row) => [...row]);
  let left = euclideanIdentity(nrows, ar);
  const zero = ar.ringZero();
  const one = ar.ringOne();

  // case 1: if a[0,0] == 0 and a[k,0] != 0 for some k, swap rows 0 and k.
  if (a[0]![0]!.isZero()) {
    let k = 0;
    while (a[k]![0]!.isZero()) {
      k += 1;
      if (k === nrows) {
        // first column is zero
        return { left, a };
      }
    }
    const swap = euclideanIdentity(nrows, ar);
    swap[0]![0] = zero;
    swap[k]![k] = zero;
    swap[0]![k] = one;
    swap[k]![0] = one.neg();
    left = swap;
    a = euclideanMatMul(left, m, ar);
  }

  // case 2: clear the column with the 2x2 unimodular matrix built from the
  // gcd of a[0,0] and a[k,0].
  for (let k = 1; k < nrows; k++) {
    const a00 = a[0]![0]!;
    const ak0 = a[k]![0]!;
    if (ak0.quo_rem(a00)[1]!.isZero()) continue; // a[k,0] in ideal(a[0,0])

    const B = a00.gcd(ak0);
    const c = euclideanInverseMod(exactQuotient(a00, B), exactQuotient(ak0, B), ar);
    const d = exactQuotient(c.mul(a00).sub(B), ak0);
    if (!c.mul(a00).sub(d.mul(ak0)).eq(B)) {
      throw new ArithmeticError('failed to clear column: c*a00 - d*ak0 != gcd');
    }
    let e: EuclideanElement;
    let f: EuclideanElement;
    if (!c.isZero()) {
      e = euclideanInverseMod(d, c, ar);
      f = exactQuotient(one.sub(d.mul(e)), c);
    } else {
      e = exactQuotient(ak0, B).neg();
      f = one;
    }
    if (!e.mul(d).add(c.mul(f)).eq(one)) {
      throw new ArithmeticError('failed to clear column: e*d + c*f != 1');
    }
    const newlmat = euclideanIdentity(nrows, ar);
    newlmat[0]![0] = c;
    newlmat[0]![k] = d.neg();
    newlmat[k]![0] = e;
    newlmat[k]![k] = f;
    a = euclideanMatMul(newlmat, a, ar);
    left = euclideanMatMul(newlmat, left, ar);
  }

  // now everything in column 0 is divisible by the pivot
  const pivot = a[0]![0]!;
  if (!pivot.isZero()) {
    for (let i = 1; i < nrows; i++) {
      const s = exactQuotient(a[i]![0]!, pivot);
      if (s.isZero()) continue;
      for (let j = 0; j < ncols; j++) {
        a[i]![j] = a[i]![j]!.sub(s.mul(a[0]![j]!));
      }
      for (let j = 0; j < nrows; j++) {
        left[i]![j] = left[i]![j]!.sub(s.mul(left[0]![j]!));
      }
    }
  }

  return { left, a };
}

/**
 * A matrix over a Euclidean base ring, carrying its shape explicitly so that
 * empty submatrices keep their column (or row) count.
 */
interface EMat {
  readonly r: number;
  readonly c: number;
  readonly a: EuclideanElement[][];
}

function emOf(a: EuclideanElement[][], r: number, c: number): EMat {
  return { r, c, a };
}

function emIdentity(n: number, ar: EuclideanArithmetic): EMat {
  return emOf(euclideanIdentity(n, ar), n, n);
}

function emMul(A: EMat, B: EMat, ar: EuclideanArithmetic): EMat {
  const zero = ar.ringZero();
  const out: EuclideanElement[][] = [];
  for (let i = 0; i < A.r; i++) {
    const row: EuclideanElement[] = [];
    for (let j = 0; j < B.c; j++) {
      let acc = zero;
      for (let t = 0; t < A.c; t++) {
        const x = A.a[i]![t]!;
        if (x.isZero()) continue;
        acc = acc.add(x.mul(B.a[t]![j]!));
      }
      row.push(acc);
    }
    out.push(row);
  }
  return emOf(out, A.r, B.c);
}

function emTranspose(M: EMat): EMat {
  const out: EuclideanElement[][] = [];
  for (let j = 0; j < M.c; j++) {
    const row: EuclideanElement[] = [];
    for (let i = 0; i < M.r; i++) {
      row.push(M.a[i]![j]!);
    }
    out.push(row);
  }
  return emOf(out, M.c, M.r);
}

/** `[[x]] (+) M`, the block diagonal sum with a 1x1 block. */
function emBlockSum1(x: EuclideanElement, M: EMat, ar: EuclideanArithmetic): EMat {
  const zero = ar.ringZero();
  const out: EuclideanElement[][] = [[x, ...new Array<EuclideanElement>(M.c).fill(zero)]];
  for (let i = 0; i < M.r; i++) {
    out.push([zero, ...M.a[i]!]);
  }
  return emOf(out, M.r + 1, M.c + 1);
}

/** The submatrix obtained by dropping the first row and the first column. */
function emDropFirst(M: EMat): EMat {
  return emOf(
    M.a.slice(1).map((row) => row.slice(1)),
    Math.max(M.r - 1, 0),
    Math.max(M.c - 1, 0)
  );
}

function emIsZero(M: EMat): boolean {
  return M.a.every((row) => row.every((e) => e.isZero()));
}

/**
 * Given a diagonal matrix `d`, return `dp, a, b` with `a*d*b = dp` and `dp`
 * diagonal with each entry dividing the next.
 *
 * @see Reference: sage/matrix/matrix2.pyx:20537 (_smith_diag)
 */
function euclideanSmithDiag(
  d: EMat,
  ar: EuclideanArithmetic
): { dp: EMat; left: EMat; right: EMat } {
  let dp = emOf(
    d.a.map((row) => [...row]),
    d.r,
    d.c
  );
  const n = Math.min(d.r, d.c);
  let left = emIdentity(d.r, ar);
  let right = emIdentity(d.c, ar);
  const one = ar.ringOne();

  for (let i = 0; i < n; i++) {
    const dii = dp.a[i]![i]!;
    if (!dii.isZero() && dii.degree() === 0) {
      // ideal(dp[i,i]) is the unit ideal: normalize the entry to 1.
      if (!dii.eq(one)) {
        const scale = exactQuotient(one, dii);
        for (let j = 0; j < d.r; j++) {
          left.a[i]![j] = left.a[i]![j]!.mul(scale);
        }
        dp.a[i]![i] = one;
      }
      continue;
    }
    for (let j = i + 1; j < n; j++) {
      const a = dp.a[i]![i]!;
      const b = dp.a[j]![j]!;
      if (b.isZero()) continue;
      if (!a.isZero() && b.quo_rem(a)[1]!.isZero()) continue; // dp[j,j] in (dp[i,i])
      const t = a.gcd(b);
      const lamb = euclideanInverseMod(exactQuotient(a, t), exactQuotient(b, t), ar);
      const mu = exactQuotient(t.sub(lamb.mul(a)), b);

      const newl = emIdentity(d.r, ar);
      newl.a[i]![i] = lamb;
      newl.a[i]![j] = one;
      newl.a[j]![i] = exactQuotient(b.neg().mul(mu), t);
      newl.a[j]![j] = exactQuotient(a, t);

      const newr = emIdentity(d.c, ar);
      newr.a[i]![i] = one;
      newr.a[i]![j] = exactQuotient(b, t).neg();
      newr.a[j]![i] = mu;
      newr.a[j]![j] = exactQuotient(lamb.mul(a), t);

      left = emMul(newl, left, ar);
      right = emMul(right, newr, ar);
      dp = emMul(emMul(newl, dp, ar), newr, ar);
    }
  }
  return { dp, left, right };
}

/**
 * One step of the Smith normal form: returns `a, b, c` with `a*m*c = b` and
 * row 0 and column 0 of `b` zero apart from `b[0,0]`.
 *
 * @see Reference: sage/matrix/matrix2.pyx:20730 (_smith_onestep)
 */
function euclideanSmithOnestep(
  m: EMat,
  ar: EuclideanArithmetic
): { left: EMat; a: EMat; right: EMat } {
  let left = emIdentity(m.r, ar);
  let right = emIdentity(m.c, ar);
  let a = emOf(
    m.a.map((row) => [...row]),
    m.r,
    m.c
  );

  if (emIsZero(m) || (m.r <= 1 && m.c <= 1)) {
    return { left, a, right };
  }

  const zero = ar.ringZero();
  const one = ar.ringOne();

  // preparation: if column 0 is zero, swap it with the first nonzero column
  let j = 0;
  while (a.a.every((row) => row[j]!.isZero())) {
    j += 1;
  }
  if (j > 0) {
    right.a[0]![0] = zero;
    right.a[j]![j] = zero;
    right.a[0]![j] = one;
    right.a[j]![0] = one.neg();
    a = emMul(a, right, ar);
  }

  const cleared = genericClearColumn(a.a, ar);
  left = emOf(cleared.left, m.r, m.r);
  a = emOf(cleared.a, m.r, m.c);

  // test whether everything to the right of the pivot in row 0 is zero
  let isdone = true;
  for (let jj = j + 1; jj < m.c; jj++) {
    if (!a.a[0]![jj]!.isZero()) {
      isdone = false;
    }
  }

  if (!isdone) {
    const sub = euclideanSmithOnestep(emTranspose(a), ar);
    left = emMul(emTranspose(sub.right), left, ar);
    a = emTranspose(sub.a);
    right = emMul(right, emTranspose(sub.left), ar);
  }

  return { left, a, right };
}

/**
 * Smith normal form over a Euclidean base ring: `u * m * v == d`.
 *
 * @see Reference: sage/matrix/matrix2.pyx:16732 (smith_form, generic branch)
 */
function euclideanSmithForm(m: EMat, ar: EuclideanArithmetic): { d: EMat; u: EMat; v: EMat } {
  if (emIsZero(m) || (m.r <= 1 && m.c <= 1)) {
    return {
      d: emOf(
        m.a.map((row) => [...row]),
        m.r,
        m.c
      ),
      u: emIdentity(m.r, ar),
      v: emIdentity(m.c, ar),
    };
  }

  const step = euclideanSmithOnestep(m, ar);
  let u = step.left;
  let v = step.right;
  const t = step.a;

  const rec = euclideanSmithForm(emDropFirst(t), ar);
  const d = emBlockSum1(t.a[0]![0]!, rec.d, ar);
  u = emMul(emBlockSum1(ar.ringOne(), rec.u, ar), u, ar);
  v = emMul(v, emBlockSum1(ar.ringOne(), rec.v, ar), ar);

  const diag = euclideanSmithDiag(d, ar);
  return {
    d: diag.dp,
    u: emMul(diag.left, u, ar),
    v: emMul(v, diag.right, ar),
  };
}

/**
 * Return `p.denominator()`: the lcm of the denominators of the *coefficients*
 * of the polynomial `p`, as a constant of the base ring `R`.
 *
 * This is **not** the denominator of `p` as a rational function (which is
 * always 1); it is the scalar that clears the coefficient denominators, e.g.
 * `(1/17*x^19 - 2/3*x + 1/3).denominator() == 51` over `QQ[x]`.  Upstream
 * falls back to `self.base_ring().one()` whenever the coefficients have no
 * `denominator` method -- which is exactly what happens over `GF(p)[x]`, so
 * this returns 1 there.
 *
 * @see Reference: sage/rings/polynomial/polynomial_element.pyx:4026 (Polynomial.denominator)
 */
function polynomialCoefficientDenominator(p: EuclideanElement): bigint {
  // polynomial_element.pyx:4104 -- `if self.degree() == -1: return one`
  if (p.isZero()) {
    return 1n;
  }
  const coeffs = (p as unknown as { coeffs?: readonly unknown[] }).coeffs;
  if (!Array.isArray(coeffs)) {
    return 1n;
  }
  let d = 1n;
  for (const c of coeffs) {
    // polynomial_element.pyx:4106 -- `x = self.coefficients()` (nonzero only)
    if (c === null || c === undefined) {
      return 1n;
    }
    if ((c as { isZero?: () => boolean }).isZero?.()) {
      continue;
    }
    const den = (c as { denominator?: unknown }).denominator;
    if (typeof den !== 'bigint') {
      // polynomial_element.pyx:4111 -- `except(AttributeError): return one`
      return 1n;
    }
    d = bigintLcm(d, den);
  }
  return d;
}

/**
 * Return `M.denominator()`: the lcm of the denominators of the entries of the
 * matrix `M`, taken over its *coordinate ring*.
 *
 * SageMath builds the basis matrix over `self.coordinate_ring()`, which is the
 * base ring `R` when every generator lies in `R^n` and the fraction field
 * `Frac(R)` otherwise (free_module.py:6690).  The two cases give genuinely
 * different denominators, and `integer_kernel` scales by whichever one applies:
 *
 * - over `R = K[x]`, `entry.denominator()` is the lcm of the denominators of
 *   the *coefficients* (an element of `K`, e.g. `51` over `QQ[x]`);
 * - over `Frac(R)`, it is the polynomial denominator of the rational function.
 *
 * @see Reference: sage/matrix/matrix2.pyx:3521 (Matrix.denominator)
 * @see Reference: sage/modules/free_module.py:2574 (basis_matrix, over coordinate_ring)
 */
function euclideanMatrixDenominator(
  lifted: FractionFieldElement[][],
  ear: EuclideanArithmetic
): EuclideanElement {
  // matrix2.pyx:3574 -- `if self.nrows() == 0 or self.ncols() == 0: return ZZ.one()`
  if (lifted.length === 0 || lifted[0]!.length === 0) {
    return ear.ringOne();
  }
  const overFractionField = lifted.some((row) => row.some((e) => !e.isIntegral()));
  if (overFractionField) {
    // Coordinate ring Frac(R): `entry.denominator()` is the polynomial
    // denominator, and the lcm is taken in R.
    let den = ear.ringOne();
    for (const row of lifted) {
      for (const e of row) {
        den = ear.denominatorLcm(den, e.den) as EuclideanElement;
      }
    }
    return den;
  }
  // Coordinate ring R = K[x]: `entry.denominator()` lives in K, so the lcm is
  // taken in K (over QQ, in ZZ) and not in R -- in R every nonzero constant is
  // a unit and the lcm would collapse to 1.
  let d = 1n;
  for (const row of lifted) {
    for (const e of row) {
      d = bigintLcm(d, polynomialCoefficientDenominator(ear.scaleToRing(e, ear.ringOne())));
    }
  }
  if (d === 1n || !ear.ring.__call__) {
    return ear.ringOne();
  }
  return ear.ring.__call__(d) as EuclideanElement;
}

/**
 * Return a basis of the left kernel `{c in R^k : c * S = 0}` over a Euclidean
 * base ring `R`, where the entries of `S` may lie in the fraction field.
 *
 * This is SageMath's `integer_kernel(R)`: multiply by the single denominator
 * `S.denominator()` (which does not change the kernel as a *set*, but does
 * change the basis the Smith form produces, because nothing downstream
 * normalizes the pivots) and take the kernel over `R`.  The kernel itself is
 * `left_kernel = transpose().right_kernel()`, computed from the Smith normal
 * form, and `right_kernel(basis='echelon')` echelonizes the Smith-form basis
 * with `_echelon_form_PID` before wrapping it in a module with
 * `already_echelonized=True`.
 *
 * @see Reference: sage/matrix/matrix2.pyx:5641 (integer_kernel)
 * @see Reference: sage/matrix/matrix2.pyx:4166 (_right_kernel_matrix_over_domain)
 * @see Reference: sage/matrix/matrix2.pyx:4975 (right_kernel), :5345 (left_kernel)
 */
function euclideanLeftKernelRows(S: unknown[][], ear: EuclideanArithmetic): EuclideanElement[][] {
  if (S.length === 0) {
    return [];
  }
  const lifted = liftRows(S, ear) as FractionFieldElement[][];
  const den = euclideanMatrixDenominator(lifted, ear);
  const cleared = lifted.map((row) =>
    row.map((e) => ear.scaleToRing(e as FractionFieldElement, den))
  );

  // left kernel of `cleared` = right kernel of its transpose
  const T = emTranspose(emOf(cleared, cleared.length, cleared[0]!.length));
  if (T.c === 0) {
    return [];
  }
  if (T.r === 0) {
    return emIdentity(T.c, ear).a;
  }
  const { d, v } = euclideanSmithForm(T, ear);
  const basis: EuclideanElement[][] = [];
  for (let i = 0; i < T.c; i++) {
    if (i >= T.r || d.a[i]![i]!.isZero()) {
      basis.push(v.a.map((row) => row[i]!));
    }
  }
  if (basis.length === 0) {
    return basis;
  }
  // `right_kernel(basis='echelon')` (matrix2.pyx:4941) sees the format string
  // 'computed-smith-form', so it runs `M.echelonize()` -- i.e.
  // `_echelon_form_PID` -- on the Smith-form basis and then wraps the rows in
  // a module with `already_echelonized=True`, which keeps them verbatim.  The
  // Smith-form rows are independent, so `_echelon_form_PID` produces no zero
  // rows; dropping them here only guards against a degenerate input.
  const { a } = echelonFormPID(basis, ear);
  return a.filter((row) => row.some((e) => !e.isZero()));
}

/**
 * Return `[left, a, pivots]` with `left * self == a`, `a` in Hermite normal
 * form and `left` invertible over the base ring.
 *
 * @see Reference: sage/matrix/matrix2.pyx:17305 (_echelon_form_PID)
 */
function echelonFormPID(
  M: EuclideanElement[][],
  ar: EuclideanArithmetic
): { left: EuclideanElement[][]; a: EuclideanElement[][]; pivots: number[] } {
  const nrows = M.length;
  const ncols = nrows === 0 ? 0 : M[0]!.length;

  if (ncols === 0 || nrows === 0) {
    return { left: euclideanIdentity(nrows, ar), a: M.map((row) => [...row]), pivots: [] };
  }
  if (nrows === 1) {
    const row = M[0]!;
    const j = row.findIndex((e) => !e.isZero());
    return {
      left: euclideanIdentity(1, ar),
      a: [[...row]],
      pivots: j === -1 ? [] : [j],
    };
  }

  const cleared = genericClearColumn(M, ar);
  let left = cleared.left;
  let a = cleared.a;
  let pivots: number[];

  if (!a[0]![0]!.isZero()) {
    const aa = a.slice(1).map((row) => row.slice(1));
    const sub = echelonFormPID(aa, ar);
    // left = block_diag(1, s) * left
    const blocked = euclideanIdentity(nrows, ar);
    for (let i = 0; i < nrows - 1; i++) {
      for (let j = 0; j < nrows - 1; j++) {
        blocked[i + 1]![j + 1] = sub.left[i]![j]!;
      }
    }
    left = euclideanMatMul(blocked, left, ar);
    a = euclideanMatMul(left, M, ar);
    pivots = [0, ...sub.pivots.map((x) => x + 1)];
  } else {
    const aa = a.map((row) => row.slice(1));
    const sub = echelonFormPID(aa, ar);
    left = euclideanMatMul(sub.left, left, ar);
    a = euclideanMatMul(left, M, ar);
    pivots = sub.pivots.map((x) => x + 1);
  }

  return { left, a, pivots };
}

/**
 * Return the fraction field of a Euclidean base ring as a {@link RingLike},
 * e.g. `QQ(x)` for `QQ[x]`.
 *
 * @see Reference: sage/rings/fraction_field.py:FractionField_1poly_field
 */
function fractionFieldOf(ring: RingLike): RingLike {
  const ar = new EuclideanArithmetic(ring);
  return {
    zero: () => ar.zero(),
    one: () => ar.one(),
    is_field: () => true,
    is_exact: () => true,
    __call__: (x: unknown) => ar.lift(x),
    toString: () => `Fraction Field of ${ring.toString?.() ?? 'ring'}`,
  };
}

/**
 * Return exact fraction field arithmetic for the given base ring.
 */
function arithmeticFor(ring: RingLike): FractionFieldArithmetic {
  let zero: unknown;
  try {
    zero = ring.zero();
  } catch {
    return new InexactArithmetic(ring);
  }

  const isFieldRing = isField(ring);

  if (typeof zero === 'bigint') {
    return new RationalArithmetic('bigint', isFieldRing);
  }
  if (typeof zero === 'number') {
    return new RationalArithmetic('number', isFieldRing);
  }
  if (zero instanceof Rational) {
    return new RationalArithmetic('rational', true);
  }
  if (typeof zero === 'object' && zero !== null) {
    const z = zero as {
      add?: unknown;
      mul?: unknown;
      div?: unknown;
      inv?: unknown;
      value?: unknown;
    };
    if (typeof z.value === 'bigint' && typeof z.add !== 'function') {
      return new RationalArithmetic('bigint', isFieldRing);
    }
    // Only fields get element arithmetic: over a ring that is neither ZZ-like
    // nor a field there is no echelon form, exactly as in SageMath, where the
    // generators of a submodule are then stored verbatim.
    if (isFieldRing && typeof z.add === 'function' && typeof z.mul === 'function') {
      return new RingElementArithmetic(ring);
    }
    // A Euclidean domain that is not a field (e.g. QQ[x]) has an exact
    // fraction field; SageMath computes its echelon forms with
    // `Matrix._echelon_form_PID` (matrix2.pyx:17305).
    if (isEuclideanRing(ring)) {
      return new EuclideanArithmetic(ring);
    }
  }

  return new InexactArithmetic(ring);
}

/**
 * Coerce the entries of a row into the given ring.
 *
 * Fractions are mapped to `num/den` in the target ring when that ring is a
 * field, which is how SageMath coerces a QQ-vector into GF(p) in
 * `change_ring`.  Over a non-field they are left alone: the coordinate ring of
 * the resulting module is then the fraction field of the base ring.
 */
function coerceRow(ring: RingLike, row: unknown[]): unknown[] {
  const ar = arithmeticFor(ring);
  return row.map((e) => {
    if (e instanceof Rational && ring.__call__) {
      if (e.isInteger()) {
        return ring.__call__(e.numerator);
      }
      if (isField(ring)) {
        return ar.div(ring.__call__(e.numerator), ring.__call__(e.denominator));
      }
    }
    return e;
  });
}

/**
 * Compute the GCD of two bigints.
 */
function bigintGcd(a: bigint, b: bigint): bigint {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    const t = b;
    b = a % b;
    a = t;
  }
  return a;
}

/**
 * Compute the LCM of two bigints.
 */
function bigintLcm(a: bigint, b: bigint): bigint {
  if (a === 0n || b === 0n) return 0n;
  const g = bigintGcd(a, b);
  const l = (a / g) * b;
  return l < 0n ? -l : l;
}

/**
 * Lift a matrix of entries into the fraction field.
 */
function liftRows(rows: unknown[][], ar: FractionFieldArithmetic): unknown[][] {
  return rows.map((row) => row.map((e) => ar.lift(e)));
}

/**
 * Lower a matrix of fraction field elements back to entries.
 */
function lowerRows(rows: unknown[][], ar: FractionFieldArithmetic): unknown[][] {
  return rows.map((row) => row.map((e) => ar.lower(e)));
}

/**
 * Reduced row echelon form of a lifted matrix over a field.
 *
 * @returns The RREF (without zero rows) and the list of pivot columns
 */
function rrefLifted(
  rows: unknown[][],
  ar: FractionFieldArithmetic
): { rows: unknown[][]; pivots: number[] } {
  const M = rows.map((row) => [...row]);
  const m = M.length;
  const n = m === 0 ? 0 : M[0]!.length;
  const pivots: number[] = [];
  let r = 0;

  for (let col = 0; col < n && r < m; col++) {
    // Find a pivot in this column
    let pivotRow = -1;
    for (let i = r; i < m; i++) {
      if (!ar.isZero(M[i]![col])) {
        pivotRow = i;
        break;
      }
    }
    if (pivotRow === -1) {
      continue;
    }
    if (pivotRow !== r) {
      [M[r], M[pivotRow]] = [M[pivotRow]!, M[r]!];
    }

    // Scale the pivot row so that the pivot is 1
    const pivot = M[r]![col];
    for (let j = col; j < n; j++) {
      M[r]![j] = ar.div(M[r]![j], pivot);
    }

    // Eliminate the column from every other row
    for (let i = 0; i < m; i++) {
      if (i === r) continue;
      const factor = M[i]![col];
      if (ar.isZero(factor)) continue;
      for (let j = col; j < n; j++) {
        M[i]![j] = ar.sub(M[i]![j], ar.mul(factor, M[r]![j]));
      }
    }

    pivots.push(col);
    r++;
  }

  return { rows: M.slice(0, r), pivots };
}

/**
 * Return the echelon form of the given rows over the base ring.
 *
 * Over ZZ this is the Hermite normal form (delegated to
 * `matrix_integer.hermite_normal_form`); over a field it is the reduced row
 * echelon form.  Zero rows are dropped, so the number of rows returned is the
 * rank of the input.
 *
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid._echelonized_basis
 */
function echelonRows(rows: unknown[][], ar: FractionFieldArithmetic): unknown[][] {
  if (rows.length === 0 || rows[0]!.length === 0) {
    return rows.map((row) => [...row]);
  }
  if (!ar.exact) {
    return rows.map((row) => [...row]);
  }

  const lifted = liftRows(rows, ar);

  if (ar.isIntegral) {
    // Clear denominators, take the Hermite normal form, restore denominators.
    let d = 1n;
    for (const row of lifted) {
      for (const e of row) {
        d = bigintLcm(d, ar.denominator(e) as bigint);
      }
    }
    const scaled: bigint[][] = lifted.map((row) =>
      row.map((e) => {
        const r = (e as Rational).mul(new Rational(d));
        if (!r.isInteger()) {
          throw new ArithmeticError('failed to clear denominators of the basis matrix');
        }
        return r.numerator;
      })
    );

    const H = hermite_normal_form(IntegerMatrixFromEntries(scaled), 'default', false, false) as {
      nrows: number;
      ncols: number;
      get: (i: number, j: number) => { value: bigint };
    };

    const out: unknown[][] = [];
    const dr = new Rational(d);
    for (let i = 0; i < H.nrows; i++) {
      const row: unknown[] = [];
      let nonzero = false;
      for (let j = 0; j < H.ncols; j++) {
        const v = H.get(i, j).value;
        if (v !== 0n) nonzero = true;
        row.push(ar.lower(new Rational(v).div(dr)));
      }
      if (nonzero) {
        out.push(row);
      }
    }
    return out;
  }

  if (ar.isEuclidean) {
    // _echelonized_basis (free_module.py:6900):
    //
    //     if basis.universe().coordinate_ring() == ambient.base_ring(): d = 1
    //     else: d = self._denominator(basis)
    //     if d != 1: basis = [x * d for x in basis]
    //     E = matrix(basis).echelon_form()
    //     if d != 1: E = E.matrix_over_field() * (~d)
    //     r = E.rank()
    //     if r < E.nrows(): E = E.matrix_from_rows(range(r))
    //
    // `_denominator` is the lcm of `vector.denominator()`, i.e. of the
    // *fraction field* denominators of the entries; when every generator is
    // already in `R^n` the coordinate ring is `R` and `d` is 1, which is what
    // the loop below produces because every such denominator is 1.  (Note
    // that the coefficient-clearing denominator `Polynomial.denominator()` is
    // deliberately NOT used here -- it is only used by `integer_kernel`, see
    // `euclideanMatrixDenominator`.)
    const ear = ar as EuclideanArithmetic;
    let d = ear.ringOne();
    for (const row of lifted) {
      for (const e of row) {
        d = ear.denominatorLcm(d, ear.denominator(e)) as EuclideanElement;
      }
    }
    const scaled: EuclideanElement[][] = lifted.map((row) =>
      row.map((e) => ear.scaleToRing(e as FractionFieldElement, d))
    );

    const { a } = echelonFormPID(scaled, ear);

    const E: unknown[][] = a.map((row) => row.map((e) => ear.lower(ear.divideByRing(e, d))));
    // `r = E.rank(); if r < E.nrows(): E = E.matrix_from_rows(range(r))` --
    // upstream keeps the *first* r rows, it does not filter out zero rows.
    const r = rankOfRows(E, ar);
    return r < E.length ? E.slice(0, r) : E;
  }

  const { rows: E } = rrefLifted(lifted, ar);
  return lowerRows(E, ar);
}

/**
 * Return the rank of the given rows over the base ring.
 */
function rankOfRows(rows: unknown[][], ar: FractionFieldArithmetic): number {
  if (rows.length === 0 || rows[0]!.length === 0) {
    return 0;
  }
  if (!ar.exact) {
    return rows.length;
  }
  return rrefLifted(liftRows(rows, ar), ar).rows.length;
}

/**
 * Return a basis of the right kernel `{x : A x = 0}` of the given rows,
 * in echelon form (SageMath's default `basis='echelon'`).
 */
function rightKernelRows(
  rows: unknown[][],
  ncols: number,
  ar: FractionFieldArithmetic
): unknown[][] {
  if (ncols === 0) {
    return [];
  }
  const lifted = rows.length === 0 ? [] : liftRows(rows, ar);
  const { rows: E, pivots } = rrefLifted(lifted, ar);

  const isPivot = new Array<boolean>(ncols).fill(false);
  for (const p of pivots) {
    isPivot[p] = true;
  }

  const kernel: unknown[][] = [];
  for (let free = 0; free < ncols; free++) {
    if (isPivot[free]) continue;
    const v: unknown[] = new Array(ncols).fill(ar.zero());
    v[free] = ar.one();
    for (let r = 0; r < pivots.length; r++) {
      v[pivots[r]!] = ar.neg(E[r]![free]);
    }
    kernel.push(v);
  }

  if (kernel.length === 0) {
    return [];
  }
  return lowerRows(rrefLifted(kernel, ar).rows, ar);
}

/**
 * Solve `x * B = v` exactly, where the rows of `B` are the basis vectors.
 *
 * @returns The (lifted) coefficient vector, or `null` if there is no solution
 */
function solveLeftLifted(
  B: unknown[][],
  v: unknown[],
  ar: FractionFieldArithmetic
): unknown[] | null {
  const n = B.length; // number of unknowns
  const m = v.length; // number of equations
  if (n === 0) {
    return v.every((e) => ar.isZero(e)) ? [] : null;
  }

  // Augmented matrix of the transposed system: row j is
  //   B[0][j] x_0 + ... + B[n-1][j] x_{n-1} = v[j]
  const A: unknown[][] = [];
  for (let j = 0; j < m; j++) {
    const row: unknown[] = [];
    for (let i = 0; i < n; i++) {
      row.push(B[i]![j]);
    }
    row.push(v[j]);
    A.push(row);
  }

  const { rows: E, pivots } = rrefLifted(A, ar);

  // Inconsistent if the augmented column is a pivot
  if (pivots.includes(n)) {
    return null;
  }

  const x: unknown[] = new Array(n).fill(ar.zero());
  for (let r = 0; r < pivots.length; r++) {
    x[pivots[r]!] = E[r]![n];
  }

  // Verify (free variables were set to zero, which is only valid if the
  // resulting vector really is a solution)
  for (let j = 0; j < m; j++) {
    let acc = ar.zero();
    for (let i = 0; i < n; i++) {
      acc = ar.add(acc, ar.mul(x[i], B[i]![j]));
    }
    if (!ar.eq(acc, v[j])) {
      return null;
    }
  }

  return x;
}

/**
 * Turn a basis of a QQ-kernel into a basis of the corresponding ZZ-kernel.
 *
 * The rows are scaled to be integral and then saturated, which is what
 * SageMath's `integer_kernel` computes.
 *
 * @see Reference: sage/matrix/matrix2.pyx:integer_kernel
 */
function integralKernelRows(rows: unknown[][], ar: FractionFieldArithmetic): unknown[][] {
  if (rows.length === 0) {
    return rows;
  }
  const cleared: bigint[][] = [];
  for (const row of rows) {
    const lifted = row.map((e) => ar.lift(e) as Rational);
    let d = 1n;
    for (const e of lifted) {
      d = bigintLcm(d, e.denominator);
    }
    cleared.push(lifted.map((e) => e.mul(new Rational(d)).numerator));
  }

  const S = matrixSaturation(IntegerMatrixFromEntries(cleared));
  const out: unknown[][] = [];
  for (let i = 0; i < S.nrows; i++) {
    const row: unknown[] = [];
    for (let j = 0; j < S.ncols; j++) {
      row.push(ar.lower(new Rational(S.get(i, j).value)));
    }
    out.push(row);
  }
  return out;
}

/**
 * Product of two matrices of lifted (fraction field) entries.
 */
function matMulLifted(A: unknown[][], B: unknown[][], ar: FractionFieldArithmetic): unknown[][] {
  const n = A.length;
  const k = B.length;
  const m = k === 0 ? 0 : B[0]!.length;
  const C: unknown[][] = [];
  for (let i = 0; i < n; i++) {
    const row: unknown[] = [];
    for (let j = 0; j < m; j++) {
      let acc = ar.zero();
      for (let t = 0; t < k; t++) {
        acc = ar.add(acc, ar.mul(A[i]![t], B[t]![j]));
      }
      row.push(acc);
    }
    C.push(row);
  }
  return C;
}

/**
 * Inverse of a square matrix of lifted entries, by Gauss-Jordan elimination
 * over the fraction field.
 *
 * @throws {ZeroDivisionError} If the matrix is singular
 */
function inverseLifted(A: unknown[][], ar: FractionFieldArithmetic): unknown[][] {
  const n = A.length;
  if (n === 0) {
    return [];
  }
  const aug: unknown[][] = A.map((row, i) => {
    const r = [...row];
    for (let j = 0; j < n; j++) {
      r.push(i === j ? ar.one() : ar.zero());
    }
    return r;
  });
  const { rows: E, pivots } = rrefLifted(aug, ar);
  if (pivots.length !== n || pivots.some((p, i) => p !== i)) {
    throw new ZeroDivisionError('matrix must be nonsingular');
  }
  return E.map((row) => row.slice(n));
}

/**
 * The Kronecker product of two matrices: each entry `x` of `A` is replaced by
 * `x * Bm`, giving an `(m*p) x (n*q)` matrix.
 *
 * @see Reference: sage/matrix/matrix2.pyx:9983 (Matrix.tensor_product)
 */
function kroneckerProduct(
  A: unknown[][],
  Bm: unknown[][],
  ar: FractionFieldArithmetic
): unknown[][] {
  const m = A.length;
  const n = m === 0 ? 0 : A[0]!.length;
  const p = Bm.length;
  const q = p === 0 ? 0 : Bm[0]!.length;
  const out: unknown[][] = [];
  for (let i = 0; i < m; i++) {
    for (let k = 0; k < p; k++) {
      const row: unknown[] = [];
      for (let j = 0; j < n; j++) {
        for (let l = 0; l < q; l++) {
          row.push(ar.mul(A[i]![j], Bm[k]![l]));
        }
      }
      out.push(row);
    }
  }
  return out;
}

/**
 * Determinant of a square matrix of lifted entries, computed exactly by
 * Gaussian elimination over the fraction field.
 */
function determinantLifted(M: unknown[][], ar: FractionFieldArithmetic): unknown {
  const n = M.length;
  if (n === 0) {
    return ar.one();
  }
  const A = M.map((row) => [...row]);
  let det = ar.one();

  for (let col = 0; col < n; col++) {
    let pivotRow = -1;
    for (let i = col; i < n; i++) {
      if (!ar.isZero(A[i]![col])) {
        pivotRow = i;
        break;
      }
    }
    if (pivotRow === -1) {
      return ar.zero();
    }
    if (pivotRow !== col) {
      [A[col], A[pivotRow]] = [A[pivotRow]!, A[col]!];
      det = ar.neg(det);
    }
    const pivot = A[col]![col];
    det = ar.mul(det, pivot);
    for (let i = col + 1; i < n; i++) {
      const factor = ar.div(A[i]![col], pivot);
      if (ar.isZero(factor)) continue;
      for (let j = col; j < n; j++) {
        A[i]![j] = ar.sub(A[i]![j], ar.mul(factor, A[col]![j]));
      }
    }
  }

  return det;
}

// ============================================================================
// Base classes
// ============================================================================

/**
 * Base class for modules with elements represented by elements of a free module.
 * @see Reference: sage/modules/free_module.py:Module_free_ambient
 */
export abstract class ModuleFreeAmbient implements FreeModuleParent {
  protected _baseRing: RingLike;
  protected _degree: number;
  protected _sparse: boolean;
  protected _innerProductMatrix: unknown;

  constructor(
    baseRing: RingLike,
    degree: number,
    sparse: boolean = false,
    innerProductMatrix?: unknown
  ) {
    if (degree < 0) {
      throw new ValueError(`degree (=${degree}) must be nonnegative`);
    }
    this._baseRing = baseRing;
    this._degree = degree;
    this._sparse = sparse;
    this._innerProductMatrix = innerProductMatrix;
  }

  /**
   * Return the degree of this free module. This is the dimension of the
   * ambient vector space in which it is embedded.
   */
  degree(): number {
    return this._degree;
  }

  /**
   * Return the base ring of this module.
   */
  baseRing(): RingLike {
    return this._baseRing;
  }

  /**
   * Return whether this module uses sparse representation.
   */
  isSparse(): boolean {
    return this._sparse;
  }

  /**
   * Return whether elements are represented exactly.
   */
  isExact(): boolean {
    if ('is_exact' in this._baseRing && typeof this._baseRing.is_exact === 'function') {
      return this._baseRing.is_exact();
    }
    return true; // Assume exact by default
  }

  /**
   * Return the inner product matrix for this module.
   */
  innerProductMatrix(): unknown {
    return this._innerProductMatrix;
  }

  /**
   * Create an element of this module from entries.
   */
  createElement(entries: unknown[]): FreeModuleElement {
    // Coerce entries to the base ring if needed.  An element of the fraction
    // field of a Euclidean base ring (e.g. QQ(x) for QQ[x]) is left alone: it
    // lives in the coordinate ring, not in the base ring.
    const coercedEntries = entries.map((e) => {
      if (e instanceof FractionFieldElement) {
        return e;
      }
      if (this._baseRing.__call__) {
        return this._baseRing.__call__(e);
      }
      return e;
    });

    if (this._sparse) {
      return new FreeModuleElementSparse(this, coercedEntries);
    }
    return new FreeModuleElementDense(this, coercedEntries);
  }

  /**
   * Return the zero vector in this module.
   */
  zeroVector(): FreeModuleElement {
    const zero = this._baseRing.zero();
    const entries = new Array(this._degree).fill(zero);

    if (this._sparse) {
      return new FreeModuleElementSparse(this, new Map());
    }
    return new FreeModuleElementDense(this, entries);
  }

  /**
   * Return the zero submodule of this module.
   */
  zeroSubmodule(): FreeModuleGeneric {
    return this.span([], this._baseRing);
  }

  /**
   * Return the R-span of gens, where R is the base_ring.
   * @param gens - List of vectors
   * @param baseRing - Optional base ring
   * @param options - Additional options
   */
  span(
    gens: FreeModuleElement[],
    baseRing?: RingLike,
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleGeneric {
    const ring = baseRing ?? this._baseRing;
    const self = this as unknown as FreeModuleGeneric;

    // The span lives in the ambient module, not in self
    // (free_module.py:1586: self._submodule_class(self.ambient_module(), ...)).
    const ambient = self.ambientModule ? self.ambientModule() : self;

    if (ring !== this._baseRing) {
      // The base ring changed: re-span in the ambient module over the new ring
      const M = ambient.changeRing(ring);
      return M.span(
        gens.map((g) => M.createElement(g.list())),
        ring,
        options
      );
    }

    const opts = {
      check: options?.check ?? true,
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    };

    if (isField(ring)) {
      return new FreeModuleSubspace(ambient as FreeModuleField, gens, opts);
    }
    if (arithmeticFor(ring).exact) {
      return new FreeModuleSubmodulePID(ambient, gens, opts);
    }
    // Over a general ring the generators are stored verbatim, exactly as in
    // SageMath's Submodule_free_ambient.
    return new FreeModuleSubmodule(ambient, gens, opts);
  }

  /**
   * Create the R-submodule of the ambient module with given generators.
   * @param gens - List of vectors or a free module
   * @param options - Additional options
   */
  submodule(
    gens: FreeModuleElement[] | FreeModuleGeneric,
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleGeneric {
    const list = Array.isArray(gens) ? gens : gens.gens();
    const V = this.span(list, this._baseRing, options);

    if (options?.check ?? true) {
      if (!V.isSubmodule(this)) {
        throw new ArithmeticError(
          `argument gens (= ${list.map((g) => g.toString()).join(', ')}) does not generate a submodule of self`
        );
      }
    }

    return V;
  }

  /**
   * Return the quotient of self by the given submodule.
   * @param sub - A submodule of self
   * @param check - Whether to check that sub is a submodule
   *
   * @see Reference: sage/modules/free_module.py:Module_free_ambient.quotient_module
   */
  quotientModule(sub: FreeModuleGeneric, check: boolean = true): FreeModuleQuotient {
    if (this._baseRing !== sub.baseRing()) {
      throw new ValueError('base rings must be the same');
    }
    // Check that sub is a valid submodule
    if (check) {
      if (!sub.isSubmodule(this as unknown as ModuleFreeAmbient)) {
        throw new ArithmeticError('sub must be a subspace of self');
      }
    }

    return new FreeModuleQuotient(this as unknown as FreeModuleGeneric, sub);
  }

  /**
   * Return whether self is a submodule of other.
   * @param other - Another module
   *
   * @see Reference: sage/modules/free_module.py:Module_free_ambient.is_submodule
   */
  isSubmodule(other: ModuleFreeAmbient): boolean {
    if (this === other) {
      return true;
    }

    if (this._baseRing !== other._baseRing) {
      return false;
    }

    if (this._degree !== other._degree) {
      return false;
    }

    const self = this as unknown as FreeModuleGeneric;
    const target = other as unknown as FreeModuleGeneric;

    if (target.rank() < self.rank()) {
      return false;
    }

    // The zero module is always a submodule
    const selfBasis = self.basis();
    if (selfBasis.length === 0) {
      return true;
    }

    const otherBasis = target.basis();
    if (otherBasis.length === 0) {
      return false;
    }

    // Solve  self.basis_matrix() = M * other.basis_matrix()  and require every
    // entry of M to lie in the base ring (free_module.py:2287).
    const ar = arithmeticFor(this._baseRing);
    if (!ar.exact) {
      throw new NotImplementedError(
        'could not determine whether this is a submodule over this base ring'
      );
    }
    const B = liftRows(
      otherBasis.map((b) => b.list()),
      ar
    );

    for (const gen of selfBasis) {
      const v = gen.list().map((e) => ar.lift(e));
      const x = solveLeftLifted(B, v, ar);
      if (x === null) {
        return false;
      }
      for (const c of x) {
        if (!ar.inBaseRing(c)) {
          return false;
        }
      }
    }
    return true;
  }
}

/**
 * Base class for all free modules.
 * @see Reference: sage/modules/free_module.py:FreeModule_generic
 */
export class FreeModuleGeneric extends ModuleFreeAmbient {
  protected _rank: number;
  protected _coordinateRing: RingLike;
  protected _basis: FreeModuleElement[] | null = null;

  constructor(
    baseRing: RingLike,
    rank: number,
    degree: number,
    sparse: boolean = false,
    coordinateRing?: RingLike,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, degree, sparse, innerProductMatrix);
    this._rank = rank;
    this._coordinateRing = coordinateRing ?? baseRing;
  }

  /**
   * Return the rank of this free module.
   */
  rank(): number {
    return this._rank;
  }

  /**
   * Return the dimension of this free module (same as rank).
   */
  dimension(): number {
    return this._rank;
  }

  /**
   * Return the codimension of this free module.
   */
  codimension(): number {
    return this._degree - this._rank;
  }

  /**
   * Return the basis of this module.
   */
  basis(): FreeModuleElement[] {
    if (this._basis !== null) {
      return this._basis;
    }

    // Default: standard basis
    const zero = this._baseRing.zero();
    const one = this._baseRing.one();
    const basisVectors: FreeModuleElement[] = [];

    for (let i = 0; i < this._rank; i++) {
      const entries: unknown[] = new Array(this._degree).fill(zero);
      entries[i] = one;
      const v = this.createElement(entries);
      v.setImmutable();
      basisVectors.push(v);
    }

    this._basis = basisVectors;
    return basisVectors;
  }

  /**
   * Return a tuple of basis elements.
   */
  gens(): FreeModuleElement[] {
    return this.basis();
  }

  /**
   * Return the i-th generator.
   * @param i - The index (default: 0)
   */
  gen(i: number = 0): FreeModuleElement {
    const b = this.basis();
    if (i < 0 || i >= b.length) {
      throw new ValueError(`generator index ${i} out of range [0, ${b.length})`);
    }
    return b[i]!;
  }

  /**
   * Return the number of basis elements.
   */
  ngens(): number {
    return this._rank;
  }

  /**
   * Return the matrix whose rows are the basis for this free module.
   * @param ring - Optional ring for the matrix
   */
  basisMatrix(ring?: RingLike): unknown {
    // For ambient modules, this is the identity matrix
    // For submodules, it's the matrix whose rows are the basis vectors
    const b = this.basis();
    const entries: unknown[][] = [];

    for (const v of b) {
      entries.push(v.list());
    }

    return entries;
  }

  /**
   * Return the echelonized basis matrix.
   *
   * The echelonized basis matrix is the row echelon form of the basis matrix.
   * For ambient modules, this is the identity matrix.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.echelonized_basis_matrix
   */
  echelonizedBasisMatrix(): unknown[][] {
    if (this.isAmbient()) {
      // For ambient modules, the echelonized basis is the identity
      const n = this._rank;
      const result: unknown[][] = [];
      const zero = this._baseRing.zero();
      const one = this._baseRing.one();

      for (let i = 0; i < n; i++) {
        const row: unknown[] = [];
        for (let j = 0; j < n; j++) {
          row.push(i === j ? one : zero);
        }
        result.push(row);
      }
      return result;
    }

    // For submodules: the Hermite normal form over ZZ, the reduced row
    // echelon form over a field.
    const basisMat = this.basisMatrix() as unknown[][];
    return echelonRows(basisMat, arithmeticFor(this._baseRing));
  }

  /**
   * Return the echelonized basis of this module.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid.echelonized_basis
   */
  echelonizedBasis(): FreeModuleElement[] {
    const ambient = this.ambientModule();
    return this.echelonizedBasisMatrix().map((row) => {
      const v = ambient.createElement(row);
      v.setImmutable();
      return v;
    });
  }

  /**
   * Return whether this module equals other, i.e. whether they have the same
   * ambient space and the same echelonized basis.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic._eq
   */
  equals(other: FreeModuleGeneric): boolean {
    if (this === other) {
      return true;
    }
    if (this._degree !== other.degree() || this._baseRing !== other.baseRing()) {
      return false;
    }
    if (this.rank() !== other.rank()) {
      return false;
    }
    const A = this.echelonizedBasisMatrix();
    const B = other.echelonizedBasisMatrix();
    if (A.length !== B.length) {
      return false;
    }
    const ar = arithmeticFor(this._baseRing);
    for (let i = 0; i < A.length; i++) {
      for (let j = 0; j < this._degree; j++) {
        if (!ar.eq(ar.lift(A[i]![j]), ar.lift(B[i]![j]))) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Return the basis matrix of this module.
   */
  matrix(): unknown {
    return this.basisMatrix();
  }

  /**
   * Return the Gram matrix B*A*B^T where A is the inner product matrix
   * and B is the basis matrix. If A is the identity (standard inner product),
   * this is just B*B^T.
   * @returns A 2D array representing the Gram matrix
   */
  gramMatrix(): unknown[][] {
    const zero = this._baseRing.zero();
    const one = this._baseRing.one();

    if (this.isAmbient()) {
      // The Gram matrix of an ambient module is its inner product matrix
      // (the identity when there is none).
      const A = this.innerProductMatrix();
      if (Array.isArray(A)) {
        return (A as unknown[][]).map((row) => [...row]);
      }
      const G: unknown[][] = [];
      for (let i = 0; i < this._degree; i++) {
        const row: unknown[] = [];
        for (let j = 0; j < this._degree; j++) {
          row.push(i === j ? one : zero);
        }
        G.push(row);
      }
      return G;
    }

    // G = B*A*B^t, where A is the inner product matrix of the ambient module
    // and B the basis matrix; the inner product of the basis vectors already
    // applies A.
    const b = this.basis();
    const n = b.length;
    const G: unknown[][] = [];
    for (let i = 0; i < n; i++) {
      const row: unknown[] = [];
      for (let j = 0; j < n; j++) {
        row.push(b[i]!.innerProduct(b[j]!));
      }
      G.push(row);
    }
    return G;
  }

  /**
   * Return the discriminant of this free module.
   *
   * This is the determinant of the Gram matrix.  When the module carries an
   * inner product matrix it is a free quadratic module, whose discriminant
   * carries the extra sign `(-1)^(rank/2)`.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.discriminant
   * @see Reference: sage/modules/free_quadratic_module.py:FreeQuadraticModule_generic.discriminant
   */
  discriminant(): unknown {
    const G = this.gramMatrix();
    const n = G.length;

    if (n === 0) {
      return this._baseRing.one();
    }

    const ar = arithmeticFor(this._baseRing);
    let det = determinantLifted(liftRows(G, ar), ar);

    if (this._innerProductMatrix !== null && this._innerProductMatrix !== undefined) {
      const r = Math.floor(this.rank() / 2);
      if (r % 2 === 1) {
        det = ar.neg(det);
      }
    }

    return ar.lower(det);
  }

  /**
   * Return the cardinality of this module.
   *
   * @returns The number of elements in this module, or Infinity if infinite.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.cardinality
   * @see Deviation: finite cardinalities are bigints (they routinely exceed
   *   2^53); the infinite cardinality is the JavaScript Infinity.
   */
  cardinality(): bigint | number {
    if (this._rank === 0) {
      return 1n;
    }

    // Check if base ring is finite
    if ('cardinality' in this._baseRing && typeof this._baseRing.cardinality === 'function') {
      const baseCard = this._baseRing.cardinality();
      if (typeof baseCard === 'bigint') {
        return baseCard ** BigInt(this._rank);
      }
      if (typeof baseCard === 'number') {
        if (!Number.isFinite(baseCard)) {
          return Number.POSITIVE_INFINITY;
        }
        return BigInt(baseCard) ** BigInt(this._rank);
      }
      if (typeof baseCard === 'object' && baseCard !== null) {
        const value = (baseCard as { value?: unknown }).value;
        if (typeof value === 'bigint') {
          return value ** BigInt(this._rank);
        }
      }
    }

    // Check if base ring has is_finite method
    if ('is_finite' in this._baseRing && typeof this._baseRing.is_finite === 'function') {
      if (this._baseRing.is_finite()) {
        // Finite but cardinality not directly available
        throw new NotImplementedError('cardinality not computable for this finite ring');
      }
    }

    // Infinite ring means infinite module (unless rank is 0)
    return Number.POSITIVE_INFINITY;
  }

  /**
   * Return whether this is an ambient module.
   */
  isAmbient(): boolean {
    return false;
  }

  /**
   * Return whether the inner product on this module is the one induced by the
   * ambient inner product.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.uses_ambient_inner_product
   */
  usesAmbientInnerProduct(): boolean {
    return true;
  }

  /**
   * Return whether this module is dense.
   */
  isDense(): boolean {
    return !this._sparse;
  }

  /**
   * Return whether the rank equals the degree.
   */
  isFull(): boolean {
    return this._rank === this._degree;
  }

  /**
   * Return whether this module is finite.
   */
  isFinite(): boolean {
    if (this._rank === 0) {
      return true;
    }
    // Check if base ring is finite
    if ('is_finite' in this._baseRing) {
      return (this._baseRing as { is_finite: () => boolean }).is_finite();
    }
    return false;
  }

  /**
   * Return whether the given basis has been specified by the user.
   */
  hasUserBasis(): boolean {
    return false;
  }

  /**
   * Return the coordinate ring of this module.
   */
  coordinateRing(): RingLike {
    return this._coordinateRing;
  }

  /**
   * Return the base field (fraction field of the base ring).
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.base_field
   */
  baseField(): RingLike {
    if (isField(this._baseRing)) {
      return this._baseRing;
    }

    // Try to get the fraction field
    if ('fraction_field' in this._baseRing && typeof this._baseRing.fraction_field === 'function') {
      return this._baseRing.fraction_field();
    }

    // For ZZ, the fraction field is QQ
    if (this._baseRing.toString?.() === 'Integer Ring') {
      return {
        zero: () => 0,
        one: () => 1,
        is_field: () => true,
        toString: () => 'Rational Field',
      };
    }

    // For a Euclidean domain such as QQ[x] the fraction field is the field of
    // fractions num/den (QQ(x)).
    if (isEuclideanRing(this._baseRing)) {
      return fractionFieldOf(this._baseRing);
    }

    throw new NotImplementedError('fraction_field not available for this base ring');
  }

  /**
   * Return the ambient module.
   */
  ambientModule(): FreeModuleGeneric {
    return this;
  }

  /**
   * Return the ambient vector space.
   *
   * This is the vector space over the fraction field that contains this module.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.ambient_vector_space
   */
  ambientVectorSpace(): FreeModuleField {
    const field = this.baseField();
    return new FreeModuleAmbientField(field, this._degree, this._sparse, this._innerProductMatrix);
  }

  /**
   * Write v in terms of the basis for self.
   *
   * The coordinates are computed exactly in the fraction field of the base
   * ring, so they may be rational even when the module is defined over ZZ
   * (SageMath returns them in `FreeModule(R.fraction_field(), rank)`).
   *
   * @param v - A vector
   * @param check - Whether to verify v is in self
   * @returns The list of coefficients c such that v = sum(c[i] * basis[i])
   * @throws {ArithmeticError} If v is not in the span of the basis
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.coordinates
   * @see Deviation: an integral coordinate is returned as a bigint and a
   *   non-integral one as a Rational; `check: false` still raises when v is
   *   outside the span, where SageMath returns a meaningless vector.
   */
  coordinates(v: FreeModuleElement, check: boolean = true): unknown[] {
    const basis = this.basis();
    const n = basis.length;

    if (v.degree() !== this._degree) {
      throw new ArithmeticError('vector is not in free module');
    }

    if (n === 0) {
      if (v.isZero()) {
        return [];
      }
      throw new ArithmeticError('vector is not in free module');
    }

    // For ambient modules the coordinates are just the entries
    if (this.isAmbient() && this._rank === this._degree) {
      return v.list();
    }

    // Solve  x * B = v  exactly over the fraction field of the base ring.
    const ar = arithmeticFor(this._baseRing);
    if (!ar.exact) {
      throw new NotImplementedError('coordinates are not implemented over this base ring');
    }
    const B = liftRows(
      basis.map((b) => b.list()),
      ar
    );
    const target = v.list().map((e) => ar.lift(e));

    // The exact solve fails exactly when v is not in the span of the basis,
    // which is the condition SageMath's `check` verifies.  We therefore always
    // raise, even when check is false.
    const x = solveLeftLifted(B, target, ar);
    if (x === null) {
      throw new ArithmeticError('vector is not in free module');
    }

    return x.map((c) => ar.lower(c));
  }

  /**
   * Return the coordinate vector for v with respect to the basis.
   * @param v - A vector
   * @param check - Whether to verify v is in self
   */
  coordinateVector(v: FreeModuleElement, check?: boolean): FreeModuleElement {
    const coords = this.coordinates(v, check);

    // Create a new free module for the coordinate space
    const coordModule = FreeModule(this._baseRing, coords.length);
    return coordModule.createElement(coords);
  }

  /**
   * Return the direct sum of self and other.
   * @param other - Another free module
   * @returns A new module of rank self.rank() + other.rank()
   */
  directSum(other: FreeModuleGeneric): FreeModuleGeneric {
    // Direct sum has rank = rank(self) + rank(other)
    // and degree = degree(self) + degree(other) (for embedded direct sum)
    // The basis vectors are (v, 0) for v in basis(self) and (0, w) for w in basis(other)

    const newRank = this._rank + other.rank();
    const newDegree = this._degree + other.degree();

    const newModule = FreeModule(this._baseRing, newDegree, { sparse: this._sparse });

    // Build basis vectors
    const zero = this._baseRing.zero();
    const basisVectors: FreeModuleElement[] = [];

    // Basis vectors from self: (v, 0)
    for (const v of this.basis()) {
      const entries: unknown[] = [];
      for (let i = 0; i < this._degree; i++) {
        entries.push(v.getItem(i));
      }
      for (let i = 0; i < other.degree(); i++) {
        entries.push(zero);
      }
      basisVectors.push(newModule.createElement(entries));
    }

    // Basis vectors from other: (0, w)
    for (const w of other.basis()) {
      const entries: unknown[] = [];
      for (let i = 0; i < this._degree; i++) {
        entries.push(zero);
      }
      for (let i = 0; i < other.degree(); i++) {
        entries.push(w.getItem(i));
      }
      basisVectors.push(newModule.createElement(entries));
    }

    return newModule.span(basisVectors);
  }

  /**
   * Return the product of this module by a scalar.
   * If M is a module with basis b_1, ..., b_n, then scale(c) has basis c*b_1, ..., c*b_n.
   * @param scalar - A scalar from the base ring
   * @returns A new module with scaled basis
   */
  scale(scalar: unknown): FreeModuleGeneric {
    const b = this.basis();
    const scaledVectors: FreeModuleElement[] = [];

    for (const v of b) {
      scaledVectors.push(v.mul(scalar));
    }

    return this.span(scaledVectors);
  }

  /**
   * Return a random element of this module.
   * @param prob - Probability each coefficient is non-zero (default: 1.0)
   * @param min - Minimum value for random coefficients (default: -10)
   * @param max - Maximum value for random coefficients (default: 10)
   */
  randomElement(prob: number = 1.0, min: number = -10, max: number = 10): FreeModuleElement {
    const entries: unknown[] = [];
    const zero = this._baseRing.zero();

    for (let i = 0; i < this._degree; i++) {
      if (Math.random() < prob) {
        // Generate random value
        if (typeof zero === 'bigint') {
          const val = BigInt(Math.floor(Math.random() * (max - min + 1)) + min);
          entries.push(val);
        } else if (typeof zero === 'number') {
          entries.push(Math.floor(Math.random() * (max - min + 1)) + min);
        } else if (this._baseRing.__call__) {
          const val = Math.floor(Math.random() * (max - min + 1)) + min;
          entries.push(this._baseRing.__call__(val));
        } else {
          entries.push(zero);
        }
      } else {
        entries.push(zero);
      }
    }

    return this.createElement(entries);
  }

  /**
   * Return whether the given vectors are linearly dependent.
   * @param vecs - A list of vectors
   */
  areLinearlyDependent(vecs: FreeModuleElement[]): boolean {
    if (vecs.length === 0) {
      return false; // Empty set is linearly independent
    }

    if (vecs.length > this._degree) {
      return true; // More vectors than the degree means dependent
    }

    // A = matrix(vecs); A.echelonize(); any zero row means dependence.
    const ar = arithmeticFor(this._baseRing);
    if (!ar.exact) {
      throw new NotImplementedError('linear dependence is not implemented over this base ring');
    }
    return (
      rankOfRows(
        vecs.map((v) => v.list()),
        ar
      ) < vecs.length
    );
  }

  /**
   * Return a linear combination of the basis vectors.
   * @param coefficients - The coefficients
   */
  linearCombinationOfBasis(coefficients: unknown[]): FreeModuleElement {
    const b = this.basis();
    if (coefficients.length !== b.length) {
      throw new ValueError(`coefficients must have length ${b.length}, got ${coefficients.length}`);
    }

    let result = this.zeroVector();
    for (let i = 0; i < b.length; i++) {
      result = result.add(b[i]!.mul(coefficients[i]));
    }
    return result;
  }

  /**
   * Return the dense version of this module.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.dense_module
   */
  denseModule(): FreeModuleGeneric {
    if (!this._sparse) {
      return this;
    }

    // Create a dense ambient module and span the same vectors
    const denseAmbient = FreeModule(this._baseRing, this._degree, { sparse: false });
    const basis = this.basis();
    const denseBasis: FreeModuleElement[] = [];

    for (const v of basis) {
      denseBasis.push(denseAmbient.createElement(v.list()));
    }

    return denseAmbient.span(denseBasis);
  }

  /**
   * Return the sparse version of this module.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.sparse_module
   */
  sparseModule(): FreeModuleGeneric {
    if (this._sparse) {
      return this;
    }

    // Create a sparse ambient module and span the same vectors
    const sparseAmbient = FreeModule(this._baseRing, this._degree, { sparse: true });
    const basis = this.basis();
    const sparseBasis: FreeModuleElement[] = [];

    for (const v of basis) {
      sparseBasis.push(sparseAmbient.createElement(v.list()));
    }

    return sparseAmbient.span(sparseBasis);
  }

  /**
   * Return the free module over `ring` obtained by coercing each element of
   * the basis of self into a vector over the fraction field of `ring`, then
   * taking the resulting module.
   *
   * @param ring - A principal ideal domain
   * @throws {TypeError} If the new ring is not a principal ideal domain
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid.change_ring
   * @see Reference: sage/modules/free_module.py:FreeModule_ambient.change_ring
   */
  changeRing(ring: RingLike): FreeModuleGeneric {
    if (ring === this._baseRing) {
      return this;
    }
    if (!isField(ring) && !isPID(ring)) {
      throw new TypeError(
        `the new ring ${ring.toString?.() ?? ring} should be a principal ideal domain`
      );
    }

    if (this.isAmbient()) {
      return FreeModule(ring, this._rank, { sparse: this._sparse });
    }

    // Re-span the basis, in the ambient module of the same degree, over R.
    const M = this.ambientModule().changeRing(ring);
    const B = this.basis().map((b) => M.createElement(coerceRow(ring, b.list())));
    if (this.hasUserBasis() && M instanceof FreeModulePID) {
      return M.spanOfBasis(B, ring);
    }
    return M.span(B, ring);
  }

  /**
   * Return an ambient free module isomorphic to this one.
   */
  nonembeddedFreeModule(): FreeModuleGeneric {
    return FreeModule(this._baseRing, this._rank, { sparse: this._sparse });
  }

  /**
   * String representation.
   */
  toString(): string {
    const ringName = this._baseRing.toString?.() ?? 'Ring';
    if (this._degree === this._rank) {
      if (isField(this._baseRing)) {
        return `Vector space of dimension ${this._rank} over ${ringName}`;
      }
      return `Free module of rank ${this._rank} over ${ringName}`;
    }
    if (isField(this._baseRing)) {
      return `Vector space of degree ${this._degree} and dimension ${this._rank} over ${ringName}`;
    }
    return `Free module of degree ${this._degree} and rank ${this._rank} over ${ringName}`;
  }
}

/**
 * Base class for free modules over an integral domain.
 * @see Reference: sage/modules/free_module.py:FreeModule_generic_domain
 */
export class FreeModuleDomain extends FreeModuleGeneric {
  constructor(
    baseRing: RingLike,
    rank: number,
    degree: number,
    sparse: boolean = false,
    coordinateRing?: RingLike,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, rank, degree, sparse, coordinateRing, innerProductMatrix);
  }

  /**
   * Return the sum of self and other.
   * @param other - Another submodule
   */
  add(other: FreeModuleGeneric): FreeModuleGeneric {
    // Combine generators from both modules
    const gens = [...this.gens(), ...other.gens()];
    return this.span(gens);
  }
}

/**
 * Base class for free modules over a PID.
 * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid
 */
export class FreeModulePID extends FreeModuleDomain {
  constructor(
    baseRing: RingLike,
    rank: number,
    degree: number,
    sparse: boolean = false,
    coordinateRing?: RingLike,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, rank, degree, sparse, coordinateRing, innerProductMatrix);
  }

  /**
   * Return the lattice index [other:self] of self in other.
   *
   * When self is contained in other, the lattice index is the usual index.
   * If the index is infinite, this returns Infinity.
   *
   * @param other - Another module
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.index_in
   * @see Deviation: an integral index is returned as a bigint and a
   *   non-integral one as a Rational.
   */
  indexIn(other: FreeModuleGeneric): unknown {
    if (this._baseRing !== other.baseRing()) {
      throw new NotImplementedError(
        'lattice index only defined for modules over the same base ring.'
      );
    }
    if (this._degree !== other.degree()) {
      throw new ArithmeticError('self and other must be embedded in the same ambient space.');
    }

    const ar = arithmeticFor(this._baseRing);

    if (ar.isField) {
      if (this.equals(other)) {
        return 1n;
      }
      if (this.isSubmodule(other as unknown as ModuleFreeAmbient)) {
        return Number.POSITIVE_INFINITY;
      }
      throw new ArithmeticError('self must be contained in the vector space spanned by other.');
    }

    // C = [other.coordinates(b) for b in self.basis()]
    const C: unknown[][] = [];
    for (const b of this.basis()) {
      C.push(other.coordinates(b).map((c) => ar.lift(c)));
    }

    if (this.rank() < other.rank()) {
      return Number.POSITIVE_INFINITY;
    }

    const det = determinantLifted(C, ar);
    // For ZZ the index is the absolute value of the determinant
    const r = det as Rational;
    const abs = ar.isIntegral && r instanceof Rational ? r.abs() : det;
    return ar.lower(abs);
  }

  /**
   * Return the intersection of self and other.
   *
   * @param other - Another module
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.intersection
   */
  intersection(other: FreeModuleGeneric): FreeModuleGeneric {
    if (this._degree !== other.degree()) {
      throw new ArithmeticError('self and other must be embedded in the same ambient space.');
    }

    // Dispense with the easy cases
    if (this.rank() === 0 || other.rank() === 0) {
      return this.zeroSubmodule();
    }
    if (this === other) {
      return this;
    }
    if (other.isSubmodule(this as unknown as ModuleFreeAmbient)) {
      return other;
    }
    if (this.isSubmodule(other as unknown as ModuleFreeAmbient)) {
      return this;
    }

    // Standard algorithm: let S be A1 stacked on A2; the vectors v in the
    // (left) kernel of S give the intersection as (v[:n]) * A1.
    const [V1, V2] =
      this.rank() <= other.rank()
        ? [this as unknown as FreeModuleGeneric, other]
        : [other, this as unknown as FreeModuleGeneric];

    const ar = arithmeticFor(this._baseRing);
    const A1 = V1.basis().map((v) => v.list());
    const A2 = V2.basis().map((v) => v.list());
    const S = [...A1, ...A2];
    const n = A1.length;

    let K: unknown[][];
    if (ar.isEuclidean) {
      // integer_kernel over a Euclidean domain: read the kernel off the
      // Hermite transformation matrix (matrix2.pyx:5641).
      K = euclideanLeftKernelRows(S, ar as EuclideanArithmetic);
    } else {
      // Left kernel of S = right kernel of S^t
      const St: unknown[][] = [];
      for (let j = 0; j < this._degree; j++) {
        St.push(S.map((row) => row[j]));
      }
      K = rightKernelRows(St, S.length, ar);

      if (!ar.isField && ar.isIntegral && K.length > 0) {
        // integer_kernel: clear denominators and saturate, so that the kernel
        // is the full ZZ-module of integral relations.
        K = integralKernelRows(K, ar);
      }
    }

    const gens: FreeModuleElement[] = [];
    const ambient = this.ambientModule();
    for (const v of K) {
      const coeffs = v.slice(0, n).map((e) => ar.lift(e));
      const entries: unknown[] = [];
      for (let j = 0; j < this._degree; j++) {
        let acc = ar.zero();
        for (let i = 0; i < n; i++) {
          acc = ar.add(acc, ar.mul(coeffs[i], ar.lift(A1[i]![j])));
        }
        entries.push(ar.lower(acc));
      }
      const w = ambient.createElement(entries);
      if (!w.isZero()) {
        gens.push(w);
      }
    }

    if (gens.length === 0) {
      return this.zeroSubmodule();
    }
    return this.span(gens);
  }

  /**
   * Return the index of this module in its saturation.
   *
   * The saturation of a submodule M of a free module F is the largest
   * submodule of F containing M with the same rank.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.index_in_saturation
   */
  indexInSaturation(): unknown {
    return this.indexIn(this.saturation());
  }

  /**
   * Return the saturated submodule of R^n that spans the same vector space.
   *
   * The saturation of a submodule M of a free module F is the largest
   * submodule S of F containing M with the same rank. Equivalently,
   * S is the intersection of the vector space span of M with F.
   *
   * For a lattice L in ZZ^n, the saturation is obtained by computing
   * the Hermite normal form of the basis matrix and dividing by the GCD
   * of the entries.
   *
   * @returns The saturation of this module
   *
   * @example
   * ```typescript
   * // Create a non-saturated lattice
   * const L = span([[9, 9, 6]], ZZ);
   * // L.saturation() returns span([[3, 3, 2]], ZZ)
   * ```
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.saturation
   */
  saturation(): FreeModuleGeneric {
    // If the base ring is a field, the module is already saturated
    if (isField(this._baseRing)) {
      return this;
    }
    if (this._rank === 0) {
      return this;
    }

    const basisMat = this.basisMatrix() as unknown[][];
    if (basisMat.length === 0) {
      return this;
    }

    const ar = arithmeticFor(this._baseRing);
    if (!ar.isIntegral) {
      throw new NotImplementedError('saturation is only implemented over ZZ');
    }

    // A, _ = self.basis_matrix()._clear_denom(); S = self.span(A.saturation())
    const lifted = liftRows(basisMat, ar);
    let d = 1n;
    for (const row of lifted) {
      for (const e of row) {
        d = bigintLcm(d, ar.denominator(e) as bigint);
      }
    }
    const cleared: bigint[][] = lifted.map((row) =>
      row.map((e) => (e as Rational).mul(new Rational(d)).numerator)
    );

    const S = matrixSaturation(IntegerMatrixFromEntries(cleared));

    const ambient = this.ambientModule();
    const saturatedVectors: FreeModuleElement[] = [];
    for (let i = 0; i < S.nrows; i++) {
      const row: unknown[] = [];
      for (let j = 0; j < S.ncols; j++) {
        row.push(ar.lower(new Rational(S.get(i, j).value)));
      }
      saturatedVectors.push(ambient.createElement(row));
    }

    const sat = ambient.span(saturatedVectors);

    // Return exactly self if it is already saturated
    return this.equals(sat) ? this : sat;
  }

  /**
   * Return the denominator of the basis matrix.
   *
   * This is the LCM of the denominators of all entries in the basis matrix
   * when expressed in the ambient space coordinates.
   *
   * @returns The denominator (LCM of all entry denominators)
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.denominator
   */
  denominator(): unknown {
    const basisMat = this.basisMatrix() as unknown[][];

    if (basisMat.length === 0) {
      return this._baseRing.one();
    }

    const ar = arithmeticFor(this._baseRing);
    if (!ar.exact) {
      return this._baseRing.one();
    }

    let d = ar.denominatorOne();
    for (const row of basisMat) {
      for (const entry of row) {
        d = ar.denominatorLcm(d, ar.denominator(ar.lift(entry)));
      }
    }

    if (ar.isEuclidean) {
      return d;
    }
    return ar.lower(new Rational(d as bigint));
  }

  /**
   * Return the free R-module with the given basis.
   * @param basis - A list of vectors
   * @param baseRing - Optional base ring
   * @param options - Additional options
   */
  spanOfBasis(
    basis: FreeModuleElement[],
    baseRing?: RingLike,
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleGeneric {
    const ring = baseRing ?? this._baseRing;
    const ambient = this.ambientModule();

    if (ring !== this._baseRing) {
      const M = ambient.changeRing(ring);
      return (M as FreeModulePID).spanOfBasis(
        basis.map((b) => M.createElement(b.list())),
        ring,
        options
      );
    }

    const opts = {
      check: options?.check ?? true,
      echelonize: false,
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    };

    if (isField(ring)) {
      return new FreeModuleSubspaceWithBasis(ambient as FreeModuleField, basis, opts);
    }
    return new FreeModuleWithBasis(ambient, basis, opts);
  }

  /**
   * Create the R-submodule with given basis.
   * @param basis - A list of linearly independent vectors
   * @param options - Additional options
   */
  submoduleWithBasis(
    basis: FreeModuleElement[],
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleGeneric {
    return this.spanOfBasis(basis, this._baseRing, options);
  }

  /**
   * Create a vector subspace of the ambient vector space.
   *
   * This creates a vector space over the fraction field of the base ring
   * that is spanned by the given generators.
   *
   * @param gens - A list of vectors
   * @param check - Whether to check vectors are in ambient space
   * @returns A vector space over the fraction field
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.vector_space_span
   */
  vectorSpaceSpan(gens: FreeModuleElement[], check?: boolean): FreeModuleField {
    const field = this.baseField();
    const ambient = new FreeModuleAmbientField(
      field,
      this._degree,
      this._sparse,
      this._innerProductMatrix
    );

    if (gens.length === 0) {
      return ambient.subspace([]);
    }

    // Convert generators to vectors in the ambient field
    const fieldGens: FreeModuleElement[] = [];
    for (const gen of gens) {
      const entries = gen.list();
      fieldGens.push(ambient.createElement(entries));
    }

    return ambient.subspace(fieldGens, { check: check ?? true });
  }

  /**
   * Create a vector subspace with given basis.
   *
   * This creates a vector space over the fraction field of the base ring
   * with the given vectors as its basis.
   *
   * @param gens - A list of vectors (must be linearly independent)
   * @param check - Whether to check vectors are in ambient space
   * @returns A vector space with the given basis
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.vector_space_span_of_basis
   */
  vectorSpaceSpanOfBasis(gens: FreeModuleElement[], check?: boolean): FreeModuleField {
    const field = this.baseField();
    const ambient = new FreeModuleAmbientField(
      field,
      this._degree,
      this._sparse,
      this._innerProductMatrix
    );

    if (gens.length === 0) {
      return ambient.subspaceWithBasis([]);
    }

    // Convert generators to vectors in the ambient field
    const fieldGens: FreeModuleElement[] = [];
    for (const gen of gens) {
      const entries = gen.list();
      fieldGens.push(ambient.createElement(entries));
    }

    return ambient.subspaceWithBasis(fieldGens, { check: check ?? true });
  }

  /**
   * Return the tensor product `self (x)_R other`.
   *
   * The result is embedded in `R^(deg(self) * deg(other))` and its basis is
   * the family of elementary tensors `e_i (x) f_j` in lexicographic order,
   * i.e. the Kronecker product of the two basis matrices.  The coordinate of
   * `e_k (x) f_l` in the ambient module is `k * deg(other) + l`.
   *
   * With `discardBasis: true` the tensor product is returned with the standard
   * basis and the Kronecker product of the two Gram matrices as inner product
   * matrix, which is SageMath's `discard_basis=True`.
   *
   * SageMath's `free_module.py` has no `tensor_product`; the tensor product of
   * two free modules embedded in `R^n` is defined in
   * `free_quadratic_module_integer_symmetric.py:1343`, which is what this
   * ports (via `Matrix.tensor_product`).
   *
   * @param other - Another free module over the same base ring
   *
   * @see Reference: sage/modules/free_quadratic_module_integer_symmetric.py:1343
   * @see Reference: sage/matrix/matrix2.pyx:9983 (Matrix.tensor_product)
   */
  tensorProduct(other: FreeModuleGeneric, options?: { discardBasis?: boolean }): FreeModuleGeneric {
    if (this._baseRing !== other.baseRing()) {
      throw new ValueError('base rings must be the same');
    }
    const ar = arithmeticFor(this._baseRing);
    if (!ar.exact) {
      throw new NotImplementedError('tensor products are not implemented over this base ring');
    }

    const sparse = this._sparse && other.isSparse();

    if (options?.discardBasis) {
      // gram_matrix = self.gram_matrix().tensor_product(other.gram_matrix())
      const G = lowerRows(
        kroneckerProduct(liftRows(this.gramMatrix(), ar), liftRows(other.gramMatrix(), ar), ar),
        ar
      );
      return FreeModule(this._baseRing, this._rank * other.rank(), {
        sparse,
        innerProductMatrix: G,
      });
    }

    const n = this._degree;
    const m = other.degree();

    // inner_product_matrix = self.inner_product_matrix() (x) other.inner_product_matrix()
    const ip1 = this._innerProductMatrix;
    const ip2 = other.innerProductMatrix();
    let innerProductMatrix: unknown;
    if (Array.isArray(ip1) || Array.isArray(ip2)) {
      const A = Array.isArray(ip1) ? (ip1 as unknown[][]) : identityRows(n, this._baseRing);
      const B = Array.isArray(ip2) ? (ip2 as unknown[][]) : identityRows(m, this._baseRing);
      innerProductMatrix = lowerRows(kroneckerProduct(liftRows(A, ar), liftRows(B, ar), ar), ar);
    }

    const ambient = FreeModule(this._baseRing, n * m, { sparse, innerProductMatrix });

    // basis_matrix = self.basis_matrix().tensor_product(other.basis_matrix())
    const rows = lowerRows(
      kroneckerProduct(
        liftRows(this.basisMatrix() as unknown[][], ar),
        liftRows(other.basisMatrix() as unknown[][], ar),
        ar
      ),
      ar
    );

    if (rows.length === 0) {
      return ambient.span([]);
    }
    const basis = rows.map((row) => ambient.createElement(row));
    return (ambient as FreeModulePID).spanOfBasis(basis, this._baseRing, { check: false });
  }
}

/**
 * The elementary tensor `v (x) w`, i.e. the Kronecker product of two vectors.
 *
 * If `b_0, ..., b_{r-1}` is the basis of `M` and `c_0, ..., c_{s-1}` that of
 * `N`, then `tensorProductVector(b_i, c_j)` is the `(i*s + j)`-th basis vector
 * of `M.tensorProduct(N)`.
 *
 * @see Reference: sage/matrix/matrix2.pyx:9983 (Matrix.tensor_product)
 */
export function tensorProductVector(
  v: FreeModuleElement,
  w: FreeModuleElement,
  parent?: FreeModuleGeneric
): FreeModuleElement {
  const a = v.list();
  const b = w.list();
  const ring = (v.parent() as FreeModuleGeneric).baseRing();
  const ar = arithmeticFor(ring);
  const entries: unknown[] = [];
  for (const x of a) {
    for (const y of b) {
      entries.push(ar.lower(ar.mul(ar.lift(x), ar.lift(y))));
    }
  }
  const M = parent ?? FreeModule(ring, a.length * b.length);
  return M.createElement(entries);
}

/**
 * The `n x n` identity matrix over the given ring, as rows.
 */
function identityRows(n: number, ring: RingLike): unknown[][] {
  const zero = ring.zero();
  const one = ring.one();
  const out: unknown[][] = [];
  for (let i = 0; i < n; i++) {
    const row: unknown[] = [];
    for (let j = 0; j < n; j++) {
      row.push(i === j ? one : zero);
    }
    out.push(row);
  }
  return out;
}

/**
 * Base class for vector spaces (free modules over a field).
 * @see Reference: sage/modules/free_module.py:FreeModule_generic_field
 */
export class FreeModuleField extends FreeModulePID {
  constructor(
    baseRing: RingLike,
    rank: number,
    degree: number,
    sparse: boolean = false,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, rank, degree, sparse, undefined, innerProductMatrix);
  }

  /**
   * Return the vector space of which this is a subspace.
   */
  vectorSpace(): FreeModuleField {
    return this;
  }

  /**
   * Return the subspace spanned by gens.
   * @param gens - A list of vectors
   * @param options - Additional options
   */
  subspace(
    gens: FreeModuleElement[],
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleField {
    return new FreeModuleSubspace(this.ambientModule() as FreeModuleField, gens, {
      check: options?.check ?? true,
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    });
  }

  /**
   * Return the subspace with given basis.
   * @param gens - A list of linearly independent vectors
   * @param options - Additional options
   */
  subspaceWithBasis(
    gens: FreeModuleElement[],
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleField {
    return new FreeModuleSubspaceWithBasis(this.ambientModule() as FreeModuleField, gens, {
      check: options?.check ?? true,
      echelonize: false,
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    });
  }

  /**
   * Return whether self is a subspace of other.
   *
   * @param other - Another vector space
   * @returns true if self is a subspace of other
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_field.is_subspace
   */
  isSubspace(other: FreeModuleField): boolean {
    // Use the underlying isSubmodule method
    return this.isSubmodule(other as unknown as ModuleFreeAmbient);
  }

  /**
   * Return the sum of self and other.
   * @param other - Another subspace
   */
  override add(other: FreeModuleField): FreeModuleField {
    const gens = [...this.gens(), ...other.gens()];
    return this.subspace(gens);
  }

  /**
   * Return the intersection of self and other.
   *
   * For vector spaces over a field, this uses the kernel-based algorithm.
   *
   * @param other - Another subspace
   * @returns The intersection of the two subspaces
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_field.intersection
   */
  override intersection(other: FreeModuleField): FreeModuleField {
    // Handle trivial cases
    if (this.dimension() === 0 || other.dimension() === 0) {
      return this.subspace([]);
    }

    if (this === other) {
      return this;
    }

    // If one is contained in the other, return the smaller one
    if (this.isSubspace(other)) {
      return this;
    }
    if (other.isSubspace(this)) {
      return other;
    }

    // Use the parent class intersection method (from FreeModulePID)
    const pidResult = super.intersection(other as unknown as FreeModuleGeneric);

    // Convert result to field subspace
    const gens = pidResult.gens();
    if (gens.length === 0) {
      return this.subspace([]);
    }

    // Create vectors in the ambient field space
    const ambient = this.ambientVectorSpace();
    const fieldGens: FreeModuleElement[] = [];

    for (const gen of gens) {
      fieldGens.push(ambient.createElement(gen.list()));
    }

    return ambient.subspace(fieldGens);
  }

  /**
   * Return the orthogonal complement of this subspace.
   *
   * The orthogonal complement is the set of all vectors v such that
   * <v, w> = 0 for all w in self.
   *
   * @returns The orthogonal complement of this subspace
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_field.orthogonal_complement
   */
  orthogonalComplement(): FreeModuleField {
    const basisMat = this.basisMatrix() as unknown[][];
    const n = this._degree;

    if (basisMat.length === 0) {
      // The complement of the zero space is the whole ambient space
      return this.ambientVectorSpace();
    }
    if (this.dimension() === n) {
      // The complement of the ambient space is zero
      return this.subspace([]);
    }

    // basis_matrix().right_kernel()
    const ar = arithmeticFor(this._baseRing);
    const kernel = rightKernelRows(basisMat, n, ar);

    const ambient = this.ambientVectorSpace();
    if (kernel.length === 0) {
      return ambient.subspace([]);
    }

    return ambient.subspace(kernel.map((row) => ambient.createElement(row)));
  }

  /**
   * Return the quotient `self/other`.
   *
   * This is the quotient *space*: an ambient vector space of dimension
   * `dim(self) - dim(other)` carrying the quotient map and a fixed section,
   * exactly as SageMath's `V/W`.
   *
   * @param other - A subspace of self
   * @returns The quotient space
   *
   * @see Reference: sage/modules/free_module.py:5239
   *   (FreeModule_generic_field.quotient_module)
   * @see Deviation: previously this returned an arbitrary complement of
   *   `other` inside `self`; that is a subspace of `self`, not the quotient,
   *   and it carried no quotient/lift maps.
   */
  quotient(other: FreeModuleField, check: boolean = true): FreeModuleQuotient {
    if (this._baseRing !== other.baseRing()) {
      throw new ValueError('base rings must be the same');
    }
    if (check && !other.isSubspace(this)) {
      throw new ArithmeticError('sub must be a subspace of self');
    }
    return new FreeModuleQuotient(this as unknown as FreeModuleGeneric, other);
  }

  /**
   * Return an iterator over subspaces of given dimension.
   *
   * Note: This is only implemented for finite fields where we can
   * enumerate all subspaces.
   *
   * @param dim - The dimension
   * @yields Subspaces of the given dimension
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_field.subspaces
   */
  *subspaces(dim: number): IterableIterator<FreeModuleField> {
    if (dim < 0 || dim > this.dimension()) {
      return;
    }

    if (dim === 0) {
      yield this.subspace([]);
      return;
    }

    if (dim === this.dimension()) {
      yield this;
      return;
    }

    // For general fields, enumerating all subspaces is complex
    // This requires enumerating Grassmannian, which needs more infrastructure
    throw new NotImplementedError(
      'subspaces iteration is only implemented for finite fields with enumerable elements'
    );
  }

  /**
   * Return a complement of this subspace.
   *
   * The complement is a subspace W such that self + W = ambient and
   * self intersection W = {0}.
   *
   * @returns A complement subspace
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_field.complement
   */
  complement(): FreeModuleField {
    // Simple cases
    if (this.dimension() === 0) {
      return this.ambientVectorSpace();
    }

    if (this.dimension() === this.ambientVectorSpace().dimension()) {
      return this.subspace([]);
    }

    // The orthogonal complement with respect to the standard inner product
    // is a complement (though not unique)
    return this.orthogonalComplement();
  }
}

// ============================================================================
// Ambient modules
// ============================================================================

/**
 * Ambient free module over a ring.
 * @see Reference: sage/modules/free_module.py:FreeModule_ambient
 */
export class FreeModuleAmbient extends FreeModuleGeneric {
  constructor(
    baseRing: RingLike,
    rank: number,
    sparse: boolean = false,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, rank, rank, sparse, undefined, innerProductMatrix);
  }

  /**
   * Return True since this is an ambient module.
   */
  override isAmbient(): boolean {
    return true;
  }
}

/**
 * Ambient free module over a PID.
 *
 * `FreeModule_ambient_pid` derives from `FreeModule_generic_pid` in SageMath,
 * so the full PID interface (span_of_basis, saturation, index_in, ...) is
 * available on ZZ^n.
 *
 * @see Reference: sage/modules/free_module.py:FreeModule_ambient_pid
 */
export class FreeModuleAmbientPID extends FreeModulePID {
  constructor(
    baseRing: RingLike,
    rank: number,
    sparse: boolean = false,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, rank, rank, sparse, undefined, innerProductMatrix);
  }

  /**
   * Return True since this is an ambient module.
   */
  override isAmbient(): boolean {
    return true;
  }
}

/**
 * Ambient vector space over a field.
 * @see Reference: sage/modules/free_module.py:FreeModule_ambient_field
 */
export class FreeModuleAmbientField extends FreeModuleField {
  constructor(
    baseRing: RingLike,
    dimension: number,
    sparse: boolean = false,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, dimension, dimension, sparse, innerProductMatrix);
  }

  /**
   * Return True since this is an ambient module.
   */
  override isAmbient(): boolean {
    return true;
  }
}

// ============================================================================
// Submodules
// ============================================================================

/**
 * Create a vector with the given entries, without coercing them.
 *
 * The entries produced by the exact linear algebra above already lie in the
 * coordinate ring of the module.
 */
function makeVector(parent: FreeModuleGeneric, entries: unknown[]): FreeModuleElement {
  const v = parent.isSparse()
    ? new FreeModuleElementSparse(parent, entries)
    : new FreeModuleElementDense(parent, entries);
  v.setImmutable();
  return v;
}

/**
 * Compute the user basis and the echelonized basis of a submodule.
 *
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid.__init__
 */
function submoduleBases(
  ambient: FreeModuleGeneric,
  basis: FreeModuleElement[],
  options?: { check?: boolean; echelonize?: boolean; alreadyEchelonized?: boolean }
): { user: unknown[][]; echelonized: unknown[][] | null } {
  const ar = arithmeticFor(ambient.baseRing());
  const rows = basis.map((b) => b.list());

  if (options?.alreadyEchelonized) {
    return { user: rows, echelonized: rows };
  }
  if (options?.echelonize) {
    const E = echelonRows(rows, ar);
    return { user: E, echelonized: E };
  }
  if ((options?.check ?? true) && ar.exact && rows.length > 0) {
    if (rankOfRows(rows, ar) !== rows.length) {
      throw new ValueError('the given basis vectors must be linearly independent');
    }
  }
  return { user: rows, echelonized: null };
}

/**
 * A submodule of a free module over a general ring.
 *
 * Over a ring that is not a PID no echelon form is available, so the
 * generators are stored verbatim, exactly as in SageMath's
 * `Submodule_free_ambient`.
 *
 * @see Reference: sage/modules/free_module.py:Submodule_free_ambient
 */
export class FreeModuleSubmodule extends FreeModuleGeneric {
  protected _ambient: FreeModuleGeneric;
  protected _userBasis: FreeModuleElement[];

  constructor(
    ambient: FreeModuleGeneric,
    gens: FreeModuleElement[],
    _options?: { check?: boolean; echelonize?: boolean; alreadyEchelonized?: boolean }
  ) {
    super(
      ambient.baseRing(),
      gens.length,
      ambient.degree(),
      ambient.isSparse(),
      undefined,
      ambient.innerProductMatrix()
    );

    this._ambient = ambient;
    this._userBasis = gens.map((g) => makeVector(this, g.list()));
    this._basis = this._userBasis;
  }

  override ambientModule(): FreeModuleGeneric {
    return this._ambient;
  }

  override hasUserBasis(): boolean {
    return true;
  }
}

/**
 * Submodule of a free module over a PID with a user-specified basis.
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid
 */
export class FreeModuleWithBasis extends FreeModulePID {
  protected _ambient: FreeModuleGeneric;
  protected _userBasis: FreeModuleElement[];
  protected _echelonizedBasisMatrix: unknown[][] | null = null;

  constructor(
    ambient: FreeModuleGeneric,
    basis: FreeModuleElement[],
    options?: { check?: boolean; echelonize?: boolean; alreadyEchelonized?: boolean }
  ) {
    const { user, echelonized } = submoduleBases(ambient, basis, options);
    super(
      ambient.baseRing(),
      user.length,
      ambient.degree(),
      ambient.isSparse(),
      undefined,
      ambient.innerProductMatrix()
    );

    this._ambient = ambient;
    this._userBasis = user.map((row) => makeVector(this, row));
    this._basis = this._userBasis;
    this._echelonizedBasisMatrix = echelonized;
  }

  /**
   * Return True since this module has a user-specified basis.
   */
  override hasUserBasis(): boolean {
    return true;
  }

  override ambientModule(): FreeModuleGeneric {
    return this._ambient;
  }

  /**
   * Return the basis matrix for self in row echelon form.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid.echelonized_basis_matrix
   */
  override echelonizedBasisMatrix(): unknown[][] {
    if (this._echelonizedBasisMatrix === null) {
      this._echelonizedBasisMatrix = echelonRows(
        this.basisMatrix() as unknown[][],
        arithmeticFor(this._baseRing)
      );
    }
    return this._echelonizedBasisMatrix;
  }
}

/**
 * An R-submodule of K^n where K is the fraction field of the PID R, given by
 * generators.  Its basis is the echelon form of the generating matrix.
 *
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_pid
 */
export class FreeModuleSubmodulePID extends FreeModuleWithBasis {
  constructor(
    ambient: FreeModuleGeneric,
    gens: FreeModuleElement[],
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ) {
    super(ambient, gens, {
      check: options?.check ?? true,
      echelonize: !(options?.alreadyEchelonized ?? false),
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    });
  }

  /**
   * Return False: the basis is the echelon form, not a user basis.
   */
  override hasUserBasis(): boolean {
    return false;
  }
}

/**
 * Subspace of a vector space with a user-specified basis.
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_field
 */
export class FreeModuleSubspaceWithBasis extends FreeModuleField {
  protected _ambient: FreeModuleField;
  protected _userBasis: FreeModuleElement[];
  protected _echelonizedBasisMatrix: unknown[][] | null = null;

  constructor(
    ambient: FreeModuleField,
    basis: FreeModuleElement[],
    options?: { check?: boolean; echelonize?: boolean; alreadyEchelonized?: boolean }
  ) {
    const { user, echelonized } = submoduleBases(ambient, basis, options);
    super(
      ambient.baseRing(),
      user.length,
      ambient.degree(),
      ambient.isSparse(),
      ambient.innerProductMatrix()
    );

    this._ambient = ambient;
    this._userBasis = user.map((row) => makeVector(this, row));
    this._basis = this._userBasis;
    this._echelonizedBasisMatrix = echelonized;
  }

  override hasUserBasis(): boolean {
    return true;
  }

  override ambientModule(): FreeModuleGeneric {
    return this._ambient;
  }

  override echelonizedBasisMatrix(): unknown[][] {
    if (this._echelonizedBasisMatrix === null) {
      this._echelonizedBasisMatrix = echelonRows(
        this.basisMatrix() as unknown[][],
        arithmeticFor(this._baseRing)
      );
    }
    return this._echelonizedBasisMatrix;
  }
}

/**
 * A subspace of a vector space, given by generators; its basis is the reduced
 * row echelon form of the generating matrix.
 *
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_field
 */
export class FreeModuleSubspace extends FreeModuleSubspaceWithBasis {
  constructor(
    ambient: FreeModuleField,
    gens: FreeModuleElement[],
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ) {
    super(ambient, gens, {
      check: options?.check ?? true,
      echelonize: !(options?.alreadyEchelonized ?? false),
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    });
  }

  /**
   * Return False: the basis is the echelon form, not a user basis.
   */
  override hasUserBasis(): boolean {
    return false;
  }
}

// ============================================================================
// Quotient Modules
// ============================================================================

/**
 * Compute SageMath's quotient and lift matrices for `V/W` over a field.
 *
 * Returns `[Q, L]`, where `Q` is the `n x m` matrix of the quotient map
 * (`n = dim V`, `m = dim V - dim W`) in terms of the basis of `V`, and `L` is
 * the `m x n` matrix of the chosen section, again in terms of the basis of
 * `V`.
 *
 * @see Reference: sage/modules/free_module.py:5366
 *   (`FreeModule_generic_field._FreeModule_generic_field__quotient_matrices`)
 */
function quotientMatrices(
  cover: FreeModuleGeneric,
  sub: FreeModuleGeneric,
  ar: FractionFieldArithmetic
): { Q: unknown[][]; L: unknown[][] } {
  // Step 1: Find bases for spaces
  const B = liftRows(sub.basisMatrix() as unknown[][], ar);
  const S = liftRows(cover.basisMatrix() as unknown[][], ar);

  const n = cover.rank();
  const m = n - sub.rank();

  // Step 2: Extend the basis B to a basis for self, by taking the pivot rows
  // of B stacked on S.
  const stacked = [...B, ...S];
  const d = cover.degree();
  const C: unknown[][] = [];
  for (let j = 0; j < d; j++) {
    C.push(stacked.map((row) => row[j]));
  }
  const { pivots: cPivots } = rrefLifted(C, ar);
  const A = cPivots.map((i) => stacked[i]!);

  // Step 3: D is the change of basis from S to A, i.e. D * A = S.
  const { pivots: P } = rrefLifted(A, ar);
  const AA = A.map((row) => P.map((j) => row[j]));
  const SS = S.map((row) => P.map((j) => row[j]));
  const D = matMulLifted(SS, inverseLifted(AA, ar), ar);

  // The quotient map takes the last m coordinates.
  const Q = D.map((row) => row.slice(n - m, n));

  // Step 4: the section map.
  const Dinv = inverseLifted(D, ar);
  const L = Dinv.slice(n - m, n);

  return { Q, L };
}

/**
 * Quotient of a free module by a submodule, `Q = V/W`.
 *
 * Over a field this is SageMath's `FreeModule_ambient_field_quotient`: an
 * ambient vector space of dimension `dim V - dim W` together with the quotient
 * map `V -> Q` and a fixed section `Q -> V`.
 *
 * Over `ZZ` it is SageMath's `FGP_Module`: the finitely generated module
 * `V/W` presented by the Smith normal form of the matrix expressing a basis of
 * `W` in terms of a basis of `V`.  Elements are coordinate vectors with
 * respect to the Smith form generators, reduced modulo the invariants.
 *
 * @see Reference: sage/modules/quotient_module.py:305
 *   (FreeModule_ambient_field_quotient)
 * @see Reference: sage/modules/fg_pid/fgp_module.py:268 (FGP_Module_class)
 * @see Deviation: SageMath has two distinct classes for the two cases; the
 *   port uses one class whose `invariants()` is empty in the field case.
 */
export class FreeModuleQuotient extends FreeModuleGeneric {
  protected _cover: FreeModuleGeneric;
  protected _submodule: FreeModuleGeneric;
  private readonly _ar: FractionFieldArithmetic;
  /** Field case: the n x m quotient matrix and the m x n lift matrix. */
  private readonly _quoMatrix: unknown[][] | null = null;
  private readonly _liftMatrix: unknown[][] | null = null;
  /** ZZ case: the Smith form data. */
  private readonly _invariantsAll: bigint[] = [];
  private readonly _nonOne: number[] = [];
  /** ZZ case: X, with `new coordinates = old coordinates * X`. */
  private readonly _toSmith: unknown[][] = [];
  /** ZZ case: the Smith form generators, in ambient coordinates. */
  private readonly _smithGensMatrix: unknown[][] = [];

  constructor(cover: FreeModuleGeneric, submodule: FreeModuleGeneric) {
    const ar = arithmeticFor(cover.baseRing());
    const n = cover.rank();
    const nW = submodule.rank();

    let rank: number;
    let degree: number;
    let quo: unknown[][] | null = null;
    let lift: unknown[][] | null = null;
    let invariantsAll: bigint[] = [];
    let nonOne: number[] = [];
    let toSmith: unknown[][] = [];
    let smithGens: unknown[][] = [];

    if (ar.isField) {
      const M = quotientMatrices(cover, submodule, ar);
      quo = M.Q;
      lift = M.L;
      rank = n - nW;
      degree = rank;
    } else if (ar.isIntegral) {
      const data = smithPresentation(cover, submodule, ar);
      invariantsAll = data.invariantsAll;
      nonOne = data.nonOne;
      toSmith = data.toSmith;
      smithGens = data.smithGens;
      degree = nonOne.length;
      rank = invariantsAll.filter((e) => e === 0n).length;
    } else {
      throw new NotImplementedError(
        'quotients of modules over rings other than fields or ZZ is not fully implemented'
      );
    }

    super(cover.baseRing(), rank, degree, cover.isSparse());

    this._cover = cover;
    this._submodule = submodule;
    this._ar = ar;
    this._quoMatrix = quo;
    this._liftMatrix = lift;
    this._invariantsAll = invariantsAll;
    this._nonOne = nonOne;
    this._toSmith = toSmith;
    this._smithGensMatrix = smithGens;
  }

  override isAmbient(): boolean {
    return this._quoMatrix !== null;
  }

  /**
   * Return the covering module `V`.
   * @see Reference: sage/modules/quotient_module.py:669 (cover)
   */
  coveringModule(): FreeModuleGeneric {
    return this._cover;
  }

  /** Alias of {@link coveringModule}, matching SageMath's `V()`. */
  V(): FreeModuleGeneric {
    return this._cover;
  }

  /**
   * Return the submodule `W` that we are quotienting by.
   * @see Reference: sage/modules/quotient_module.py:683 (relations)
   */
  relations(): FreeModuleGeneric {
    return this._submodule;
  }

  /** Alias of {@link relations}, matching SageMath's `W()`. */
  W(): FreeModuleGeneric {
    return this._submodule;
  }

  /**
   * Return the invariants of this quotient: the diagonal entries of the Smith
   * normal form of the relative matrix, padded with zeros, excluding ones.
   *
   * Over a field this is the empty list (the quotient is a vector space).
   *
   * @see Reference: sage/modules/fg_pid/fgp_module.py:952 (invariants)
   */
  invariants(includeOnes: boolean = false): bigint[] {
    if (includeOnes) {
      return [...this._invariantsAll];
    }
    return this._invariantsAll.filter((e) => e !== 1n);
  }

  /**
   * Return the generators of this quotient.
   *
   * Over a field these are the standard basis vectors of the quotient space;
   * over ZZ they are the Smith form generators, whose orders are the
   * invariants.
   *
   * @see Reference: sage/modules/fg_pid/fgp_module.py:1016 (smith_form_gens)
   */
  override basis(): FreeModuleElement[] {
    if (this._basis !== null) {
      return this._basis;
    }
    const zero = this._baseRing.zero();
    const one = this._baseRing.one();
    const out: FreeModuleElement[] = [];
    for (let i = 0; i < this._degree; i++) {
      const entries: unknown[] = new Array(this._degree).fill(zero);
      entries[i] = one;
      out.push(makeVector(this, entries));
    }
    this._basis = out;
    return out;
  }

  override ngens(): number {
    return this._degree;
  }

  /**
   * Return the lift of an element of the quotient to the covering module `V`.
   *
   * This is a fixed section of the quotient map, so `project(lift(x)) == x`
   * for every `x`.
   *
   * @param v - An element of the quotient
   * @returns A representative in the ambient module of `V`
   *
   * @see Reference: sage/modules/quotient_module.py:650 (lift)
   * @see Reference: sage/modules/fg_pid/fgp_element.py:lift
   */
  lift(v: FreeModuleElement): FreeModuleElement {
    const x = v.list();
    if (x.length !== this._degree) {
      throw new ArithmeticError(`vector must have degree ${this._degree}`);
    }
    const ar = this._ar;
    const ambient = this._cover.ambientModule();

    if (this._liftMatrix !== null) {
      // coordinates in V = x * L, then the corresponding element of V.
      const coeffs = matMulLifted([x.map((e) => ar.lift(e))], this._liftMatrix, ar)[0]!;
      return this.combineCover(coeffs);
    }

    // ZZ: sum over the Smith form generators.
    const full: unknown[] = new Array(this._invariantsAll.length).fill(ar.zero());
    for (let j = 0; j < this._nonOne.length; j++) {
      full[this._nonOne[j]!] = ar.lift(x[j]);
    }
    const row = matMulLifted([full], this._smithGensMatrix, ar)[0]!;
    const entries = row.map((e) => ar.lower(e));
    const w = ambient.createElement(entries);
    w.setImmutable();
    return w;
  }

  /** `coeffs * basis(V)`, as an element of the ambient module of `V`. */
  private combineCover(coeffs: unknown[]): FreeModuleElement {
    const ar = this._ar;
    const S = liftRows(this._cover.basisMatrix() as unknown[][], ar);
    const row = matMulLifted([coeffs], S, ar)[0] ?? new Array(this._cover.degree()).fill(ar.zero());
    const w = this._cover.ambientModule().createElement(row.map((e) => ar.lower(e)));
    w.setImmutable();
    return w;
  }

  /**
   * Return the image of `v` in the quotient.
   *
   * `v` must lie in the covering module `V`; the result is the coordinate
   * vector of the coset `v + W`, which is zero exactly when `v` lies in `W`.
   *
   * @param v - An element of the covering module
   *
   * @see Reference: sage/modules/quotient_module.py:604 (quotient_map)
   * @see Reference: sage/modules/fg_pid/fgp_module.py:1221 (coordinate_vector)
   */
  project(v: FreeModuleElement): FreeModuleElement {
    const ar = this._ar;
    // c = V.coordinate_vector(v); raises if v is not in V.
    const c = this._cover.coordinates(v).map((e) => ar.lift(e));

    if (this._quoMatrix !== null) {
      const row = matMulLifted([c], this._quoMatrix, ar)[0] ?? [];
      const w = makeVector(
        this,
        row.map((e) => ar.lower(e))
      );
      return w;
    }

    // ZZ: rewrite in the Smith basis and reduce modulo the invariants.
    const nc = matMulLifted([c], this._toSmith, ar)[0] ?? [];
    const entries: unknown[] = [];
    for (const i of this._nonOne) {
      const e = ar.lower(nc[i]) as bigint;
      const inv = this._invariantsAll[i]!;
      entries.push(inv === 0n ? e : ((e % inv) + inv) % inv);
    }
    return makeVector(this, entries);
  }

  /**
   * Return the number of elements of this quotient.
   *
   * @see Reference: sage/modules/fg_pid/fgp_module.py:1729 (cardinality)
   */
  override cardinality(): bigint | number {
    if (this._quoMatrix !== null) {
      return super.cardinality();
    }
    let c = 1n;
    for (const e of this.invariants()) {
      if (e === 0n) {
        return Number.POSITIVE_INFINITY;
      }
      c *= e;
    }
    return c;
  }

  override toString(): string {
    const ringName = this._baseRing.toString?.() ?? 'Ring';
    if (this._quoMatrix !== null) {
      const kind = this._sparse ? 'Sparse vector' : 'Vector';
      return (
        `${kind} space quotient V/W of dimension ${this._rank} over ${ringName} where\n` +
        `V: ${this._cover.toString()}\nW: ${this._submodule.toString()}`
      );
    }
    const inv = this.invariants();
    const body = inv.length === 1 ? `(${inv[0]})` : `(${inv.join(', ')})`;
    return `Finitely generated module V/W over ${ringName} with invariants ${body}`;
  }
}

/**
 * Compute the Smith normal form presentation of `V/W` over ZZ.
 *
 * @see Reference: sage/modules/fg_pid/fgp_module.py:890 (_relative_matrix)
 * @see Reference: sage/modules/fg_pid/fgp_module.py:915 (_smith_form)
 * @see Reference: sage/modules/fg_pid/fgp_module.py:1016 (smith_form_gens)
 */
function smithPresentation(
  cover: FreeModuleGeneric,
  sub: FreeModuleGeneric,
  ar: FractionFieldArithmetic
): {
  invariantsAll: bigint[];
  nonOne: number[];
  toSmith: unknown[][];
  smithGens: unknown[][];
} {
  const n = cover.rank();

  // _relative_matrix: each basis vector of W in terms of the basis of V.
  const A: bigint[][] = [];
  for (const b of sub.basis()) {
    const coords = cover.coordinates(b);
    A.push(
      coords.map((c) => {
        const r = ar.lift(c) as Rational;
        if (!r.isInteger()) {
          throw new ArithmeticError('sub must be a submodule of self');
        }
        return r.numerator;
      })
    );
  }

  if (n === 0) {
    return { invariantsAll: [], nonOne: [], toSmith: [], smithGens: [] };
  }

  let D: bigint[][];
  let X: bigint[][];
  if (A.length === 0) {
    // W = 0: the Smith form is empty and X is the identity.
    D = [];
    X = identityBigint(n);
  } else {
    const [Dm, , Xm] = smith_form_integer(IntegerMatrixFromEntries(A), true) as [
      { nrows: number; ncols: number; get: (i: number, j: number) => { value: bigint } },
      unknown,
      { nrows: number; ncols: number; get: (i: number, j: number) => { value: bigint } },
    ];
    D = matrixToBigint(Dm);
    X = matrixToBigint(Xm);
  }

  // invariants: the diagonal of D, padded with zeros to length n.
  const invariantsAll: bigint[] = [];
  for (let i = 0; i < D.length; i++) {
    const e = D[i]![i] ?? 0n;
    invariantsAll.push(e < 0n ? -e : e);
  }
  while (invariantsAll.length < n) {
    invariantsAll.push(0n);
  }

  const nonOne: number[] = [];
  for (let i = 0; i < n; i++) {
    if (invariantsAll[i] !== 1n) {
      nonOne.push(i);
    }
  }

  // Y = X^{-1}; Z = Y * basis_matrix(V) expresses the new basis in ambient
  // coordinates.
  const Xl = X.map((row) => row.map((e) => ar.lift(e)));
  const Y = inverseLifted(Xl, ar);
  const S = liftRows(cover.basisMatrix() as unknown[][], ar);
  const smithGensAll = matMulLifted(Y, S, ar);

  return {
    invariantsAll,
    nonOne,
    toSmith: Xl,
    smithGens: smithGensAll,
  };
}

/** The n x n identity matrix over ZZ. */
function identityBigint(n: number): bigint[][] {
  const I: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < n; j++) {
      row.push(i === j ? 1n : 0n);
    }
    I.push(row);
  }
  return I;
}

/** Convert an IntegerMatrix to a plain bigint matrix. */
function matrixToBigint(M: {
  nrows: number;
  ncols: number;
  get: (i: number, j: number) => { value: bigint };
}): bigint[][] {
  const out: bigint[][] = [];
  for (let i = 0; i < M.nrows; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < M.ncols; j++) {
      row.push(M.get(i, j).value);
    }
    out.push(row);
  }
  return out;
}
