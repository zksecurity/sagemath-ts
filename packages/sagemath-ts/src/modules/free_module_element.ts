/**
 * @module sage/modules/free_module_element
 * @description Free module elements (vectors)
 * @see Reference: sage/modules/free_module_element.py
 */

import { ArithmeticError, NotImplementedError, ValueError } from '../errors.js';
import { Rational } from '../rings/rational.js';

/**
 * Raised when a sequence index is out of range.
 *
 * Mirrors Python's builtin :exc:`IndexError`, which SageMath raises for
 * out-of-range vector indices (``free_module_element.pyx:1910``).
 */
export class IndexError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'IndexError';
  }
}

/**
 * Ring-like interface for the base ring of a free module.
 */
export interface RingLike {
  zero(): unknown;
  one(): unknown;
  __call__?(x: unknown): unknown;
  is_field?(): boolean;
  is_exact?(): boolean;
  toString?(): string;
}

/**
 * Parent module interface.
 */
export interface FreeModuleParent {
  degree(): number;
  baseRing(): RingLike;
  isSparse(): boolean;
  innerProductMatrix?(): unknown;
  isAmbient?(): boolean;
  usesAmbientInnerProduct?(): boolean;
  ambientModule?(): FreeModuleParent;
  coordinateVector?(v: FreeModuleElement, check?: boolean): FreeModuleElement;
}

/**
 * An element of a free module.
 * @see Reference: sage/modules/free_module_element.py:FreeModuleElement
 */
export class FreeModuleElement {
  protected _parent: FreeModuleParent;
  protected _entries: unknown[];
  protected _degree: number;
  protected _isMutable: boolean = true;

  /**
   * Create a free module element.
   * @param parent - The parent free module
   * @param entries - The entries of the vector
   */
  constructor(parent: FreeModuleParent, entries: unknown[]) {
    this._parent = parent;
    this._degree = parent.degree();
    this._isMutable = true;

    // Validate and copy entries
    if (entries.length !== this._degree) {
      throw new ValueError(`entries must have length ${this._degree}, got ${entries.length}`);
    }
    this._entries = [...entries];
  }

  /**
   * Return the parent module of this element.
   */
  parent(): FreeModuleParent {
    return this._parent;
  }

  /**
   * Return the i-th entry of this vector.
   *
   * Negative indices count from the end of the vector, exactly as in Python.
   *
   * @param i - The index
   * @throws {IndexError} If the (wrapped) index is out of range
   * @see Reference: sage/modules/free_module_element.pyx:__getitem__
   */
  getItem(i: number): unknown {
    let n = i;
    if (n < 0) {
      n += this._degree;
    }
    if (n < 0 || n >= this._degree) {
      throw new IndexError('vector index out of range');
    }
    return this._entries[n];
  }

  /**
   * Set the i-th entry of this vector.
   *
   * Negative indices count from the end of the vector, exactly as in Python.
   *
   * @param i - The index
   * @param value - The value to set
   * @throws {IndexError} If the (wrapped) index is out of range
   * @see Reference: sage/modules/free_module_element.pyx:__setitem__
   */
  setItem(i: number, value: unknown): void {
    if (!this._isMutable) {
      throw new ValueError('vector is immutable; please change a copy instead (use copy())');
    }
    let n = i;
    if (n < 0) {
      n += this._degree;
    }
    if (n < 0 || n >= this._degree) {
      throw new IndexError('vector index out of range');
    }
    this._entries[n] = value;
  }

  /**
   * Return the entries of this vector as a list.
   */
  list(): unknown[] {
    return [...this._entries];
  }

  /**
   * Return the degree (number of entries) of this vector.
   */
  degree(): number {
    return this._degree;
  }

  /**
   * Return the sum of this vector and other.
   * @param other - Another vector
   */
  add(other: FreeModuleElement): FreeModuleElement {
    if (this._degree !== other._degree) {
      throw new ArithmeticError(
        `vectors must have same degree: ${this._degree} vs ${other._degree}`
      );
    }

    const result: unknown[] = [];
    for (let i = 0; i < this._degree; i++) {
      result.push(addElements(this._entries[i], other._entries[i]));
    }
    return new FreeModuleElement(this._parent, result);
  }

  /**
   * Return the difference of this vector and other.
   * @param other - Another vector
   */
  sub(other: FreeModuleElement): FreeModuleElement {
    if (this._degree !== other._degree) {
      throw new ArithmeticError(
        `vectors must have same degree: ${this._degree} vs ${other._degree}`
      );
    }

    const result: unknown[] = [];
    for (let i = 0; i < this._degree; i++) {
      result.push(subElements(this._entries[i], other._entries[i]));
    }
    return new FreeModuleElement(this._parent, result);
  }

  /**
   * Return the negation of this vector.
   */
  neg(): FreeModuleElement {
    const result: unknown[] = [];
    for (let i = 0; i < this._degree; i++) {
      result.push(negElement(this._entries[i]));
    }
    return new FreeModuleElement(this._parent, result);
  }

  /**
   * Return the scalar multiple of this vector.
   * @param scalar - A scalar
   */
  mul(scalar: unknown): FreeModuleElement {
    const result: unknown[] = [];
    for (let i = 0; i < this._degree; i++) {
      result.push(mulElements(scalar, this._entries[i]));
    }
    return new FreeModuleElement(this._parent, result);
  }

  /**
   * Return the dot product of this vector with other.
   * @param other - Another vector
   */
  dotProduct(other: FreeModuleElement): unknown {
    if (this._degree !== other._degree) {
      throw new ArithmeticError(
        `vectors must have same degree: ${this._degree} vs ${other._degree}`
      );
    }

    if (this._degree === 0) {
      return this._parent.baseRing().zero();
    }

    let sum = mulElements(this._entries[0], other._entries[0]);
    for (let i = 1; i < this._degree; i++) {
      sum = addElements(sum, mulElements(this._entries[i], other._entries[i]));
    }
    return sum;
  }

  /**
   * Return the inner product of this vector with other.
   * Uses the inner product matrix of the parent module.
   * @param other - Another vector
   */
  innerProduct(other: FreeModuleElement): unknown {
    const M = this._parent;

    // (x)^t A y, where A is the inner product matrix of the ambient module.
    // If the module carries no inner product matrix (or the identity), the
    // inner product is the dot product.
    const isAmbient = M.isAmbient?.() ?? true;
    if (isAmbient) {
      const A = innerProductRows(M.innerProductMatrix?.(), this._degree);
      if (A === null) {
        return this.dotProduct(other);
      }
      return bilinearForm(A, this.list(), other.list(), M.baseRing());
    }

    // Submodules use the inner product induced by the ambient module unless
    // they were given one of their own.
    if (M.usesAmbientInnerProduct?.() ?? true) {
      const ambient = M.ambientModule?.() ?? M;
      const A = innerProductRows(ambient.innerProductMatrix?.(), this._degree);
      if (A === null) {
        return this.dotProduct(other);
      }
      return bilinearForm(A, this.list(), other.list(), M.baseRing());
    }

    const A = innerProductRows(M.innerProductMatrix?.(), this._degree);
    if (A === null) {
      return this.dotProduct(other);
    }
    if (typeof M.coordinateVector !== 'function') {
      throw new NotImplementedError('inner product requires coordinate_vector on the parent');
    }
    const v = M.coordinateVector(this).list();
    const w = M.coordinateVector(other).list();
    return bilinearForm(A, v, w, M.baseRing());
  }

  /**
   * Return the pairwise product of this vector with other.
   * @param other - Another vector
   */
  pairwiseProduct(other: FreeModuleElement): FreeModuleElement {
    if (this._degree !== other._degree) {
      throw new ArithmeticError(
        `vectors must have same degree: ${this._degree} vs ${other._degree}`
      );
    }

    const result: unknown[] = [];
    for (let i = 0; i < this._degree; i++) {
      result.push(mulElements(this._entries[i], other._entries[i]));
    }
    return new FreeModuleElement(this._parent, result);
  }

  /**
   * Return the cross product of this vector with other.
   *
   * Only defined for vectors of length 3 (via the quaternions) or 7 (via the
   * octonions).
   *
   * @param other - Another vector of the same (3 or 7) degree
   * @throws {TypeError} If the degrees are not both 3 or both 7
   * @see Reference: sage/modules/free_module_element.pyx:cross_product
   */
  crossProduct(other: FreeModuleElement): FreeModuleElement {
    const l = this._entries;
    const r = other._entries;

    if (l.length === 3 && r.length === 3) {
      // cross product: [l1*r2 - l2*r1, l2*r0 - l0*r2, l0*r1 - l1*r0]
      const result = [
        subElements(mulElements(l[1], r[2]), mulElements(l[2], r[1])),
        subElements(mulElements(l[2], r[0]), mulElements(l[0], r[2])),
        subElements(mulElements(l[0], r[1]), mulElements(l[1], r[0])),
      ];
      return new FreeModuleElement(this._parent, result);
    }

    if (l.length === 7 && r.length === 7) {
      // Seven dimensional cross product, via the octonions.
      const term = (i: number, j: number): unknown =>
        subElements(mulElements(l[i], r[j]), mulElements(l[j], r[i]));
      const three = (a: [number, number], b: [number, number], c: [number, number]): unknown =>
        addElements(addElements(term(a[0], a[1]), term(b[0], b[1])), term(c[0], c[1]));

      const result = [
        three([1, 3], [2, 6], [4, 5]),
        three([2, 4], [3, 0], [5, 6]),
        three([3, 5], [4, 1], [6, 0]),
        three([4, 6], [5, 2], [0, 1]),
        three([5, 0], [6, 3], [1, 2]),
        three([6, 1], [0, 4], [2, 3]),
        three([0, 2], [1, 5], [3, 4]),
      ];
      return new FreeModuleElement(this._parent, result);
    }

    throw new TypeError(
      `Cross product only defined for vectors of length three or seven, not (${l.length} and ${r.length})`
    );
  }

  /**
   * Return the p-norm of this vector.
   *
   * The norm is computed exactly whenever the p-th root of
   * `sum |x_i|^p` is rational: the result is then a bigint (when it is an
   * integer) or a {@link Rational}.  When the root is irrational SageMath
   * returns a symbolic expression; this port has no symbolic ring and
   * returns the double-precision value instead.
   *
   * @param p - The norm parameter (default: 2), `Infinity` for the max norm
   * @see Reference: sage/modules/free_module_element.pyx:norm
   * @see Deviation: irrational norms are returned as doubles
   */
  norm(p: number = 2): unknown {
    const exact = this.normExact(p);
    if (exact === null) {
      return this.normNumeric(p);
    }
    return lowerRational(exact, this._entries[0]);
  }

  /**
   * Return the p-norm as an exact rational, or `null` when it is irrational.
   */
  private normExact(p: number): Rational | null {
    const abs = this._entries.map((e) => absExact(e));

    if (p === Number.POSITIVE_INFINITY) {
      if (this._degree === 0) {
        throw new ValueError('max() arg is an empty sequence');
      }
      let best = abs[0]!;
      for (const a of abs) {
        if (a.gt(best)) {
          best = a;
        }
      }
      return best;
    }

    if (p < 1) {
      throw new ValueError(`${p} is not greater than or equal to 1`);
    }

    if (!Number.isInteger(p)) {
      return null;
    }

    // s = sum |x_i|^p, computed exactly
    let s = Rational.zero();
    for (const a of abs) {
      s = s.add(a.pow(BigInt(p)));
    }
    if (p === 1) {
      return s;
    }
    if (s.is_nth_power(BigInt(p))) {
      return s.nth_root(BigInt(p));
    }
    return null;
  }

  /**
   * Return the p-norm in double precision (used when it is irrational).
   */
  private normNumeric(p: number): number {
    let s = 0;
    for (const e of this._entries) {
      s += absExact(e).toNumber() ** p;
    }
    return s ** (1 / p);
  }

  /**
   * Return the p-norm of this vector.
   *
   * Alias for {@link norm}; SageMath spells this `norm(p)`.
   *
   * @param p - The norm parameter (default: 2)
   */
  pNorm(p: number = 2): unknown {
    return this.norm(p);
  }

  /**
   * Return this vector divided by its p-norm.
   *
   * Note that normalizing a vector changes its base ring, exactly as in
   * SageMath: the result is a vector of {@link Rational}s (or of doubles when
   * the norm is irrational).
   *
   * @param p - The norm parameter (default: 2)
   * @see Reference: sage/modules/free_module_element.pyx:normalized
   * @see Deviation: the entries are Rationals (or doubles when the norm is
   *   irrational), since this port has no symbolic ring.
   */
  normalized(p: number = 2): FreeModuleElement {
    const nr = this.normExact(p);

    if (nr === null) {
      // Irrational norm: fall back to double precision (see norm()).
      const n = this.normNumeric(p);
      if (n === 0) {
        throw new ArithmeticError('cannot normalize the zero vector');
      }
      const result = this._entries.map((e) => absExactSigned(e).toNumber() / n);
      return new FreeModuleElement(numericParent(this._degree), result);
    }

    if (nr.isZero()) {
      throw new ArithmeticError('cannot normalize the zero vector');
    }

    const result = this._entries.map((e) => absExactSigned(e).div(nr));
    return new FreeModuleElement(rationalParent(this._degree), result);
  }

  /**
   * Return whether this vector is zero.
   */
  isZero(): boolean {
    const zero = this._parent.baseRing().zero();
    for (let i = 0; i < this._degree; i++) {
      if (!elementsEqual(this._entries[i], zero)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Return whether this vector is mutable.
   */
  isMutable(): boolean {
    return this._isMutable;
  }

  /**
   * Make this vector immutable.
   */
  setImmutable(): void {
    this._isMutable = false;
  }

  /**
   * Return a copy of this vector.
   */
  copy(): FreeModuleElement {
    const result = new FreeModuleElement(this._parent, [...this._entries]);
    return result;
  }

  /**
   * Return the support of this vector (indices of non-zero entries).
   */
  support(): number[] {
    const zero = this._parent.baseRing().zero();
    const indices: number[] = [];
    for (let i = 0; i < this._degree; i++) {
      if (!elementsEqual(this._entries[i], zero)) {
        indices.push(i);
      }
    }
    return indices;
  }

  /**
   * Return the number of non-zero entries.
   */
  hammingWeight(): number {
    return this.support().length;
  }

  /**
   * Return a string representation.
   */
  toString(): string {
    const entries = this._entries.map((e) => String(e)).join(', ');
    return `(${entries})`;
  }

  /**
   * Return whether this vector equals other.
   * @param other - Another vector
   */
  equals(other: FreeModuleElement): boolean {
    if (this._degree !== other._degree) {
      return false;
    }
    for (let i = 0; i < this._degree; i++) {
      if (!elementsEqual(this._entries[i], other._entries[i])) {
        return false;
      }
    }
    return true;
  }
}

/**
 * A dense element of a free module.
 * @see Reference: sage/modules/free_module_element.py:FreeModuleElement_generic_dense
 */
export class FreeModuleElementDense extends FreeModuleElement {
  constructor(parent: FreeModuleParent, entries: unknown[]) {
    super(parent, entries);
  }
}

/**
 * A sparse element of a free module.
 * @see Reference: sage/modules/free_module_element.py:FreeModuleElement_generic_sparse
 */
export class FreeModuleElementSparse extends FreeModuleElement {
  private _sparseEntries: Map<number, unknown>;

  constructor(parent: FreeModuleParent, entries: Map<number, unknown> | unknown[]) {
    // Convert to dense for the parent constructor
    const degree = parent.degree();
    const zero = parent.baseRing().zero();
    const denseEntries: unknown[] = new Array(degree).fill(zero);

    if (entries instanceof Map) {
      for (const [idx, val] of entries) {
        if (idx >= 0 && idx < degree) {
          denseEntries[idx] = val;
        }
      }
    } else {
      for (let i = 0; i < Math.min(entries.length, degree); i++) {
        denseEntries[i] = entries[i];
      }
    }

    super(parent, denseEntries);

    // Store sparse representation
    this._sparseEntries = new Map();
    for (let i = 0; i < degree; i++) {
      if (!elementsEqual(denseEntries[i], zero)) {
        this._sparseEntries.set(i, denseEntries[i]);
      }
    }
  }

  /**
   * Return a dictionary mapping indices to non-zero values.
   */
  dict(): Map<number, unknown> {
    return new Map(this._sparseEntries);
  }

  /**
   * Override support to use sparse representation.
   */
  override support(): number[] {
    return Array.from(this._sparseEntries.keys()).sort((a, b) => a - b);
  }

  /**
   * Override hammingWeight to use sparse representation.
   */
  override hammingWeight(): number {
    return this._sparseEntries.size;
  }
}

// ============================================================================
// Helper functions for element arithmetic
// ============================================================================

/**
 * Convert a ring element to an exact rational, preserving its sign.
 */
function absExactSigned(e: unknown): Rational {
  if (typeof e === 'bigint') {
    return new Rational(e);
  }
  if (typeof e === 'number') {
    return Rational.from(e);
  }
  if (e instanceof Rational) {
    return e;
  }
  if (typeof e === 'object' && e !== null) {
    const value = (e as { value?: unknown }).value;
    if (typeof value === 'bigint') {
      return new Rational(value);
    }
    const num = (e as { numerator?: unknown }).numerator;
    const den = (e as { denominator?: unknown }).denominator;
    if (typeof num === 'bigint' && typeof den === 'bigint') {
      return new Rational(num, den);
    }
  }
  return Rational.from(String(e));
}

/**
 * Return |e| as an exact rational.
 */
function absExact(e: unknown): Rational {
  return absExactSigned(e).abs();
}

/**
 * Convert an exact rational back to the representation used by the entries.
 *
 * Integers are returned as bigints when the entries are bigints, doubles when
 * the entries are JavaScript numbers, and {@link Rational} otherwise.
 */
function lowerRational(r: Rational, sample: unknown): unknown {
  if (typeof sample === 'number') {
    return r.toNumber();
  }
  if (r.isInteger()) {
    return r.numerator;
  }
  return r;
}

/**
 * A parent for vectors with double-precision entries.
 */
function numericParent(degree: number): FreeModuleParent {
  return {
    degree: () => degree,
    baseRing: () => ({
      zero: () => 0,
      one: () => 1,
      is_field: () => true,
      is_exact: () => false,
      toString: () => 'Real Double Field',
    }),
    isSparse: () => false,
  };
}

/**
 * A parent for vectors with exact rational entries.
 */
function rationalParent(degree: number): FreeModuleParent {
  return {
    degree: () => degree,
    baseRing: () => ({
      zero: () => Rational.zero(),
      one: () => Rational.one(),
      is_field: () => true,
      is_exact: () => true,
      toString: () => 'Rational Field',
    }),
    isSparse: () => false,
  };
}

/**
 * Normalize an inner product matrix to a list of rows.
 *
 * Returns `null` when the module has no inner product matrix, or when the
 * inner product matrix is given as the scalar 1 (in which case the inner
 * product is the dot product, cf. `_inner_product_is_dot_product`).
 */
function innerProductRows(ipm: unknown, degree: number): unknown[][] | null {
  if (ipm === null || ipm === undefined) {
    return null;
  }
  // inner_product_matrix=1 means the identity matrix
  if (typeof ipm === 'number' || typeof ipm === 'bigint') {
    return null;
  }

  let rows: unknown[][];
  if (Array.isArray(ipm)) {
    rows = ipm as unknown[][];
  } else if (
    typeof ipm === 'object' &&
    typeof (ipm as { get?: unknown }).get === 'function' &&
    typeof (ipm as { nrows?: unknown }).nrows === 'number'
  ) {
    const M = ipm as { nrows: number; ncols: number; get: (i: number, j: number) => unknown };
    rows = [];
    for (let i = 0; i < M.nrows; i++) {
      const row: unknown[] = [];
      for (let j = 0; j < M.ncols; j++) {
        const entry = M.get(i, j);
        const value = (entry as { value?: unknown } | null)?.value;
        row.push(typeof value === 'bigint' ? value : entry);
      }
      rows.push(row);
    }
  } else {
    return null;
  }

  if (rows.length !== degree) {
    throw new ArithmeticError(`inner product matrix must have ${degree} rows, got ${rows.length}`);
  }
  return rows;
}

/**
 * Compute `x^t * A * y`.
 */
function bilinearForm(A: unknown[][], x: unknown[], y: unknown[], ring: RingLike): unknown {
  let total: unknown = null;
  for (let i = 0; i < x.length; i++) {
    for (let j = 0; j < y.length; j++) {
      const a = A[i]?.[j];
      if (a === undefined) {
        throw new ArithmeticError('inner product matrix has the wrong shape');
      }
      const term = mulElements(mulElements(x[i], a), y[j]);
      total = total === null ? term : addElements(total, term);
    }
  }
  return total === null ? ring.zero() : total;
}

/**
 * Add two ring elements.
 */
function addElements(a: unknown, b: unknown): unknown {
  if (typeof a === 'bigint' && typeof b === 'bigint') {
    return a + b;
  }
  if (typeof a === 'number' && typeof b === 'number') {
    return a + b;
  }
  if (
    typeof a === 'object' &&
    a !== null &&
    'add' in a &&
    typeof (a as { add: unknown }).add === 'function'
  ) {
    return (a as { add: (x: unknown) => unknown }).add(b);
  }
  throw new ArithmeticError(`cannot add ${typeof a} and ${typeof b}`);
}

/**
 * Subtract two ring elements.
 */
function subElements(a: unknown, b: unknown): unknown {
  if (typeof a === 'bigint' && typeof b === 'bigint') {
    return a - b;
  }
  if (typeof a === 'number' && typeof b === 'number') {
    return a - b;
  }
  if (
    typeof a === 'object' &&
    a !== null &&
    'sub' in a &&
    typeof (a as { sub: unknown }).sub === 'function'
  ) {
    return (a as { sub: (x: unknown) => unknown }).sub(b);
  }
  throw new ArithmeticError(`cannot subtract ${typeof a} and ${typeof b}`);
}

/**
 * Negate a ring element.
 */
function negElement(a: unknown): unknown {
  if (typeof a === 'bigint') {
    return -a;
  }
  if (typeof a === 'number') {
    return -a;
  }
  if (
    typeof a === 'object' &&
    a !== null &&
    'neg' in a &&
    typeof (a as { neg: unknown }).neg === 'function'
  ) {
    return (a as { neg: () => unknown }).neg();
  }
  throw new ArithmeticError(`cannot negate ${typeof a}`);
}

/**
 * Multiply two ring elements.
 */
function mulElements(a: unknown, b: unknown): unknown {
  if (typeof a === 'bigint' && typeof b === 'bigint') {
    return a * b;
  }
  if (typeof a === 'number' && typeof b === 'number') {
    return a * b;
  }
  if (typeof a === 'bigint' && typeof b === 'number') {
    return a * BigInt(b);
  }
  if (typeof a === 'number' && typeof b === 'bigint') {
    return BigInt(a) * b;
  }
  if (
    typeof a === 'object' &&
    a !== null &&
    'mul' in a &&
    typeof (a as { mul: unknown }).mul === 'function'
  ) {
    return (a as { mul: (x: unknown) => unknown }).mul(b);
  }
  if (
    typeof b === 'object' &&
    b !== null &&
    'mul' in b &&
    typeof (b as { mul: unknown }).mul === 'function'
  ) {
    return (b as { mul: (x: unknown) => unknown }).mul(a);
  }
  throw new ArithmeticError(`cannot multiply ${typeof a} and ${typeof b}`);
}

/**
 * Check if two ring elements are equal.
 */
function elementsEqual(a: unknown, b: unknown): boolean {
  if (a === b) {
    return true;
  }
  if (typeof a === 'bigint' && typeof b === 'bigint') {
    return a === b;
  }
  if (typeof a === 'number' && typeof b === 'number') {
    return a === b;
  }
  if (
    typeof a === 'object' &&
    a !== null &&
    'eq' in a &&
    typeof (a as { eq: unknown }).eq === 'function'
  ) {
    return (a as { eq: (x: unknown) => boolean }).eq(b);
  }
  if (
    typeof a === 'object' &&
    a !== null &&
    'equals' in a &&
    typeof (a as { equals: unknown }).equals === 'function'
  ) {
    return (a as { equals: (x: unknown) => boolean }).equals(b);
  }
  return false;
}

// ============================================================================
// Factory functions
// ============================================================================

/**
 * Create a vector from the given entries.
 * @param entries - The entries of the vector
 * @param ring - Optional base ring
 * @see Reference: sage/modules/free_module_element.py:vector
 */
export function vector(entries: unknown[], ring?: RingLike): FreeModuleElement {
  // Infer ring from entries if not provided
  const effectiveRing: RingLike = ring ?? {
    zero: () => {
      if (entries.length === 0) return 0n;
      const first = entries[0];
      if (typeof first === 'bigint') return 0n;
      if (typeof first === 'number') return 0;
      if (typeof first === 'object' && first !== null && 'parent' in first) {
        const parent = (first as { parent?: () => { zero?: () => unknown } }).parent?.();
        if (parent && 'zero' in parent && typeof parent.zero === 'function') {
          return parent.zero();
        }
      }
      return 0n;
    },
    one: () => {
      if (entries.length === 0) return 1n;
      const first = entries[0];
      if (typeof first === 'bigint') return 1n;
      if (typeof first === 'number') return 1;
      return 1n;
    },
  };

  const parent: FreeModuleParent = {
    degree: () => entries.length,
    baseRing: () => effectiveRing,
    isSparse: () => false,
  };

  return new FreeModuleElementDense(parent, entries);
}

/**
 * Create a zero vector of the given degree.
 * @param ring - The base ring
 * @param degree - The number of entries
 */
export function zeroVector(ring: RingLike, degree: number | bigint): FreeModuleElement {
  const degreeNum = typeof degree === 'bigint' ? Number(degree) : degree;
  if (degreeNum < 0) {
    throw new ValueError(`degree (=${degreeNum}) must be nonnegative`);
  }

  const zero = ring.zero();
  const entries = new Array(degreeNum).fill(zero);

  const parent: FreeModuleParent = {
    degree: () => degreeNum,
    baseRing: () => ring,
    isSparse: () => false,
  };

  return new FreeModuleElementDense(parent, entries);
}

/**
 * Create a random vector.
 * @param ring - The base ring
 * @param degree - The number of entries
 */
export function randomVector(ring: RingLike, degree: number | bigint): FreeModuleElement {
  const degreeNum = typeof degree === 'bigint' ? Number(degree) : degree;
  if (degreeNum < 0) {
    throw new ValueError(`degree (=${degreeNum}) must be nonnegative`);
  }

  // For random vectors, we need a ring with random_element method
  if (!ring || typeof (ring as { random_element?: unknown }).random_element !== 'function') {
    throw new NotImplementedError('random_vector requires a ring with random_element() method');
  }

  const randomElement = (ring as unknown as { random_element: () => unknown }).random_element;
  const entries: unknown[] = [];
  for (let i = 0; i < degreeNum; i++) {
    entries.push(randomElement.call(ring));
  }

  const parent: FreeModuleParent = {
    degree: () => degreeNum,
    baseRing: () => ring,
    isSparse: () => false,
  };

  return new FreeModuleElementDense(parent, entries);
}
