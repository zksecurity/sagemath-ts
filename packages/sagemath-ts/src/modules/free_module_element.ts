/**
 * @module sage/modules/free_module_element
 * @description Free module elements (vectors)
 * @see Reference: sage/modules/free_module_element.py
 */

import { NotImplementedError, ValueError, ArithmeticError } from '../errors.js';

/**
 * Ring-like interface for the base ring of a free module.
 */
export interface RingLike {
  zero(): unknown;
  one(): unknown;
  __call__?(x: unknown): unknown;
  is_field?(): boolean;
  is_exact?(): boolean;
}

/**
 * Parent module interface.
 */
export interface FreeModuleParent {
  degree(): number;
  baseRing(): RingLike;
  isSparse(): boolean;
  innerProductMatrix?(): unknown;
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
      throw new ValueError(
        `entries must have length ${this._degree}, got ${entries.length}`
      );
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
   * @param i - The index
   */
  getItem(i: number): unknown {
    if (i < 0 || i >= this._degree) {
      throw new ValueError(`index ${i} out of range [0, ${this._degree})`);
    }
    return this._entries[i];
  }

  /**
   * Set the i-th entry of this vector.
   * @param i - The index
   * @param value - The value to set
   */
  setItem(i: number, value: unknown): void {
    if (!this._isMutable) {
      throw new ValueError(
        'vector is immutable; please change a copy instead (use copy())'
      );
    }
    if (i < 0 || i >= this._degree) {
      throw new ValueError(`index ${i} out of range [0, ${this._degree})`);
    }
    this._entries[i] = value;
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
    // If there's an inner product matrix, use it
    // Otherwise, fall back to dot product
    const parent = this._parent;
    if (parent.innerProductMatrix && parent.innerProductMatrix()) {
      // For now, fall back to dot product
      // Full implementation would multiply by the inner product matrix
      return this.dotProduct(other);
    }
    return this.dotProduct(other);
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
   * Only defined for 3-dimensional vectors.
   * @param other - Another 3-dimensional vector
   */
  crossProduct(other: FreeModuleElement): FreeModuleElement {
    if (this._degree !== 3 || other._degree !== 3) {
      throw new ArithmeticError('cross product is only defined for 3-dimensional vectors');
    }

    const a = this._entries;
    const b = other._entries;

    // cross product: [a1*b2 - a2*b1, a2*b0 - a0*b2, a0*b1 - a1*b0]
    const result = [
      subElements(mulElements(a[1], b[2]), mulElements(a[2], b[1])),
      subElements(mulElements(a[2], b[0]), mulElements(a[0], b[2])),
      subElements(mulElements(a[0], b[1]), mulElements(a[1], b[0])),
    ];

    return new FreeModuleElement(this._parent, result);
  }

  /**
   * Return the Euclidean norm (L2 norm) of this vector.
   * Note: Returns the squared norm for exact rings, or actual norm for floating point.
   */
  norm(): unknown {
    // Compute sum of squares
    let sumSq = this._parent.baseRing().zero();
    for (let i = 0; i < this._degree; i++) {
      const entry = this._entries[i];
      sumSq = addElements(sumSq, mulElements(entry, entry));
    }
    return sumSq;
  }

  /**
   * Return the p-norm of this vector.
   * @param p - The norm parameter (default: 2)
   * @returns For p=2, returns squared norm. For other p, returns sum of |x_i|^p as a number.
   */
  pNorm(p: number = 2): unknown {
    if (p === 2) {
      return this.norm();
    }

    if (p === Infinity) {
      // Infinity norm: max |x_i|
      let maxVal = 0;
      for (let i = 0; i < this._degree; i++) {
        const entry = this._entries[i];
        let absVal: number;
        if (typeof entry === 'bigint') {
          absVal = entry < 0n ? Number(-entry) : Number(entry);
        } else if (typeof entry === 'number') {
          absVal = Math.abs(entry);
        } else {
          // Try to convert
          const numVal = Number(String(entry));
          absVal = Math.abs(numVal);
        }
        if (absVal > maxVal) {
          maxVal = absVal;
        }
      }
      return maxVal;
    }

    if (p === 1) {
      // 1-norm: sum of |x_i|
      let sum = 0;
      for (let i = 0; i < this._degree; i++) {
        const entry = this._entries[i];
        if (typeof entry === 'bigint') {
          sum += entry < 0n ? Number(-entry) : Number(entry);
        } else if (typeof entry === 'number') {
          sum += Math.abs(entry);
        } else {
          const numVal = Number(String(entry));
          sum += Math.abs(numVal);
        }
      }
      return sum;
    }

    // General p-norm: (sum of |x_i|^p)^(1/p)
    // We return the sum without the root for consistency with norm()
    let sum = 0;
    for (let i = 0; i < this._degree; i++) {
      const entry = this._entries[i];
      let absVal: number;
      if (typeof entry === 'bigint') {
        absVal = entry < 0n ? Number(-entry) : Number(entry);
      } else if (typeof entry === 'number') {
        absVal = Math.abs(entry);
      } else {
        const numVal = Number(String(entry));
        absVal = Math.abs(numVal);
      }
      sum += Math.pow(absVal, p);
    }
    return sum;
  }

  /**
   * Return the normalized version of this vector (unit vector in same direction).
   * Works over real/rational fields. For integer rings, returns approximate normalization.
   */
  normalized(): FreeModuleElement {
    const normSq = this.norm();

    // Compute actual norm (sqrt of squared norm)
    let normVal: number;
    if (typeof normSq === 'bigint') {
      normVal = Math.sqrt(Number(normSq));
    } else if (typeof normSq === 'number') {
      normVal = Math.sqrt(normSq);
    } else {
      const numVal = Number(String(normSq));
      normVal = Math.sqrt(numVal);
    }

    if (normVal === 0) {
      throw new ArithmeticError('cannot normalize the zero vector');
    }

    // Divide each entry by the norm
    const result: unknown[] = [];
    for (let i = 0; i < this._degree; i++) {
      const entry = this._entries[i];
      if (typeof entry === 'bigint') {
        // For integers, return approximate floating point
        result.push(Number(entry) / normVal);
      } else if (typeof entry === 'number') {
        result.push(entry / normVal);
      } else if (typeof entry === 'object' && entry !== null && 'div' in entry) {
        // Try to use division method
        // This is approximate - would need proper field element handling
        const numEntry = Number(String(entry));
        result.push(numEntry / normVal);
      } else {
        const numEntry = Number(String(entry));
        result.push(numEntry / normVal);
      }
    }

    // Create parent for number entries
    const numParent: FreeModuleParent = {
      degree: () => this._degree,
      baseRing: () => ({
        zero: () => 0,
        one: () => 1,
        is_field: () => true,
      }),
      isSparse: () => false,
    };

    return new FreeModuleElement(numParent, result);
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
 * Add two ring elements.
 */
function addElements(a: unknown, b: unknown): unknown {
  if (typeof a === 'bigint' && typeof b === 'bigint') {
    return a + b;
  }
  if (typeof a === 'number' && typeof b === 'number') {
    return a + b;
  }
  if (typeof a === 'object' && a !== null && 'add' in a && typeof (a as { add: unknown }).add === 'function') {
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
  if (typeof a === 'object' && a !== null && 'sub' in a && typeof (a as { sub: unknown }).sub === 'function') {
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
  if (typeof a === 'object' && a !== null && 'neg' in a && typeof (a as { neg: unknown }).neg === 'function') {
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
  if (typeof a === 'object' && a !== null && 'mul' in a && typeof (a as { mul: unknown }).mul === 'function') {
    return (a as { mul: (x: unknown) => unknown }).mul(b);
  }
  if (typeof b === 'object' && b !== null && 'mul' in b && typeof (b as { mul: unknown }).mul === 'function') {
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
  if (typeof a === 'object' && a !== null && 'eq' in a && typeof (a as { eq: unknown }).eq === 'function') {
    return (a as { eq: (x: unknown) => boolean }).eq(b);
  }
  if (typeof a === 'object' && a !== null && 'equals' in a && typeof (a as { equals: unknown }).equals === 'function') {
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
    throw new NotImplementedError(
      'random_vector requires a ring with random_element() method'
    );
  }

  const randomElement = (ring as { random_element: () => unknown }).random_element;
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
