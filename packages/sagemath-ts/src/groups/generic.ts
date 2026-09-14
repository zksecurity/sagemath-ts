/**
 * @module sage/groups/generic
 * @description Miscellaneous generic group functions
 *
 * A collection of functions implementing generic algorithms in arbitrary
 * groups, including additive and multiplicative groups.
 *
 * Port of: sage/groups/generic.py
 * Reference: reference/sage/src/sage/groups/generic.py
 */

import {
  CRT_list,
  type Factorization,
  factor,
  is_prime,
  isqrt,
  valuation,
} from '../arith/misc.js';
import { ArithmeticError, ValueError, ZeroDivisionError } from '../errors.js';
import { current_randstate } from '../misc/randstate.js';
import { IntegerModRing } from '../rings/finite_rings/integer_mod_ring.js';
import type { Integer } from '../rings/integer_ring.js';
import type { Rational } from '../rings/rational.js';
import { RealDoubleElement } from '../rings/real_double.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';

/**
 * Group operation names accepted as the 'operation' parameter.
 */
export const MULTIPLICATION_NAMES = ['multiplication', 'times', 'product', '*'] as const;
export const ADDITION_NAMES = ['addition', 'plus', 'sum', '+'] as const;

export type OperationType =
  | (typeof MULTIPLICATION_NAMES)[number]
  | (typeof ADDITION_NAMES)[number]
  | 'other';

const STANDARD_OP_ERROR =
  "in order to specify custom identity/inverse/op, operation must be 'other'";

function isMultiplicative(operation: OperationType): boolean {
  return MULTIPLICATION_NAMES.includes(operation as (typeof MULTIPLICATION_NAMES)[number]);
}

function isAdditive(operation: OperationType): boolean {
  return ADDITION_NAMES.includes(operation as (typeof ADDITION_NAMES)[number]);
}

function validateGroupOps<T>(
  operation: OperationType,
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): void {
  if (
    (isMultiplicative(operation) || isAdditive(operation)) &&
    (identity !== undefined || inverse !== undefined || op !== undefined)
  ) {
    throw new ValueError(STANDARD_OP_ERROR);
  }
  if (
    !isMultiplicative(operation) &&
    !isAdditive(operation) &&
    (identity === undefined || inverse === undefined || op === undefined)
  ) {
    throw new ValueError(
      'identity, inverse and operation must all be specified when operation is neither addition nor multiplication'
    );
  }
}

/**
 * Interface for group elements that support generic group operations.
 */
export interface GroupElement {
  eq(other: this): boolean;
}

/**
 * Interface for multiplicative group elements.
 */
export interface MultiplicativeGroupElement extends GroupElement {
  mul(other: MultiplicativeGroupElement): MultiplicativeGroupElement;
  inv(): MultiplicativeGroupElement;
  pow(n: bigint): MultiplicativeGroupElement;
  isOne(): boolean;
}

/**
 * Interface for additive group elements.
 */
export interface AdditiveGroupElement extends GroupElement {
  add(other: AdditiveGroupElement): AdditiveGroupElement;
  neg(): AdditiveGroupElement;
  mul(n: bigint): AdditiveGroupElement; // scalar multiplication
  isZero(): boolean;
}

/**
 * Generic group operations configuration.
 */
export interface GroupOps<T> {
  identity: T;
  inverse: (x: T) => T;
  op: (x: T, y: T) => T;
  power: (x: T, n: bigint) => T;
  isIdentity: (x: T) => boolean;
}

function elementsEqual<T>(a: T, b: T): boolean {
  if (a !== null && typeof a === 'object' && 'eq' in (a as object)) {
    const eqFn = (a as unknown as { eq?: (other: T) => boolean }).eq;
    if (typeof eqFn === 'function') {
      return eqFn.call(a, b);
    }
  }
  return a === b;
}

type GroupParent = { one?: () => unknown; zero?: () => unknown };

/** Sage parent methods map to parent() or a stored parent in the port. */
function getGroupParent(sample: unknown): GroupParent | undefined {
  const parent = (sample as { parent?: GroupParent | (() => GroupParent) }).parent;
  return typeof parent === 'function' ? parent.call(sample) : parent;
}

/** Resolve a parent identity or an existing host element's neutral power/action. */
function standardGroupIdentity<T>(sample: T, multiplicative: boolean): T | undefined {
  const parent = getGroupParent(sample);
  const identity = multiplicative ? parent?.one?.() : parent?.zero?.();
  if (identity !== undefined) return identity as T;
  if (typeof sample === 'bigint') return (multiplicative ? 1n : 0n) as T;
  if (typeof sample === 'number') return (multiplicative ? 1 : 0) as T;
  const element = sample as { pow?: (n: bigint) => T; mul?: (n: bigint) => T };
  return multiplicative ? element.pow?.(0n) : element.mul?.(0n);
}

/** Sage's operator.inv is exposed as __invert__ or the existing inv adapter. */
function invertGroupElement<T>(value: T): T {
  const element = value as { __invert__?: () => T; inv?: () => T };
  const inverse = element.__invert__ ?? element.inv;
  if (typeof inverse === 'function') return inverse.call(value);
  if (typeof value === 'number') return (1 / value) as T;
  throw new ValueError('Cannot invert element');
}

/** Accept Sage predicate names as well as existing camel-case host adapters. */
function isStandardIdentity<T>(value: T, multiplicative: boolean): boolean {
  const element = value as {
    is_one?: () => boolean;
    isOne?: () => boolean;
    is_zero?: () => boolean;
    isZero?: () => boolean;
  };
  const predicate = multiplicative
    ? (element.is_one ?? element.isOne)
    : (element.is_zero ?? element.isZero);
  if (typeof predicate === 'function') return predicate.call(value);
  const identity = standardGroupIdentity(value, multiplicative);
  return identity !== undefined && elementsEqual(value, identity);
}

function deriveStandardGroupOps<T extends GroupElement>(
  sample: T,
  operation: OperationType
): GroupOps<T> {
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);
  if (!isMult && !isAdd) {
    throw new ValueError(
      'identity, inverse and operation must all be specified when operation is neither addition nor multiplication'
    );
  }

  const identity = standardGroupIdentity(sample, isMult);

  if (identity === undefined) {
    throw new ValueError('identity could not be determined for standard operation');
  }

  const op = isMult
    ? (x: T, y: T) => {
        if ((x as unknown as MultiplicativeGroupElement).mul) {
          return (x as unknown as MultiplicativeGroupElement).mul(
            y as unknown as MultiplicativeGroupElement
          ) as unknown as T;
        }
        if (typeof x === 'bigint') {
          return ((x as unknown as bigint) * (y as unknown as bigint)) as unknown as T;
        }
        if (typeof x === 'number') {
          return ((x as unknown as number) * (y as unknown as number)) as unknown as T;
        }
        throw new ValueError('Cannot multiply elements');
      }
    : (x: T, y: T) => {
        if ((x as unknown as AdditiveGroupElement).add) {
          return (x as unknown as AdditiveGroupElement).add(
            y as unknown as AdditiveGroupElement
          ) as unknown as T;
        }
        if (typeof x === 'bigint') {
          return ((x as unknown as bigint) + (y as unknown as bigint)) as unknown as T;
        }
        if (typeof x === 'number') {
          return ((x as unknown as number) + (y as unknown as number)) as unknown as T;
        }
        throw new ValueError('Cannot add elements');
      };

  const inverse = isMult
    ? (x: T) => invertGroupElement(x)
    : (x: T) => {
        if ((x as unknown as AdditiveGroupElement).neg) {
          return (x as unknown as AdditiveGroupElement).neg() as unknown as T;
        }
        if (typeof x === 'bigint') {
          return -(x as unknown as bigint) as unknown as T;
        }
        if (typeof x === 'number') {
          return -(x as unknown as number) as unknown as T;
        }
        throw new ValueError('Cannot negate element');
      };

  const isIdentity = isMult
    ? (x: T) => {
        if ((x as unknown as MultiplicativeGroupElement).isOne) {
          return isStandardIdentity(x, true);
        }
        return elementsEqual(x, identity as T);
      }
    : (x: T) => {
        if ((x as unknown as AdditiveGroupElement).isZero) {
          return isStandardIdentity(x, false);
        }
        return elementsEqual(x, identity as T);
      };

  return {
    identity,
    inverse,
    op,
    power: (x: T, n: bigint) => multiple(x, n, 'other', identity, inverse, op),
    isIdentity,
  };
}

/**
 * Parse group operation specification.
 *
 * @param operation - Operation type or 'other' for custom operations
 * @param identity - Identity element (required for 'other')
 * @param inverse - Inverse function (required for 'other')
 * @param op - Binary operation function (required for 'other')
 * @param sample - Element whose parent supplies a standard identity
 * @returns Normalized group operations
 * @see Deviation: Generic group iterator and parent adapters
 */
export function parseGroupOps(
  operation: Exclude<OperationType, 'other'>,
  identity: undefined,
  inverse: undefined,
  op: undefined,
  sample: Integer
): GroupOps<Integer | Rational>;
export function parseGroupOps<T extends GroupElement>(
  operation: OperationType,
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T,
  sample?: T
): GroupOps<T>;
export function parseGroupOps<T extends GroupElement>(
  operation: OperationType,
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T,
  sample?: T
): GroupOps<T> {
  if (isMultiplicative(operation) || isAdditive(operation)) {
    validateGroupOps(operation, identity, inverse, op);
    if (sample === undefined) {
      throw new ValueError('identity could not be determined for standard operation');
    }
    return deriveStandardGroupOps(sample, operation);
  }

  // Custom operation
  if (identity === undefined || inverse === undefined || op === undefined) {
    throw new ValueError(
      "identity, inverse and operation must all be specified when operation is neither addition nor multiplication"
    );
  }

  return {
    identity,
    inverse,
    op,
    power: (x: T, n: bigint) => multiple(x, n, operation, identity, inverse, op),
    isIdentity: (x: T) => x.eq(identity),
  };
}

/** _power_func keeps native powers for normalized standard operations. */
function groupPowerFunction<T extends GroupElement>(
  operation: OperationType,
  ops: GroupOps<T>
): (x: T, n: bigint) => T {
  if (isMultiplicative(operation)) {
    return (x, n) => (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
  }
  if (isAdditive(operation)) {
    return (x, n) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
  }
  return ops.power;
}

/**
 * Compute n*a or a^n using binary algorithm.
 *
 * @param a - Group element
 * @param n - Integer (can be negative)
 * @param operation - Type of group operation
 * @param identity - Identity element (for custom operations)
 * @param inverse - Inverse function (for custom operations)
 * @param op - Binary operation (for custom operations)
 * @returns n*a (additive) or a^n (multiplicative)
 *
 * @example
 * ```typescript
 * // Multiplicative group
 * const x = Mod(3n, 7n);
 * const result = multiple(x, 5n, '*'); // 3^5 mod 7 = 5
 *
 * // Additive group (elliptic curves)
 * const P = E.point(...);
 * const result = multiple(P, 5n, '+'); // 5*P
 * ```
 */
export function multiple(
  a: Integer,
  n: IntegerLike,
  operation?: (typeof MULTIPLICATION_NAMES)[number],
  identity?: undefined,
  inverse?: undefined,
  op?: undefined
): Integer | Rational;
export function multiple<T extends GroupElement>(
  a: T,
  n: IntegerLike,
  operation?: OperationType,
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): T;
export function multiple<T extends GroupElement>(
  a: T,
  n: IntegerLike,
  operation: OperationType = '*',
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): T {
  let nBig = toBigInt(n);
  validateGroupOps(operation, identity, inverse, op);
  // generic.py:multiple always uses its binary operation schedule.
  const ops = parseGroupOps(operation, identity, inverse, op, a);

  if (nBig === 0n) {
    return ops.identity;
  }

  if (nBig < 0n) {
    nBig = -nBig;
    a = ops.inverse(a);
  }

  if (nBig === 1n) {
    return a;
  }

  // generic.py:342–379: idempotence and the native binary schedule.
  const aa = ops.op(a, a);
  if (elementsEqual(aa, a)) return a;
  if (nBig === 2n) return aa;
  if (nBig === 3n) return ops.op(aa, a);
  if (nBig === 4n) return ops.op(aa, aa);
  const m = nBig & 1n;
  nBig >>= 1n;
  let apow = aa;
  while ((nBig & 1n) === 0n) {
    apow = ops.op(apow, apow);
    nBig >>= 1n;
  }
  let result = apow;
  nBig >>= 1n;
  if (m) result = ops.op(result, a);
  while (nBig !== 0n) {
    apow = ops.op(apow, apow);
    if (nBig & 1n) result = ops.op(result, apow);
    nBig >>= 1n;
  }
  return result;
}

/**
 * Baby-step giant-step algorithm for discrete logarithm.
 *
 * Solves n*a = b (additive) or a^n = b (multiplicative) with lb <= n <= ub.
 *
 * @param a - Base element
 * @param b - Target element
 * @param bounds - Tuple [lb, ub] with 0 <= lb <= ub
 * @param operation - Type of group operation
 * @param identity - Identity element (for custom operations)
 * @param inverse - Inverse function (for custom operations)
 * @param op - Binary operation (for custom operations)
 * @returns n such that a^n = b (or n*a = b)
 * @throws {ValueError} If no such n exists
 *
 * @example
 * ```typescript
 * const b = Mod(2n, 37n);
 * const a = b.pow(20n);
 * const x = bsgs(b, a, [0n, 36n], '*'); // returns 20n
 * ```
 *
 * @see Deviation: Generic Group API and DLP
 */
export function bsgs<T extends GroupElement>(
  a: T,
  b: T,
  bounds: [IntegerLike, IntegerLike],
  operation: OperationType = '*',
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): bigint {
  const ops = parseGroupOps(operation, identity, inverse, op, a);
  const { power, op: multiply, inverse: invert, identity: identityElem } = ops;
  const isId = (x: T) => elementsEqual(identityElem, x);
  const [lbInput, ubInput] = bounds;
  const lb = toBigInt(lbInput);
  const ub = toBigInt(ubInput);

  if (lb < 0n || ub < lb) {
    throw new ValueError('bsgs() requires 0<=lb<=ub');
  }

  const range = 1n + ub - lb;

  // Handle identity base specially
  if (elementsEqual(a, identityElem) && !elementsEqual(b, identityElem)) {
    throw new ValueError('no solution in bsgs()');
  }

  // Compute b^(-1) * a^lb
  const bInv = invert(b);
  const aLb = power(a, lb);
  const c = multiply(bInv, aLb);

  if (range < 30n) {
    // Sage uses a simple linear search for small ranges (generic.py:625-633).
    // Without it the giant-step loop can return an n outside [lb, ub].
    let d = c;
    for (let i0 = 0n; i0 < range; i0++) {
      const i = lb + i0;
      if (isId(d)) {
        // identity == b^(-1)*a^i, so return i
        return i;
      }
      d = multiply(a, d);
    }
    throw new ValueError('no solution in bsgs()');
  }

  // Baby-step giant-step
  const m = isqrt(range) + 1n;

  // Baby steps: compute a^i for i = 0, 1, ..., m-1
  // Store in hash table: value -> index
  const table = new Map<string, Array<{ element: T; index: bigint }>>();
  let d = c;

  for (let i = 0n; i < m; i++) {
    if (isId(d)) {
      return lb + i;
    }
    // Use string representation as hash key
    const key = elementToString(d);
    const bucket = table.get(key) ?? [];
    // Python dict resolves hash collisions using identity, then equality.
    const entry = bucket.find(({ element }) => element === d || elementsEqual(element, d));
    if (entry) entry.index = lb + i;
    else bucket.push({ element: d, index: lb + i });
    table.set(key, bucket);
    d = multiply(d, a);
  }

  // Sage reuses the final baby-step state to form a^(-m).
  const giantFactor = multiply(c, invert(d));
  let giant = identityElem;

  for (let k = 0n; k < m; k++) {
    const key = elementToString(giant);
    const entry = table
      .get(key)
      ?.find(({ element }) => element === giant || elementsEqual(element, giant));
    if (entry !== undefined) {
      return k * m + entry.index;
    }
    giant = multiply(giantFactor, giant);
  }

  throw new ValueError(`log of ${b} to the base ${a} does not exist in (${lb}, ${ub})`);
}

/** Empty Python format specification in generic.py's f-string errors. */
function formatGroupElement(value: unknown): string {
  if (value instanceof RealDoubleElement) {
    // RDF.__format__ delegates to float, unlike its Sage _repr_.
    if (Number.isNaN(value.value)) return 'nan';
    if (value.value === Infinity) return 'inf';
    if (value.value === -Infinity) return '-inf';
  }
  return String(value);
}

/**
 * Convert a group element to a string for hashing.
 */
function elementToString<T>(elem: T): string {
  if (elem === null || elem === undefined) {
    return String(elem);
  }
  // Handle IntegerMod and similar types
  if (typeof (elem as unknown as { value: bigint }).value === 'bigint') {
    return String((elem as unknown as { value: bigint }).value);
  }
  // Handle objects with toString
  if (typeof (elem as unknown as { toString(): string }).toString === 'function') {
    return (elem as unknown as { toString(): string }).toString();
  }
  return String(elem);
}

/**
 * Core of Sage's `discrete_log`: Pohlig-Hellman with per-digit BSGS.
 *
 * This is a direct transcription of the loop in `sage/groups/generic.py:1077-1119`
 * (the `bounds is None` case). In particular it replicates Sage's repair of an
 * `ord` that is a proper multiple of the order of `base`, and reduces the CRT
 * result modulo the repaired `ord`.
 *
 * @internal
 */
function _discrete_log_core<T extends GroupElement>(
  a: T,
  base: T,
  ordIn: bigint,
  factorization: Factorization,
  ops: {
    power: (x: T, n: bigint) => T;
    multiply: (x: T, y: T) => T;
  },
  operation: OperationType,
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): bigint {
  const { power, multiply } = ops;
  let ord = ordIn;
  // Drop the unit factor (-1) if the factorization carries one.
  const f = factorization.filter(([p]) => p > 0n);

  const l: bigint[] = new Array(f.length).fill(0n);
  const mods: bigint[] = [];
  let runningMod = 1n;
  let i = -1;

  for (let idx = 0; idx < f.length; idx++) {
    i = idx;
    const pi = f[idx]![0];
    let ri = f[idx]![1];
    let gamma = power(base, ord / pi);
    // Pohlig-Hellman does not work with an incorrect order, and the caller
    // might have provided a proper multiple of the order of base.
    while (elementsEqual(gamma, power(gamma, 0n)) && ri > 0n) {
      ord = ord / pi;
      ri -= 1n;
      gamma = power(base, ord / pi);
    }
    const bound = ord - 1n;
    let runningBound = bound < pi ** ri - 1n ? bound : pi ** ri - 1n;
    let j = -1n;
    for (let jj = 0n; jj < ri; jj++) {
      j = jj;
      const tempBound = runningBound < pi - 1n ? runningBound : pi - 1n;
      const h = power(multiply(a, power(base, -l[idx]!)), ord / pi ** (jj + 1n));
      const c = bsgs(gamma, h, [0n, tempBound], operation, identity, inverse, op);
      l[idx] = l[idx]! + c * pi ** jj;
      runningBound = runningBound / pi;
      runningMod = runningMod * pi;
      if (runningMod > bound) break;
    }
    mods.push(pi ** (j + 1n));
    // we have log % runningMod; if log < runningMod we have the value of log
    if (runningMod > bound) break;
  }

  const residues = l.slice(0, i + 1);
  if (residues.length === 0) {
    return 0n;
  }
  const crt = CRT_list(residues, mods);
  return ord === 0n ? crt : ((crt % ord) + ord) % ord;
}

/**
 * Pohlig-Hellman algorithm for discrete logarithm.
 *
 * Reduces the discrete log problem to subproblems in prime power order subgroups,
 * then combines results using CRT.
 *
 * @param a - Target element (find x such that base^x = a)
 * @param base - Base element
 * @param n - Order of base (or a multiple of the order)
 * @param factorization - Factorization of n (computed if not provided)
 * @param operation - Type of group operation
 * @param identity - Identity element (for custom operations)
 * @param inverse - Inverse function (for custom operations)
 * @param op - Binary operation (for custom operations)
 * @returns x such that base^x = a
 *
 * @example
 * ```typescript
 * const base = Mod(2n, 37n);
 * const a = base.pow(20n);
 * const x = pohlig_hellman(a, base, 36n); // returns 20n
 * ```
 */
export function pohlig_hellman<T extends GroupElement>(
  a: T,
  base: T,
  n: IntegerLike,
  factorization?: Factorization,
  operation: OperationType = '*',
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): bigint {
  const nBig = toBigInt(n);
  const parsed = parseGroupOps(operation, identity, inverse, op, a);
  const power = groupPowerFunction(operation, parsed);
  const multiply = parsed.op;

  // Get factorization of n
  const factors = factorization ?? factor(nBig);

  return _discrete_log_core(
    a,
    base,
    nBig,
    factors,
    { power, multiply },
    'other',
    parsed.identity,
    parsed.inverse,
    parsed.op
  );
}

/**
 * Generic discrete logarithm function.
 *
 * Find x such that base^x = a (multiplicative) or x*base = a (additive).
 *
 * Uses Pohlig-Hellman combined with baby-step giant-step for efficiency.
 *
 * @param a - Target element
 * @param base - Base element
 * @param ord - A multiple of the order of base (computed from the element if
 *   not provided). Passing e.g. the cardinality of an elliptic curve is the
 *   documented usage; the Pohlig-Hellman loop repairs an over-large `ord`.
 * @param operation - Type of group operation
 * @param identity - Identity element (for custom operations)
 * @param inverse - Inverse function (for custom operations)
 * @param op - Binary operation (for custom operations)
 * @returns x such that base^x = a
 *
 * @example
 * ```typescript
 * // In (Z/37Z)*
 * const b = Mod(2n, 37n);
 * const a = b.pow(20n);
 * const x = discrete_log(a, b, 36n, '*'); // returns 20n
 *
 * // In an elliptic curve group
 * const P = E.generator();
 * const Q = P.mul(123n);
 * const x = discrete_log(Q, P, E.order(), '+'); // returns 123n
 * ```
 *
 * @see Reference: sage/groups/generic.py:discrete_log (line 823)
 * @see Deviation: Generic Group API and DLP
 */
export function discrete_log<T extends GroupElement>(
  a: T,
  base: T,
  ord?: IntegerLike,
  operation: OperationType = '*',
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): bigint {
  const parsed = parseGroupOps(operation, identity, inverse, op, a);
  const power = groupPowerFunction(operation, parsed);
  const multiply = parsed.op;
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  // Convert ord to bigint if provided
  let ordBig: bigint | undefined;
  if (ord !== undefined) {
    ordBig = toBigInt(ord);
  }

  // If order not provided, try to get it from the element (Sage: _ord_from_op)
  if (ordBig === undefined) {
    // Try to get order from base element
    if (isMult) {
      const elem = base as unknown as { multiplicative_order?(): bigint };
      if (typeof elem.multiplicative_order === 'function') {
        ordBig = elem.multiplicative_order();
      }
    } else if (isAdd) {
      const elem = base as unknown as { additive_order?(): bigint };
      if (typeof elem.additive_order === 'function') {
        ordBig = elem.additive_order();
      }
    }
  }

  if (ordBig === undefined) {
    throw new ValueError('ord must be specified when element order cannot be determined');
  }

  const ordFixed = ordBig;
  try {
    // base is the identity but a is not: no solution
    if (elementsEqual(base, power(base, 0n)) && !elementsEqual(a, base)) {
      throw new ValueError('no solution');
    }

    const result = _discrete_log_core(
      a,
      base,
      ordFixed,
      factor(ordFixed),
      { power, multiply },
      'other',
      parsed.identity,
      parsed.inverse,
      parsed.op
    );

    if (!elementsEqual(power(base, result), a)) {
      throw new ValueError('verification failed');
    }
    return result;
  } catch (e) {
    if (e instanceof ValueError) {
      throw new ValueError(
        `no discrete log of ${formatGroupElement(a)} found to base ${formatGroupElement(base)}`
      );
    }
    throw e;
  }
}

/**
 * Find the order of a group element given a multiple of its order.
 *
 * @param a - Group element
 * @param orderMultiple - A known multiple of the order (i.e., a^multiple = identity)
 * @param factorization - Factorization of multiple (computed if not provided)
 * @param operation - Type of group operation (default `'+'`, as in SageMath)
 * @param identity - Identity element (for custom operations)
 * @param inverse - Inverse function (for custom operations)
 * @param op - Binary operation (for custom operations)
 * @param options - `plist` (prime factors of the multiple, kept for
 *   compatibility) and `check` (default `true`): verify that `orderMultiple`
 *   really is a multiple of the order
 * @returns The exact order of a
 *
 * @example
 * ```typescript
 * const a = Mod(5n, 37n);
 * // We know 5^36 = 1 mod 37 (by Fermat's little theorem)
 * const order = order_from_multiple(a, 36n, undefined, '*');
 * console.log(order); // The multiplicative order of 5 mod 37
 * ```
 *
 * @see Reference: sage/groups/generic.py:order_from_multiple (line 1328)
 * @see Deviation: Generic Group API and DLP
 */
export function order_from_multiple<T extends GroupElement>(
  a: T,
  orderMultiple: IntegerLike,
  factorization?: Factorization,
  operation: OperationType = '+',
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T,
  options?: { plist?: IntegerLike[]; check?: boolean }
): bigint {
  const parsed = parseGroupOps(operation, identity, inverse, op, a);
  return orderFromMultipleParsed(
    a,
    orderMultiple,
    factorization,
    parsed.identity,
    groupPowerFunction(operation, parsed),
    options
  );
}

/** Reduce an order multiple using the identity and power already parsed by the caller. */
function orderFromMultipleParsed<T extends GroupElement>(
  a: T,
  orderMultiple: IntegerLike,
  factorization: Factorization | undefined,
  identity: T,
  power: (x: T, n: bigint) => T,
  options?: { plist?: IntegerLike[]; check?: boolean }
): bigint {
  const check = options?.check ?? true;
  const isId = (x: T) => elementsEqual(x, identity);

  // Identity has order 1
  if (isId(a)) {
    return 1n;
  }

  const orderMultipleBig = toBigInt(orderMultiple);
  if (check && !isId(power(a, orderMultipleBig))) {
    throw new ValueError(
      `The order of P(=${formatGroupElement(a)}) does not divide ${orderMultipleBig}`
    );
  }

  // Get factorization
  let factors: Factorization;
  if (factorization?.length) {
    factors = factorization;
  } else if (options?.plist?.length) {
    // Sage calls M.valuation(p); preserve its validation and GMP removal path.
    const valued = options.plist.map((p) => {
      const prime = toBigInt(p);
      return [prime, valuation(orderMultipleBig, prime)] as const;
    });
    if (valued[0]![1] === 'Infinity') {
      // M=0: every valuation is +Infinity. For a single base the native
      // helper repeatedly applies that base until it reaches the identity.
      if (valued.length === 1) {
        const prime = valued[0]![0];
        let Q = a,
          exponent = 0n;
        while (!isId(Q)) {
          Q = power(Q, prime);
          exponent++;
        }
        return prime ** exponent;
      }
      // The multi-factor helper cannot form its infinite cost/product.
      // infinity.py:400,452 and expression.pyx:685 define these diagnostics.
      for (const [prime] of valued) {
        if (prime < 0n) throw new TypeError('Python infinity cannot have complex phase.');
        if (prime <= 1n) {
          const error = new ArithmeticError(
            prime === 0n
              ? 'cannot add infinity to minus infinity'
              : 'cannot multiply infinity by zero'
          );
          error.name = 'SignError';
          throw error;
        }
      }
      throw new TypeError(
        "unsupported operand parent(s) for ^: 'The Infinity Ring' and 'The Infinity Ring'"
      );
    }
    factors = valued.map(([prime, exponent]) => [prime, exponent as bigint]);
  } else {
    factors = factor(orderMultipleBig);
  }

  // Filter out sign factor and convert to list form for the helper
  const L: Array<[bigint, bigint]> = [];
  for (const [p, e] of factors) {
    if (p !== -1n) {
      L.push([p, e]);
    }
  }

  // Special case: M itself is prime (or prime power with single factor)
  if (L.length === 1 && L[0]![0] === orderMultipleBig && L[0]![1] === 1n) {
    return orderMultipleBig;
  }

  // Compute total cost S = sum of e * log(p) for all factors
  // This represents the approximate "cost" of multiplications
  const totalCost = L.reduce((sum, [p, e]) => sum + Number(e) * Math.log(Number(p)), 0);

  /**
   * Internal recursive helper to minimize group operations.
   *
   * Uses cost-aware splitting to balance the work between left and right halves
   * of the factor list. The cost of a factor p^e is approximately e * log(p),
   * which represents the number of group operations needed.
   *
   * @param Q - Current group element
   * @param factorList - List of (prime, exponent) tuples to process
   * @param S - Sum of costs for factors in factorList
   * @returns The order contribution from these factors
   */
  function _order_from_multiple_helper(
    Q: T,
    factorList: Array<[bigint, bigint]>,
    S: number
  ): bigint {
    const l = factorList.length;

    if (l === 1) {
      // Base case: single prime factor
      // Determine the power of p dividing the order
      const [p, e] = factorList[0]!;
      let e0 = 0n;

      // Efficiency improvement: avoid the last multiplication by p.
      // For example, if M itself is prime, the code used to compute M*P
      // twice (unless P=0), now it does it once.
      while (!isId(Q) && e0 < e - 1n) {
        Q = power(Q, p);
        e0 += 1n;
      }
      if (!isId(Q)) {
        e0 += 1n;
      }
      return p ** e0;
    } else {
      // Recursive case: split the list to balance costs
      // Try to find k such that sum of costs for L[:k] is closest to S/2
      let sumLeft = 0;
      let k = 0;
      for (k = 0; k < l; k++) {
        const [p, e] = factorList[k]!;
        // Cost of p^e is approximately e * log(p)
        const v = Number(e) * Math.log(Number(p));
        // Check if adding this factor would take us farther from S/2
        if (Math.abs(sumLeft + v - S / 2) > Math.abs(sumLeft - S / 2)) {
          break;
        }
        sumLeft += v;
      }

      // Ensure we make progress (avoid empty splits)
      if (k <= 0 || k >= l) {
        k = Math.floor(l / 2);
      }

      const L1 = factorList.slice(0, k);
      const L2 = factorList.slice(k);

      // Compute product of p^e for all factors in L2
      let productL2 = 1n;
      for (const [p, e] of L2) {
        productL2 *= p ** e;
      }

      // Recursive calls:
      // First, compute order contribution from L1 factors
      // by multiplying Q by the product of L2 factors
      const o1 = _order_from_multiple_helper(power(Q, productL2), L1, sumLeft);

      // Then compute order contribution from L2 factors
      // by multiplying Q by o1 (the order from L1)
      const o2 = _order_from_multiple_helper(power(Q, o1), L2, S - sumLeft);

      return o1 * o2;
    }
  }

  return _order_from_multiple_helper(a, L, totalCost);
}

/**
 * Find the order of a group element given only upper and lower bounds
 * for a multiple of the order (e.g., bounds on the order of the group).
 *
 * Uses BSGS to find n with lb <= n <= ub such that n*P = identity,
 * then calls order_from_multiple() to get the exact order.
 *
 * @param P - A group element
 * @param bounds - A 2-tuple (lb, ub) such that m*P = identity for some m
 *   with lb <= m <= ub. If undefined, gradually increasing bounds will be
 *   tried (may loop infinitely if the element has no torsion).
 * @param d - Optional positive integer; only m which are multiples of d
 *   will be considered
 * @param operation - Type of group operation ('+', '*', or 'other')
 * @param identity - Identity element (for custom operations)
 * @param inverse - Inverse function (for custom operations)
 * @param op - Binary operation (for custom operations)
 * @returns The exact order of P
 * @throws {ValueError} If no suitable n found in the given bounds
 *
 * @example
 * ```typescript
 * // In GF(5^5)*, find the order of an element
 * // The group order is 5^5 - 1 = 3124
 * const b = Mod(3n, 3125n); // Example element
 * const order = order_from_bounds(b, [625n, 3125n], undefined, '*');
 *
 * // Without bounds - automatically increases search range
 * const order2 = order_from_bounds(b, undefined, undefined, '*');
 *
 * // With divisibility constraint
 * const order3 = order_from_bounds(b, [1n, 3125n], 7n, '*');
 * // Will only find orders that are multiples of 7
 * ```
 *
 * @see Reference: sage/groups/generic.py:order_from_bounds (lines 1476-1563)
 */
export function order_from_bounds<T extends GroupElement>(
  P: T,
  bounds?: [IntegerLike, IntegerLike],
  d?: IntegerLike,
  operation: OperationType = '+',
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): bigint {
  const parsed = parseGroupOps(operation, identity, inverse, op, P);
  return orderFromBoundsParsed(P, bounds, d, parsed, groupPowerFunction(operation, parsed));
}

/** Keep one parsed parent identity throughout exponential bounds and order reduction. */
function orderFromBoundsParsed<T extends GroupElement>(
  P: T,
  bounds: [IntegerLike, IntegerLike] | undefined,
  d: IntegerLike | undefined,
  parsed: GroupOps<T>,
  nativePower: (x: T, n: bigint) => T
): bigint {
  const { power, identity: identityElem } = parsed;

  // Handle bounds=undefined case: gradually increase bounds
  if (bounds === undefined) {
    let lb = 1n;
    let ub = 256n;
    while (true) {
      try {
        return orderFromBoundsParsed(P, [lb, ub], d, parsed, nativePower);
      } catch (e) {
        if (e instanceof ValueError) {
          lb = ub + 1n;
          ub *= 16n;
        } else {
          throw e;
        }
      }
    }
  }

  // Parse bounds
  const [lbInput, ubInput] = bounds;
  let lb = toBigInt(lbInput);
  let ub = toBigInt(ubInput);

  // Handle d parameter
  let Q = P;
  const dBig = d !== undefined ? toBigInt(d) : 1n;

  if (dBig > 1n) {
    // Q = d*P
    Q = power(P, dBig);
    // Adjust bounds: divide by d with ceiling/floor
    // We need to find m such that lb <= d*m <= ub
    // So ceiling(lb/d) <= m <= floor(ub/d)
    // Sage divides Integer bounds exactly before applying ceil/floor.
    lb = lb / dBig + (lb % dBig > 0n ? 1n : 0n);
    ub = ub / dBig - (ub % dBig < 0n ? 1n : 0n);
  }

  // Use bsgs to find n = d*m with lb <= n <= ub and n*P = identity
  const m = bsgs(Q, identityElem, [lb, ub], 'other', identityElem, parsed.inverse, parsed.op);
  const n = dBig * m;

  // Use order_from_multiple to find exact order
  return orderFromMultipleParsed(P, n, undefined, identityElem, nativePower, { check: false });
}

/**
 * Find a multiple of the order of a group element.
 *
 * This is useful when the group order is not known exactly, but we can
 * compute powers and test for identity.
 *
 * @param a - Group element
 * @param operation - Type of group operation
 * @param identity - Identity element (for custom operations)
 * @param inverse - Inverse function (for custom operations)
 * @param op - Binary operation (for custom operations)
 * @param maxIterations - Maximum number of iterations (default: 2^20)
 * @returns Some m with a^m = identity
 * @throws {ValueError} If no multiple found within maxIterations
 *
 * @example
 * ```typescript
 * const a = Mod(5n, 37n);
 * const mult = multiple_of_order(a, '*');
 * // mult is some value m where 5^m = 1 mod 37
 * ```
 */
export function multiple_of_order<T extends GroupElement>(
  a: T,
  operation: OperationType = '*',
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T,
  maxIterations: IntegerLike = 1n << 20n
): bigint {
  const maxIterationsBig = toBigInt(maxIterations);
  validateGroupOps(operation, identity, inverse, op);

  // Use order_from_bounds with undefined bounds, which implements an efficient
  // exponential search strategy using BSGS. This is O(sqrt(order)) per range
  // with exponentially increasing ranges, much better than O(order) linear search.
  //
  // The algorithm:
  // 1. Try ranges [1, 256], [257, 4096], [4097, 65536], ... (16x growth)
  // 2. Within each range, use baby-step giant-step which is O(sqrt(range))
  // 3. Once a multiple m is found, use order_from_multiple to get exact order
  //
  // For an element of order n:
  // - Linear search (old): O(n) group operations
  // - BSGS-based (new): O(sqrt(n)) group operations
  //
  // The returned value is the exact order, which is always a valid multiple.
  //
  // Note: maxIterations is kept for API compatibility but order_from_bounds
  // will keep searching until it finds the order (or forever for infinite order).
  try {
    return order_from_bounds(a, undefined, undefined, operation, identity, inverse, op);
  } catch (e) {
    if (e instanceof ValueError) {
      throw new ValueError(
        `could not find multiple of order within ${maxIterationsBig} iterations: ` +
          'element may have infinite order'
      );
    }
    throw e;
  }
}

/**
 * Test if a group element `P` has order exactly equal to a given positive
 * integer `n`.
 *
 * In some cases order *testing* can be much faster than *computing* the order
 * using {@link order_from_multiple}; this uses Sage's divide-and-conquer
 * recursion over the factorization of `n`.
 *
 * @param P - Group element
 * @param n - Proposed order, or its factorization
 * @param operation - Type of group operation (default `'+'`, as in SageMath)
 * @returns true if the order of P is exactly n
 *
 * @example
 * ```typescript
 * const a = Mod(2n, 7n);
 * has_order(a, 3n, '*'); // true (2^3 = 8 = 1 mod 7)
 * has_order(a, 6n, '*'); // false
 * ```
 *
 * @see Reference: sage/groups/generic.py:has_order (line 1566)
 * @see Deviation: Generic Group API and DLP
 */
export function has_order<T extends GroupElement>(
  P: T,
  n: IntegerLike | Factorization,
  operation: OperationType = '+'
): boolean {
  let fn: Array<[bigint, bigint]>;
  if (Array.isArray(n)) {
    fn = n.map(([p, e]) => [toBigInt(p), toBigInt(e)] as [bigint, bigint]);
  } else {
    const nBig = toBigInt(n);
    if (nBig <= 0n) {
      return false;
    }
    fn = factor(nBig).filter(([p]) => p > 0n);
  }

  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let mult: (x: T, k: bigint) => T;
  let isId: (x: T) => boolean;

  if (isAdd) {
    mult = (x: T, k: bigint) => multiple(x, k, '+');
    isId = (x: T) => isStandardIdentity(x, false);
  } else if (isMult) {
    mult = (x: T, k: bigint) => multiple(x, k, '*');
    isId = (x: T) => isStandardIdentity(x, true);
  } else {
    throw new ValueError('unknown group operation');
  }

  const _rec = (Q: T, factors: Array<[bigint, bigint]>): boolean => {
    if (factors.length === 0) {
      return isId(Q);
    }

    if (factors.length === 1) {
      const [p, k] = factors[0]!;
      let R = Q;
      for (let i = 0n; i < k; i++) {
        if (isId(R)) {
          return false;
        }
        R = mult(R, p);
      }
      return isId(R);
    }

    const fl = factors.filter((_, idx) => idx % 2 === 0);
    const fr = factors.filter((_, idx) => idx % 2 === 1);
    let left = 1n;
    for (const [p, k] of fl) left *= p ** k;
    let right = 1n;
    for (const [p, k] of fr) right *= p ** k;
    const L = mult(Q, right);
    const R = mult(Q, left);
    return _rec(L, fl) && _rec(R, fr);
  };

  return _rec(P, fn);
}

/**
 * Shallow copy adapter for Python copy.copy in the multiples constructor.
 * @see Deviation: Generic group iterator and parent adapters
 */
function copyGroupElement<T>(value: T): T {
  if (value === null || typeof value !== 'object') return value;
  const hooks = value as { __copy__?: () => T; copy?: () => T };
  if (typeof hooks.__copy__ === 'function') return hooks.__copy__();
  if (typeof hooks.copy === 'function') return hooks.copy();
  return Object.create(Object.getPrototypeOf(value), Object.getOwnPropertyDescriptors(value)) as T;
}

/**
 * Iterator for multiples/powers of a group element.
 *
 * For additive groups: yields P0, P0+P, P0+2P, ..., P0+(n-1)P
 * For multiplicative groups: yields P0, P0*P, P0*P^2, ..., P0*P^(n-1)
 *
 * @param P - Step element
 * @param n - Number of multiples to generate
 * @param P0 - Offset (default: identity)
 * @param indexed - If true, yield `[i, element]` pairs instead of elements
 * @param operation - Type of group operation (default `'+'`, as in SageMath)
 * @param op - Binary operation (required when operation is neither '+' nor '*')
 * @returns Generator yielding `P0 + i*P` (or `[i, P0 + i*P]` when indexed)
 *
 * @see Reference: sage/groups/generic.py:multiples (line 388)
 * @see Deviation: Generic group iterator and parent adapters
 */
export function multiples<T extends GroupElement>(
  P: T,
  n: IntegerLike,
  P0: T | undefined,
  indexed: true,
  operation?: OperationType,
  op?: (x: T, y: T) => T
): Generator<[bigint, T], void, undefined>;
export function multiples<T extends GroupElement>(
  P: T,
  n: IntegerLike,
  P0?: T,
  indexed?: false,
  operation?: OperationType,
  op?: (x: T, y: T) => T
): Generator<T, void, undefined>;
export function multiples<T extends GroupElement>(
  P: T,
  n: IntegerLike,
  P0?: T,
  indexed = false,
  operation: OperationType = '+',
  op?: (x: T, y: T) => T
): Generator<T | [bigint, T], void, undefined> {
  const nBig = toBigInt(n);
  if (nBig < 0n) {
    throw new ValueError('n cannot be negative in multiples');
  }

  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let multiply: (x: T, y: T) => T;
  let start: T;

  if (isMult) {
    multiply = (x: T, y: T) =>
      (x as unknown as MultiplicativeGroupElement).mul(
        y as unknown as MultiplicativeGroupElement
      ) as unknown as T;
    const parent = P0 === undefined ? getGroupParent(P) : undefined;
    start =
      P0 !== undefined
        ? P0
        : ((parent && typeof parent.one === 'function'
            ? parent.one()
            : (P as unknown as MultiplicativeGroupElement).pow(0n)) as T);
  } else if (isAdd) {
    multiply = (x: T, y: T) =>
      (x as unknown as AdditiveGroupElement).add(
        y as unknown as AdditiveGroupElement
      ) as unknown as T;
    const parent = P0 === undefined ? getGroupParent(P) : undefined;
    start =
      P0 !== undefined
        ? P0
        : ((parent && typeof parent.zero === 'function'
            ? parent.zero()
            : (P as unknown as AdditiveGroupElement).mul(0n)) as T);
  } else {
    if (P0 === undefined) {
      throw new ValueError(
        'P0 must be supplied when operation is neither addition nor multiplication'
      );
    }
    if (op === undefined) {
      throw new ValueError(
        'op() must both be supplied when operation is neither addition nor multiplication'
      );
    }
    multiply = op;
    start = P0;
  }

  // generic.py:480–505: copy at construction and advance before returning.
  const step = copyGroupElement(P);
  let current = copyGroupElement(start);
  if (step == null || current == null) throw new ValueError('P and Q must not be None');
  let index = 0n;
  let closed = false;
  // Retain the existing Generator interface and inherited iterator helpers.
  // The native state machine survives a callback error; a generator body cannot.
  const iterator = (function* (): Generator<T | [bigint, T], void, undefined> {})();
  iterator.next = () => {
    if (closed || index >= nBig) return { value: undefined, done: true };
    const i = index++,
      value = current;
    current = multiply(current, step);
    return { value: indexed ? [i, value] : value, done: false };
  };
  iterator.return = (value) => {
    closed = true;
    return { value, done: true };
  };
  iterator.throw = (error) => {
    closed = true;
    throw error;
  };
  return iterator;
}

/**
 * Pollard Lambda (Kangaroo) algorithm for computing discrete logarithm
 * when the logarithm is known to lie in a bounded interval [lb, ub].
 *
 * Uses only O(log(ub-lb)) memory and O(sqrt(ub-lb)) time.
 *
 * @param a - Target element (find x such that base^x = a)
 * @param base - Base element
 * @param bounds - Tuple [lb, ub] with 0 <= lb <= ub
 * @param operation - Type of group operation
 * @param identity - Identity element (for custom operations)
 * @param inverse - Inverse function (for custom operations)
 * @param op - Binary operation (for custom operations)
 * @param hashFunction - Custom hash function (defaults to elementToString-based hash)
 * @returns x such that base^x = a with lb <= x <= ub
 * @throws {ValueError} If no such x exists or algorithm fails to find solution
 *
 * @example
 * ```typescript
 * // In (Z/37Z)*
 * const base = Mod(2n, 37n);
 * const a = base.pow(20n);
 * // We know the log is between 15 and 25
 * const x = discrete_log_lambda(a, base, [15n, 25n], '*'); // returns 20n
 * ```
 *
 * @see SageMath reference: sage/groups/generic.py discrete_log_lambda
 */
export function discrete_log_lambda<T extends GroupElement>(
  a: T,
  base: T,
  bounds: [IntegerLike, IntegerLike],
  operation: OperationType = '*',
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T,
  hashFunction?: (x: T) => bigint
): bigint {
  const parsed = parseGroupOps(operation, identity, inverse, op, a);
  const power = groupPowerFunction(operation, parsed);
  const multiply = parsed.op;
  const [lbInput, ubInput] = bounds;
  const lb = toBigInt(lbInput);
  const ub = toBigInt(ubInput);

  if (lb < 0n || ub < lb) {
    throw new ValueError('discrete_log_lambda() requires 0<=lb<=ub');
  }

  const mutable = 'set_immutable' in Object(base);
  const freeze = (x: T) => {
    if (mutable) (x as unknown as { set_immutable(): void }).set_immutable();
  };

  // Default hash function
  const hash =
    hashFunction ??
    ((x: T) => {
      const s = elementToString(x);
      // Simple hash: sum of char codes with position weighting
      let h = 0n;
      for (let i = 0; i < s.length; i++) {
        h = (h * 31n + BigInt(s.charCodeAt(i))) % (1n << 62n);
      }
      return h;
    });

  const hashIndex = (x: T, k: bigint): number => {
    const h = hash(x);
    if (k === 0n) throw new ZeroDivisionError('integer modulo by zero');
    return Number(((h % k) + k) % k);
  };

  const width = ub - lb;
  const N = isqrt(width) + 1n;
  const randstate = current_randstate();

  const M = new Map<number, [bigint, T]>();

  // Retry loop to handle random walk failures
  for (let attempt = 0; attempt < 10; attempt++) {
    // Set up random walk function
    // We need k values r_i and corresponding base^r_i
    // k should be chosen such that 2^k >= N
    let k = 0n;
    while (1n << k < N) {
      k++;
    }
    const kInt = Number(k);

    // generic.py:1214: prandom uses the state's CPython stream, not GMP.
    // Random step sizes r_i, Sage: randrange(1, N)
    const maxR = N - 1n;
    for (let i = 0; i < kInt; i++) {
      const r = randstate.python_random().randint(1n, maxR);
      const e = power(base, r);
      M.set(i, [r, e]);
    }

    // First random walk: "tame" kangaroo
    // Starts at base^ub, takes N steps
    let H = power(base, ub);
    let c = ub;
    for (let i = 0n; i < N; i++) {
      freeze(H);
      const hashIdx = hashIndex(H, k);
      const [r, e] = M.get(hashIdx)!;
      H = multiply(H, e);
      c = c + r;
    }

    // H is now the "trap": H = base^c where c = ub + (sum of steps)
    freeze(H);
    const mem = H;

    // Second random walk: "wild" kangaroo
    // Starts at a, tries to reach the same position as tame kangaroo
    H = a;
    let d = 0n;

    // Walk until we either find the trap or overshoot
    while (c - d >= lb) {
      freeze(H);
      // Check if we've found the trap
      if (ub >= c - d && elementsEqual(H, mem)) {
        // H = a * base^d = base^(c) = mem
        // So a = base^(c - d)
        return c - d;
      }

      const hashIdx = hashIndex(H, k);
      const [r, e] = M.get(hashIdx)!;
      H = multiply(H, e);
      d = d + r;
    }
  }

  throw new ValueError('Pollard Lambda failed to find a log');
}

/**
 * DJB2 hash function for group elements (used by Pollard rho).
 * @internal
 */
function rhoHash<T>(elem: T): number {
  const str = elementToString(elem);
  let hash = 5381;
  for (let i = 0; i < str.length; i++) {
    hash = ((hash << 5) + hash) ^ str.charCodeAt(i);
    hash = hash >>> 0; // Convert to unsigned 32-bit
  }
  return hash;
}

/**
 * Pollard Rho algorithm for computing discrete logarithm in a cyclic group
 * of prime order.
 *
 * Uses Teske's space-efficient collision detection algorithm.
 * Falls back to BSGS for small orders (isqrt(ord) < 20).
 *
 * @param a - Target element (find x such that base^x = a)
 * @param base - Base element (generator of a prime order group)
 * @param ord - Prime order of the group (required)
 * @param operation - Type of group operation
 * @param identity - Identity element (for custom operations)
 * @param inverse - Inverse function (for custom operations)
 * @param op - Binary operation (for custom operations)
 * @returns x such that base^x = a
 * @throws {ValueError} If ord is not prime or if no solution exists
 *
 * @example
 * ```typescript
 * // In (Z/37Z)* (order 36 is not prime, so we need a prime order subgroup)
 * // For a prime order p, this works directly:
 * const base = Mod(4n, 1019n); // order 509, a prime
 * const target = base.pow(250n);
 * const x = discrete_log_rho(target, base, 509n, '*'); // 250n
 * ```
 *
 * @see SageMath reference: sage/groups/generic.py discrete_log_rho
 * @see Deviation: Generic Group API and DLP
 */
export function discrete_log_rho<T extends GroupElement>(
  a: T,
  base: T,
  ord: IntegerLike,
  operation: OperationType = '*',
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): bigint {
  const parsed = parseGroupOps(operation, identity, inverse, op, a);
  const power = groupPowerFunction(operation, parsed);
  const multiply = parsed.op;
  const ordBig = toBigInt(ord);

  // Check that order is prime
  if (!is_prime(ordBig)) {
    throw new ValueError('for Pollard rho algorithm the order of the group must be prime');
  }

  const mutable = 'set_immutable' in Object(base);
  const freeze = (x: T) => {
    if (mutable) (x as unknown as { set_immutable(): void }).set_immutable();
  };

  const isqrtord = isqrt(ordBig);

  // SageMath parameters
  const partitionSize = 20;
  const memorySize = 4;

  // Fall back to BSGS for small orders
  if (isqrtord < BigInt(partitionSize)) {
    return bsgs(base, a, [0n, ordBig], 'other', parsed.identity, parsed.inverse, parsed.op);
  }

  // Reset bound: 8 * sqrt(ord) iterations per attempt
  const resetBound = 8n * isqrtord;

  // Modular arithmetic on exponents
  const exponentRing = new IntegerModRing(ordBig); // ordBig is a positive prime here.
  const I = (x: bigint) => exponentRing.__call__(x);

  // Outer loop for multiple attempts (to avoid infinite loops)
  for (let s = 0; s < 10; s++) {
    // Setup random walk function
    // m[i], n[i] are random exponents for partition i
    // M[i] = base^m[i] * a^n[i] is the multiplier for partition i
    const m: bigint[] = [];
    const n: bigint[] = [];
    const M: T[] = [];

    // Sage draws all m values, then all n values, from IntegerModRing.
    for (let i = 0; i < partitionSize; i++) m.push(exponentRing.random_element().value);
    for (let i = 0; i < partitionSize; i++) n.push(exponentRing.random_element().value);
    for (let i = 0; i < partitionSize; i++) {
      M.push(multiply(power(base, m[i]!), power(a, n[i]!)));
    }

    // Initial random point: x = base^ax
    let ax = exponentRing.random_element();
    let x = power(base, ax.value);
    freeze(x);
    let bx = I(0n);

    type Entry = { element: T; exponents: [bigint, bigint] };
    // Keep each stored element so colliding string keys can be removed separately.
    const sigma: Array<[number, { key: string; entry: Entry } | null]> = [];
    for (let i = 0; i < memorySize; i++) {
      sigma.push([0, null]);
    }

    // H is the hash table storing (ax, bx) for each element seen
    const H = new Map<string, Entry[]>();

    let i0 = 0;
    let nextsigma = 0;

    for (let i = 0; i < Number(resetBound); i++) {
      // Random walk step
      const hashIdx = rhoHash(x) % partitionSize;
      x = multiply(M[hashIdx]!, x);
      ax = ax.add(I(m[hashIdx]!));
      bx = bx.add(I(n[hashIdx]!));
      freeze(x);

      // Get string key for element
      const xKey = elementToString(x);

      // Look for collisions
      const stored = H.get(xKey)?.find(({ element }) => element === x || elementsEqual(element, x));
      if (stored !== undefined) {
        const [ay, by] = stored.exponents;
        const bxVal = bx.value;
        if (bxVal === by) {
          // Same beta, can't solve - break and restart
          break;
        } else {
          // Found collision: base^ax * a^bx = base^ay * a^by
          // base^(ay - ax) = a^(bx - by)
          // log(a) = (ay - ax) / (bx - by) mod ord
          const ayMinusAx = I(ay - ax.value);
          const bxMinusBy = I(bxVal - by);
          const result = ayMinusAx.mul(bxMinusBy.inv()).value;

          // Verify
          if (elementsEqual(power(base, result), a)) {
            return result;
          } else {
            // Wrong result, break and restart
            break;
          }
        }
      }

      // Should we remember this value?
      if (i >= nextsigma) {
        // Remove old value from H
        const old = sigma[i0]![1];
        if (old !== null) {
          const bucket = H.get(old.key)!.filter((entry) => entry !== old.entry);
          if (bucket.length) H.set(old.key, bucket);
          else H.delete(old.key);
        }

        // Store current value
        const entry: Entry = { element: x, exponents: [ax.value, bx.value] };
        sigma[i0] = [i, { key: xKey, entry }];
        i0 = (i0 + 1) % memorySize;

        // Next storage time: 3 * oldest stored iteration
        // This spreads out storage points geometrically
        nextsigma = 3 * sigma[i0]![0];

        const bucket = H.get(xKey) ?? [];
        bucket.push(entry);
        H.set(xKey, bucket);
      }
    }
  }

  throw new ValueError('Pollard rho algorithm failed to find a logarithm');
}

/**
 * Alias for discrete_log for backward compatibility.
 */
export const discrete_log_generic = discrete_log;
