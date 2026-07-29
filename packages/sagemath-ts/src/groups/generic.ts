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
  integer_ceil,
  integer_floor,
  is_prime,
  isqrt,
} from '../arith/misc.js';
import { ValueError } from '../errors.js';
import { current_randstate } from '../misc/randstate.js';
import { Mod } from '../rings/finite_rings/integer_mod.js';
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

function assertNoCustomOps<T>(
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
}

/**
 * Interface for group elements that support generic group operations.
 */
export interface GroupElement {
  eq(other: GroupElement): boolean;
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

function deriveStandardGroupOps<T extends GroupElement>(
  sample: T,
  operation: OperationType
): GroupOps<T> {
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);
  if (!isMult && !isAdd) {
    throw new ValueError(
      "identity, inverse and operation must all be specified when operation is 'other'"
    );
  }

  const parent = (sample as unknown as { parent?: { one?: () => unknown; zero?: () => unknown } })
    .parent;
  let identity: T | undefined;

  if (isMult && parent && typeof parent.one === 'function') {
    identity = parent.one() as T;
  } else if (isAdd && parent && typeof parent.zero === 'function') {
    identity = parent.zero() as T;
  }

  if (identity === undefined) {
    if (typeof sample === 'bigint') {
      identity = (isMult ? 1n : 0n) as unknown as T;
    } else if (typeof sample === 'number') {
      identity = (isMult ? 1 : 0) as unknown as T;
    }
  }

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
    ? (x: T) => {
        if ((x as unknown as MultiplicativeGroupElement).inv) {
          return (x as unknown as MultiplicativeGroupElement).inv() as unknown as T;
        }
        if (typeof x === 'number') {
          return (1 / (x as unknown as number)) as unknown as T;
        }
        throw new ValueError('Cannot invert element');
      }
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
          return (x as unknown as MultiplicativeGroupElement).isOne();
        }
        return elementsEqual(x, identity as T);
      }
    : (x: T) => {
        if ((x as unknown as AdditiveGroupElement).isZero) {
          return (x as unknown as AdditiveGroupElement).isZero();
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
 * @returns Normalized group operations
 */
export function parseGroupOps<T extends GroupElement>(
  operation: OperationType,
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T,
  sample?: T
): GroupOps<T> {
  if (isMultiplicative(operation) || isAdditive(operation)) {
    assertNoCustomOps(operation, identity, inverse, op);
    if (sample === undefined) {
      throw new ValueError('identity could not be determined for standard operation');
    }
    return deriveStandardGroupOps(sample, operation);
  }

  // Custom operation
  if (identity === undefined || inverse === undefined || op === undefined) {
    throw new ValueError(
      "identity, inverse and operation must all be specified when operation is 'other'"
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
export function multiple<T extends GroupElement>(
  a: T,
  n: IntegerLike,
  operation: OperationType = '*',
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): T {
  let nBig = toBigInt(n);
  assertNoCustomOps(operation, identity, inverse, op);
  // For multiplicative groups, use native pow if available
  if (isMultiplicative(operation)) {
    const elem = a as MultiplicativeGroupElement;
    if (typeof elem.pow === 'function') {
      return elem.pow(nBig) as unknown as T;
    }
  }

  // For additive groups, use native scalar multiplication if available
  if (isAdditive(operation)) {
    const elem = a as AdditiveGroupElement;
    if (typeof elem.mul === 'function' && elem.mul.length === 1) {
      return elem.mul(nBig) as unknown as T;
    }
  }

  // Fall back to generic binary algorithm
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

  // Binary exponentiation
  let result = ops.identity;
  let base = a;

  while (nBig > 0n) {
    if ((nBig & 1n) === 1n) {
      result = ops.op(result, base);
    }
    base = ops.op(base, base);
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
  assertNoCustomOps(operation, identity, inverse, op);
  const [lbInput, ubInput] = bounds;
  const lb = toBigInt(lbInput);
  const ub = toBigInt(ubInput);

  if (lb < 0n || ub < lb) {
    throw new ValueError('bsgs() requires 0<=lb<=ub');
  }

  const range = 1n + ub - lb;

  // For multiplicative groups, use native pow directly
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  // Define group operations based on type
  let power: (x: T, n: bigint) => T;
  let multiply: (x: T, y: T) => T;
  let invert: (x: T) => T;
  let isId: (x: T) => boolean;

  if (isMult) {
    power = (x: T, n: bigint) =>
      (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    multiply = (x: T, y: T) =>
      (x as unknown as MultiplicativeGroupElement).mul(
        y as unknown as MultiplicativeGroupElement
      ) as unknown as T;
    invert = (x: T) => (x as unknown as MultiplicativeGroupElement).inv() as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    multiply = (x: T, y: T) =>
      (x as unknown as AdditiveGroupElement).add(
        y as unknown as AdditiveGroupElement
      ) as unknown as T;
    invert = (x: T) => (x as unknown as AdditiveGroupElement).neg() as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError(
        "identity, inverse and operation must all be specified when operation is 'other'"
      );
    }
    power = (x: T, n: bigint) => multiple(x, n, operation, identity, inverse, op);
    multiply = op;
    invert = inverse;
    isId = (x: T) => x.eq(identity);
  }

  // Handle identity base specially
  if (isId(a) && !isId(b)) {
    throw new ValueError('no solution in bsgs()');
  }

  // Compute b^(-1) * a^lb
  const aLb = power(a, lb);
  const bInv = invert(b);
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
  const table = new Map<string, bigint>();
  let d = c;

  for (let i = 0n; i < m; i++) {
    if (isId(d)) {
      return lb + i;
    }
    // Use string representation as hash key
    const key = elementToString(d);
    table.set(key, lb + i);
    d = multiply(d, a);
  }

  // Giant steps: we have d = c * a^m = b^(-1) * a^(lb+m)
  // We need a^(-m) for giant steps
  const am = power(a, m);
  const aInvM = invert(am);

  // Start from identity and multiply by a^(-m) each step
  // We're looking for: identity * a^(-km) = c * a^j for some j in table
  // Which means: a^(-km) = b^(-1) * a^(lb+j)
  // So: b = a^(km + lb + j)
  let giant = power(a, 0n); // Start at identity (a^0 = 1)

  for (let k = 0n; k < m; k++) {
    const key = elementToString(giant);
    const j = table.get(key);
    if (j !== undefined) {
      return k * m + j;
    }
    giant = multiply(giant, aInvM);
  }

  throw new ValueError(`log of ${b} to the base ${a} does not exist in (${lb}, ${ub})`);
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
    invert: (x: T) => T;
    isId: (x: T) => boolean;
  },
  operation: OperationType,
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): bigint {
  const { power, multiply, invert, isId } = ops;
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
    while (isId(gamma) && ri > 0n) {
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
      const h = power(multiply(a, invert(power(base, l[idx]!))), ord / pi ** (jj + 1n));
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
  assertNoCustomOps(operation, identity, inverse, op);
  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let power: (x: T, n: bigint) => T;
  let multiply: (x: T, y: T) => T;
  let invert: (x: T) => T;
  let isId: (x: T) => boolean;

  if (isMult) {
    power = (x: T, n: bigint) =>
      (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    multiply = (x: T, y: T) =>
      (x as unknown as MultiplicativeGroupElement).mul(
        y as unknown as MultiplicativeGroupElement
      ) as unknown as T;
    invert = (x: T) => (x as unknown as MultiplicativeGroupElement).inv() as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    multiply = (x: T, y: T) =>
      (x as unknown as AdditiveGroupElement).add(
        y as unknown as AdditiveGroupElement
      ) as unknown as T;
    invert = (x: T) => (x as unknown as AdditiveGroupElement).neg() as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError(
        "identity, inverse and operation must all be specified when operation is 'other'"
      );
    }
    power = (x: T, n: bigint) => multiple(x, n, operation, identity, inverse, op);
    multiply = op;
    invert = inverse;
    isId = (x: T) => x.eq(identity);
  }

  // Get factorization of n
  const factors = factorization ?? factor(nBig);

  return _discrete_log_core(
    a,
    base,
    nBig,
    factors,
    { power, multiply, invert, isId },
    operation,
    identity,
    inverse,
    op
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
  assertNoCustomOps(operation, identity, inverse, op);
  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let power: (x: T, n: bigint) => T;
  let multiply: (x: T, y: T) => T;
  let invert: (x: T) => T;
  let isId: (x: T) => boolean;

  if (isMult) {
    power = (x: T, n: bigint) =>
      (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    multiply = (x: T, y: T) =>
      (x as unknown as MultiplicativeGroupElement).mul(
        y as unknown as MultiplicativeGroupElement
      ) as unknown as T;
    invert = (x: T) => (x as unknown as MultiplicativeGroupElement).inv() as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    multiply = (x: T, y: T) =>
      (x as unknown as AdditiveGroupElement).add(
        y as unknown as AdditiveGroupElement
      ) as unknown as T;
    invert = (x: T) => (x as unknown as AdditiveGroupElement).neg() as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError(
        "identity, inverse and operation must all be specified when operation is 'other'"
      );
    }
    power = (x: T, n: bigint) => multiple(x, n, operation, identity, inverse, op);
    multiply = op;
    invert = inverse;
    isId = (x: T) => x.eq(identity);
  }

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
    if (isId(base) && !elementsEqual(a, base)) {
      throw new ValueError('no solution');
    }

    const result = _discrete_log_core(
      a,
      base,
      ordFixed,
      factor(ordFixed),
      { power, multiply, invert, isId },
      operation,
      identity,
      inverse,
      op
    );

    if (!elementsEqual(power(base, result), a)) {
      throw new ValueError('verification failed');
    }
    return result;
  } catch (e) {
    if (e instanceof ValueError) {
      throw new ValueError(`no discrete log of ${a} found to base ${base}`);
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
  const orderMultipleBig = toBigInt(orderMultiple);
  const check = options?.check ?? true;
  assertNoCustomOps(operation, identity, inverse, op);
  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let power: (x: T, n: bigint) => T;
  let isId: (x: T) => boolean;

  if (isMult) {
    power = (x: T, n: bigint) =>
      (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError(
        "identity, inverse and operation must all be specified when operation is 'other'"
      );
    }
    power = (x: T, n: bigint) => multiple(x, n, operation, identity, inverse, op);
    isId = (x: T) => x.eq(identity);
  }

  // Identity has order 1
  if (isId(a)) {
    return 1n;
  }

  if (check && !isId(power(a, orderMultipleBig))) {
    throw new ValueError(`The order of P(=${a}) does not divide ${orderMultipleBig}`);
  }

  // Get factorization
  let factors: Factorization;
  if (factorization) {
    factors = factorization;
  } else if (options?.plist) {
    factors = options.plist.map((p) => {
      const _p = toBigInt(p);
      let e = 0n;
      let m = orderMultipleBig;
      while (m % _p === 0n) {
        m /= _p;
        e += 1n;
      }
      return [_p, e] as [bigint, bigint];
    });
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
  assertNoCustomOps(operation, identity, inverse, op);

  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let power: (x: T, n: bigint) => T;
  let getIdentity: () => T;

  if (isMult) {
    power = (x: T, n: bigint) =>
      (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    getIdentity = () => (P as unknown as MultiplicativeGroupElement).pow(0n) as unknown as T;
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    getIdentity = () => (P as unknown as AdditiveGroupElement).mul(0n) as unknown as T;
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError(
        "identity, inverse and operation must all be specified when operation is 'other'"
      );
    }
    power = (x: T, n: bigint) => multiple(x, n, operation, identity, inverse, op);
    getIdentity = () => identity;
  }

  // Handle bounds=undefined case: gradually increase bounds
  if (bounds === undefined) {
    let lb = 1n;
    let ub = 256n;
    while (true) {
      try {
        return order_from_bounds(P, [lb, ub], d, operation, identity, inverse, op);
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

  // Get identity for bsgs
  const identityElem = identity ?? getIdentity();

  // Handle d parameter
  let Q = P;
  const dBig = d !== undefined ? toBigInt(d) : 1n;

  if (dBig > 1n) {
    // Q = d*P
    Q = power(P, dBig);
    // Adjust bounds: divide by d with ceiling/floor
    // We need to find m such that lb <= d*m <= ub
    // So ceiling(lb/d) <= m <= floor(ub/d)
    const lbNum = Number(lb);
    const ubNum = Number(ub);
    const dNum = Number(dBig);

    // For large numbers, use bigint division with ceiling/floor semantics
    if (
      lb > Number.MAX_SAFE_INTEGER ||
      ub > Number.MAX_SAFE_INTEGER ||
      dBig > Number.MAX_SAFE_INTEGER
    ) {
      // ceiling(lb/d) = (lb + d - 1) / d for positive integers
      lb = (lb + dBig - 1n) / dBig;
      // floor(ub/d) = ub / d (integer division)
      ub = ub / dBig;
    } else {
      lb = integer_ceil(lbNum / dNum);
      ub = integer_floor(ubNum / dNum);
    }
  }

  // Use bsgs to find n = d*m with lb <= n <= ub and n*P = identity
  const m = bsgs(Q, identityElem, [lb, ub], operation, identity, inverse, op);
  const n = dBig * m;

  // Use order_from_multiple to find exact order
  return order_from_multiple(P, n, undefined, operation, identity, inverse, op);
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
  assertNoCustomOps(operation, identity, inverse, op);

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
    mult = (x: T, k: bigint) => (x as unknown as AdditiveGroupElement).mul(k) as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else if (isMult) {
    mult = (x: T, k: bigint) => (x as unknown as MultiplicativeGroupElement).pow(k) as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
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
    return _rec(mult(Q, right), fl) && _rec(mult(Q, left), fr);
  };

  return _rec(P, fn);
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
 * @see Reference: sage/groups/generic.py:multiples (line 355)
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
export function* multiples<T extends GroupElement>(
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
    start = P0 ?? ((P as unknown as MultiplicativeGroupElement).pow(0n) as unknown as T);
  } else if (isAdd) {
    multiply = (x: T, y: T) =>
      (x as unknown as AdditiveGroupElement).add(
        y as unknown as AdditiveGroupElement
      ) as unknown as T;
    start = P0 ?? ((P as unknown as AdditiveGroupElement).mul(0n) as unknown as T);
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

  let current = start;

  for (let i = 0n; i < nBig; i++) {
    yield indexed ? [i, current] : current;
    current = multiply(current, P);
  }
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
  assertNoCustomOps(operation, identity, inverse, op);
  const [lbInput, ubInput] = bounds;
  const lb = toBigInt(lbInput);
  const ub = toBigInt(ubInput);

  if (lb < 0n || ub < lb) {
    throw new ValueError('discrete_log_lambda() requires 0 <= lb <= ub');
  }

  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let power: (x: T, n: bigint) => T;
  let multiply: (x: T, y: T) => T;

  if (isMult) {
    power = (x: T, n: bigint) =>
      (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    multiply = (x: T, y: T) =>
      (x as unknown as MultiplicativeGroupElement).mul(
        y as unknown as MultiplicativeGroupElement
      ) as unknown as T;
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    multiply = (x: T, y: T) =>
      (x as unknown as AdditiveGroupElement).add(
        y as unknown as AdditiveGroupElement
      ) as unknown as T;
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError(
        "identity, inverse and operation must all be specified when operation is 'other'"
      );
    }
    power = (x: T, n: bigint) => multiple(x, n, operation, identity, inverse, op);
    multiply = op;
  }

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

  const width = ub - lb;
  const N = isqrt(width) + 1n;
  const randstate = current_randstate();

  // Retry loop to handle random walk failures
  for (let attempt = 0; attempt < 10; attempt++) {
    // Set up random walk function
    // We need k values r_i and corresponding base^r_i
    // k should be chosen such that 2^k >= N
    let k = 0n;
    while (1n << k < N) {
      k++;
    }
    // Ensure k >= 1 to have at least one step option
    if (k < 1n) k = 1n;
    const kInt = Number(k);

    // Random step sizes r_i, Sage: randrange(1, N)
    const M: Map<number, [bigint, T]> = new Map();
    const maxR = N > 1n ? N - 1n : 1n;
    for (let i = 0; i < kInt; i++) {
      const r = randstate.randint(1n, maxR);
      const e = power(base, r);
      M.set(i, [r, e]);
    }

    // First random walk: "tame" kangaroo
    // Starts at base^ub, takes N steps
    let H = power(base, ub);
    let c = ub;
    for (let i = 0n; i < N; i++) {
      const hashIdx = Number(hash(H) % k);
      const entry = M.get(hashIdx);
      const [r, e] = entry ?? [1n, power(base, 1n)];
      H = multiply(H, e);
      c = c + r;
    }

    // H is now the "trap": H = base^c where c = ub + (sum of steps)
    const mem = H;

    // Second random walk: "wild" kangaroo
    // Starts at a, tries to reach the same position as tame kangaroo
    H = a;
    let d = 0n;

    // Walk until we either find the trap or overshoot
    while (c - d >= lb) {
      // Check if we've found the trap
      if (ub >= c - d && elementsEqual(H, mem)) {
        // H = a * base^d = base^(c) = mem
        // So a = base^(c - d)
        return c - d;
      }

      const hashIdx = Number(hash(H) % k);
      const entry = M.get(hashIdx);
      const [r, e] = entry ?? [1n, power(base, 1n)];
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
 * const p = 1019n; // prime
 * const base = Mod(2n, p);
 * const target = base.pow(500n);
 * const x = discrete_log_rho(target, base, p - 1n, '*'); // p-1 must be prime!
 *
 * // For prime p where p-1 = 2q with q prime, use the order-q subgroup:
 * // const generator = base.pow(2n); // has order q
 * // const x = discrete_log_rho(target, generator, q, '*');
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
  assertNoCustomOps(operation, identity, inverse, op);
  const ordBig = toBigInt(ord);

  if (ordBig <= 0n) {
    throw new ValueError('discrete_log_rho() requires positive order');
  }

  // Check that order is prime
  if (!is_prime(ordBig)) {
    throw new ValueError('for Pollard rho algorithm the order of the group must be prime');
  }

  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let power: (x: T, n: bigint) => T;
  let multiply: (x: T, y: T) => T;
  let isId: (x: T) => boolean;

  if (isMult) {
    power = (x: T, n: bigint) =>
      (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    multiply = (x: T, y: T) =>
      (x as unknown as MultiplicativeGroupElement).mul(
        y as unknown as MultiplicativeGroupElement
      ) as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    multiply = (x: T, y: T) =>
      (x as unknown as AdditiveGroupElement).add(
        y as unknown as AdditiveGroupElement
      ) as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError(
        "identity, inverse and operation must all be specified when operation is 'other'"
      );
    }
    power = (x: T, n: bigint) => multiple(x, n, operation, identity, inverse, op);
    multiply = op;
    isId = (x: T) => x.eq(identity);
  }

  // Handle identity case
  if (isId(a)) {
    return 0n;
  }

  // Handle base being identity
  if (isId(base)) {
    throw new ValueError('base must not be the identity');
  }

  // Verify a is in the group generated by base
  // (a^ord should be identity)
  if (!isId(power(a, ordBig))) {
    throw new ValueError('target element is not in the group generated by base');
  }

  const isqrtord = isqrt(ordBig);

  // SageMath parameters
  const partitionSize = 20;
  const memorySize = 4;

  // Fall back to BSGS for small orders
  if (isqrtord < BigInt(partitionSize)) {
    return bsgs(base, a, [0n, ordBig - 1n], operation, identity, inverse, op);
  }

  // Reset bound: 8 * sqrt(ord) iterations per attempt
  const resetBound = 8n * isqrtord;

  // Modular arithmetic on exponents
  const I = (x: bigint) => Mod(x, ordBig);

  // Get random state
  const randstate = current_randstate();

  // Outer loop for multiple attempts (to avoid infinite loops)
  for (let s = 0; s < 10; s++) {
    // Setup random walk function
    // m[i], n[i] are random exponents for partition i
    // M[i] = base^m[i] * a^n[i] is the multiplier for partition i
    const m: bigint[] = [];
    const n: bigint[] = [];
    const M: T[] = [];

    for (let i = 0; i < partitionSize; i++) {
      const mi = randstate.randint(0n, ordBig - 1n);
      const ni = randstate.randint(0n, ordBig - 1n);
      m.push(mi);
      n.push(ni);
      M.push(multiply(power(base, mi), power(a, ni)));
    }

    // Initial random point: x = base^ax
    let ax = I(randstate.randint(0n, ordBig - 1n));
    let x = power(base, ax.value);
    let bx = I(0n);

    // Sigma array tracks when we stored values (iteration, element key)
    const sigma: Array<[number, string | null]> = [];
    for (let i = 0; i < memorySize; i++) {
      sigma.push([0, null]);
    }

    // H is the hash table storing (ax, bx) for each element seen
    const H = new Map<string, [bigint, bigint]>();

    let i0 = 0;
    let nextsigma = 0;

    for (let i = 0; i < Number(resetBound); i++) {
      // Random walk step
      const hashIdx = rhoHash(x) % partitionSize;
      x = multiply(M[hashIdx]!, x);
      ax = ax.add(I(m[hashIdx]!));
      bx = bx.add(I(n[hashIdx]!));

      // Get string key for element
      const xKey = elementToString(x);

      // Look for collisions
      const stored = H.get(xKey);
      if (stored !== undefined) {
        const [ay, by] = stored;
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
        const oldKey = sigma[i0]![1];
        if (oldKey !== null) {
          H.delete(oldKey);
        }

        // Store current value
        sigma[i0] = [i, xKey];
        i0 = (i0 + 1) % memorySize;

        // Next storage time: 3 * oldest stored iteration
        // This spreads out storage points geometrically
        nextsigma = 3 * sigma[i0]![0];

        H.set(xKey, [ax.value, bx.value]);
      }
    }
  }

  throw new ValueError('Pollard Rho failed to find discrete log');
}

/**
 * Alias for discrete_log for backward compatibility.
 */
export const discrete_log_generic = discrete_log;
