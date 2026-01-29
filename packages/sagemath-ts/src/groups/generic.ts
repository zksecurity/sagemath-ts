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

import { CRT_list, factor, is_prime, isqrt, type Factorization } from '../arith/misc.js';
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
  if ((isMultiplicative(operation) || isAdditive(operation)) && (identity !== undefined || inverse !== undefined || op !== undefined)) {
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

function deriveStandardGroupOps<T extends GroupElement>(sample: T, operation: OperationType): GroupOps<T> {
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);
  if (!isMult && !isAdd) {
    throw new ValueError("identity, inverse and operation must all be specified when operation is 'other'");
  }

  const parent = (sample as unknown as { parent?: { one?: () => unknown; zero?: () => unknown } }).parent;
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
    ? ((x: T, y: T) => {
        if ((x as unknown as MultiplicativeGroupElement).mul) {
          return (x as unknown as MultiplicativeGroupElement).mul(y as unknown as MultiplicativeGroupElement) as unknown as T;
        }
        if (typeof x === 'bigint') {
          return ((x as unknown as bigint) * (y as unknown as bigint)) as unknown as T;
        }
        if (typeof x === 'number') {
          return ((x as unknown as number) * (y as unknown as number)) as unknown as T;
        }
        throw new ValueError('Cannot multiply elements');
      })
    : ((x: T, y: T) => {
        if ((x as unknown as AdditiveGroupElement).add) {
          return (x as unknown as AdditiveGroupElement).add(y as unknown as AdditiveGroupElement) as unknown as T;
        }
        if (typeof x === 'bigint') {
          return ((x as unknown as bigint) + (y as unknown as bigint)) as unknown as T;
        }
        if (typeof x === 'number') {
          return ((x as unknown as number) + (y as unknown as number)) as unknown as T;
        }
        throw new ValueError('Cannot add elements');
      });

  const inverse = isMult
    ? ((x: T) => {
        if ((x as unknown as MultiplicativeGroupElement).inv) {
          return (x as unknown as MultiplicativeGroupElement).inv() as unknown as T;
        }
        if (typeof x === 'number') {
          return (1 / (x as unknown as number)) as unknown as T;
        }
        throw new ValueError('Cannot invert element');
      })
    : ((x: T) => {
        if ((x as unknown as AdditiveGroupElement).neg) {
          return (x as unknown as AdditiveGroupElement).neg() as unknown as T;
        }
        if (typeof x === 'bigint') {
          return (-(x as unknown as bigint)) as unknown as T;
        }
        if (typeof x === 'number') {
          return (-(x as unknown as number)) as unknown as T;
        }
        throw new ValueError('Cannot negate element');
      });

  const isIdentity = isMult
    ? ((x: T) => {
        if ((x as unknown as MultiplicativeGroupElement).isOne) {
          return (x as unknown as MultiplicativeGroupElement).isOne();
        }
        return elementsEqual(x, identity as T);
      })
    : ((x: T) => {
        if ((x as unknown as AdditiveGroupElement).isZero) {
          return (x as unknown as AdditiveGroupElement).isZero();
        }
        return elementsEqual(x, identity as T);
      });

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
    throw new ValueError('bsgs() requires 0 <= lb <= ub');
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
    power = (x: T, n: bigint) => (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    multiply = (x: T, y: T) => (x as unknown as MultiplicativeGroupElement).mul(y as unknown as MultiplicativeGroupElement) as unknown as T;
    invert = (x: T) => (x as unknown as MultiplicativeGroupElement).inv() as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    multiply = (x: T, y: T) => (x as unknown as AdditiveGroupElement).add(y as unknown as AdditiveGroupElement) as unknown as T;
    invert = (x: T) => (x as unknown as AdditiveGroupElement).neg() as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError("identity, inverse and operation must all be specified when operation is 'other'");
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
  let c = multiply(bInv, aLb);

  // For small ranges, use simple linear search
  if (range < 30n) {
    let d = c;
    for (let i = lb; i <= ub; i++) {
      if (isId(d)) {
        return i;
      }
      d = multiply(d, a);
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

  throw new ValueError(`log of ${b} to the base ${a} does not exist in [${lb}, ${ub}]`);
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
    power = (x: T, n: bigint) => (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    multiply = (x: T, y: T) => (x as unknown as MultiplicativeGroupElement).mul(y as unknown as MultiplicativeGroupElement) as unknown as T;
    invert = (x: T) => (x as unknown as MultiplicativeGroupElement).inv() as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    multiply = (x: T, y: T) => (x as unknown as AdditiveGroupElement).add(y as unknown as AdditiveGroupElement) as unknown as T;
    invert = (x: T) => (x as unknown as AdditiveGroupElement).neg() as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError("identity, inverse and operation must all be specified when operation is 'other'");
    }
    power = (x: T, n: bigint) => multiple(x, n, operation, identity, inverse, op);
    multiply = op;
    invert = inverse;
    isId = (x: T) => x.eq(identity);
  }

  // Get factorization of n
  const factors = factorization ?? factor(nBig);

  // Handle trivial cases
  if (isId(a)) {
    return 0n;
  }

  // Filter out sign factor
  const primeFactors = factors.filter(([p]) => p > 0n);

  if (primeFactors.length === 0) {
    // n = 1, so a must be identity
    if (isId(a)) {
      return 0n;
    }
    throw new ValueError(`no discrete log: ${a} is not the identity but order is 1`);
  }

  // Solve DLP in each prime power subgroup
  const residues: bigint[] = [];
  const moduli: bigint[] = [];

  for (const [p, e] of primeFactors) {
    const pe = p ** e;
    const cofactor = nBig / pe;

    // gamma = base^cofactor has order dividing p^e
    const gamma = power(base, cofactor);

    // h = a^cofactor - target in the p^e subgroup
    const h = power(a, cofactor);

    // Now solve gamma^x = h for x in [0, p^e)
    // Use the standard Pohlig-Hellman lifting for prime powers

    // gamma_0 = gamma^(p^(e-1)) has order p
    const gamma_0 = power(gamma, p ** (e - 1n));

    let x_i = 0n;
    let h_k = h; // h_k = h * gamma^(-x_k) where x_k = sum of digits found so far

    for (let k = 0n; k < e; k++) {
      // Compute h_k^(p^(e-1-k))
      const target = power(h_k, p ** (e - 1n - k));

      // Solve gamma_0^d_k = target where d_k is the k-th digit in base p
      let d_k: bigint;
      try {
        d_k = bsgs(gamma_0, target, [0n, p - 1n], operation, identity, inverse, op);
      } catch {
        d_k = 0n;
      }

      x_i += d_k * p ** k;

      // Update h_k = h_k * gamma^(-d_k * p^k) = h * gamma^(-x_i)
      if (k < e - 1n) {
        const adjustment = power(gamma, d_k * p ** k);
        h_k = multiply(h_k, invert(adjustment));
      }
    }

    residues.push(x_i);
    moduli.push(pe);
  }

  // Combine using CRT
  return CRT_list(residues, moduli) % nBig;
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
 * @param ord - Order of base (computed if not provided)
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
  let isId: (x: T) => boolean;

  if (isMult) {
    power = (x: T, n: bigint) => (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError("identity, inverse and operation must all be specified when operation is 'other'");
    }
    power = (x: T, n: bigint) => multiple(x, n, operation, identity, inverse, op);
    isId = (x: T) => x.eq(identity);
  }

  // Convert ord to bigint if provided
  let ordBig: bigint | undefined;
  if (ord !== undefined) {
    ordBig = toBigInt(ord);
  }

  // If order not provided, try to get it from the element
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

  // Handle identity case
  if (isId(a)) {
    return 0n;
  }

  // Handle base being identity
  if (isId(base)) {
    if (!isId(a)) {
      throw new ValueError(`no discrete log of ${a} found to base ${base}`);
    }
    return 0n;
  }

  // Use Pohlig-Hellman for composite orders, BSGS for prime orders
  const factors = factor(ordBig);

  // Filter out -1 factor
  const primeFactors = factors.filter(([p]) => p > 0n);

  // If order has only one prime factor (possibly with multiplicity), use BSGS directly
  if (primeFactors.length === 1) {
    return bsgs(base, a, [0n, ordBig - 1n], operation, identity, inverse, op);
  }

  // Use Pohlig-Hellman
  const result = pohlig_hellman(a, base, ordBig, factors, operation, identity, inverse, op);

  // Verify result
  if (!power(base, result).eq(a)) {
    throw new ValueError(`no discrete log of ${a} found to base ${base}`);
  }

  return result;
}

/**
 * Find the order of a group element given a multiple of its order.
 *
 * @param a - Group element
 * @param orderMultiple - A known multiple of the order (i.e., a^multiple = identity)
 * @param factorization - Factorization of multiple (computed if not provided)
 * @param operation - Type of group operation
 * @param identity - Identity element (for custom operations)
 * @param inverse - Inverse function (for custom operations)
 * @param op - Binary operation (for custom operations)
 * @returns The exact order of a
 *
 * @example
 * ```typescript
 * const a = Mod(5n, 37n);
 * // We know 5^36 = 1 mod 37 (by Fermat's little theorem)
 * const order = order_from_multiple(a, 36n, undefined, '*');
 * console.log(order); // The multiplicative order of 5 mod 37
 * ```
 */
export function order_from_multiple<T extends GroupElement>(
  a: T,
  orderMultiple: IntegerLike,
  factorization?: Factorization,
  operation: OperationType = '*',
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): bigint {
  const orderMultipleBig = toBigInt(orderMultiple);
  assertNoCustomOps(operation, identity, inverse, op);
  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let power: (x: T, n: bigint) => T;
  let isId: (x: T) => boolean;

  if (isMult) {
    power = (x: T, n: bigint) => (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError("identity, inverse and operation must all be specified when operation is 'other'");
    }
    power = (x: T, n: bigint) => multiple(x, n, operation, identity, inverse, op);
    isId = (x: T) => x.eq(identity);
  }

  // Identity has order 1
  if (isId(a)) {
    return 1n;
  }

  // Get factorization
  const factors = factorization ?? factor(orderMultipleBig);

  // The order divides the multiple, so we try removing prime factors
  let order = orderMultipleBig;

  for (const [p, e] of factors) {
    if (p === -1n) continue; // Skip sign factor

    // Try to divide out as many copies of p as possible
    for (let i = 0n; i < e; i++) {
      const testOrder = order / p;
      const powered = power(a, testOrder);

      if (isId(powered)) {
        order = testOrder;
      } else {
        break;
      }
    }
  }

  return order;
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
  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let multiply: (x: T, y: T) => T;
  let isId: (x: T) => boolean;

  if (isMult) {
    multiply = (x: T, y: T) => (x as unknown as MultiplicativeGroupElement).mul(y as unknown as MultiplicativeGroupElement) as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
  } else if (isAdd) {
    multiply = (x: T, y: T) => (x as unknown as AdditiveGroupElement).add(y as unknown as AdditiveGroupElement) as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError("identity, inverse and operation must all be specified when operation is 'other'");
    }
    multiply = op;
    isId = (x: T) => x.eq(identity);
  }

  // Identity has order 1
  if (isId(a)) {
    return 1n;
  }

  // We look for a^m = identity by computing powers of a
  let current = a;
  let m = 1n;

  while (m < maxIterationsBig) {
    current = multiply(current, a);
    m++;

    if (isId(current)) {
      return m;
    }
  }

  throw new ValueError(`could not find multiple of order within ${maxIterationsBig} iterations`);
}

/**
 * Check if an element has exactly the given order.
 *
 * @param a - Group element
 * @param n - Proposed order
 * @param operation - Type of group operation
 * @returns true if the order of a is exactly n
 *
 * @example
 * ```typescript
 * const a = Mod(2n, 7n);
 * has_order(a, 3n, '*'); // true (2^3 = 8 = 1 mod 7)
 * has_order(a, 6n, '*'); // false
 * ```
 */
export function has_order<T extends GroupElement>(
  a: T,
  n: IntegerLike,
  operation: OperationType = '*'
): boolean {
  const nBig = toBigInt(n);
  if (nBig <= 0n) {
    return false;
  }

  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let power: (x: T, n: bigint) => T;
  let isId: (x: T) => boolean;

  if (isMult) {
    power = (x: T, n: bigint) => (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else {
    // Custom operation not supported without full parameters
    return false;
  }

  // a^n must be identity
  const powered = power(a, nBig);
  if (!isId(powered)) {
    return false;
  }

  // For all prime factors p of n, a^(n/p) must NOT be identity
  const factors = factor(nBig);

  for (const [p] of factors) {
    if (p === -1n) continue;

    const testPower = power(a, nBig / p);
    if (isId(testPower)) {
      return false;
    }
  }

  return true;
}

/**
 * Iterator for multiples/powers of a group element.
 *
 * For additive groups: yields P0, P0+P, P0+2P, ..., P0+(n-1)P
 * For multiplicative groups: yields P0, P0*P, P0*P^2, ..., P0*P^(n-1)
 *
 * @param P - Step element
 * @param n - Number of multiples to generate
 * @param P0 - Starting element (default: identity)
 * @param operation - Type of group operation
 * @param identity - Identity element (for custom operations)
 * @param inverse - Inverse function (for custom operations)
 * @param op - Binary operation (for custom operations)
 * @returns Generator yielding [index, element] pairs
 */
export function* multiples<T extends GroupElement>(
  P: T,
  n: IntegerLike,
  P0?: T,
  operation: OperationType = '*',
  identity?: T,
  inverse?: (x: T) => T,
  op?: (x: T, y: T) => T
): Generator<[bigint, T], void, undefined> {
  const nBig = toBigInt(n);
  assertNoCustomOps(operation, identity, inverse, op);
  if (nBig < 0n) {
    throw new ValueError('n cannot be negative in multiples');
  }

  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let multiply: (x: T, y: T) => T;
  let getIdentity: () => T;

  if (isMult) {
    multiply = (x: T, y: T) => (x as unknown as MultiplicativeGroupElement).mul(y as unknown as MultiplicativeGroupElement) as unknown as T;
    getIdentity = () => (P as unknown as MultiplicativeGroupElement).pow(0n) as unknown as T;
  } else if (isAdd) {
    multiply = (x: T, y: T) => (x as unknown as AdditiveGroupElement).add(y as unknown as AdditiveGroupElement) as unknown as T;
    getIdentity = () => (P as unknown as AdditiveGroupElement).mul(0n) as unknown as T;
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError("identity, inverse and operation must all be specified when operation is 'other'");
    }
    multiply = op;
    getIdentity = () => identity;
  }

  let current = P0 ?? getIdentity();

  for (let i = 0n; i < nBig; i++) {
    yield [i, current];
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
    power = (x: T, n: bigint) => (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    multiply = (x: T, y: T) => (x as unknown as MultiplicativeGroupElement).mul(y as unknown as MultiplicativeGroupElement) as unknown as T;
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    multiply = (x: T, y: T) => (x as unknown as AdditiveGroupElement).add(y as unknown as AdditiveGroupElement) as unknown as T;
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError("identity, inverse and operation must all be specified when operation is 'other'");
    }
    power = (x: T, n: bigint) => multiple(x, n, operation, identity, inverse, op);
    multiply = op;
  }

  // Default hash function
  const hash = hashFunction ?? ((x: T) => {
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

  // Retry loop to handle random walk failures
  for (let attempt = 0; attempt < 10; attempt++) {
    // Set up random walk function
    // We need k values r_i and corresponding base^r_i
    // k should be chosen such that 2^k >= N
    let k = 0n;
    while ((1n << k) < N) {
      k++;
    }
    // Ensure k >= 1 to have at least one step option
    if (k < 1n) k = 1n;
    const kInt = Number(k);

    // Random step sizes r_i in [1, N)
    // Use deterministic "random" based on attempt number for reproducibility
    const M: Map<number, [bigint, T]> = new Map();
    const maxR = N > 1n ? N - 1n : 1n;
    for (let i = 0; i < kInt; i++) {
      // Pseudo-random step size based on i and attempt
      const r = 1n + ((BigInt(i + 1) * 17n + BigInt(attempt) * 23n + 7n) % maxR);
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
    throw new ValueError('discrete_log_rho() requires prime order group');
  }

  // Define group operations based on type
  const isMult = isMultiplicative(operation);
  const isAdd = isAdditive(operation);

  let power: (x: T, n: bigint) => T;
  let multiply: (x: T, y: T) => T;
  let isId: (x: T) => boolean;

  if (isMult) {
    power = (x: T, n: bigint) => (x as unknown as MultiplicativeGroupElement).pow(n) as unknown as T;
    multiply = (x: T, y: T) => (x as unknown as MultiplicativeGroupElement).mul(y as unknown as MultiplicativeGroupElement) as unknown as T;
    isId = (x: T) => (x as unknown as MultiplicativeGroupElement).isOne();
  } else if (isAdd) {
    power = (x: T, n: bigint) => (x as unknown as AdditiveGroupElement).mul(n) as unknown as T;
    multiply = (x: T, y: T) => (x as unknown as AdditiveGroupElement).add(y as unknown as AdditiveGroupElement) as unknown as T;
    isId = (x: T) => (x as unknown as AdditiveGroupElement).isZero();
  } else {
    if (identity === undefined || inverse === undefined || op === undefined) {
      throw new ValueError("identity, inverse and operation must all be specified when operation is 'other'");
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
