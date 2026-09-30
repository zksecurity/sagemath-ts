/** Balanced products, powering and order schedules from basemath/bb_group.c. */
import { Z_factor } from './ifactor.js';
import { PariError } from './errors.js';
export function gen_product<T>(values: readonly T[], multiply: (a: T, b: T) => T): T | bigint {
  if (!values.length) return 1n;
  if (values.length === 1) return values[0]!;
  // producttree_scheme splits sizes first, then evaluates each complete level.
  // Depth-first recursion has the same grouping but a different callback order.
  const sizes = producttree_scheme(values.length);
  let products: T[] = [], offset = 0;
  for (const size of sizes) {
    products.push(size === 1 ? values[offset]! : multiply(values[offset]!, values[offset + 1]!));
    offset += size;
  }
  while (products.length > 1) {
    const next: T[] = [];
    for (let i = 0; i < products.length; i += 2)
      next.push(multiply(products[i]!, products[i + 1]!));
    products = next;
  }
  return products[0]!;
}

/** Nonzero unsigned-word exponent; preserves the original operation order. */
export function gen_powu_i<T>(
  x: T,
  n: bigint,
  square: (x: T) => T,
  multiply: (x: T, y: T) => T
): T {
  if (n <= 0n || n >= 1n << 64n)
    throw new RangeError('gen_powu_i requires a nonzero unsigned word');
  if (n < 512n) {
    let y = x;
    for (const bit of n.toString(2).slice(1)) {
      y = square(y);
      if (bit === '1') y = multiply(y, x);
    }
    return y;
  }
  let width = n < 1n << 25n ? 2 : 3;
  const table = [x],
    x2 = square(x);
  for (let i = 1; i < 1 << (width - 1); i++) table.push(multiply(table[i - 1]!, x2));
  let top = n.toString(2).length - 1;
  let result: T | undefined;
  while (top >= 0) {
    width = Math.min(width, top + 1);
    const word = Number((n >> BigInt(top + 1 - width)) & ((1n << BigInt(width)) - 1n));
    let trailing = 0;
    while (!(word & (1 << trailing))) trailing++;
    top -= width;
    const factor = table[word >> (trailing + 1)]!;
    if (result === undefined) result = factor;
    else {
      for (let i = 0; i < width - trailing; i++) result = square(result);
      result = multiply(result, factor);
    }
    for (let i = 0; i < trailing; i++) result = square(result);
    while (top >= 0 && !(n & (1n << BigInt(top)))) {
      result = square(result);
      top--;
    }
  }
  return result!;
}

/** Nonzero integer exponent magnitude, with PARI's arbitrary-precision windows.
 * @see Deviation: PARI polynomial GCD adapters
 */
export function gen_pow_i<T>(x: T, n: bigint, square: (x: T) => T, multiply: (x: T, y: T) => T): T {
  if (n === 0n) throw new RangeError('gen_pow_i requires a nonzero exponent');
  if (n < 0n) n = -n;
  if (n < 1n << 64n) return gen_powu_i(x, n, square, multiply);
  let top = n.toString(2).length - 1;
  const width = top <= 64 ? 3 : top <= 160 ? 4 : top <= 384 ? 5 : top <= 896 ? 6 : 7;
  const table = [x],
    x2 = square(x);
  for (let i = 1; i < 1 << (width - 1); i++) table.push(multiply(table[i - 1]!, x2));
  let result: T | undefined;
  while (top >= 0) {
    const length = Math.min(width, top + 1),
      window = Number((n >> BigInt(top + 1 - length)) & ((1n << BigInt(length)) - 1n));
    let trailing = 0;
    while ((window & (1 << trailing)) === 0) trailing++;
    top -= length;
    const term = table[window >> (trailing + 1)]!;
    if (result === undefined) result = term;
    else {
      for (let i = 0; i < length - trailing; i++) result = square(result);
      result = multiply(result, term);
    }
    for (let i = 0; i < trailing; i++) result = square(result);
    while (top >= 0 && ((n >> BigInt(top)) & 1n) === 0n) {
      result = square(result);
      top--;
    }
  }
  return result!;
}

/** Fused left-to-right powering for a nonzero exponent magnitude.
 * @see Deviation: PARI modular square-root adapters
 */
export function gen_pow_fold<T>(x: T, n: bigint, square: (x: T) => T, multiplySquare: (x: T) => T): T {
  if (n === 0n) throw new RangeError('gen_pow_fold requires a nonzero exponent');
  if (n < 0n) n = -n;
  let result = x;
  for (const bit of n.toString(2).slice(1))
    result = bit === '1' ? multiplySquare(result) : square(result);
  return result;
}

/** Native producttree_scheme (bb_group.c:332), left child gets an odd remainder. */
export function producttree_scheme(n: number): number[] {
  if (n <= 2) return [n];
  let sizes = [n];
  for (let level = 0; level < Math.floor(Math.log2(n - 1)); level++) {
    const next: number[] = [];
    for (const size of sizes) next.push(Math.ceil(size / 2), Math.floor(size / 2));
    sizes = next;
  }
  return sizes;
}


export type GroupOrderFactors = readonly (readonly [bigint, bigint])[];
export type GroupOrder = bigint | GroupOrderFactors | readonly [bigint, GroupOrderFactors];
/** bb_group.c:560–578/668–713. The supplied order must annihilate the element.
 * Factor rows are [prime, exponent]; a vector is [order, factor rows].
 * Native recursion tests identity and powers prime factors without verifying the bound.
 * @see Deviation: PARI generic and extension-curve order adapters
 */
export function gen_order<T>(
  a: T,
  order: GroupOrder | null,
  power: (a: T, n: bigint) => T,
  isIdentity: (a: T) => boolean
): bigint {
  if (order === null)
    throw new PariError(`incorrect type in gen_order [missing order] (${typeof a === 'bigint' ? 't_INT' : 't_VEC'}).`);
  const vector = Array.isArray(order) && typeof order[0] === 'bigint';
  const kind = typeof order === 'bigint' ? 't_INT' : vector ? 't_VEC' : 't_MAT';
  const invalid = () => { throw new PariError(`incorrect type in generic discrete logarithm (order factorization) (${kind}).`); };
  let o: bigint;
  let factors: GroupOrderFactors;
  if (typeof order === 'bigint') {
    if (order <= 0n) return invalid();
    o = order;
    factors = Z_factor(o);
  } else {
    if (vector && order.length !== 2) return invalid();
    factors = vector ? (order[1] as GroupOrderFactors) : order as GroupOrderFactors;
    if (!Array.isArray(factors) || factors.some(row => !Array.isArray(row) || row.length !== 2 ||
      typeof row[0] !== 'bigint' || typeof row[1] !== 'bigint' || row[0] <= 0n || row[1] <= 0n)) return invalid();
    o = vector ? (order[0] as bigint) :
      factors.reduce((v, [p, e]) => v * p ** e, 1n);
    if (o <= 0n) return invalid();
  }
  if (o === 1n) return 1n;
  const rec = (a: T, o: bigint, lo: number, hi: number): bigint => {
    if (isIdentity(a)) return 1n;
    if (lo === hi) {
      const [p, e] = factors[lo]!;
      if (e >= 1n << 63n) throw new PariError('overflow in t_INT-->long assignment.');
      let b = a;
      for (let i = 0n; i < e; i++) {
        if (isIdentity(b)) return p ** i;
        b = power(b, p);
      }
      return p ** e;
    }
    const mid = Math.floor((lo + hi) / 2);
    let cofactor = 1n;
    for (let i = lo; i <= mid; i++) cofactor *= factors[i]![0] ** factors[i]![1];
    const right = rec(power(a, cofactor), o / cofactor, mid + 1, hi);
    const left = rec(power(a, right), o / right, lo, mid);
    return right * left;
  };
  return rec(a, o, 0, factors.length - 1);
}
