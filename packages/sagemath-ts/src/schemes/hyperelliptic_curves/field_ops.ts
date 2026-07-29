/**
 * @module sage/schemes/hyperelliptic_curves/field_ops
 * @description Base-ring helper operations used by the hyperelliptic curve port.
 *
 * This file has **no upstream counterpart**.  In SageMath every base ring is a
 * `Parent` and every element a `RingElement`, so `K.characteristic()`,
 * `a.is_square()`, `a.sqrt(all=True)`, `iter(K)` and `sorted(...)` are uniformly
 * available.  In this port the base rings are plain TypeScript classes
 * (`FiniteFieldPrime`, `FiniteFieldExtension`, `RationalField`, ...) that do not
 * share a single interface, so the operations Sage takes for granted are
 * collected here and dispatched structurally.
 *
 * Everything in this file is an implementation detail of
 * `schemes/hyperelliptic_curves/`; the semantics mirror Sage exactly:
 *
 * - `characteristic_of`  -> `K.characteristic()`
 * - `cardinality_of`     -> `K.cardinality()` (`null` when `K` is infinite)
 * - `degree_of`          -> `K.degree()` (absolute degree over the prime field)
 * - `iterate_field`      -> `iter(K)`
 * - `is_square_of`       -> `a.is_square()`
 * - `sqrt_all_of`        -> `a.sqrt(all=True, extend=False)`
 * - `compare_elements`   -> the ordering used by Python's `sorted`
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import type { CoefficientRing, RingElement } from '../../rings/polynomial/polynomial_element.js';

/**
 * A base ring for a hyperelliptic curve.
 *
 * Structurally this is Sage's `Parent`: the operations actually used are looked
 * up at run time, so `FiniteFieldPrime`, `FiniteFieldExtension`, `PrimeField`
 * and `RationalField` all qualify without a shared base class.
 */
export interface HyperellipticBaseRing<C extends RingElement> extends CoefficientRing<C> {
  toString(): string;
}

/** Structural view of the members we probe for on a base ring. */
interface RingProbe {
  characteristic?: bigint | (() => bigint | number);
  cardinality?: () => bigint | number;
  order?: bigint | (() => bigint | number | string);
  degree?: number | (() => bigint | number);
  elements?: () => Iterable<unknown>;
  list?: () => unknown[];
  is_field?: () => boolean;
  [Symbol.iterator]?: () => Iterator<unknown>;
}

/** Structural view of the members we probe for on a ring element. */
interface ElementProbe {
  is_square?: () => boolean;
  sqrt?: (opts?: unknown) => unknown;
  div?: (other: unknown) => unknown;
  inv?: () => unknown;
  pow?: (n: bigint | number) => unknown;
  toBigInt?: () => bigint;
  lift?: () => bigint;
  integer_representation?: () => bigint;
  cmp?: (other: unknown) => number;
  value?: bigint;
}

/**
 * Return the characteristic of `K`.
 *
 * Mirrors `K.characteristic()`.  `FiniteFieldPrime` and `FiniteFieldExtension`
 * expose it as a data property, `RationalField` as a method.
 */
export function characteristic_of<C extends RingElement>(K: HyperellipticBaseRing<C>): bigint {
  const probe = K as unknown as RingProbe;
  const c = probe.characteristic;
  if (typeof c === 'bigint') {
    return c;
  }
  if (typeof c === 'function') {
    return BigInt((c as () => bigint | number).call(K));
  }
  throw new NotImplementedError(
    `SAGE_NOT_IMPLEMENTED: characteristic of base ring ${K} is not available`
  );
}

/**
 * Return `K.cardinality()`, or `null` when `K` is infinite.
 */
export function cardinality_of<C extends RingElement>(K: HyperellipticBaseRing<C>): bigint | null {
  const probe = K as unknown as RingProbe;
  if (typeof probe.cardinality === 'function') {
    return BigInt(probe.cardinality.call(K));
  }
  if (typeof probe.order === 'bigint') {
    return probe.order;
  }
  if (characteristic_of(K) === 0n) {
    return null;
  }
  throw new NotImplementedError(
    `SAGE_NOT_IMPLEMENTED: cardinality of base ring ${K} is not available`
  );
}

/**
 * Return `K.degree()`: the degree of `K` over its prime field.
 *
 * Prime fields (and `QQ`) have degree 1.
 */
export function degree_of<C extends RingElement>(K: HyperellipticBaseRing<C>): number {
  const probe = K as unknown as RingProbe;
  if (typeof probe.degree === 'number') {
    return probe.degree;
  }
  if (typeof probe.degree === 'function') {
    return Number((probe.degree as () => bigint | number).call(K));
  }
  return 1;
}

/** Whether `K` is a finite field. */
export function is_finite_field<C extends RingElement>(K: HyperellipticBaseRing<C>): boolean {
  return characteristic_of(K) !== 0n;
}

/**
 * Iterate over all elements of a finite `K`, mirroring `for x in K`.
 */
export function iterate_field<C extends RingElement>(K: HyperellipticBaseRing<C>): Iterable<C> {
  const probe = K as unknown as RingProbe;
  if (typeof probe[Symbol.iterator] === 'function') {
    return K as unknown as Iterable<C>;
  }
  if (typeof probe.elements === 'function') {
    return probe.elements.call(K) as Iterable<C>;
  }
  if (typeof probe.list === 'function') {
    return probe.list.call(K) as C[];
  }
  throw new NotImplementedError(`SAGE_NOT_IMPLEMENTED: base ring ${K} is not iterable`);
}

/** `a / b` in the base ring. */
export function div_elements<C extends RingElement>(a: C, b: C): C {
  const probe = a as unknown as ElementProbe;
  if (typeof probe.div === 'function') {
    return probe.div.call(a, b) as C;
  }
  const bp = b as unknown as ElementProbe;
  if (typeof bp.inv === 'function') {
    return a.mul(bp.inv.call(b) as C) as C;
  }
  throw new ValueError('base ring does not support division');
}

/** `a ** n` in the base ring (n >= 0). */
export function pow_element<C extends RingElement>(a: C, n: bigint): C {
  const probe = a as unknown as ElementProbe;
  if (typeof probe.pow === 'function') {
    return probe.pow.call(a, n) as C;
  }
  // Square-and-multiply fallback.
  if (n < 0n) {
    throw new ValueError('negative exponent requires an inv() method');
  }
  let result: C | null = null;
  let base = a;
  let e = n;
  while (e > 0n) {
    if (e & 1n) {
      result = result === null ? base : (result.mul(base) as C);
    }
    base = base.mul(base) as C;
    e >>= 1n;
  }
  if (result === null) {
    throw new ValueError('pow_element: cannot build one() without a pow method');
  }
  return result;
}

/**
 * Mirror `a.is_square()`.
 *
 * Over a finite field of odd order this is the Euler criterion; in
 * characteristic 2 the Frobenius is a bijection so every element is a square.
 */
export function is_square_of<C extends RingElement>(K: HyperellipticBaseRing<C>, a: C): boolean {
  const probe = a as unknown as ElementProbe;
  if (typeof probe.is_square === 'function') {
    return probe.is_square.call(a);
  }
  const q = cardinality_of(K);
  if (q === null) {
    throw new NotImplementedError(
      `SAGE_NOT_IMPLEMENTED: is_square over the infinite base ring ${K}`
    );
  }
  if (a.isZero()) {
    return true;
  }
  if (q % 2n === 0n) {
    return true;
  }
  return pow_element(a, (q - 1n) / 2n).eq(K.one() as C);
}

/**
 * Mirror `a.sqrt(all=True, extend=False)`: all square roots of `a` in `K`,
 * sorted with {@link compare_elements} so that the result is deterministic
 * exactly as Sage's `ys.sort()` makes `lift_x` deterministic.
 */
export function sqrt_all_of<C extends RingElement>(K: HyperellipticBaseRing<C>, a: C): C[] {
  const roots = sqrt_all_unsorted(K, a);
  roots.sort((u, v) => compare_elements(u, v));
  return roots;
}

function sqrt_all_unsorted<C extends RingElement>(K: HyperellipticBaseRing<C>, a: C): C[] {
  const q = cardinality_of(K);

  if (q === null) {
    // Infinite ring (QQ, ...): use the element's own sqrt.
    const probe = a as unknown as ElementProbe;
    if (typeof probe.sqrt === 'function') {
      if (!(typeof probe.is_square === 'function' && probe.is_square.call(a))) {
        return [];
      }
      const res = probe.sqrt.call(a, { all: true, extend: false });
      return Array.isArray(res) ? (res as C[]) : [res as C];
    }
    throw new NotImplementedError(`SAGE_NOT_IMPLEMENTED: square roots over the base ring ${K}`);
  }

  if (a.isZero()) {
    return [K.zero() as C];
  }

  if (q % 2n === 0n) {
    // Characteristic 2: x -> x^2 is a bijection, the unique square root of a is
    // a^(q/2).
    return [pow_element(a, q / 2n)];
  }

  if (!is_square_of(K, a)) {
    return [];
  }

  const probe = a as unknown as ElementProbe;
  const r =
    typeof probe.sqrt === 'function' ? (probe.sqrt.call(a) as C) : finite_field_sqrt(K, a, q);
  const minusR = r.neg() as C;
  if (r.eq(minusR)) {
    return [r];
  }
  return [r, minusR];
}

/** Cache of a quadratic non-residue per base ring, keyed by identity. */
const nonResidueCache = new WeakMap<object, RingElement>();

/**
 * Tonelli–Shanks over an arbitrary finite field of odd order.
 *
 * Used only when the element class does not provide its own `sqrt()`.
 */
function finite_field_sqrt<C extends RingElement>(K: HyperellipticBaseRing<C>, a: C, q: bigint): C {
  const one = K.one() as C;

  // q - 1 = 2^s * t with t odd
  let t = q - 1n;
  let s = 0n;
  while (t % 2n === 0n) {
    t /= 2n;
    s += 1n;
  }

  if (s === 1n) {
    // q = 3 mod 4
    return pow_element(a, (q + 1n) / 4n);
  }

  // Find a quadratic non-residue (cached per base ring: the search is O(1)
  // expected but scans the field iterator, which is expensive to redo).
  let z = nonResidueCache.get(K as object) as C | undefined;
  if (z === undefined) {
    for (const c of iterate_field(K)) {
      if (c.isZero()) {
        continue;
      }
      if (!pow_element(c, (q - 1n) / 2n).eq(one)) {
        z = c;
        break;
      }
    }
    if (z === undefined) {
      throw new ValueError('no quadratic non-residue found; base ring is not a field');
    }
    nonResidueCache.set(K as object, z);
  }

  let m = s;
  let c = pow_element(z, t);
  let tt = pow_element(a, t);
  let r = pow_element(a, (t + 1n) / 2n);

  while (!tt.eq(one)) {
    let i = 0n;
    let temp = tt;
    while (!temp.eq(one)) {
      temp = temp.mul(temp) as C;
      i += 1n;
      if (i === m) {
        throw new ValueError('element is not a square');
      }
    }
    let b = c;
    for (let j = 0n; j < m - i - 1n; j++) {
      b = b.mul(b) as C;
    }
    m = i;
    c = b.mul(b) as C;
    tt = tt.mul(c) as C;
    r = r.mul(b) as C;
  }

  return r;
}

/**
 * Total order on base-ring elements matching Python's `sorted`.
 *
 * - Finite fields: by the integer representation (for `GF(p)` the residue in
 *   `[0, p)`, for `GF(p^n)` the base-`p` digits of the coefficient vector),
 *   which is the order Sage's `sorted` uses on `FiniteField` elements.
 * - `QQ`: numeric order.
 */
export function compare_elements<C extends RingElement>(a: C, b: C): number {
  const pa = a as unknown as ElementProbe;
  const pb = b as unknown as ElementProbe;

  if (typeof pa.cmp === 'function') {
    return pa.cmp.call(a, b);
  }

  const ia = element_to_bigint(a);
  const ib = element_to_bigint(b);
  if (ia !== null && ib !== null) {
    return ia < ib ? -1 : ia > ib ? 1 : 0;
  }

  void pb;
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: ordering of base ring elements without cmp/integer representation'
  );
}

/**
 * Order a list of roots the way Sage's `Polynomial.roots()` does.
 *
 * `roots()` reads the linear factors off `self.factor()`, and
 * `Factorization.sort` (`sage/structure/factorization.py:671-741`) orders them
 * by `(degree, exponent, prime)`.  For the simple roots this module deals with
 * all the linear factors have degree 1 and exponent 1, so the order is the
 * order of the polynomials `x - r`, i.e. the order of the constant terms `-r`.
 */
export function sort_roots_like_sage<C extends RingElement>(roots: C[]): C[] {
  return roots.slice().sort((u, v) => compare_elements(u.neg() as C, v.neg() as C));
}

/** The canonical non-negative integer attached to a finite-field element. */
export function element_to_bigint<C extends RingElement>(a: C): bigint | null {
  const probe = a as unknown as ElementProbe;
  if (typeof probe.integer_representation === 'function') {
    return probe.integer_representation.call(a);
  }
  if (typeof probe.toBigInt === 'function') {
    return probe.toBigInt.call(a);
  }
  if (typeof probe.lift === 'function') {
    return probe.lift.call(a);
  }
  if (typeof probe.value === 'bigint') {
    return probe.value;
  }
  return null;
}

/**
 * A field embedding `K -> L` of finite fields, mirroring `Hom(K, L)[0]`
 * (`hyperelliptic_finite_field.py:1311`), and the canonical coercion used by
 * `change_ring` (`hyperelliptic_generic.py:113-143`).
 *
 * For a prime `K` this is the canonical map on residues.  Otherwise
 * `K = GF(p)[a]/(m(a))` and any root `r` of `m` in `L` induces an embedding
 * `sum c_i a^i -> sum c_i r^i`; Sage picks the first element of `Hom(K, L)`,
 * this port picks the first root found while iterating `L`.  The number of
 * points of a curve over `L` does not depend on the choice, since two
 * embeddings differ by an automorphism of `L`.
 */
export function field_embedding<C extends RingElement>(
  K: HyperellipticBaseRing<C>,
  L: HyperellipticBaseRing<RingElement>
): (c: C) => RingElement {
  if (K === (L as unknown as HyperellipticBaseRing<C>)) {
    return (cc: C): RingElement => cc;
  }

  if (degree_of(K) === 1) {
    // A prime base ring has a unique ring map into `L`, so whatever coercion
    // `L` provides is the right one.  (For a base ring of higher degree it is
    // *not*: re-reading the coefficient vector of an element of `GF(p^m)` in
    // `GF(p^(mn))` is not a field homomorphism, so that path is never taken.)
    let directWorks = false;
    try {
      directWorks = L.__call__(K.one()).eq(L.one() as never) && L.__call__(K.zero()).isZero();
    } catch {
      directWorks = false;
    }
    if (directWorks) {
      return (cc: C): RingElement => L.__call__(cc);
    }
    return (cc: C): RingElement => {
      const v = element_to_bigint(cc);
      if (v === null) {
        throw new NotImplementedError(
          `SAGE_NOT_IMPLEMENTED: embedding of the base ring ${K} into ${L}`
        );
      }
      return L.__call__(v);
    };
  }

  const probeRing = K as unknown as {
    modulus?: { degree(): number; getCoeff(i: number): { toBigInt?: () => bigint } };
  };
  const modulus = probeRing.modulus;
  if (modulus === undefined) {
    throw new NotImplementedError(
      `SAGE_NOT_IMPLEMENTED: embedding of the base field ${K} into an extension`
    );
  }

  const md = modulus.degree();
  const modCoeffs: RingElement[] = [];
  for (let i = 0; i <= md; i++) {
    const c = modulus.getCoeff(i);
    const v = c.toBigInt?.();
    if (v === undefined) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: lifting the modulus of the base field to its extension'
      );
    }
    modCoeffs.push(L.__call__(v));
  }

  // Find a root of the modulus in L by Horner evaluation.
  let root: RingElement | null = null;
  for (const y of iterate_field(L)) {
    let acc = modCoeffs[md]!;
    for (let i = md - 1; i >= 0; i--) {
      acc = acc.mul(y).add(modCoeffs[i]!);
    }
    if (acc.isZero()) {
      root = y;
      break;
    }
  }
  if (root === null) {
    throw new ValueError(
      `the modulus of ${K} has no root in the extension; not a valid field extension`
    );
  }

  const probeElement = (c: C): bigint[] => {
    const el = c as unknown as { coefficients?: () => { toBigInt(): bigint }[] };
    if (typeof el.coefficients !== 'function') {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: coefficient vector of an element of the base field'
      );
    }
    return el.coefficients().map((x) => x.toBigInt());
  };

  return (cc: C): RingElement => {
    const coeffs = probeElement(cc);
    let acc = L.zero();
    let pw = L.one();
    for (const v of coeffs) {
      acc = acc.add(pw.mul(L.__call__(v)));
      pw = pw.mul(root);
    }
    return acc;
  };
}

/**
 * Whether the absolute trace `Tr_{K/F_p}(a)` vanishes.
 *
 * Mirrors `a.trace() == 0` for a finite field element; for a prime field the
 * absolute trace is the element itself.
 */
export function absolute_trace_is_zero<C extends RingElement>(
  K: HyperellipticBaseRing<C>,
  a: C
): boolean {
  const probe = a as unknown as { trace?: () => { isZero(): boolean } };
  if (typeof probe.trace === 'function') {
    return probe.trace.call(a).isZero();
  }
  if (degree_of(K) !== 1) {
    throw new NotImplementedError(`SAGE_NOT_IMPLEMENTED: absolute trace over the base ring ${K}`);
  }
  return a.isZero();
}

/**
 * `n * a` in the base ring, for a plain integer `n` (Sage coerces `n` into the
 * ring and multiplies).
 */
export function int_times<C extends RingElement>(K: HyperellipticBaseRing<C>, n: bigint, a: C): C {
  return (K.__call__(n) as C).mul(a) as C;
}
