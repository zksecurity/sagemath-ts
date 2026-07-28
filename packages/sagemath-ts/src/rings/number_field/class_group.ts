/**
 * @module sage/rings/number_field/class_group
 * @description Class groups of number fields
 *
 * Port of: sage/rings/number_field/class_group.py
 * Reference: reference/sage/src/sage/rings/number_field/class_group.py
 *
 * NOTE: Full class group computation requires PARI's bnfinit.
 * This implementation provides the structure and basic operations.
 * For full PARI compatibility, see DEVIATIONS.md.
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import type { NumberField } from './number_field.js';
import type { NumberFieldIdeal } from './number_field_ideal.js';

/**
 * The class group of a number field.
 *
 * The class group is the quotient of the group of fractional ideals by
 * the subgroup of principal ideals.
 *
 * Class group computation is one of the fundamental problems in
 * algebraic number theory. The structure is:
 *   Cl(K) ≅ Z/d_1 × Z/d_2 × ... × Z/d_r
 * where d_i | d_{i+1} (elementary divisors).
 *
 * The class number h(K) = |Cl(K)| = d_1 * d_2 * ... * d_r.
 *
 * @see Reference: sage/rings/number_field/class_group.py:ClassGroup
 */
export class ClassGroup {
  private readonly _number_field: NumberField;
  private readonly _invariants: readonly bigint[];
  private readonly _gens: readonly ClassGroupElement[];

  constructor(
    number_field: NumberField,
    invariants: bigint[] = [],
    gens: ClassGroupElement[] = []
  ) {
    this._number_field = number_field;
    this._invariants = Object.freeze([...invariants]);
    this._gens = Object.freeze([...gens]);
  }

  /**
   * Return the number field this class group belongs to.
   */
  number_field(): NumberField {
    return this._number_field;
  }

  /**
   * Return the class number (order of the class group).
   *
   * The class number h(K) = |Cl(K)| is a fundamental invariant of the number field.
   * h(K) = 1 means the ring of integers is a unique factorization domain.
   *
   * @see Reference: sage/rings/number_field/class_group.py:order
   */
  order(): bigint {
    if (this._invariants.length === 0) {
      return 1n;
    }
    return this._invariants.reduce((a, b) => a * b, 1n);
  }

  /**
   * Return the class number.
   *
   * Alias for order().
   *
   * @see Reference: sage/rings/number_field/number_field.py:class_number
   */
  class_number(): bigint {
    return this.order();
  }

  /**
   * Return the exponent of the class group.
   *
   * The exponent is the smallest positive integer e such that
   * x^e = 1 for all x in the class group. It equals the largest invariant.
   */
  exponent(): bigint {
    if (this._invariants.length === 0) {
      return 1n;
    }
    // PARI/Sage list the invariants in decreasing order (d_{i+1} | d_i), so the
    // exponent is the largest one, whichever order the caller supplied.
    return this._invariants.reduce((a, b) => (a > b ? a : b), 1n);
  }

  /**
   * Return the invariants (elementary divisors) of the class group.
   *
   * These are the d_i such that Cl(K) ≅ Z/d_1 × Z/d_2 × ... × Z/d_r
   * with d_i | d_{i+1}.
   */
  invariants(): bigint[] {
    return [...this._invariants];
  }

  /**
   * Return the generators of the class group.
   *
   * These are ideal classes that generate the class group.
   * Each generator has order equal to the corresponding invariant.
   *
   * @see Reference: sage/rings/number_field/class_group.py:gens
   */
  gens(): ClassGroupElement[] {
    return [...this._gens];
  }

  /**
   * Return the generators of the class group.
   *
   * Alias for gens(). Returns ideal class generators that generate the class group.
   *
   * @see Reference: sage/rings/number_field/class_group.py:gens
   */
  class_group_generators(): ClassGroupElement[] {
    return this.gens();
  }

  /**
   * Return the number of generators.
   */
  ngens(): number {
    return this._gens.length;
  }

  /**
   * Return the i-th generator.
   */
  gen(i: number): ClassGroupElement {
    if (i < 0 || i >= this._gens.length) {
      throw new ValueError(`generator index ${i} out of range`);
    }
    return this._gens[i]!;
  }

  /**
   * Return the identity element.
   */
  identity(): ClassGroupElement {
    // The identity is represented by zero exponents for each invariant
    return new ClassGroupElement(
      this,
      Array(this._invariants.length)
        .fill(null)
        .map(() => 0n)
    );
  }

  /**
   * Return the class of an ideal.
   *
   * Maps an ideal to its equivalence class in the class group.
   * Two ideals are equivalent if their quotient is a principal ideal.
   *
   * @see Reference: sage/rings/number_field/class_group.py:__call__
   */
  __call__(ideal: NumberFieldIdeal): ClassGroupElement {
    // For trivial class group (class number 1), everything is principal
    if (this.order() === 1n) {
      return this.identity();
    }

    // For principal ideals (1 generator), return identity
    if (ideal.ngens() === 1) {
      return this.identity();
    }

    // For quadratic fields with known class groups, we can compute the class
    const K = this._number_field;
    if (K.degree() === 2) {
      const norm = ideal.norm();
      const invariants = this.invariants();

      // For cyclic class groups, find the exponent
      if (invariants.length === 1) {
        const h = invariants[0]!;

        // The class of an ideal is determined by the prime factorization
        // For small class numbers, we can search
        if (h <= 10n) {
          // Check each class
          for (let e = 0n; e < h; e++) {
            const classElem = new ClassGroupElement(this, [e]);
            // Would need to compare representative ideals
            // For now, return based on norm congruence
          }
        }
      }
    }

    throw new NotImplementedError('__call__ requires PARI bnfisprincipal for general ideals');
  }

  /**
   * Return a random element of the class group.
   *
   * Generates random exponents for each generator.
   */
  random_element(): ClassGroupElement {
    const invariants = this._invariants;
    const exponents: bigint[] = [];

    for (const d of invariants) {
      // Random exponent in [0, d)
      const randomExp = BigInt(Math.floor(Math.random() * Number(d)));
      exponents.push(randomExp);
    }

    return new ClassGroupElement(this, exponents);
  }

  /**
   * Check if the class group is trivial (class number is 1).
   */
  is_trivial(): boolean {
    return this.order() === 1n;
  }

  /**
   * Check if the class group is cyclic.
   */
  is_cyclic(): boolean {
    return this._invariants.length <= 1;
  }

  /**
   * Return a list of all elements.
   *
   * WARNING: This can be very expensive for large class groups.
   */
  list(): ClassGroupElement[] {
    // Generate all combinations of exponents
    const elements: ClassGroupElement[] = [];
    const numInvariants = this._invariants.length;

    if (numInvariants === 0) {
      return [this.identity()];
    }

    const generateElements = (idx: number, exponents: bigint[]): void => {
      if (idx === numInvariants) {
        elements.push(new ClassGroupElement(this, [...exponents]));
        return;
      }

      const order = this._invariants[idx]!;
      for (let e = 0n; e < order; e++) {
        exponents[idx] = e;
        generateElements(idx + 1, exponents);
      }
    };

    generateElements(0, Array(numInvariants).fill(0n));
    return elements;
  }

  /**
   * Iterate over all elements.
   */
  *[Symbol.iterator](): Iterator<ClassGroupElement> {
    for (const element of this.list()) {
      yield element;
    }
  }

  toString(): string {
    const order = this.order();
    if (order === 1n) {
      return `Trivial class group of ${this._number_field}`;
    }

    const structure = this._invariants.map((d) => `C${d}`).join(' x ');
    return `Class group of order ${order} with structure ${structure} of ${this._number_field}`;
  }
}

/**
 * An element of a class group.
 *
 * Represented by exponents with respect to the generators.
 *
 * @see Reference: sage/rings/number_field/class_group.py:FractionalIdealClass
 */
export class ClassGroupElement {
  private readonly _parent: ClassGroup;
  private readonly _exponents: readonly bigint[];
  private _ideal?: NumberFieldIdeal;

  constructor(parent: ClassGroup, exponents: bigint[], ideal?: NumberFieldIdeal) {
    this._parent = parent;
    // Reduce exponents modulo the invariants
    const invariants = parent.invariants();
    const normalizedExps = exponents.map((e, i) => {
      const d = invariants[i];
      if (d === undefined || d === 0n) return e;
      return ((e % d) + d) % d;
    });
    this._exponents = Object.freeze(normalizedExps);
    this._ideal = ideal;
  }

  /**
   * Return the parent class group.
   */
  parent(): ClassGroup {
    return this._parent;
  }

  /**
   * Return a representative ideal for this class.
   */
  ideal(): NumberFieldIdeal {
    if (!this._ideal) {
      throw new NotImplementedError('ideal computation requires PARI');
    }
    return this._ideal;
  }

  /**
   * Return the exponent vector with respect to generators.
   */
  exponents(): bigint[] {
    return [...this._exponents];
  }

  /**
   * Return the order of this element in the class group.
   *
   * The order is the smallest positive integer n such that this^n = identity.
   */
  order(): bigint {
    // The order divides the exponent of the group
    const groupOrder = this._parent.order();
    const invariants = this._parent.invariants();

    if (this.is_trivial()) {
      return 1n;
    }

    // Compute order as LCM of orders in each cyclic factor
    let order = 1n;
    for (let i = 0; i < this._exponents.length; i++) {
      const e = this._exponents[i]!;
      const d = invariants[i]!;
      if (e !== 0n) {
        const g = gcdBigInt(e, d);
        const localOrder = d / g;
        order = lcmBigInt(order, localOrder);
      }
    }

    return order;
  }

  /**
   * Check if this is the trivial class.
   */
  is_trivial(): boolean {
    return this._exponents.every((e) => e === 0n);
  }

  /**
   * Check if this is the trivial class (alias for is_trivial).
   */
  is_principal(): boolean {
    return this.is_trivial();
  }

  /**
   * Multiply two class group elements.
   */
  mul(other: ClassGroupElement): ClassGroupElement {
    if (this._parent !== other._parent) {
      throw new ValueError('elements must be in the same class group');
    }
    const newExponents = this._exponents.map((e, i) => e + other._exponents[i]!);
    return new ClassGroupElement(this._parent, newExponents);
  }

  /**
   * Return the inverse.
   */
  inv(): ClassGroupElement {
    const newExponents = this._exponents.map((e) => -e);
    return new ClassGroupElement(this._parent, newExponents);
  }

  /**
   * Return self raised to power n.
   */
  pow(n: bigint): ClassGroupElement {
    const newExponents = this._exponents.map((e) => e * n);
    return new ClassGroupElement(this._parent, newExponents);
  }

  /**
   * Check equality.
   */
  eq(other: ClassGroupElement): boolean {
    if (this._parent !== other._parent) {
      return false;
    }
    for (let i = 0; i < this._exponents.length; i++) {
      if (this._exponents[i] !== other._exponents[i]) {
        return false;
      }
    }
    return true;
  }

  toString(): string {
    if (this.is_trivial()) {
      return 'Trivial principal fractional ideal class';
    }
    const exps = this._exponents.map((e, i) => `g${i}^${e}`).join(' * ');
    return `Fractional ideal class (${exps})`;
  }
}

/**
 * The narrow (or strict) class group of a number field.
 *
 * The narrow class group uses totally positive principal ideals
 * instead of all principal ideals.
 *
 * @see Reference: sage/rings/number_field/class_group.py:NarrowClassGroup
 */
export class NarrowClassGroup extends ClassGroup {
  // The narrow class group considers only totally positive elements
  // as generators of principal ideals
}

/**
 * The ray class group modulo a modulus m.
 *
 * @see Reference: sage/rings/number_field/class_group.py:RayClassGroup
 */
export class RayClassGroup extends ClassGroup {
  private readonly _modulus: unknown;

  constructor(
    number_field: NumberField,
    modulus: unknown,
    invariants: bigint[] = [],
    gens: ClassGroupElement[] = []
  ) {
    super(number_field, invariants, gens);
    this._modulus = modulus;
  }

  /**
   * Return the modulus.
   */
  modulus(): unknown {
    return this._modulus;
  }
}

// Helper functions

function gcdBigInt(a: bigint, b: bigint): bigint {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    [a, b] = [b, a % b];
  }
  return a;
}

function lcmBigInt(a: bigint, b: bigint): bigint {
  if (a === 0n || b === 0n) return 0n;
  return (a * b) / gcdBigInt(a, b);
}

// ---------------------------------------------------------------------------
// Class groups of quadratic fields, via binary quadratic forms
// ---------------------------------------------------------------------------
//
// SageMath answers `K.class_group()` with PARI's `bnfinit`, and
// `K.class_number()` with `bnfinit`/`qfbclassno`.  `parigp-ts` has neither, so
// the form class group of the (fundamental) discriminant is computed here:
// for a quadratic field the form class group of disc(K) *is* the class group
// (Gauss), and for D > 0 the ordinary class group is the quotient of the
// narrow one by the class of the negative principal form.
//
// This enumerates reduced forms, which is exponential in log|D| where PARI's
// Shanks/McCurley machinery is subexponential; a size guard keeps that honest.
//
// @see Deviation: Quadratic Class Groups Computed From Binary Quadratic Forms

/** Largest |D| for which the reduced-form enumeration is attempted. */
const CLASS_GROUP_DISC_BOUND = 2_000_000n;

type Form = readonly [bigint, bigint, bigint];

function babs(x: bigint): bigint {
  return x < 0n ? -x : x;
}

function isqrtBig(n: bigint): bigint {
  if (n < 0n) throw new ValueError('isqrt of a negative number');
  if (n < 2n) return n;
  let x = n;
  let y = (x + 1n) / 2n;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }
  return x;
}

function xgcdBig(a: bigint, b: bigint): [bigint, bigint, bigint] {
  let [oldR, r] = [a, b];
  let [oldS, s] = [1n, 0n];
  let [oldT, t] = [0n, 1n];
  while (r !== 0n) {
    const q = oldR / r;
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
    [oldT, t] = [t, oldT - q * t];
  }
  if (oldR < 0n) return [-oldR, -oldS, -oldT];
  return [oldR, oldS, oldT];
}

/** Reduce a positive definite form (Gauss). */
function reduceDefinite(f: Form): Form {
  const D = f[1] * f[1] - 4n * f[0] * f[2];
  let [a, b, c] = f;
  for (let guard = 0; guard < 100000; guard++) {
    if (b > a || b <= -a) {
      // Normalise b into (-a, a].
      const m = 2n * a;
      let r = ((b % m) + m) % m; // in [0, 2a)
      if (r > a) r -= m; // in (-a, a]
      b = r;
      c = (b * b - D) / (4n * a);
    }
    if (a > c || (a === c && b < 0n)) {
      [a, b, c] = [c, -b, a];
      continue;
    }
    return [a, b, c];
  }
  throw new ValueError('definite reduction did not terminate');
}

/** Reduce an indefinite form by repeated application of rho (Cohen 5.6.5). */
function rhoIndefinite(f: Form, sqrtD: bigint): Form {
  const [a, b, c] = f;
  const D = b * b - 4n * a * c;
  const cAbs = babs(c);
  const m = 2n * cAbs;
  let r = ((-b % m) + m) % m; // in [0, 2|c|)
  if (cAbs >= sqrtD) {
    // unique r = -b (mod 2c) with -|c| < r <= |c|
    if (r > cAbs) r -= m;
  } else {
    // unique r = -b (mod 2c) with sqrt(D) - 2|c| < r <= sqrt(D)
    const lo = sqrtD - m;
    while (r <= lo) r += m;
    while (r > lo + m) r -= m;
  }
  const newC = (r * r - D) / (4n * c);
  return [c, r, newC];
}

/**
 * `(a, b, c)` is reduced when `|sqrt(D) - 2|a|| < b < sqrt(D)`, equivalently
 * `sqrt(D) - b < 2|a| < sqrt(D) + b` with `b > 0`.  `D` is not a square, so
 * these translate exactly into integer comparisons with `isqrt(D)`.
 */
function isReducedIndefinite(f: Form, sqrtD: bigint): boolean {
  const [a, b] = f;
  if (b <= 0n || b > sqrtD) return false;
  const a2 = 2n * babs(a);
  return a2 + b > sqrtD && a2 - b <= sqrtD;
}

function reduceIndefinite(f: Form, sqrtD: bigint): Form {
  let cur = f;
  for (let i = 0; i < 10000; i++) {
    if (isReducedIndefinite(cur, sqrtD)) return cur;
    cur = rhoIndefinite(cur, sqrtD);
  }
  throw new ValueError('indefinite reduction did not terminate');
}

/**
 * Dirichlet composition of two forms of the same discriminant
 * (Cohen, Algorithm 5.4.7).  The result is not reduced.
 */
function composeForms(f1: Form, f2: Form): Form {
  let [a1, b1, c1] = f1;
  let [a2, b2, c2] = f2;
  void c1;
  if (a1 > a2) {
    [a1, b1, c1, a2, b2, c2] = [a2, b2, c2, a1, b1, c1];
  }
  const s = (b1 + b2) / 2n;
  const n = b2 - s;
  let y1: bigint;
  let d: bigint;
  if (a2 % a1 === 0n) {
    y1 = 0n;
    d = a1;
  } else {
    const [g, u] = xgcdBig(a2, a1);
    y1 = u;
    d = g;
  }
  let x2: bigint;
  let y2: bigint;
  let d1: bigint;
  if (s % d === 0n) {
    y2 = -1n;
    x2 = 0n;
    d1 = d;
  } else {
    const [g, u, v] = xgcdBig(s, d);
    x2 = u;
    y2 = -v;
    d1 = g;
  }
  const v1 = a1 / d1;
  const v2 = a2 / d1;
  let r = (y1 * y2 * n - x2 * c2) % v1;
  if (r < 0n) r += babs(v1);
  const a3 = v1 * v2;
  const b3 = b2 + 2n * v2 * r;
  const c3 = (c2 * d1 + r * (b2 + v2 * r)) / v1;
  return [a3, b3, c3];
}

function formKey(f: Form): string {
  return `${f[0]},${f[1]},${f[2]}`;
}

/** Enumerate the reduced primitive forms of a (negative) discriminant. */
function reducedDefiniteForms(D: bigint): Form[] {
  const out: Form[] = [];
  const bound = isqrtBig(-D / 3n) + 1n;
  for (let a = 1n; a <= bound; a++) {
    for (let b = -a + 1n; b <= a; b++) {
      const t = b * b - D;
      if (t % (4n * a) !== 0n) continue;
      const c = t / (4n * a);
      if (c < a) continue;
      if (gcdBigInt(gcdBigInt(a, b), c) !== 1n) continue;
      if ((a === c || b === a) && b < 0n) continue;
      out.push([a, b, c]);
    }
  }
  return out;
}

/** Enumerate the reduced primitive forms of a positive non-square discriminant. */
function reducedIndefiniteForms(D: bigint): Form[] {
  const sqrtD = isqrtBig(D);
  const out: Form[] = [];
  for (let b = 1n; b <= sqrtD; b++) {
    if ((D - b * b) % 4n !== 0n) continue;
    const prod = (b * b - D) / 4n; // = a*c
    if (prod === 0n) continue;
    const hi = (sqrtD + b) / 2n + 1n;
    for (let aAbs = 1n; aAbs <= hi; aAbs++) {
      if (prod % aAbs !== 0n) continue;
      for (const a of [aAbs, -aAbs]) {
        const c = prod / a;
        if (gcdBigInt(gcdBigInt(babs(a), b), babs(c)) !== 1n) continue;
        const f: Form = [a, b, c];
        if (isReducedIndefinite(f, sqrtD)) out.push(f);
      }
    }
  }
  return out;
}

interface FormClassGroup {
  /** Canonical representatives of the classes. */
  classes: Form[];
  /** Multiply two classes (given by their canonical representatives). */
  compose: (x: Form, y: Form) => Form;
  /** The identity class. */
  identity: Form;
}

function principalForm(D: bigint): Form {
  if (D % 4n === 0n) return [1n, 0n, -D / 4n];
  return [1n, 1n, (1n - D) / 4n];
}

/** Build the form class group of the discriminant D (narrow class group when D > 0). */
function formClassGroup(D: bigint): FormClassGroup {
  if (babs(D) > CLASS_GROUP_DISC_BOUND) {
    throw new NotImplementedError(
      `SAGE_NOT_IMPLEMENTED: class group of discriminant ${D} requires PARI bnfinit/quadclassunit`
    );
  }
  if (D < 0n) {
    const forms = reducedDefiniteForms(D);
    const canonical = (f: Form) => reduceDefinite(f);
    return {
      classes: forms,
      compose: (x, y) => canonical(composeForms(x, y)),
      identity: canonical(principalForm(D)),
    };
  }
  const sqrtD = isqrtBig(D);
  if (sqrtD * sqrtD === D) {
    throw new ValueError('discriminant must not be a perfect square');
  }
  const forms = reducedIndefiniteForms(D);
  // Group the reduced forms into rho-cycles; a class is a cycle.
  const cycleOf = new Map<string, Form>();
  const reps: Form[] = [];
  for (const f of forms) {
    if (cycleOf.has(formKey(f))) continue;
    const cycle: Form[] = [];
    let cur = f;
    do {
      cycle.push(cur);
      cur = rhoIndefinite(cur, sqrtD);
    } while (formKey(cur) !== formKey(f) && cycle.length < 100000);
    let rep = cycle[0]!;
    for (const g of cycle) {
      if (
        g[0] < rep[0] ||
        (g[0] === rep[0] && (g[1] < rep[1] || (g[1] === rep[1] && g[2] < rep[2])))
      ) {
        rep = g;
      }
    }
    for (const g of cycle) cycleOf.set(formKey(g), rep);
    reps.push(rep);
  }
  const canonical = (f: Form): Form => {
    const red = reduceIndefinite(f, sqrtD);
    const rep = cycleOf.get(formKey(red));
    if (rep === undefined) {
      throw new ValueError(`reduced form ${formKey(red)} is not in the enumerated list`);
    }
    return rep;
  };
  return {
    classes: reps,
    compose: (x, y) => canonical(composeForms(x, y)),
    identity: canonical(principalForm(D)),
  };
}

/**
 * Elementary divisors of a finite abelian group given by its element list and
 * multiplication, in PARI/Sage order (decreasing, each dividing the previous).
 */
function abelianInvariants<T>(
  elements: T[],
  key: (x: T) => string,
  compose: (x: T, y: T) => T,
  identity: T
): bigint[] {
  const h = BigInt(elements.length);
  if (h === 1n) return [];
  const idKey = key(identity);
  const orderOf = (x: T): bigint => {
    let cur = x;
    let n = 1n;
    while (key(cur) !== idKey) {
      cur = compose(cur, x);
      n++;
      if (n > h) throw new ValueError('element order exceeds the group order');
    }
    return n;
  };
  const orders = elements.map(orderOf);
  // Prime factorisation of h.
  const primes: bigint[] = [];
  let m = h;
  for (let p = 2n; p * p <= m; p++) {
    if (m % p === 0n) {
      primes.push(p);
      while (m % p === 0n) m /= p;
    }
  }
  if (m > 1n) primes.push(m);

  // For each prime, the p-primary type from the counts |G[p^k]|.
  const parts: bigint[][] = [];
  for (const p of primes) {
    let pk = p;
    const counts: bigint[] = [];
    for (;;) {
      let c = 0n;
      for (const o of orders) {
        if (pk % o === 0n) c++;
      }
      counts.push(c);
      if (counts.length >= 2 && counts[counts.length - 1] === counts[counts.length - 2]) break;
      pk *= p;
    }
    // r_k = log_p(|G[p^k]| / |G[p^(k-1)]|) is the number of parts >= k.
    const rank: number[] = [];
    let prev = 1n;
    for (const c of counts) {
      let ratio = c / prev;
      let r = 0;
      while (ratio > 1n) {
        ratio /= p;
        r++;
      }
      rank.push(r);
      prev = c;
    }
    const maxRank = rank[0] ?? 0;
    const type: bigint[] = [];
    for (let i = 0; i < maxRank; i++) {
      let k = 0;
      while (k < rank.length && rank[k]! > i) k++;
      type.push(p ** BigInt(k));
    }
    parts.push(type);
  }
  // Combine: the j-th largest elementary divisor is the product over primes of
  // the j-th largest p-part.
  let r = 0;
  for (const t of parts) r = Math.max(r, t.length);
  const invariants: bigint[] = [];
  for (let j = 0; j < r; j++) {
    let d = 1n;
    for (const t of parts) {
      if (j < t.length) d *= t[j]!;
    }
    invariants.push(d);
  }
  return invariants;
}

/**
 * Elementary divisors of the class group of the quadratic field of
 * (fundamental) discriminant `D`, in Sage/PARI order.
 */
export function quadraticClassGroupInvariants(D: bigint): bigint[] {
  const G = formClassGroup(D);
  if (D < 0n) {
    return abelianInvariants(G.classes, formKey, G.compose, G.identity);
  }
  // D > 0: the form class group is the *narrow* class group Cl^+.  The
  // ordinary class group is Cl^+ / <[-f0]>, where f0 is the principal form.
  const f0 = principalForm(D);
  const negIdentity = G.compose([-f0[0], f0[1], -f0[2]], G.identity);
  if (formKey(negIdentity) === formKey(G.identity)) {
    // The fundamental unit has norm -1: Cl = Cl^+.
    return abelianInvariants(G.classes, formKey, G.compose, G.identity);
  }
  // Quotient by the order-2 subgroup {1, [-f0]}.
  const cosetRep = new Map<string, string>();
  const reps: Form[] = [];
  for (const c of G.classes) {
    if (cosetRep.has(formKey(c))) continue;
    const other = G.compose(c, negIdentity);
    cosetRep.set(formKey(c), formKey(c));
    cosetRep.set(formKey(other), formKey(c));
    reps.push(c);
  }
  const quotientCompose = (x: Form, y: Form): Form => {
    const prod = G.compose(x, y);
    const repKey = cosetRep.get(formKey(prod));
    if (repKey === undefined) throw new ValueError('coset representative not found');
    const rep = reps.find((f) => formKey(f) === repKey);
    if (rep === undefined) throw new ValueError('coset representative not found');
    return rep;
  };
  const qId = quotientCompose(G.identity, G.identity);
  return abelianInvariants(reps, formKey, quotientCompose, qId);
}

/**
 * The class number of the quadratic field of (fundamental) discriminant `D`.
 *
 * SageMath calls PARI's `qfbclassno`; this counts classes of binary quadratic
 * forms, which is Gauss's theorem in the same form.
 */
export function quadraticClassNumber(D: bigint): bigint {
  return quadraticClassGroupInvariants(D).reduce((a, b) => a * b, 1n);
}
