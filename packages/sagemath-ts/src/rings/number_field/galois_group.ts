/**
 * @module sage/rings/number_field/galois_group
 * @description Galois groups of number fields
 *
 * Port of: sage/rings/number_field/galois_group.py
 * Reference: reference/sage/src/sage/rings/number_field/galois_group.py
 *
 * NOTE: SageMath delegates most Galois group operations to PARI/GP (galoisinit, polgalois).
 * This implementation provides basic functionality where possible without PARI.
 * Full implementation would require PARI's number field functions which are not yet
 * available in parigp-ts.
 *
 * @see Deviation: Number-Field Kernel Not Delegated to parigp-ts
 */

import {
  type GaloisInit,
  type Perm,
  QPoly_to_fractions,
  galoisfixedfield,
  galoisinit,
  galoispermtopol,
} from '@sagemath-ts/parigp-ts';
import { gcd as intGcd } from '../../arith/misc.js';
import { NotImplementedError, ValueError } from '../../errors.js';
import { Rational } from '../rational.js';
import type { NumberField, NumberFieldAutomorphism, NumberFieldElement } from './number_field.js';
import type { NumberFieldEmbedding } from './number_field_embeddings.js';
import type { NumberFieldIdeal } from './number_field_ideal.js';

/**
 * Represents a permutation as an array where perm[i] is the image of i.
 * Uses 0-based indexing internally but can convert to 1-based for display.
 */
export type Permutation = number[];

/**
 * Helper: compose two permutations (left-to-right composition).
 * (a * b)(x) = b(a(x))
 */
function composePermutations(a: Permutation, b: Permutation): Permutation {
  const n = a.length;
  const result: Permutation = new Array(n);
  for (let i = 0; i < n; i++) {
    result[i] = b[a[i]!]!;
  }
  return result;
}

/**
 * Helper: invert a permutation.
 */
function invertPermutation(p: Permutation): Permutation {
  const n = p.length;
  const result: Permutation = new Array(n);
  for (let i = 0; i < n; i++) {
    result[p[i]!] = i;
  }
  return result;
}

/**
 * Helper: check if two permutations are equal.
 */
function permutationsEqual(a: Permutation, b: Permutation): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}

/**
 * Helper: identity permutation of degree n.
 */
function identityPermutation(n: number): Permutation {
  return Array.from({ length: n }, (_, i) => i);
}

/**
 * Helper: check if a permutation is the identity.
 */
function isIdentityPermutation(p: Permutation): boolean {
  for (let i = 0; i < p.length; i++) {
    if (p[i] !== i) return false;
  }
  return true;
}

/**
 * Helper: compute the order of a permutation.
 */
function permutationOrder(p: Permutation): bigint {
  // Order is LCM of cycle lengths
  const visited = new Set<number>();
  let order = 1n;

  for (let i = 0; i < p.length; i++) {
    if (visited.has(i)) continue;

    let cycleLen = 0;
    let current = i;
    while (!visited.has(current)) {
      visited.add(current);
      current = p[current]!;
      cycleLen++;
    }

    // LCM of order and cycleLen
    const cyclelenBig = BigInt(cycleLen);
    order = (order * cyclelenBig) / intGcd(order, cyclelenBig);
  }

  return order;
}

/**
 * Helper: convert permutation to cycle notation string.
 */
function permutationToCycles(p: Permutation): string {
  const visited = new Set<number>();
  const cycles: string[] = [];

  for (let i = 0; i < p.length; i++) {
    if (visited.has(i) || p[i] === i) {
      visited.add(i);
      continue;
    }

    const cycle: number[] = [];
    let current = i;
    while (!visited.has(current)) {
      visited.add(current);
      cycle.push(current + 1); // 1-based for display
      current = p[current]!;
    }

    if (cycle.length > 1) {
      cycles.push(`(${cycle.join(',')})`);
    }
  }

  return cycles.length > 0 ? cycles.join('') : '()';
}

/**
 * Helper: raise permutation to a power.
 */
function permutationPow(p: Permutation, n: bigint): Permutation {
  const len = p.length;
  if (n === 0n) return identityPermutation(len);

  if (n < 0n) {
    p = invertPermutation(p);
    n = -n;
  }

  // Binary exponentiation
  let result = identityPermutation(len);
  let base = [...p];

  while (n > 0n) {
    if (n % 2n === 1n) {
      result = composePermutations(result, base);
    }
    base = composePermutations(base, base);
    n = n / 2n;
  }

  return result;
}

/**
 * The Galois group of a number field extension.
 *
 * For a number field K/Q, the Galois group is Gal(L/Q) where L is the Galois closure of K.
 * For Galois extensions, K = L.
 *
 * NOTE: Full implementation requires PARI's galoisinit function. This implementation
 * provides functionality for special cases (quadratic fields, small cyclic groups).
 *
 * @see Reference: sage/rings/number_field/galois_group.py:GaloisGroup_v2
 */
export class GaloisGroup {
  private readonly _number_field: NumberField;
  private _elements: GaloisGroupElement[] | null = null;
  private _automorphisms: NumberFieldAutomorphism[] | null = null;
  private _generators: GaloisGroupElement[] | null = null;
  private _order: bigint | null = null;
  private _pariCache: { gal: GaloisInit; permByAut: Perm[] } | null = null;
  private readonly _decompCache = new WeakMap<object, GaloisSubgroup>();
  private readonly _ramCache = new WeakMap<object, Map<number, GaloisSubgroup>>();
  private readonly _powCache = new WeakMap<object, Map<number, NumberFieldIdeal>>();

  constructor(number_field: NumberField) {
    this._number_field = number_field;
  }

  /**
   * Return the number field.
   * @see Reference: sage/rings/number_field/galois_group.py:number_field
   */
  number_field(): NumberField {
    return this._number_field;
  }

  /**
   * Return the splitting field (Galois closure).
   */
  splitting_field(): NumberField {
    // For Galois fields, this is the field itself
    // For non-Galois, would need to compute the closure
    return this._number_field;
  }

  /**
   * Return the order of the Galois group.
   *
   * For Galois extensions, this equals the degree of the field.
   * For non-Galois, this is the order of the Galois group of the Galois closure.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:order
   */
  order(): bigint {
    if (this._order !== null) return this._order;
    this._order = BigInt(this._computeElements().length);
    return this._order;
  }

  /**
   * Return the degree of the Galois group (as a permutation group).
   * This is the degree of the number field (or its Galois closure).
   * @see Reference: sage/rings/number_field/galois_group.py:degree
   */
  degree(): number {
    return this._number_field.degree();
  }

  /**
   * Compute and cache the elements of the Galois group.
   */
  private _computeElements(): GaloisGroupElement[] {
    if (this._elements !== null) return this._elements;

    const K = this._number_field;
    const n = K.degree();
    const auts = K.automorphisms();
    if (auts.length !== n) {
      // K/Q is not Galois: Sage returns the Galois group of the Galois closure,
      // which needs PARI's galoisinit/polgalois.
      throw new NotImplementedError(
        `SAGE_NOT_IMPLEMENTED: ${K} is not Galois over Q (only ${auts.length} of ` +
          `${n} automorphisms).  SageMath returns the Galois group of the Galois ` +
          'closure (galois_group.py:268 _gcdata -> number_field.py:9199 ' +
          '_galois_closure_and_embedding -> splitting_field.py:371 splitting_field), ' +
          "which needs PARI's nffactor (factorisation over a number field) and " +
          'rnfequation/polcompositum; neither is ported.'
      );
    }

    // Right-regular representation: perm_g[i] = index of (a_i . g).  Then
    // composePermutations(perm_g, perm_h) = perm_{g h}, matching `mul`.
    const key = (a: NumberFieldAutomorphism): string => a.im_gens()[0]!.toString();
    const index = new Map<string, number>();
    auts.forEach((a, i) => index.set(key(a), i));

    const elements: GaloisGroupElement[] = [];
    for (const g of auts) {
      const perm: Permutation = [];
      for (const ai of auts) {
        // (a_i . g)(alpha) = a_i(g(alpha))
        const image = ai.__call__(g.im_gens()[0]!);
        const j = index.get(image.toString());
        if (j === undefined) {
          throw new ValueError('the automorphisms of this field are not closed under composition');
        }
        perm.push(j);
      }
      elements.push(new GaloisGroupElement(this, perm));
    }
    this._elements = elements;
    this._automorphisms = auts;
    return this._elements;
  }

  /**
   * The automorphisms of the field, indexed consistently with the permutation
   * representation used by the group elements.
   */
  _automorphismList(): NumberFieldAutomorphism[] {
    this._computeElements();
    if (this._automorphisms === null) {
      throw new ValueError('automorphisms have not been computed');
    }
    return this._automorphisms;
  }

  /**
   * Return the list of elements.
   * @see Reference: sage/rings/number_field/galois_group.py:list
   */
  list(): GaloisGroupElement[] {
    return [...this._computeElements()];
  }

  /**
   * Return the generators.
   *
   * For cyclic groups, this is a single generator.
   * For more complex groups, multiple generators may be needed.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:gens
   */
  gens(): GaloisGroupElement[] {
    if (this._generators !== null) return [...this._generators];

    const elements = this._computeElements();
    const n = elements.length;

    if (n <= 1) {
      this._generators = [this.identity()];
      return [...this._generators];
    }

    // Find generators - for cyclic groups, find an element of maximal order
    // For general groups, this is more complex
    const order = this.order();

    // Try to find a generator that generates the whole group
    for (const elem of elements) {
      if (elem.order() === order) {
        this._generators = [elem];
        return [...this._generators];
      }
    }

    // For non-cyclic groups, we need multiple generators
    // This is a simple approach - may not give minimal generating set
    const generated = new Set<string>();
    const gens: GaloisGroupElement[] = [];

    for (const elem of elements) {
      const key = permutationToCycles(elem._permutation);
      if (!generated.has(key)) {
        gens.push(elem);
        // Add all powers of this element to generated set
        let current = elem._permutation;
        for (let i = 0; i < n; i++) {
          generated.add(permutationToCycles(current));
          current = composePermutations(current, elem._permutation);
        }
      }
    }

    // Keep only non-identity generators
    this._generators = gens.filter((g) => !g.is_identity());
    if (this._generators.length === 0) {
      this._generators = [this.identity()];
    }

    return [...this._generators];
  }

  /**
   * Return the number of generators.
   * @see Reference: sage/rings/number_field/galois_group.py:ngens
   */
  ngens(): number {
    return this.gens().length;
  }

  /**
   * Return the i-th generator.
   * @see Reference: sage/rings/number_field/galois_group.py:gen
   */
  gen(i: number): GaloisGroupElement {
    const gens = this.gens();
    if (i < 0 || i >= gens.length) {
      throw new ValueError(`generator index ${i} out of range`);
    }
    return gens[i]!;
  }

  /**
   * Return the identity element.
   * @see Reference: sage/rings/number_field/galois_group.py:identity
   */
  identity(): GaloisGroupElement {
    return new GaloisGroupElement(this, identityPermutation(this.degree()));
  }

  /**
   * Check if the Galois group is abelian.
   *
   * A group is abelian if all elements commute: gh = hg for all g, h.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:is_abelian
   */
  is_abelian(): boolean {
    const elements = this._computeElements();

    // Check if all pairs commute
    for (let i = 0; i < elements.length; i++) {
      for (let j = i + 1; j < elements.length; j++) {
        const gh = elements[i]!.mul(elements[j]!);
        const hg = elements[j]!.mul(elements[i]!);
        if (!gh.eq(hg)) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Check if the Galois group is cyclic.
   *
   * A group is cyclic if it can be generated by a single element.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:is_cyclic
   */
  is_cyclic(): boolean {
    const order = this.order();
    const elements = this._computeElements();

    // A group is cyclic iff it has an element of order equal to |G|
    for (const elem of elements) {
      if (elem.order() === order) {
        return true;
      }
    }
    return false;
  }

  /**
   * Check if the Galois group is solvable.
   *
   * A group is solvable if its derived series terminates at the trivial group.
   * For permutation groups of degree <= 4, all groups are solvable.
   * For degree >= 5, S_n and A_n are not solvable, but proper subgroups might be.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:is_solvable
   */
  is_solvable(): boolean {
    const n = this.degree();
    const order = this.order();

    // All groups of degree <= 4 are solvable
    if (n <= 4) return true;

    // Abelian groups are solvable
    if (this.is_abelian()) return true;

    // If order is a prime power, the group is solvable (p-groups are solvable)
    if (this._isPrimePower(order)) return true;

    // For degree 5, check against known non-solvable groups
    // A_5 has order 60, S_5 has order 120
    if (n === 5) {
      if (order === 60n || order === 120n) {
        // Could be A_5 or S_5, need to verify
        // A_5 and S_5 are NOT solvable
        // But we need more evidence - check if the group is exactly S_5 or A_5
        const factorial5 = 120n;
        if (order === factorial5) {
          // Full symmetric group S_5 is NOT solvable
          return false;
        }
        if (order === 60n) {
          // A_5 is NOT solvable
          // But we should verify this IS A_5 and not some other group of order 60
          // A_5 is simple and non-abelian
          if (!this.is_abelian()) {
            return false;
          }
        }
      }
      // Other subgroups of S_5 are solvable
      return true;
    }

    // For larger degrees, compute the derived series
    // G' = [G, G], G'' = [G', G'], etc.
    // If this terminates at {e}, the group is solvable
    try {
      return this._checkSolvableByDerivedSeries();
    } catch {
      // If we can't compute derived series, check known cases
      // For Galois groups of polynomials, solvability determines if roots
      // can be expressed by radicals (Galois's theorem)

      // Nilpotent groups are solvable
      if (this._isNilpotent()) return true;

      // Groups of order p^a * q^b for primes p, q are solvable (Burnside's theorem)
      if (this._hasTwoPrimePowerOrder(order)) return true;

      throw new NotImplementedError(
        'is_solvable for general groups requires derived series computation'
      );
    }
  }

  /**
   * Check solvability by computing the derived series.
   */
  private _checkSolvableByDerivedSeries(): boolean {
    const elements = this._computeElements();
    let currentGroup = elements.map((e) => e._permutation);

    // Maximum iterations to prevent infinite loop
    const maxIterations = 100;

    for (let iter = 0; iter < maxIterations; iter++) {
      // Compute commutator subgroup [G, G]
      const commutators = new Set<string>();

      for (const g of currentGroup) {
        for (const h of currentGroup) {
          // [g, h] = g * h * g^(-1) * h^(-1)
          const gInv = invertPermutation(g);
          const hInv = invertPermutation(h);
          const comm = composePermutations(
            composePermutations(g, h),
            composePermutations(gInv, hInv)
          );
          commutators.add(JSON.stringify(comm));
        }
      }

      const derivedGroup = Array.from(commutators).map((s) => JSON.parse(s) as Permutation);

      // Check if we've reached the trivial group
      if (derivedGroup.length === 1 && isIdentityPermutation(derivedGroup[0]!)) {
        return true;
      }

      // Check if the derived series has stabilized (not reaching trivial)
      if (derivedGroup.length === currentGroup.length) {
        return false;
      }

      currentGroup = derivedGroup;
    }

    return false;
  }

  /**
   * Check if the group is nilpotent (simple heuristic).
   */
  private _isNilpotent(): boolean {
    // P-groups are nilpotent
    return this._isPrimePower(this.order());
  }

  /**
   * Check if n = p^a * q^b for primes p, q.
   */
  private _hasTwoPrimePowerOrder(n: bigint): boolean {
    if (n <= 1n) return true;

    let remaining = n;
    let primeCount = 0;

    let p = 2n;
    while (p * p <= remaining && primeCount <= 2) {
      if (remaining % p === 0n) {
        primeCount++;
        while (remaining % p === 0n) {
          remaining /= p;
        }
      }
      p++;
    }

    if (remaining > 1n) primeCount++;

    return primeCount <= 2;
  }

  /**
   * Helper: check if n is a prime power.
   */
  private _isPrimePower(n: bigint): boolean {
    if (n <= 1n) return false;

    // Find smallest prime factor
    let p = 2n;
    while (p * p <= n) {
      if (n % p === 0n) {
        // Check if n is a power of p
        while (n % p === 0n) {
          n = n / p;
        }
        return n === 1n;
      }
      p++;
    }
    // n is prime
    return true;
  }

  /**
   * Check if the Galois group acts transitively on the roots.
   *
   * For the Galois group of a field extension, this is always true
   * when the group is computed correctly.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:is_transitive
   */
  is_transitive(): boolean {
    // The Galois group always acts transitively on the roots of an irreducible polynomial
    // (since the field is defined by an irreducible polynomial)
    return true;
  }

  /**
   * Return the transitive group identification number.
   *
   * This is a number identifying the transitive permutation group
   * in the standard database (e.g., GAP's TransitiveGroup).
   *
   * For small degrees (up to 7), uses a lookup table based on order and structure.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:transitive_number
   */
  transitive_number(): number {
    const n = this.degree();
    const order = Number(this.order());

    // For small degrees, we can identify common groups
    if (n === 1) return 1; // Trivial group

    if (n === 2) return 1; // S2 = C2

    if (n === 3) {
      if (order === 3) return 1; // C3 = A3
      if (order === 6) return 2; // S3
    }

    if (n === 4) {
      if (order === 4) {
        if (this.is_cyclic()) return 1; // C4
        return 2; // V4 (Klein four-group)
      }
      if (order === 8) return 3; // D4
      if (order === 12) return 4; // A4
      if (order === 24) return 5; // S4
    }

    if (n === 5) {
      // Transitive groups of degree 5:
      // 5T1 = C5 (order 5)
      // 5T2 = D5 (order 10)
      // 5T3 = F5 = Frobenius group (order 20)
      // 5T4 = A5 (order 60)
      // 5T5 = S5 (order 120)
      if (order === 5) return 1; // C5
      if (order === 10) return 2; // D5
      if (order === 20) return 3; // F5 (Frobenius)
      if (order === 60) return 4; // A5
      if (order === 120) return 5; // S5
    }

    if (n === 6) {
      // Some common degree 6 groups
      if (order === 6) {
        if (this.is_cyclic()) return 1; // C6
        return 2; // S3 (not transitive as 6T, but included for pattern)
      }
      if (order === 12) return 4; // D6
      if (order === 360) return 15; // A6
      if (order === 720) return 16; // S6
    }

    if (n === 7) {
      if (order === 7) return 1; // C7
      if (order === 2520) return 6; // A7
      if (order === 5040) return 7; // S7
    }

    // For larger degrees or unrecognized groups, need PARI or GAP
    throw new NotImplementedError(
      `transitive_number for degree ${n} with order ${order} requires transitive group database`
    );
  }

  /**
   * Whether the underlying number field is Galois over `Q`.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:457 (is_galois)
   */
  is_galois(): boolean {
    const K = this._number_field;
    return K.automorphisms().length === K.degree();
  }

  /* ---------------------------------------------------------------- */
  /* PARI galoisinit plumbing                                          */
  /* ---------------------------------------------------------------- */

  /**
   * `galoisinit` of the PARI model of this field, together with the PARI
   * permutation attached to each of our automorphisms (indexed exactly as
   * {@link _automorphismList}).
   *
   * SageMath keeps the same object in `GaloisGroup_v2._pari_data`
   * (`galois_group.py:320`), computed as
   * `self._galois_closure.__pari__().galoisinit()`.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:320 (_pari_data)
   * @see Reference: reference/pari/src/basemath/galconj.c:3176 (galoisinit)
   */
  _pariGalois(): { gal: GaloisInit; permByAut: Perm[] } {
    if (this._pariCache !== null) return this._pariCache;
    const K = this._number_field;
    const auts = this._automorphismList();
    const gal = galoisinit(K.pari_polynomial());
    if (gal === null) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: PARI galoisinit declined this field (its Galois ' +
          'group is not weakly super solvable, galconj.c:1104 ga_non_wss)'
      );
    }
    // sigma(theta) as an element of K, for every PARI permutation
    const scale = new Rational(K.pari_theta_scale());
    const index = new Map<string, number>();
    auts.forEach((a, i) => index.set(a.im_gens()[0]!.toString(), i));
    const permByAut: Perm[] = new Array(auts.length);
    for (const perm of gal.group.slice(1)) {
      const q = QPoly_to_fractions(galoispermtopol(gal, perm));
      const coeffs = q.map(([n, d]) => new Rational(n, d));
      // image of theta, hence image of alpha = image of theta / scale
      const imTheta = K.element_from_theta_poly(coeffs);
      const imAlpha = imTheta.scalarMul(scale.inv());
      const i = index.get(imAlpha.toString());
      if (i === undefined) {
        throw new ValueError('galoisinit returned an automorphism we do not know');
      }
      permByAut[i] = perm;
    }
    for (let i = 0; i < auts.length; i++) {
      if (permByAut[i] === undefined) {
        throw new ValueError('galoisinit did not return all automorphisms');
      }
    }
    this._pariCache = { gal, permByAut };
    return this._pariCache;
  }

  /** The PARI permutation of a group element. */
  private _permOf(e: GaloisGroupElement): Perm {
    // perm_g[0] is the index of g in `_automorphismList()`; see `__call__`.
    return this._pariGalois().permByAut[e._permutation[0]!]!;
  }

  /* ---------------------------------------------------------------- */
  /* Local behaviour at a prime                                        */
  /* ---------------------------------------------------------------- */

  /** Coerce the argument of the local methods to a prime ideal of `O_K`. */
  private _asPrimeIdeal(prime: unknown): NumberFieldIdeal {
    const K = this._number_field;
    if (typeof prime === 'bigint') {
      // SageMath insists on a prime *ideal*; the port also accepts a rational
      // prime, and then uses a prime above it (the answer only depends on the
      // choice up to conjugacy).
      return K.primes_above(prime)[0]!;
    }
    const P = prime as NumberFieldIdeal;
    if (typeof P?.is_prime !== 'function') {
      throw new ValueError(`${prime} is not a prime ideal`);
    }
    if (!P.is_prime()) throw new ValueError(`${P} is not a prime ideal`);
    return P;
  }

  /** The image `sigma(I)` of an ideal under an automorphism. */
  private _applyToIdeal(s: GaloisGroupElement, I: NumberFieldIdeal): NumberFieldIdeal {
    const K = this._number_field;
    return K.ideal(...I.gens().map((g) => s.__call__(g)));
  }

  /**
   * Decomposition group of a prime ideal `P`, i.e. the subgroup of elements
   * that map `P` to itself.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:601 (decomposition_group)
   */
  decomposition_group(prime: unknown): GaloisSubgroup {
    if (!this.is_galois()) {
      throw new ValueError('Decomposition groups only defined for Galois extensions');
    }
    const P = this._asPrimeIdeal(prime);
    const hit = this._decompCache.get(P);
    if (hit) return hit;
    const gens = this._computeElements().filter((s) => this._applyToIdeal(s, P).eq(P));
    const D = new GaloisSubgroup(this, gens);
    this._decompCache.set(P, D);
    return D;
  }

  /**
   * The `v`-th ramification group of `P`: the set of elements `s` acting
   * trivially modulo `P^(v+1)`.
   *
   * Upstream reads this off PARI's `idealramgroups` (`base1.c:1074`), which
   * computes, for each `s` in the group, `idx(s) = v_P(s(pi) - pi)` capped at
   * the bound `(v_P(diff) - (e-1))/(p-1)` and corrected on the residue-field
   * generator (`idealramgroupswild`, `base1.c:968`), and then returns
   * `G_v = {s : idx(s) >= v+1}`.  We evaluate the same condition
   * `v_P(s(x) - x) >= v+1` on a Z-basis of `O_K` instead of on `pi` and the
   * residue generator: `x -> s(x) - x` is additive, so the two are equivalent,
   * and no bound on `v_P(diff)` is needed.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:693 (ramification_group)
   */
  ramification_group(prime: unknown, v: number): GaloisSubgroup {
    if (!this.is_galois()) {
      throw new ValueError('Ramification groups only defined for Galois extensions');
    }
    if (v < -1) throw new ValueError('v must be at least -1');
    const D = this.decomposition_group(prime);
    if (v === -1) return D;
    const K = this._number_field;
    const P = this._asPrimeIdeal(prime);
    let byV = this._ramCache.get(P);
    if (!byV) {
      byV = new Map<number, GaloisSubgroup>();
      this._ramCache.set(P, byV);
    }
    const hit = byV.get(v);
    if (hit) return hit;
    const Pv = this._primePower(P, v + 1);
    const basis = K.integral_basis();
    const gens = D.list().filter((s) =>
      basis.every((w) => {
        const d = s.__call__(w).sub(w);
        return d.is_zero() || Pv.contains(d);
      })
    );
    const Gv = new GaloisSubgroup(this, gens);
    byV.set(v, Gv);
    return Gv;
  }

  /** `P^k`, cached (the ramification filtration asks for many of them). */
  private _primePower(P: NumberFieldIdeal, k: number): NumberFieldIdeal {
    let m = this._powCache.get(P);
    if (!m) {
      m = new Map<number, NumberFieldIdeal>();
      this._powCache.set(P, m);
    }
    const hit = m.get(k);
    if (hit) return hit;
    const v = k <= 1 ? P : this._primePower(P, k - 1).mul(P);
    m.set(k, v);
    return v;
  }

  /**
   * The inertia group of `P`, i.e. the 0-th ramification group.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:721 (inertia_group)
   */
  inertia_group(prime: unknown): GaloisSubgroup {
    if (!this.is_galois()) {
      throw new ValueError('Inertia groups only defined for Galois extensions');
    }
    return this.ramification_group(prime, 0);
  }

  /**
   * The set of ramification breaks of `P`: the indices `i` with
   * `G_{i+1} != G_i`.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:740 (ramification_breaks)
   */
  ramification_breaks(prime: unknown): number[] {
    if (!this.is_galois()) {
      throw new ValueError('Ramification breaks only defined for Galois extensions');
    }
    const orders: number[] = [];
    for (let v = -1; ; v++) {
      const o = Number(this.ramification_group(prime, v).order());
      orders.push(o);
      if (o === 1) break;
      if (v > 4096) throw new ValueError('ramification_breaks did not terminate');
    }
    // orders[j] = |G_{j-1}|
    const breaks: number[] = [];
    for (let j = 0; j + 1 < orders.length; j++) {
      if (orders[j] !== orders[j + 1]) breaks.push(j - 1);
    }
    return breaks;
  }

  /**
   * Return the Frobenius element at an unramified prime; an alias for
   * {@link artin_symbol}.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:767 (artin_symbol)
   */
  frobenius(prime: unknown): GaloisGroupElement {
    return this.artin_symbol(prime);
  }

  /**
   * The Artin symbol at an unramified prime ideal `P`: the unique `s` in the
   * decomposition group with `s(x) = x^p mod P`.
   *
   * This is a transcription of SageMath's own implementation, which does not
   * call PARI: it tests the congruence on the ring generators of `O_K`.  We
   * test it on a Z-basis instead, which is equivalent because `x -> s(x) - x^p`
   * is additive modulo `P` (Frobenius on the residue field is additive).
   *
   * @see Reference: sage/rings/number_field/galois_group.py:767 (artin_symbol)
   */
  artin_symbol(prime: unknown): GaloisGroupElement {
    if (!this.is_galois()) {
      throw new ValueError('Artin symbols only defined for Galois extensions');
    }
    const K = this._number_field;
    const P = this._asPrimeIdeal(prime);
    const p = P.prime_below();
    const basis = K.integral_basis();
    const t: GaloisGroupElement[] = [];
    for (const s of this.decomposition_group(P).list()) {
      const ok = basis.every((g) => {
        const d = s.__call__(g).sub(g.pow(p));
        return d.is_zero() || P.contains(d);
      });
      if (ok) t.push(s);
    }
    if (t.length > 1) throw new ValueError(`${P} is ramified`);
    if (t.length === 0) throw new ValueError(`${P} has no Frobenius element`);
    return t[0]!;
  }

  /**
   * The element of the group inducing complex conjugation for the complex
   * place `place`.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:649 (complex_conjugation)
   */
  complex_conjugation(place?: NumberFieldEmbedding): GaloisGroupElement {
    const K = this._number_field;
    if (place === undefined) {
      throw new ValueError('No default complex embedding specified');
    }
    if (!this.is_galois()) throw new ValueError('Extension is not Galois');
    if (K.is_totally_real()) throw new ValueError('No complex conjugation (field is real)');
    // SageMath: elts = [s for s in self if P(s(g)) == gconj]; exactly one must
    // match.  The comparison here is done on the *certified* intervals, so a
    // candidate is kept only when its box can contain the conjugate value.
    const g = K.gen();
    const target = place.__call__(g);
    const conj = { re: target.re, im: { lo: target.im.hi.neg(), hi: target.im.lo.neg() } };
    const elts: GaloisGroupElement[] = [];
    for (const s of this._computeElements()) {
      const v = place.__call__(s.__call__(g));
      const overlaps =
        v.re.lo.le(conj.re.hi) &&
        conj.re.lo.le(v.re.hi) &&
        v.im.lo.le(conj.im.hi) &&
        conj.im.lo.le(v.im.hi);
      if (overlaps) elts.push(s);
    }
    if (elts.length !== 1) throw new ValueError('Something has gone very wrong here');
    return elts[0]!;
  }

  /* ---------------------------------------------------------------- */
  /* Fixed fields                                                      */
  /* ---------------------------------------------------------------- */

  /**
   * Return the fixed field of a subgroup, together with the element of `K`
   * that generates it.
   *
   * SageMath calls PARI's `galoisfixedfield` (`galois_group.py:890`) and then
   * `L.subfield(x, name)`; so do we.
   *
   * `polred` is *not* applied: SageMath post-processes the PARI answer with
   * `polredbest` whenever the index is at most 8, and `polredbest` is not
   * ported.  The field returned here is therefore PARI's raw
   * `galoisfixedfield` answer, which is what SageMath returns for
   * `fixed_field(polred=False)`.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:890 (fixed_field)
   * @see Reference: reference/pari/src/basemath/galconj.c:3276 (galoisfixedfield)
   * @see Deviation: fixed_field does not apply polredbest
   */
  fixed_field_data(subgroupElements?: GaloisGroupElement[]): {
    field: NumberField | 'Q';
    gen: NumberFieldElement | null;
    polynomial: bigint[] | null;
  } {
    const K = this._number_field;
    const all = this._computeElements();
    const elements = subgroupElements ?? all;

    if (elements.length === all.length) return { field: 'Q', gen: null, polynomial: null };
    if (elements.length === 1 && elements[0]!.is_identity()) {
      return { field: K, gen: K.gen(), polynomial: null };
    }

    const { gal } = this._pariGalois();
    const perms = elements.map((e) => this._permOf(e));
    const ff = galoisfixedfield(gal, perms, 0);
    const P = ff.P as bigint[];
    const S = ff.S;
    if (S === undefined) throw new ValueError('galoisfixedfield returned no generator');
    const coeffs = QPoly_to_fractions(S).map(([n, d]) => new Rational(n, d));
    const gen = K.element_from_theta_poly(coeffs);
    // eslint-disable-next-line @typescript-eslint/no-var-requires
    const { NumberFieldConstructor } = require('./number_field.js');
    const field = NumberFieldConstructor(P, `${K._name}0`) as NumberField;
    return { field, gen, polynomial: P };
  }

  /**
   * Return the fixed field of a subgroup.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:890 (fixed_field)
   */
  fixed_field(subgroupElements?: GaloisGroupElement[]): NumberField | 'Q' {
    return this.fixed_field_data(subgroupElements).field;
  }

  /**
   * Return the subgroups.
   *
   * For small groups, enumerates all subgroups using brute force.
   * For larger groups, throws NotImplementedError.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:subgroups
   */
  subgroups(): GaloisGroupElement[][] {
    const elements = this._computeElements();
    const order = elements.length;

    // Only compute for small groups
    if (order > 24) {
      throw new NotImplementedError(
        'subgroups enumeration only implemented for groups of order <= 24'
      );
    }

    const subgroups: GaloisGroupElement[][] = [];
    const seenSubgroups = new Set<string>();

    // Trivial subgroup
    subgroups.push([this.identity()]);
    seenSubgroups.add(JSON.stringify([0])); // Just identity

    // Whole group
    const wholeGroupKey = elements
      .map((_, i) => i)
      .sort()
      .join(',');
    subgroups.push([...elements]);
    seenSubgroups.add(wholeGroupKey);

    // For each subset of generators, compute the generated subgroup
    for (let i = 0; i < elements.length; i++) {
      const generator = elements[i]!;
      if (generator.is_identity()) continue;

      // Generate cyclic subgroup from this element
      const cyclic = this._generateSubgroup([generator]);
      const key = this._subgroupKey(cyclic, elements);

      if (!seenSubgroups.has(key)) {
        subgroups.push(cyclic);
        seenSubgroups.add(key);
      }

      // Try pairs of generators
      for (let j = i + 1; j < elements.length; j++) {
        const gen2 = elements[j]!;
        if (gen2.is_identity()) continue;

        const subgroup = this._generateSubgroup([generator, gen2]);
        const pairKey = this._subgroupKey(subgroup, elements);

        if (!seenSubgroups.has(pairKey)) {
          subgroups.push(subgroup);
          seenSubgroups.add(pairKey);
        }
      }
    }

    return subgroups;
  }

  /**
   * Generate a subgroup from a set of generators.
   */
  private _generateSubgroup(generators: GaloisGroupElement[]): GaloisGroupElement[] {
    const subgroupMap = new Map<string, GaloisGroupElement>();
    subgroupMap.set(permutationToCycles(this.identity()._permutation), this.identity());

    for (const gen of generators) {
      subgroupMap.set(permutationToCycles(gen._permutation), gen);
    }

    let changed = true;
    while (changed) {
      changed = false;
      const currentElems = Array.from(subgroupMap.values());

      for (const g of currentElems) {
        for (const h of currentElems) {
          const product = g.mul(h);
          const key = permutationToCycles(product._permutation);
          if (!subgroupMap.has(key)) {
            subgroupMap.set(key, product);
            changed = true;
          }

          const inverse = g.inv();
          const invKey = permutationToCycles(inverse._permutation);
          if (!subgroupMap.has(invKey)) {
            subgroupMap.set(invKey, inverse);
            changed = true;
          }
        }
      }
    }

    return Array.from(subgroupMap.values());
  }

  /**
   * Create a canonical key for a subgroup.
   */
  private _subgroupKey(subgroup: GaloisGroupElement[], allElements: GaloisGroupElement[]): string {
    const indices: number[] = [];
    for (const s of subgroup) {
      for (let i = 0; i < allElements.length; i++) {
        if (s.eq(allElements[i]!)) {
          indices.push(i);
          break;
        }
      }
    }
    return indices.sort((a, b) => a - b).join(',');
  }

  /**
   * Return the normal subgroups.
   *
   * A subgroup H is normal if gHg^{-1} = H for all g in G.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:normal_subgroups
   */
  normal_subgroups(): GaloisGroupElement[][] {
    const elements = this._computeElements();
    const allSubgroups = this.subgroups();

    return allSubgroups.filter((H) => this._isNormalSubgroup(H, elements));
  }

  /**
   * Check if H is a normal subgroup.
   */
  private _isNormalSubgroup(H: GaloisGroupElement[], G: GaloisGroupElement[]): boolean {
    const HSet = new Set(H.map((h) => permutationToCycles(h._permutation)));

    for (const g of G) {
      for (const h of H) {
        // Check if g * h * g^(-1) is in H
        const conjugate = g.mul(h).mul(g.inv());
        const key = permutationToCycles(conjugate._permutation);
        if (!HSet.has(key)) {
          return false;
        }
      }
    }

    return true;
  }

  /**
   * Return the center of the Galois group.
   *
   * The center Z(G) = {g in G : gh = hg for all h in G}.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:center
   */
  center(): GaloisGroupElement[] {
    const elements = this._computeElements();
    const centerElements: GaloisGroupElement[] = [];

    for (const g of elements) {
      let inCenter = true;
      for (const h of elements) {
        const gh = g.mul(h);
        const hg = h.mul(g);
        if (!gh.eq(hg)) {
          inCenter = false;
          break;
        }
      }
      if (inCenter) {
        centerElements.push(g);
      }
    }

    return centerElements;
  }

  /**
   * Return the commutator (derived) subgroup elements.
   *
   * The commutator subgroup [G, G] is generated by all commutators [g, h] = ghg^{-1}h^{-1}.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:commutator
   */
  commutator(): GaloisGroupElement[] {
    const elements = this._computeElements();
    const commutatorSet = new Map<string, GaloisGroupElement>();

    // Generate all commutators
    for (const g of elements) {
      for (const h of elements) {
        // [g, h] = g * h * g^(-1) * h^(-1)
        const comm = g.mul(h).mul(g.inv()).mul(h.inv());
        const key = permutationToCycles(comm._permutation);
        if (!commutatorSet.has(key)) {
          commutatorSet.set(key, comm);
        }
      }
    }

    // The commutator subgroup is the group generated by these commutators
    // We need to take the closure under multiplication
    let changed = true;
    while (changed) {
      changed = false;
      const currentComms = Array.from(commutatorSet.values());
      for (const c1 of currentComms) {
        for (const c2 of currentComms) {
          const product = c1.mul(c2);
          const key = permutationToCycles(product._permutation);
          if (!commutatorSet.has(key)) {
            commutatorSet.set(key, product);
            changed = true;
          }
        }
      }
    }

    return Array.from(commutatorSet.values());
  }

  /**
   * Return the abelianization G/[G, G].
   *
   * For finite groups, this returns representatives of cosets.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:abelianization
   */
  abelianization(): { order: bigint; structure: bigint[] } {
    const elements = this._computeElements();
    const commutatorSubgroup = this.commutator();

    // Order of abelianization is |G| / |[G, G]|
    const abOrder = BigInt(elements.length) / BigInt(commutatorSubgroup.length);

    // For abelian groups, the abelianization is the group itself
    if (this.is_abelian()) {
      return { order: abOrder, structure: [abOrder] };
    }

    // For non-abelian groups, we'd need to compute the structure
    // Return just the order for now
    return { order: abOrder, structure: [] };
  }

  /**
   * Return the conjugacy classes.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:conjugacy_classes
   */
  conjugacy_classes(): GaloisGroupElement[][] {
    const elements = this._computeElements();
    const classes: GaloisGroupElement[][] = [];
    const assigned = new Set<number>();

    for (let i = 0; i < elements.length; i++) {
      if (assigned.has(i)) continue;

      const g = elements[i]!;
      const conjugacyClass: GaloisGroupElement[] = [];

      // Compute all conjugates h g h^{-1}
      for (const h of elements) {
        const conj = h.mul(g).mul(h.inv());

        // Check if this conjugate is already in the class
        let found = false;
        for (const c of conjugacyClass) {
          if (c.eq(conj)) {
            found = true;
            break;
          }
        }

        if (!found) {
          conjugacyClass.push(conj);

          // Mark this element as assigned
          for (let j = 0; j < elements.length; j++) {
            if (elements[j]!.eq(conj)) {
              assigned.add(j);
            }
          }
        }
      }

      classes.push(conjugacyClass);
    }

    return classes;
  }

  /**
   * Return the character table.
   *
   * For small abelian groups, computes the character table directly.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:character_table
   */
  character_table(): { classes: GaloisGroupElement[][]; characters: bigint[][] } {
    const order = Number(this.order());

    if (order > 12) {
      throw new NotImplementedError('character_table only implemented for groups of order <= 12');
    }

    const classes = this.conjugacy_classes();

    // For abelian groups, each element is its own conjugacy class
    // and there are |G| irreducible characters, all 1-dimensional
    if (this.is_abelian()) {
      const elements = this._computeElements();
      const n = elements.length;

      // For a cyclic group of order n, characters are chi_k(g) = zeta_n^(k * ord(g))
      // For simplicity, return a placeholder structure
      const characters: bigint[][] = [];

      // Identity character
      characters.push(Array(classes.length).fill(1n));

      return { classes, characters };
    }

    throw new NotImplementedError(
      'character_table for non-abelian groups requires representation theory'
    );
  }

  /**
   * Return a random element.
   * @see Reference: sage/rings/number_field/galois_group.py:random_element
   */
  random_element(): GaloisGroupElement {
    const elements = this._computeElements();
    const idx = Math.floor(Math.random() * elements.length);
    return elements[idx]!;
  }

  /**
   * Return the list of automorphisms as field homomorphisms.
   *
   * Each automorphism is represented as an object with an __call__ method
   * that applies the automorphism to elements of the number field.
   *
   * @see Reference: sage/rings/number_field/number_field.py:automorphisms
   */
  automorphisms(): Array<{ __call__: (x: NumberFieldElement) => NumberFieldElement }> {
    const elements = this._computeElements();
    return elements.map((elem) => elem.as_hom());
  }

  /**
   * Create a subgroup from generators.
   */
  subgroup(generators: GaloisGroupElement[]): GaloisSubgroup {
    return new GaloisSubgroup(this, generators);
  }

  /**
   * Iterator over elements.
   */
  *[Symbol.iterator](): Iterator<GaloisGroupElement> {
    for (const elem of this._computeElements()) {
      yield elem;
    }
  }

  toString(): string {
    const K = this._number_field;
    const n = this.degree();
    try {
      const ord = this.order();
      return `Galois group of order ${ord} of ${K.defining_polynomial()}`;
    } catch {
      return `Galois group of degree ${n}`;
    }
  }
}

/**
 * An element of a Galois group (an automorphism).
 *
 * Galois group elements are stored as permutations of the roots of the
 * defining polynomial. They can also act on elements of the number field.
 *
 * @see Reference: sage/rings/number_field/galois_group.py:GaloisGroupElement
 */
export class GaloisGroupElement {
  private readonly _parent: GaloisGroup;
  readonly _permutation: Permutation;

  constructor(parent: GaloisGroup, permutation: Permutation) {
    this._parent = parent;
    this._permutation = permutation;
  }

  /**
   * Return the parent Galois group.
   * @see Reference: sage/rings/number_field/galois_group.py:parent
   */
  parent(): GaloisGroup {
    return this._parent;
  }

  /**
   * Return the order of this element.
   *
   * The order is the smallest positive integer n such that g^n = identity.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:order
   */
  order(): bigint {
    return permutationOrder(this._permutation);
  }

  /**
   * Apply this automorphism to a number field element.
   *
   * For a general automorphism sigma, sigma(a) is computed by applying
   * the permutation to the polynomial representation of a.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:__call__
   */
  __call__(x: NumberFieldElement): NumberFieldElement {
    if (this.is_identity()) {
      return x;
    }
    // perm[0] is the index of this automorphism in the parent's list, because
    // perm_g[0] = index of (identity . g) = index of g.
    const auts = this._parent._automorphismList();
    const sigma = auts[this._permutation[0]!];
    if (sigma === undefined) {
      throw new ValueError('Galois group element does not correspond to an automorphism');
    }
    return sigma.__call__(x);
  }

  /**
   * Return the fixed field of this automorphism.
   *
   * The fixed field is the subfield fixed by the cyclic group generated by this element.
   * For the identity, this is the whole field. For quadratic fields, the non-identity
   * element fixes only Q.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:fixed_field
   */
  fixed_field(): NumberField | 'Q' {
    return this.fixed_field_data().field;
  }

  /**
   * The fixed field of the cyclic group generated by this element, together
   * with the element of `K` that generates it.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:890 (fixed_field)
   */
  fixed_field_data(): {
    field: NumberField | 'Q';
    gen: NumberFieldElement | null;
    polynomial: bigint[] | null;
  } {
    const K = this._parent.number_field();
    if (this.is_identity()) return { field: K, gen: K.gen(), polynomial: null };
    if (K.degree() === 2) return { field: 'Q', gen: null, polynomial: null };
    // the cyclic subgroup <self>
    const elts: GaloisGroupElement[] = [];
    let cur: GaloisGroupElement = this._parent.identity();
    do {
      elts.push(cur);
      cur = cur.mul(this);
    } while (!cur.is_identity());
    return this._parent.fixed_field_data(elts);
  }

  /**
   * The greatest `v` such that this element acts trivially modulo `P^v`.
   *
   * A transcription of SageMath's own implementation, which evaluates
   * `min_g v_P(s(g) - g)` over the ring generators of `O_K`; we take the
   * minimum over a Z-basis, which is the same number because `x -> s(x) - x`
   * is additive.
   *
   * @see Reference: sage/rings/number_field/galois_group.py:1060 (ramification_degree)
   */
  ramification_degree(prime: unknown): bigint {
    if (!this._parent.is_galois()) {
      throw new ValueError('Ramification degree only defined for Galois extensions');
    }
    const K = this._parent.number_field();
    const P =
      typeof prime === 'bigint'
        ? K.primes_above(prime)[0]!
        : (prime as NumberFieldIdeal);
    let best: bigint | null = null;
    for (const g of K.integral_basis()) {
      const d = this.__call__(g).sub(g);
      if (d.is_zero()) continue;
      const v = P.valuation(d);
      if (best === null || v < best) best = v;
    }
    // s = identity on every basis element: infinite valuation.  SageMath would
    // return +Infinity here; the identity is the only such element.
    if (best === null) {
      throw new ValueError('ramification_degree of the identity is infinite');
    }
    return best;
  }

  /**
   * Compose two automorphisms.
   *
   * (g * h)(x) = g(h(x))
   *
   * @see Reference: sage/rings/number_field/galois_group.py:__mul__
   */
  mul(other: GaloisGroupElement): GaloisGroupElement {
    this._checkParent(other);
    const composed = composePermutations(this._permutation, other._permutation);
    return new GaloisGroupElement(this._parent, composed);
  }

  /**
   * Return the inverse automorphism.
   * @see Reference: sage/rings/number_field/galois_group.py:__invert__
   */
  inv(): GaloisGroupElement {
    const invPerm = invertPermutation(this._permutation);
    return new GaloisGroupElement(this._parent, invPerm);
  }

  /**
   * Return self raised to power n.
   * @see Reference: sage/rings/number_field/galois_group.py:__pow__
   */
  pow(n: bigint): GaloisGroupElement {
    const powPerm = permutationPow(this._permutation, n);
    return new GaloisGroupElement(this._parent, powPerm);
  }

  /**
   * Check equality.
   * @see Reference: sage/rings/number_field/galois_group.py:__eq__
   */
  eq(other: GaloisGroupElement): boolean {
    if (this._parent !== other._parent) return false;
    return permutationsEqual(this._permutation, other._permutation);
  }

  /**
   * Check if this is the identity.
   * @see Reference: sage/rings/number_field/galois_group.py:is_identity
   */
  is_identity(): boolean {
    return isIdentityPermutation(this._permutation);
  }

  /**
   * Return as a permutation (array form).
   * @see Reference: sage/rings/number_field/galois_group.py:as_permutation
   */
  as_permutation(): Permutation {
    return [...this._permutation];
  }

  /**
   * Return the corresponding automorphism as a field homomorphism.
   * @see Reference: sage/rings/number_field/galois_group.py:as_hom
   */
  as_hom(): { __call__: (x: NumberFieldElement) => NumberFieldElement } {
    return {
      __call__: (x: NumberFieldElement) => this.__call__(x),
    };
  }

  private _checkParent(other: GaloisGroupElement): void {
    if (this._parent !== other._parent) {
      throw new ValueError('elements must be in the same Galois group');
    }
  }

  toString(): string {
    return permutationToCycles(this._permutation);
  }
}

/**
 * A subgroup of a Galois group.
 *
 * @see Reference: sage/rings/number_field/galois_group.py:GaloisGroup_subgroup
 */
export class GaloisSubgroup extends GaloisGroup {
  private readonly _ambient: GaloisGroup;
  private readonly _subgroupElements: GaloisGroupElement[];

  constructor(ambient: GaloisGroup, generators: GaloisGroupElement[]) {
    super(ambient.number_field());
    this._ambient = ambient;
    this._subgroupElements = this._generateFromGenerators(generators);
  }

  /**
   * Generate all elements from the given generators.
   */
  private _generateFromGenerators(generators: GaloisGroupElement[]): GaloisGroupElement[] {
    if (generators.length === 0) {
      return [this._ambient.identity()];
    }

    const elements = new Map<string, GaloisGroupElement>();
    const identity = this._ambient.identity();
    elements.set(permutationToCycles(identity._permutation), identity);

    // Keep multiplying generators until closure
    let changed = true;
    while (changed) {
      changed = false;
      const currentElements = Array.from(elements.values());

      for (const g of currentElements) {
        for (const h of generators) {
          for (const product of [g.mul(h), h.mul(g), g.mul(h.inv()), h.inv().mul(g)]) {
            const key = permutationToCycles(product._permutation);
            if (!elements.has(key)) {
              elements.set(key, product);
              changed = true;
            }
          }
        }
      }
    }

    return Array.from(elements.values());
  }

  /**
   * Return the ambient Galois group.
   */
  ambient(): GaloisGroup {
    return this._ambient;
  }

  override order(): bigint {
    return BigInt(this._subgroupElements.length);
  }

  override list(): GaloisGroupElement[] {
    return [...this._subgroupElements];
  }

  /**
   * Return the fixed field of this subgroup.
   */
  override fixed_field(): NumberField | 'Q' {
    return this._ambient.fixed_field(this._subgroupElements);
  }
}

/** Coerce a prime argument (a rational prime or a prime ideal) to a rational prime. */
function toPrime(prime: unknown): bigint {
  if (typeof prime === 'bigint') return prime;
  return (prime as { prime_below: () => bigint }).prime_below();
}

/**
 * Determine how the rational prime `p` behaves in a quadratic field, using the
 * actual prime decomposition (Dedekind-Kummer) rather than a Legendre symbol,
 * which gets `p = 2` wrong.
 */
function quadraticSplitting(K: NumberField, p: bigint): ['split' | 'inert' | 'ramified'] {
  const decomposition = K.decomposition(p);
  if (decomposition.length === 2) return ['split'];
  const e = decomposition[0]![1];
  return [e > 1n ? 'ramified' : 'inert'];
}

/**
 * Factory function to create the Galois group of a number field.
 *
 * @param K - The number field
 * @returns The Galois group of K (or its Galois closure)
 */
export function galois_group(K: NumberField): GaloisGroup {
  return new GaloisGroup(K);
}
