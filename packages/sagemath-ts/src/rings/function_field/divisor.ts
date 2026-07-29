/**
 * @module sage/rings/function_field/divisor
 * @description Divisors of function fields
 *
 * Port of: sage/rings/function_field/divisor.py
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import { divide_constants } from './constant_field.js';
import type { ConstantFieldElement } from './constant_field.js';
import type { FunctionFieldElement_rational } from './element_rational.js';
import type { FunctionField } from './function_field.js';
import type { FunctionFieldPlace } from './place.js';

/** Entry of a divisor: a place together with its (nonzero) multiplicity. */
type DivisorEntry<C extends ConstantFieldElement> = [FunctionFieldPlace<C>, bigint];

/**
 * Construct a divisor from the data.
 *
 * @see Reference: sage/rings/function_field/divisor.py:78 (divisor)
 */
export function divisor<C extends ConstantFieldElement>(
  field: FunctionField<C>,
  data: Iterable<DivisorEntry<C>>
): FunctionFieldDivisor<C> {
  const divisor_group = field.divisor_group();
  return new FunctionFieldDivisor(divisor_group, data);
}

/**
 * Construct a prime divisor from the place.
 *
 * @see Reference: sage/rings/function_field/divisor.py:103 (prime_divisor)
 */
export function prime_divisor<C extends ConstantFieldElement>(
  field: FunctionField<C>,
  place: FunctionFieldPlace<C>,
  m: bigint | number = 1
): FunctionFieldDivisor<C> {
  const divisor_group = field.divisor_group();
  return new FunctionFieldDivisor(divisor_group, [[place, BigInt(m)]]);
}

/**
 * Divisors of function fields.
 *
 * @see Reference: sage/rings/function_field/divisor.py:129 (FunctionFieldDivisor)
 */
export class FunctionFieldDivisor<C extends ConstantFieldElement> {
  readonly _parent: DivisorGroup<C>;
  /** Place key -> (place, multiplicity).  Zero multiplicities are never stored. */
  readonly _data: Map<string, DivisorEntry<C>>;

  private _functionSpaceCache:
    | [Array<FunctionFieldElement_rational<C>>, (f: FunctionFieldElement_rational<C>) => C[]]
    | null = null;
  private _differentialSpaceCache:
    | [Array<FunctionFieldElement_rational<C>>, (f: FunctionFieldElement_rational<C>) => C[]]
    | null = null;

  constructor(parent: DivisorGroup<C>, data: Iterable<DivisorEntry<C>>) {
    this._parent = parent;
    // Upstream's `FunctionFieldDivisor.__init__` stores the dictionary as
    // given, zero multiplicities included: `p.divisor(0)` really does print as
    // `0*Place (...)` and has `p` in its support.  Only `_add_` and `_rmul_`
    // drop zeros.  Verified against SageMath 10.3.
    this._data = new Map();
    for (const [place, m] of data) {
      this._data.set(place._key(), [place, m]);
    }
  }

  parent(): DivisorGroup<C> {
    return this._parent;
  }

  /** The function field this divisor lives on. */
  function_field(): FunctionField<C> {
    return this._parent._field;
  }

  private _entries(): Array<DivisorEntry<C>> {
    const entries = [...this._data.values()];
    entries.sort((a, b) => a[0].cmp(b[0]));
    return entries;
  }

  /**
   * Return a string representation of the divisor.
   *
   * @see Reference: sage/rings/function_field/divisor.py:182 (_format)
   */
  _format(mul: string, cr: string): string {
    const plus = ' + ';
    const minus = ' - ';

    const entries = this._entries();
    if (entries.length === 0) {
      return '0';
    }

    const [p0, m0] = entries[0]!;
    let r: string;
    if (m0 === 1n) {
      r = p0.toString();
    } else if (m0 === -1n) {
      r = `- ${p0}`;
    } else {
      r = `${m0}${mul}${p0}`;
    }
    for (let i = 1; i < entries.length; i++) {
      const [p, m] = entries[i]!;
      if (m === 1n) {
        r += cr + plus + p.toString();
      } else if (m === -1n) {
        r += cr + minus + p.toString();
      } else if (m > 0n) {
        r += `${cr}${plus}${m}${mul}${p}`;
      } else if (m < 0n) {
        r += `${cr}${minus}${-m}${mul}${p}`;
      }
    }
    return r;
  }

  /**
   * @see Reference: sage/rings/function_field/divisor.py:232 (_repr_)
   */
  _repr_(split: boolean = true): string {
    return this._format('*', split ? '\n' : '');
  }

  toString(): string {
    return this._repr_();
  }

  /**
   * Compare this divisor with ``other`` lexicographically, viewed as lists of
   * (place, multiplicity) pairs read from the largest place downwards.
   *
   * @see Reference: sage/rings/function_field/divisor.py:270 (_richcmp_)
   */
  cmp(other: FunctionFieldDivisor<C>): number {
    const s = this._entries();
    const o = other._entries();
    let si = s.length - 1;
    let oi = o.length - 1;
    while (si >= 0 && oi >= 0) {
      const sk = s[si]!;
      const ok = o[oi]!;
      const c = sk[0].cmp(ok[0]);
      if (c === 0) {
        if (sk[1] === ok[1]) {
          si--;
          oi--;
          continue;
        }
        return sk[1] < ok[1] ? -1 : 1;
      }
      return c;
    }
    const sl = si + 1;
    const ol = oi + 1;
    return sl < ol ? -1 : sl > ol ? 1 : 0;
  }

  eq(other: FunctionFieldDivisor<C>): boolean {
    return this.cmp(other) === 0;
  }

  /**
   * Return the additive inverse of the divisor.
   *
   * @see Reference: sage/rings/function_field/divisor.py:311 (_neg_)
   */
  neg(): FunctionFieldDivisor<C> {
    return new FunctionFieldDivisor(
      this._parent,
      this._entries().map(([p, m]) => [p, -m] as DivisorEntry<C>)
    );
  }

  /**
   * Add the divisor to the other divisor (or to a place).
   *
   * @see Reference: sage/rings/function_field/divisor.py:338 (_add_)
   */
  add(other: FunctionFieldDivisor<C> | FunctionFieldPlace<C>): FunctionFieldDivisor<C> {
    const rhs = other instanceof FunctionFieldDivisor ? other : other.divisor();
    const data = new Map<string, DivisorEntry<C>>();
    for (const [key, [p, m]] of this._data) {
      data.set(key, [p, m]);
    }
    for (const [key, [p, m]] of rhs._data) {
      const cur = data.get(key);
      const sum = (cur ? cur[1] : 0n) + m;
      if (sum === 0n) {
        data.delete(key);
      } else {
        data.set(key, [p, sum]);
      }
    }
    return new FunctionFieldDivisor(this._parent, data.values());
  }

  /** Return ``this - other``. */
  sub(other: FunctionFieldDivisor<C> | FunctionFieldPlace<C>): FunctionFieldDivisor<C> {
    const rhs = other instanceof FunctionFieldDivisor ? other : other.divisor();
    return this.add(rhs.neg());
  }

  /**
   * Multiply the integer ``i`` into the divisor.
   *
   * @see Reference: sage/rings/function_field/divisor.py:366 (_rmul_)
   */
  scalar_mul(i: bigint | number): FunctionFieldDivisor<C> {
    const k = BigInt(i);
    return new FunctionFieldDivisor(
      this._parent,
      this._entries()
        .map(([p, m]) => [p, k * m] as DivisorEntry<C>)
        .filter(([, m]) => m !== 0n)
    );
  }

  /**
   * Return the dictionary representing the divisor.
   *
   * @see Reference: sage/rings/function_field/divisor.py:391 (dict)
   */
  dict(): Map<string, DivisorEntry<C>> {
    return this._data;
  }

  /**
   * Return the list of place and multiplicity pairs of the divisor.
   *
   * @see Reference: sage/rings/function_field/divisor.py:409 (list)
   */
  list(): Array<DivisorEntry<C>> {
    return this._entries();
  }

  /**
   * Return the support of the divisor.
   *
   * @see Reference: sage/rings/function_field/divisor.py:427 (support)
   */
  support(): Array<FunctionFieldPlace<C>> {
    return this._entries().map(([p]) => p);
  }

  /**
   * Return the multiplicity of the divisor at the place.
   *
   * @see Reference: sage/rings/function_field/divisor.py:445 (multiplicity)
   */
  multiplicity(place: FunctionFieldPlace<C>): bigint {
    const e = this._data.get(place._key());
    return e ? e[1] : 0n;
  }

  /** Alias of {@link multiplicity}, as in SageMath. */
  valuation(place: FunctionFieldPlace<C>): bigint {
    return this.multiplicity(place);
  }

  /**
   * Return ``true`` if this divisor has nonnegative multiplicity at all places.
   *
   * @see Reference: sage/rings/function_field/divisor.py:470 (is_effective)
   */
  is_effective(): boolean {
    for (const [, m] of this._data.values()) {
      if (m < 0n) {
        return false;
      }
    }
    return true;
  }

  /**
   * Return the numerator (positive) part of the divisor.
   *
   * @see Reference: sage/rings/function_field/divisor.py:490 (numerator)
   */
  numerator(): FunctionFieldDivisor<C> {
    return new FunctionFieldDivisor(
      this._parent,
      this._entries().filter(([, m]) => m > 0n)
    );
  }

  /**
   * Return the denominator part of the divisor: the negative of its negative
   * part.
   *
   * @see Reference: sage/rings/function_field/divisor.py:514 (denominator)
   */
  denominator(): FunctionFieldDivisor<C> {
    return new FunctionFieldDivisor(
      this._parent,
      this._entries()
        .filter(([, m]) => m < 0n)
        .map(([p, m]) => [p, -m] as DivisorEntry<C>)
    );
  }

  /**
   * Return the degree of the divisor.
   *
   * @see Reference: sage/rings/function_field/divisor.py:539 (degree)
   */
  degree(): bigint {
    let s = 0n;
    for (const [p, m] of this._entries()) {
      s += p.degree() * m;
    }
    return s;
  }

  /**
   * Return the dimension of the Riemann-Roch space of the divisor.
   *
   * @see Reference: sage/rings/function_field/divisor.py:554 (dimension)
   */
  dimension(): bigint {
    return BigInt(this.basis_function_space().length);
  }

  /**
   * Return a basis of the Riemann-Roch space of the divisor.
   *
   * @see Reference: sage/rings/function_field/divisor.py:572 (basis_function_space)
   */
  basis_function_space(): Array<FunctionFieldElement_rational<C>> {
    return this._function_space()[0];
  }

  /**
   * Return an (echelon) basis and coordinates function for the Riemann-Roch
   * space of the divisor.
   *
   * @see Reference: sage/rings/function_field/divisor.py:632 (_function_space)
   */
  _function_space(): [
    Array<FunctionFieldElement_rational<C>>,
    (f: FunctionFieldElement_rational<C>) => C[],
  ] {
    if (this._functionSpaceCache === null) {
      this._functionSpaceCache = this._echelon_basis(this._basis());
    }
    return this._functionSpaceCache;
  }

  /**
   * Return the vector space of the Riemann-Roch space of the divisor.
   *
   * SageMath returns a `VectorSpace` object plus two morphisms; we return the
   * dimension plus the two maps, because this port has no `VectorSpace` type.
   *
   * @see Reference: sage/rings/function_field/divisor.py:589 (function_space)
   * @see Deviation: function-field Riemann-Roch space returned as maps
   */
  function_space(): [
    number,
    (v: C[]) => FunctionFieldElement_rational<C>,
    (f: FunctionFieldElement_rational<C>) => C[],
  ] {
    const F = this._parent._field;
    const [basis, coordinates] = this._function_space();
    const n = basis.length;

    const from_V = (v: C[]): FunctionFieldElement_rational<C> => {
      let s = F.zero() as FunctionFieldElement_rational<C>;
      for (let i = 0; i < n; i++) {
        s = s.add(basis[i]!.scalar_mul(v[i]!)) as FunctionFieldElement_rational<C>;
      }
      return s;
    };
    return [n, from_V, coordinates];
  }

  /**
   * Return an (echelon) basis and coordinates function for the differential
   * space `\Omega(D)` of the divisor.
   *
   * @see Reference: sage/rings/function_field/divisor.py:736 (_differential_space)
   */
  _differential_space(): [
    Array<FunctionFieldElement_rational<C>>,
    (f: FunctionFieldElement_rational<C>) => C[],
  ] {
    if (this._differentialSpaceCache === null) {
      const F = this._parent._field;
      const x = F.base_field().gen();
      const d = x.divisor_of_poles().scalar_mul(-2n).add(F.different()).sub(this);
      this._differentialSpaceCache = this._echelon_basis(d._basis());
    }
    return this._differentialSpaceCache;
  }

  /**
   * Return a basis of the space of differentials `\Omega(D)`.
   *
   * @see Reference: sage/rings/function_field/divisor.py:664 (basis_differential_space)
   */
  basis_differential_space(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: FunctionFieldDivisor.basis_differential_space ' +
        '(needs sage/rings/function_field/differential.py); use _differential_space() ' +
        'for the underlying function basis'
    );
  }

  /**
   * Return a basis of the Riemann-Roch space of the divisor.
   *
   * This implements Hess' algorithm 6.1 in [Hes2002]_.  SageMath's version is
   * written for a function field of arbitrary degree `n` over its rational
   * base field; this port only has rational function fields, where `n = 1` and
   * the weak-Popov transformation of step 3 is vacuous (a `1 x 1` matrix has no
   * pivot conflicts).  Steps 1, 2, 2.5 and 4 are carried out exactly as
   * upstream.
   *
   * @see Reference: sage/rings/function_field/divisor.py:771 (_basis)
   */
  _basis(): Array<FunctionFieldElement_rational<C>> {
    const F = this._parent._field;
    const n = Number(F.degree());
    if (n !== 1) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: FunctionFieldDivisor._basis for function fields of degree > 1 ' +
          '(needs the weak Popov form step of Hess algorithm 6.1)'
      );
    }
    const O = F.maximal_order();
    const Oinf = F.maximal_order_infinite();

    // Step 1
    let I = O.ideal(F.one());
    let J = Oinf.ideal(F.one());
    for (const [p, m] of this._entries()) {
      if (p.is_infinite_place()) {
        J = J.mul(p.prime_ideal().pow(-m));
      } else {
        I = I.mul(p.prime_ideal().pow(-m));
      }
    }

    // Step 2: M = C * B^{-1} is the 1x1 matrix [i/j].
    const i = I.gens_over_base()[0] as FunctionFieldElement_rational<C>;
    const j = J.gens_over_base()[0] as FunctionFieldElement_rational<C>;
    const M00 = i.div(j);

    // Step 2.5
    const den = M00.denominator();
    const num = M00.numerator();
    const ideg = num.degree();

    // Step 4
    const x = F.gen();
    const basis: Array<FunctionFieldElement_rational<C>> = [];
    let xk = F.one();
    for (let k = 0; k <= den.degree() - ideg; k++) {
      basis.push(xk.mul(i) as FunctionFieldElement_rational<C>);
      xk = xk.mul(x);
    }
    return basis;
  }

  /**
   * Compute an echelonized basis of the subspace generated by ``basis`` over
   * the constant field `k`, together with a coordinates function.
   *
   * @see Reference: sage/rings/function_field/divisor.py:882 (_echelon_basis)
   */
  _echelon_basis(
    basis: Array<FunctionFieldElement_rational<C>>
  ): [Array<FunctionFieldElement_rational<C>>, (f: FunctionFieldElement_rational<C>) => C[]] {
    const F = this._parent._field;
    const k = F.constant_base_field();
    const m = basis.length;

    // The rational function field is one-dimensional over itself, so the
    // "vectors" of upstream's ``free_module`` are the elements themselves.
    const vbasis: Array<FunctionFieldElement_rational<C>> = [...basis];

    const pivot = (v: FunctionFieldElement_rational<C>): [number, number] | null => {
      if (v.is_zero()) {
        return null;
      }
      return [0, v.numerator().degree() - v.denominator().degree()];
    };
    const greater = (v: [number, number], w: [number, number]): boolean =>
      v[0] < w[0] || (v[0] === w[0] && v[1] > w[1]);
    const pkey = (p: [number, number]): string => `${p[0]},${p[1]}`;

    // collate rows by their pivot position
    const pivot_rows = new Map<string, [[number, number], number[]]>();
    const addRow = (p: [number, number], idx: number): void => {
      const key = pkey(p);
      const cur = pivot_rows.get(key);
      if (cur) {
        cur[1].push(idx);
      } else {
        pivot_rows.set(key, [p, [idx]]);
      }
    };
    for (let idx = 0; idx < m; idx++) {
      const p = pivot(vbasis[idx]!);
      if (p !== null) {
        addRow(p, idx);
      }
    }

    // leading coefficient of the (single) component, as in upstream:
    //   c = v.numerator().lc() / v.denominator().lc()
    const lead = (v: FunctionFieldElement_rational<C>): C =>
      divide_constants(v.numerator().leading_coefficient(), v.denominator().leading_coefficient());

    const nbasis: Array<FunctionFieldElement_rational<C>> = [];
    const npivots: Array<[number, number]> = [];
    while (pivot_rows.size > 0) {
      const pivots = [...pivot_rows.values()].map(([p]) => p);
      let head = pivots[0]!;
      for (const p of pivots.slice(1)) {
        if (!greater(head, p)) {
          head = p;
        }
      }
      const rows = pivot_rows.get(pkey(head))![1];
      if (rows.length > 1) {
        const r = rows[0]!;
        const cr = lead(vbasis[r]!);
        for (const idx of rows.slice(1)) {
          const ci = lead(vbasis[idx]!);
          const factor = divide_constants(ci, cr);
          vbasis[idx] = vbasis[idx]!.sub(
            vbasis[r]!.scalar_mul(factor)
          ) as FunctionFieldElement_rational<C>;
          const p = pivot(vbasis[idx]!);
          if (p !== null) {
            addRow(p, idx);
          }
        }
      }
      nbasis.push(vbasis[rows[0]!]!);
      npivots.push(head);
      pivot_rows.delete(pkey(head));
    }

    const coordinates = (f: FunctionFieldElement_rational<C>): C[] => {
      let v = f;
      const coords: C[] = [];
      for (let idx = 0; idx < m; idx++) {
        coords.push(k.zero());
      }
      while (!v.is_zero()) {
        const p = pivot(v)!;
        const ind = npivots.findIndex((q) => q[0] === p[0] && q[1] === p[1]);
        if (ind < 0) {
          throw new ValueError('element is not in the Riemann-Roch space');
        }
        const w = nbasis[ind]!;
        const c = divide_constants(lead(v), lead(w));
        v = v.sub(w.scalar_mul(c)) as FunctionFieldElement_rational<C>;
        coords[ind] = c;
      }
      return coords;
    };

    return [nbasis, coordinates];
  }
}

/**
 * Groups of divisors of function fields.
 *
 * @see Reference: sage/rings/function_field/divisor.py:985 (DivisorGroup)
 */
export class DivisorGroup<C extends ConstantFieldElement> {
  readonly _field: FunctionField<C>;

  constructor(field: FunctionField<C>) {
    this._field = field;
  }

  /**
   * @see Reference: sage/rings/function_field/divisor.py:1017 (_repr_)
   */
  _repr_(): string {
    return `Divisor group of ${this._field}`;
  }

  toString(): string {
    return this._repr_();
  }

  /**
   * Construct a divisor from ``x``.
   *
   * @see Reference: sage/rings/function_field/divisor.py:1030 (_element_constructor_)
   */
  __call__(x: unknown): FunctionFieldDivisor<C> {
    if (x === 0 || x === 0n) {
      return new FunctionFieldDivisor(this, []);
    }
    throw new ValueError(`cannot construct a divisor from ${x}`);
  }

  /** Return the zero divisor. */
  zero(): FunctionFieldDivisor<C> {
    return new FunctionFieldDivisor(this, []);
  }

  /**
   * Return the function field to which the divisor group is attached.
   *
   * @see Reference: sage/rings/function_field/divisor.py:1068 (function_field)
   */
  function_field(): FunctionField<C> {
    return this._field;
  }
}
