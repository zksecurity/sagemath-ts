/**
 * @module sage/schemes/hyperelliptic_curves/invariants
 * @description Compute invariants of quintics and sextics via 'Ueberschiebung'
 *
 * Port of: `sage/schemes/hyperelliptic_curves/invariants.py`
 *
 * AUTHOR (upstream): Nick Alexander
 *
 * ## Representation
 *
 * Upstream works with multivariate polynomials and expands the differential
 * operator `(f g)_k` symbolically in `QQ[dfdx, dfdy, dgdx, dgdy]`
 * (`invariants.py:44-138`).  Every object that appears in `ubs` is a
 * *homogeneous binary form*, and transvectants of homogeneous forms are again
 * homogeneous, so this port represents a form of degree `d` densely as
 *
 * ```
 * F = sum_{i=0}^{d} c[i] * x^i * y^(d-i)
 * ```
 *
 * and expands `(fx gy - fy gx)^k` with the binomial theorem, which is exactly
 * what `differential_operator` + `diffsymb` compute:
 *
 * ```
 * (f g)_k = const * sum_{i=0}^{k} (-1)^i C(k,i)
 *              (d^k f / dx^{k-i} dy^i) (d^k g / dx^i dy^{k-i})
 * ```
 *
 * with `const = (m-k)! (n-k)! / (m! n!)`, `n = max(deg f, k)`,
 * `m = max(deg g, k)`.
 */

import { binomial, factorial } from '../../arith/misc.js';
import { NotImplementedError, ValueError } from '../../errors.js';
import type { Polynomial, RingElement } from '../../rings/polynomial/polynomial_element.js';
import { type HyperellipticBaseRing, characteristic_of, div_elements } from './field_ops.js';

/**
 * A homogeneous binary form of degree `d`:
 * `sum_{i=0}^{d} c[i] x^i y^(d-i)`.
 *
 * `d` is the *nominal* degree; a form may be identically zero while keeping a
 * positive nominal degree, exactly as in the symbolic computation upstream.
 */
export interface BinaryForm<C extends RingElement> {
  readonly d: number;
  readonly c: readonly C[];
}

function zero_form<C extends RingElement>(K: HyperellipticBaseRing<C>, d: number): BinaryForm<C> {
  const c: C[] = [];
  for (let i = 0; i <= d; i++) {
    c.push(K.zero() as C);
  }
  return { d, c };
}

/** `n * a` in the base ring, via double-and-add on the element itself. */
function int_times<C extends RingElement>(K: HyperellipticBaseRing<C>, n: bigint, a: C): C {
  return (K.__call__(n) as C).mul(a) as C;
}

/**
 * Differentiate `F` `xtimes` times with respect to `x` and `ytimes` times with
 * respect to `y`.
 *
 * Port of `invariants.py:18-41` (`diffxy`).
 */
export function diffxy<C extends RingElement>(
  K: HyperellipticBaseRing<C>,
  F: BinaryForm<C>,
  xtimes: number,
  ytimes: number
): BinaryForm<C> {
  let H = F;
  for (let i = 0; i < xtimes; i++) {
    H = derivative_x(K, H);
  }
  for (let j = 0; j < ytimes; j++) {
    H = derivative_y(K, H);
  }
  return H;
}

function derivative_x<C extends RingElement>(
  K: HyperellipticBaseRing<C>,
  F: BinaryForm<C>
): BinaryForm<C> {
  if (F.d <= 0) {
    return zero_form(K, 0);
  }
  const c: C[] = [];
  for (let j = 0; j <= F.d - 1; j++) {
    c.push(int_times(K, BigInt(j + 1), F.c[j + 1]!));
  }
  return { d: F.d - 1, c };
}

function derivative_y<C extends RingElement>(
  K: HyperellipticBaseRing<C>,
  F: BinaryForm<C>
): BinaryForm<C> {
  if (F.d <= 0) {
    return zero_form(K, 0);
  }
  const c: C[] = [];
  for (let j = 0; j <= F.d - 1; j++) {
    c.push(int_times(K, BigInt(F.d - j), F.c[j]!));
  }
  return { d: F.d - 1, c };
}

function mul_forms<C extends RingElement>(
  K: HyperellipticBaseRing<C>,
  A: BinaryForm<C>,
  B: BinaryForm<C>
): BinaryForm<C> {
  const d = A.d + B.d;
  const c: C[] = [];
  for (let k = 0; k <= d; k++) {
    let s = K.zero() as C;
    for (let i = Math.max(0, k - B.d); i <= Math.min(A.d, k); i++) {
      s = s.add(A.c[i]!.mul(B.c[k - i]!) as C) as C;
    }
    c.push(s);
  }
  return { d, c };
}

function add_scaled<C extends RingElement>(
  K: HyperellipticBaseRing<C>,
  A: BinaryForm<C>,
  B: BinaryForm<C>,
  factor: bigint
): BinaryForm<C> {
  if (A.d !== B.d) {
    throw new ValueError(`cannot add binary forms of degrees ${A.d} and ${B.d}`);
  }
  const c: C[] = [];
  for (let i = 0; i <= A.d; i++) {
    c.push(A.c[i]!.add(int_times(K, factor, B.c[i]!)) as C);
  }
  return { d: A.d, c };
}

/**
 * Return the transvectant `(f g)_k` of Mestre, p 315 [Mes1991].
 *
 * Port of `invariants.py:113-138` (`Ueberschiebung`), with the symbolic
 * differential operator expanded by the binomial theorem.
 *
 * `differential_operator` (`invariants.py:44-78`) and `diffsymb`
 * (`invariants.py:81-110`) are not exposed separately: they only exist
 * upstream to build and then apply this expansion.
 */
export function Ueberschiebung<C extends RingElement>(
  K: HyperellipticBaseRing<C>,
  f: BinaryForm<C>,
  g: BinaryForm<C>,
  k: number
): BinaryForm<C> {
  const n = Math.max(f.d, k);
  const m = Math.max(g.d, k);

  // const = (m-k)! (n-k)! / (m! n!)
  const num = factorial(BigInt(m - k)) * factorial(BigInt(n - k));
  const den = factorial(BigInt(m)) * factorial(BigInt(n));
  const konst = div_elements(K.__call__(num) as C, K.__call__(den) as C);

  const outDeg = f.d + g.d - 2 * k;
  let res = zero_form(K, Math.max(outDeg, 0));

  for (let i = 0; i <= k; i++) {
    const a = diffxy(K, f, k - i, i);
    const b = diffxy(K, g, i, k - i);
    const prod = mul_forms(K, a, b);
    if (prod.d !== res.d) {
      // Only possible when a derivative bottomed out at degree 0; the product
      // is then identically zero and contributes nothing.
      continue;
    }
    const coeff = binomial(BigInt(k), BigInt(i)) * (i % 2 === 0 ? 1n : -1n);
    res = add_scaled(K, res, prod, coeff);
  }

  const c = res.c.map((x) => x.mul(konst) as C);
  return { d: res.d, c };
}

/** The dictionary of forms computed by `ubs` (`invariants.py:141-205`). */
export interface UbsResult<C extends RingElement> {
  f: BinaryForm<C>;
  i: BinaryForm<C>;
  Delta: BinaryForm<C>;
  y1: BinaryForm<C>;
  y2: BinaryForm<C>;
  y3: BinaryForm<C>;
  A: C;
  B: C;
  C: C;
  D: C;
}

/**
 * Turn a univariate polynomial of degree at most 6 into the corresponding
 * binary sextic `sum_{i=0}^{6} f[i] x^i y^(6-i)` (`invariants.py:190-193`).
 */
export function sextic_form<C extends RingElement>(f: Polynomial<C>): BinaryForm<C> {
  const K = f.parent.base_ring as HyperellipticBaseRing<C>;
  if (f.degree() > 6) {
    throw new ValueError(`the argument must have degree at most 6, got ${f.degree()}`);
  }
  const c: C[] = [];
  for (let i = 0; i <= 6; i++) {
    c.push(f.getCoeff(i));
  }
  void K;
  return { d: 6, c };
}

/**
 * Given a sextic form `f`, return the invariants of Mestre, p 317 [Mes1991].
 *
 * Port of `invariants.py:141-205` (`ubs`).
 */
export function ubs<C extends RingElement>(
  K: HyperellipticBaseRing<C>,
  form: BinaryForm<C>
): UbsResult<C> {
  const ub = (a: BinaryForm<C>, b: BinaryForm<C>, k: number) => Ueberschiebung(K, a, b, k);

  const f = form;
  const i = ub(f, f, 4);
  const Delta = ub(i, i, 2);
  const y1 = ub(f, i, 4);
  const y2 = ub(i, y1, 2);
  const y3 = ub(i, y2, 2);
  const A = ub(f, f, 6);
  const B = ub(i, i, 4);
  const Cc = ub(i, Delta, 4);
  const D = ub(y3, y1, 2);

  const constant = (F: BinaryForm<C>): C => {
    if (F.d !== 0) {
      throw new ValueError(`expected a constant invariant, got degree ${F.d}`);
    }
    return F.c[0]!;
  };

  return {
    f,
    i,
    Delta,
    y1,
    y2,
    y3,
    A: constant(A),
    B: constant(B),
    C: constant(Cc),
    D: constant(D),
  };
}

/**
 * Convert Clebsch invariants `A, B, C, D` to Igusa invariants
 * `I_2, I_4, I_6, I_10`.
 *
 * Port of `invariants.py:208-231`.
 *
 * @param K - the base ring; required because this port's ring elements do not
 *   carry a reference to their parent (see DESIGN.md).
 */
export function clebsch_to_igusa<C extends RingElement>(
  K: HyperellipticBaseRing<C>,
  A: C,
  B: C,
  Cc: C,
  D: C
): [C, C, C, C] {
  const t = (n: bigint, a: C): C => int_times(K, n, a);
  const A2 = A.mul(A) as C;
  const A3 = A2.mul(A) as C;
  const A5 = A3.mul(A2) as C;
  const B2 = B.mul(B) as C;

  const I2 = t(-120n, A);
  const I4 = t(-720n, A2).add(t(6750n, B)) as C;
  const I6 = t(8640n, A3)
    .add(t(-108000n, A.mul(B) as C))
    .add(t(202500n, Cc)) as C;
  const I10 = t(-62208n, A5)
    .add(t(972000n, A3.mul(B) as C))
    .add(t(1620000n, A2.mul(Cc) as C))
    .add(t(-3037500n, A.mul(B2) as C))
    .add(t(-6075000n, B.mul(Cc) as C))
    .add(t(-4556250n, D)) as C;
  return [I2, I4, I6, I10];
}

/**
 * Convert Igusa invariants `I_2, I_4, I_6, I_10` to Clebsch invariants.
 *
 * Port of `invariants.py:234-257`.
 */
export function igusa_to_clebsch<C extends RingElement>(
  K: HyperellipticBaseRing<C>,
  I2: C,
  I4: C,
  I6: C,
  I10: C
): [C, C, C, C] {
  const t = (n: bigint, a: C): C => int_times(K, n, a);
  const over = (a: C, n: bigint): C => div_elements(a, K.__call__(n) as C);

  const I2_2 = I2.mul(I2) as C;
  const I2_3 = I2_2.mul(I2) as C;
  const I2_5 = I2_3.mul(I2_2) as C;

  const A = over(I2.neg() as C, 120n);
  // B = -(-I2^2 - 20*I4)/135000 = (I2^2 + 20*I4)/135000
  const B = over(I2_2.add(t(20n, I4)) as C, 135000n);
  const C1 = over(
    (I2_3.add(t(80n, I2.mul(I4) as C)).add(t(-600n, I6)) as C).neg() as C,
    121500000n
  );
  const D = over(
    (
      t(9n, I2_5)
        .add(t(700n, I2_3.mul(I4) as C))
        .add(t(-3600n, I2_2.mul(I6) as C))
        .add(t(-12400n, I2.mul(I4).mul(I4) as C))
        .add(t(48000n, I4.mul(I6) as C))
        .add(t(10800000n, I10)) as C
    ).neg() as C,
    49207500000000n
  );
  return [A, B, C1, D];
}

/**
 * Given a sextic form `f`, return the Clebsch invariants `(A, B, C, D)` of
 * Mestre, p 317 [Mes1991].
 *
 * Port of `invariants.py:260-290`.
 */
export function clebsch_invariants<C extends RingElement>(f: Polynomial<C>): [C, C, C, C] {
  const K = f.parent.base_ring as HyperellipticBaseRing<C>;
  const ch = characteristic_of(K);
  if (ch === 2n || ch === 3n || ch === 5n) {
    throw new NotImplementedError(
      'Invariants of binary sextics/genus 2 hyperelliptic curves not implemented in characteristics 2, 3, and 5'
    );
  }
  const U = ubs(K, sextic_form(f));
  return [U.A, U.B, U.C, U.D];
}

/**
 * Given a sextic form `f`, return the Igusa-Clebsch invariants
 * `I_2, I_4, I_6, I_10` [IJ1960].
 *
 * Port of `invariants.py:293-331`.
 */
export function igusa_clebsch_invariants<C extends RingElement>(f: Polynomial<C>): [C, C, C, C] {
  const K = f.parent.base_ring as HyperellipticBaseRing<C>;
  const [A, B, Cc, D] = clebsch_invariants(f);
  return clebsch_to_igusa(K, A, B, Cc, D);
}

/**
 * The three absolute Igusa invariants used by van Wamelen [Wam1999].
 *
 * Port of `invariants.py:334-370`.
 */
export function absolute_igusa_invariants_wamelen<C extends RingElement>(
  f: Polynomial<C>
): [C, C, C] {
  const [I2, I4, I6, I10] = igusa_clebsch_invariants(f);
  const I2_2 = I2.mul(I2) as C;
  const I2_3 = I2_2.mul(I2) as C;
  const I2_5 = I2_3.mul(I2_2) as C;
  return [
    div_elements(I2_5, I10),
    div_elements(I2_3.mul(I4) as C, I10),
    div_elements(I2_2.mul(I6) as C, I10),
  ];
}

/**
 * The three absolute Igusa invariants used by Kohel [KohECHIDNA].
 *
 * Port of `invariants.py:373-409`.
 */
export function absolute_igusa_invariants_kohel<C extends RingElement>(
  f: Polynomial<C>
): [C, C, C] {
  const [I2, I4, I6, I10] = igusa_clebsch_invariants(f);
  const I2_2 = I2.mul(I2) as C;
  const I2_3 = I2_2.mul(I2) as C;
  return [
    div_elements(I4.mul(I6) as C, I10),
    div_elements(I2_3.mul(I4) as C, I10),
    div_elements(I2_2.mul(I6) as C, I10),
  ];
}
