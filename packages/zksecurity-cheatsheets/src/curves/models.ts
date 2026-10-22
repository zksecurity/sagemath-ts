import { inverse_mod } from '../../../sagemath-ts/src/arith/misc.js';
import { ValueError } from '../../../sagemath-ts/src/errors.js';
import type { EllipticCurvePoint } from '../../../sagemath-ts/src/schemes/elliptic_curves/index.js';
import type { MontgomeryCurveParams, TwistedEdwardsCurveParams } from './params.js';
import { type PrimeFieldCurve, createPrimeFieldCurve } from './weierstrass.js';

export interface MontgomeryPoint {
  u: bigint;
  v: bigint;
}

export interface EdwardsPoint {
  x: bigint;
  y: bigint;
}

/** `B v² = u³ + A u² + u`. `null` stands for the point at infinity. */
export interface MontgomeryCurve extends PrimeFieldCurve {
  model: 'montgomery';
  A: bigint;
  B: bigint;
  isOnCurve(P: MontgomeryPoint): boolean;
  toWeierstrass(P: MontgomeryPoint | null): EllipticCurvePoint;
  fromWeierstrass(P: EllipticCurvePoint): MontgomeryPoint | null;
  /** A point with u-coordinate `u mod p`, or `null` when u lies on the quadratic twist. */
  liftU(u: bigint): EllipticCurvePoint | null;
}

/** `a x² + y² = 1 + d x² y²`. Its Montgomery form shares `curve`. */
export interface TwistedEdwardsCurve extends PrimeFieldCurve {
  model: 'twisted-edwards';
  a: bigint;
  d: bigint;
  montgomery: MontgomeryCurve;
  isOnCurve(P: EdwardsPoint): boolean;
  toWeierstrass(P: EdwardsPoint): EllipticCurvePoint;
  fromWeierstrass(P: EllipticCurvePoint): EdwardsPoint;
  toMontgomery(P: EdwardsPoint): MontgomeryPoint | null;
  fromMontgomery(P: MontgomeryPoint | null): EdwardsPoint;
}

function modArith(p: bigint, name: string) {
  if (p <= 3n) throw new ValueError(`${name}: the characteristic must exceed 3`);
  const m = (a: bigint) => ((a % p) + p) % p;
  const div = (a: bigint, b: bigint) => m(a * inverse_mod(m(b), p));
  return { m, div };
}

function sameCurve(P: EllipticCurvePoint, base: PrimeFieldCurve): boolean {
  const E = P.curve;
  return (
    E === base.curve ||
    (E.field.characteristic === base.field.characteristic &&
      E.a.value === base.curve.a.value &&
      E.b.value === base.curve.b.value)
  );
}

/**
 * Montgomery parameters with their short-Weierstrass model
 * `y² = x³ + (3 - A²)/(3B²) x + (2A³ - 9A)/(27B³)` under `x = (3u + A)/(3B)`, `y = v/B`.
 * `params.order` is `#E`; the result's `order` is the generator's, `#E / cofactor`.
 */
export function createMontgomeryCurve(params: MontgomeryCurveParams): MontgomeryCurve {
  const p = params.field.p;
  const { m, div } = modArith(p, params.name);
  const A = m(params.A);
  const B = m(params.B);
  if (B === 0n || m(A * A - 4n) === 0n) {
    throw new ValueError(`${params.name}: singular Montgomery curve`);
  }
  if (params.order % params.cofactor !== 0n) {
    throw new ValueError(`${params.name}: the cofactor does not divide the curve order`);
  }

  const isOnCurve = ({ u, v }: MontgomeryPoint) => m(B * v * v - u * (u * u + A * u + 1n)) === 0n;
  const check = (P: MontgomeryPoint) => {
    if (!isOnCurve(P)) throw new ValueError(`(${P.u}, ${P.v}) is not on ${params.name}`);
  };
  const toXY = ({ u, v }: MontgomeryPoint) => ({ x: div(3n * u + A, 3n * B), y: div(v, B) });

  check(params.generator);
  const base = createPrimeFieldCurve({
    name: params.name,
    field: { p },
    a: div(3n - A * A, 3n * B * B),
    b: div(2n * A ** 3n - 9n * A, 27n * B ** 3n),
    order: params.order / params.cofactor,
    cofactor: params.cofactor,
    generator: toXY(params.generator),
  });

  const toWeierstrass = (P: MontgomeryPoint | null): EllipticCurvePoint => {
    if (P === null) return base.curve.zero();
    const Q = { u: m(P.u), v: m(P.v) };
    check(Q);
    const { x, y } = toXY(Q);
    return base.curve.point(x, y);
  };

  return {
    ...base,
    model: 'montgomery',
    A,
    B,
    isOnCurve: (P) => isOnCurve({ u: m(P.u), v: m(P.v) }),
    toWeierstrass,
    fromWeierstrass(P) {
      if (!sameCurve(P, base)) throw new ValueError(`P is not a point of ${params.name}'s curve`);
      if (P.isZero()) return null;
      return { u: m(B * P.x!.value - div(A, 3n)), v: m(B * P.y!.value) };
    },
    liftU(u) {
      const U = m(u);
      const rhs = base.field.__call__(div(U * (U * U + A * U + 1n), B));
      if (!rhs.is_square()) return null;
      return toWeierstrass({ u: U, v: rhs.sqrt({ extend: false }).value });
    },
  };
}

/**
 * Twisted Edwards parameters with Montgomery and short-Weierstrass models.
 *
 * The Montgomery form is `A = 2(a + d)/(a - d)`, `B = 4/(a - d)` under
 * `u = (1 + y)/(1 - y)`, `v = u/x`. Passing `montgomery` targets that curve instead
 * (for ed25519, curve25519): `v = s·u/x` with `s` fixed by the two generators.
 *
 * Exceptional points: `(0, 1)` is the identity, `(0, -1)` maps to `(0, 0)`. The other
 * 2-torsion points `(u, 0)` exist when `d/a` is a square, and `u = -1` when `d` is a
 * square; both sit at infinity on the Edwards model, so `fromMontgomery` throws.
 */
export function createTwistedEdwardsCurve(
  params: TwistedEdwardsCurveParams,
  montgomery?: MontgomeryCurve
): TwistedEdwardsCurve {
  const p = params.field.p;
  const { m, div } = modArith(p, params.name);
  const a = m(params.a);
  const d = m(params.d);
  if (a === 0n || d === 0n || a === d) {
    throw new ValueError(`${params.name}: singular twisted Edwards curve`);
  }
  if (params.order !== params.cofactor * params.subgroupOrder) {
    throw new ValueError(`${params.name}: order must equal cofactor * subgroupOrder`);
  }

  const isOnCurve = ({ x, y }: EdwardsPoint) =>
    m(a * x * x + y * y - 1n - d * x * x * y * y) === 0n;
  const check = (P: EdwardsPoint) => {
    if (!isOnCurve(P)) throw new ValueError(`(${P.x}, ${P.y}) is not on ${params.name}`);
  };
  const toMontgomeryScaled = ({ x, y }: EdwardsPoint, s: bigint): MontgomeryPoint | null => {
    if (x === 0n) return y === 1n ? null : { u: 0n, v: 0n };
    const u = div(1n + y, 1n - y);
    return { u, v: div(s * u, x) };
  };

  const G = { x: m(params.generator.x), y: m(params.generator.y) };
  check(G);
  const canonical = toMontgomeryScaled(G, 1n);
  if (canonical === null) throw new ValueError(`${params.name}: the generator is the identity`);
  const A = div(2n * (a + d), a - d);
  const B = div(4n, a - d);

  let s = 1n;
  let mont: MontgomeryCurve;
  if (montgomery === undefined) {
    mont = createMontgomeryCurve({
      name: params.name,
      field: { p },
      A,
      B,
      order: params.order,
      cofactor: params.cofactor,
      generator: canonical,
    });
  } else {
    const target = montgomery.fromWeierstrass(montgomery.generator);
    if (
      montgomery.field.characteristic !== p ||
      montgomery.A !== A ||
      montgomery.order !== params.subgroupOrder ||
      target === null ||
      target.u !== canonical.u
    ) {
      throw new ValueError(`${params.name} does not correspond to ${montgomery.name}`);
    }
    s = div(target.v, canonical.v);
    if (m(s * s * montgomery.B) !== B) {
      throw new ValueError(`${params.name} does not correspond to ${montgomery.name}`);
    }
    mont = montgomery;
  }

  const toMontgomery = (P: EdwardsPoint): MontgomeryPoint | null => {
    const Q = { x: m(P.x), y: m(P.y) };
    check(Q);
    return toMontgomeryScaled(Q, s);
  };
  const fromMontgomery = (P: MontgomeryPoint | null): EdwardsPoint => {
    if (P === null) return { x: 0n, y: 1n };
    const u = m(P.u);
    const v = m(P.v);
    if (!mont.isOnCurve({ u, v })) throw new ValueError(`(${u}, ${v}) is not on ${mont.name}`);
    if (v === 0n) {
      if (u === 0n) return { x: 0n, y: p - 1n };
      throw new ValueError(`(${u}, 0) has order 2 and lies at infinity on ${params.name}`);
    }
    if (u === p - 1n) {
      throw new ValueError(`(${u}, ${v}) has u = -1 and lies at infinity on ${params.name}`);
    }
    return { x: div(s * u, v), y: div(u - 1n, u + 1n) };
  };

  return {
    name: params.name,
    model: 'twisted-edwards',
    field: mont.field,
    curve: mont.curve,
    generator: mont.toWeierstrass(toMontgomeryScaled(G, s)),
    order: params.subgroupOrder,
    cofactor: params.cofactor,
    a,
    d,
    montgomery: mont,
    isOnCurve: (P) => isOnCurve({ x: m(P.x), y: m(P.y) }),
    toWeierstrass: (P) => mont.toWeierstrass(toMontgomery(P)),
    fromWeierstrass: (P) => fromMontgomery(mont.fromWeierstrass(P)),
    toMontgomery,
    fromMontgomery,
  };
}
