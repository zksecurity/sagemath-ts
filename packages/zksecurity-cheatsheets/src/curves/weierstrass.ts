import { is_prime, isqrt } from '../../../sagemath-ts/src/arith/misc.js';
import { ValueError } from '../../../sagemath-ts/src/errors.js';
// `PrimeField` is `GF`'s return type for prime orders but is not re-exported from any
// index, so it has to be imported from the defining module.
import type { PrimeField } from '../../../sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
import { GF } from '../../../sagemath-ts/src/rings/finite_rings/index.js';
import {
  EllipticCurve,
  type EllipticCurveFiniteField,
  type EllipticCurvePoint,
} from '../../../sagemath-ts/src/schemes/elliptic_curves/index.js';
import type { PrimeFieldCurveParams } from './params.js';

export type CurveModel = 'short-weierstrass' | 'montgomery' | 'twisted-edwards';

/**
 * A curve over a prime field, as a short-Weierstrass sagemath-ts curve. `model` is the
 * form the curve is specified in; Montgomery and Edwards curves carry maps to it.
 */
export interface PrimeFieldCurve {
  name: string;
  model: CurveModel;
  field: PrimeField;
  curve: EllipticCurveFiniteField;
  /** Generates the subgroup of prime order `order`. */
  generator: EllipticCurvePoint;
  order: bigint;
  /** `#E(Fp) = cofactor * order`. */
  cofactor: bigint;
}

export interface WeierstrassCurve extends PrimeFieldCurve {
  model: 'short-weierstrass';
}

/**
 * The cofactor `h` when `order` alone fixes it: if `order > 4√p`, exactly one multiple
 * of it lies in the Hasse interval `[p + 1 - 2√p, p + 1 + 2√p]`. Returns `undefined`
 * otherwise. The caller must know `order` is prime and divides `#E`.
 */
export function hasseCofactor(p: bigint, order: bigint): bigint | undefined {
  const s = isqrt(4n * p) + 1n;
  if (order <= 2n * s) return undefined;
  const h = (p + 1n + s) / order;
  return h * order >= p + 1n - s ? h : undefined;
}

/**
 * Build the field, curve and generator for short-Weierstrass parameters.
 *
 * A missing `cofactor` is derived with `hasseCofactor`; when the order does not fix it,
 * this throws rather than count points. When `order` is prime, kills the generator and
 * fixes the cofactor, it is seeded into the point and curve order caches, so `.order()`
 * and `.cardinality()` return it instead of running SEA.
 */
export function createPrimeFieldCurve(params: PrimeFieldCurveParams): WeierstrassCurve {
  // `GF` is `GFExtended`: its second parameter is a generator *name*, not an options
  // object, and it returns `PrimeField | FiniteFieldExtension`. Every curve here is over
  // a prime field, so narrow to that branch for the constructors below.
  const field = GF(params.field.p) as PrimeField;
  // The port carries two prime-field classes: `PrimeField` (finite_field_extension.ts),
  // which is what `GF` returns, and `FiniteFieldPrime` (finite_field_constructor.ts),
  // which is what `EllipticCurve` declares. They are runtime-compatible for everything
  // used here — `FiniteFieldPrime` merely declares four extra members — but not
  // assignable, so the cast stands in until the two classes are reconciled upstream.
  const curve = EllipticCurve(field as unknown as Parameters<typeof EllipticCurve>[0], [
    params.a,
    params.b,
  ]);
  const generator = curve.point(params.generator.x, params.generator.y);
  const cofactor = params.cofactor ?? hasseCofactor(params.field.p, params.order);
  if (cofactor === undefined) {
    throw new ValueError(`${params.name}: the order does not determine the cofactor; pass it`);
  }
  // Parameters under audit may be wrong, so seed the caches only when the values are proven.
  if (is_prime(params.order) && generator.mul(params.order).isZero()) {
    generator.setOrder(params.order, false);
    if (hasseCofactor(params.field.p, params.order) === cofactor) {
      curve.set_order(cofactor * params.order, false);
    }
  }

  return {
    name: params.name,
    model: 'short-weierstrass',
    field,
    curve,
    generator,
    order: params.order,
    cofactor,
  };
}
