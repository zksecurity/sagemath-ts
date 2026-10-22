import { ValueError } from '../../../sagemath-ts/src/errors.js';
import type {
  FiniteFieldElement,
  FiniteFieldExtension,
} from '../../../sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
import { GF } from '../../../sagemath-ts/src/rings/finite_rings/index.js';
import type { EllipticCurvePoint as GenericPoint } from '../../../sagemath-ts/src/schemes/elliptic_curves/ell_point.js';
import {
  EllipticCurve,
  type EllipticCurveGeneric,
  type EllipticCurvePoint,
  ate_pairing,
} from '../../../sagemath-ts/src/schemes/elliptic_curves/index.js';
import { namedCurve } from './named.js';
import { bls12_381, bn254 } from './params.js';
import type { WeierstrassCurve } from './weierstrass.js';

export type Fp2Point = GenericPoint<FiniteFieldElement>;
export type Fp12Point = GenericPoint<FiniteFieldElement>;

export interface G2Curve {
  name: string;
  model: 'short-weierstrass';
  /** Fp[u]/(u² + 1). */
  field: FiniteFieldExtension;
  curve: EllipticCurveGeneric<FiniteFieldElement>;
  generator: Fp2Point;
  order: bigint;
  cofactor: bigint;
  /** `D`: `y² = x³ + b/ξ`, untwisted by `(x·w², y·w³)`. `M`: `y² = x³ + b·ξ`, by `(x/w², y/w³)`. */
  twist: 'D' | 'M';
}

export interface PairingCurve {
  name: string;
  g1: WeierstrassCurve;
  g2: G2Curve;
  Fp2: FiniteFieldExtension;
  /** Fp[w]/(w¹² - 2c·w⁶ + c² + 1), so that w⁶ = ξ = c + u. */
  Fp12: FiniteFieldExtension;
  xi: FiniteFieldElement;
  /** `y² = x³ + b` over Fp12; `liftG1` and `untwist` land here. */
  E12: EllipticCurveGeneric<FiniteFieldElement>;
  embeddingDegree: 12;
  /** Trace of Frobenius of E(Fp), `p + 1 - #E(Fp)`, as `ate_pairing` takes it. */
  trace: bigint;
  fp2(c0: bigint, c1: bigint): FiniteFieldElement;
  /** A G2-curve point from `[c0, c1]` coordinates meaning `c0 + c1·u`. */
  g2Point(x: [bigint, bigint], y: [bigint, bigint]): Fp2Point;
  /** The field embedding `u ↦ w⁶ - c`. */
  fp2ToFp12(a: FiniteFieldElement): FiniteFieldElement;
  liftG1(P: EllipticCurvePoint): Fp12Point;
  untwist(Q: Fp2Point): Fp12Point;
  /**
   * Sage's reduced ate pairing `ate_pairing(liftG1(P), untwist(Q), r, 12, trace)`.
   * Bilinear and non-degenerate on G1 × G2, but not the optimal ate pairing of other
   * libraries: compare relations between values, not values.
   */
  pairing(P: EllipticCurvePoint, Q: Fp2Point): FiniteFieldElement;
}

export type PairingCurveName = 'bn254' | 'bls12_381';

const specs = {
  bn254: { params: bn254, c: 9n, twist: 'D' },
  bls12_381: { params: bls12_381, c: 1n, twist: 'M' },
} as const;

function createPairingCurve(name: PairingCurveName): PairingCurve {
  const { params, c, twist } = specs[name];
  const g1 = namedCurve(name);
  const p = params.field.p;
  const r = params.r;

  const Fp2 = GF(p ** 2n, { name: 'u', modulus: [1n, 0n, 1n] }) as FiniteFieldExtension;
  const fp2 = (c0: bigint, c1: bigint) => Fp2.__call__([c0, c1]);
  const xi = fp2(c, 1n);
  const b = Fp2.__call__(params.b);
  const E2 = EllipticCurve<FiniteFieldElement>(Fp2, [
    0n,
    0n,
    0n,
    0n,
    twist === 'D' ? b.div(xi) : b.mul(xi),
  ]);
  const g2Point = (x: [bigint, bigint], y: [bigint, bigint]) => E2.point([fp2(...x), fp2(...y)]);

  const Fp12 = GF(p ** 12n, {
    name: 'w',
    modulus: [c * c + 1n, 0n, 0n, 0n, 0n, 0n, p - 2n * c, 0n, 0n, 0n, 0n, 0n, 1n],
  }) as FiniteFieldExtension;
  const E12 = EllipticCurve<FiniteFieldElement>(Fp12, [0n, 0n, 0n, 0n, params.b]);
  const w = Fp12.gen();
  const [wx, wy] = twist === 'D' ? [w.pow(2n), w.pow(3n)] : [w.pow(2n).inv(), w.pow(3n).inv()];

  const fp2ToFp12 = (a: FiniteFieldElement) => {
    const [c0, c1] = Fp2.__call__(a).coefficients();
    const a0 = c0?.value ?? 0n;
    const a1 = c1?.value ?? 0n;
    return Fp12.__call__([a0 - c * a1 + c * p, 0n, 0n, 0n, 0n, 0n, a1]);
  };

  const g1Curve = g1.curve;
  const liftG1 = (P: EllipticCurvePoint) => {
    const E = P.curve;
    if (
      E !== g1Curve &&
      (E.field.characteristic !== p || E.a.value !== 0n || E.b.value !== params.b)
    ) {
      throw new ValueError(`P is not a point of ${g1.name}`);
    }
    if (P.isZero()) return E12.zero();
    return E12.point([Fp12.__call__(P.x!.value), Fp12.__call__(P.y!.value)]);
  };
  const untwist = (Q: Fp2Point) => {
    if (Q.curve !== E2)
      throw new ValueError(`Q is not a point of ${name}'s g2.curve; build it with g2Point`);
    if (Q.is_zero()) return E12.zero();
    return E12.point([fp2ToFp12(Q.x()).mul(wx), fp2ToFp12(Q.y()).mul(wy)]);
  };
  const trace = p + 1n - g1.cofactor * g1.order;
  const g2Cofactor = params.g2.cofactor!;
  const g2Generator = g2Point(params.g2.x, params.g2.y);
  g2Generator.setOrder(r, false);
  E2._order = g2Cofactor * r;

  return {
    name: params.name,
    g1,
    g2: {
      name: `${params.name}-g2`,
      model: 'short-weierstrass',
      field: Fp2,
      curve: E2,
      generator: g2Generator,
      order: r,
      cofactor: g2Cofactor,
      twist,
    },
    Fp2,
    Fp12,
    xi,
    E12,
    embeddingDegree: 12,
    trace,
    fp2,
    g2Point,
    fp2ToFp12,
    liftG1,
    untwist,
    pairing: (P, Q) => ate_pairing(liftG1(P), untwist(Q), r, 12, trace),
  };
}

const cache = new Map<PairingCurveName, PairingCurve>();

/** BN254 or BLS12-381 with G1, G2 and the Fp12 curve both map into for pairings. */
export function pairingCurve(name: PairingCurveName): PairingCurve {
  if (!Object.hasOwn(specs, name)) {
    throw new ValueError(`unknown pairing curve ${name}; expected bn254 or bls12_381`);
  }
  let curve = cache.get(name);
  if (curve === undefined) {
    curve = createPairingCurve(name);
    cache.set(name, curve);
  }
  return curve;
}
