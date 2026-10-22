import { ValueError } from '../../../sagemath-ts/src/errors.js';
import {
  type MontgomeryCurve,
  type TwistedEdwardsCurve,
  createMontgomeryCurve,
  createTwistedEdwardsCurve,
} from './models.js';
import {
  babyJubJub,
  bls12_377_g1,
  bls12_381_g1,
  bn254_g1,
  curve25519,
  ed25519,
  grumpkin,
  jubjub,
  p256,
  p384,
  p521,
  pallas,
  secp256k1,
  vesta,
} from './params.js';
import { type WeierstrassCurve, createPrimeFieldCurve } from './weierstrass.js';

export interface NamedCurves {
  secp256k1: WeierstrassCurve;
  p256: WeierstrassCurve;
  p384: WeierstrassCurve;
  p521: WeierstrassCurve;
  /** G1 of BN254; `pairingCurve('bn254')` adds G2. */
  bn254: WeierstrassCurve;
  /** G1 of BLS12-381; `pairingCurve('bls12_381')` adds G2. */
  bls12_381: WeierstrassCurve;
  /** G1 of BLS12-377. */
  bls12_377: WeierstrassCurve;
  grumpkin: WeierstrassCurve;
  pallas: WeierstrassCurve;
  vesta: WeierstrassCurve;
  curve25519: MontgomeryCurve;
  /** Shares `curve` with curve25519 (RFC 7748 §4.1 map). */
  ed25519: TwistedEdwardsCurve;
  babyJubJub: TwistedEdwardsCurve;
  jubjub: TwistedEdwardsCurve;
}

export type NamedCurveName = keyof NamedCurves;
export type NamedCurve = NamedCurves[NamedCurveName];

const builders: { [N in NamedCurveName]: () => NamedCurves[N] } = {
  secp256k1: () => createPrimeFieldCurve(secp256k1),
  p256: () => createPrimeFieldCurve(p256),
  p384: () => createPrimeFieldCurve(p384),
  p521: () => createPrimeFieldCurve(p521),
  bn254: () => createPrimeFieldCurve(bn254_g1),
  bls12_381: () => createPrimeFieldCurve(bls12_381_g1),
  bls12_377: () => createPrimeFieldCurve(bls12_377_g1),
  grumpkin: () => createPrimeFieldCurve(grumpkin),
  pallas: () => createPrimeFieldCurve(pallas),
  vesta: () => createPrimeFieldCurve(vesta),
  curve25519: () => createMontgomeryCurve(curve25519),
  ed25519: () => createTwistedEdwardsCurve(ed25519, namedCurve('curve25519')),
  babyJubJub: () => createTwistedEdwardsCurve(babyJubJub),
  jubjub: () => createTwistedEdwardsCurve(jubjub),
};

export const namedCurveNames = Object.keys(builders) as NamedCurveName[];

const cache = new Map<NamedCurveName, NamedCurve>();

/** A ready-to-use named curve, built once. Orders and cofactors are hardcoded. */
export function namedCurve<N extends NamedCurveName>(name: N): NamedCurves[N] {
  if (!Object.hasOwn(builders, name)) {
    throw new ValueError(`unknown curve ${name}; expected one of ${namedCurveNames.join(', ')}`);
  }
  let curve = cache.get(name);
  if (curve === undefined) {
    curve = builders[name]();
    cache.set(name, curve);
  }
  return curve as NamedCurves[N];
}
