/**
 * Fast affine short-Weierstrass arithmetic over a prime field, plus ECDSA signing
 * and key verification. Used by the ECDSA nonce attacks to plant weak signatures
 * in tests and to check recovered private keys against the public key.
 *
 * This is plain bigint arithmetic: not constant time, educational only.
 */

import type { PrimeFieldCurveParams } from '../curves.js';

export interface EcdsaCurve {
  /** Field characteristic. */
  p: bigint;
  /** Curve coefficient a in y^2 = x^3 + a x + b. */
  a: bigint;
  /** Curve coefficient b. */
  b: bigint;
  /** Order of the generator's subgroup. */
  n: bigint;
  /** Generator x-coordinate. */
  Gx: bigint;
  /** Generator y-coordinate. */
  Gy: bigint;
}

export type AffinePoint = { x: bigint; y: bigint } | null;

export interface EcdsaSignature {
  r: bigint;
  s: bigint;
  /** Message hash reduced to an integer (truncated to bitlength(n) by the caller). */
  h: bigint;
}

/** Build an {@link EcdsaCurve} from a cheatsheet prime-field curve definition. */
export function ecdsaCurveFromParams(params: PrimeFieldCurveParams): EcdsaCurve {
  return {
    p: params.field.p,
    a: params.a,
    b: params.b,
    n: params.order,
    Gx: params.generator.x,
    Gy: params.generator.y,
  };
}

export function mod(a: bigint, m: bigint): bigint {
  const r = a % m;
  return r < 0n ? r + m : r;
}

/** Modular inverse via the extended Euclidean algorithm. Throws if not invertible. */
export function modInv(a: bigint, m: bigint): bigint {
  let [oldR, r] = [mod(a, m), m];
  let [oldS, s] = [1n, 0n];
  while (r !== 0n) {
    const q = oldR / r;
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
  }
  if (oldR !== 1n) {
    throw new Error(`not invertible: gcd(${a}, ${m}) = ${oldR}`);
  }
  return mod(oldS, m);
}

export function pointDouble(curve: EcdsaCurve, P: AffinePoint): AffinePoint {
  if (P === null) return null;
  if (P.y === 0n) return null;
  const { p, a } = curve;
  const lambda = mod((3n * P.x * P.x + a) * modInv(2n * P.y, p), p);
  const x = mod(lambda * lambda - 2n * P.x, p);
  const y = mod(lambda * (P.x - x) - P.y, p);
  return { x, y };
}

export function pointAdd(curve: EcdsaCurve, P: AffinePoint, Q: AffinePoint): AffinePoint {
  if (P === null) return Q;
  if (Q === null) return P;
  const { p } = curve;
  if (P.x === Q.x) {
    if (mod(P.y + Q.y, p) === 0n) return null;
    return pointDouble(curve, P);
  }
  const lambda = mod((Q.y - P.y) * modInv(Q.x - P.x, p), p);
  const x = mod(lambda * lambda - P.x - Q.x, p);
  const y = mod(lambda * (P.x - x) - P.y, p);
  return { x, y };
}

export function scalarMul(curve: EcdsaCurve, k: bigint, P: AffinePoint): AffinePoint {
  let scalar = mod(k, curve.n);
  let result: AffinePoint = null;
  let addend: AffinePoint = P;
  while (scalar > 0n) {
    if (scalar & 1n) result = pointAdd(curve, result, addend);
    addend = pointDouble(curve, addend);
    scalar >>= 1n;
  }
  return result;
}

export function pointEq(P: AffinePoint, Q: AffinePoint): boolean {
  if (P === null || Q === null) return P === Q;
  return P.x === Q.x && P.y === Q.y;
}

export function generator(curve: EcdsaCurve): AffinePoint {
  return { x: curve.Gx, y: curve.Gy };
}

/** Public key d*G. */
export function publicKey(curve: EcdsaCurve, d: bigint): AffinePoint {
  return scalarMul(curve, d, generator(curve));
}

/** True if d is the private key for the given public key. */
export function verifyPrivateKey(curve: EcdsaCurve, d: bigint, pub: AffinePoint): boolean {
  const dd = mod(d, curve.n);
  if (dd === 0n) return false;
  return pointEq(scalarMul(curve, dd, generator(curve)), pub);
}

/**
 * Produce an ECDSA signature with a caller-chosen nonce. Returns null when the
 * nonce yields r = 0 or s = 0 (the caller should pick another nonce), mirroring
 * how a real signer rejects such nonces.
 */
export function signWithNonce(
  curve: EcdsaCurve,
  d: bigint,
  k: bigint,
  h: bigint
): EcdsaSignature | null {
  const { n } = curve;
  const kk = mod(k, n);
  if (kk === 0n) return null;
  const R = scalarMul(curve, kk, generator(curve));
  if (R === null) return null;
  const r = mod(R.x, n);
  if (r === 0n) return null;
  const s = mod(modInv(kk, n) * (mod(h, n) + r * mod(d, n)), n);
  if (s === 0n) return null;
  return { r, s, h: mod(h, n) };
}
