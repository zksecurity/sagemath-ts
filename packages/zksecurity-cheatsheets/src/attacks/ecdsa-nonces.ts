/**
 * Recover ECDSA private keys from signatures produced with weak nonces.
 *
 * Covers algebraic attacks (nonce reuse, known linear relations, nonces that are
 * affine in the private key) and hidden-number-problem lattice attacks (known
 * zero most/least significant bits, an unknown shared prefix or postfix, and
 * nonces from a linear congruential generator such as java.util.Random).
 *
 * Inputs are the curve, the public key, and a list of signatures (r, s, h). All
 * recovered keys are verified against the public key before being returned, so a
 * returned key is always correct and a null means the attack did not succeed.
 *
 * Educational only: plain bigint arithmetic, not constant time, not audited.
 */

import {
  type IntegerMatrix,
  IntegerMatrixFromEntries,
  LLL,
} from '../../../sagemath-ts/src/matrix/matrix_integer.js';
import { BKZ } from '../../../sagemath-ts/src/modules/bkz.js';
import {
  type AffinePoint,
  type EcdsaCurve,
  type EcdsaSignature,
  mod,
  modInv,
  verifyPrivateKey,
} from './ecdsa-ec.js';

export type { EcdsaCurve, EcdsaSignature, AffinePoint };

/** One hidden-number-problem row: a small value `k = (u + t*x) mod n` with `0 <= k < bound`. */
export interface HnpInstance {
  t: bigint;
  u: bigint;
}

export interface LatticeOptions {
  /** Lattice reduction to use. BKZ is stronger but slower; default LLL. */
  reduction?: 'LLL' | 'BKZ';
  /** BKZ block size (ignored for LLL). */
  blockSize?: number;
}

/** Linear congruential generator `state' = (a*state + c) mod m`. */
export interface Lcg {
  a: bigint;
  c: bigint;
  m: bigint;
}

/** java.util.Random: 48-bit LCG. Its `nextLong`/seed stream advances by this recurrence. */
export const JAVA_LCG: Lcg = { a: 0x5deece66dn, c: 0xbn, m: 1n << 48n };

function bitLength(x: bigint): number {
  return x <= 0n ? 0 : x.toString(2).length;
}

function reduceLattice(rows: bigint[][], options?: LatticeOptions): bigint[][] {
  const matrix = IntegerMatrixFromEntries(rows);
  const reduced =
    options?.reduction === 'BKZ'
      ? BKZ(matrix, { blockSize: options.blockSize ?? Math.min(20, rows.length), delta: 0.99 })
      : (LLL(matrix, 0.99) as IntegerMatrix);
  const out: bigint[][] = [];
  for (let i = 0; i < reduced.nrows; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < reduced.ncols; j++) row.push(reduced.get(i, j).value);
    out.push(row);
  }
  return out;
}

/**
 * Solve a hidden-number-problem: given rows with `k_i = (u_i + t_i*x) mod n`,
 * `0 <= k_i < bound`, return candidate values of the hidden `x` (mod n). The
 * caller verifies each candidate against whatever ground truth it has.
 *
 * Uses the recentered Nguyen-Shparlinski embedding lattice, reduced with LLL (or
 * BKZ). The embedding scales the modular columns by `n` so that the target
 * short vector carries the hidden number in a dedicated coordinate.
 */
export function solveHnp(
  n: bigint,
  instances: HnpInstance[],
  bound: bigint,
  options?: LatticeOptions
): bigint[] {
  const m = instances.length;
  if (m === 0) return [];
  const half = bound >> 1n;
  const dim = m + 2;
  const K = n;

  const rows: bigint[][] = [];
  for (let i = 0; i < m; i++) {
    const row = new Array<bigint>(dim).fill(0n);
    row[i] = K * n;
    rows.push(row);
  }
  const rowT = new Array<bigint>(dim).fill(0n);
  const rowC = new Array<bigint>(dim).fill(0n);
  for (let i = 0; i < m; i++) {
    rowT[i] = K * mod(instances[i]!.t, n);
    rowC[i] = K * mod(instances[i]!.u - half, n);
  }
  rowT[m] = 1n;
  rowC[m + 1] = K * bound;
  rows.push(rowT, rowC);

  const reduced = reduceLattice(rows, options);

  const candidates = new Set<bigint>();
  const add = (x: bigint): void => {
    const v = mod(x, n);
    if (v !== 0n) candidates.add(v);
  };
  const tInv0 = modInv(mod(instances[0]!.t, n), n);
  for (const row of reduced) {
    const dCol = row[m]!;
    if (dCol !== 0n) {
      add(dCol);
      add(-dCol);
    }
    // Backup: reconstruct x from the first instance's recovered error term.
    const scaled = row[0]!;
    if (scaled % K === 0n) {
      const e0 = scaled / K;
      for (const k0 of [e0 + half, -e0 + half]) {
        if (k0 >= 0n && k0 < bound) {
          add((k0 - instances[0]!.u) * tInv0);
        }
      }
    }
  }
  return [...candidates];
}

function firstVerifying(curve: EcdsaCurve, pub: AffinePoint, candidates: bigint[]): bigint | null {
  for (const d of candidates) {
    if (verifyPrivateKey(curve, d, pub)) return d;
    const neg = mod(-d, curve.n);
    if (verifyPrivateKey(curve, neg, pub)) return neg;
  }
  return null;
}

// =============================================================================
// Algebraic attacks
// =============================================================================

/**
 * Recover the key when two signatures reused the same nonce (identical `r`).
 * Scans every pair sharing an `r` and returns the first recovered key.
 */
export function recoverFromReusedNonce(
  curve: EcdsaCurve,
  pub: AffinePoint,
  sigs: EcdsaSignature[]
): bigint | null {
  const { n } = curve;
  const byR = new Map<string, EcdsaSignature[]>();
  for (const sig of sigs) {
    const key = sig.r.toString();
    const bucket = byR.get(key);
    if (bucket) bucket.push(sig);
    else byR.set(key, [sig]);
  }
  for (const bucket of byR.values()) {
    for (let i = 0; i < bucket.length; i++) {
      for (let j = i + 1; j < bucket.length; j++) {
        const a = bucket[i]!;
        const b = bucket[j]!;
        const sDiff = mod(a.s - b.s, n);
        if (sDiff === 0n) continue;
        const k = mod((a.h - b.h) * modInv(sDiff, n), n);
        if (k === 0n) continue;
        const d = mod(mod(a.s * k - a.h, n) * modInv(a.r, n), n);
        if (verifyPrivateKey(curve, d, pub)) return d;
      }
    }
  }
  return null;
}

/**
 * Recover the key when the nonces are affine in one shared unknown scalar `t`:
 * `k_i = alpha_i + beta_i * t`, with `alpha_i`, `beta_i` known. Nonce reuse is the
 * case `alpha_i = 0, beta_i = 1`; a known offset between two nonces is `beta_i = 1`
 * with differing `alpha_i`. Any two signatures give a 2x2 linear system in
 * `(t, d)` solved mod n.
 */
export function recoverFromLinearNonceRelations(
  curve: EcdsaCurve,
  pub: AffinePoint,
  sigs: EcdsaSignature[],
  relations: { alpha: bigint; beta: bigint }[]
): bigint | null {
  const { n } = curve;
  if (sigs.length !== relations.length) {
    throw new Error('sigs and relations must have the same length');
  }
  // s_i*(alpha_i + beta_i*t) = h_i + r_i*d  =>  (s_i*beta_i)*t - r_i*d = h_i - s_i*alpha_i
  for (let i = 0; i < sigs.length; i++) {
    for (let j = i + 1; j < sigs.length; j++) {
      const a1 = mod(sigs[i]!.s * relations[i]!.beta, n);
      const b1 = mod(-sigs[i]!.r, n);
      const c1 = mod(sigs[i]!.h - sigs[i]!.s * relations[i]!.alpha, n);
      const a2 = mod(sigs[j]!.s * relations[j]!.beta, n);
      const b2 = mod(-sigs[j]!.r, n);
      const c2 = mod(sigs[j]!.h - sigs[j]!.s * relations[j]!.alpha, n);
      const det = mod(a1 * b2 - a2 * b1, n);
      if (det === 0n) continue;
      const detInv = modInv(det, n);
      const d = mod((a1 * c2 - a2 * c1) * detInv, n);
      if (verifyPrivateKey(curve, d, pub)) return d;
    }
  }
  return null;
}

/**
 * Recover the key when each nonce is affine in the private key itself:
 * `k_i = a_i + b_i * d`, with `a_i`, `b_i` known. A single signature suffices:
 * `s_i*(a_i + b_i*d) = h_i + r_i*d  =>  d = (h_i - s_i*a_i)/(s_i*b_i - r_i)`.
 */
export function recoverFromKeyDependentNonce(
  curve: EcdsaCurve,
  pub: AffinePoint,
  sigs: EcdsaSignature[],
  relations: { a: bigint; b: bigint }[]
): bigint | null {
  const { n } = curve;
  if (sigs.length !== relations.length) {
    throw new Error('sigs and relations must have the same length');
  }
  for (let i = 0; i < sigs.length; i++) {
    const denom = mod(sigs[i]!.s * relations[i]!.b - sigs[i]!.r, n);
    if (denom === 0n) continue;
    const d = mod((sigs[i]!.h - sigs[i]!.s * relations[i]!.a) * modInv(denom, n), n);
    if (verifyPrivateKey(curve, d, pub)) return d;
  }
  return null;
}

// =============================================================================
// Hidden-number-problem lattice attacks
// =============================================================================

/** Map a signature to `k = (u + t*d) mod n` with `t = r/s`, `u = h/s`. */
function hnpCoefficients(curve: EcdsaCurve, sig: EcdsaSignature): HnpInstance {
  const sInv = modInv(sig.s, curve.n);
  return { t: mod(sInv * sig.r, curve.n), u: mod(sInv * sig.h, curve.n) };
}

/**
 * Recover the key when the top `msbBits` bits of every nonce are zero (short
 * nonces). The nonce bound is `2^(bitlen(n) - msbBits)`.
 */
export function recoverWithKnownMsbNonce(
  curve: EcdsaCurve,
  pub: AffinePoint,
  sigs: EcdsaSignature[],
  msbBits: number,
  options?: LatticeOptions
): bigint | null {
  const { n } = curve;
  const bound = 1n << BigInt(bitLength(n) - msbBits);
  const instances = sigs.map((sig) => hnpCoefficients(curve, sig));
  return firstVerifying(curve, pub, solveHnp(n, instances, bound, options));
}

/**
 * Recover the key when the low `lsbBits` bits of every nonce are a known
 * constant `knownLow` (default zero). Dividing the relation by `2^lsbBits`
 * turns it into the short-nonce problem.
 */
export function recoverWithKnownLsbNonce(
  curve: EcdsaCurve,
  pub: AffinePoint,
  sigs: EcdsaSignature[],
  lsbBits: number,
  knownLow = 0n,
  options?: LatticeOptions
): bigint | null {
  const { n } = curve;
  const inv = modInv(1n << BigInt(lsbBits), n);
  const bound = 1n << BigInt(bitLength(n) - lsbBits);
  const instances = sigs.map((sig) => {
    const base = hnpCoefficients(curve, sig);
    return { t: mod(inv * base.t, n), u: mod(inv * (base.u - knownLow), n) };
  });
  return firstVerifying(curve, pub, solveHnp(n, instances, bound, options));
}

/**
 * Recover the key when all nonces share the same unknown high bits (a common
 * prefix), so only the low `varyingBits` differ. Differences of nonces are then
 * bounded by `2^varyingBits`, which is a hidden-number problem on the pairwise
 * differences against the first signature.
 */
export function recoverWithSharedNoncePrefix(
  curve: EcdsaCurve,
  pub: AffinePoint,
  sigs: EcdsaSignature[],
  varyingBits: number,
  options?: LatticeOptions
): bigint | null {
  const { n } = curve;
  const base = sigs.map((sig) => hnpCoefficients(curve, sig));
  const bound = 1n << BigInt(varyingBits + 1);
  const half = 1n << BigInt(varyingBits);
  const instances: HnpInstance[] = [];
  for (let i = 1; i < base.length; i++) {
    instances.push({
      t: mod(base[i]!.t - base[0]!.t, n),
      u: mod(base[i]!.u - base[0]!.u + half, n),
    });
  }
  return firstVerifying(curve, pub, solveHnp(n, instances, bound, options));
}

/**
 * Recover the key when all nonces share the same unknown low `sharedLowBits`
 * (a common postfix). The high parts `k_i >> sharedLowBits` are small and their
 * pairwise differences give a hidden-number problem after dividing by
 * `2^sharedLowBits`.
 */
export function recoverWithSharedNoncePostfix(
  curve: EcdsaCurve,
  pub: AffinePoint,
  sigs: EcdsaSignature[],
  sharedLowBits: number,
  options?: LatticeOptions
): bigint | null {
  const { n } = curve;
  const base = sigs.map((sig) => hnpCoefficients(curve, sig));
  const inv = modInv(1n << BigInt(sharedLowBits), n);
  const highBits = bitLength(n) - sharedLowBits;
  const bound = 1n << BigInt(highBits + 1);
  const half = 1n << BigInt(highBits);
  const instances: HnpInstance[] = [];
  for (let i = 1; i < base.length; i++) {
    instances.push({
      t: mod(inv * (base[i]!.t - base[0]!.t), n),
      u: mod(inv * (base[i]!.u - base[0]!.u) + half, n),
    });
  }
  return firstVerifying(curve, pub, solveHnp(n, instances, bound, options));
}

/**
 * Recover the key when the nonces are consecutive outputs of a linear
 * congruential generator (`state' = (a*state + c) mod lcg.m`), each output used
 * directly as a (short) nonce. The LCG recurrence relates consecutive nonces by
 * `k_{i+1} = a*k_i + c - m*w_i` with a small quotient `w_i < a`; eliminating the
 * private key gives a hidden-number problem in the quotients.
 *
 * `sigs` must be in generation order. Needs at least four signatures.
 */
export function recoverFromLcgNonces(
  curve: EcdsaCurve,
  pub: AffinePoint,
  sigs: EcdsaSignature[],
  lcg: Lcg,
  options?: LatticeOptions
): bigint | null {
  const { n } = curve;
  const { a, c, m } = lcg;
  const base = sigs.map((sig) => hnpCoefficients(curve, sig));
  // For each consecutive pair, d = const_i + coef_i * w_i (mod n), with w_i < a.
  const consts: bigint[] = [];
  const coefs: bigint[] = [];
  for (let i = 0; i + 1 < base.length; i++) {
    const dt = mod(base[i + 1]!.t - a * base[i]!.t, n);
    if (dt === 0n) return null;
    const invDt = modInv(dt, n);
    consts.push(mod(-invDt * (base[i + 1]!.u - a * base[i]!.u - c), n));
    coefs.push(mod(-invDt * m, n));
  }
  if (coefs.length < 2) return null;

  // Hidden number is w_0; every other w_i is a small residue tied to it.
  const coef0Inv0 = coefs[0]!;
  const instances: HnpInstance[] = [];
  for (let i = 1; i < coefs.length; i++) {
    const invCoef = modInv(coefs[i]!, n);
    instances.push({
      t: mod(invCoef * coef0Inv0, n),
      u: mod(-invCoef * (consts[i]! - consts[0]!), n),
    });
  }
  const candidatesW0 = solveHnp(n, instances, a, options);
  const keys = candidatesW0.map((w0) => mod(consts[0]! + coefs[0]! * w0, n));
  return firstVerifying(curve, pub, keys);
}
