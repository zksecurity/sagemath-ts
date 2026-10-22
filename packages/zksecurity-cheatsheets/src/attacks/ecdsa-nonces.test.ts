import { describe, expect, test } from 'bun:test';
import { p256, secp256k1 } from '../curves.js';
import {
  type EcdsaCurve,
  type EcdsaSignature,
  ecdsaCurveFromParams,
  mod,
  publicKey,
  signWithNonce,
} from './ecdsa-ec.js';
import {
  JAVA_LCG,
  recoverFromKeyDependentNonce,
  recoverFromLcgNonces,
  recoverFromLinearNonceRelations,
  recoverFromReusedNonce,
  recoverWithKnownLsbNonce,
  recoverWithKnownMsbNonce,
  recoverWithSharedNoncePostfix,
  recoverWithSharedNoncePrefix,
} from './ecdsa-nonces.js';

const SECP256K1: EcdsaCurve = ecdsaCurveFromParams(secp256k1);
const P256: EcdsaCurve = ecdsaCurveFromParams(p256);

const MASK64 = (1n << 64n) - 1n;

function makeRng(seed: bigint): (bits: number) => bigint {
  let state = seed & MASK64;
  const next64 = (): bigint => {
    state = (state + 0x9e3779b97f4a7c15n) & MASK64;
    let z = state;
    z = ((z ^ (z >> 30n)) * 0xbf58476d1ce4e5b9n) & MASK64;
    z = ((z ^ (z >> 27n)) * 0x94d049bb133111ebn) & MASK64;
    z = (z ^ (z >> 31n)) & MASK64;
    return z;
  };
  return (bits: number): bigint => {
    let r = 0n;
    let got = 0;
    while (got < bits) {
      r = (r << 64n) | next64();
      got += 64;
    }
    return r & ((1n << BigInt(bits)) - 1n);
  };
}

/** Sign with a chosen nonce, retrying nothing (the caller controls k). */
function sign(curve: EcdsaCurve, d: bigint, k: bigint, h: bigint): EcdsaSignature {
  const sig = signWithNonce(curve, d, k, h);
  if (sig === null) throw new Error('degenerate nonce in test');
  return sig;
}

describe('ECDSA algebraic nonce attacks', () => {
  test('nonce reuse recovers the key; distinct nonces do not', () => {
    const rng = makeRng(1n);
    for (const curve of [SECP256K1, P256]) {
      const d = (rng(256) % (curve.n - 1n)) + 1n;
      const pub = publicKey(curve, d);
      const k = (rng(256) % (curve.n - 1n)) + 1n;
      const s1 = sign(curve, d, k, rng(256) % curve.n);
      const s2 = sign(curve, d, k, rng(256) % curve.n);
      expect(s1.r).toBe(s2.r);
      expect(recoverFromReusedNonce(curve, pub, [s1, s2])).toBe(d);

      const distinct = [
        sign(curve, d, (rng(256) % (curve.n - 1n)) + 1n, rng(256) % curve.n),
        sign(curve, d, (rng(256) % (curve.n - 1n)) + 1n, rng(256) % curve.n),
      ];
      expect(recoverFromReusedNonce(curve, pub, distinct)).toBeNull();
    }
  });

  test('known linear relation between nonces (shared t with offset)', () => {
    const rng = makeRng(2n);
    const curve = SECP256K1;
    const d = (rng(256) % (curve.n - 1n)) + 1n;
    const pub = publicKey(curve, d);
    const t = (rng(256) % (curve.n - 1n)) + 1n;
    const offsets = [0n, 12345n, mod(-99n, curve.n)];
    const relations = offsets.map((o) => ({ alpha: o, beta: 1n }));
    const sigs = offsets.map((o) => sign(curve, d, mod(t + o, curve.n), rng(256) % curve.n));
    expect(recoverFromLinearNonceRelations(curve, pub, sigs, relations)).toBe(d);
  });

  test('nonce affine in the private key, single signature', () => {
    const rng = makeRng(3n);
    const curve = P256;
    const d = (rng(256) % (curve.n - 1n)) + 1n;
    const pub = publicKey(curve, d);
    const a = 7n;
    const b = 3n;
    const k = mod(a + b * d, curve.n);
    const sig = sign(curve, d, k, rng(256) % curve.n);
    expect(recoverFromKeyDependentNonce(curve, pub, [sig], [{ a, b }])).toBe(d);
  });
});

describe('ECDSA HNP lattice attacks', () => {
  function biasedMsbKey(curve: EcdsaCurve, seed: bigint, nonceBits: number, count: number) {
    const rng = makeRng(seed);
    const d = (rng(256) % (curve.n - 1n)) + 1n;
    const pub = publicKey(curve, d);
    const sigs: EcdsaSignature[] = [];
    while (sigs.length < count) {
      const k = rng(nonceBits);
      if (k === 0n) continue;
      sigs.push(sign(curve, d, k, rng(256) % curve.n));
    }
    return { d, pub, sigs };
  }

  test('known-zero MSBs (short nonces) on secp256k1 and P-256', () => {
    const msbBits = 120;
    const nonceBits = 256 - msbBits;
    for (const curve of [SECP256K1, P256]) {
      const { d, pub, sigs } = biasedMsbKey(curve, 10n, nonceBits, 8);
      expect(recoverWithKnownMsbNonce(curve, pub, sigs, msbBits)).toBe(d);
    }
  });

  test('no false positive: full-entropy nonces are not recoverable as short', () => {
    const curve = SECP256K1;
    const rng = makeRng(11n);
    const d = (rng(256) % (curve.n - 1n)) + 1n;
    const pub = publicKey(curve, d);
    const sigs: EcdsaSignature[] = [];
    while (sigs.length < 8) {
      const k = (rng(256) % (curve.n - 1n)) + 1n;
      sigs.push(sign(curve, d, k, rng(256) % curve.n));
    }
    expect(recoverWithKnownMsbNonce(curve, pub, sigs, 120)).toBeNull();
  });

  test('known-zero LSBs', () => {
    const curve = SECP256K1;
    const lsbBits = 120;
    const rng = makeRng(20n);
    const d = (rng(256) % (curve.n - 1n)) + 1n;
    const pub = publicKey(curve, d);
    const sigs: EcdsaSignature[] = [];
    while (sigs.length < 8) {
      const high = rng(256 - lsbBits);
      const k = high << BigInt(lsbBits);
      if (k === 0n) continue;
      sigs.push(sign(curve, d, mod(k, curve.n), rng(256) % curve.n));
    }
    expect(recoverWithKnownLsbNonce(curve, pub, sigs, lsbBits)).toBe(d);
  });

  test('shared unknown nonce prefix (common high bits)', () => {
    const curve = SECP256K1;
    const varyingBits = 120;
    const rng = makeRng(30n);
    const d = (rng(256) % (curve.n - 1n)) + 1n;
    const pub = publicKey(curve, d);
    const prefix = rng(256 - varyingBits) << BigInt(varyingBits);
    const sigs: EcdsaSignature[] = [];
    while (sigs.length < 9) {
      const k = mod(prefix + rng(varyingBits), curve.n);
      if (k === 0n) continue;
      sigs.push(sign(curve, d, k, rng(256) % curve.n));
    }
    expect(recoverWithSharedNoncePrefix(curve, pub, sigs, varyingBits)).toBe(d);
  });

  test('shared unknown nonce postfix (common low bits)', () => {
    const curve = SECP256K1;
    const sharedLowBits = 150;
    const highBits = 256 - sharedLowBits;
    const rng = makeRng(40n);
    const d = (rng(256) % (curve.n - 1n)) + 1n;
    const pub = publicKey(curve, d);
    const postfix = rng(sharedLowBits);
    const sigs: EcdsaSignature[] = [];
    while (sigs.length < 9) {
      const high = rng(highBits);
      const k = mod((high << BigInt(sharedLowBits)) + postfix, curve.n);
      if (k === 0n) continue;
      sigs.push(sign(curve, d, k, rng(256) % curve.n));
    }
    expect(recoverWithSharedNoncePostfix(curve, pub, sigs, sharedLowBits)).toBe(d);
  });
});

describe('ECDSA LCG nonce attack', () => {
  test('java.util.Random (48-bit LCG) consecutive nonces', () => {
    const curve = SECP256K1;
    const rng = makeRng(50n);
    const d = (rng(256) % (curve.n - 1n)) + 1n;
    const pub = publicKey(curve, d);
    let state = rng(48);
    const sigs: EcdsaSignature[] = [];
    for (let i = 0; i < 6; i++) {
      const k = state === 0n ? 1n : state;
      sigs.push(sign(curve, d, k, rng(256) % curve.n));
      state = mod(JAVA_LCG.a * state + JAVA_LCG.c, JAVA_LCG.m);
    }
    expect(recoverFromLcgNonces(curve, pub, sigs, JAVA_LCG)).toBe(d);
  });
});
