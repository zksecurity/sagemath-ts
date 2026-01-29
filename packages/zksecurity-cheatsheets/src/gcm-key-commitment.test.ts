import { describe, test, expect } from 'bun:test';
import { randomBytes } from 'node:crypto';
import {
  GF2_128,
  GCM_POLYNOMIAL,
  ghash,
  computeH,
  computeJ,
  aesGcmEncrypt,
  aesGcmDecrypt,
  demonstrateKeyCommitmentAttack,
  multiKeyCollision,
} from './gcm-key-commitment.js';

describe('GF(2^128)', () => {
  test('zero and one', () => {
    const zero = GF2_128.zero();
    const one = GF2_128.one();
    expect(zero.isZero()).toBe(true);
    expect(one.isZero()).toBe(false);
    expect(zero.add(one).eq(one)).toBe(true);
  });

  test('addition is XOR', () => {
    const a = GF2_128.fromBigInt(0x123456789abcdef0n);
    const b = GF2_128.fromBigInt(0xfedcba9876543210n);
    const sum = a.add(b);
    expect(sum.toBigInt()).toBe(0x123456789abcdef0n ^ 0xfedcba9876543210n);
  });

  test('self-addition is zero', () => {
    const a = GF2_128.fromBigInt(0xdeadbeefcafebabebeefcafedeadbeefn);
    expect(a.add(a).isZero()).toBe(true);
  });

  test('multiplication by one', () => {
    const a = GF2_128.fromBigInt(0x123456789abcdef0fedcba9876543210n);
    const one = GF2_128.one();
    expect(a.mul(one).eq(a)).toBe(true);
  });

  test('multiplication by zero', () => {
    const a = GF2_128.fromBigInt(0x123456789abcdef0fedcba9876543210n);
    const zero = GF2_128.zero();
    expect(a.mul(zero).isZero()).toBe(true);
  });

  test('inverse', () => {
    const a = GF2_128.fromBigInt(0xdeadbeefcafebabebeefcafedeadbeefn);
    const aInv = a.inv();
    const product = a.mul(aInv);
    expect(product.eq(GF2_128.one())).toBe(true);
  });

  test('bytes roundtrip', () => {
    const original = randomBytes(16);
    const elem = GF2_128.fromBytes(original);
    const recovered = elem.toBytes();
    expect(recovered.equals(original)).toBe(true);
  });

  test('power', () => {
    const a = GF2_128.fromBigInt(0xabcdef0123456789n);
    const a2 = a.pow(2n);
    const aSquared = a.sqr();
    expect(a2.eq(aSquared)).toBe(true);

    const a4 = a.pow(4n);
    expect(a4.eq(a2.sqr())).toBe(true);
  });
});

describe('GHASH', () => {
  test('empty inputs', () => {
    const H = GF2_128.fromBigInt(0x12345678n);
    const result = ghash(H, Buffer.alloc(0), Buffer.alloc(0));
    // Empty input should give just the length block contribution
    // L = 0||0, so GHASH = L*H = 0
    expect(result.isZero()).toBe(true);
  });

  test('single block ciphertext', () => {
    const H = GF2_128.fromBigInt(0xabcdef0123456789fedcba9876543210n);
    const ct = randomBytes(16);
    const result = ghash(H, Buffer.alloc(0), ct);

    // Manual computation: C*H^2 + L*H where L encodes len(ct)=128 bits
    const C = GF2_128.fromBytes(ct);
    const lenBlock = Buffer.alloc(16);
    lenBlock.writeBigUInt64BE(0n, 0);
    lenBlock.writeBigUInt64BE(128n, 8);
    const L = GF2_128.fromBytes(lenBlock);

    const expected = C.mul(H.sqr()).add(L.mul(H));
    expect(result.eq(expected)).toBe(true);
  });
});

describe('AES-GCM helpers', () => {
  test('encrypt and decrypt roundtrip', () => {
    const key = randomBytes(32);
    const iv = randomBytes(12);
    const plaintext = Buffer.from('Hello, World!!! (16 bytes)');
    const aad = Buffer.from('additional data');

    const { ciphertext, tag } = aesGcmEncrypt(key, iv, plaintext, aad);
    const decrypted = aesGcmDecrypt(key, iv, ciphertext, tag, aad);

    expect(decrypted.equals(plaintext)).toBe(true);
  });

  test('computeH is deterministic', () => {
    const key = randomBytes(32);
    const H1 = computeH(key);
    const H2 = computeH(key);
    expect(H1.eq(H2)).toBe(true);
  });

  test('computeJ is deterministic', () => {
    const key = randomBytes(32);
    const iv = randomBytes(12);
    const J1 = computeJ(key, iv);
    const J2 = computeJ(key, iv);
    expect(J1.eq(J2)).toBe(true);
  });
});

describe('Key Commitment Attack', () => {
  test('basic attack produces different plaintexts', () => {
    const result = demonstrateKeyCommitmentAttack();

    // Keys should be different
    expect(result.key1.equals(result.key2)).toBe(false);

    // Plaintexts should be different (with overwhelming probability)
    expect(result.plaintext1.equals(result.plaintext2)).toBe(false);

    // Both should decrypt with their respective keys
    const decrypted1 = aesGcmDecrypt(
      result.key1,
      result.iv,
      result.ciphertext,
      result.tag,
      result.aad
    );
    const decrypted2 = aesGcmDecrypt(
      result.key2,
      result.iv,
      result.ciphertext,
      result.tag,
      result.aad
    );

    expect(decrypted1.equals(result.plaintext1)).toBe(true);
    expect(decrypted2.equals(result.plaintext2)).toBe(true);
  });

  test('attack with fixed keys', () => {
    // Use deterministic keys for reproducibility
    const key1 = Buffer.alloc(32);
    const key2 = Buffer.alloc(32);
    key1.fill(0x01);
    key2.fill(0x02);

    const result = demonstrateKeyCommitmentAttack({ key1, key2 });

    expect(result.key1.equals(key1)).toBe(true);
    expect(result.key2.equals(key2)).toBe(true);

    // Verify tag is valid for both keys
    const decrypted1 = aesGcmDecrypt(
      key1,
      result.iv,
      result.ciphertext,
      result.tag,
      result.aad
    );
    const decrypted2 = aesGcmDecrypt(
      key2,
      result.iv,
      result.ciphertext,
      result.tag,
      result.aad
    );

    // Plaintexts should differ
    expect(decrypted1.equals(decrypted2)).toBe(false);
  });

  test('attack works with custom IV', () => {
    const iv = Buffer.from('custom_iv_12'); // 12 bytes
    const result = demonstrateKeyCommitmentAttack({ iv });

    expect(result.iv.equals(iv)).toBe(true);
  });
});

describe('Multi-Key Collision', () => {
  test('two keys collision', () => {
    const keys = [randomBytes(32), randomBytes(32)];
    const iv = randomBytes(12);

    const { ciphertext, tag, plaintexts } = multiKeyCollision(keys, iv);

    // Verify both keys produce valid decryption
    for (let i = 0; i < keys.length; i++) {
      const decrypted = aesGcmDecrypt(keys[i]!, iv, ciphertext, tag);
      expect(decrypted.equals(plaintexts[i]!)).toBe(true);
    }

    // Plaintexts should be different
    expect(plaintexts[0]!.equals(plaintexts[1]!)).toBe(false);
  });

  test('three keys collision', () => {
    const keys = [randomBytes(32), randomBytes(32), randomBytes(32)];
    const iv = randomBytes(12);

    const { ciphertext, tag, plaintexts } = multiKeyCollision(keys, iv);

    // Verify all keys produce valid decryption
    for (let i = 0; i < keys.length; i++) {
      const decrypted = aesGcmDecrypt(keys[i]!, iv, ciphertext, tag);
      expect(decrypted.equals(plaintexts[i]!)).toBe(true);
    }

    // All plaintexts should be distinct (with overwhelming probability)
    expect(plaintexts[0]!.equals(plaintexts[1]!)).toBe(false);
    expect(plaintexts[0]!.equals(plaintexts[2]!)).toBe(false);
    expect(plaintexts[1]!.equals(plaintexts[2]!)).toBe(false);
  });
});

describe('Attack verification', () => {
  test('GHASH values match expected tag relationship', () => {
    const result = demonstrateKeyCommitmentAttack();

    // Manually compute GHASH and verify tag equation
    const ghash1 = ghash(result.H1, result.aad, result.ciphertext);
    const ghash2 = ghash(result.H2, result.aad, result.ciphertext);

    const J1 = computeJ(result.key1, result.iv);
    const J2 = computeJ(result.key2, result.iv);

    const tag1 = ghash1.add(J1);
    const tag2 = ghash2.add(J2);

    // Tags should be equal
    expect(tag1.eq(tag2)).toBe(true);

    // Tag should match the computed tag
    const expectedTag = GF2_128.fromBytes(result.tag);
    expect(tag1.eq(expectedTag)).toBe(true);
  });
});
