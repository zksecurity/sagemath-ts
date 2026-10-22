import { describe, expect, test } from 'bun:test';
import { createHash, generateKeyPairSync } from 'node:crypto';
import { is_prime, next_prime } from '../../../sagemath-ts/src/arith/misc.js';
import {
  batchGcd,
  factorWithNMinus1Divisor,
  findSharedNMinus1Factors,
  findSharedPrimes,
} from './rsa-batch.js';

function random(bits: number, seed: string): bigint {
  let hex = '';
  for (let i = 0; hex.length * 4 < bits; i++) {
    hex += createHash('sha256').update(`${seed}/${i}`).digest('hex');
  }
  return (BigInt(`0x${hex}`) >> BigInt(hex.length * 4 - bits)) | (1n << BigInt(bits - 1));
}

function randomPrime(bits: number, seed: string): bigint {
  return next_prime(random(bits, seed));
}

function opensslModulus(): bigint {
  const { publicKey } = generateKeyPairSync('rsa', { modulusLength: 2048 });
  const { n } = publicKey.export({ format: 'jwk' });
  return BigInt(`0x${Buffer.from(n!, 'base64url').toString('hex')}`);
}

const SOUND = [opensslModulus(), opensslModulus(), opensslModulus()];

describe('batchGcd', () => {
  test('returns gcd of each value with the product of the others', () => {
    expect(batchGcd([6n, 10n, 35n, 11n])).toEqual([2n, 10n, 5n, 1n]);
  });
});

describe('findSharedPrimes', () => {
  const p = randomPrime(1024, 'shared-p');
  const q0 = randomPrime(1024, 'shared-q0');
  const q1 = randomPrime(1024, 'shared-q1');

  test('factors moduli sharing primes, including one whose primes are both shared', () => {
    const moduli = [p * q0, p * q1, SOUND[0]!, p * q0, q0 * q1];
    const findings = findSharedPrimes(moduli);
    expect(findings.map((f) => [f.check, f.keys, f.shared])).toEqual([
      ['duplicate-modulus', [0, 3], undefined],
      ['shared-prime', [0, 1, 3], p],
      ['shared-prime', [0, 3, 4], q0],
      ['shared-prime', [1, 4], q1],
    ]);
    for (const { key, p: a, q: b } of findings.flatMap((f) => f.factored ?? [])) {
      expect(a * b).toBe(moduli[key]!);
      expect(a > 1n && a <= b).toBe(true);
    }
  });

  test('finds nothing among OpenSSL keys', () => {
    expect(findSharedPrimes(SOUND)).toEqual([]);
  });
});

describe('findSharedNMinus1Factors', () => {
  const r = randomPrime(600, 'pm1-r');
  /** A prime 1 + r·x with an even 424-bit x. */
  const primeOneModR = (seed: string) => {
    for (let i = 0; ; i++) {
      const candidate = 1n + r * (random(424, `${seed}/${i}`) << 1n);
      if (is_prime(candidate)) return candidate;
    }
  };
  const p1 = primeOneModR('p1');
  const q1 = primeOneModR('q1');
  const p2 = primeOneModR('p2');
  const q2 = primeOneModR('q2');

  test('finds r and factors both keys', () => {
    const moduli = [p1 * q1, SOUND[0]!, p2 * q2];
    const findings = findSharedNMinus1Factors(moduli);
    expect(findings).toHaveLength(1);
    expect(findings[0]).toMatchObject({
      check: 'shared-n-minus-1',
      severity: 'broken',
      keys: [0, 2],
    });
    expect(findings[0]!.shared! % r).toBe(0n);
    expect(findings[0]!.factored).toEqual([
      { key: 0, p: p1 < q1 ? p1 : q1, q: p1 < q1 ? q1 : p1 },
      { key: 2, p: p2 < q2 ? p2 : q2, q: p2 < q2 ? q2 : p2 },
    ]);
  });

  test('factorWithNMinus1Divisor needs the right r', () => {
    expect(factorWithNMinus1Divisor(p1 * q1, r)).not.toBeNull();
    expect(factorWithNMinus1Divisor(SOUND[0]!, r)).toBeNull();
  });

  test('finds nothing among OpenSSL keys', () => {
    expect(findSharedNMinus1Factors(SOUND)).toEqual([]);
  });
});
