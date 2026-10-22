import { describe, expect, test } from 'bun:test';
import { createHash, generateKeyPairSync } from 'node:crypto';
import {
  gcd,
  inverse_mod,
  is_prime,
  next_prime,
  power_mod,
} from '../../../sagemath-ts/src/arith/misc.js';
import {
  type RsaPublicKey,
  checkRsaKey,
  checkRsaKeys,
  fermatFactor,
  pollardPm1,
  rationalApproximation,
  rationalApproximationFactor,
  rocaDiscreteLog,
  smallFactors,
  wienerAttack,
} from './rsa-keys.js';
import { primesUpTo } from './rsa-math.js';

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

/** The first prime of the form f(i), i = 0, 1, 2, ... */
function firstPrime(f: (i: number) => bigint): bigint {
  for (let i = 0; ; i++) {
    const candidate = f(i);
    if (is_prime(candidate)) return candidate;
  }
}

function opensslKey(bits: number, e = 65537): RsaPublicKey {
  const { publicKey } = generateKeyPairSync('rsa', { modulusLength: bits, publicExponent: e });
  const { n } = publicKey.export({ format: 'jwk' });
  return { n: BigInt(`0x${Buffer.from(n!, 'base64url').toString('hex')}`), e: BigInt(e) };
}

const SOUND = opensslKey(2048);
const sorted = (a: bigint, b: bigint): [bigint, bigint] => (a < b ? [a, b] : [b, a]);

describe('fermatFactor', () => {
  test('factors primes about 2^518 apart', () => {
    const p = randomPrime(1024, 'fermat-p');
    const q = next_prime(p + random(518, 'fermat-gap'));
    expect(fermatFactor(p * q)).toEqual(sorted(p, q));
    const finding = checkRsaKey({ n: p * q, e: 65537n }).find((f) => f.check === 'fermat');
    expect(finding).toMatchObject({ severity: 'broken', p: sorted(p, q)[0], q: sorted(p, q)[1] });
  });

  test('factors primes sharing their top 256 and bottom 512 bits only with lowBits', () => {
    const high = random(256, 'shared-high') << 768n;
    const low = random(512, 'shared-low') | 1n;
    const p = firstPrime((i) => high + (random(256, `shared-p/${i}`) << 512n) + low);
    const q = firstPrime((i) => high + (random(256, `shared-q/${i}`) << 512n) + low);
    const n = p * q;
    expect(fermatFactor(n, { maxIterations: 1 << 16 })).toBeNull();
    expect(fermatFactor(n, { lowBits: 512, maxIterations: 16 })).toEqual(sorted(p, q));
    const finding = checkRsaKey({ n, e: 65537n }).find((f) => f.check === 'fermat-shared-bits');
    expect(finding).toMatchObject({ severity: 'broken', p: sorted(p, q)[0], q: sorted(p, q)[1] });
  });

  test('finds nothing on an OpenSSL key', () => {
    expect(fermatFactor(SOUND.n)).toBeNull();
    expect(fermatFactor(SOUND.n, { lowBits: 64 })).toBeNull();
  });
});

describe('smallFactors', () => {
  test('finds a small prime factor', () => {
    const n = 1009n * randomPrime(2038, 'small-q');
    expect(smallFactors(n)).toEqual([1009n]);
    const finding = checkRsaKey({ n, e: 65537n }).find((f) => f.check === 'small-factors');
    expect(finding).toMatchObject({ severity: 'broken', p: 1009n, q: n / 1009n });
  });

  test('finds nothing on an OpenSSL key', () => {
    expect(smallFactors(SOUND.n)).toEqual([]);
  });
});

describe('wienerAttack', () => {
  const p = randomPrime(1024, 'wiener-p');
  const q = randomPrime(1024, 'wiener-q');
  const phi = (p - 1n) * (q - 1n);
  const coprimeExponent = (bits: number, seed: string) => {
    let d = random(bits, seed) | 1n;
    while (gcd(d, phi) !== 1n) d += 2n;
    return d;
  };

  test('recovers a 500-bit private exponent', () => {
    const d = coprimeExponent(500, 'wiener-d');
    const e = inverse_mod(d, phi);
    expect(wienerAttack(p * q, e)).toEqual({ d, p: sorted(p, q)[0]!, q: sorted(p, q)[1]! });
    const finding = checkRsaKey({ n: p * q, e }).find((f) => f.check === 'wiener');
    expect(finding).toMatchObject({ severity: 'broken', d });
  });

  test('finds nothing for a full-size private exponent or e = 65537', () => {
    const e = inverse_mod(coprimeExponent(2000, 'wiener-big-d'), phi);
    expect(wienerAttack(p * q, e)).toBeNull();
    expect(wienerAttack(SOUND.n, SOUND.e)).toBeNull();
  });
});

describe('rationalApproximation', () => {
  const repeated = (byte: number) => BigInt(`0x${byte.toString(16).repeat(128)}`);

  test('factors primes with repeating byte patterns', () => {
    const p = next_prime(repeated(0xc1));
    const q = next_prime(repeated(0xc7));
    const n = p * q;
    expect(fermatFactor(n, { maxIterations: 1 << 16 })).toBeNull();
    const approx = rationalApproximation(n);
    expect(approx?.precisionBits).toBeGreaterThan(900);
    expect(rationalApproximationFactor(n, approx!)).toEqual(sorted(p, q));
    const finding = checkRsaKey({ n, e: 65537n }).find((f) => f.check === 'rational-approximation');
    expect(finding).toMatchObject({ severity: 'broken', p: sorted(p, q)[0], q: sorted(p, q)[1] });
  });

  test('factors primes near 3·2^1022 and 2^1023 with 300 random low bits', () => {
    const p = next_prime((3n << 1022n) + random(300, 'rational-p'));
    const q = next_prime((1n << 1023n) + random(300, 'rational-q'));
    const approx = rationalApproximation(p * q);
    expect(approx).toMatchObject({ numerator: 3n, denominator: 4n, exponent: 2047 });
    expect(rationalApproximationFactor(p * q, approx!)).toEqual(sorted(p, q));
  });

  test('finds nothing on an OpenSSL key', () => {
    expect(rationalApproximation(SOUND.n)).toBeNull();
  });
});

describe('pollardPm1', () => {
  const primes = primesUpTo(100_000);
  /** A prime p of at least 1024 bits with p − 1 = 2·cofactor·(primes below 10^5). */
  const smoothPrime = (seed: string, cofactor: bigint) =>
    firstPrime((attempt) => {
      let m = 2n * cofactor;
      for (let i = 0; m < 1n << 1023n; i++) {
        m *= primes[Number(random(32, `${seed}/${attempt}/${i}`) % BigInt(primes.length))]!;
      }
      return m + 1n;
    });
  const q = randomPrime(1024, 'pm1-q');

  test('stage 1 splits n when p − 1 is 10^5-smooth', () => {
    const p = smoothPrime('pm1-stage1', 1n);
    expect(pollardPm1(p * q, { bound: 100_000, bound2: 100_000 })).toEqual(sorted(p, q));
    const finding = checkRsaKey({ n: p * q, e: 65537n }).find((f) => f.check === 'pollard-p-1');
    expect(finding).toMatchObject({ severity: 'broken', p: sorted(p, q)[0], q: sorted(p, q)[1] });
  });

  test('stage 2 handles one prime of p − 1 between the bounds', () => {
    const p = smoothPrime('pm1-stage2', next_prime(700_000n));
    expect(pollardPm1(p * q, { bound: 100_000, bound2: 100_000 })).toBeNull();
    expect(pollardPm1(p * q)).toEqual(sorted(p, q));
  });

  test('finds nothing on an OpenSSL key', () => {
    expect(pollardPm1(SOUND.n)).toBeNull();
  });
});

describe('rocaDiscreteLog', () => {
  test('fingerprints a 2048-bit key built like Infineon RSALib', () => {
    // RSALib's 2048-bit keys use M = the product of the first 126 primes.
    const M = primesUpTo(701).reduce((a, b) => a * b, 1n);
    const kBits = 1024 - M.toString(2).length;
    const rocaPrime = (seed: string) =>
      firstPrime(
        (i) =>
          random(kBits, `${seed}/k/${i}`) * M + power_mod(65537n, random(62, `${seed}/a/${i}`), M)
      );
    const n = rocaPrime('roca-p') * rocaPrime('roca-q');
    const c = rocaDiscreteLog(n);
    const M39 = primesUpTo(167).reduce((a, b) => a * b, 1n);
    expect(c).not.toBeNull();
    expect(power_mod(65537n, c!, M39)).toBe(n % M39);
    expect(checkRsaKey({ n, e: 65537n }).find((f) => f.check === 'roca')?.severity).toBe(
      'vulnerable'
    );
  });

  test('finds nothing on an OpenSSL key', () => {
    expect(rocaDiscreteLog(SOUND.n)).toBeNull();
  });
});

describe('checkRsaKey', () => {
  test('reports nothing for a 2048-bit OpenSSL key', () => {
    expect(checkRsaKey(SOUND)).toEqual([]);
  });

  test('warns about a 1024-bit key with e = 3', () => {
    const findings = checkRsaKey(opensslKey(1024, 3));
    expect(findings.map((f) => [f.check, f.severity])).toEqual([
      ['size', 'warning'],
      ['exponent', 'warning'],
    ]);
  });

  test('runs only the checks it is given', () => {
    expect(checkRsaKey({ n: 1009n * 1013n, e: 65537n }, { checks: [] })).toEqual([]);
  });

  test('checkRsaKeys combines per-key and batch findings', () => {
    const p = randomPrime(1024, 'batch-p');
    const keys = [
      { n: p * randomPrime(1024, 'batch-q0'), e: 65537n },
      SOUND,
      { n: p * randomPrime(1024, 'batch-q1'), e: 3n },
    ];
    const report = checkRsaKeys(keys);
    expect(report.keys.map((fs) => fs.map((f) => f.check))).toEqual([[], [], ['exponent']]);
    expect(report.batch.map((f) => [f.check, f.keys])).toEqual([['shared-prime', [0, 2]]]);
  });
});
