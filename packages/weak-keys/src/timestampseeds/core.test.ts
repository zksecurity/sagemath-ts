import { describe, expect, test } from 'bun:test';
import { randomInt } from 'node:crypto';
import { readFileSync } from 'node:fs';
import fc from 'fast-check';
import {
  CURVES,
  type Curve,
  PROFILES,
  type Profile,
  candidate,
  draws_per_key,
  parse_time,
  scalar,
  scalar_from_words,
  stream_at,
} from './core.js';
import { json, parse_json } from './json.js';
import { RNGS, type Rng, type Word, words } from './rngs.js';
import { fingerprint } from './targets.js';

const take = (rng: Rng, seed: bigint, count: number) => {
  const stream = words(rng, seed);
  return Array.from({ length: count }, () => stream.next().value![0]);
};
function block(value: bigint): Iterator<Word> {
  const buffer = Buffer.from(value.toString(16).padStart(64, '0'), 'hex');
  return Array.from({ length: 8 }, (_, i): Word => [BigInt(buffer.readUInt32BE(i * 4)), 4])[
    Symbol.iterator
  ]();
}

describe('timestamp RNG/profile regressions', () => {
  test('native known-answer words and distinct MT seed expansion', () => {
    for (const [rng, seed, expected] of [
      ['glibc-rand', 1n, [1804289383n, 846930886n, 1681692777n, 1714636915n, 1957747793n]],
      ['ansi-lcg31', 1n, [1103527590n, 377401575n, 662824084n]],
      ['python-mt19937', 1n, [577090037n]],
      ['java-random', 0n, [3139482720n, 3571011896n, 1033096058n, 2604232894n, 2737687030n]],
      ['msvc-rand', 1n, [41n, 18467n, 6334n, 26500n, 19169n]],
      ['cpp-mt19937', 5489n, [3499211612n, 581869302n, 3890346734n, 3586334585n, 545404204n]],
    ] as const)
      expect(take(rng, seed, expected.length)).toEqual([...expected]);
  });
  test('invalid seeds and algorithms fail', () => {
    expect(() => words('lrand48', -1n).next()).toThrow('nonnegative');
    expect(() => words('lrand48', 1 as unknown as bigint).next()).toThrow('nonnegative');
    expect(() => words('invalid' as Rng, 1n).next()).toThrow('Unknown RNG');
    for (const profile of [
      'python-getrandbits-reject',
      'python-randbytes-reject',
      'java-nextbytes-reject',
    ] as Profile[])
      expect(() => draws_per_key('lrand48', profile)).toThrow('requires --rng');
    expect(() => stream_at('lrand48', 1n, 'concat-be-reject', -1)).toThrow('Skipped draws');
    expect(() => stream_at('lrand48', 1n, 'concat-be-reject', 0, 1000)).toThrow('Key index');
  });
  test('word endianness differs from whole-buffer endianness', () => {
    const bytes = Buffer.from(Array.from({ length: 32 }, (_, i) => i + 1));
    const value = BigInt(`0x${bytes.toString('hex')}`);
    const reversed = BigInt(`0x${Buffer.from(bytes).reverse().toString('hex')}`);
    expect(scalar_from_words(block(value), 'buffer-le-reject', 1n << 256n)).toBe(reversed);
    expect(scalar_from_words(block(value), 'concat-le-reject', 1n << 256n)).not.toBe(reversed);
  });
  test('low-byte consumes 32 words and six-byte words discard the unused tail', () => {
    const stream = Array.from({ length: 33 }, (_, i): Word => [BigInt(0x123400 + i), 4])[
      Symbol.iterator
    ]();
    expect(scalar_from_words(stream, 'low-byte-reject', 1n << 256n)).toBe(
      BigInt(`0x${Buffer.from(Array.from({ length: 32 }, (_, i) => i)).toString('hex')}`)
    );
    expect(stream.next().value).toEqual([0x123420n, 4]);
    const words6 = Array.from(
      { length: 7 },
      (_, i): Word => [
        BigInt(
          `0x${Buffer.from(Array.from({ length: 6 }, (_, j) => i * 6 + j + 1)).toString('hex')}`
        ),
        6,
      ]
    );
    for (const profile of ['concat-be-reject', 'concat-le-reject'] as Profile[]) {
      const iterator = words6[Symbol.iterator]();
      const bytes = Buffer.concat(
        words6.slice(0, 6).map(([word]) => {
          const b = Buffer.from(word.toString(16).padStart(12, '0'), 'hex');
          return profile === 'concat-le-reject' ? b.reverse() : b;
        })
      ).subarray(0, 32);
      expect(scalar_from_words(iterator, profile, 1n << 256n)).toBe(
        BigInt(`0x${bytes.toString('hex')}`)
      );
      expect(iterator.next().value).toEqual(words6[6]);
    }
  });
  test('rejection consumes exactly one block, including zero and subgroup order', () => {
    for (const { order } of Object.values(CURVES))
      for (const value of [0n, 1n, order - 1n, order, order + 1n, (1n << 256n) - 1n]) {
        expect(scalar_from_words(block(value), 'concat-be-reject', order)).toBe(
          value > 0n && value < order ? value : null
        );
        expect(scalar_from_words(block(value), 'concat-be-mod', order)).toBe(
          (value % (order - 1n)) + 1n
        );
      }
    const stream = [
      ...Array.from({ length: 8 }, (): Word => [0n, 4]),
      ...Array.from({ length: 7 }, (): Word => [0n, 4]),
      [1n, 4] as Word,
      ...Array.from({ length: 7 }, (): Word => [0n, 4]),
      [2n, 4] as Word,
    ][Symbol.iterator]();
    expect(scalar_from_words(stream, 'concat-be-reject', 2n)).toBeNull();
    expect(scalar_from_words(stream, 'concat-be-reject', 2n)).toBe(1n);
    expect(scalar_from_words(stream, 'concat-be-reject', 2n)).toBeNull();
    expect(stream.next().done).toBe(true);
  });
  test('existing upstream public fixtures retain points and hashes, including >64-bit seeds', () => {
    const fixtures = parse_json(
      readFileSync(
        new URL('../../timestampseeds/tests/fixtures/public_keys.json', import.meta.url),
        'utf8'
      )
    ) as {
      cases: {
        rng: Rng;
        seed: number | bigint;
        profile: Profile;
        curve: Curve;
        skip_draws: number;
        key_index: number;
        public: null | { sec1_compressed_hex: string; spki_sha256: string };
      }[];
    };
    for (const row of fixtures.cases) {
      const result = candidate(
        BigInt(row.seed),
        'seconds',
        row.rng,
        row.profile,
        row.curve,
        row.skip_draws,
        row.key_index
      );
      if (row.public === null) expect(result).toBeNull();
      else {
        expect(result).toMatchObject(row.public);
        expect(fingerprint(Buffer.from(result!.spki_der_base64, 'base64'))).toBe(
          result!.spki_sha256
        );
      }
    }
  });
});

test('exact timestamp boundaries and timezone formats', () => {
  for (const [stamp, unit, expected] of [
    ['1970-01-01T01:00:01+01:00', 'seconds', 1n],
    ['1970-01-01T00:00:00.001Z', 'milliseconds', 1n],
    ['1970-01-01T00:00:01.0000000Z', 'seconds', 1n],
    ['1970-01-01T00:00:00.001000000Z', 'milliseconds', 1n],
    ['1969-12-31T23:59:59Z', 'seconds', -1n],
    ['18446744073709551658', 'seconds', 18446744073709551658n],
    ['19700101T000001Z', 'seconds', 1n],
    ['1970-W01-4T00:00:01Z', 'seconds', 1n],
    ['1970-01-01T00:00:31+00:00:30', 'seconds', 1n],
  ] as const)
    expect(parse_time(stamp, unit)).toBe(expected);
  for (const [stamp, unit] of [
    ['1970-01-01T00:00:01.0000001Z', 'seconds'],
    ['1970-01-01T00:00:00.0010001Z', 'milliseconds'],
    ['1970-01-01T00:00:01,0000001Z', 'seconds'],
    ['1970-01-01T00:00:01+00:00:00.0000001', 'seconds'],
    ['1970-01-01T00:00:01+00:00:00.5', 'seconds'],
    ['1970-01-01T00:00:00', 'seconds'],
    ['1970-01-01T00:00:00.001Z', 'seconds'],
    ['2024-02-30T00:00:00Z', 'seconds'],
    ['1970-01-01T00:00:00+24:00', 'seconds'],
  ] as const)
    expect(() => parse_time(stamp, unit)).toThrow();
});

test('fresh seed properties: truncation, stream consumption, profiles, exact JSON and UTC units', () => {
  const replay = Number(process.env.WEAK_KEYS_SEED ?? randomInt(0, 0x7fffffff));
  console.info(`timestamp core generator=v1 replay seed=${replay}`);
  fc.assert(
    fc.property(
      fc.bigInt({ min: 0n, max: (1n << 96n) - 1n }),
      fc.integer({ min: 0, max: 15 }),
      fc.integer({ min: 0, max: 4 }),
      (seed, skip, index) => {
        expect(parse_json(json({ seed }))).toEqual({
          seed: seed <= BigInt(Number.MAX_SAFE_INTEGER) ? Number(seed) : seed,
        });
        for (const rng of RNGS) {
          const bits = rng === 'java-random' ? 48n : 32n;
          if (rng !== 'python-mt19937')
            expect(take(rng, seed, 20)).toEqual(take(rng, seed + (1n << bits), 20));
          for (const profile of PROFILES) {
            if (
              (profile.startsWith('python-') && rng !== 'python-mt19937') ||
              (profile.startsWith('java-') && rng !== 'java-random')
            )
              continue;
            const continuous = stream_at(rng, seed, profile, skip);
            let value: bigint | null = null;
            for (let i = 0; i <= index; i++)
              value = scalar_from_words(continuous, profile, CURVES.p256.order);
            expect(scalar(rng, seed, profile, CURVES.p256.order, skip, index)).toBe(value);
          }
        }
        const epoch = Number(seed % 2000000000n);
        const iso = new Date(epoch * 1000).toISOString();
        expect(parse_time(iso, 'seconds')).toBe(BigInt(epoch));
        expect(parse_time(iso, 'milliseconds')).toBe(BigInt(epoch) * 1000n);
      }
    ),
    { seed: replay, numRuns: 25 }
  );
});
