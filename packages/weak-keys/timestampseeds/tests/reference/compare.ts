/** Fresh generated inputs, live native output, replay seed printed on every run. */
import assert from 'node:assert/strict';
import { randomInt } from 'node:crypto';
import { readFileSync } from 'node:fs';
import fc from 'fast-check';
import { PROFILES, candidate } from '../../../src/timestampseeds/core.js';
import { json, parse_json } from '../../../src/timestampseeds/json.js';
import { RNGS, type Rng, words } from '../../../src/timestampseeds/rngs.js';
import { build, public_reference } from './build_fixtures.js';
import { References } from './references.js';

const replay = Number(process.env.WEAK_KEYS_SEED ?? randomInt(0, 0x7fffffff));
console.log(`timestamp native generator=v1 replay seed=${replay}`);
const windows = process.argv.includes('--windows');
if (windows && process.platform !== 'win32')
  throw new Error('Microsoft CRT reference requires Windows');
const reference = new References(windows);
let comparisons = 0;
try {
  const names = windows
    ? (['msvc-rand'] as Rng[])
    : RNGS.filter((name) => process.platform === 'linux' || name !== 'glibc-rand');
  const seeds = [
    0n,
    1n,
    42n,
    1700000000n,
    0x80000000n,
    0xffffffffn,
    0x100000001n,
    (1n << 48n) + 42n,
    (1n << 64n) + 42n,
    ...fc.sample(fc.bigInt({ min: 0n, max: (1n << 96n) - 1n }), { seed: replay, numRuns: 8 }),
  ];
  for (const rng of names)
    for (const seed of seeds) {
      const count = rng.includes('mt19937') ? 1300 : 100;
      const stream = words(rng, seed);
      assert.deepEqual(
        Array.from({ length: count }, () => stream.next().value![0]),
        reference.words(rng, seed, count),
        `${rng} seed=${seed}`
      );
      comparisons += count;
    }
  if (!windows) {
    for (const [row, seed] of seeds.slice(-4).entries())
      for (const rng of names)
        for (const profile of PROFILES) {
          if (
            (profile.startsWith('python-') && rng !== 'python-mt19937') ||
            (profile.startsWith('java-') && rng !== 'java-random')
          )
            continue;
          const skip = row + 1,
            index = row % 3;
          const value = reference.value(rng, seed, profile, skip, index);
          for (const curve of ['p256', 'secp256k1'] as const) {
            const expected = public_reference(curve, value, profile);
            const actual = candidate(seed, 'seconds', rng, profile, curve, skip, index);
            assert.deepEqual(
              actual && {
                sec1_compressed_hex: actual.sec1_compressed_hex,
                spki_sha256: actual.spki_sha256,
              },
              expected,
              `${rng}/${profile}/${curve} seed=${seed}`
            );
            comparisons++;
          }
        }
    if (process.platform === 'linux') {
      const original = parse_json(
        readFileSync(new URL('../fixtures/public_keys.json', import.meta.url), 'utf8')
      ) as { cases: unknown[] };
      assert.equal(json(build(reference).cases), json(original.cases));
      comparisons += original.cases.length;
    }
  }
  console.log(
    `${comparisons} live native comparisons passed; ${windows ? 'Microsoft CRT' : process.platform === 'linux' ? 'all eight RNGs and original fixtures' : 'seven RNGs; glibc job runs on Linux'}.`
  );
} finally {
  reference.close();
}
