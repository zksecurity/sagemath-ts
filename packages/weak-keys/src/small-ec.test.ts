import { expect, test } from 'bun:test';
import { spawnSync } from 'node:child_process';
import { randomInt } from 'node:crypto';
import fc from 'fast-check';
import { p256, p384, secp256k1 } from '../../zksecurity-cheatsheets/src/curves/params.js';
import { buildSmallScalarDictionary, hasSmallEcScalar, pointKey } from './small-ec.js';

const curves = [p256, p384, secp256k1] as const;
const nativeNames = ['prime256v1', 'secp384r1', 'secp256k1'] as const;
function nativePoint(index: number, scalar: bigint) {
  const script = `const {createECDH}=require('node:crypto');const e=createECDH(process.argv[1]);e.setPrivateKey(Buffer.from(process.argv[2],'hex'));process.stdout.write(e.getPublicKey('hex','uncompressed'));`;
  const width = index === 1 ? 96 : 64;
  const run = spawnSync(
    'node',
    ['-e', script, nativeNames[index]!, scalar.toString(16).padStart(width, '0')],
    { encoding: 'utf8' }
  );
  if (run.error || run.status !== 0) throw new Error('Independent native EC oracle unavailable');
  const point = run.stdout;
  expect(point.length).toBe(2 + 2 * width);
  expect(point.slice(0, 2)).toBe('04');
  return {
    x: BigInt(`0x${point.slice(2, 2 + width)}`),
    y: BigInt(`0x${point.slice(2 + width)}`),
  };
}

for (let index = 0; index < curves.length; index++) {
  const curve = curves[index]!;
  test(`${curve.name}: exact inclusive signed boundary and outside range`, () => {
    const dictionary = buildSmallScalarDictionary(curve, 65536);
    expect(dictionary.size).toBe(131072);
    for (const scalar of [1n, 65536n]) {
      expect(dictionary.has(pointKey(nativePoint(index, scalar)))).toBe(true);
      expect(dictionary.has(pointKey(nativePoint(index, curve.order - scalar)))).toBe(true);
    }
    expect(dictionary.has(pointKey(nativePoint(index, 65537n)))).toBe(false);
    expect(dictionary.has(pointKey(nativePoint(index, curve.order - 65537n)))).toBe(false);
  });
}

test('fresh signed native points agree across curves and bounded intervals', () => {
  const seed = process.env.WEAK_KEYS_SEED ? Number(process.env.WEAK_KEYS_SEED) : randomInt(2 ** 31);
  console.info(`small-ec generator=v1 replay seed=${seed}`);
  fc.assert(
    fc.property(
      fc.integer({ min: 0, max: 2 }),
      fc.integer({ min: 1, max: 16 }),
      fc.integer({ min: 1, max: 131072 }),
      fc.boolean(),
      (index, bits, scalar, negated) => {
        const curve = curves[index]!;
        const point = nativePoint(index, negated ? curve.order - BigInt(scalar) : BigInt(scalar));
        expect(hasSmallEcScalar(curve, point, 2 ** bits)).toBe(scalar <= 2 ** bits);
      }
    ),
    { seed, numRuns: 32 }
  );
});

test('cache is keyed by exact scope and all curve parameters', () => {
  const point = nativePoint(0, 37n);
  expect(hasSmallEcScalar(p256, point, 64)).toBe(true);
  expect(hasSmallEcScalar(p256, point, 32)).toBe(false);
  expect(hasSmallEcScalar(p256, point, 37)).toBe(true);
  expect(hasSmallEcScalar(p384, nativePoint(1, 37n), 32)).toBe(false);
  expect(hasSmallEcScalar(p256, point, 37)).toBe(true);
  const modified = { ...p256, generator: nativePoint(0, 2n) };
  expect(hasSmallEcScalar(modified, point, 37)).toBe(false);
  expect(hasSmallEcScalar(p256, point, 37)).toBe(true);
});

test('non-power-of-two dictionaries and invalid scopes', () => {
  expect(buildSmallScalarDictionary(p256, 37).size).toBe(74);
  for (const bound of [0, -1, 1.5, 65537, Number.NaN, Number.POSITIVE_INFINITY])
    expect(() => buildSmallScalarDictionary(p256, bound)).toThrow();
  expect(() => buildSmallScalarDictionary({ ...p256, order: 100n }, 50)).toThrow();
  expect(() => buildSmallScalarDictionary({ ...p256, generator: { x: 0n, y: 0n } }, 4)).toThrow();
});
