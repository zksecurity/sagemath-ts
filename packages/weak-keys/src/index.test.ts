import { describe, expect, test } from 'bun:test';
import { spawnSync } from 'node:child_process';
import { randomInt } from 'node:crypto';
import { fileURLToPath } from 'node:url';
import fc from 'fast-check';
import { CHECKS, type Check, screenPublicKeys } from './index.js';

const rsa = (n: bigint | string, e: bigint | string = 3n) => ({ type: 'rsa', n, e });
const status = (input: unknown, check: Check, options = {}) =>
  screenPublicKeys([input], { checks: [check], ...options }).records[0]!.checks.find(
    (r) => r.check === check
  )!;
function ec(scalar: number, curve = 'p256') {
  const reference = fileURLToPath(new URL('../tests/reference/ec.cjs', import.meta.url));
  const native = spawnSync('node', [reference, curve, String(scalar)], { encoding: 'utf8' });
  expect(native.status).toBe(0);
  return { type: 'ec', curve, point: native.stdout };
}

describe('bounded public-only screening', () => {
  test('prime status is probable, composite outcome exact; no certified-prime claim', () => {
    expect(status(rsa(101n), 'rsa-prime-modulus')).toMatchObject({
      status: 'matched',
      confidence: 'probable',
    });
    for (const n of [561n, 1105n, 341550071728321n])
      expect(status(rsa(n), 'rsa-prime-modulus')).toMatchObject({
        status: 'no-match',
        confidence: 'exact',
      });
  });
  test('general perfect powers include composite bases and odd powers', () => {
    for (const n of [9n, 27n, 64n, 216n, 81n, 15625n])
      expect(status(rsa(n), 'rsa-perfect-power').status).toBe('matched');
    for (const n of [15n, 21n, 35n, 101n])
      expect(status(rsa(n), 'rsa-perfect-power').status).toBe('no-match');
  });
  test('Fermat is opt-in, verifies factors, and never serializes them', () => {
    const input = {
      ...rsa(10007n * 10009n, 65537n),
      label: 'CONFIDENTIAL_SENTINEL',
      p: '10007',
      pattern: { word: 'SECRET_WORD' },
    };
    const report = screenPublicKeys([input], { checks: ['rsa-fermat'], fermatIterations: 1 });
    expect(report.records[0]!.checks.find((r) => r.check === 'rsa-fermat')).toEqual({
      check: 'rsa-fermat',
      status: 'matched',
      confidence: 'exact',
      scope: { iterations: 1 },
    });
    const serialized = JSON.stringify(report);
    for (const forbidden of [
      '10007',
      '10009',
      '100160063',
      'CONFIDENTIAL_SENTINEL',
      'SECRET_WORD',
      'pattern',
      'message',
      '"scalar":',
      'fingerprint',
    ])
      expect(serialized).not.toContain(forbidden);
    expect(
      screenPublicKeys([input]).records[0]!.checks.find((r) => r.check === 'rsa-fermat')!.status
    ).toBe('not-run');
  });
  test('declared checked scope differs from empty or exhausted budgets', () => {
    expect(status(rsa(101n), 'rsa-fermat', { fermatIterations: 1 })).toMatchObject({
      status: 'no-match',
      confidence: 'bounded',
      scope: { iterations: 1 },
    });
    expect(status(rsa(101n), 'rsa-fermat', { fermatIterations: 0 }).status).toBe(
      'budget-exhausted'
    );
    expect(status(rsa(101n), 'rsa-fermat', { maxFermatSteps: 0 }).status).toBe('budget-exhausted');
    const report = screenPublicKeys([rsa(101n), rsa(103n)], {
      checks: ['rsa-fermat'],
      fermatIterations: 2,
      maxFermatSteps: 3,
    });
    expect(report.records.map((r) => r.checks[2]!.scope)).toEqual([
      { iterations: 2 },
      { iterations: 1 },
    ]);
    expect(status(rsa(101n), 'rsa-prime-modulus', { maxKeys: 0 }).status).toBe('budget-exhausted');
  });
  test('malformed and unsupported records do not reach inappropriate math', () => {
    for (const input of [null, [], {}, rsa(0n), rsa(1n), rsa(10n, 2n), rsa('SECRET'), rsa(7n, 7n)])
      expect(status(input, 'rsa-prime-modulus').status).toBe('malformed');
    expect(status(rsa(257n), 'rsa-prime-modulus', { maxModulusBits: 8 }).status).toBe(
      'unsupported'
    );
    expect(status({ type: 'ed25519' }, 'rsa-prime-modulus').status).toBe('unsupported');
    expect(
      status({ type: 'ec', curve: 'UNTRUSTED_CURVE', point: 'SECRET' }, 'ec-small-scalar').status
    ).toBe('unsupported');
    expect(status({ type: 'ec', curve: 'p256', point: '00' }, 'ec-small-scalar').status).toBe(
      'malformed'
    );
  });
  test('all numeric budgets are finite, integral and hard-capped', () => {
    for (const maxKeys of [-1, Number.NaN, Number.POSITIVE_INFINITY, 0.5, 257])
      expect(() => screenPublicKeys([], { maxKeys })).toThrow();
    expect(() => screenPublicKeys([], { checks: ['unknown' as Check] })).toThrow();
    expect(() => screenPublicKeys([], { fermatIterations: 4097 })).toThrow();
    expect(() => screenPublicKeys([], { ecTableBits: 17 })).toThrow();
    expect(() => screenPublicKeys([], { ecTableBits: 0 })).toThrow();
  });
  test('native OpenSSL points match bounded EC screens; cached table cannot expand range', () => {
    const input = ec(271);
    const found = status(input, 'ec-small-scalar', { ecTableBits: 10 });
    expect(found).toMatchObject({
      status: 'matched',
      confidence: 'exact',
      scope: { maxScalar: 1024 },
    });
    expect(status(input, 'ec-small-scalar', { ecTableBits: 8 })).toMatchObject({
      status: 'no-match',
      scope: { maxScalar: 256 },
    });
    expect(status(ec(256), 'ec-small-scalar', { ecTableBits: 8 }).status).toBe('matched');
    expect(status(ec(257), 'ec-small-scalar', { ecTableBits: 8 }).status).toBe('no-match');
    const positive = ec(17);
    const field = 2n ** 256n - 2n ** 224n + 2n ** 192n + 2n ** 96n - 1n;
    const y = BigInt('0x' + positive.point.slice(66));
    const negative = {
      ...positive,
      point: positive.point.slice(0, 66) + (field - y).toString(16).padStart(64, '0'),
    };
    expect(status(negative, 'ec-small-scalar').status).toBe('matched');
    const offCurve = {
      ...positive,
      point: positive.point.slice(0, 66) + (y + 1n).toString(16).padStart(64, '0'),
    };
    expect(status(offCurve, 'ec-small-scalar').status).toBe('malformed');
    expect(status(ec(17, 'secp256k1'), 'ec-small-scalar').status).toBe('matched');
    expect(status(ec(17, 'p384'), 'ec-small-scalar').status).toBe('matched');
    const serialized = JSON.stringify(
      screenPublicKeys([input], { checks: ['ec-small-scalar'], ecTableBits: 10 })
    );
    expect(serialized).not.toContain('271');
    expect(serialized).not.toContain(input.point);
    expect(serialized).not.toContain('pattern');
    expect(serialized).not.toContain('"scalar":');
    expect(status(input, 'ec-small-scalar', { maxEcKeys: 0 }).status).toBe('budget-exhausted');
  });
  test('reports have a strict allowlist, even with attacker-controlled metadata', () => {
    fc.assert(
      fc.property(fc.jsonValue(), (metadata) => {
        const report = screenPublicKeys([{ ...rsa(101n), metadata }], { checks: CHECKS });
        expect(Object.keys(report).sort()).toEqual(['records', 'schemaVersion']);
        for (const record of report.records) {
          expect(Object.keys(record).sort()).toEqual(['checks', 'index']);
          for (const check of record.checks)
            expect(Object.keys(check).sort()).toEqual(
              check.scope
                ? ['check', 'confidence', 'scope', 'status']
                : ['check', 'confidence', 'status']
            );
        }
      }),
      { numRuns: 50 }
    );
  });
});

test('fresh live TypeScript properties: prime and general perfect-power adapters', () => {
  const seed = randomInt(0, 0x7fffffff);
  console.info(`weak-keys live TypeScript replay seed=${seed}`);
  const reference = fileURLToPath(new URL('../tests/reference/structure.ts', import.meta.url));
  fc.assert(
    fc.property(
      fc.array(fc.integer({ min: 5, max: 10000 }), { minLength: 40, maxLength: 60 }),
      (values) => {
        const native = spawnSync(process.execPath, [reference], {
          input: JSON.stringify(values),
          encoding: 'utf8',
        });
        expect(native.status).toBe(0); // Missing native reference is a failure, never a silent skip.
        const expected = JSON.parse(native.stdout) as Array<{ prime: boolean; power: boolean }>;
        const report = screenPublicKeys(
          values.map((n) => rsa(BigInt(n))),
          { maxKeys: 64 }
        );
        report.records.forEach((record, i) => {
          expect(record.checks[0]!.status === 'matched').toBe(expected[i]!.prime);
          expect(record.checks[1]!.status === 'matched').toBe(expected[i]!.power);
        });
      }
    ),
    { seed, numRuns: 5 }
  );
});
