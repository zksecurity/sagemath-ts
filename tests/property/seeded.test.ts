import { expect, test } from 'bun:test';
import { materializeSuite, generateArgument, generateRecipe, derivedSeed } from './seeded.js';
import { MersenneTwister } from './typescript/mersenne-twister.js';

const suite = { module: 'demo', cases: [
  { function: 'f', argGenerators: ['randomBigint(-100000000000000000000, 100000000000000000000)', 'randomList(randomBigint(-100, 100), 5)'] },
  { function: 'regression', rows: [[42, '9007199254740993']] as [number, string][] },
] };

test('replay regenerates identical inputs and another seed explores new inputs', () => {
  const a = materializeSuite(suite, 123, 10);
  expect(a).toEqual(materializeSuite(suite, 123, 10));
  expect(a).not.toEqual(materializeSuite(suite, 124, 10));
  expect(a.cases[1]).toEqual(suite.cases[1]);
  expect(JSON.parse(JSON.stringify(a))).toEqual(a);
  expect(derivedSeed(123, 'a')).not.toBe(derivedSeed(123, 'b'));
});

test('input generation rejects invalid domains instead of silently fabricating a prime', () => {
  const rng = new MersenneTwister(42);
  expect(() => generateArgument('randomPrime(14, 16)', rng)).toThrow('No prime');
  expect(() => generateArgument('randomBigint(2, 1)', rng)).toThrow('Empty');
  expect(() => generateArgument('randomList(randomBigint(0, 2), -1)', rng)).toThrow('Invalid');
  expect(() => materializeSuite(suite, NaN, 10)).toThrow('seed');
  expect(() => materializeSuite(suite, 1, 0)).toThrow('runs');
});

test('recipes preserve matrix shape and positive scaling, byte bounds, and integer precision', () => {
  const recipe = [
    { scale: { values: ['9007199254740993', 0, 0, '9007199254740993'], range: [1, 7] as [number, number] } },
    { bytes: 32 }, { list: { length: 6, range: [-8, 8] as [number, number] } },
    { integer: ['9007199254740993', '9007199254741093'] as [string, string] },
    { text: 'expression' as const },
  ];
  for (const seed of [1, 123, 456, 0xffffffff]) {
    const [matrix, bytes, coefficients, big, text] = generateRecipe(recipe, new MersenneTwister(seed));
    expect(matrix).toHaveLength(4);
    const m = matrix as bigint[];
    expect(m[0]).toBe(m[3]); expect(m[1]).toBe(0n); expect(m[2]).toBe(0n);
    expect(m[0]! / 9007199254740993n).toBeGreaterThanOrEqual(1n);
    expect(m[0]! % 9007199254740993n).toBe(0n);
    expect(bytes).toHaveLength(32);
    expect((bytes as bigint[]).every(v => v >= 0n && v <= 255n)).toBe(true);
    expect(coefficients).toHaveLength(6);
    expect((coefficients as bigint[]).every(v => v >= -8n && v <= 8n)).toBe(true);
    expect(big as bigint).toBeGreaterThanOrEqual(9007199254740993n);
    expect(String.fromCharCode(...(text as bigint[]).map(Number))).toMatch(/^\(\d+ \+ \d+\) \* \d+$/);
  }
  expect(() => generateRecipe([{ bytes: -1 }], new MersenneTwister(1))).toThrow('length');
  expect(() => generateRecipe([{ integer: [2, 1] }], new MersenneTwister(1))).toThrow('interval');
});

test('modular polynomial divisor recipes always have a nonzero leading coefficient', () => {
  for (let seed = 0; seed < 100; seed++) {
    const [divisor] = generateRecipe([{ list: { length: 3, range: [0, 6], nonzeroLast: true } }], new MersenneTwister(seed));
    expect((divisor as bigint[])[2]).toBeGreaterThan(0n);
    expect((divisor as bigint[])[2]).toBeLessThan(7n);
  }
  expect(() => generateRecipe([{ list: { length: 1, range: [0, 0], nonzeroLast: true } }], new MersenneTwister(1))).toThrow('Nonzero');
});

test('retained failure seeds replay identical payloads alongside fresh trials without duplicating rows', () => {
  const original = { module: 'demo', cases: [{ function: 'f', argGenerators: ['randomBigint(0, 1000000000)'] }] };
  const failing = materializeSuite(original, 123, 3).cases[0] as { rows: unknown[] };
  const retained = { module: 'demo', cases: [{ ...original.cases[0]!, replaySeeds: [{ seed: 123, runs: 3, reason: 'regression' }] }] };
  for (const seed of [123, 456, 789]) {
    const rows = (materializeSuite(retained, seed, 3).cases[0] as { rows: unknown[] }).rows;
    for (const row of failing.rows) expect(rows).toContainEqual(row);
    expect(rows.length).toBe(seed === 123 ? 3 : 6);
  }
});

test('large-integer regression recipes reconstruct shifted values without stored decimal blobs', () => {
  const [values] = generateRecipe([{ shifted: [
    { coefficient: 1, shift: 50100, offset: 0 },
    { coefficient: 1, shift: 50097, offset: 3 },
  ] }], new MersenneTwister(1));
  expect(values).toEqual([1n << 50100n, (1n << 50097n) + 3n]);
  expect(() => generateRecipe([{ shifted: [{ coefficient: 1, shift: -1, offset: 0 }] }], new MersenneTwister(1))).toThrow('shift');
});
