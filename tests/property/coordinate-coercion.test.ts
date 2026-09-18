import { expect, test } from 'bun:test';
import { readFileSync } from 'node:fs';
import { compareResults, runPythonTests, runTypeScriptTests } from './compare.js';
import { freshSeed, materializeSuite } from './seeded.js';
import type { CaseSuite } from './case-format.js';

test('curve coordinate conversion and GF(2) constructors match the original implementation', async () => {
  const seed =
    process.env.SAGEMATH_TEST_SEED === undefined
      ? freshSeed()
      : Number(process.env.SAGEMATH_TEST_SEED);
  const runs = Number(process.env.SAGEMATH_TEST_RUNS ?? 300);
  const source = JSON.parse(
    readFileSync(new URL('./cases/ec_advanced.cases.json', import.meta.url), 'utf8')
  ) as CaseSuite;
  source.cases = source.cases.filter((c) => c.function === 'ec_coordinate_coercion');
  expect(source.cases.length).toBeGreaterThan(0);
  const input = JSON.stringify(materializeSuite(source, seed, runs));
  console.log(
    `Property seed: ${seed}; replay: SAGEMATH_TEST_SEED=${seed} SAGEMATH_TEST_RUNS=${runs} bun test tests/property/coordinate-coercion.test.ts`
  );
  const original = await runPythonTests(input);
  const port = await runTypeScriptTests(input);
  const result = compareResults(original, port);
  console.log(
    `Total: ${result.total}, Passed: ${result.passed}, Failed: ${result.failed}, Errors: ${result.errors}`
  );
  if (result.failed || result.errors) {
    throw new Error(JSON.stringify(result.results.filter((r) => !r.match).slice(0, 3)));
  }
  expect(result.passed).toBeGreaterThanOrEqual(runs + 9);
}, 120_000);
