import { expect, test } from 'bun:test';
import { readFileSync } from 'node:fs';
import { compareResults, runPythonTests, runTypeScriptTests } from './compare.js';
import { freshSeed, materializeSuite } from './seeded.js';
import type { CaseSuite } from './case-format.js';

test('curve coordinate coercion, promotion and base change match the original implementation', async () => {
  const seed =
    process.env.SAGEMATH_TEST_SEED === undefined
      ? freshSeed()
      : Number(process.env.SAGEMATH_TEST_SEED);
  const runs = Number(process.env.SAGEMATH_TEST_RUNS ?? 300);
  const source = JSON.parse(
    readFileSync(new URL('./cases/ec_advanced.cases.json', import.meta.url), 'utf8')
  ) as CaseSuite;
  source.cases = source.cases.filter((c) =>
    ['ec_coordinate_coercion', 'ec_curve_base_change', 'ec_lift_extension', 'ec_finite_coordinates', 'ec_finite_coordinate_roots'].includes(c.function)
  );
  expect(source.cases.length).toBeGreaterThan(0);
  const suite = materializeSuite(source, seed, runs);
  const expected = suite.cases.reduce((n, c) => n + (c.rows?.length ?? 0), 0);
  const input = JSON.stringify(suite);
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
  expect(result.passed).toBe(expected);
}, 120_000);

test('number-field polynomial variable names and scalar zero methods match Sage', async () => {
  const seed =
    process.env.SAGEMATH_TEST_SEED === undefined
      ? freshSeed()
      : Number(process.env.SAGEMATH_TEST_SEED);
  const runs = Number(process.env.SAGEMATH_TEST_RUNS ?? 100);
  const source = JSON.parse(
    readFileSync(new URL('./cases/number_fields.cases.json', import.meta.url), 'utf8')
  ) as CaseSuite;
  source.cases = source.cases.filter((c) => c.function === 'nf_named_polynomial');
  const suite = materializeSuite(source, seed, runs);
  const input = JSON.stringify(suite);
  console.log(
    `Property seed: ${seed}; replay: SAGEMATH_TEST_SEED=${seed} SAGEMATH_TEST_RUNS=${runs} bun test tests/property/coordinate-coercion.test.ts`
  );
  const result = compareResults(await runPythonTests(input), await runTypeScriptTests(input));
  console.log(
    `Total: ${result.total}, Passed: ${result.passed}, Failed: ${result.failed}, Errors: ${result.errors}`
  );
  if (result.failed || result.errors)
    throw new Error(JSON.stringify(result.results.filter((r) => !r.match).slice(0, 3)));
  expect(result.passed).toBe(runs + 3);
}, 120_000);
