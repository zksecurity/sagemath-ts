import { expect, test } from 'bun:test';
import { readFileSync } from 'node:fs';
import { compareResults, runPythonTests, runTypeScriptTests } from './compare.js';
import { freshSeed, materializeSuite } from './seeded.js';
import type { CaseSuite } from './case-format.js';

test('elliptic torsion callers match bundled Sage', async () => {
  const seed = Number(process.env.SAGEMATH_TEST_SEED ?? freshSeed());
  const runs = Number(process.env.SAGEMATH_TEST_RUNS ?? 60);
  const source = JSON.parse(
    readFileSync(new URL('./cases/ec_advanced.cases.json', import.meta.url), 'utf8')
  ) as CaseSuite;
  source.cases = source.cases.filter((c) =>
    ['ec_primary_torsion', 'ec_division_points', 'ec_scalar_order'].includes(c.function)
  );
  const suite = materializeSuite(source, seed, runs),
    input = JSON.stringify(suite);
  console.log(
    `Property seed: ${seed}; replay: SAGEMATH_TEST_SEED=${seed} SAGEMATH_TEST_RUNS=${runs} bun test tests/property/elliptic-torsion.test.ts`
  );
  const result = compareResults(await runPythonTests(input), await runTypeScriptTests(input));
  console.log(
    `Total: ${result.total}, Passed: ${result.passed}, Failed: ${result.failed}, Errors: ${result.errors}`
  );
  if (result.failed || result.errors)
    throw new Error(JSON.stringify(result.results.filter((r) => !r.match).slice(0, 3)));
  expect(result.passed).toBe(suite.cases.reduce((n, c) => n + (c.rows?.length ?? 0), 0));
}, 120_000);

test('torsion linear-relation dependency matches bundled Sage', async () => {
  const seed = Number(process.env.SAGEMATH_TEST_SEED ?? freshSeed());
  const runs = Number(process.env.SAGEMATH_TEST_RUNS ?? 120);
  const source = JSON.parse(
    readFileSync(new URL('./cases/groups_modn.cases.json', import.meta.url), 'utf8')
  ) as CaseSuite;
  source.cases = source.cases.filter((c) =>
    ['gg_linear_relation', 'im_additive_order'].includes(c.function)
  );
  const suite = materializeSuite(source, seed, runs),
    input = JSON.stringify(suite);
  console.log(`Property seed: ${seed}; runs: ${runs}; linear_relation`);
  const result = compareResults(await runPythonTests(input), await runTypeScriptTests(input));
  console.log(
    `Total: ${result.total}, Passed: ${result.passed}, Failed: ${result.failed}, Errors: ${result.errors}`
  );
  if (result.failed || result.errors)
    throw new Error(JSON.stringify(result.results.filter((r) => !r.match).slice(0, 3)));
  expect(result.passed).toBe(suite.cases.reduce((n, c) => n + (c.rows?.length ?? 0), 0));
}, 120_000);

test('division polynomial zero and negative indices match Sage', async () => {
  const source = JSON.parse(
    readFileSync(new URL('./cases/elliptic_curves.cases.json', import.meta.url), 'utf8')
  ) as CaseSuite;
  source.cases = source.cases.filter(
    (c) => c.function === 'division_polynomial' && Array.isArray(c.rows)
  );
  const input = JSON.stringify(source);
  const result = compareResults(await runPythonTests(input), await runTypeScriptTests(input));
  if (result.failed || result.errors)
    throw new Error(JSON.stringify(result.results.filter((r) => !r.match).slice(0, 3)));
  expect(result.passed).toBe(4);
}, 120_000);
