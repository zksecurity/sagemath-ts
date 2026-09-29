import { expect, test } from 'bun:test';
import { readFileSync } from 'node:fs';
import { compareResults, runPythonTests, runTypeScriptTests } from './compare.js';
import { freshSeed, materializeSuite } from './seeded.js';
import type { CaseSuite } from './case-format.js';

test('general curve PARI models match native callers', async () => {
  const seed = Number(process.env.SAGEMATH_TEST_SEED ?? freshSeed());
  const runs = Number(process.env.SAGEMATH_TEST_RUNS ?? 60);
  const source = JSON.parse(
    readFileSync(new URL('./cases/ec_advanced.cases.json', import.meta.url), 'utf8')
  ) as CaseSuite;
  source.cases = source.cases.filter((c) =>
    [
      'ec_pari_scalar',
      'ec_pari_order',
      'pari_elliptic_order',
      'pari_elliptic_model',
      'ec_pari_transformed',
      'pari_elliptic_transformed',
    ].includes(c.function)
  );
  const suite = materializeSuite(source, seed, runs),
    input = JSON.stringify(suite);
  console.log(
    `Property seed: ${seed}; replay: SAGEMATH_TEST_SEED=${seed} SAGEMATH_TEST_RUNS=${runs} bun test tests/property/elliptic-pari-model.test.ts`
  );
  const result = compareResults(await runPythonTests(input), await runTypeScriptTests(input));
  console.log(
    `Total: ${result.total}, Passed: ${result.passed}, Failed: ${result.failed}, Errors: ${result.errors}`
  );
  if (result.failed || result.errors)
    throw new Error(JSON.stringify(result.results.filter((r) => !r.match).slice(0, 3)));
  expect(result.passed).toBe(suite.cases.reduce((n, c) => n + (c.rows?.length ?? 0), 0));
}, 120_000);
