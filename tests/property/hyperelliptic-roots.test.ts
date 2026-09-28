import { expect, test } from 'bun:test';
import { readFileSync } from 'node:fs';
import { compareResults, runPythonTests, runTypeScriptTests } from './compare.js';
import { freshSeed, materializeSuite } from './seeded.js';
import type { CaseSuite } from './case-format.js';

test('hyperelliptic root callers match bundled Sage dependencies', async () => {
  const seed = Number(process.env.SAGEMATH_TEST_SEED ?? freshSeed());
  const runs = Number(process.env.SAGEMATH_TEST_RUNS ?? 200);
  const source = JSON.parse(
    readFileSync(new URL('./cases/hyperelliptic.cases.json', import.meta.url), 'utf8')
  ) as CaseSuite;
  source.cases = source.cases.filter((c) => c.function === 'hyp_root_callers');
  const suite = materializeSuite(source, seed, runs),
    input = JSON.stringify(suite);
  console.log(
    `Property seed: ${seed}; replay: SAGEMATH_TEST_SEED=${seed} SAGEMATH_TEST_RUNS=${runs} bun test tests/property/hyperelliptic-roots.test.ts`
  );
  const result = compareResults(await runPythonTests(input), await runTypeScriptTests(input));
  console.log(
    `Total: ${result.total}, Passed: ${result.passed}, Failed: ${result.failed}, Errors: ${result.errors}`
  );
  if (result.failed || result.errors)
    throw new Error(JSON.stringify(result.results.filter((r) => !r.match).slice(0, 3)));
  expect(result.passed).toBe(suite.cases.reduce((n, c) => n + (c.rows?.length ?? 0), 0));
}, 120_000);
