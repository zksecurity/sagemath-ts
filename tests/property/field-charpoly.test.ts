import { expect, test } from 'bun:test';
import { readFileSync } from 'node:fs';
import { compareResults, runPythonTests, runTypeScriptTests } from './compare.js';
import { freshSeed, materializeSuite } from './seeded.js';
import type { CaseSuite } from './case-format.js';

test('finite-field characteristic polynomials and bivariate resultants match PARI', async () => {
  const seed = Number(process.env.SAGEMATH_TEST_SEED ?? freshSeed());
  const runs = Number(process.env.SAGEMATH_TEST_RUNS ?? 300);
  const source = JSON.parse(
    readFileSync(new URL('./cases/finite_fields.cases.json', import.meta.url), 'utf8')
  ) as CaseSuite;
  source.cases = source.cases.filter((c) =>
    [
      'pari_field_charpoly',
      'pari_field_charpoly_scalar',
      'pari_bivariate_resultant',
      'ff_extension_norm',
    ].includes(c.function)
  );
  const suite = materializeSuite(source, seed, runs),
    input = JSON.stringify(suite);
  console.log(
    `Property seed: ${seed}; replay: SAGEMATH_TEST_SEED=${seed} SAGEMATH_TEST_RUNS=${runs} bun test tests/property/field-charpoly.test.ts`
  );
  const result = compareResults(await runPythonTests(input), await runTypeScriptTests(input));
  console.log(
    `Total: ${result.total}, Passed: ${result.passed}, Failed: ${result.failed}, Errors: ${result.errors}`
  );
  if (result.failed || result.errors)
    throw new Error(JSON.stringify(result.results.filter((r) => !r.match).slice(0, 3)));
  expect(result.passed).toBe(suite.cases.reduce((n, c) => n + (c.rows?.length ?? 0), 0));
}, 120_000);
