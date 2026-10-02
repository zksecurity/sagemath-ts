import { expect, test } from 'bun:test';
import { readFileSync } from 'node:fs';
import { compareResults, runPythonTests, runTypeScriptTests } from './compare.js';
import { freshSeed, materializeSuite } from './seeded.js';
import type { CaseSuite } from './case-format.js';

test('polynomial precision lifting matches bundled PARI', async () => {
  const seed = Number(process.env.SAGEMATH_TEST_SEED ?? freshSeed());
  // Precision recipes share one compact case entry; retain multiple fresh
  // visits per constructor variant despite removing duplicate case wrappers.
  const runs = Number(process.env.SAGEMATH_TEST_RUNS ?? 200);
  const source = JSON.parse(
    readFileSync(new URL('./cases/function_fields.cases.json', import.meta.url), 'utf8')
  ) as CaseSuite;
  source.cases = source.cases.filter((c) =>
    ['ff_pari_zp_precision', 'ff_pari_zp_lift', 'ff_pari_zp_binary_linear', 'ff_parent'].includes(c.function)
  );
  const suite = materializeSuite(source, seed, runs),
    input = JSON.stringify(suite);
  console.log(
    `Property seed: ${seed}; replay: SAGEMATH_TEST_SEED=${seed} SAGEMATH_TEST_RUNS=${runs} bun test tests/property/padic-precision.test.ts`
  );
  const result = compareResults(await runPythonTests(input), await runTypeScriptTests(input));
  console.log(
    `Total: ${result.total}, Passed: ${result.passed}, Failed: ${result.failed}, Errors: ${result.errors}`
  );
  if (result.failed || result.errors)
    throw new Error(JSON.stringify(result.results.filter((r) => !r.match).slice(0, 3)));
  expect(result.passed).toBe(suite.cases.reduce((n, c) => n + (c.rows?.length ?? 0), 0));
}, 180_000);


test('binary linear lift corrects the bundled PARI zero-exponent crash', async () => {
  const input=JSON.stringify({module:'function_fields',cases:[{function:'ff_pari_zp_precision',
    rows:[[268484677,21,[391],[1],[1],2,13]]}]});
  // The exact native body is guarded at its illegal generic-power call.
  const native=await runPythonTests(input);
  expect(native).toHaveLength(1);
  expect(native[0]!.errorType).toBe('PariError');
  expect(native[0]!.error).toBe('bug in zero exponent in Fp_pow2n, please report.');
  const { Zp_sqrtnlift }=await import('../../packages/parigp-ts/src/Zp.js');
  expect(Zp_sqrtnlift(783n,1n,1n,2n,13)).toBe(783n);
},30_000);
