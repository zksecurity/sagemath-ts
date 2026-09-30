import { expect, test } from 'bun:test';
import { readFileSync } from 'node:fs';
import { compareResults, runPythonTests, runTypeScriptTests } from './compare.js';
import { freshSeed, materializeSuite } from './seeded.js';
import type { CaseSuite } from './case-format.js';
import { F3xq_ellcardj, setrand } from '@sagemath-ts/parigp-ts';

test('ternary count corrects the bundled PARI random twist-sign bug', async () => {
  const T=[1,2,0,1],a=[0,0,2],b=[2,1,2];
  const input=JSON.stringify({module:'ec_advanced',cases:[{function:'pari_word_extension_card',rows:
    [[1,3,T,a,b,1],[2,3,T,a,b,4]]}]});
  const native=await runPythonTests(input);
  expect(native).toHaveLength(2);
  expect(native.every(r=>r.error===null)).toBe(true);
  // Independent enumeration is checked below. Upstream returns the twist at seed 4.
  expect(native.map(r=>JSON.parse(r.result!).value[0])).toEqual(['19','37']);
  for(const seed of [1n,4n]) {
    setrand(seed);
    expect(F3xq_ellcardj(a.map(BigInt),b.map(BigInt),T.map(BigInt),27n,3)).toBe(19n);
  }
},30_000);

test('ternary supersingular counts match independent enumeration and native RNG', async () => {
  const seed=Number(process.env.SAGEMATH_TEST_SEED??freshSeed());
  const runs=Number(process.env.SAGEMATH_TEST_RUNS??60);
  const source=JSON.parse(readFileSync(new URL('./cases/ec_advanced.cases.json',import.meta.url),'utf8')) as CaseSuite;
  source.cases=source.cases.filter(c=>c.function==='pari_ternary_supersingular_card');
  const suite=materializeSuite(source,seed,runs),input=JSON.stringify(suite);
  console.log(`Property seed: ${seed}; replay: SAGEMATH_TEST_SEED=${seed} SAGEMATH_TEST_RUNS=${runs} bun test tests/property/elliptic-ternary-card.test.ts`);
  const result=compareResults(await runPythonTests(input),await runTypeScriptTests(input));
  console.log(`Total: ${result.total}, Passed: ${result.passed}, Failed: ${result.failed}, Errors: ${result.errors}`);
  if(result.failed||result.errors)throw new Error(JSON.stringify(result.results.filter(r=>!r.match).slice(0,3)));
  expect(result.passed).toBe(suite.cases.reduce((n,c)=>n+(c.rows?.length??0),0));
},120_000);
