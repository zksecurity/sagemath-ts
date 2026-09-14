import { readFileSync } from 'node:fs';
import { join, relative, dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { runPythonTests } from './compare.js';
import { materializeSuite, freshSeed, derivedSeed, GENERATOR_VERSION } from './seeded.js';
const root = resolve(import.meta.dir, '../..');
const catalog = JSON.parse(readFileSync(join(import.meta.dir, 'native-suites.json'), 'utf8'));
const processSeed = process.env.SAGEMATH_TEST_SEED === undefined ? freshSeed() : Number(process.env.SAGEMATH_TEST_SEED);
const processRuns = process.env.SAGEMATH_TEST_RUNS === undefined ? 8 : Number(process.env.SAGEMATH_TEST_RUNS);
let announced = false;

/** Generate fresh payloads in constrained constructor domains, plus named regressions. */
export function nativeInputs(name, seed = processSeed, runs = processRuns) {
  materializeSuite({ module: 'validation', cases: [] }, seed, runs);
  const descriptor = catalog[name];
  if (!descriptor) throw new Error(`Unregistered native comparison suite: ${name}`);
  const suites = descriptor.sources.map((source, index) => materializeSuite({
    module: source.area,
    cases: [{ function: source.function, recipes: source.recipes }],
  }, derivedSeed(seed, `${name}:${index}`), runs));
  for (const regression of descriptor.regressions ?? []) suites.push({
    module: regression.area,
    cases: regression.recipe
      ? materializeSuite({ module: regression.area, cases: [{ function: regression.function, recipes: [regression.recipe] }] }, regression.seed, 1).cases
      : [{ function: regression.function, rows: [regression.row] }],
    regressionId: regression.id,
  });
  return suites;
}

/** Expected values are produced by the original implementation on every run. */
export async function nativeFixtures(testURL, formerFixture) {
  const name = relative(root, resolve(dirname(fileURLToPath(testURL)), formerFixture));
  if (!announced) {
    console.error(`Native property seed: ${processSeed}; ${GENERATOR_VERSION}; replay with SAGEMATH_TEST_SEED=${processSeed} SAGEMATH_TEST_RUNS=${processRuns} bun test <same test files>`);
    announced = true;
  }
  const suites = nativeInputs(name), descriptor = catalog[name], output = [];
  for (const suite of suites) {
    const result = await runPythonTests(JSON.stringify(suite));
    const rows = suite.cases[0].rows;
    if (result.length !== rows.length) throw new Error(`Native reference omitted results for ${name}`);
    for (let i = 0; i < result.length; i++) {
      const r = result[i], row = rows[i];
      const args = row.slice(1).map(a => Array.isArray(a) ? `[${a.join(', ')}]` : String(a));
      if (r.function !== suite.cases[0].function || r.seed !== row[0] || JSON.stringify(r.args) !== JSON.stringify(args))
        throw new Error(`Native reference returned a result for the wrong input in ${name}`);
      if (r.errorType === undefined) throw new Error('Native reference omitted exception class');
      if (descriptor.successOnly && (r.result === null || r.error !== null || r.errorType !== null))
        throw new Error(`Generated native success-domain input errors for ${name}: ${r.errorType}: ${r.error}`);
      if (/^(Unknown module|Unknown function):/.test(r.error ?? ''))
        throw new Error(`Native reference dispatch failed for ${name}: ${r.error}`);
      if (suite.regressionId !== undefined) r.regressionId = suite.regressionId;
      output.push(r);
    }
  }
  return output;
}

/** Adapt the original small-fixture test shapes to live generated native results. */
export async function legacyNativeFixtures(testURL, suiteID) {
  const rows = await nativeFixtures(testURL, suiteID);
  return rows.map((row, index) => ({
    probe: row.regressionId === undefined ? index : Number(row.regressionId),
    function: row.function,
    kind: row.function === 'pari_word_linear' ? 'word'
      : row.function === 'pari_integer_linear' ? 'integer'
      : row.function === 'pari_lll_heuristic' ? 'heuristic' : 'proved',
    group: row.regressionId === '162' ? 'extra' : row.function === 'pari_lll_wrapper' ? 'main' : 'stage',
    args: row.args.map(a => a.startsWith('[')
      ? a.slice(1, -1).split(',').map(v => v.trim()).filter(Boolean)
      : Number.isSafeInteger(Number(a)) ? Number(a) : a),
    ...(row.error === null ? { result: row.result, expected: JSON.parse(row.result) }
      : { error: row.error, errorType: row.errorType }),
  }));
}
