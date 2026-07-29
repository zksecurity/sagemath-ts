#!/usr/bin/env bun
/**
 * Stable fast/slow partitions for the unit and live differential suites.
 *
 * The slow lists are intentionally explicit. New tests and new property areas
 * enter the fast tier by default, so adding coverage cannot silently make the
 * routinely-run tier smaller.
 */

import { basename, join } from 'node:path';

type Suite = 'unit' | 'property';
type Tier = 'fast' | 'slow';

const projectRoot = join(import.meta.dir, '..');

const slowUnitFiles = new Set([
  'packages/parigp-ts/src/buch.test.ts',
  'packages/parigp-ts/src/elliptic/advanced.test.ts',
  'packages/parigp-ts/src/elliptic/ellsea.test.ts',
  'packages/parigp-ts/src/galconj.test.ts',
  'packages/parigp-ts/src/mpqs.test.ts',
  'packages/parigp-ts/src/polmodular.test.ts',
  'packages/parigp-ts/src/qfb.test.ts',
  'packages/sagemath-ts/src/coding/bch_code.test.ts',
  'packages/sagemath-ts/src/coding/goppa_code.test.ts',
  'packages/sagemath-ts/src/coding/reed_muller_code.test.ts',
  'packages/sagemath-ts/src/coding/reed_solomon.test.ts',
  'packages/sagemath-ts/src/rings/number_field/galois_group.test.ts',
  'packages/sagemath-ts/src/rings/number_field/number_field_embeddings.test.ts',
  'packages/sagemath-ts/src/rings/number_field/number_field_ideal.test.ts',
  'packages/sagemath-ts/src/rings/polynomial/polynomial_factorization.test.ts',
  'packages/sagemath-ts/src/schemes/elliptic_curves/ell_finite_field.test.ts',
  'packages/sagemath-ts/src/schemes/elliptic_curves/ell_point.test.ts',
  'packages/sagemath-ts/src/schemes/elliptic_curves/ell_torsion.test.ts',
  'packages/sagemath-ts/src/schemes/elliptic_curves/isogeny_class.test.ts',
  'packages/sagemath-ts/src/schemes/hyperelliptic_curves/hyperelliptic_finite_field.test.ts',
  'packages/sagemath-ts/src/schemes/hyperelliptic_curves/jacobian_morphism.test.ts',
  'packages/sagemath-ts/src/stats/distributions/discrete_gaussian_integer.test.ts',
  'packages/sagemath-ts/src/stats/distributions/discrete_gaussian_lattice.test.ts',
]);

const slowPropertyAreas = new Set([
  'coding_crypto',
  'ec_advanced',
  'groups_modn',
  'lattices',
  'matrix_extended',
  'mpfr',
  'padics_series',
  'rand_stats',
]);

function selected<T>(all: T[], slow: Set<T>, tier: Tier): T[] {
  return all.filter((item) => (tier === 'slow' ? slow.has(item) : !slow.has(item)));
}

function assertPartition<T>(all: T[], slow: Set<T>, label: string): void {
  const actual = new Set(all);
  const missing = [...slow].filter((item) => !actual.has(item));
  if (missing.length !== 0) {
    throw new Error(`Unknown ${label} in slow tier: ${missing.join(', ')}`);
  }
  const fast = selected(all, slow, 'fast');
  const slowSelected = selected(all, slow, 'slow');
  if (fast.length + slowSelected.length !== all.length) {
    throw new Error(`${label} tier partition is not exhaustive`);
  }
}

async function unitFiles(): Promise<string[]> {
  const patterns = ['packages/**/*.test.ts', 'tests/**/*.test.ts', 'tutorial/**/*.test.ts'];
  const files = patterns
    .flatMap((pattern) => [
      ...new Bun.Glob(pattern).scanSync({ cwd: projectRoot, onlyFiles: true }),
    ])
    .sort();
  assertPartition(files, slowUnitFiles, 'unit test file');
  return files;
}

async function propertyAreas(): Promise<string[]> {
  const glob = new Bun.Glob('tests/property/cases/*.cases.json');
  const areas = [...glob.scanSync({ cwd: projectRoot, onlyFiles: true })]
    .map((path) => basename(path, '.cases.json'))
    .sort();
  assertPartition(areas, slowPropertyAreas, 'property area');
  return areas;
}

async function main(): Promise<void> {
  const [suiteArg, tierArg, ...forwarded] = process.argv.slice(2);
  if (
    (suiteArg !== 'unit' && suiteArg !== 'property') ||
    (tierArg !== 'fast' && tierArg !== 'slow')
  ) {
    console.error(
      'usage: bun tests/run-test-tier.ts <unit|property> <fast|slow> [--list|suite options]'
    );
    process.exit(2);
  }

  const suite = suiteArg as Suite;
  const tier = tierArg as Tier;
  const all = suite === 'unit' ? await unitFiles() : await propertyAreas();
  const slow = suite === 'unit' ? slowUnitFiles : slowPropertyAreas;
  const chosen = selected(all, slow, tier);

  if (forwarded.includes('--list')) {
    console.log(chosen.join('\n'));
    return;
  }

  console.log(`Running ${chosen.length}/${all.length} ${suite} ${tier}-tier entries`);
  const command =
    suite === 'unit'
      ? ['bun', 'test', ...chosen, ...forwarded]
      : ['bun', 'tests/property/compare.ts', '--cases', chosen.join(','), ...forwarded];
  const child = Bun.spawn(command, {
    cwd: projectRoot,
    stdin: 'inherit',
    stdout: 'inherit',
    stderr: 'inherit',
  });
  process.exit(await child.exited);
}

await main();
