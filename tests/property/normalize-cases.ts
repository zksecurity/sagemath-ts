#!/usr/bin/env bun
/**
 * Rewrite property-test case files into the canonical compact form described
 * in `case-format.ts`, or verify that they already are.
 *
 * Usage:
 *   bun tests/property/normalize-cases.ts            # rewrite every cases/*.cases.json
 *   bun tests/property/normalize-cases.ts arith ...  # rewrite the named areas only
 *   bun tests/property/normalize-cases.ts --check    # exit 1 if any file is not canonical
 *
 * Run this after adding cases in the verbose `fixedValue` form; the
 * `case-format.test.ts` unit test fails until every checked-in file is canonical.
 */

import {
  closeSync,
  openSync,
  readdirSync,
  renameSync,
  writeFileSync,
  statSync,
  unlinkSync,
} from 'node:fs';
import { basename, join } from 'node:path';
import {
  type CaseSuite,
  countCases,
  isCanonicalText,
  normalizeSuite,
  serializeSuite,
} from './case-format.js';

import { readText, isCaseFile, caseArea, writeGzipChunks } from './storage.js';

const CASES_DIR = join(import.meta.dir, 'cases');

function caseFiles(areas: string[]): string[] {
  const all = readdirSync(CASES_DIR).filter(isCaseFile).sort();
  if (areas.length === 0) return all;
  const missing = areas.filter((a) => !all.some((f) => caseArea(f) === a));
  if (missing.length !== 0) throw new Error(`Unknown area(s): ${missing.join(', ')}`);
  return areas.map((a) => all.find((f) => caseArea(f) === a)!);
}

export async function normalizeFile(
  path: string
): Promise<{ before: number; after: number; cases: number }> {
  const text = readText(path);
  const before = statSync(path).size;
  const suite = normalizeSuite(JSON.parse(text) as CaseSuite);
  const expected = caseArea(basename(path));
  if (suite.module !== expected) {
    throw new Error(
      `${path}: "module" is ${JSON.stringify(suite.module)}, expected ${JSON.stringify(expected)}`
    );
  }
  const tmp = `${path}.tmp`;
  if (path.endsWith('.gz') || Buffer.byteLength(text) >= 1_000_000) {
    const output = path.endsWith('.gz') ? path : `${path}.gz`;
    await writeGzipChunks(tmp, serializeSuite(suite));
    renameSync(tmp, output);
    if (output !== path) unlinkSync(path);
    return { before, after: statSync(output).size, cases: countCases(suite) };
  }
  const fd = openSync(tmp, 'w');
  let after = 0;
  try {
    for (const chunk of serializeSuite(suite)) {
      writeFileSync(fd, chunk);
      after += Buffer.byteLength(chunk);
    }
  } finally {
    closeSync(fd);
  }
  renameSync(tmp, path);
  return { before, after, cases: countCases(suite) };
}

async function main(): Promise<void> {
  const argv = process.argv.slice(2);
  const check = argv.includes('--check');
  const files = caseFiles(argv.filter((a) => !a.startsWith('--')));
  let dirty = 0;
  for (const file of files) {
    const path = join(CASES_DIR, file);
    if (check) {
      if (!isCanonicalText(readText(path))) {
        dirty++;
        console.error(`not canonical: ${file}`);
      }
      continue;
    }
    const { before, after, cases } = await normalizeFile(path);
    const mb = (n: number) => (n / 1e6).toFixed(1);
    console.log(`${file}: ${mb(before)} MB -> ${mb(after)} MB (${cases} cases)`);
  }
  if (check) {
    if (dirty !== 0) {
      console.error(`${dirty} file(s) need \`bun tests/property/normalize-cases.ts\``);
      process.exit(1);
    }
    console.log(`${files.length} case file(s) canonical`);
  }
}

if (import.meta.main) await main();
