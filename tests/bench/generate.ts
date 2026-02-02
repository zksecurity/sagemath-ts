#!/usr/bin/env bun
/**
 * Benchmark input generator.
 *
 * Reads benchmark case definitions from tests/bench/cases/*.bench.json
 * and writes deterministic inputs to tests/bench/inputs/*.bench.inputs.json.
 */

import { mkdirSync, readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { basename, join, relative } from 'node:path';
import { is_prime, next_prime } from '../../packages/sagemath-ts/src/index.js';
import { MersenneTwister } from '../property/typescript/mersenne-twister.js';

interface BenchCaseDefinition {
  id?: string;
  function: string;
  seeds?: number[];
  samples?: number;
  argGenerators?: string[];
  warmup?: number;
  iterations?: number;
  repeats?: number;
}

interface BenchSuiteDefinition {
  module: string;
  cases: BenchCaseDefinition[];
}

interface BenchCaseInputs {
  id: string;
  function: string;
  inputs: (string | string[])[][];
  warmup: number;
  iterations: number;
  repeats: number;
}

interface BenchSuiteInputs {
  module: string;
  generatedAt: string;
  sourceFile: string;
  cases: BenchCaseInputs[];
}

const DEFAULTS = {
  samples: 50,
  warmup: 1000,
  iterations: 10000,
  repeats: 5,
} as const;

class SeededRandom {
  private mt: MersenneTwister;

  constructor(seed: number) {
    this.mt = new MersenneTwister(seed);
  }

  integer(min: bigint, max: bigint): bigint {
    return this.mt.randint(min, max);
  }

  prime(min: bigint, max: bigint): bigint {
    const maxAttempts = 10000;
    for (let i = 0; i < maxAttempts; i++) {
      const candidate = this.integer(min, max);
      if (is_prime(candidate)) {
        return candidate;
      }
    }
    return next_prime(min - 1n);
  }
}

function generateArg(generatorSpec: string, random: SeededRandom): bigint | bigint[] {
  if (generatorSpec.startsWith('randomBigint(')) {
    const params = generatorSpec.slice(13, -1);
    const parts = params.split(',');
    const minVal = BigInt(parts[0]!.trim());
    const maxVal = BigInt(parts[1]!.trim());
    return random.integer(minVal, maxVal);
  }

  if (generatorSpec.startsWith('randomPrime(')) {
    const params = generatorSpec.slice(12, -1);
    const parts = params.split(',');
    const minVal = BigInt(parts[0]!.trim());
    const maxVal = BigInt(parts[1]!.trim());
    return random.prime(minVal, maxVal);
  }

  if (generatorSpec.startsWith('randomBigintBits(')) {
    const params = generatorSpec.slice(17, -1);
    const bits = Number.parseInt(params.trim(), 10);
    if (!Number.isFinite(bits) || bits < 1) {
      throw new Error(`invalid bit size for randomBigintBits: ${params}`);
    }
    const minVal = 1n << BigInt(bits - 1);
    const maxVal = (1n << BigInt(bits)) - 1n;
    return random.integer(minVal, maxVal);
  }

  if (generatorSpec.startsWith('fixedValue(')) {
    const valueStr = generatorSpec.slice(11, -1);
    return BigInt(valueStr.trim());
  }

  if (generatorSpec.startsWith('randomList(')) {
    const params = generatorSpec.slice(11, -1);
    const lastComma = params.lastIndexOf(',');
    const innerGenerator = params.slice(0, lastComma).trim();
    const length = Number.parseInt(params.slice(lastComma + 1).trim(), 10);
    const result: bigint[] = [];
    for (let i = 0; i < length; i++) {
      const val = generateArg(innerGenerator, random);
      if (Array.isArray(val)) {
        result.push(...val);
      } else {
        result.push(val);
      }
    }
    return result;
  }

  return BigInt(generatorSpec);
}

function formatArg(arg: bigint | bigint[]): string | string[] {
  if (Array.isArray(arg)) {
    return arg.map((value) => value.toString());
  }
  return arg.toString();
}

function parseArgs() {
  const args = process.argv.slice(2);
  let caseFilter: string | null = null;
  let casesDir = join(process.cwd(), 'tests/bench/cases');
  let outDir = join(process.cwd(), 'tests/bench/inputs');

  for (let i = 0; i < args.length; i++) {
    if (args[i] === '--case' && i + 1 < args.length) {
      caseFilter = args[i + 1] || null;
      i++;
      continue;
    }
    if (args[i] === '--cases-dir' && i + 1 < args.length) {
      casesDir = args[i + 1]!;
      i++;
      continue;
    }
    if (args[i] === '--out-dir' && i + 1 < args.length) {
      outDir = args[i + 1]!;
      i++;
      continue;
    }
  }

  return { caseFilter, casesDir, outDir };
}

function resolveCaseFiles(casesDir: string, caseFilter: string | null): string[] {
  const all = readdirSync(casesDir)
    .filter((name) => name.endsWith('.bench.json'))
    .map((name) => join(casesDir, name));

  if (!caseFilter) {
    return all;
  }

  const normalized = caseFilter.endsWith('.bench.json') ? caseFilter : `${caseFilter}.bench.json`;
  const matches = all.filter(
    (file) => basename(file) === normalized || file.endsWith(`/${normalized}`)
  );
  if (matches.length === 0) {
    throw new Error(`No benchmark case file matched: ${caseFilter}`);
  }
  return matches;
}

function generateInputsForSuite(suite: BenchSuiteDefinition, sourceFile: string): BenchSuiteInputs {
  const moduleName = suite.module || 'arith';
  const casesOut: BenchCaseInputs[] = [];

  suite.cases.forEach((benchCase, index) => {
    const seeds = benchCase.seeds ?? [42];
    const samples = benchCase.samples ?? DEFAULTS.samples;
    const argGenerators = benchCase.argGenerators ?? [];
    const warmup = benchCase.warmup ?? DEFAULTS.warmup;
    const iterations = benchCase.iterations ?? DEFAULTS.iterations;
    const repeats = benchCase.repeats ?? DEFAULTS.repeats;
    const inputs: (string | string[])[][] = [];

    for (const seed of seeds) {
      const random = new SeededRandom(seed);
      for (let s = 0; s < samples; s++) {
        const args = argGenerators.map((gen) => generateArg(gen, random));
        inputs.push(args.map((arg) => formatArg(arg)) as (string | string[])[]);
      }
    }

    casesOut.push({
      id: benchCase.id ?? `${benchCase.function}-${index + 1}`,
      function: benchCase.function,
      inputs,
      warmup,
      iterations,
      repeats,
    });
  });

  return {
    module: moduleName,
    generatedAt: new Date().toISOString(),
    sourceFile: relative(process.cwd(), sourceFile),
    cases: casesOut,
  };
}

function main() {
  const { caseFilter, casesDir, outDir } = parseArgs();
  const caseFiles = resolveCaseFiles(casesDir, caseFilter);

  mkdirSync(outDir, { recursive: true });

  for (const file of caseFiles) {
    const raw = readFileSync(file, 'utf8');
    const suite = JSON.parse(raw) as BenchSuiteDefinition;
    const output = generateInputsForSuite(suite, file);
    const outName = basename(file).replace('.bench.json', '.bench.inputs.json');
    const outPath = join(outDir, outName);
    writeFileSync(outPath, JSON.stringify(output, null, 2));
    process.stdout.write(`Generated ${outPath}\n`);
  }
}

main();
