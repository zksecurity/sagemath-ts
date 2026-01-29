#!/usr/bin/env bun
/**
 * Benchmark runner for sagemath-ts.
 *
 * Usage:
 *   bun tests/bench/typescript/runner.ts --input tests/bench/inputs/arith.bench.inputs.json
 *   cat tests/bench/inputs/arith.bench.inputs.json | bun tests/bench/typescript/runner.ts
 */

import {
  gcd,
  lcm,
  xgcd,
  factor,
  is_prime,
  is_prime_power,
  next_prime,
  previous_prime,
  euler_phi,
  radical,
  moebius,
  kronecker_symbol,
  legendre_symbol,
  jacobi_symbol,
  power_mod,
  inverse_mod,
  crt,
  isqrt,
  is_square,
  is_squarefree,
  divisors,
  number_of_divisors,
  sigma,
  prime_range,
  squarefree_part,
  prime_factors,
  valuation,
  arith,
} from '../../../packages/sagemath-ts/src/index.js';
import { readFileSync } from 'node:fs';
import {
  bls12_381_g1,
  createPrimeFieldCurve,
  p256,
  secp256k1,
} from '../../../packages/zksecurity-cheatsheets/src/index.js';

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
  sourceFile?: string;
  cases: BenchCaseInputs[];
}

interface BenchCaseResult {
  id: string;
  function: string;
  inputs: number;
  iterations: number;
  warmup: number;
  repeats: number;
  total_ns: string[];
  per_call_ns: number[];
  p50_ns: number | null;
  p90_ns: number | null;
  error: string | null;
}

interface BenchSuiteResult {
  module: string;
  inputFile: string | null;
  runtime: string;
  timestamp: string;
  results: BenchCaseResult[];
}

const secp256k1Context = createPrimeFieldCurve(secp256k1);
const p256Context = createPrimeFieldCurve(p256);
const bls12_381_g1Context = createPrimeFieldCurve(bls12_381_g1);

function parseArg(arg: string | string[]): bigint | bigint[] {
  if (Array.isArray(arg)) {
    return arg.map((value) => BigInt(value));
  }
  return BigInt(arg);
}

function percentile(values: number[], p: number): number | null {
  if (values.length === 0) {
    return null;
  }
  const sorted = [...values].sort((a, b) => a - b);
  const index = (sorted.length - 1) * p;
  const lower = Math.floor(index);
  const upper = Math.ceil(index);
  if (lower === upper) {
    return sorted[lower]!;
  }
  const weight = index - lower;
  return sorted[lower]! * (1 - weight) + sorted[upper]! * weight;
}

function executeFunction(
  module: string,
  functionName: string,
  args: (bigint | bigint[])[]
): unknown {
  const functionMap: Record<string, Record<string, (...args: any[]) => unknown>> = {
    arith: {
      gcd: (a: bigint, b: bigint) => gcd(a, b),
      lcm: (a: bigint, b: bigint) => lcm(a, b),
      xgcd: (a: bigint, b: bigint) => xgcd(a, b),
      factor: (n: bigint) => factor(n),
      is_prime: (n: bigint) => is_prime(n),
      is_prime_power: (n: bigint) => is_prime_power(n),
      next_prime: (n: bigint) => next_prime(n),
      previous_prime: (n: bigint) => previous_prime(n),
      euler_phi: (n: bigint) => euler_phi(n),
      radical: (n: bigint) => radical(n),
      moebius: (n: bigint) => moebius(n),
      kronecker_symbol: (a: bigint, n: bigint) => kronecker_symbol(a, n),
      legendre_symbol: (a: bigint, p: bigint) => legendre_symbol(a, p),
      jacobi_symbol: (a: bigint, n: bigint) => jacobi_symbol(a, n),
      power_mod: (a: bigint, n: bigint, m: bigint) => power_mod(a, n, m),
      inverse_mod: (a: bigint, m: bigint) => inverse_mod(a, m),
      crt: (a: bigint, b: bigint, m: bigint, n: bigint) => crt(a, b, m, n),
      isqrt: (n: bigint) => isqrt(n),
      is_square: (n: bigint) => is_square(n),
      is_squarefree: (n: bigint) => is_squarefree(n),
      divisors: (n: bigint) => divisors(n),
      number_of_divisors: (n: bigint) => number_of_divisors(n),
      sigma: (n: bigint, k?: bigint) => sigma(n, k),
      prime_range: (start: bigint, stop?: bigint) =>
        stop !== undefined ? prime_range(start, stop) : prime_range(start),
      trial_division: (n: bigint, bound?: bigint) =>
        bound !== undefined ? arith.trial_division(n, bound) : arith.trial_division(n),
      squarefree_part: (n: bigint) => squarefree_part(n),
      prime_factors: (n: bigint) => prime_factors(n),
      valuation: (n: bigint, p: bigint) => valuation(n, p),
    },
    bench_crypto: {
      rsa_2048_pow: (m: bigint, e: bigint, n: bigint) => power_mod(m % n, e, n),
      secp256k1_mul: (k: bigint) => secp256k1Context.generator.mul(k),
      p256_mul: (k: bigint) => p256Context.generator.mul(k),
      bls12_381_g1_mul: (k: bigint) => bls12_381_g1Context.generator.mul(k),
    },
  };

  if (!(module in functionMap)) {
    throw new Error(`Unknown module: ${module}`);
  }

  if (!(functionName in functionMap[module]!)) {
    throw new Error(`Unknown function: ${module}.${functionName}`);
  }

  const func = functionMap[module]![functionName]!;
  return func(...args);
}

function runCase(moduleName: string, benchCase: BenchCaseInputs): BenchCaseResult {
  const inputs = benchCase.inputs.map((args) => args.map((arg) => parseArg(arg)));
  const warmup = benchCase.warmup ?? 0;
  const iterations = benchCase.iterations ?? 1;
  const repeats = benchCase.repeats ?? 1;

  let error: string | null = null;
  let totalNs: string[] = [];
  let perCallNs: number[] = [];

  try {
    for (let i = 0; i < warmup; i++) {
      const args = inputs[i % inputs.length]!;
      executeFunction(moduleName, benchCase.function, args);
    }

    for (let r = 0; r < repeats; r++) {
      const start = process.hrtime.bigint();
      for (let i = 0; i < iterations; i++) {
        const args = inputs[i % inputs.length]!;
        executeFunction(moduleName, benchCase.function, args);
      }
      const end = process.hrtime.bigint();
      const elapsed = end - start;
      totalNs.push(elapsed.toString());
      perCallNs.push(Number(elapsed) / iterations);
    }
  } catch (e) {
    error = e instanceof Error ? e.message : String(e);
    totalNs = [];
    perCallNs = [];
  }

  return {
    id: benchCase.id,
    function: benchCase.function,
    inputs: inputs.length,
    iterations,
    warmup,
    repeats,
    total_ns: totalNs,
    per_call_ns: perCallNs,
    p50_ns: percentile(perCallNs, 0.5),
    p90_ns: percentile(perCallNs, 0.9),
    error,
  };
}

async function readInput(inputPath: string | null): Promise<string> {
  if (inputPath) {
    return readFileSync(inputPath, 'utf8');
  }

  const chunks: string[] = [];
  const reader = Bun.stdin.stream().getReader();

  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(new TextDecoder().decode(value));
  }

  return chunks.join('');
}

function parseArgs() {
  const args = process.argv.slice(2);
  let inputPath: string | null = null;

  for (let i = 0; i < args.length; i++) {
    if (args[i] === '--input' && i + 1 < args.length) {
      inputPath = args[i + 1]!;
      i++;
    }
  }

  return { inputPath };
}

async function main() {
  const { inputPath } = parseArgs();
  const inputData = await readInput(inputPath);

  if (!inputData.trim()) {
    console.error('Error: No input provided. Pass benchmark inputs via --input or stdin.');
    process.exit(1);
  }

  const suite = JSON.parse(inputData) as BenchSuiteInputs;
  const results = suite.cases.map((benchCase) => runCase(suite.module, benchCase));

  const output: BenchSuiteResult = {
    module: suite.module,
    inputFile: inputPath,
    runtime: `bun ${Bun.version}`,
    timestamp: new Date().toISOString(),
    results,
  };

  console.log(JSON.stringify(output, null, 2));
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
