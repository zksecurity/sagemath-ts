#!/usr/bin/env bun
/**
 * Property test runner for sagemath-ts.
 * Executes operations with deterministic seeding and outputs results.
 *
 * This script reads test case definitions from stdin (JSON format) and
 * executes them using our TypeScript implementation, outputting results
 * in a format that can be compared with the SageMath implementation.
 *
 * Usage:
 *   bun runner.ts < test_cases.json
 *   bun runner.ts --seed 42 < test_cases.json
 *
 * Output format:
 *   {"function": "gcd", "args": [12, 8], "result": "4", "seed": 42}
 *
 * Areas
 * -----
 * This runner contains NO area-specific code. Every area lives in its own module
 * under `areas/` and is discovered by globbing that directory, so several agents
 * can add areas in parallel without ever editing this file. See
 * tests/property/README.md.
 */

import { existsSync, readdirSync } from 'node:fs';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import {
  type Factorization,
  formatFactorization,
  is_prime,
  next_prime,
} from '../../../packages/sagemath-ts/src/index.js';
import { MersenneTwister } from './mersenne-twister.js';

/** Directory holding the per-area dispatch modules. */
const AREAS_DIR = join(import.meta.dir, 'areas');

/**
 * Seeded pseudo-random number generator using Mersenne Twister.
 * This matches Python's random module exactly for the same seed.
 */
class SeededRandom {
  private mt: MersenneTwister;

  constructor(seed: number) {
    this.mt = new MersenneTwister(seed);
  }

  /**
   * Generate a random integer in range [min, max] inclusive
   */
  integer(min: bigint, max: bigint): bigint {
    return this.mt.randint(min, max);
  }

  /**
   * Generate a random prime in range [min, max]
   */
  prime(min: bigint, max: bigint): bigint {
    // Simple approach: generate random numbers until we find a prime
    // This is not efficient but works for testing with small ranges
    const maxAttempts = 10000;
    for (let i = 0; i < maxAttempts; i++) {
      const candidate = this.integer(min, max);
      if (is_prime(candidate)) {
        return candidate;
      }
    }
    // Fallback: find next prime from min
    const p = next_prime(min - 1n);
    while (p > max) {
      // If no prime in range, just return the next prime after min
      return p;
    }
    return p;
  }
}

/**
 * Test case definition from JSON
 */
interface TestCase {
  module?: string;
  function: string;
  seeds?: number[];
  argGenerators?: string[];
}

/**
 * Test suite definition from JSON
 */
interface TestSuite {
  module: string;
  cases: TestCase[];
}

/**
 * Test result output
 */
interface TestResult {
  function: string;
  args: string[];
  result: string | null;
  error: string | null;
  seed: number;
}

/**
 * A discovered area: its dispatch table plus any per-function formatter overrides.
 */
interface Area {
  functions: Record<string, (...args: never[]) => unknown>;
  formatters: Record<string, (result: unknown) => string>;
}

/**
 * Return the names of the area modules present in `areas/`.
 *
 * An area is any `areas/<area>.ts` whose basename does not start with `_`
 * (underscore-prefixed modules are shared helpers, not areas). The name must
 * match the `module` field of `cases/<area>.cases.json`.
 */
function discoverAreaNames(): string[] {
  if (!existsSync(AREAS_DIR)) {
    return [];
  }
  return readdirSync(AREAS_DIR)
    .filter(
      (f) =>
        f.endsWith('.ts') && !f.startsWith('_') && !f.endsWith('.d.ts') && !f.endsWith('.test.ts')
    )
    .map((f) => f.slice(0, -'.ts'.length))
    .sort();
}

/**
 * Import `areas/<name>.ts` and return its dispatch tables.
 *
 * Every area module must export a `functions` object mapping the `function`
 * names used in `cases/<name>.cases.json` to callables. It may optionally
 * export a `formatters` object mapping a function name to a
 * `(result) => string` callable, for results the generic `formatResult` below
 * cannot render.
 *
 * Areas are imported lazily and one at a time, so a broken area module can
 * never take down the areas owned by other agents.
 */
async function loadArea(name: string): Promise<Area> {
  const url = pathToFileURL(join(AREAS_DIR, `${name}.ts`)).href;
  const mod = (await import(url)) as {
    functions?: unknown;
    formatters?: unknown;
  };

  if (typeof mod.functions !== 'object' || mod.functions === null) {
    throw new Error(
      `Area module areas/${name}.ts must export a \`functions\` object mapping function names to callables`
    );
  }

  return {
    functions: mod.functions as Record<string, (...args: never[]) => unknown>,
    formatters: (mod.formatters ?? {}) as Record<string, (result: unknown) => string>,
  };
}

/**
 * Generate an argument based on the generator specification.
 */
function generateArg(generatorSpec: string, random: SeededRandom): bigint | bigint[] {
  if (generatorSpec.startsWith('randomBigint(')) {
    // Parse: randomBigint(min, max)
    const params = generatorSpec.slice(13, -1);
    const parts = params.split(',');
    const minVal = BigInt(parts[0]!.trim());
    const maxVal = BigInt(parts[1]!.trim());
    return random.integer(minVal, maxVal);
  }

  if (generatorSpec.startsWith('randomPrime(')) {
    // Parse: randomPrime(min, max)
    const params = generatorSpec.slice(12, -1);
    const parts = params.split(',');
    const minVal = BigInt(parts[0]!.trim());
    const maxVal = BigInt(parts[1]!.trim());
    return random.prime(minVal, maxVal);
  }

  if (generatorSpec.startsWith('fixedValue(')) {
    // Parse: fixedValue(value) - value can be integer or array
    const valueStr = generatorSpec.slice(11, -1).trim();
    if (valueStr.startsWith('[') && valueStr.endsWith(']')) {
      // Parse array: [1, 2, 3]
      const inner = valueStr.slice(1, -1).trim();
      if (!inner) {
        return [];
      }
      const parts = inner.split(',');
      return parts.map((p) => BigInt(p.trim()));
    }
    return BigInt(valueStr);
  }

  if (generatorSpec.startsWith('randomList(')) {
    // Parse: randomList(generator, length)
    const params = generatorSpec.slice(11, -1);
    const lastComma = params.lastIndexOf(',');
    const innerGenerator = params.slice(0, lastComma).trim();
    const length = Number.parseInt(params.slice(lastComma + 1).trim());
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

  // Assume it's a literal value
  return BigInt(generatorSpec);
}

/**
 * Format a result to match SageMath's string representation.
 *
 * This is the generic fallback used by every area. If your area produces a
 * shape this cannot render, either return an already-formatted string from the
 * area function or export a `formatters` entry from the area module — do not
 * add another special case here (it would make this file a merge hotspot).
 */
function formatResult(result: unknown, functionName?: string): string {
  // Handle null
  if (result === null) {
    return 'null';
  }

  // Special case: factor always returns a Factorization
  if (functionName === 'factor' && Array.isArray(result)) {
    return formatFactorization(result as Factorization);
  }

  // Special case: poly_factor returns [(coeffs, mult), ...]
  if (
    (functionName === 'poly_factor_ff' || functionName === 'poly_factor') &&
    Array.isArray(result)
  ) {
    // Format polynomial factorization as list of (factor_coeffs, multiplicity) pairs
    const formatted = result.map((item) => {
      if (Array.isArray(item) && item.length === 2 && Array.isArray(item[0])) {
        const [coeffs, mult] = item as [bigint[], number];
        const coeffsStr = '[' + coeffs.map((c) => c.toString()).join(', ') + ']';
        return `(${coeffsStr}, ${mult})`;
      }
      return formatResult(item);
    });
    return '[' + formatted.join(', ') + ']';
  }

  // Special case: poly_quo_rem returns (quotient, remainder) as tuple
  if (
    (functionName === 'poly_quo_rem_ff' || functionName === 'poly_quo_rem') &&
    Array.isArray(result)
  ) {
    if (result.length === 2 && Array.isArray(result[0]) && Array.isArray(result[1])) {
      const [q, r] = result as [bigint[], bigint[]];
      const qStr = '[' + q.map((c) => c.toString()).join(', ') + ']';
      const rStr = '[' + r.map((c) => c.toString()).join(', ') + ']';
      return `(${qStr}, ${rStr})`;
    }
  }

  if (typeof result === 'boolean') {
    return result ? 'True' : 'False';
  }

  if (Array.isArray(result)) {
    // Check if it's a 2D array (matrix) for HNF/LLL results
    if (
      result.length > 0 &&
      Array.isArray(result[0]) &&
      result[0].length > 0 &&
      (functionName?.startsWith('hnf_') || functionName?.startsWith('lll_'))
    ) {
      // It's a matrix result - format as list of lists
      const formattedRows = result.map(
        (row: Array<unknown>) => '[' + row.map((x) => String(x)).join(', ') + ']'
      );
      return '[' + formattedRows.join(', ') + ']';
    }

    // Check if it looks like a factorization (array of [prime, exp] tuples)
    // Must have length 2 inner arrays where both elements are bigint
    // But exclude matrix operations that return 2D arrays
    if (
      result.length > 0 &&
      Array.isArray(result[0]) &&
      result[0].length === 2 &&
      typeof result[0][0] === 'bigint' &&
      typeof result[0][1] === 'bigint' &&
      // Exclude matrix functions
      !functionName?.startsWith('hnf_') &&
      !functionName?.startsWith('lll_') &&
      !functionName?.startsWith('snf_')
    ) {
      // It's a Factorization
      return formatFactorization(result as Factorization);
    }

    // Check if it's a tuple-like result from xgcd (always 3 elements)
    // Only format as tuple if this is specifically xgcd function
    if (
      functionName === 'xgcd' &&
      result.length === 3 &&
      result.every((x) => typeof x === 'bigint')
    ) {
      return `(${result.map((x) => x!.toString()).join(', ')})`;
    }

    // Check if it's a signature result (tuple of 2 numbers)
    if (
      functionName === 'quadratic_signature' &&
      result.length === 2 &&
      result.every((x) => typeof x === 'number')
    ) {
      return `(${result[0]}, ${result[1]})`;
    }

    // Regular array
    return '[' + result.map((x) => formatResult(x)).join(', ') + ']';
  }

  if (typeof result === 'bigint') {
    return result.toString();
  }

  return String(result);
}

/**
 * Execute a function with the given arguments.
 */
function executeFunction(
  module: string,
  area: Area | undefined,
  functionName: string,
  args: (bigint | bigint[])[]
): unknown {
  if (area === undefined) {
    throw new Error(`Unknown module: ${module}`);
  }

  if (!(functionName in area.functions)) {
    throw new Error(`Unknown function: ${module}.${functionName}`);
  }

  const func = area.functions[functionName]! as (...a: unknown[]) => unknown;
  return func(...args);
}

/**
 * Run a single test case with the given seed.
 */
function runTestCase(
  testCase: TestCase,
  seed: number,
  module: string,
  area: Area | undefined
): TestResult {
  const random = new SeededRandom(seed);
  const functionName = testCase.function;
  const argGenerators = testCase.argGenerators || [];

  // Generate arguments
  const args: (bigint | bigint[])[] = [];
  for (const gen of argGenerators) {
    const arg = generateArg(gen, random);
    args.push(arg);
  }

  // Execute function
  let formattedResult: string | null = null;
  let error: string | null = null;

  try {
    const result = executeFunction(module, area, functionName, args);
    const override = area?.formatters[functionName];
    formattedResult = override ? override(result) : formatResult(result, functionName);
  } catch (e) {
    error = e instanceof Error ? e.message : String(e);
  }

  return {
    function: functionName,
    args: args.map((a) => (Array.isArray(a) ? `[${a.join(', ')}]` : a.toString())),
    result: formattedResult,
    error,
    seed,
  };
}

/**
 * Run all test cases in a test suite.
 */
async function runTestSuite(testSuite: TestSuite): Promise<TestResult[]> {
  const module = testSuite.module || 'arith';
  const cases = testSuite.cases || [];
  const results: TestResult[] = [];

  // Load only the area this suite needs, so a broken sibling area module can
  // never break an unrelated suite.
  const area = discoverAreaNames().includes(module) ? await loadArea(module) : undefined;

  for (const testCase of cases) {
    const seeds = testCase.seeds || [42];

    for (const seed of seeds) {
      const result = runTestCase(testCase, seed, module, area);
      results.push(result);
    }
  }

  return results;
}

/**
 * Main entry point
 */
async function main() {
  // Parse command line arguments
  let seedOverride: number | null = null;
  const args = process.argv.slice(2);

  for (let i = 0; i < args.length; i++) {
    if (args[i] === '--seed' && i + 1 < args.length) {
      seedOverride = Number.parseInt(args[i + 1]!);
      break;
    }
  }

  // Read test cases from stdin
  const inputChunks: string[] = [];
  const reader = Bun.stdin.stream().getReader();

  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    inputChunks.push(new TextDecoder().decode(value));
  }

  const inputData = inputChunks.join('');

  if (!inputData.trim()) {
    console.error('Error: No input provided. Pass test cases as JSON via stdin.');
    process.exit(1);
  }

  let testSuite: TestSuite;
  try {
    testSuite = JSON.parse(inputData);
  } catch (e) {
    console.error(`Error parsing JSON: ${e}`);
    process.exit(1);
  }

  // Override seeds if specified
  if (seedOverride !== null) {
    for (const testCase of testSuite.cases || []) {
      testCase.seeds = [seedOverride];
    }
  }

  // Run tests
  const results = await runTestSuite(testSuite);

  // Output results as JSON
  console.log(JSON.stringify(results, null, 2));
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
