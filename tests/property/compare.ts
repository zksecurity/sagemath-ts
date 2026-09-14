#!/usr/bin/env bun
import { materializeSuite, freshSeed, GENERATOR_VERSION } from './seeded.js';
import { writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { createHash } from 'node:crypto';
import { readText, isCaseFile, caseArea } from './storage.js';
/**
 * Live differential runner. Generate shared inputs from a fresh seed, execute
 * both original and port, and retain failing inputs only. Use --seed for an
 * exact generator replay or --replay with the failure artifact printed by the run.
 */

import { spawn } from 'node:child_process';
import { existsSync } from 'node:fs';
import { mkdtemp, readdir } from 'node:fs/promises';
import { join } from 'node:path';

// Colors for terminal output
const colors = {
  reset: '\x1b[0m',
  red: '\x1b[31m',
  green: '\x1b[32m',
  yellow: '\x1b[33m',
  blue: '\x1b[34m',
  magenta: '\x1b[35m',
  cyan: '\x1b[36m',
  gray: '\x1b[90m',
};

/**
 * Test result from either Python or TypeScript runner
 */
export interface TestResult {
  function: string;
  args: string[];
  result: string | null;
  error: string | null;
  errorType?: string | null;
  seed: number;
}

/**
 * Comparison result for a single test case
 */
interface ComparisonResult {
  function: string;
  args: string[];
  seed: number;
  pythonResult: string | null;
  typescriptResult: string | null;
  pythonError: string | null;
  typescriptError: string | null;
  match: boolean;
}

/**
 * Summary of comparison results
 */
interface ComparisonSummary {
  total: number;
  passed: number;
  failed: number;
  errors: number;
  results: ComparisonResult[];
}

// Paths
const SCRIPT_DIR = import.meta.dir;
const PROJECT_ROOT = join(SCRIPT_DIR, '../..');
const CASES_DIR = join(SCRIPT_DIR, 'cases');
const PYTHON_RUNNER = join(SCRIPT_DIR, 'python/runner.py');
const TS_RUNNER = join(SCRIPT_DIR, 'typescript/runner.ts');
const PYTHON_AREAS_DIR = join(SCRIPT_DIR, 'python/areas');
const TS_AREAS_DIR = join(SCRIPT_DIR, 'typescript/areas');

/**
 * Run a command and capture output
 */
async function runCommand(
  command: string,
  args: string[],
  input?: string
): Promise<{ stdout: string; stderr: string; exitCode: number }> {
  const timeoutMs = Number(process.env.SAGEMATH_TEST_TIMEOUT_MS ?? 120_000);
  if (!Number.isSafeInteger(timeoutMs) || timeoutMs < 1 || timeoutMs > 2_147_483_647)
    throw new Error('SAGEMATH_TEST_TIMEOUT_MS must be a positive timer-safe integer');
  return new Promise((resolve) => {
    const proc = spawn(command, args, {
      cwd: PROJECT_ROOT,
      detached: process.platform !== 'win32',
      stdio: ['pipe', 'pipe', 'pipe'],
    });

    const timer = setTimeout(() => {
      stderr += `\nReference/test process timed out after ${timeoutMs}ms`;
      try {
        if (process.platform !== 'win32' && proc.pid) process.kill(-proc.pid, 'SIGKILL');
        else proc.kill('SIGKILL');
      } catch {} // The process may already have exited.
    }, timeoutMs);
    let stdout = '';
    let stderr = '';

    proc.stdout.on('data', (data) => {
      stdout += data.toString();
    });

    proc.stderr.on('data', (data) => {
      stderr += data.toString();
    });

    proc.stdin.on('error', () => {}); // Early process exits are reported below.
    proc.stdin.end(input);

    proc.on('close', (code) => {
      clearTimeout(timer);
      resolve({
        stdout,
        stderr,
        exitCode: code ?? 1,
      });
    });

    proc.on('error', (err) => {
      clearTimeout(timer);
      resolve({
        stdout,
        stderr: err.message,
        exitCode: 1,
      });
    });
  });
}

/**
 * Run Python/SageMath tests
 */
export async function runPythonTests(testCaseJson: string): Promise<TestResult[]> {
  const result = await runCommand('sage', [PYTHON_RUNNER], testCaseJson);

  if (result.exitCode !== 0) {
    console.error(`${colors.red}SageMath runner failed:${colors.reset}\n${result.stderr}`);
    throw new Error('SageMath runner failed');
  }

  try {
    return JSON.parse(result.stdout);
  } catch (e) {
    console.error(`${colors.red}Failed to parse SageMath output:${colors.reset}\n${result.stdout}`);
    throw e;
  }
}

/**
 * Run TypeScript tests
 */
export async function runTypeScriptTests(testCaseJson: string): Promise<TestResult[]> {
  const result = await runCommand('bun', ['run', TS_RUNNER], testCaseJson);

  if (result.exitCode !== 0) {
    console.error(`${colors.red}TypeScript runner failed:${colors.reset}\n${result.stderr}`);
    throw new Error('TypeScript runner failed');
  }

  try {
    return JSON.parse(result.stdout);
  } catch (e) {
    console.error(
      `${colors.red}Failed to parse TypeScript output:${colors.reset}\n${result.stdout}`
    );
    throw e;
  }
}

/**
 * Compare results from Python and TypeScript
 */
export function compareResults(
  pythonResults: TestResult[],
  tsResults: TestResult[]
): ComparisonSummary {
  const results: ComparisonResult[] = [];
  let passed = 0;
  let failed = 0;
  let errors = 0;

  // Create a map for easier lookup
  const tsMap = new Map<string, TestResult[]>();
  for (const r of tsResults) {
    const key = JSON.stringify([r.function, r.args, r.seed]);
    const bucket = tsMap.get(key) ?? [];
    bucket.push(r);
    tsMap.set(key, bucket);
  }

  for (const pyResult of pythonResults) {
    const key = JSON.stringify([pyResult.function, pyResult.args, pyResult.seed]);
    const tsResult = tsMap.get(key)?.shift();

    if (!tsResult) {
      results.push({
        function: pyResult.function,
        args: pyResult.args,
        seed: pyResult.seed,
        pythonResult: pyResult.result,
        typescriptResult: null,
        pythonError: pyResult.error,
        typescriptError: 'Test not found in TypeScript results',
        match: false,
      });
      errors++;
      continue;
    }

    // Empty messages (notably AssertionError) still represent failures.
    const pyHasError = pyResult.error !== null || pyResult.errorType != null;
    const tsHasError = tsResult.error !== null || tsResult.errorType != null;
    if (pyHasError && tsHasError) {
      const dispatchFailure =
        /^(Unknown module|Unknown function):/.test(pyResult.error ?? '') ||
        /^(Unknown module|Unknown function):/.test(tsResult.error ?? '');
      results.push({
        function: pyResult.function,
        args: pyResult.args,
        seed: pyResult.seed,
        pythonResult: pyResult.result,
        typescriptResult: tsResult.result,
        pythonError: pyResult.error,
        typescriptError: tsResult.error,
        match:
          !dispatchFailure &&
          !!pyResult.errorType &&
          pyResult.errorType === tsResult.errorType &&
          pyResult.error === tsResult.error,
      });
      if (dispatchFailure) {
        errors++;
      } else if (!pyResult.errorType || !tsResult.errorType) {
        // Legacy transcripts have no exception class: regenerate them.
        errors++;
      } else if (pyResult.errorType === tsResult.errorType && pyResult.error === tsResult.error) {
        passed++;
      } else {
        failed++;
      }
      continue;
    }

    // One has error, one doesn't
    if (pyHasError || tsHasError) {
      results.push({
        function: pyResult.function,
        args: pyResult.args,
        seed: pyResult.seed,
        pythonResult: pyResult.result,
        typescriptResult: tsResult.result,
        pythonError: pyResult.error,
        typescriptError: tsResult.error,
        match: false,
      });
      failed++;
      continue;
    }

    // Compare results
    const match = pyResult.result === tsResult.result;
    results.push({
      function: pyResult.function,
      args: pyResult.args,
      seed: pyResult.seed,
      pythonResult: pyResult.result,
      typescriptResult: tsResult.result,
      pythonError: null,
      typescriptError: null,
      match,
    });

    if (match) {
      passed++;
    } else {
      failed++;
    }
  }

  for (const bucket of tsMap.values()) for (const extra of bucket) {
    errors++;
    results.push({ function: extra.function, args: extra.args, seed: extra.seed,
      pythonResult: null, typescriptResult: extra.result,
      pythonError: 'Unexpected TypeScript result without a reference input',
      typescriptError: extra.error, match: false });
  }
  if (!pythonResults.length && !tsResults.length) errors++;

  return {
    total: passed + failed + errors,
    passed,
    failed,
    errors,
    results,
  };
}

/**
 * Print comparison summary
 */
function printSummary(caseName: string, summary: ComparisonSummary, verbose: boolean): void {
  console.log(`\n${colors.blue}=== ${caseName} ===${colors.reset}`);
  console.log(
    `Total: ${summary.total}, ` +
      `${colors.green}Passed: ${summary.passed}${colors.reset}, ` +
      `${colors.red}Failed: ${summary.failed}${colors.reset}, ` +
      `${colors.yellow}Errors: ${summary.errors}${colors.reset}`
  );

  // Show failed/error cases
  const failures = summary.results.filter((r) => !r.match);

  if (failures.length > 0 || verbose) {
    console.log('');

    for (const result of summary.results) {
      if (result.match && !verbose) continue;

      const status = result.match
        ? `${colors.green}PASS${colors.reset}`
        : `${colors.red}FAIL${colors.reset}`;

      console.log(
        `  ${status} ${result.function}(${result.args.join(', ')}) [seed=${result.seed}]`
      );

      if (!result.match || verbose) {
        if (result.pythonError) {
          console.log(`    ${colors.gray}Python error: ${result.pythonError}${colors.reset}`);
        } else {
          console.log(`    ${colors.gray}Python:     ${result.pythonResult}${colors.reset}`);
        }

        if (result.typescriptError) {
          console.log(
            `    ${colors.gray}TypeScript error: ${result.typescriptError}${colors.reset}`
          );
        } else {
          console.log(`    ${colors.gray}TypeScript: ${result.typescriptResult}${colors.reset}`);
        }
      }
    }
  }
}

/**
 * Check if SageMath is available
 */
async function checkSageMath(): Promise<boolean> {
  try {
    const result = await runCommand('sage', ['--version']);
    return result.exitCode === 0;
  } catch {
    return false;
  }
}

/**
 * Load test case files
 */
async function loadTestCases(
  caseFilters?: Set<string>
): Promise<{ name: string; content: string }[]> {
  const cases: { name: string; content: string }[] = [];

  if (!existsSync(CASES_DIR)) {
    console.error(`${colors.red}Test cases directory not found: ${CASES_DIR}${colors.reset}`);
    return cases;
  }

  const files = await readdir(CASES_DIR);

  for (const file of files) {
    if (!isCaseFile(file)) continue;

    const name = caseArea(file);

    if (caseFilters && !caseFilters.has(name)) continue;
    if (
      !existsSync(join(PYTHON_AREAS_DIR, `${name}.py`)) ||
      !existsSync(join(TS_AREAS_DIR, `${name}.ts`))
    ) {
      throw new Error(
        `Property area '${name}' needs both python/areas/${name}.py and typescript/areas/${name}.ts`
      );
    }

    const content = readText(join(CASES_DIR, file));
    cases.push({ name, content });
  }

  return cases;
}

/** Live comparisons retain only failing inputs, never successful transcripts. */
async function main(): Promise<void> {
  const argv = process.argv.slice(2);
  let seed = freshSeed(), runs = 25, verbose = false;
  let caseFilters: Set<string> | undefined;
  let replay: string | undefined;
  for (let i = 0; i < argv.length; i++) {
    const value = () => {
      if (!argv[i + 1]) throw new Error(`Missing value for ${argv[i]}`);
      return argv[++i]!;
    };
    switch (argv[i]) {
      case '--case': caseFilters = new Set([value()]); break;
      case '--cases': caseFilters = new Set(value().split(',')); break;
      case '--seed': seed = Number(value()); break;
      case '--runs': runs = Number(value()); break;
      case '--replay': replay = value(); break;
      case '--verbose': case '-v': verbose = true; break;
      case '--help': case '-h':
        console.log('Live property comparisons: --case AREA | --cases A,B --seed UINT32 --runs N --replay FAILURE.json --verbose');
        console.log('A fresh seed is chosen by default. Inputs are generated once for both runtimes. Only failures are saved.');
        return;
      case '--generate': case '--typescript-only':
        throw new Error('Stored-output modes have been retired. Run live comparisons with SageMath installed.');
      default: throw new Error(`Unknown argument: ${argv[i]}`);
    }
  }
  // Validate options even for an area containing only explicit regressions.
  materializeSuite({ module: 'validation', cases: [] }, seed, runs);
  if (!await checkSageMath()) throw new Error('SageMath is required for live comparative tests. No cached-output fallback is used.');
  let testCases = replay
    ? [{ name: 'replay', content: JSON.stringify(JSON.parse(readText(replay)).suite) }]
    : await loadTestCases(caseFilters);
  if (!testCases.length) throw new Error('No property areas selected');
  if (caseFilters) for (const area of caseFilters)
    if (!testCases.some(c => c.name === area)) throw new Error(`Unknown property area: ${area}`);
  console.log(`Property seed: ${seed}; runs per generator: ${runs}; generator: ${GENERATOR_VERSION}`);
  let passed = 0, failed = 0, errors = 0;
  for (const source of testCases) {
    const suite = replay ? JSON.parse(source.content) : materializeSuite(JSON.parse(source.content), seed, runs);
    const input = JSON.stringify(suite);
    console.log(`Running live comparisons for ${suite.module}...`);
    try {
      const python = await runPythonTests(input);
      const typescript = await runTypeScriptTests(input);
      const summary = compareResults(python, typescript);
      printSummary(suite.module, summary, verbose);
      passed += summary.passed; failed += summary.failed; errors += summary.errors;
      if (summary.failed || summary.errors) {
        const keys = new Set(summary.results.filter(r => !r.match).map(r => JSON.stringify([r.function, r.args, r.seed])));
        const cases = suite.cases.map((c: any) => ({ function: c.function, rows: c.rows.filter((r: any[]) => {
          const args = r.slice(1).map(a => Array.isArray(a) ? `[${a.join(', ')}]` : String(a));
          return keys.has(JSON.stringify([c.function, args, r[0]]));
        }) })).filter((c: any) => c.rows.length);
        const directory = await mkdtemp(join(tmpdir(), 'sagemath-property-failure-'));
        const path = join(directory, `${suite.module}.json`);
        const revision = await runCommand('git', ['rev-parse', 'HEAD']);
        const sage = await runCommand('sage', ['--version']);
        await writeFile(path, JSON.stringify({ seed, runs, generator: GENERATOR_VERSION,
          sourceHash: createHash('sha256').update(source.content).digest('hex'),
          revision: revision.stdout.trim(), reference: sage.stdout.trim(),
          suite: { module: suite.module, cases },
        }, null, 2) + '\n');
        console.error(`Replay failing inputs: bun tests/property/compare.ts --replay ${path}`);
      }
    } catch (e) {
      errors++;
      const directory = await mkdtemp(join(tmpdir(), 'sagemath-property-failure-'));
      const path = join(directory, `${suite.module}.json`);
      await writeFile(path, JSON.stringify({ seed, runs, generator: GENERATOR_VERSION,
        sourceHash: createHash('sha256').update(source.content).digest('hex'), suite,
        error: e instanceof Error ? e.message : String(e),
      }) + '\n');
      console.error(`${suite.module}: ${e instanceof Error ? e.message : e}`);
      console.error(`Replay failed process inputs: bun tests/property/compare.ts --replay ${path}`);
    }
  }
  console.log(`Total: ${passed + failed + errors}, Passed: ${passed}, Failed: ${failed}, Errors: ${errors}`);
  if (failed || errors) process.exitCode = 1;
}

if (import.meta.main) {
  main().catch(e => { console.error(e); process.exitCode = 1; });
}
