#!/usr/bin/env bun
/**
 * Property test comparison harness.
 *
 * Runs both Python/SageMath and TypeScript implementations with the same
 * test cases and compares the results to ensure behavioral equivalence.
 *
 * Usage:
 *   bun compare.ts                           # Run all comparison tests
 *   bun compare.ts --case arith              # Run specific module
 *   bun compare.ts --generate                # Generate expected results only
 *   bun compare.ts --typescript-only         # Run TypeScript tests only
 *   bun compare.ts --verbose                 # Verbose output
 */

import { spawn } from 'node:child_process';
import { existsSync } from 'node:fs';
import { mkdir, readFile, readdir, writeFile } from 'node:fs/promises';
import { basename, join } from 'node:path';

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
interface TestResult {
  function: string;
  args: string[];
  result: string | null;
  error: string | null;
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
const TRANSCRIPTS_DIR = join(SCRIPT_DIR, 'transcripts');
const PYTHON_RUNNER = join(SCRIPT_DIR, 'python/runner.py');
const TS_RUNNER = join(SCRIPT_DIR, 'typescript/runner.ts');

/**
 * Run a command and capture output
 */
async function runCommand(
  command: string,
  args: string[],
  input?: string
): Promise<{ stdout: string; stderr: string; exitCode: number }> {
  return new Promise((resolve) => {
    const proc = spawn(command, args, {
      cwd: PROJECT_ROOT,
      stdio: ['pipe', 'pipe', 'pipe'],
    });

    let stdout = '';
    let stderr = '';

    proc.stdout.on('data', (data) => {
      stdout += data.toString();
    });

    proc.stderr.on('data', (data) => {
      stderr += data.toString();
    });

    if (input) {
      proc.stdin.write(input);
      proc.stdin.end();
    }

    proc.on('close', (code) => {
      resolve({
        stdout,
        stderr,
        exitCode: code ?? 1,
      });
    });

    proc.on('error', (err) => {
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
async function runPythonTests(testCaseJson: string): Promise<TestResult[]> {
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
async function runTypeScriptTests(testCaseJson: string): Promise<TestResult[]> {
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
function compareResults(pythonResults: TestResult[], tsResults: TestResult[]): ComparisonSummary {
  const results: ComparisonResult[] = [];
  let passed = 0;
  let failed = 0;
  let errors = 0;

  // Create a map for easier lookup
  const tsMap = new Map<string, TestResult>();
  for (const r of tsResults) {
    const key = `${r.function}:${r.args.join(',')}:${r.seed}`;
    tsMap.set(key, r);
  }

  for (const pyResult of pythonResults) {
    const key = `${pyResult.function}:${pyResult.args.join(',')}:${pyResult.seed}`;
    const tsResult = tsMap.get(key);

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

    // Both have errors
    if (pyResult.error && tsResult.error) {
      // Consider it a match if both errored (behavior matches)
      results.push({
        function: pyResult.function,
        args: pyResult.args,
        seed: pyResult.seed,
        pythonResult: pyResult.result,
        typescriptResult: tsResult.result,
        pythonError: pyResult.error,
        typescriptError: tsResult.error,
        match: true, // Both errored
      });
      passed++;
      continue;
    }

    // One has error, one doesn't
    if (pyResult.error || tsResult.error) {
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

  return {
    total: pythonResults.length,
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
async function loadTestCases(caseFilter?: string): Promise<{ name: string; content: string }[]> {
  const cases: { name: string; content: string }[] = [];

  if (!existsSync(CASES_DIR)) {
    console.error(`${colors.red}Test cases directory not found: ${CASES_DIR}${colors.reset}`);
    return cases;
  }

  const files = await readdir(CASES_DIR);

  for (const file of files) {
    if (!file.endsWith('.cases.json')) continue;

    const name = basename(file, '.cases.json');

    if (caseFilter && name !== caseFilter) continue;

    const content = await readFile(join(CASES_DIR, file), 'utf-8');
    cases.push({ name, content });
  }

  return cases;
}

/**
 * Save transcripts for debugging
 */
async function saveTranscripts(
  caseName: string,
  pythonResults: TestResult[],
  tsResults: TestResult[]
): Promise<void> {
  const pyDir = join(TRANSCRIPTS_DIR, 'python');
  const tsDir = join(TRANSCRIPTS_DIR, 'typescript');

  await mkdir(pyDir, { recursive: true });
  await mkdir(tsDir, { recursive: true });

  await writeFile(join(pyDir, `${caseName}.json`), JSON.stringify(pythonResults, null, 2));

  await writeFile(join(tsDir, `${caseName}.json`), JSON.stringify(tsResults, null, 2));
}

/**
 * Main entry point
 */
async function main(): Promise<void> {
  // Parse arguments
  const args = process.argv.slice(2);
  let caseFilter: string | undefined;
  let generateOnly = false;
  let typescriptOnly = false;
  let verbose = false;

  for (let i = 0; i < args.length; i++) {
    switch (args[i]) {
      case '--case':
        caseFilter = args[++i];
        break;
      case '--generate':
        generateOnly = true;
        break;
      case '--typescript-only':
        typescriptOnly = true;
        break;
      case '--verbose':
      case '-v':
        verbose = true;
        break;
      case '--help':
      case '-h':
        console.log(`
Property Test Comparison Tool

Usage:
  bun compare.ts [options]

Options:
  --case <name>       Run only tests for specified module (e.g., arith)
  --generate          Generate expected results from SageMath only
  --typescript-only   Run TypeScript tests only (compare with saved transcripts)
  --verbose, -v       Show all test results, not just failures
  --help, -h          Show this help message

Examples:
  bun compare.ts                      # Run all comparison tests
  bun compare.ts --case arith         # Run arith module tests only
  bun compare.ts --typescript-only    # Run TS tests against saved transcripts
        `);
        return;
    }
  }

  console.log(`${colors.blue}========================================${colors.reset}`);
  console.log(`${colors.blue}  Property Test Comparison${colors.reset}`);
  console.log(`${colors.blue}  SageMath vs sagemath-ts${colors.reset}`);
  console.log(`${colors.blue}========================================${colors.reset}`);

  // Check for SageMath if needed
  if (!typescriptOnly) {
    const hasSage = await checkSageMath();
    if (!hasSage) {
      console.error(`\n${colors.yellow}Warning: SageMath not found.${colors.reset}`);
      console.error('Install SageMath to run comparison tests, or use --typescript-only.\n');
      console.error('Install instructions:');
      console.error('  macOS: brew install sage');
      console.error('  Ubuntu: sudo apt-get install sagemath\n');

      if (!generateOnly) {
        typescriptOnly = true;
        console.log('Falling back to TypeScript-only mode.\n');
      } else {
        process.exit(1);
      }
    }
  }

  // Load test cases
  const testCases = await loadTestCases(caseFilter);

  if (testCases.length === 0) {
    console.error(`\n${colors.red}No test cases found.${colors.reset}`);
    console.error(`Create test case files in: ${CASES_DIR}`);
    process.exit(1);
  }

  console.log(
    `\nFound ${testCases.length} test case file(s): ${testCases.map((c) => c.name).join(', ')}`
  );

  // Run tests
  let totalPassed = 0;
  let totalFailed = 0;
  let totalErrors = 0;

  for (const testCase of testCases) {
    let pythonResults: TestResult[] = [];
    let tsResults: TestResult[] = [];

    // Run Python/SageMath tests
    if (!typescriptOnly) {
      try {
        console.log(
          `\n${colors.cyan}Running SageMath tests for ${testCase.name}...${colors.reset}`
        );
        pythonResults = await runPythonTests(testCase.content);
        console.log(`  ${colors.green}Completed ${pythonResults.length} test(s)${colors.reset}`);
      } catch (e) {
        console.error(`  ${colors.red}Failed to run SageMath tests${colors.reset}`);
        totalErrors++;
        continue;
      }
    } else {
      // Load from saved transcripts
      const transcriptPath = join(TRANSCRIPTS_DIR, 'python', `${testCase.name}.json`);
      if (existsSync(transcriptPath)) {
        const content = await readFile(transcriptPath, 'utf-8');
        pythonResults = JSON.parse(content);
        console.log(
          `\n${colors.cyan}Loaded ${pythonResults.length} SageMath result(s) from transcript${colors.reset}`
        );
      } else {
        console.error(
          `\n${colors.yellow}No SageMath transcript found for ${testCase.name}${colors.reset}`
        );
        console.error('Run without --typescript-only first to generate transcripts.');
        continue;
      }
    }

    if (generateOnly) {
      // Save transcripts and continue
      await saveTranscripts(testCase.name, pythonResults, []);
      console.log(`  ${colors.green}Saved Python transcript${colors.reset}`);
      continue;
    }

    // Run TypeScript tests
    try {
      console.log(`${colors.cyan}Running TypeScript tests for ${testCase.name}...${colors.reset}`);
      tsResults = await runTypeScriptTests(testCase.content);
      console.log(`  ${colors.green}Completed ${tsResults.length} test(s)${colors.reset}`);
    } catch (e) {
      console.error(`  ${colors.red}Failed to run TypeScript tests${colors.reset}`);
      totalErrors++;
      continue;
    }

    // Save transcripts
    await saveTranscripts(testCase.name, pythonResults, tsResults);

    // Compare results
    const summary = compareResults(pythonResults, tsResults);
    printSummary(testCase.name, summary, verbose);

    totalPassed += summary.passed;
    totalFailed += summary.failed;
    totalErrors += summary.errors;
  }

  // Final summary
  console.log(`\n${colors.blue}========================================${colors.reset}`);
  console.log(`${colors.blue}  Final Summary${colors.reset}`);
  console.log(`${colors.blue}========================================${colors.reset}`);
  console.log(
    `Total: ${totalPassed + totalFailed + totalErrors}, ` +
      `${colors.green}Passed: ${totalPassed}${colors.reset}, ` +
      `${colors.red}Failed: ${totalFailed}${colors.reset}, ` +
      `${colors.yellow}Errors: ${totalErrors}${colors.reset}`
  );

  // Exit with error if any failures
  if (totalFailed > 0 || totalErrors > 0) {
    process.exit(1);
  }
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
