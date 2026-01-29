/**
 * Utilities for property-based testing with fast-check.
 *
 * Provides:
 * - Custom arbitraries for mathematical types (bigint ranges, primes, etc.)
 * - Automatic regression test saving when counterexamples are found
 * - SageMath comparison helpers
 */

import * as fc from 'fast-check';
import { writeFileSync, readFileSync, existsSync, mkdirSync } from 'node:fs';
import { join, dirname } from 'node:path';

// ============================================
// Regression Test Management
// ============================================

export interface RegressionEntry {
  id: string;
  date: string;
  module: string;
  property: string;
  counterexample: string; // JSON stringified
  notes?: string;
}

const REGRESSIONS_DIR = join(import.meta.dir, '..', '..', 'regressions');
const REGISTRY_FILE = join(REGRESSIONS_DIR, 'registry.json');

/**
 * Load existing regression entries
 */
export function loadRegressions(): RegressionEntry[] {
  if (!existsSync(REGISTRY_FILE)) {
    return [];
  }
  return JSON.parse(readFileSync(REGISTRY_FILE, 'utf-8'));
}

/**
 * Save a new regression when a property test fails.
 *
 * @example
 * fc.assert(fc.property(...), {
 *   reporter: createRegressionReporter('arith', 'gcd commutativity')
 * });
 */
export function saveRegression(
  module: string,
  property: string,
  counterexample: unknown[],
  notes?: string
): string {
  // Ensure directory exists
  if (!existsSync(REGRESSIONS_DIR)) {
    mkdirSync(REGRESSIONS_DIR, { recursive: true });
  }

  const regressions = loadRegressions();

  const entry: RegressionEntry = {
    id: `${module}-${Date.now()}-${Math.random().toString(36).slice(2, 8)}`,
    date: new Date().toISOString(),
    module,
    property,
    counterexample: JSON.stringify(counterexample, (_, v) =>
      typeof v === 'bigint' ? `${v}n` : v
    ),
    notes,
  };

  regressions.push(entry);
  writeFileSync(REGISTRY_FILE, JSON.stringify(regressions, null, 2));

  console.log(`\n📝 Regression saved: ${entry.id}`);
  console.log(`   Module: ${module}`);
  console.log(`   Property: ${property}`);
  console.log(`   Counterexample: ${entry.counterexample}`);

  return entry.id;
}

/**
 * Create a reporter that auto-saves regressions on failure.
 */
export function createRegressionReporter(module: string, property: string) {
  return (result: fc.RunDetails<unknown[]>) => {
    if (result.failed && result.counterexample) {
      saveRegression(module, property, result.counterexample);
    }
  };
}

/**
 * Run a property test with automatic regression saving.
 * When a property fails, the counterexample is saved before the error is thrown.
 */
export function assertProperty<T extends unknown[]>(
  module: string,
  propertyName: string,
  arbitrary: fc.Arbitrary<T>,
  predicate: (...args: T) => boolean | void,
  options?: fc.Parameters<T>
): void {
  try {
    fc.assert(
      fc.property(arbitrary, (...args: T) => predicate(...args)),
      options
    );
  } catch (error) {
    // Extract counterexample from error and save regression
    if (error instanceof Error && error.message.includes('Counterexample:')) {
      const match = error.message.match(/Counterexample: \[(.*?)\]/s);
      if (match) {
        saveRegression(module, propertyName, [match[1]], error.message);
      }
    }
    // Re-throw so the test fails
    throw error;
  }
}

// ============================================
// Custom Arbitraries for Mathematical Types
// ============================================

/**
 * Arbitrary for positive bigints in a range.
 */
export function positiveBigInt(
  max: bigint = 10n ** 15n
): fc.Arbitrary<bigint> {
  return fc.bigInt({ min: 1n, max });
}

/**
 * Arbitrary for non-negative bigints in a range.
 */
export function nonNegativeBigInt(
  max: bigint = 10n ** 15n
): fc.Arbitrary<bigint> {
  return fc.bigInt({ min: 0n, max });
}

/**
 * Arbitrary for bigints that could be negative.
 */
export function anyBigInt(
  absMax: bigint = 10n ** 15n
): fc.Arbitrary<bigint> {
  return fc.bigInt({ min: -absMax, max: absMax });
}

/**
 * Arbitrary for small primes (useful for finite field testing).
 */
export function smallPrime(): fc.Arbitrary<bigint> {
  const primes = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n, 41n, 43n, 47n];
  return fc.constantFrom(...primes);
}

/**
 * Arbitrary for medium primes (for more extensive testing).
 */
export function mediumPrime(): fc.Arbitrary<bigint> {
  const primes = [
    2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n, 41n, 43n, 47n,
    53n, 59n, 61n, 67n, 71n, 73n, 79n, 83n, 89n, 97n, 101n, 103n, 107n, 109n, 113n,
    127n, 131n, 137n, 139n, 149n, 151n, 157n, 163n, 167n, 173n, 179n, 181n, 191n, 193n, 197n, 199n,
  ];
  return fc.constantFrom(...primes);
}

/**
 * Arbitrary for cryptographic-size primes.
 * These are fixed known primes for reproducibility.
 */
export function largePrime(): fc.Arbitrary<bigint> {
  const primes = [
    // Small Mersenne primes
    2n ** 13n - 1n,  // 8191
    2n ** 17n - 1n,  // 131071
    2n ** 19n - 1n,  // 524287
    // Other known primes
    65537n,          // F4 Fermat prime
    104729n,         // 10000th prime
    1000003n,        // Just over 10^6
    15485863n,       // 1 millionth prime
  ];
  return fc.constantFrom(...primes);
}

/**
 * Arbitrary for coprime pairs.
 */
export function coprimePair(
  max: bigint = 10n ** 10n
): fc.Arbitrary<[bigint, bigint]> {
  return fc
    .tuple(
      fc.bigInt({ min: 1n, max }),
      fc.bigInt({ min: 1n, max })
    )
    .filter(([a, b]) => gcdSimple(a, b) === 1n);
}

// Simple GCD for filtering (doesn't use the library being tested)
function gcdSimple(a: bigint, b: bigint): bigint {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    [a, b] = [b, a % b];
  }
  return a;
}

/**
 * Arbitrary for polynomial coefficients (list of field elements).
 */
export function polynomialCoeffs(
  p: bigint,
  maxDegree: number = 10
): fc.Arbitrary<bigint[]> {
  return fc
    .array(fc.bigInt({ min: 0n, max: p - 1n }), { minLength: 1, maxLength: maxDegree + 1 })
    .filter((coeffs) => coeffs.length === 0 || coeffs[coeffs.length - 1] !== 0n);
}

/**
 * Arbitrary for a point on an elliptic curve (as [x, y] or 'infinity').
 */
export function curvePoint(
  p: bigint
): fc.Arbitrary<[bigint, bigint] | 'infinity'> {
  return fc.oneof(
    fc.constant('infinity' as const),
    fc.tuple(
      fc.bigInt({ min: 0n, max: p - 1n }),
      fc.bigInt({ min: 0n, max: p - 1n })
    )
  );
}

// ============================================
// SageMath Comparison Helpers
// ============================================

/**
 * Format a value for comparison with SageMath output.
 */
export function formatForSage(value: unknown): string {
  if (typeof value === 'bigint') {
    return value.toString();
  }
  if (typeof value === 'boolean') {
    return value ? 'True' : 'False';
  }
  if (Array.isArray(value)) {
    // Check if it's a factorization (array of [prime, exp] pairs)
    if (value.length > 0 && Array.isArray(value[0]) && value[0].length === 2) {
      // Format as SageMath factorization: "2^2 * 3 * 5^2"
      return value
        .map(([p, e]: [bigint, bigint]) =>
          e === 1n ? p.toString() : `${p}^${e}`
        )
        .join(' * ') || '1';
    }
    // Regular list
    return `[${value.map(formatForSage).join(', ')}]`;
  }
  if (value === null) {
    return 'None';
  }
  return String(value);
}

/**
 * Parse SageMath output into TypeScript value.
 */
export function parseFromSage(sageOutput: string): unknown {
  const trimmed = sageOutput.trim();

  // Boolean
  if (trimmed === 'True') return true;
  if (trimmed === 'False') return false;

  // None/null
  if (trimmed === 'None') return null;

  // Integer
  if (/^-?\d+$/.test(trimmed)) {
    return BigInt(trimmed);
  }

  // List
  if (trimmed.startsWith('[') && trimmed.endsWith(']')) {
    // This is simplified; real parsing would need proper bracket matching
    const inner = trimmed.slice(1, -1);
    if (inner === '') return [];
    // Simple split for basic cases
    return inner.split(',').map((s) => parseFromSage(s.trim()));
  }

  // Tuple (convert to array)
  if (trimmed.startsWith('(') && trimmed.endsWith(')')) {
    const inner = trimmed.slice(1, -1);
    if (inner === '') return [];
    return inner.split(',').map((s) => parseFromSage(s.trim()));
  }

  // Factorization "2^2 * 3 * 5"
  if (/^\d+(\^\d+)?(\s*\*\s*\d+(\^\d+)?)*$/.test(trimmed)) {
    const factors = trimmed.split('*').map((f) => f.trim());
    return factors.map((f) => {
      if (f.includes('^')) {
        const [p, e] = f.split('^');
        return [BigInt(p), BigInt(e)];
      }
      return [BigInt(f), 1n];
    });
  }

  return trimmed;
}

// ============================================
// Test Configuration
// ============================================

/**
 * Default fast-check configuration for this project.
 */
export const defaultFcParams: fc.Parameters<unknown[]> = {
  numRuns: 100,           // Run 100 iterations per property
  seed: 42,               // Deterministic seed for reproducibility
  endOnFailure: false,    // Continue to find minimal counterexample
  verbose: false,         // Set true for debugging
};

/**
 * Thorough fast-check configuration (for CI or release testing).
 */
export const thoroughFcParams: fc.Parameters<unknown[]> = {
  numRuns: 1000,
  seed: 42,
  endOnFailure: false,
  verbose: false,
};
