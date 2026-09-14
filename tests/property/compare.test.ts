import { expect, test } from 'bun:test';
import { compareResults } from './compare.js';

const result = (errorType: string | null, error: string | null = 'bad input') => ({
  function: 'regression',
  args: ['0'],
  seed: 1,
  result: null,
  error,
  errorType,
});

test('different exceptions do not establish behavioral equivalence', () => {
  const comparison = compareResults([result('ValueError')], [result('RangeError')]);
  expect(comparison.passed).toBe(0);
  expect(comparison.failed).toBe(1);
});

test('empty exception messages remain errors, not successful null results', () => {
  const comparison = compareResults([result('AssertionError', '')], [result(null, null)]);
  expect(comparison.passed).toBe(0);
  expect(comparison.failed).toBe(1);
});

test('matching exceptions compare their class', () => {
  expect(compareResults([result('ValueError')], [result('ValueError')]).passed).toBe(1);
});

test('unknown dispatch and legacy exception transcripts cannot pass', () => {
  expect(
    compareResults(
      [result('ValueError', 'Unknown function: missing')],
      [result('ValueError', 'Unknown function: missing')]
    ).errors
  ).toBe(1);
  const legacy = { ...result(null), errorType: undefined };
  expect(compareResults([legacy], [legacy]).passed).toBe(0);
});

test('different original exception messages fail comparison', () => {
  const comparison = compareResults(
    [result('TypeError', 'cannot write image')],
    [result('TypeError', 'cannot create image')]
  );
  expect(comparison.passed).toBe(0);
  expect(comparison.failed).toBe(1);
  expect(comparison.results[0]!.match).toBe(false);
});

test('live comparisons reject missing, extra, and empty result streams', () => {
  const r = result(null, null);
  expect(compareResults([r, r], [r]).errors).toBe(1);
  expect(compareResults([r], [r, r]).errors).toBe(1);
  expect(compareResults([], []).errors).toBe(1);
  expect(compareResults([r, r], [r, r]).passed).toBe(2);
});

test('argument boundaries cannot collide in result lookup', () => {
  const a = { ...result(null, null), args: ['1,2', '3'] };
  const b = { ...result(null, null), args: ['1', '2,3'] };
  expect(compareResults([a], [b]).errors).toBe(2);
});
