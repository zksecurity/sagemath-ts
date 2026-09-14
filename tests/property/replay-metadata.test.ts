import { expect, test } from 'bun:test';
import * as fc from 'fast-check';
import { replayMetadata, defaultFcParams, thoroughFcParams, createRegressionReporter } from './typescript/utils/property-test-utils.js';

test('structured fast-check replay preserves a nested BigInt counterexample', () => {
  const property = fc.property(fc.array(fc.bigInt({ min: 9007199254740993n, max: 9007199254741993n }), { minLength: 2 }), () => false);
  const failed = fc.check(property, { seed: 123, numRuns: 10 });
  const metadata = replayMetadata(failed);
  const replayed = fc.check(property, { seed: metadata.seed, path: metadata.path, endOnFailure: true });
  expect(replayed.failed).toBe(true);
  expect(replayed.counterexample).toEqual(failed.counterexample);
  expect(metadata.generatorVersion).toBe(`fast-check@${fc.__version}`);
  expect(defaultFcParams.seed).toBeUndefined();
  expect(thoroughFcParams.seed).toBeUndefined();
});

test('custom regression reporter cannot turn an exhausted property into a pass', () => {
  // No counterexample is available here, so the reporter has nothing to persist.
  const exhausted = fc.check(fc.property(fc.integer(), () => { fc.pre(false); }), { numRuns: 1, maxSkipsPerRun: 1 });
  expect(exhausted.failed).toBe(true);
  expect(exhausted.counterexample).toBeNull();
  expect(() => createRegressionReporter('test', 'exhausted')(exhausted)).toThrow();
});
