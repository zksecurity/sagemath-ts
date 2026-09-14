import { expect, test } from 'bun:test';
import { F2Ms_ker, F2Ms_colelim } from './F2v.js';
import { setrand, getrand, pari_rand } from './random.js';
import { mpqsInternals } from './mpqs.js';

function seeded(seed: bigint, run: () => void): void {
  const saved = getrand();
  try {
    setrand(seed);
    run();
  } finally {
    setrand(saved);
  }
}

test('native sparse-kernel dispatch switches after row 640', () =>
  seeded(1n, () => {
    const state = getrand();
    expect(F2Ms_ker([[], [], []], 640)).toEqual([1n, 2n, 4n]);
    expect(getrand()).toBe(state);
    expect(F2Ms_ker([[], [], []], 641)).toEqual([4n, 6n, 5n]);
    expect(getrand()).not.toBe(state);
  }));

test('native singleton elimination returns surviving original column indices', () => {
  expect(F2Ms_colelim([[1], [2], [2], [], [3, 4], [4]], 641)).toEqual([2, 3, 4]);
});

test('native block Lanczos can return a proper kernel subspace', () =>
  seeded(1n, () => {
    expect(
      F2Ms_ker(
        Array.from({ length: 130 }, () => []),
        641
      )
    ).toHaveLength(64);
  }));

test('native singular-block retries consume thirteen complete 65-word starts', () =>
  seeded(3n, () => {
    const cycle = Array.from({ length: 65 }, (_, i) => [i + 1, ((i + 1) % 65) + 1]);
    // Bundled PARI: twelve failed starts, then success, with this seed and matrix.
    for (let i = 0; i < 13 * 65; i++) pari_rand();
    const expectedState = getrand();
    setrand(3n);
    expect(F2Ms_ker(cycle, 641)).toEqual([(1n << 65n) - 1n]);
    expect(getrand()).toBe(expectedState);
  }));

test('MPQS preserves one-based bit positions across word boundaries', () =>
  seeded(3n, () => {
    const cycle = Array.from({ length: 65 }, (_, i) => [i + 1, ((i + 1) % 65) + 1]);
    expect(mpqsInternals.F2Ms_ker(cycle, 641).map((v) => [...v])).toEqual([
      [0xfffffffe, 0xffffffff, 3],
    ]);
  }));
