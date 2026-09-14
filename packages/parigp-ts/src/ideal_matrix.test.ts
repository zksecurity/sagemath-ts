/** Compact regressions supplement the 1,494 direct bundled-native comparisons. */
import { expect, test } from 'bun:test';
import { integerMatrixInverse, rationalLift, wordMatrixAdjoint } from './_matrix_inverse.js';
import { hnf_divscale } from './hnf_snf.js';
import { halfgcdii } from './kernel/none/halfgcd.js';

test('word adjoint retains the singular rank-one result', () => {
  expect(
    wordMatrixAdjoint(
      [
        [1n, 2n],
        [2n, 4n],
      ],
      101n
    )
  ).toEqual([
    [4n, 99n],
    [99n, 1n],
  ]);
  expect(
    wordMatrixAdjoint(
      [
        [0n, 0n],
        [0n, 0n],
      ],
      101n
    )
  ).toEqual([
    [0n, 0n],
    [0n, 0n],
  ]);
});

test('integer inverse preserves primitive content and native singular sentinels', () => {
  expect(
    integerMatrixInverse([
      [2n, 4n],
      [0n, 6n],
    ])
  ).toEqual([
    [
      [3n, -2n],
      [0n, 1n],
    ],
    6n,
  ]);
  expect(
    integerMatrixInverse([
      [1n, 2n],
      [2n, 4n],
    ])
  ).toBeNull();
  const zero = Array.from({ length: 3 }, () => [0n, 0n, 0n]);
  expect(integerMatrixInverse(zero)).toEqual([zero, 0n]);
});

test('half-GCD matches the native equal-input swap, including zero', () => {
  expect(halfgcdii(0n, 0n)).toEqual([
    [
      [0n, 1n],
      [1n, 0n],
    ],
    [0n, 0n],
  ]);
  const [M, V] = halfgcdii(-(2n ** 8192n + 17n), 2n ** 8190n + 31n);
  const a = -(2n ** 8192n + 17n),
    b = 2n ** 8190n + 31n;
  expect(V).toEqual([M[0][0] * a + M[0][1] * b, M[1][0] * a + M[1][1] * b]);
  expect(M[0][0] * M[1][1] - M[0][1] * M[1][0]).toBeOneOf([-1n, 1n]);
});

test('rational lift and HNF division retain exact fractions', () => {
  expect(rationalLift(34n, 101n, 7n)).toEqual([1n, 3n]);
  expect(rationalLift(0n, 1n, 0n)).toEqual([0n, 1n]);
  expect(
    hnf_divscale(
      [
        [2n, 0n],
        [1n, 3n],
      ],
      [
        [4n, 0n],
        [5n, 3n],
      ],
      3n
    )
  ).toEqual([
    [6n, 0n],
    [6n, 3n],
  ]);
  expect(() => hnf_divscale([[2n]], [[1n]], 1n)).toThrow('integral solution');
  expect(() => hnf_divscale([[1n]], [], 1n)).toThrow('same dimension');
});
