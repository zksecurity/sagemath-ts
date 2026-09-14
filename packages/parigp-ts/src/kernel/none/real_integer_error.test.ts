import { expect, test } from 'bun:test';
import { gcvtoi, itor, real_0_bit, type MpReal } from '../../qfb.js';
import { gcvtoi as buchGcvtoi } from '../../buch.js';

test('exact cancellation retains finite native real accuracy', () => {
  const x = itor(1n, 64);
  expect(gcvtoi(x)).toEqual([1n, -63]);
  expect(buchGcvtoi(x)).toEqual({ z: 1n, e: -63 });
});

test('zero conversion error follows allocation and the integral-bit branch', () => {
  const x: MpReal = { s: 0, e: 192, p: 128, m: 0n };
  expect(gcvtoi(real_0_bit(0))).toEqual([0n, 1]);
  expect(buchGcvtoi(real_0_bit(0))).toEqual({ z: 0n, e: 1 });
  expect(gcvtoi(x)).toEqual([0n, 65]);
  expect(buchGcvtoi(x)).toEqual({ z: 0n, e: 65 });
});
