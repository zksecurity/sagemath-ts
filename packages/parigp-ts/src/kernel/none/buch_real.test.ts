import { expect, test } from 'bun:test';
import * as buch from '../../buch.js';
import { PariError } from '../../errors.js';
import { itor, truncr, gcvtoi } from '../../qfb.js';

test('Buchmann constructors distinguish canonical and allocated native zeros', () => {
  expect(buch.real_0_bit(10)).toEqual({ s: 0, e: 10, m: 0n, p: 0 });
  expect(buch.real_0(64)).toEqual({ s: 0, e: -64, m: 0n, p: 0 });
  expect(buch.itor(0n, 200)).toEqual({ s: 0, e: -256, m: 0n, p: 256 });
  expect(buch.setprec(buch.real_0_bit(1), 65)).toEqual({ s: 0, e: -128, m: 0n, p: 128 });
  expect(buch.real_1(200)).toEqual({ s: 1, e: 0, m: 1n << 255n, p: 256 });
});

test('Buchmann arithmetic retains native precision and cancellation accuracy', () => {
  const a = itor(1n, 64), b = itor(1n, 128);
  expect(buch.addrr(a, b)).toEqual({ s: 1, e: 1, m: 1n << 63n, p: 64 });
  expect(buch.subrr(a, a)).toEqual({ s: 0, e: -63, m: 0n, p: 0 });
  expect(buch.mulrr(a, b)).toEqual(a);
  expect(buch.divrr(a, b)).toEqual(a);
  expect(buch.mulir(0n, { ...a, e: 3 })).toEqual({ s: 0, e: -61, m: 0n, p: 0 });
});

test('Buchmann square and multiplication preserve distinct native truncated paths', () => {
  const x = { s: 1, e: 0, p: 832, m: (5n << 829n) - 3n };
  expect(buch.sqrr(x)).toEqual({ s: 1, e: 0, p: 832, m: (25n << 827n) - 7n });
  expect(buch.mulrr(x, x)).toEqual(buch.sqrr(x));
  expect(buch.mulrr(x, { ...x })).toEqual({ s: 1, e: 0, p: 832, m: (25n << 827n) - 8n });
});

test('Buchmann word operations retain integers beyond binary64 precision', () => {
  const n = (1n << 53n) + 1n;
  expect(buch.mulur(n, itor(1n, 128))).toEqual({ s: 1, e: 53, p: 128, m: n << 74n });
  expect(buch.divru(itor(n, 128), n)).toEqual(itor(1n, 128));
});

test('native truncation checks precision while conversion with error metadata does not', () => {
  const x = { ...itor(1n, 64), e: 64 };
  for (const truncate of [truncr, buch.truncr]) {
    expect(() => truncate(x)).toThrow(PariError);
    expect(() => truncate(x)).toThrow('precision too low in truncr (precision loss in truncation)');
    expect(truncate({ ...x, e: 63 })).toBe(1n << 63n);
    expect(truncate({ ...x, s: 0 })).toBe(0n);
  }
  expect(gcvtoi(x)).toEqual([1n << 64n, 1]);
  expect(buch.gcvtoi(x)).toEqual({ z: 1n << 64n, e: 1 });
});
