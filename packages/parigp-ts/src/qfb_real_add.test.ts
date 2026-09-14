import { expect, test } from 'bun:test';
import { addir, addrr, itor, real_0_bit, subrr } from './qfb.js';

test('native real addition preserves word precision through a carry', () => {
  expect(addrr(itor(1n, 64), itor(1n, 64))).toEqual({
    s: 1, e: 1, m: 1n << 63n, p: 64,
  });
});

test('native cancellation retains the remaining full mantissa word', () => {
  const x = { s: 1 as const, e: 0, p: 128, m: (1n << 127n) + 1n };
  expect(subrr(x, itor(1n, 128))).toEqual({
    s: 1, e: -127, m: 1n << 63n, p: 64,
  });
  expect(subrr(x, x)).toEqual(real_0_bit(-127));
});

test('native zero accuracy selects the integer conversion precision', () => {
  expect(addir(1n, real_0_bit(-65))).toEqual(itor(1n, 128));
  expect(addir(1n, real_0_bit(0))).toEqual(real_0_bit(0));
});

test('native whole-word extension uses a deterministic zero guard', () => {
  // Direct C addrr_sign with zero padding gives this frame. A padding word >=2
  // changes the native last bit: see the documented upstream out-of-range read.
  const x = { s: 1 as const, e: 0, p: 64, m: (1n << 63n) + ((1n << 63n) - 1n) / 3n };
  const y = { s: 1 as const, e: 64, p: 192, m: (1n << 192n) - 2n };
  expect(addrr(x, y)).toEqual({
    s: 1, e: 65, m: 170141183460469231737836218407120622932n, p: 128,
  });
});
