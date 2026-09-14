import { expect, test } from 'bun:test';
import { itor, shiftr } from './qfb.js';
import { mpfactr } from './trans2.js';
import { mpexp } from './trans1.js';
import { gen_product, gen_powu_i } from './bb_group.js';
import { quadratic_prec_mask } from './Zp.js';

test('PARI factorial returns its native real mantissa and a wide exponent', () => {
  expect(mpfactr(5n)).toEqual({ s: 1, e: 6n, m: 17293822569102704640n, p: 64 });
  expect(mpfactr(1n << 54n)).toEqual({ s: 1, e: 946788236117799971n, m: 17930630872012024344n, p: 64 });
  expect(() => mpfactr(1n << 56n)).toThrow('overflow in expo()');
  expect(() => mpfactr(-1n)).toThrow('domain error in factorial: argument < 0');
});

test('native Newton reciprocal retains the original last bit near exp(0)', () => {
  for (const [n, correction] of [[1n, 0n], [3n, 8n]]) {
    expect(mpexp(shiftr(itor(-n!, 4096), -4096))).toEqual({
      s: 1, e: -1n, p: 8192, m: (1n << 8192n) - (n! << 4096n) + correction!,
    });
  }
});

test('native product evaluates complete levels in the original callback order', () => {
  const trace: string[] = [];
  expect(gen_product([1n, 2n, 3n, 4n, 5n, 6n], (a, b) => {
    trace.push(`${a},${b}`);
    return a * b;
  })).toBe(720n);
  expect(trace).toEqual(['1,2', '4,5', '2,3', '20,6', '6,120']);
  expect(gen_product([], (a: bigint, b: bigint) => a * b)).toBe(1n);
});

test('native kernel preconditions reject values outside their C domains', () => {
  expect(() => mpfactr(1n, 65)).toThrow(RangeError);
  expect(() => mpfactr(1n << 63n)).toThrow(RangeError);
  expect(() => gen_powu_i(1n, 0n, (a) => a * a, (a, b) => a * b)).toThrow(RangeError);
  expect(() => quadratic_prec_mask(1)).toThrow(RangeError);
});
