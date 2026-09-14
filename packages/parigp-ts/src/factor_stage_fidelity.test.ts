import { expect, test } from 'bun:test';
import { PariError } from './errors.js';
import { pollardbrent, squfof, Z_pollardbrent } from './ifactor.js';

test('Pollard-Brent retains round budgets beyond the 32-bit shift boundary', () => {
  expect(Z_pollardbrent(1022117n, 2 ** 26, 0)).toEqual([1009n, 1013n]);
  expect(Z_pollardbrent(1022117n, 2 ** 40, 0)).toEqual([1009n, 1013n]);
});

test('Pollard-Brent seed selection adds native integers exactly', () => {
  expect(Z_pollardbrent(100160063n, 14, 2 ** 53)).toEqual([10007n, 10009n]);
  expect(Z_pollardbrent(1022117n, 14, 2 ** 60)).toEqual([1009n, 1013n]);
});

test('Pollard-Brent retry exhaustion preserves the native error', () => {
  let caught: unknown;
  try {
    Z_pollardbrent(299n, 14, 3);
  } catch (e) {
    caught = e;
  }
  expect(caught?.constructor).toBe(PariError);
  expect((caught as Error).message).toBe('bug in , please report.');
});

test('Pollard-Brent uses the native low-limb gate', () => {
  expect(pollardbrent(7876759719473978540033n)).toBeNull();
  expect(pollardbrent(29809938423118930378753n)).toEqual([11n, 2709994402101720943523n]);
  expect(pollardbrent(1022117n)).toBeNull();
});

test('SQUFOF declines at the native 64-bit build threshold', () => {
  expect(squfof(101000606000909n)).toBeNull();
});

test('the SQUFOF factor-list adapter includes the residual cofactor', () => {
  expect(squfof(122089n)).toEqual([11n, 11n, 1009n]);
});
