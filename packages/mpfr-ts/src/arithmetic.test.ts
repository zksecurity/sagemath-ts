import { expect, test } from 'bun:test';
import { mpfr_init2, mpfr_set_str, mpfr_add, mpfr_mul, mpfr_get_d } from './index.js';

test('native MPFR addition preserves precision, aliases and subtraction midpoint status', () => {
  const x = mpfr_init2(127),
    y = mpfr_init2(128);
  mpfr_set_str(x, '1e-1000');
  mpfr_set_str(y, '-0.1');
  expect(mpfr_add(x, x, y)).toBe(2);
  expect(x).toEqual({
    precision: 127,
    kind: 'finite',
    sign: -1,
    mantissa: 136112946768375385385349842972707284582n,
    exponent: -3,
  });
});

test('native MPFR binary operations use destination precision and singular signs', () => {
  const x = mpfr_init2(100),
    y = mpfr_init2(100),
    z = mpfr_init2(2);
  mpfr_set_str(x, '1.5');
  mpfr_set_str(y, '1.5');
  expect(mpfr_mul(z, x, y)).toBe(-1);
  expect(mpfr_get_d(z)).toBe(2);
  expect(mpfr_add(z, x, y)).toBe(0);
  expect(mpfr_get_d(z)).toBe(3);
  mpfr_set_str(x, '-0');
  mpfr_set_str(y, '-0');
  expect(mpfr_add(z, x, y)).toBe(0);
  expect(Object.is(mpfr_get_d(z), -0)).toBe(true);
  mpfr_set_str(y, 'inf');
  expect(mpfr_mul(z, x, y)).toBe(0);
  expect(z.kind).toBe('nan');
});
