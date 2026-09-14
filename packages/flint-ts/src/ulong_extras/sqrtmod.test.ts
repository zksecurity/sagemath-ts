import { expect, test } from 'bun:test';
import { _n_jacobi_unsigned, n_jacobi_unsigned, n_jacobi } from './jacobi.js';
import { n_is_square } from './is_square.js';
import { n_preinvert_limb } from './preinvert_limb.js';
import { n_powmod2_ui_preinv } from './powmod2_ui_preinv.js';
import { n_sqrtmod } from './sqrtmod.js';

test('FLINT word roots retain native representative choices across algorithm branches', () => {
  expect(n_sqrtmod(100n,101n)).toBe(10n);
  for(const p of [607n,613n,617n,65537n,18446744073709551557n]) {
    for(const x of [0n,1n,2n,5n,17n,p-1n]) {
      const a=x*x%p,r=n_sqrtmod(a,p);
      expect(r*r%p).toBe(a);
    }
  }
  expect(n_sqrtmod(3n,7n)).toBe(0n);
  expect(n_sqrtmod(2n,625n)).toBe(0n);
});

test('word Jacobi parity, unsigned lifts and exact square boundaries', () => {
  expect(n_jacobi_unsigned(2n,7n)).toBe(1);
  expect(_n_jacobi_unsigned(2n,7n,2)).toBe(-1);
  expect(n_jacobi(-1n,7n)).toBe(-1);
  const square=4294967295n**2n;
  expect(n_is_square(square)).toBe(true);
  expect(n_is_square(square-1n)).toBe(false);
  expect(n_is_square(square+1n)).toBe(false);
  const p=18446744073709551557n;
  expect(n_powmod2_ui_preinv(p-1n,1n<<63n,p,n_preinvert_limb(p))).toBe(1n);
});

test('word kernels reject inputs outside native argument preconditions', () => {
  expect(() => n_sqrtmod(1n,1n)).toThrow(RangeError);
  expect(() => n_jacobi_unsigned(1n,2n)).toThrow(RangeError);
  expect(() => n_jacobi(1n<<63n,3n)).toThrow(RangeError);
  expect(() => n_is_square(-1n)).toThrow(RangeError);
  expect(() => n_preinvert_limb(0n)).toThrow(RangeError);
  expect(() => n_powmod2_ui_preinv(1n,-1n,3n,0n)).toThrow(RangeError);
});
