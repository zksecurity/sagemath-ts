import { expect, test } from 'bun:test';
import {
  _fmpz_poly_gcd,
  _fmpz_poly_gcd_heuristic,
  _fmpz_poly_gcd_modular,
  _fmpz_poly_gcd_subresultant,
} from './gcd.js';

// Exact expectations are also compared against the FLINT C entry points in
// tests/property/{python,typescript}/areas/polynomial_ops.*.
for (const kernel of [
  _fmpz_poly_gcd,
  _fmpz_poly_gcd_subresultant,
  _fmpz_poly_gcd_heuristic,
  _fmpz_poly_gcd_modular,
]) {
  test(`${kernel.name}: content, sign and immutable coefficient buffers`, () => {
    const a = Object.freeze([-6n, 0n, 6n, 0n]);
    const b = Object.freeze([9n, -9n]);
    expect(kernel(a, b)).toEqual([-3n, 3n]);
    expect(kernel([], [-6n, -12n])).toEqual([6n, 12n]);
    expect(kernel([], [])).toEqual([]);
    expect(a).toEqual([-6n, 0n, 6n, 0n]);
  });
}
