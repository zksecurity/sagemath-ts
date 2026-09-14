import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
const fixtures = await loadLiveNative(import.meta.url, "./finite_polynomial_roots_state.native.json");
import { getrand, setrand } from '@sagemath-ts/parigp-ts/src/random.js';
import { PrimeField } from '../finite_rings/finite_field_extension.js';
import { FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import { GF2 } from '../finite_rings/gf2.js';
import { PolynomialRing } from './polynomial_ring.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';

for (const row of fixtures) {
  test(`bundled PARI finite polynomial root effects ${row.seed}`, () => {
    const kind = BigInt(row.args[0]!),
      p = BigInt(row.args[1]!);
    const text = row.args[2]!.slice(1, -1).trim();
    const coefficients = text ? text.split(',').map((c) => BigInt(c.trim())) : [];
    const saved = getrand();
    try {
      setrand(1n);
      const B = kind === 0n ? new PrimeField(p) : kind === 1n ? new FiniteFieldPrime(p) : GF2;
      const R = new PolynomialRing(B as unknown as CoefficientRing<RingElement>, 'x');
      const f = R.__call__(coefficients),
        before = f.coeffs.map(String);
      const value =
        row.function === 'finite_polynomial_distinct_roots_state'
          ? f
              .roots({ multiplicities: false })
              .map((r) => [String(r), (r as unknown as { parent: unknown }).parent === B])
          : f
              .roots()
              .map(([r, m]) => [
                String(r),
                String(m),
                (r as unknown as { parent: unknown }).parent === B,
              ]);
      const result = JSON.stringify({ value: [value, String(getrand())] });
      expect({ result, error: null, errorType: null }).toEqual({
        result: row.result,
        error: row.error,
        errorType: row.errorType,
      });
      expect(f.coeffs.map(String)).toEqual(before);
    } finally {
      setrand(saved);
    }
  });
}
