import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
const fixtures = await loadLiveNative(import.meta.url, "./finite_polynomial_roots.native.json");
import { GFpn, PrimeField, FiniteFieldElement } from '../finite_rings/finite_field_extension.js';
import { FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import { GF2 } from '../finite_rings/gf2.js';
import { PolynomialRing } from './polynomial_ring.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';

for (const row of fixtures) {
  test(`Sage ordered finite polynomial ${row.function} ${row.seed}`, () => {
    const kind = BigInt(row.args[0]!),
      p = BigInt(row.args[1]!),
      degree = Number(row.args[2]!);
    const parse = (s: string) =>
      s.slice(1, -1).trim()
        ? s
            .slice(1, -1)
            .split(',')
            .map((c) => BigInt(c.trim()))
        : [];
    const modulus = parse(row.args[3]!),
      coefficients = parse(row.args[4]!);
    let result: string;
    try {
      // GFpn's array adapter omits the implicit leading one.
      const B =
        kind === 0n
          ? new PrimeField(p)
          : kind === 1n
            ? new FiniteFieldPrime(p)
            : kind === 2n
              ? GF2
              : GFpn(p, degree, modulus.slice(0, -1).map(Number), 'a');
      const R = new PolynomialRing(B as unknown as CoefficientRing<RingElement>, 'x');
      const f = R.__call__(
        coefficients.map((c) =>
          kind === 3n ? (B as ReturnType<typeof GFpn>).fromInteger(c) : B.__call__(c)
        )
      );
      const before = f.coeffs.map(String);
      const value =
        row.function === 'finite_polynomial_factors'
          ? f
              .factor()
              .map(([g, e]) => [
                g.coeffs.map((c) =>
                  c instanceof FiniteFieldElement ? String(c.integer_representation()) : String(c)
                ),
                String(e),
              ])
          : f
              .roots()
              .map(([r, m]) => [
                r instanceof FiniteFieldElement ? String(r.integer_representation()) : String(r),
                String(m),
                (r as unknown as { parent: unknown }).parent === B,
              ]);
      expect(f.coeffs.map(String)).toEqual(before);
      result = JSON.stringify({ value });
    } catch (e) {
      result = JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
    }
    expect({ result, error: null, errorType: null }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
  });
}
