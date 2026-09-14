import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
const fixtures = await loadLiveNative(import.meta.url, "./polynomial_distinct_roots.native.json");
import { Integer } from '../integer_ring.js';
import { Rational } from '../rational.js';
import { QQ } from '../rational_field.js';
import { GFpn, PrimeField, FiniteFieldElement } from '../finite_rings/finite_field_extension.js';
import { FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import { GF2 } from '../finite_rings/gf2.js';
import { PolynomialRing } from './polynomial_ring.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';

for (const row of fixtures) {
  test(`bundled Sage distinct polynomial roots ${row.seed}`, () => {
    const kind = BigInt(row.args[0]!);
    const parse = (s: string) =>
      s.slice(1, -1).trim()
        ? s
            .slice(1, -1)
            .split(',')
            .map((c) => BigInt(c.trim()))
        : [];
    let result: string;
    try {
      const finite = row.function === 'finite_polynomial_distinct_roots';
      const integerBase = {
        zero: () => new Integer(0n),
        one: () => new Integer(1n),
        __call__: (c: unknown) => new Integer(c as bigint),
        is_field: () => false,
        toString: () => 'Integer Ring',
      };
      let B: unknown, coefficients: unknown[];
      if (finite) {
        const p = BigInt(row.args[1]!),
          n = Number(row.args[2]!),
          modulus = parse(row.args[3]!);
        const base =
          kind === 0n
            ? new PrimeField(p)
            : kind === 1n
              ? new FiniteFieldPrime(p)
              : kind === 2n
                ? GF2
                : GFpn(p, n, modulus.slice(0, -1).map(Number), 'a');
        B = base;
        coefficients = parse(row.args[4]!).map((c) =>
          kind === 3n ? (base as ReturnType<typeof GFpn>).fromInteger(c) : base.__call__(c)
        );
      } else {
        B = kind === 0n ? integerBase : QQ;
        coefficients = parse(row.args[1]!).map((c) =>
          kind === 0n ? c : new Rational(c, BigInt(row.args[2]!))
        );
      }
      const R = new PolynomialRing(B as CoefficientRing<RingElement>, 'x'),
        f = R.__call__(coefficients);
      const before = f.coeffs.map(String);
      const roots = f.roots({ multiplicities: false });
      expect(f.coeffs.map(String)).toEqual(before);
      result = JSON.stringify({
        value: roots.map((r) =>
          finite
            ? [
                r instanceof FiniteFieldElement ? String(r.integer_representation()) : String(r),
                (r as unknown as { parent: unknown }).parent === B,
              ]
            : [
                String(r),
                r instanceof Integer
                  ? 'Integer Ring'
                  : r instanceof Rational
                    ? 'Rational Field'
                    : typeof r,
              ]
        ),
      });
    } catch (e) {
      result = JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
    }
    expect({ result, error: null, errorType: null }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
  }, 30000);
}
