import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
const fixtures = await loadLiveNative(import.meta.url, "./number_field_automorphism_order.native.json");
import { NumberField, RationalPolynomial } from './number_field.js';
import { Rational } from '../rational.js';

for (const row of fixtures) {
  test(`Sage automorphism order ${row.seed}`, () => {
    const coefficients = row.args[0]!.slice(1, -1)
      .split(',')
      .map((c) => BigInt(c.trim()));
    const denominator = BigInt(row.args[1]!);
    const K = new NumberField(
      new RationalPolynomial(coefficients.map((c) => new Rational(c, denominator))),
      'a'
    );
    const automorphisms = K.automorphisms();
    const result = JSON.stringify(
      automorphisms.map((sigma) => sigma.__call__(K.gen()).list().map(String))
    );
    expect({ result, error: null, errorType: null }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
    expect(automorphisms[0]!.is_identity()).toBe(true);
  }, 30000);
}
