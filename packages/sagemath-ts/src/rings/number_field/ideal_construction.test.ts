import { expect, test } from 'bun:test';
import { idealintersect } from '@sagemath-ts/parigp-ts/src/base4.js';
import { randomi, setrand } from '@sagemath-ts/parigp-ts/src/random.js';
import { Integer } from '../integer_ring.js';
import { Rational } from '../rational.js';
import { NumberField, RationalPolynomial } from './number_field.js';
import { NumberFieldFractionalIdeal, NumberFieldIdeal } from './number_field_ideal.js';
import fixtures from './ideal_intersection.native.json' with { type: 'json' };

const field = (constant = 2n, name = 'a') =>
  new NumberField(RationalPolynomial.fromBigInts([constant, 0n, 1n]), name);

// Exact Sage observations are also exercised by nf_ideal_inputs/factory/cross_ideal.
test('ideal factories preserve fractional class, identity, and original generators', () => {
  const K = field(),
    I = K.ideal([new Integer(2n), new Rational(3n, 2n)]);
  expect(I).toBeInstanceOf(NumberFieldFractionalIdeal);
  expect(I.gens().map(String)).toEqual(['2', '3/2']);
  expect(K.ideal(I)).toBe(I);
  expect(K.fractional_ideal(I)).toBe(I);
  expect(I.add(K.ideal(3n))).toBeInstanceOf(NumberFieldFractionalIdeal);
  expect(
    I.intersection([new Rational(3n, 2n), 3n])
      .norm()
      .toString()
  ).toBe('9/4');
  expect(I.intersection([]).is_zero()).toBe(true);
});

test('zero factories and direct fractional constructors retain different errors', () => {
  const K = field(),
    zero = K.ideal([]);
  expect(zero).not.toBeInstanceOf(NumberFieldFractionalIdeal);
  expect(zero.gens().map(String)).toEqual(['0']);
  expect(K.ideal(0n)).not.toBe(zero);
  expect(zero.toString()).toBe('Ideal (0) of Number Field in a with defining polynomial x^2 + 2');
  const nonzero = 'gens must have a nonzero element (zero ideal is not a fractional ideal)';
  const empty = 'gens must have length at least 1 (zero ideal is not a fractional ideal)';
  expect(() => K.fractional_ideal(0n)).toThrow(nonzero);
  expect(() => new NumberFieldFractionalIdeal(K, [])).toThrow(empty);
  expect(() => new NumberFieldFractionalIdeal(K, [[]])).toThrow(nonzero);
  expect(() => new NumberFieldIdeal(K, [[]])).toThrow(empty);
});

test('base ideals are not accepted by the fractional-ideal identity shortcut', () => {
  const K = field(),
    base = new NumberFieldIdeal(K, [3n]);
  const message =
    'unable to convert Ideal (3) of Number Field in a with defining polynomial x^2 + 2 to Number Field in a with defining polynomial x^2 + 2';
  expect(() => K.ideal(base)).toThrow(message);
  expect(() => K.ideal(2n).intersection(base)).toThrow(message);
});

test('abstract foreign fields coerce rational constants but require nonrational embeddings', () => {
  const K = field(),
    L = field(3n, 'b');
  expect(K.__call__(L.__call__(new Rational(3n, 2n))).toString()).toBe('3/2');
  expect(K.ideal(L.ideal(2n)).norm().toString()).toBe('4');
  expect(() => K.ideal(L.gen())).toThrow(
    'No compatible natural embeddings found for Number Field in a with defining polynomial x^2 + 2 and Number Field in b with defining polynomial x^2 + 3'
  );
});

test('native ideal intersection normalizes redundant denominators without mutating inputs', () => {
  const nf = field(-2n)._pari_ideal_data();
  for (const fixture of fixtures) {
    const [n, seed, iDen, jDen, flatI, flatJ] = fixture.args as [
      number,
      number,
      string,
      string,
      string[],
      string[],
    ];
    const columns = (flat: string[]) =>
      Array.from({ length: flat.length / n }, (_, j) =>
        Array.from({ length: n }, (_, i) => BigInt(flat[i * (flat.length / n) + j]!))
      );
    const I = columns(flatI),
      J = columns(flatJ),
      original = structuredClone([I, J]);
    setrand(BigInt(seed));
    const [H, denominator] = idealintersect(nf, I, BigInt(iDen), J, BigInt(jDen));
    expect(
      JSON.stringify([H, denominator, randomi(1n << 128n)], (_, v) =>
        typeof v === 'bigint' ? String(v) : v
      )
    ).toBe(fixture.expected);
    expect([I, J]).toEqual(original);
  }
  expect(idealintersect(nf, [], 7n, [], 3n)).toEqual([[], 1n]);
  expect(() => idealintersect(nf, [], 0n, [], 1n)).toThrow(RangeError);
  expect(() => idealintersect(nf, [[1n]], 1n, [], 1n)).toThrow(RangeError);
});
