import { expect, test } from 'bun:test';
import { legacyNativeFixtures } from '../../../tests/property/native-live.mjs';
const fixtures = await legacyNativeFixtures(import.meta.url, './lll.fixtures.json');
import { fplll_fast } from './lll.js';

// Generate constrained inputs and execute unmodified PARI lll.c on each run.
for (const fixture of fixtures) {
  test(`PARI fast LLL status, basis and transform: native probe ${fixture.probe}`, () => {
    const [, m, n, dn, dd, en, ed, keep, track, flat] = fixture.args as [
      number,
      number,
      number,
      number,
      number,
      number,
      number,
      number,
      number,
      string[],
    ];
    const B = Array.from({ length: n }, (_, j) =>
      Array.from({ length: m }, (_, i) => BigInt(flat[i * n + j]!))
    );
    const original = B.map((c) => [...c]);
    const result = fplll_fast(B, dn / dd, en / ed, !!keep, !!track);
    expect(
      JSON.parse(JSON.stringify(result, (_, v) => (typeof v === 'bigint' ? String(v) : v)))
    ).toEqual(fixture.expected);
    expect(B).toEqual(original);
    const [, reduced, U] = result;
    if (U) {
      const product = U.map((column) =>
        Array.from({ length: m }, (_, i) =>
          column.reduce((sum, u, j) => sum + original[j]![i]! * u, 0n)
        )
      );
      expect(product).toEqual(reduced);
    }
  });
}

test('fast LLL rejects matrices without the native nonempty rectangular shape', () => {
  for (const B of [[], [[]], [[1n], []]]) {
    expect(() => fplll_fast(B)).toThrow('fplll_fast requires a nonempty rectangular basis');
  }
});
