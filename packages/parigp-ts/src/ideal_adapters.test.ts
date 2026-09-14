import { expect, test } from 'bun:test';
import { diviiround } from './gen3.js';
import { ZM_hnfcenter } from './hnf_snf.js';
import { nativeFixtures } from '../../../tests/property/native-live.mjs';
const fixtures = (await nativeFixtures(import.meta.url, './ideal_adapters.native.json')).map(row => ({
  ...row, args: row.args.map(a => a.startsWith('[') ? a.slice(1, -1).split(',').map(v => v.trim()).filter(Boolean) : a),
  expected: row.error === null ? { result: row.result! } : { error: row.error },
}));

for (const [i, fixture] of fixtures.entries()) {
  test(`original PARI rounding and centered HNF ${i}`, () => {
    const operation = () => {
      if (fixture.function === 'pari_diviiround')
        return diviiround(BigInt(fixture.args[0] as string), BigInt(fixture.args[1] as string));
      const n = Number(fixture.args[0]),
        flat = fixture.args[1] as string[];
      const input = Array.from({ length: n }, (_, j) =>
        Array.from({ length: n }, (_, row) => BigInt(flat[row * n + j]!))
      );
      const copy = structuredClone(input),
        result = ZM_hnfcenter(input);
      expect(input).toEqual(copy);
      return result;
    };
    if ('error' in fixture.expected) {
      expect(operation).toThrow(fixture.expected.error!);
    } else {
      expect(JSON.stringify(operation(), (_, v) => (typeof v === 'bigint' ? String(v) : v))).toBe(
        fixture.expected.result
      );
    }
  });
}

test('centered HNF shape and diagonal adapter boundaries', () => {
  expect(() => ZM_hnfcenter([[0n]])).toThrow(RangeError);
  expect(() => ZM_hnfcenter([[1n, 0n]])).toThrow(RangeError);
});
