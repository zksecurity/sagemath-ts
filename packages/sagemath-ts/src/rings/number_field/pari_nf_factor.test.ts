import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
const fixtures = await loadLiveNative(import.meta.url, "./pari_nf_factor.native.json");
import { zpFactorSquarefree, zpIsIrreducibleOverQ } from './pari_nf.js';
import { getrand, setrand } from '@sagemath-ts/parigp-ts/src/random.js';

for (const row of fixtures) {
  test(`bundled PARI legacy factor ${row.seed}`, () => {
    const op = BigInt(row.args[0]!),
      raw = row.args[1]!.slice(1, -1).trim(),
      coefficients = raw ? raw.split(',').map((c) => BigInt(c.trim())) : [],
      before = [...coefficients],
      saved = getrand();
    try {
      setrand(1n);
      const value =
        op === 1n ? zpIsIrreducibleOverQ(coefficients) : zpFactorSquarefree(coefficients);
      const result = JSON.stringify([value, getrand()], (_, v) =>
        typeof v === 'bigint' ? String(v) : v
      );
      expect({ result, error: null, errorType: null }).toEqual({
        result: row.result,
        error: row.error,
        errorType: row.errorType,
      });
      expect(coefficients).toEqual(before);
    } finally {
      setrand(saved);
    }
  }, 30000);
}
