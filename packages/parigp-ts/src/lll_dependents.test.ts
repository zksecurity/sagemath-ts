import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
import { lllgramint } from './qfrep.js';
import { ZM_hnflll, ZM_snf_group } from './buch.js';
import { algdep } from './bibli1.js';
const fixtures = (await loadLiveNative(import.meta.url, "./lll_dependents.native.json.gz")).map(row => ({
  ...row, op: Number(row.args[0]), n: Number(row.args[1]), flag: Number(row.args[2]),
  flat: row.args[3]!.slice(1, -1).split(',').map(v => v.trim()).filter(Boolean),
}));
const relations = await loadLiveNative(import.meta.url, "./algdep_lll.native.json");

for (const [i, f] of fixtures.entries())
  test(`bundled PARI Gram, Hermite and Smith reduction ${i}`, () => {
    const flat = f.flat.map(BigInt),
      rows = f.n ? (flat.length - (f.op === 3 ? 2 : 0)) / f.n : 0;
    const M = Array.from({ length: f.n }, (_, j) =>
      Array.from({ length: rows }, (_, r) => flat[r * f.n + j]!)
    );
    const B = [[], ...M.map((c) => [0n, ...c])],
      saved = structuredClone([M, B]);
    const strip = (A: bigint[][] | null) => A?.slice(1).map((c) => c.slice(1)) ?? null;
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      let value: unknown;
      if (f.op === 0) value = lllgramint(M);
      else if (f.op === 3)
        value = lllgramint(M, { n: flat[flat.length - 2]!, d: flat[flat.length - 1]! });
      else if (f.op === 4) value = strip(ZM_hnflll(B, false, f.flag !== 0).H);
      else if (f.op === 1) {
        const r = ZM_hnflll(B, true, f.flag !== 0);
        value = [strip(r.H), strip(r.B)];
      } else {
        const r = ZM_snf_group(B);
        value = [r.D.slice(1), strip(r.Ui)];
      }
      result = JSON.stringify(value, (_, v) => (typeof v === 'bigint' ? String(v) : v));
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: f.result,
      error: f.error,
      errorType: f.errorType,
    });
    expect([M, B]).toEqual(saved);
  });

for (const [i, f] of relations.entries())
  test(`bundled PARI relation basis and double input ${i}`, () => {
    const bits = new DataView(new ArrayBuffer(8));
    bits.setBigUint64(0, BigInt(f.args[0]!));
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = '[' + algdep(bits.getFloat64(0), BigInt(f.args[1]!)).join(', ') + ']';
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: f.result,
      error: f.error,
      errorType: f.errorType,
    });
  });
