import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
const fixtures = await loadLiveNative(import.meta.url, "./integer_mod_root_helpers.native.json");
import { IntegerModRing, Zmod } from './integer_mod_ring.js';
import { IntegerMod } from './integer_mod.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import type { CoefficientRing, RingElement } from '../polynomial/polynomial_element.js';

for (const row of fixtures)
  test(`bundled Sage modular root helper ${row.seed}`, () => {
    const values = row.args.map((s) =>
      s.startsWith('[')
        ? s.slice(1, -1).trim()
          ? s
              .slice(1, -1)
              .split(',')
              .map((c) => BigInt(c.trim()))
          : []
        : BigInt(s)
    );
    let result: string;
    const normalize = (v: unknown): unknown =>
      Array.isArray(v) ? v.map(normalize) : typeof v === 'boolean' || v === null ? v : String(v);
    try {
      let value: unknown;
      const n = values[0] as bigint;
      if (row.function === 'mi_residue_root_lift') {
        const e = values[1] as bigint,
          K = Zmod(n ** e) as IntegerModRing,
          R = new PolynomialRing(K as unknown as CoefficientRing<IntegerMod & RingElement>, 'x');
        const f = R.__call__(values[2] as bigint[]),
          r = (Zmod(n) as IntegerModRing).__call__(values[3] as bigint);
        const before = f.coeffs.map(String);
        value = IntegerModRing._lift_residue_field_root(n, e, f, f.derivative(), r).map((v) => [
          String(v),
          String(v.parent),
          v.parent === K,
        ]);
        expect(f.coeffs.map(String)).toEqual(before);
      } else {
        const K = values[1] ? new IntegerModRing(n) : (Zmod(n) as IntegerModRing);
        if (row.function === 'mi_field') {
          const F = K.field();
          value = [
            String(F),
            F.order,
            F.characteristic,
            F.degree,
            F === K.field(),
            (values[2] as bigint[]).map((v) => [String(F.__call__(v)), F.__call__(v).parent === F]),
          ];
        } else {
          const f = K.factored_order();
          value = [f, f === K.factored_order()];
        }
      }
      result = JSON.stringify({ value: normalize(value) });
    } catch (e) {
      result = JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
    }
    expect({ result, error: null, errorType: null }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
  }, 30000);
