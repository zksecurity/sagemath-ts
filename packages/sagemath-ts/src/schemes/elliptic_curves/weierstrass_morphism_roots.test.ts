import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
import { execFileSync } from 'node:child_process';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import { GFpn } from '../../rings/finite_rings/finite_field_extension.js';
import { QQ } from '../../rings/rational_field.js';
import { EllipticCurve } from './constructor.js';
import { _isomorphisms, baseWI, WeierstrassIsomorphism } from './weierstrass_morphism.js';
const native = await loadLiveNative(import.meta.url, "./weierstrass_morphism_roots.native.json");

// The generic curve interface does not structurally model Rational; these
// native comparisons exercise the existing runtime QQ/finite-field adapters.
type Any = any;
function evaluate(name: string, args: Any[]): string {
  try {
    if (name === 'wm_isomorphism_argument_errors') {
      const K = GF(5n),
        E = EllipticCurve(K, [0n, 1n]);
      const values: Any[] = [E, null, 7, [], {}];
      return JSON.stringify({
        value: [..._isomorphisms(values[Number(args[0])], values[Number(args[1])])].length,
      });
    }
    const [pRaw, dRaw, m, l, r, t] = args;
    const p = BigInt(pRaw),
      degree = Number(dRaw);
    const K: Any = degree > 1 ? GFpn(p, degree, m.slice(0, -1).map(Number), 'a') : p ? GF(p) : QQ;
    const decode = (v: Any) => (degree > 1 ? K.fromInteger(BigInt(v)) : K.__call__(BigInt(v)));
    const encode = (v: Any) => (degree > 1 ? String(v.integer_representation()) : String(v));
    const E = EllipticCurve(K, l.map(decode));
    const F = EllipticCurve(
      K,
      t.length
        ? new baseWI(...(t.map(decode) as [Any, Any, Any, Any])).call(E.a_invariants())
        : r.map(decode)
    );
    const before = [E.a_invariants().map(encode), F.a_invariants().map(encode)];
    const rows = [..._isomorphisms(E, F)];
    const value = rows.map((row) => [
      row.map(encode),
      row.every((v) => !p || v.parent === K),
      new baseWI(...row).call(E.a_invariants()).every((v, i) => v.eq(F.a_invariants()[i]!)),
    ]);
    let first;
    try {
      first = { value: new WeierstrassIsomorphism(E, null, F).tuple().map(encode) };
    } catch (e) {
      first = { error: (e as Error).name, message: (e as Error).message };
    }
    if (
      JSON.stringify(before) !==
      JSON.stringify([E.a_invariants().map(encode), F.a_invariants().map(encode)])
    ) {
      throw Error('isomorphism enumeration mutated an input curve');
    }
    return JSON.stringify({ value: [value, first] });
  } catch (e) {
    return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
  }
}

for (const row of native) {
  test(`native Weierstrass polynomial roots ${row.function} ${row.seed}`, () => {
    // Native transcript arrays contain decimal integer literals. Preserve every
    // digit before JSON parsing, including the 61/89-bit transformation inputs.
    const args = row.args.map((arg) =>
      arg.startsWith('[') ? JSON.parse(arg.replace(/-?\d+/g, '"$&"')) : arg
    );
    if (row.seed >= 28600000) {
      // Guard large-field callers against accidentally enumerating their fields.
      const script = `
        import { GF } from '${import.meta.dir}/../../rings/finite_rings/finite_field_constructor.ts';
        import { GFpn } from '${import.meta.dir}/../../rings/finite_rings/finite_field_extension.ts';
        import { QQ } from '${import.meta.dir}/../../rings/rational_field.ts';
        import { EllipticCurve } from '${import.meta.dir}/constructor.ts';
        import { _isomorphisms, baseWI, WeierstrassIsomorphism } from '${import.meta.dir}/weierstrass_morphism.ts';
        console.log((${evaluate.toString()})(${JSON.stringify(row.function)}, ${JSON.stringify(args)}));`;
      const output = execFileSync(process.execPath, ['--eval', script], {
        encoding: 'utf8',
        timeout: 10_000,
      });
      expect(output.trim()).toBe(row.result);
    } else {
      expect(evaluate(row.function, args)).toBe(row.result);
    }
  }, 15_000);
}
