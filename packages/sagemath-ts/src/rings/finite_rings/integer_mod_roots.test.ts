import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";


import { execFileSync } from 'node:child_process';
import { expect, test } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./integer_mod_roots.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];
import { IntegerMod } from './integer_mod.js';
import { IntegerModRing, Zmod } from './integer_mod_ring.js';
import type { FiniteFieldElement } from './finite_field_prime.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import type { CoefficientRing, RingElement } from '../polynomial/polynomial_element.js';

for (const row of fixtures)
  test(`bundled Sage modular polynomial roots ${row.seed}`, () => {
    const K = Zmod(BigInt(row.args[0]!)) as IntegerModRing;
    const R = new PolynomialRing(K as unknown as CoefficientRing<IntegerMod & RingElement>, 'x');
    const text = row.args[1]!.slice(1, -1).trim(),
      coefficients = text ? text.split(',').map((c) => BigInt(c.trim())) : [];
    const f = R.__call__(coefficients),
      before = f.coeffs.map(String);
    const multiplicities = BigInt(row.args[2]!) !== 0n;
    let result: string;
    if (K.order > 1_000_000_000_000n && !multiplicities) {
      // Composite moduli exceed 10^12: enumeration cannot complete this watch.
      const program = `import {Zmod} from ${JSON.stringify(new URL('./integer_mod_ring.ts', import.meta.url).pathname)};
      import {PolynomialRing} from ${JSON.stringify(new URL('../polynomial/polynomial_ring.ts', import.meta.url).pathname)};
      const K=Zmod(${row.args[0]}n),f=new PolynomialRing(K,'x').__call__(${JSON.stringify(coefficients.map(String))}.map(BigInt));
      const before=f.coeffs.map(String);
      const value=f.roots({multiplicities:false}).map(r=>[String(r),String(r.parent),r.parent===K]);
      if(JSON.stringify(before)!==JSON.stringify(f.coeffs.map(String)))throw Error('input mutated');
      console.log(JSON.stringify({value}));`;
      result = execFileSync(process.execPath, ['-e', program], {
        encoding: 'utf8',
        timeout: 10000,
      }).trim();
    } else {
      try {
        const hook = row.function === 'mi_modular_roots_hook' ? BigInt(row.args[3]!) : -1n;
        const roots =
          hook < 0
            ? multiplicities
              ? f.roots()
              : f.roots({ multiplicities: false })
            : K._roots_univariate_polynomial(f, {
                multiplicities,
                ring: hook === 0n ? null : hook === 1n ? K : (Zmod(K.order + 1n) as IntegerModRing),
                algorithm: 'ignored',
              });
        const value = multiplicities
          ? (roots as Array<[FiniteFieldElement, number]>).map(([r, m]) => [
              String(r),
              String(m),
              String(r.parent),
              r.parent === K.field(),
            ])
          : (roots as IntegerMod[]).map((r) => [String(r), String(r.parent), r.parent === K]);
        result = JSON.stringify({ value });
      } catch (e) {
        result = JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
      }
    }
    expect({ result, error: null, errorType: null }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
    expect(f.coeffs.map(String)).toEqual(before);
  }, 30000);
