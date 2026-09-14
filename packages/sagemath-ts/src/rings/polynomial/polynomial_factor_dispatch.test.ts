import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
import { execFileSync } from 'node:child_process';
const fixtures = await loadLiveNative(import.meta.url, "./polynomial_factor_dispatch.native.json");
import { Integer } from '../integer_ring.js';
import { Rational } from '../rational.js';
import { QQ } from '../rational_field.js';
import { PrimeField } from '../finite_rings/finite_field_extension.js';
import { GF2 } from '../finite_rings/gf2.js';
import { PolynomialRing } from './polynomial_ring.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';

for (const row of fixtures) {
  test(`Sage polynomial factor dispatch ${row.seed}`, () => {
    const op = BigInt(row.args[0]!),
      kind = BigInt(row.args[1]!);
    const text = row.args[2]!.slice(1, -1).trim();
    const coefficients = text ? text.split(',').map((c) => BigInt(c.trim())) : [];
    const denominator = BigInt(row.args[3]!);
    const original = [...coefficients];
    const run = (): string => {
      try {
        const integerBase = {
          zero: () => new Integer(0n),
          one: () => new Integer(1n),
          __call__: (x: unknown) => new Integer(x as bigint),
          is_field: () => false,
          toString: () => 'Integer Ring',
        };
        const B =
          kind === 0n
            ? integerBase
            : kind === 1n
              ? QQ
              : kind === 2n
                ? new PrimeField(7n)
                : kind === 3n
                  ? GF2
                  : new PrimeField((1n << 89n) - 1n);
        const R = new PolynomialRing(B as unknown as CoefficientRing<RingElement>, 'x');
        const f = R.__call__(
          coefficients.map((c) => (denominator === 1n ? c : new Rational(c, denominator)))
        );
        const value =
          op === 1n
            ? f.is_irreducible()
            : f.factor().map(([g, e]) => [g.coeffs.map(String), String(e)]);
        return JSON.stringify({ value });
      } catch (e) {
        return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
      }
    };
    let result: string;
    if (coefficients.length > 101) {
      // Bound regressions to the obsolete high-degree recombination driver.
      // Timeout is failure, never an expected factorization result.
      const program = `import {Integer} from ${JSON.stringify(new URL('../integer_ring.ts', import.meta.url).pathname)};
        import {QQ} from ${JSON.stringify(new URL('../rational_field.ts', import.meta.url).pathname)};
        import {Rational} from ${JSON.stringify(new URL('../rational.ts', import.meta.url).pathname)};
        import {PrimeField} from ${JSON.stringify(new URL('../finite_rings/finite_field_extension.ts', import.meta.url).pathname)};
        import {GF2} from ${JSON.stringify(new URL('../finite_rings/gf2.ts', import.meta.url).pathname)};
        import {PolynomialRing} from ${JSON.stringify(new URL('./polynomial_ring.ts', import.meta.url).pathname)};
        const op=${op}n, kind=${kind}n, coefficients=${JSON.stringify(coefficients.map(String))}.map(BigInt), denominator=${denominator}n;
        console.log((${run.toString()})());`;
      result = execFileSync(process.execPath, ['-e', program], {
        encoding: 'utf8',
        timeout: 10000,
      }).trim();
    } else result = run();
    expect({ result, error: null, errorType: null }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
    expect(coefficients).toEqual(original);
  }, 30000);
}
