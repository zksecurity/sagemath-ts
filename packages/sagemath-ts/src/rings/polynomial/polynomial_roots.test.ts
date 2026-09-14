import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
import { execFileSync } from 'node:child_process';
const fixtures = await loadLiveNative(import.meta.url, "./polynomial_roots.native.json");
import { Integer } from '../integer_ring.js';
import { Rational } from '../rational.js';
import { QQ } from '../rational_field.js';
import { PolynomialRing } from './polynomial_ring.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';

for (const row of fixtures) {
  test(`Sage polynomial roots ${row.seed}`, () => {
    const kind = BigInt(row.args[0]!);
    const text = row.args[1]!.slice(1, -1).trim();
    const coefficients = text ? text.split(',').map((c) => BigInt(c.trim())) : [];
    const denominator = BigInt(row.args[2]!);
    const before = [...coefficients];
    const run = (): string => {
      try {
        const integerBase = {
          zero: () => new Integer(0n),
          one: () => new Integer(1n),
          __call__: (x: unknown) => new Integer(x as bigint),
          is_field: () => false,
          toString: () => 'Integer Ring',
        };
        const R = new PolynomialRing(
          (kind === 0n ? integerBase : QQ) as unknown as CoefficientRing<RingElement>,
          'x'
        );
        const f = R.__call__(
          coefficients.map((c) => (kind === 0n ? c : new Rational(c, denominator)))
        );
        return JSON.stringify({
          value: f
            .roots()
            .map(([r, m]) => [
              String(r),
              String(m),
              r instanceof Integer
                ? 'Integer Ring'
                : r instanceof Rational
                  ? 'Rational Field'
                  : typeof r,
            ]),
        });
      } catch (e) {
        return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
      }
    };
    let result: string;
    if (coefficients.some((c) => c > 1n << 200n || c < -(1n << 200n)) && coefficients.length < 10) {
      // A regression to integer factoring can stall on these 234-bit semiprimes.
      // The watchdog is a test failure, never an accepted oracle result.
      const program = `import {Integer} from ${JSON.stringify(new URL('../integer_ring.ts', import.meta.url).pathname)};
        import {QQ} from ${JSON.stringify(new URL('../rational_field.ts', import.meta.url).pathname)};
        import {Rational} from ${JSON.stringify(new URL('../rational.ts', import.meta.url).pathname)};
        import {PolynomialRing} from ${JSON.stringify(new URL('./polynomial_ring.ts', import.meta.url).pathname)};
        const kind=${kind}n, coefficients=${JSON.stringify(coefficients.map(String))}.map(BigInt), denominator=${denominator}n;
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
    expect(coefficients).toEqual(before);
  }, 30000);
}
