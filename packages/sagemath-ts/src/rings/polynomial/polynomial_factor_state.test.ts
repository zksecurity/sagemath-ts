import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
import { execFileSync } from 'node:child_process';
const fixtures = await loadLiveNative(import.meta.url, "./polynomial_factor_state.native.json");
import { getrand, setrand } from '@sagemath-ts/parigp-ts/src/random.js';
import { Integer } from '../integer_ring.js';
import { Rational } from '../rational.js';
import { QQ } from '../rational_field.js';
import { PolynomialRing } from './polynomial_ring.js';
import { Polynomial, type CoefficientRing, type RingElement } from './polynomial_element.js';

for (const row of fixtures) {
  test(`bundled polynomial factor state ${row.seed}`, () => {
    const saved = getrand();
    try {
      setrand(1n);
      const op = BigInt(row.args[0]!),
        kind = BigInt(row.args[1]!),
        denominator = BigInt(row.args[3]!);
      const text = row.args[2]!.slice(1, -1).trim();
      const coefficients = text ? text.split(',').map((c) => BigInt(c.trim())) : [];
      const original = [...coefficients];
      const execute = (): string => {
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
          coefficients.map((c) => (denominator === 1n ? c : new Rational(c, denominator)))
        );
        const before = JSON.stringify(f.coeffs.map(String));
        let value: unknown;
        if (op === 0n)
          value = [f.factor().map(([g, e]) => [g.coeffs.map(String), String(e)]), getrand()];
        else
          value = [
            [f.is_irreducible(), getrand()],
            [f.is_irreducible(), getrand()],
            [new Polynomial([...f.coeffs], f.parent).is_irreducible(), getrand()],
          ];
        if (JSON.stringify(f.coeffs.map(String)) !== before)
          throw new Error('factorization mutated input coefficients');
        return JSON.stringify(value, (_, v) => (typeof v === 'bigint' ? String(v) : v));
      };
      let result: string;
      if (coefficients.length > 101) {
        // Regressions to the legacy recombination route must fail within a bound.
        const program = `import {getrand,setrand} from '@sagemath-ts/parigp-ts/src/random.js';
          import {Integer} from ${JSON.stringify(new URL('../integer_ring.ts', import.meta.url).pathname)};
          import {QQ} from ${JSON.stringify(new URL('../rational_field.ts', import.meta.url).pathname)};
          import {Rational} from ${JSON.stringify(new URL('../rational.ts', import.meta.url).pathname)};
          import {PolynomialRing} from ${JSON.stringify(new URL('./polynomial_ring.ts', import.meta.url).pathname)};
          import {Polynomial} from ${JSON.stringify(new URL('./polynomial_element.ts', import.meta.url).pathname)};
          const op=${op}n, kind=${kind}n, coefficients=${JSON.stringify(coefficients.map(String))}.map(BigInt), denominator=${denominator}n;
          setrand(1n);console.log((${execute.toString()})());`;
        result = execFileSync(process.execPath, ['-e', program], {
          encoding: 'utf8',
          timeout: 10000,
        }).trim();
      } else result = execute();
      expect({ result, error: null, errorType: null }).toEqual({
        result: row.result,
        error: row.error,
        errorType: row.errorType,
      });
      expect(coefficients).toEqual(original);
    } finally {
      setrand(saved);
    }
  }, 30000);
}
