import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { expect, test } from 'bun:test';
import { FiniteFieldExtension as RootField, PrimeField } from './finite_field_extension.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import { setrand as rootSetrand, getrand as rootGetrand } from '@sagemath-ts/parigp-ts';
const native = await loadLiveNative(import.meta.url, "./extension_is_square.native.json");
type Any = any;
const functions: Record<string, (...args: Any[]) => string> = {};
functions.ff_extension_is_square = (p: bigint, T: bigint[], a: bigint[], seed: bigint) => {
  const K = new PrimeField(p),
    R = new PolynomialRing(K, 'x');
  const F = new RootField(p, T.length - 1, R.__call__(T.map((c) => K.__call__(c))), 'a');
  const x = F.__call__(R.__call__(a.map((c) => K.__call__(c))));
  rootSetrand(seed);
  return JSON.stringify({ value: x.is_square(), state: String(rootGetrand()) });
};

for (const [i, row] of native.entries())
  test(`native extension square predicate ${i}`, () => {
    const raw = row.args.map((arg) =>
      arg.startsWith('[') ? JSON.parse(arg.replace(/-?\d+/g, '"$&"')) : arg
    );
    const decode = (v: Any): Any => (Array.isArray(v) ? v.map(decode) : BigInt(v));
    expect(functions.ff_extension_is_square!(...raw.map(decode))).toBe(row.result);
  });
