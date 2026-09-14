import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";


import { expect, test } from 'bun:test';
import { execFileSync } from 'node:child_process';
import { FF_issquareall, PariType } from '@sagemath-ts/parigp-ts';
import {
  setrand as rootSetrand,
  getrand as rootGetrand,
} from '@sagemath-ts/parigp-ts';
import {
  FiniteFieldExtension as RootField,
  FiniteFieldElement as ExtensionElement,
  PrimeField,
} from './finite_field_extension.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
const native = (await loadLiveNative(import.meta.url, "./extension_square_root.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];
type Any = any;
const functions: Record<string, (...args: Any[]) => string> = {};
function invoke(x: object, name: string, ...args: unknown[]): unknown {
  return (x as Record<string, (...values: unknown[]) => unknown>)[name]!.apply(x, args);
}
functions.pari_ff_square_root = (p: bigint, T: bigint[], a: bigint[], seed: bigint) => {
  rootSetrand(seed);
  const result = FF_issquareall({
    type: PariType.t_FFELT,
    p,
    degree: T.length - 1,
    value: a,
    definingPoly: T,
  });
  return JSON.stringify({
    value:
      result === null
        ? null
        : (typeof result.value === 'bigint' ? [result.value] : result.value).map(String),
    state: String(rootGetrand()),
  });
};
functions.ff_extension_sqrt = (p: bigint, T: bigint[], a: bigint[], seed: bigint, mode: bigint) => {
  const K = new PrimeField(p),
    R = new PolynomialRing(K, 'x');
  const F = new RootField(p, T.length - 1, R.__call__(T.map((c) => K.__call__(c))), 'a');
  const x = F.__call__(R.__call__(a.map((c) => K.__call__(c))));
  const options =
    mode === 0n
      ? undefined
      : mode === 1n
        ? { extend: false }
        : mode === 2n
          ? { extend: true }
          : mode === 3n
            ? { all: true }
            : mode === 4n
              ? { extend: false, all: true }
              : mode === 5n
                ? { extend: true, all: true }
                : { unknown: true };
  rootSetrand(seed);
  try {
    const result = invoke(x, 'sqrt', ...(options === undefined ? [] : [options]));
    const encode = (v: ExtensionElement) => v.lift.coeffs.map((c) => String(c.value));
    return JSON.stringify({
      value: Array.isArray(result) ? result.map(encode) : encode(result as ExtensionElement),
      state: String(rootGetrand()),
    });
  } catch (e) {
    return JSON.stringify({
      error: (e as Error).name,
      message: (e as Error).message,
      state: String(rootGetrand()),
    });
  }
};

for (const [i, row] of native.entries()) {
  test(`native extension scalar root ${i} ${row.function}`, () => {
    const args = row.args.map((arg) =>
      arg.startsWith('[') ? JSON.parse(arg.replace(/-?\d+/g, '"$&"')) : arg
    );
    const decode = (v: Any): Any => (Array.isArray(v) ? v.map(decode) : BigInt(v));
    expect(functions[row.function]!(...args.map(decode))).toBe(row.result);
  });
}
// A separate process bounds extension searches and retains exact native roots
// and states, including binary degrees spanning more than two packed words.
for (const p of ['3', '18446744073709551557', '618970019642690137449562111', '2']) {
  const row = native.find(
    (r) =>
      r.function === 'ff_extension_sqrt' &&
      r.args[0] === p &&
      r.args[4] === '3' &&
      r.args[2] !== '[]' &&
      r.args[2] !== '[1]' &&
      (p !== '2' || JSON.parse(r.args[1]).length > 100)
  )!;
  test(`watched native extension sqrt ${p}`, () => {
    const args = row.args.map((arg) =>
      arg.startsWith('[') ? JSON.parse(arg.replace(/-?\d+/g, '"$&"')) : arg
    );
    const script = `
      import { FiniteFieldExtension as RootField, PrimeField } from '${import.meta.dir}/finite_field_extension.ts';
      import { PolynomialRing } from '${import.meta.dir}/../polynomial/polynomial_ring.ts';
      import { setrand as rootSetrand, getrand as rootGetrand } from '${import.meta.dir}/../../../../parigp-ts/src/random.ts';
      const invoke=${invoke.toString()};
      const decode=v=>Array.isArray(v)?v.map(decode):BigInt(v);
      console.log((${functions.ff_extension_sqrt!.toString()})(...${JSON.stringify(args)}.map(decode)));`;
    expect(
      execFileSync(process.execPath, ['--eval', script], {
        encoding: 'utf8',
        timeout: 10000,
      }).trim()
    ).toBe(row.result);
  }, 15000);
}

test('scalar roots retain input ownership and the original parent', () => {
  const K = new PrimeField(3n),
    R = new PolynomialRing(K, 'x');
  const F = new RootField(3n, 2, R.__call__([1n, 0n, 1n].map((c) => K.__call__(c))), 'a');
  const x = F.__call__(2n),
    before = x.lift.coeffs.map(String),
    state = rootGetrand();
  try {
    rootSetrand(17n);
    const roots = x.sqrt({ all: true });
    expect(roots).toHaveLength(2);
    expect(roots.every((r) => r.parent === F && r.pow(2n).eq(x))).toBe(true);
    expect(x.lift.coeffs.map(String)).toEqual(before);
  } finally {
    rootSetrand(state);
  }
});
