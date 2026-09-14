import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';
const cases = await loadLiveNative(import.meta.url, "./prime_valuation.native.json");
import { NumberField, RationalPolynomial } from '../../rings/number_field/number_field.js';
import { Rational } from '../../rings/rational.js';
import { primedec_end } from '@sagemath-ts/parigp-ts/src/base2.js';
import { idealval } from '@sagemath-ts/parigp-ts/src/base4.js';
const encode = (v: unknown) => JSON.stringify(v, (_, x) => (typeof x === 'bigint' ? String(x) : x));
const compare = (A: bigint[][], B: bigint[][]) => {
  for (let j = 0; j < A.length; j++)
    for (let i = 0; i < A[j]!.length; i++) {
      if (A[j]![i]! < B[j]![i]!) return -1;
      if (A[j]![i]! > B[j]![i]!) return 1;
    }
  return 0;
};
function data(cs: bigint[], p: bigint) {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    nf = K._pari_ideal_data(),
    H = K.decomposition(p)
      .map(([P]) => (P as any)._computeHNF().entries as bigint[][])
      .sort(compare);
  return { K, nf, H, primes: primedec_end(nf, H, p) };
}
export function pari_prime_valuation_data(cs: bigint[], p: bigint): string {
  const { H, primes } = data(cs, p);
  return encode(primes.map((P, i) => [H[i], [P.p, P.generator, P.e, P.f, P.tau]]));
}
export function pari_ideal_valuation(
  cs: bigint[],
  p: bigint,
  flat: bigint[],
  d: bigint,
  count: bigint
): string {
  const { K, nf, H, primes } = data(cs, p),
    n = K.degree(),
    gens = Array.from({ length: Number(count) }, (_, i) =>
      K.__call__(flat.slice(i * n, (i + 1) * n).map((c) => new Rational(c, d)))
    ),
    I = K.ideal(gens),
    h = I.is_zero() ? { entries: [], denominator: 1n } : (I as any)._computeHNF();
  return encode(
    primes.map((P, i) => {
      const v = idealval(nf, h.entries, h.denominator, P);
      return [H[i], v === 'Infinity' ? '+oo' : v];
    })
  );
}

import { ZC_nfval } from '@sagemath-ts/parigp-ts/src/base3.js';
export function pari_prime_element_valuation(cs: bigint[], p: bigint, flat: bigint[]): string {
  const { H, primes } = data(cs, p);
  return encode(primes.map((P, i) => [H[i], ZC_nfval(flat, P)]));
}

const functions = { pari_prime_valuation_data, pari_ideal_valuation, pari_prime_element_valuation };
const arg = (s: string): bigint | bigint[] =>
  s.startsWith('[') ? JSON.parse(s.replace(/-?\d+/g, '"$&"')).map(BigInt) : BigInt(s);
for (const [i, r] of cases.entries())
  test('native prime_valuation ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = (functions as Record<string, (...args: any[]) => string>)[r.function]!(
        ...r.args.map(arg)
      );
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: r.result,
      error: r.error,
      errorType: r.errorType,
    });
  });

import resource from './valuation.resource.native.json' with { type: 'json' };
test('native zero-vector nontermination has an explicit low-level boundary', () => {
  const { primes } = data([-3n, 0n, 1n], 7n);
  expect(resource.native.error).toBe('bundled PARI element timed out; no comparison result');
  expect(() => ZC_nfval([0n, 0n], primes[0]!)).toThrow(resource.port.error!);
});
test('valuation preserves its integral inputs and rejects malformed adapter buffers', () => {
  const { nf, primes } = data([-2n, 0n, 1n], 2n),
    P = primes[0]!,
    H = [
      [4n, 0n],
      [0n, 4n],
    ],
    x = [8n, 8n];
  const saved = structuredClone({ nf, P, H, x });
  expect(idealval(nf, H, 8n, P)).toBe(-2n);
  expect(ZC_nfval(x, P)).toBe(6n);
  expect({ nf, P, H, x }).toEqual(saved);
  expect(() => idealval(nf, H, 0n, P)).toThrow('idealval requires a positive denominator');
  expect(() => idealval(nf, [[1n]], 1n, P)).toThrow('idealval requires a square ideal HNF');
});
