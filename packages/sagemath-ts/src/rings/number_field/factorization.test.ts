import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const cases = (await loadLiveNative(import.meta.url, "./factorization.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { NumberField, RationalPolynomial } from '../../rings/number_field/number_field.js';
import { NumberFieldIdeal } from '../../rings/number_field/number_field_ideal.js';
import { Rational } from '../../rings/rational.js';
import { hnf } from '../../rings/number_field/pari_nf.js';
import { lcm } from '../../arith/misc.js';
function lattice(I: NumberFieldIdeal) {
  if (I.is_zero()) return [];
  const rows = I.zk_basis().map((x) => x.list());
  let d = 1n;
  for (const row of rows) for (const c of row) d = lcm(d, c.denominator);
  return hnf(
    rows.map((row) => row.map((c) => c.numerator * (d / c.denominator))),
    I.number_field().degree()
  ).map((row) => row.map((c) => String(new Rational(c, d))));
}

export function nf_ideal_factor(
  op: bigint,
  cs: bigint[],
  flat: bigint[],
  d: bigint,
  count: bigint
): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    n = K.degree(),
    gens = Array.from({ length: Number(count) }, (_, i) =>
      K.__call__(flat.slice(i * n, (i + 1) * n).map((c) => new Rational(c, d)))
    ),
    I = K.ideal(gens);
  const F = op === 0n || op === 3n ? I.factor() : op === 1n ? K.factor(I) : K.factor(gens);
  const value = F.map(([P, e]) => [
    op === 3n
      ? (P as any)._computeHNF().entries.map((col: bigint[]) => col.map(String))
      : lattice(P),
    String(e),
  ]);
  if (op !== 3n)
    value.sort((a, b) =>
      JSON.stringify(a[0]) < JSON.stringify(b[0])
        ? -1
        : JSON.stringify(a[0]) > JSON.stringify(b[0])
          ? 1
          : 0
    );
  return JSON.stringify(value);
}

import { NumberFieldFractionalIdeal } from '../../rings/number_field/number_field_ideal.js';
export function nf_factor_properties(
  cs: bigint[],
  flat: bigint[],
  d: bigint,
  count: bigint
): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    n = K.degree(),
    I = K.ideal(
      Array.from({ length: Number(count) }, (_, i) =>
        K.__call__(flat.slice(i * n, (i + 1) * n).map((c) => new Rational(c, d)))
      )
    );
  const F = I.factor();
  let product = K.ideal(1n);
  const rows = [];
  for (const [P, e] of F) {
    product = product.mul(P.pow(e));
    rows.push([
      lattice(P),
      String(e),
      P.is_prime(),
      P instanceof NumberFieldFractionalIdeal,
      String(P.ramification_index()),
      String(P.residue_class_degree()),
      String(I.valuation(P)),
    ]);
  }
  rows.sort((a, b) =>
    JSON.stringify(a[0]) < JSON.stringify(b[0])
      ? -1
      : JSON.stringify(a[0]) > JSON.stringify(b[0])
        ? 1
        : 0
  );
  return JSON.stringify([F === I.factor(), F === K.factor(I), lattice(product), rows]);
}

const functions = { nf_ideal_factor, nf_factor_properties };
const arg = (s: string): bigint | bigint[] =>
  s.startsWith('[')
    ? s.slice(1, -1).trim()
      ? s
          .slice(1, -1)
          .split(',')
          .map((v) => BigInt(v.trim()))
      : []
    : BigInt(s);
for (const [i, row] of cases.filter((row) => Object.hasOwn(functions, row.function)).entries())
  test('bundled/Sage full ideal ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = (functions as any)[row.function](...row.args.map(arg));
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
  });
