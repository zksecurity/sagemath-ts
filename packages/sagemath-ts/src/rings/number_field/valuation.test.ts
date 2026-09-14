import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';
const cases = await loadLiveNative(import.meta.url, "./valuation.native.json");
import { NumberField, RationalPolynomial } from '../../rings/number_field/number_field.js';
import { NumberFieldIdeal } from '../../rings/number_field/number_field_ideal.js';
import { Rational } from '../../rings/rational.js';
import { ValueError } from '../../errors.js';
import { Integer } from '../../rings/integer_ring.js';
import { hnf } from '../../rings/number_field/pari_nf.js';
import { lcm } from '../../arith/misc.js';
function lattice(I: NumberFieldIdeal): string[][] {
  if (I.is_zero()) return [];
  const rows = I.zk_basis().map((x) => x.list());
  let d = 1n;
  for (const row of rows) for (const c of row) d = lcm(d, c.denominator);
  return hnf(
    rows.map((row) => row.map((c) => c.numerator * (d / c.denominator))),
    I.number_field().degree()
  ).map((row) => row.map((c) => String(new Rational(c, d))));
}
export function nf_ideal_valuation(
  cs: bigint[],
  p: bigint,
  flat: bigint[],
  d: bigint,
  count: bigint,
  base: bigint
): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    n = K.degree(),
    gens = Array.from({ length: Number(count) }, (_, i) =>
      K.__call__(flat.slice(i * n, (i + 1) * n).map((c) => new Rational(c, d)))
    ),
    I = base ? new NumberFieldIdeal(K, gens) : K.ideal(gens);
  const out = K.decomposition(p).map(([P]) => {
    const v = I.valuation(P as any);
    return [lattice(P), v === 'Infinity' ? '+Infinity' : String(v)];
  });
  out.sort((a, b) =>
    JSON.stringify(a[0]) < JSON.stringify(b[0])
      ? -1
      : JSON.stringify(a[0]) > JSON.stringify(b[0])
        ? 1
        : 0
  );
  return JSON.stringify(out);
}
export function nf_ideal_valuation_input(
  cs: bigint[],
  receiver: bigint,
  kind: bigint,
  n: bigint,
  d: bigint
): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    L = new NumberField(RationalPolynomial.fromBigInts(cs), 'b'),
    q = new Rational(n, d);
  const I = [
    K.ideal(0n),
    K.ideal(1n),
    K.ideal(6n),
    K.ideal(K.gen().add(1n)),
    new NumberFieldIdeal(K, [K.zero()]),
    new NumberFieldIdeal(K, [K.__call__(6n)]),
  ][Number(receiver)]!;
  const values = [
    n,
    new Integer(n),
    q,
    Number(n) / Number(d),
    K.__call__(q),
    K.gen().add(q),
    [],
    [q],
    Array(K.degree()).fill(q),
    K.ideal(q),
    K.ideal(0n),
    new NumberFieldIdeal(K, [K.__call__(q)]),
    L.__call__(q),
    L.gen().add(q),
    L.ideal(q),
    L.decomposition(2n)[0]![0],
    '2',
    [K.gen(), q],
  ];
  try {
    const v = I.valuation(values[Number(kind)] as any);
    return v === 'Infinity' ? '+Infinity' : String(v);
  } catch (error) {
    if (
      error instanceof ValueError &&
      (kind === 8n || kind === 17n) &&
      error.message.startsWith('p (= ')
    ) {
      const P = K.ideal(values[Number(kind)] as any);
      throw new ValueError(
        'p (= ' +
          JSON.stringify(lattice(P)) +
          ')' +
          error.message.slice(error.message.lastIndexOf(')') + 1)
      );
    }
    throw error;
  }
}

export function nf_prime_below(cs: bigint[], p: bigint): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    out = K.decomposition(p).map(([P]) => [lattice(P), String(P.prime_below())]);
  out.sort((a, b) =>
    JSON.stringify(a[0]) < JSON.stringify(b[0])
      ? -1
      : JSON.stringify(a[0]) > JSON.stringify(b[0])
        ? 1
        : 0
  );
  return JSON.stringify(out);
}

const functions = { nf_ideal_valuation, nf_ideal_valuation_input, nf_prime_below };
const arg = (s: string): bigint | bigint[] =>
  s.startsWith('[') ? JSON.parse(s.replace(/-?\d+/g, '"$&"')).map(BigInt) : BigInt(s);
for (const [i, r] of cases.entries())
  test('native valuation ' + i, () => {
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

import { execFileSync } from 'node:child_process';
import hangs from './valuation.hangs.native.json' with { type: 'json' };
test('former fractional and large-prime hangs complete with native results', () => {
  const program = `import {NumberField,RationalPolynomial}from ${JSON.stringify(import.meta.dir + '/number_field.js')};
 const K=new NumberField(RationalPolynomial.fromBigInts([-2n,0n,1n]),'a'),P=K.decomposition(2n)[0][0],x=K.gen().add(1n).div(K.__call__(2n));
 let rejected;try{P.valuation(x);}catch(e){rejected={error:e.message,errorType:e.name};}
 console.log(JSON.stringify({rejected,valuation:String(K.ideal(x).valuation(P)),prime:String(K.ideal(18446744073709551629n).prime_below())}));`;
  const actual = JSON.parse(
    execFileSync(process.execPath, ['-e', program], { encoding: 'utf8', timeout: 5000 })
  );
  expect(actual).toEqual({
    rejected: hangs.native.original_argument,
    valuation: hangs.native.ideal_valuation.result,
    prime: '18446744073709551629',
  });
}, 10000);
