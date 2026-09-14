import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';
const fixtures = await loadLiveNative(import.meta.url, "./ideal_classes.native.json");
import { NumberField, RationalPolynomial } from './number_field.js';
import { NumberFieldIdeal } from './number_field_ideal.js';
import { Rational } from '../rational.js';
import { Integer } from '../integer_ring.js';
import { hnf } from './pari_nf.js';
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
function nf_ideal_class_method(
  cs: bigint[],
  factory: bigint,
  n: bigint,
  d: bigint,
  op: bigint
): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    q = new Rational(n, d),
    I = factory ? K.ideal(q) : new NumberFieldIdeal(K, [K.__call__(q)]);
  let r: unknown;
  if (op === 0n) r = I.inverse();
  else if (op === 1n) r = I.numerator();
  else if (op === 2n) r = I.denominator();
  else if (op === 3n) r = I.divides(K.ideal(6n));
  else if (op === 4n) r = I.is_coprime(K.ideal(3n));
  else if (op === 5n) r = I.pow(-1n);
  else if (op === 6n) r = I.pow(new Integer(-1n));
  else if (op === 7n) r = I.div(K.ideal(2n));
  else if (op === 8n) r = I.div(2n);
  else if (op === 9n) r = I.pow(2n);
  else if (op === 10n) r = I.factor().map(([P, e]) => [String(P.norm()), String(e)]);
  else if (op === 11n) r = I.ramification_index();
  else if (op === 12n) r = I.residue_class_degree();
  else if (op === 13n) {
    const F = I.residue_field();
    r = [String(F.characteristic), String(F.order), String(F.degree)];
  } else if (op === 14n) r = I.is_maximal();
  else if (op === 15n) r = I.is_integral();
  else if (op === 16n) r = I.pow(0n);
  else r = I.pow(1n);
  if (r instanceof NumberFieldIdeal) r = [r.constructor.name, lattice(r)];
  else if (typeof r === 'bigint') r = String(r);
  return JSON.stringify(r);
}
function nf_prime_ideal_class(cs: bigint[], p: bigint, op: bigint): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    out: unknown[][] = [];
  for (const [I] of K.decomposition(p)) {
    let r: unknown;
    if (op === 0n) r = I.constructor.name;
    else if (op === 1n) {
      const inv = I.inverse();
      r = [inv.constructor.name, String(inv.norm())];
    } else if (op === 2n) r = [String(I.ramification_index()), String(I.residue_class_degree())];
    else if (op === 3n) r = I.divides(p);
    else if (op === 4n) r = I.is_coprime(p);
    else r = [lattice(I.numerator()), lattice(I.denominator())];
    out.push([lattice(I), r]);
  }
  out.sort((a, b) =>
    JSON.stringify(a[0]) < JSON.stringify(b[0])
      ? -1
      : JSON.stringify(a[0]) > JSON.stringify(b[0])
        ? 1
        : 0
  );
  return JSON.stringify(out);
}

const operations = { nf_ideal_class_method, nf_prime_ideal_class };
const argument = (s: string) =>
  s.startsWith('[')
    ? s.slice(1, -1).trim()
      ? s
          .slice(1, -1)
          .split(',')
          .map((v) => BigInt(v.trim()))
      : []
    : BigInt(s);
for (const [i, f] of fixtures.entries())
  test(`original ideal class and prime method ${i}`, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = (operations[f.function as keyof typeof operations] as (...args: any[]) => string)(
        ...f.args.map(argument)
      );
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: f.result,
      error: f.error,
      errorType: f.errorType,
    });
  });
