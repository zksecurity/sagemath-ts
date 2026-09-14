import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';
const cases = await loadLiveNative(import.meta.url, "./trace_norm.native.json");
import { NumberField, RationalPolynomial } from '../../rings/number_field/number_field.js';
import { Rational } from '../../rings/rational.js';
export function nf_trace_norm(op: bigint, cs: bigint[], flat: bigint[], d: bigint): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a', undefined, false),
    a = K.__call__(flat.map((c) => new Rational(c, d)));
  return String(
    op === 0n
      ? a.trace()
      : op === 1n
        ? a.norm()
        : op === 2n
          ? a.absolute_norm()
          : op === 3n
            ? a.absolute_trace()
            : op === 4n
              ? a.relative_norm()
              : a.relative_trace()
  );
}

const arg = (s: string): bigint | bigint[] =>
  s.startsWith('[')
    ? s.slice(1, -1).trim()
      ? s
          .slice(1, -1)
          .split(',')
          .map((v) => BigInt(v.trim()))
      : []
    : BigInt(s);
for (const [i, row] of cases.entries())
  test('original trace/norm trace_norm ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = (nf_trace_norm as any)(...row.args.map(arg));
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

import { execFileSync } from 'node:child_process';
import bounds from './trace_norm.bounds.native.json' with { type: 'json' };
for (const row of bounds.native)
  test('native value without quartic trace/norm work ' +
    row.op +
    ':' +
    row.degree +
    ':' +
    row.kind, () => {
    // x^n-2 is irreducible by Eisenstein. Disable both constructors' redundant
    // irreducibility check to isolate the operation being compared.
    const program = `import {NumberField,RationalPolynomial}from ${JSON.stringify(import.meta.dir + '/number_field.js')};
  const n=${row.degree},K=new NumberField(RationalPolynomial.fromBigInts([-2n,...Array(n-1).fill(0n),1n]),'a',undefined,false),x=${row.kind}===0?K.zero():${row.kind}===1?K.__call__(3n):K.gen().add(1n);
  console.log(String(${row.op}===0?x.trace():x.norm()));`;
    expect(
      execFileSync(process.execPath, ['-e', program], { encoding: 'utf8', timeout: 3000 }).trim()
    ).toBe(row.result);
  });
