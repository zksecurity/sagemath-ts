import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { expect, test } from 'bun:test';
import { execFileSync } from 'node:child_process';
import {
  FF_issquareall,
  FpX_resultant,
  Flx_resultant,
  FpXQ_norm,
  Flxq_norm,
  FpXQ_issquare,
  Flxq_issquare,
  FF_norm,
  FF_issquare,
  PariType,
} from './index.js';
import { setrand as rootSetrand, getrand as rootGetrand } from './random.js';
const native = (await loadLiveNative(import.meta.url, "./field_predicates.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];
type Any = any;
const functions: Record<string, (...args: Any[]) => string> = {};
const decode = (v: Any): Any => (Array.isArray(v) ? v.map(decode) : BigInt(v));
const raw = (args: string[]) =>
  args.map((arg) => (arg.startsWith('[') ? JSON.parse(arg.replace(/-?\d+/g, '"$&"')) : arg));
functions.pari_field_predicates = (
  op: bigint,
  p: bigint,
  a: bigint[],
  T: bigint[],
  seed: bigint
) => {
  rootSetrand(seed);
  const x = {
    type: PariType.t_FFELT as PariType.t_FFELT,
    p,
    degree: T.length - 1,
    value: a,
    definingPoly: T,
  };
  const value =
    op === 0n
      ? FpX_resultant(a, T, p)
      : op === 1n
        ? Flx_resultant(a, T, p)
        : op === 2n
          ? FpXQ_norm(a, T, p)
          : op === 3n
            ? Flxq_norm(a, T, p)
            : op === 4n
              ? FpXQ_issquare(a, T, p)
              : op === 5n
                ? Flxq_issquare(a, T, p)
                : op === 6n
                  ? FF_norm(x)
                  : FF_issquare(x);
  return JSON.stringify({
    value: typeof value === 'boolean' ? value : String(value),
    state: String(rootGetrand()),
  });
};

for (const [i, row] of native.entries())
  test(`native field resultant/norm/predicate ${i}`, () => {
    expect(functions[row.function]!(...raw(row.args).map(decode))).toBe(row.result);
  });
for (const p of ['3', '2305843009213693951', '618970019642690137449562111']) {
  const row = native.find(
    (r) =>
      r.args[0] === '0' &&
      r.args[1] === p &&
      JSON.parse(r.args[2]).length > (p === '3' ? 852 : p === '2305843009213693951' ? 1874 : 388) &&
      JSON.parse(r.result).value !== '0'
  )!;
  test(`watched native half-resultant ${p}`, () => {
    const script = `
      import {FpX_resultant,Flx_resultant,FpXQ_norm,Flxq_norm,FpXQ_issquare,Flxq_issquare,FF_norm,FF_issquare,PariType} from '${import.meta.dir}/index.ts';
      import {setrand as rootSetrand,getrand as rootGetrand} from '${import.meta.dir}/random.ts';
      const decode=${decode.toString()};
      console.log((${functions.pari_field_predicates!.toString()})(...${JSON.stringify(raw(row.args))}.map(decode)));`;
    expect(
      execFileSync(process.execPath, ['--eval', script], {
        encoding: 'utf8',
        timeout: 10000,
      }).trim()
    ).toBe(row.result);
  }, 15000);
}
test('half-resultant calls preserve both polynomial inputs', () => {
  const row = native.find(
    (r) => r.args[0] === '0' && r.args[1] === '3' && JSON.parse(r.args[2]).length > 852
  )!;
  const args = raw(row.args).map(decode),
    snapshot = structuredClone(args);
  expect(functions.pari_field_predicates!(...args)).toBe(row.result);
  expect(args).toEqual(snapshot);
});

functions.pari_field_record_scalar=(op:bigint,p:bigint,a:bigint,seed:bigint)=>{
  rootSetrand(seed);
  const x={type:PariType.t_FFELT as PariType.t_FFELT,p,degree:1,value:a};
  if(op===0n){
    const root=FF_issquareall(x);
    return JSON.stringify({value:root===null?null:(typeof root.value==='bigint'?[root.value]:root.value).map(String),state:String(rootGetrand())});
  }
  const value=op===1n?FF_norm(x):FF_issquare(x);
  return JSON.stringify({value:typeof value==='boolean'?value:String(value),state:String(rootGetrand())});
};
