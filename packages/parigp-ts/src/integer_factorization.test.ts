import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import {test,expect}from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./integer_factorization.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];
import * as q from './QX_factor.js';
import { setrand, getrand } from './random.js';
export function pari_qx_factor(
  op: bigint,
  p: bigint,
  a: bigint,
  b: bigint,
  k: bigint,
  B: bigint,
  P: bigint[],
  F: bigint[],
  M: bigint[]
): string {
  const unpack = (f: bigint[]) => {
    const out: bigint[][] = [];
    for (let i = 0; i < f.length; ) {
      const n = Number(f[i++]!);
      out.push(f.slice(i, i + n));
      i += n;
    }
    return out;
  };
  const saved = getrand();
  try {
    setrand(1n);
    let x: unknown;
    if (op === 0n) x = q.ZX_gcd_all(P, F);
    else if (op === 1n) x = q.ZX_squff(P);
    else if (op === 2n) x = q.pick_prime(P, Number(k));
    else if (op === 3n) x = q.ZX_DDF_max(P, Number(k));
    else if (op === 4n) x = q.ZX_DDF(P);
    else if (op === 5n || op === 12n) {
      const f = op === 5n ? q.ZX_factor(P) : q.QX_factor([P, B]);
      x = [f.map(([p]) => p), f.map(([, e]) => e)];
    } else if (op === 6n) x = q.ZX_is_irred(P) ? 1 : 0;
    else if (op === 7n) x = q.combine_factors(P, unpack(F), p, Number(k));
    else if (op === 8n) x = q.LLL_cmbf(P, unpack(F), p, p ** a, B, Number(a), Number(b));
    else if (op === 9n) x = q.chk_factors(P, unpack(M), B, unpack(F), p ** a);
    else if (op === 10n) x = q.chk_factors_get(B || null, unpack(F), M, null, p);
    return JSON.stringify([x, getrand()], (_, v) =>
      typeof v === 'bigint' || typeof v === 'number' ? String(v) : v
    );
  } finally {
    setrand(saved);
  }
}

const invoke:Record<string,Function>={pari_qx_factor};

const arg=(s:string):bigint|bigint[]=>s.startsWith('[')?s.slice(1,-1).trim()?s.slice(1,-1).split(',').map(v=>BigInt(v.trim())):[]:BigInt(s);
for(const[i,row]of fixtures.entries())test('native integer factorization '+row.function+' '+i,()=>{
 let result:string|null=null,error:string|null=null,errorType:string|null=null;
 try{result=(invoke[row.function] as any)(...row.args.map(arg));}catch(e){error=(e as Error).message;errorType=(e as Error).name;}
 expect({result,error,errorType}).toEqual({result:row.result,error:row.error,errorType:row.errorType});
},10000);

test('integer factors and selected modular products own their coefficient arrays',()=>{
 const saved=getrand();try {const P=[-3n,0n,0n,0n,3n],before=P.slice(),f=q.ZX_factor(P);f[0]![0][0]=99n;expect(P).toEqual(before);
 const modular=[[1n,1n],[4n,1n]],v=q.chk_factors_get(null,modular,[1n,0n],null,5n)as bigint[];
 v[0]=98n;expect(modular).toEqual([[1n,1n],[4n,1n]]);
 }finally{setrand(saved);}
});
