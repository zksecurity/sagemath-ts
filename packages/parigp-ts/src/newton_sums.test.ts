import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import {test,expect}from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./newton_sums.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { polsym, polsym_gen } from './polarit2.js';
export function pari_newton_sums(
  op: bigint,
  n: bigint,
  N: bigint,
  use: bigint,
  P: bigint[],
  c: bigint[]
): string {
  let result: unknown;
  if (op) result = polsym(P, Number(n));
  else if (N) result = polsym_gen(P, use ? c : null, Number(n), null, N);
  else {
    const y = use
      ? Array.from(
          { length: c.length / 2 },
          (_, i) => [c[2 * i]!, c[2 * i + 1]!] as [bigint, bigint]
        )
      : null;
    result = polsym_gen(P, y, Number(n), null, null);
  }
  return JSON.stringify(result, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}

const invoke:Record<string,Function>={pari_newton_sums};

const arg=(s:string):bigint|bigint[]=>s.startsWith('[')?s.slice(1,-1).trim()?s.slice(1,-1).split(',').map(v=>BigInt(v.trim())):[]:BigInt(s);
for(const[i,row]of fixtures.entries())test('native integer factorization '+row.function+' '+i,()=>{
 let result:string|null=null,error:string|null=null,errorType:string|null=null;
 try{result=(invoke[row.function] as any)(...row.args.map(arg));}catch(e){error=(e as Error).message;errorType=(e as Error).name;}
 expect({result,error,errorType}).toEqual({result:row.result,error:row.error,errorType:row.errorType});
},10000);

test('Newton sum output preserves polynomial and supplied prefix storage',()=>{
 const P=[1n,0n,2n],prefix:[bigint,bigint][]=[[2n,1n],[0n,1n]],before=structuredClone([P,prefix]);
 const r=polsym_gen(P,prefix,4,null,null);r[0]![0]=99n;expect([P,prefix]).toEqual(before);
 const cached=[2n,0n],s=polsym_gen(P,cached,4,null,7n);s[0]=98n;expect(cached).toEqual([2n,0n]);
});
