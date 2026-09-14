import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import {test,expect} from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./bounded_recombination.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { ZX_divides_i, ZX_divides, cmbf_precs, cmbf } from './QX_factor.js';
import { FpXV_prod } from './FpX.js';
import { centermodii } from './polarit2.js';
import { cmbf_maxK } from './nffactor.js';
export function pari_recombination(op: bigint, p: bigint, a: bigint, b: bigint, k: bigint, B: bigint, x: bigint[], f: bigint[]): string {
  const factors=()=>{const out:bigint[][]=[];for(let i=0;i<f.length;){const n=Number(f[i++]!);out.push(f.slice(i,i+n));i+=n;}return out;};
  let result:unknown;
  if(op===0n)result=ZX_divides_i(x,f,B<0n?null:B);
  else if(op===6n)result=ZX_divides(x,f);
  else if(op===1n)result=cmbf_precs(p,x[0]!,B);
  else if(op===2n){const r=cmbf(x,factors(),B,p,Number(a),Number(b),Number(k));result=[r[0],r[1],r[2],r[3]?1:0];}
  else if(op===3n)result=centermodii(x[0]!,p,a?B:null);
  else if(op===4n)result=FpXV_prod(factors(),p);
  else result=cmbf_maxK(Number(a));
  return JSON.stringify(result,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}

const arg=(s:string):bigint|bigint[]=>s.startsWith('[')?s.slice(1,-1).trim()?s.slice(1,-1).split(',').map(v=>BigInt(v.trim())):[]:BigInt(s);
for(const[i,row]of fixtures.entries())test('bundled PARI bounded recombination '+i,()=>{
 let result:string|null=null,error:string|null=null,errorType:string|null=null;
 try{result=(pari_recombination as any)(...row.args.map(arg));}catch(e){error=(e as Error).message;errorType=(e as Error).name;}
 expect({result,error,errorType}).toEqual({result:row.result,error:row.error,errorType:row.errorType});
},10000);
test('bounded recombination and division preserve supplied storage',()=>{
 const pol=[-1n,0n,1n],factors=[[1n,1n],[124n,1n]], before=structuredClone([pol,factors]);
 const result=cmbf(pol,factors,3n,5n,3,1,1);
 expect([pol,factors]).toEqual(before);
 result[0][0]![0]=99n;result[1][0]![0]![0]=98n;
 expect([pol,factors]).toEqual(before);
 const divisor=[1n,1n],q=ZX_divides(pol,divisor)!;q[0]=97n;
 expect(pol).toEqual([-1n,0n,1n]);expect(divisor).toEqual([1n,1n]);
});
test('singleton modular products return a copy without normalization',()=>{
 const f=[-5n,0n,9n],product=FpXV_prod([f],3n) as bigint[];
 expect(product).toEqual(f);product[0]=101n;expect(f).toEqual([-5n,0n,9n]);
});
