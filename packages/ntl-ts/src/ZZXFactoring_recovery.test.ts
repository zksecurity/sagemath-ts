import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import {test,expect} from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./ZZXFactoring_recovery.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import * as q from './ZZXFactoring.js';
export function ntl_factor_recovery(op:bigint,v:bigint[],F:bigint[],D:bigint[],flat:bigint[],W:bigint[]):string{
 let x:unknown;
 if(op===0n)x=q.PolyEval(F,v[0]!);
 else if(op===1n)x=q.RootBound(F);
 else if(op===2n){const m=Number(v[0]),n=Number(v[1]),C=Number(v[2]),r=Number(v[3]),d=Number(v[4]),M=Array.from({length:m},(_,i)=>flat.slice(i*n,(i+1)*n));x=q.CutAway(D,M,C!,r!,d!)}
 else{const w:bigint[][]=[];for(let i=0,k=1;i<Number(W[0]);i++){const n=Number(W[k++]!);w.push(W.slice(k,k+n));k+=n}x=q.AdditionalLifting(v[0]!,Number(v[1]),w,v[2]!,Number(v[3]),F,Boolean(v[4]))}
 return JSON.stringify(x,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}

const arg=(s:string):bigint|bigint[]=>s.startsWith('[')?s.slice(1,-1).trim()?s.slice(1,-1).split(',').map(v=>BigInt(v.trim())):[]:BigInt(s);
for(const[i,row]of fixtures.entries())test('native NTL factor recovery '+i,()=>{
 let result:string|null=null,error:string|null=null,errorType:string|null=null;
 try{result=(ntl_factor_recovery as Function)(...row.args.map(arg));}catch(e){error=(e as Error).message;errorType=(e as Error).name;}
 expect({result,error,errorType}).toEqual({result:row.result,error:row.error,errorType:row.errorType});
},10000);
test('NTL factor row removal preserves determinant and matrix inputs',()=>{
 const D=[1n,1n,1n,100n],M=[[1n,0n,0n],[0n,1n,0n],[0n,0n,10n]],copy=M.map(r=>r.slice());const B=q.CutAway(D,M,1,2,1);
 expect(B).toEqual([[1n,0n],[0n,1n]]);B[0]![0]=99n;expect(M).toEqual(copy);expect(D).toEqual([1n,1n,1n,100n]);
});
test('NTL additional lifting returns independent factors and ignores the old modulus slot',()=>{
 const W=[[1n,1n],[2n,1n]],F=[7n,3n,1n];const next=q.AdditionalLifting(-19n,1,W,5n,4,F,false);
 expect(next).toEqual([625n,4,[[281n,1n],[347n,1n]]]);next[2][0]![0]=99n;expect(W).toEqual([[1n,1n],[2n,1n]]);expect(F).toEqual([7n,3n,1n]);
});
