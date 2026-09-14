import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import {test,expect}from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./general_hnf.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];
import { ZM_hnf, ZM_hnfall_i, ZM_hnfall, hnfall } from './hnf_snf.js';
export function pari_general_hnf(op:bigint,u:bigint,remove:bigint,m:bigint,n:bigint,flat:bigint[]):string {
 const A=Array.from({length:Number(n)},(_,j)=>Array.from({length:Number(m)},(_,i)=>flat[i*Number(n)+j]!));
 const r=op===0n?ZM_hnf(A):op===4n?hnfall(A):op===2n?ZM_hnfall_i(A,!!u,Number(remove)as 0|1|2):ZM_hnfall(A,!!u,Number(remove)as 0|1|2);
 return JSON.stringify(r,(_,v)=>typeof v==='bigint'?String(v):v);
}

const arg=(s:string):bigint|bigint[]=>s.startsWith('[')?s.slice(1,-1).trim()?s.slice(1,-1).split(',').map(v=>BigInt(v.trim())):[]:BigInt(s);
for(const[i,row]of fixtures.entries())test('bundled PARI general Hermite '+i,()=>{
 let result:string|null=null,error:string|null=null,errorType:string|null=null;
 try{result=(pari_general_hnf as any)(...row.args.map(arg));}catch(e){error=(e as Error).message;errorType=(e as Error).name;}
 expect({result,error,errorType}).toEqual({result:row.result,error:row.error,errorType:row.errorType});
},10000);
test('general Hermite output owns its matrix and transformation columns',()=>{
 const A=[[1n,2n],[2n,4n],[1n,0n]],before=structuredClone(A);
 for(const[H,U]of [ZM_hnfall(A,true,0),ZM_hnfall_i(A,true,1),ZM_hnfall(A,true,2),hnfall(A)]) {
   H[0]![0]=123n;U![0]![0]=124n;expect(A).toEqual(before);
 }
 const H=ZM_hnf(A);H[0]![0]=125n;expect(A).toEqual(before);
});
