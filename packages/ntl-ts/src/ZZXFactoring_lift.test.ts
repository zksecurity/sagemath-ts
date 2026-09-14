import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import {test,expect} from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./ZZXFactoring_lift.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import * as q from './ZZXFactoring.js';
import * as z from './ZZX1.js';
export function ntl_multilift(op:bigint,p:bigint,e:bigint,F:bigint[],A:bigint[],B:bigint[]):string {
 let x:unknown;
 if(op===0n)x=z.mul(A,B);
 else if(op===1n)x=z.sqr(A);
 else{const a:bigint[][]=[];for(let i=0,k=1;i<Number(A[0]);i++){const n=Number(A[k++]!);a.push(A.slice(k,k+n));k+=n}x=q.MultiLift(a,F,Number(e),p)}
 return JSON.stringify(x,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}

const arg=(s:string):bigint|bigint[]=>s.startsWith('[')?s.slice(1,-1).trim()?s.slice(1,-1).split(',').map(v=>BigInt(v.trim())):[]:BigInt(s);
for(const[i,row]of fixtures.entries())test('native NTL integer products and lifting '+i,()=>{
 let result:string|null=null,error:string|null=null,errorType:string|null=null;
 try{result=(ntl_multilift as Function)(...row.args.map(arg));}catch(e){error=(e as Error).message;errorType=(e as Error).name;}
 expect({result,error,errorType}).toEqual({result:row.result,error:row.error,errorType:row.errorType});
},10000);
test('integer products preserve normalization and alias-square input ownership',()=>{
 const a=[-2n,0n,3n,0n],b=[5n,-7n,1n];const p=z.mul(a,b),s=z.mul(a,a);
 expect(p).toEqual([-10n,14n,13n,-21n,3n]);expect(s).toEqual(z.sqr(a));expect(s).toEqual([4n,0n,-12n,0n,9n]);
 p[0]=999n;s[0]=999n;expect(a).toEqual([-2n,0n,3n,0n]);expect(b).toEqual([5n,-7n,1n]);
});
test('multifactor lifting preserves inputs and owns every returned factor',()=>{
 const a=[[1n,1n],[2n,1n]],f=[7n,3n,1n];const one=q.MultiLift(a,f,1,5n),lifted=q.MultiLift(a,f,4,5n);
 expect(one).toEqual(a);expect(lifted).toEqual([[281n,1n],[347n,1n]]);one[0]![0]=99n;lifted[1]![0]=99n;
 expect(a).toEqual([[1n,1n],[2n,1n]]);expect(f).toEqual([7n,3n,1n]);
});
