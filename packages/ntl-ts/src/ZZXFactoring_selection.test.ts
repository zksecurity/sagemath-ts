import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import {test,expect} from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./ZZXFactoring_selection.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import * as q from './ZZXFactoring.js';
export function ntl_factor_selection(op:bigint,p:bigint,A:bigint[],B:bigint[],C:bigint[],packed:bigint[]):string{
 const W:bigint[][]=[];for(let i=0,k=1;i<Number(packed[0]);i++){const n=Number(packed[k++]!);W.push(packed.slice(k,k+n));k+=n;}
 const a=A.map(Number),b=B.map(Number);let x:unknown;
 switch(Number(op)){
 case 0:x=q.inplace_rev(A);break;
 case 1:x=q.RecordPattern(a,W.map((f,i)=>[f,b[i]!]as const),p);break;
 case 2:x=q.NumFactors(a);break;
 case 3:x=q.CalcPossibleDegrees(a);break;
 case 4:x=q.CalcPossibleDegrees(W,b[0]!,p);break;
 case 5:x=q.ConstTermTest(W,a,B[0]!,B[1]!,C,b[2]!,p);break;
 case 6:x=q.BalCopy(A,p);break;
 case 7:x=q.mul(W,p);break;
 case 8:x=q.mul(W,p,a);break;
 case 9:x=q.InvMul(W,a,p);break;
 case 10:x=q.RemoveFactors(W,a,p);break;
 case 11:x=q.unpack(A[0]!,b[0]!);break;
 case 12:x=q.SubPattern(a,b);break;
 default:throw new Error('unknown NTL selection operation');
 }
 return JSON.stringify(x,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}

const arg=(s:string):bigint|bigint[]=>s.startsWith('[')?s.slice(1,-1).trim()?s.slice(1,-1).split(',').map(v=>BigInt(v.trim())):[]:BigInt(s);
for(const[i,row]of fixtures.entries())test('native NTL factor selection '+i,()=>{
 let result:string|null=null,error:string|null=null,errorType:string|null=null;
 try{result=(ntl_factor_selection as Function)(...row.args.map(arg));}catch(e){error=(e as Error).message;errorType=(e as Error).name;}
 expect({result,error,errorType}).toEqual({result:row.result,error:row.error,errorType:row.errorType});
},10000);
test('NTL factor selection leaves source vectors and nested coefficients independent',()=>{
 const W=[[1n,1n],[2n,1n],[3n,1n]],copy=W.map(f=>f.slice());const kept=q.RemoveFactors(W,[1],101n);kept[0]![0]=99n;
 expect(W).toEqual(copy);expect(q.mul(W,101n)).toEqual([6n,11n,6n,1n]);expect(W).toEqual(copy);
 const pat=[1,2,3];const updated=q.SubPattern(pat,[0,1,2]);updated[0]=99;expect(pat).toEqual([1,2,3]);
});
test('NTL factor constant-term cache retains unused suffix and returns a copy',()=>{
 const prod=[7n,8n,9n];const out=q.ConstTermTest([[2n,1n],[3n,1n]],[0,1],12n,1n,prod,0,101n);
 expect(out).toEqual([1,[2n,6n,9n],1]);out[1][0]=99n;expect(prod).toEqual([7n,8n,9n]);
});
