import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import {test,expect} from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./ZZXFactoring_traces.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import * as q from './ZZXFactoring.js';
export function ntl_factor_traces(op:bigint,v:bigint[],f:bigint[],Tr:bigint[],C:bigint[],pb:bigint[],A:bigint[],B:bigint[]):string {
 const p=v[0]!,[d,d1,n]=v.slice(1,4).map(Number),root=v[4]!,bit_delta=Number(v[5]),lc=v[6]!,P=v[7]!,pd=v[8]!,delta=Number(v[9]);
 const unpack=(a:bigint[],columns:number)=>Array.from({length:Number(a[0])},(_,i)=>a.slice(1+i*columns,1+(i+1)*columns));
 let x:unknown;
 if(op===0n)x=q.ComputeTrace(Tr,f,d!,P);
 else if(op===1n)x=q.ChopTraces(C,Tr,d!,pb,pd,P,lc);
 else if(op===2n)x=q.DenseChopTraces(C,Tr,d!,d1!,root,pd,P,lc,unpack(A,d!));
 else if(op===3n)x=q.Compute_pb(C.map(Number),pb,p!,d!,root,n!);
 else if(op===4n)x=q.Compute_pdelta(delta,pd,p!,bit_delta);
 else if(op===5n)x=q.BuildReductionMatrix(n!,d!,pd,unpack(A,d!),unpack(B,n!));
 else if(op===6n)x=q.Compute_pb_eff(p!,d!,root,n!,bit_delta);
 else if(op===7n)x=q.d1_val(bit_delta,n!,d1!);
 return JSON.stringify(x,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}

const arg=(s:string):bigint|bigint[]=>s.startsWith('[')?s.slice(1,-1).trim()?s.slice(1,-1).split(',').map(v=>BigInt(v.trim())):[]:BigInt(s);
for(const[i,row]of fixtures.entries())test('native NTL factor traces '+i,()=>{
 let result:string|null=null,error:string|null=null,errorType:string|null=null;
 try{result=(ntl_factor_traces as Function)(...row.args.map(arg));}catch(e){error=(e as Error).message;errorType=(e as Error).name;}
 expect({result,error,errorType}).toEqual({result:row.result,error:row.error,errorType:row.errorType});
});
test('NTL trace updates copy input vectors and retain unused tails',()=>{
 const traces=[0n,0n,91n],copy=traces.slice(),f=[-2n,0n,1n];
 const out=q.ComputeTrace(traces,f,2,125n);expect(out).toEqual([0n,4n,91n]);out[0]=99n;expect(traces).toEqual(copy);
 const C=[0n,0n,99n],chopped=q.ChopTraces(C,[4n,8n],2,[5n,5n],5n,125n,2n);
 expect(chopped).toEqual([2n,1n,99n]);chopped[2]=0n;expect(C).toEqual([0n,0n,99n]);
 const b=[3,99],pb=[27n,99n],next=q.Compute_pb(b,pb,3,2,2n,4);
 expect(next).toEqual([[3,4],[27n,81n]]);expect(b).toEqual([3,99]);expect(pb).toEqual([27n,99n]);
});
test('NTL factor lattice results do not alias either input matrix',()=>{
 const A=[[1n,2n],[3n,4n]],B=[[1n,0n],[1n,-1n]],result=q.BuildReductionMatrix(2,2,9n,A,B);
 expect(result).toEqual([[[2n,0n,1n,2n],[2n,-2n,-2n,-2n],[0n,0n,9n,0n],[0n,0n,0n,9n]],2]);
 result[0][0]![0]=999n;expect(A).toEqual([[1n,2n],[3n,4n]]);expect(B).toEqual([[1n,0n],[1n,-1n]]);
});
