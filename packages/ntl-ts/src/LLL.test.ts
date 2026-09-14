import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import {test,expect} from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./LLL.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import * as q from './LLL.js';
import {_ntl_gexteucl} from './lip.js';
export function ntl_lll(op:bigint,mm:bigint,nn:bigint,aa:bigint,bb:bigint,u:bigint,flat:bigint[],y:bigint[],initial:bigint[]):string {
 const m=Number(mm),n=Number(nn),a=Number(aa),b=Number(bb),M=Array.from({length:m},(_,i)=>flat.slice(i*n,(i+1)*n));
 const x=op===0n?q.LLL(M,a,b,Boolean(u)):op===1n?q.LLL_plus(M,a,b,Boolean(u)):op===2n?q.image(M,Boolean(u)):op===3n?q.LatticeSolve(M,y,a,initial,n):_ntl_gexteucl(y[0]!,y[1]!);
 return JSON.stringify(x,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}

const arg=(s:string):bigint|bigint[]=>s.startsWith('[')?s.slice(1,-1).trim()?s.slice(1,-1).split(',').map(v=>BigInt(v.trim())):[]:BigInt(s);
for(const[i,row]of fixtures.entries())test('native NTL lattice and extended GCD '+i,()=>{
 let result:string|null=null,error:string|null=null,errorType:string|null=null;
 try{result=(ntl_lll as Function)(...row.args.map(arg));}catch(e){error=(e as Error).message;errorType=(e as Error).name;}
 expect({result,error,errorType}).toEqual({result:row.result,error:row.error,errorType:row.errorType});
},10000);
test('NTL lattice outputs own their basis, transformation and determinant vectors',()=>{
 const B=[[1n,1n,1n],[-1n,0n,2n],[3n,5n,6n]],copy=B.map(r=>r.slice());
 const [rank,D,M,U]=q.LLL_plus(B,3,4,true);expect(rank).toBe(3);expect(D).toEqual([1n,1n,2n,9n]);
 M[0]![0]=99n;U![0]![0]=99n;D[0]=99n;expect(B).toEqual(copy);
 expect(q.LLL_plus(B)[1]).toEqual([1n,1n,2n,9n]);expect(q.LLL_plus(B)[3]).toBeNull();
});
test('NTL image preserves the input and keeps dependent zero rows at the front',()=>{
 const B=[[2n,0n],[0n,3n],[1n,1n]],copy=B.map(r=>r.slice());const [rank,det,M,U]=q.image(B,true);
 expect(rank).toBe(2);expect(det).toBe(1n);expect(M).toEqual([[0n,0n],[1n,0n],[0n,1n]]);
 M[0]![0]=99n;U![0]![0]=99n;expect(B).toEqual(copy);
});
test('NTL lattice solve copies successful and unchanged unsuccessful outputs',()=>{
 const A=[[2n,0n],[0n,3n],[1n,1n]],y=[1n,2n],initial=[19n];const solved=q.LatticeSolve(A,y,0,initial);
 expect(solved).toEqual([1,[1n,1n,-1n]]);solved[1][0]=99n;expect(y).toEqual([1n,2n]);expect(A).toEqual([[2n,0n],[0n,3n],[1n,1n]]);
 const absent=q.LatticeSolve([[2n,0n],[0n,2n]],[1n,0n],0,initial);expect(absent).toEqual([0,[19n]]);absent[1][0]=99n;expect(initial).toEqual([19n]);
});
