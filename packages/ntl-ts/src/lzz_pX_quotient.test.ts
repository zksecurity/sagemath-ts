import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import {test,expect} from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./lzz_pX_quotient.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import * as q from './lzz_pX.js';
export function ntl_word_quotient(op:bigint,p:bigint,initialized:bigint,exponent:bigint,scalar:bigint,a:bigint[],b:bigint[],f:bigint[]):string{
 let x:unknown;
 if(op===4n)x=q.MulByXMod(a,f,p);
 else if(op===5n)x=q.InvMod(a,f,p);
 else if(op===6n)x=q.InvModStatus(a,f,p);
 else {const F=new q.zz_pXModulus(initialized?f:null,p);switch(Number(op)){
 case 0:x=[F.n,F.f];break;
 case 11:x=[F.modCrossover===45?1:F.modCrossover===90?2:3,F.modCrossover,Number(F.reciprocal.length>0)];break;
 case 10:{const before=q.rem(a,F);q.build(F,b);x=[before,q.rem(a,F),F.f,F.n];break;}
 case 1:x=q.rem(a,F);break;
 case 2:x=q.MulMod(a,b,F);break;
 case 3:x=q.SqrMod(a,F);break;
 case 7:x=q.PowerXMod(exponent,F);break;
 case 8:x=q.PowerXPlusAMod(scalar,exponent,F);break;
 case 9:x=q.PowerMod(a,exponent,F);break;
 default:throw new Error('unknown NTL word quotient operation');
 }}return JSON.stringify(x,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}

const arg=(s:string):bigint|bigint[]=>s.startsWith('[')?s.slice(1,-1).trim()?s.slice(1,-1).split(',').map(v=>BigInt(v.trim())):[]:BigInt(s);
for(const[i,row]of fixtures.entries())test('native NTL word polynomial quotient '+i,()=>{
 let result:string|null=null,error:string|null=null,errorType:string|null=null;
 try{result=(ntl_word_quotient as Function)(...row.args.map(arg));}catch(e){error=(e as Error).message;errorType=(e as Error).name;}
 expect({result,error,errorType}).toEqual({result:row.result,error:row.error,errorType:row.errorType});
},10000);
test('NTL word quotient construction and returned coefficients remain independent',()=>{
 const f=[1n,0n,1n],F=new q.zz_pXModulus(f,5n);f[0]=2n;expect(F.f).toEqual([1n,0n,1n]);
 const a=[1n,2n],r=q.rem(a,F);r[0]=99n;expect(a).toEqual([1n,2n]);expect(q.PowerXMod(3n,F)).toEqual([0n,4n]);
});
test('NTL plain modulus rebuild retains but does not use the previous reciprocal',()=>{
 const f=[1n,...Array<bigint>(96).fill(0n),1n],F=new q.zz_pXModulus(f,17n);expect(F.reciprocal.length).toBe(96);
 const previous=F.reciprocal.slice();q.build(F,[-1n,1n]);expect(F.reciprocal).toEqual(previous);expect(F.f).toEqual([16n,1n]);expect(q.rem([1n,2n,3n],F)).toEqual([6n]);
});
