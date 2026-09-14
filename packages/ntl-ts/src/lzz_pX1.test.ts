import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import {test,expect} from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./lzz_pX1.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import {zz_pXModulus} from './lzz_pX.js';
import * as q from './lzz_pX1.js';
import {TraceMap,PowerCompose} from './lzz_pXFactoring.js';
export function ntl_word_composition(op:bigint,p:bigint,param:bigint[],packed:bigint[]):string{
 const w:bigint[][]=[];for(let i=0,k=1;i<Number(packed[0]);i++){const n=Number(packed[k++]!);w.push(packed.slice(k,k+n));k+=n;}
 let x:unknown;const m=Number(param[0]),d=Number(param[1]);
 if(op===8n)x=q.GCD(w[2]!,w[3]!,p);
 else{const F=new zz_pXModulus(param[2]?w[0]!:null,p),H=new q.zz_pXNewArgument();
 switch(Number(op)){
 case 0:x=q.CompMod(w[2]!,w[1]!,F);break;
 case 1:if(param[3])q.build(H,w[1]!,F,m);x=q.CompMod(w[2]!,H,F);break;
 case 2:q.build(H,w[1]!,F,m);x=[H.mat,H.poly];break;
 case 3:{if(param[3])q.build(H,w[1]!,F,m);const G=new zz_pXModulus(w[5]!,p);q.reduce(H,G);x=[[H.mat,H.poly],q.CompMod(w[2]!,H,G)];break;}
 case 4:x=q.Comp2Mod(w[2]!,w[3]!,w[1]!,F);break;
 case 5:x=q.Comp3Mod(w[2]!,w[3]!,w[4]!,w[1]!,F);break;
 case 6:x=TraceMap(w[2]!,d,F,w[1]!);break;
 case 7:x=PowerCompose(w[1]!,d,F);break;
 default:throw new Error('unknown NTL composition operation');
 }}return JSON.stringify(x,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}

const arg=(s:string):bigint|bigint[]=>s.startsWith('[')?s.slice(1,-1).trim()?s.slice(1,-1).split(',').map(v=>BigInt(v.trim())):[]:BigInt(s);
for(const[i,row]of fixtures.entries())test('native NTL word composition '+i,()=>{
 let result:string|null=null,error:string|null=null,errorType:string|null=null;
 try{result=(ntl_word_composition as Function)(...row.args.map(arg));}catch(e){error=(e as Error).message;errorType=(e as Error).name;}
 expect({result,error,errorType}).toEqual({result:row.result,error:row.error,errorType:row.errorType});
},10000);
test('NTL composition cache and outputs preserve input polynomials',()=>{
 const F=new zz_pXModulus([1n,0n,1n],5n),h=[1n,1n],g=[1n,2n,3n],H=new q.zz_pXNewArgument();q.build(H,h,F,2);
 const x=q.CompMod(g,H,F);expect(x).toEqual([3n,3n]);x[0]=99n;expect(g).toEqual([1n,2n,3n]);expect(h).toEqual([1n,1n]);expect(H.mat).toEqual([[1n,0n],[1n,1n]]);
});
test('NTL composition reduction replaces the old cache independently',()=>{
 const F=new zz_pXModulus([1n,0n,1n],5n),H=new q.zz_pXNewArgument();q.build(H,[1n,1n],F,2);const old=H.mat;
 const G=new zz_pXModulus([-1n,1n],5n);q.reduce(H,G);old[0]![0]=99n;expect(H.mat).toEqual([[1n],[2n]]);expect(H.poly).toEqual([4n]);expect(q.CompMod([1n,2n,3n],H,G)).toEqual([2n]);
});
test('NTL word GCD returns an independent normalized polynomial',()=>{
 const a=[-2n,0n,2n],b=[-3n,3n],g=q.GCD(a,b,5n);expect(g).toEqual([4n,1n]);g[0]=99n;expect(a).toEqual([-2n,0n,2n]);expect(b).toEqual([-3n,3n]);
});
