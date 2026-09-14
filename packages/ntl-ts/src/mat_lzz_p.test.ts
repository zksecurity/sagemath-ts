import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import {test,expect} from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./mat_lzz_p.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import {mul} from './mat_lzz_p.js';
export function ntl_word_matrix(op:bigint,p:bigint,dims:bigint[],seed:bigint,a:bigint[],b:bigint[]):string{
 if(p<=1n||p>=1n<<60n)mul([],[],p);
 const m=Number(dims[0]),k=Number(dims[1]),k2=Number(dims[2]),n=Number(dims[3]);
 if(m<0||k2<0)throw new Error('SetDims: bad args');
 let state=seed;const next=()=>{state=(state*6364136223846793005n+1442695040888963407n)&((1n<<64n)-1n);const v=((state^(state>>32n))>>1n)%p;return state&1n?-v:v;};
 const A=Array.from({length:m},(_,i)=>Array.from({length:k},(_,j)=>op?next():a[i*k+j]!));
 const B=Array.from({length:k2},(_,i)=>Array.from({length:n},(_,j)=>op?next():b[i*n+j]!));
 const C=mul(A,B,p,n,k);
 const x=op?[String(m),String(n),C.map(row=>row.map(v=>v.toString(16).padStart(Math.ceil((p-1n).toString(2).length/4),'0')).join('')).join('')]:C;
 return JSON.stringify(x,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}

const arg=(s:string):bigint|bigint[]=>s.startsWith('[')?s.slice(1,-1).trim()?s.slice(1,-1).split(',').map(v=>BigInt(v.trim())):[]:BigInt(s);
for(const[i,row]of fixtures.entries())test('native NTL word matrix '+i,()=>{
 let result:string|null=null,error:string|null=null,errorType:string|null=null;
 try{result=(ntl_word_matrix as Function)(...row.args.map(arg));}catch(e){error=(e as Error).message;errorType=(e as Error).name;}
 expect({result,error,errorType}).toEqual({result:row.result,error:row.error,errorType:row.errorType});
},10000);
test('NTL word matrix multiplication preserves both input matrices',()=>{
 const A=[[1n,2n],[3n,4n]],B=[[5n,6n],[7n,8n]],C=mul(A,B,101n);expect(C).toEqual([[19n,22n],[43n,50n]]);C[0]![0]=99n;expect(A).toEqual([[1n,2n],[3n,4n]]);expect(B).toEqual([[5n,6n],[7n,8n]]);
});
test('NTL word matrix multiplication handles shared operands without mutation',()=>{
 const A=[[1n,2n],[3n,4n]];expect(mul(A,A,5n)).toEqual([[2n,0n],[0n,2n]]);expect(A).toEqual([[1n,2n],[3n,4n]]);
});
