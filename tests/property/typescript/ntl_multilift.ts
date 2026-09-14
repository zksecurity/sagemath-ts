import * as q from '../../../packages/ntl-ts/src/ZZXFactoring.js';
import * as z from '../../../packages/ntl-ts/src/ZZX1.js';
export function ntl_multilift(op:bigint,p:bigint,e:bigint,F:bigint[],A:bigint[],B:bigint[]):string {
 let x:unknown;
 if(op===0n)x=z.mul(A,B);
 else if(op===1n)x=z.sqr(A);
 else{const a:bigint[][]=[];for(let i=0,k=1;i<Number(A[0]);i++){const n=Number(A[k++]!);a.push(A.slice(k,k+n));k+=n}x=q.MultiLift(a,F,Number(e),p)}
 return JSON.stringify(x,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}
