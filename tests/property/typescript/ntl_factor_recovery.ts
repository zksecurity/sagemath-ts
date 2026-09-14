import * as q from '../../../packages/ntl-ts/src/ZZXFactoring.js';
export function ntl_factor_recovery(op:bigint,v:bigint[],F:bigint[],D:bigint[],flat:bigint[],W:bigint[]):string{
 let x:unknown;
 if(op===0n)x=q.PolyEval(F,v[0]!);
 else if(op===1n)x=q.RootBound(F);
 else if(op===2n){const m=Number(v[0]),n=Number(v[1]),C=Number(v[2]),r=Number(v[3]),d=Number(v[4]),M=Array.from({length:m},(_,i)=>flat.slice(i*n,(i+1)*n));x=q.CutAway(D,M,C!,r!,d!)}
 else{const w:bigint[][]=[];for(let i=0,k=1;i<Number(W[0]);i++){const n=Number(W[k++]!);w.push(W.slice(k,k+n));k+=n}x=q.AdditionalLifting(v[0]!,Number(v[1]),w,v[2]!,Number(v[3]),F,Boolean(v[4]))}
 return JSON.stringify(x,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}
