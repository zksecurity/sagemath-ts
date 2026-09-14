import * as q from '../../../packages/ntl-ts/src/ZZXFactoring.js';
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
