import { ZM_hnf, ZM_hnfall_i, ZM_hnfall, hnfall } from '../../../packages/parigp-ts/src/hnf_snf.js';
export function pari_general_hnf(op:bigint,u:bigint,remove:bigint,m:bigint,n:bigint,flat:bigint[]):string {
 const A=Array.from({length:Number(n)},(_,j)=>Array.from({length:Number(m)},(_,i)=>flat[i*Number(n)+j]!));
 const r=op===0n?ZM_hnf(A):op===4n?hnfall(A):op===2n?ZM_hnfall_i(A,!!u,Number(remove)as 0|1|2):ZM_hnfall(A,!!u,Number(remove)as 0|1|2);
 return JSON.stringify(r,(_,v)=>typeof v==='bigint'?String(v):v);
}
