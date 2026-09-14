import * as q from '../../../packages/ntl-ts/src/ZZXFactoring.js';
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
