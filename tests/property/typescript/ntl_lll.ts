import * as q from '../../../packages/ntl-ts/src/LLL.js';
import {_ntl_gexteucl} from '../../../packages/ntl-ts/src/lip.js';
export function ntl_lll(op:bigint,mm:bigint,nn:bigint,aa:bigint,bb:bigint,u:bigint,flat:bigint[],y:bigint[],initial:bigint[]):string {
 const m=Number(mm),n=Number(nn),a=Number(aa),b=Number(bb),M=Array.from({length:m},(_,i)=>flat.slice(i*n,(i+1)*n));
 const x=op===0n?q.LLL(M,a,b,Boolean(u)):op===1n?q.LLL_plus(M,a,b,Boolean(u)):op===2n?q.image(M,Boolean(u)):op===3n?q.LatticeSolve(M,y,a,initial,n):_ntl_gexteucl(y[0]!,y[1]!);
 return JSON.stringify(x,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}
