import { ZX_divides_i, ZX_divides, cmbf_precs, cmbf } from '../../../packages/parigp-ts/src/QX_factor.js';
import { FpXV_prod } from '../../../packages/parigp-ts/src/FpX.js';
import { centermodii } from '../../../packages/parigp-ts/src/polarit2.js';
import { cmbf_maxK } from '../../../packages/parigp-ts/src/nffactor.js';
export function pari_recombination(op: bigint, p: bigint, a: bigint, b: bigint, k: bigint, B: bigint, x: bigint[], f: bigint[]): string {
  const factors=()=>{const out:bigint[][]=[];for(let i=0;i<f.length;){const n=Number(f[i++]!);out.push(f.slice(i,i+n));i+=n;}return out;};
  let result:unknown;
  if(op===0n)result=ZX_divides_i(x,f,B<0n?null:B);
  else if(op===6n)result=ZX_divides(x,f);
  else if(op===1n)result=cmbf_precs(p,x[0]!,B);
  else if(op===2n){const r=cmbf(x,factors(),B,p,Number(a),Number(b),Number(k));result=[r[0],r[1],r[2],r[3]?1:0];}
  else if(op===3n)result=centermodii(x[0]!,p,a?B:null);
  else if(op===4n)result=FpXV_prod(factors(),p);
  else result=cmbf_maxK(Number(a));
  return JSON.stringify(result,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}
