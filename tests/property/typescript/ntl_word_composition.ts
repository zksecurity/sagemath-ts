import {zz_pXModulus} from '../../../packages/ntl-ts/src/lzz_pX.js';
import * as q from '../../../packages/ntl-ts/src/lzz_pX1.js';
import {TraceMap,PowerCompose} from '../../../packages/ntl-ts/src/lzz_pXFactoring.js';
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
