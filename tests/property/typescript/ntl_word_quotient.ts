import * as q from '../../../packages/ntl-ts/src/lzz_pX.js';
export function ntl_word_quotient(op:bigint,p:bigint,initialized:bigint,exponent:bigint,scalar:bigint,a:bigint[],b:bigint[],f:bigint[]):string{
 let x:unknown;
 if(op===4n)x=q.MulByXMod(a,f,p);
 else if(op===5n)x=q.InvMod(a,f,p);
 else if(op===6n)x=q.InvModStatus(a,f,p);
 else {const F=new q.zz_pXModulus(initialized?f:null,p);switch(Number(op)){
 case 0:x=[F.n,F.f];break;
 case 11:x=[F.modCrossover===45?1:F.modCrossover===90?2:3,F.modCrossover,Number(F.reciprocal.length>0)];break;
 case 10:{const before=q.rem(a,F);q.build(F,b);x=[before,q.rem(a,F),F.f,F.n];break;}
 case 1:x=q.rem(a,F);break;
 case 2:x=q.MulMod(a,b,F);break;
 case 3:x=q.SqrMod(a,F);break;
 case 7:x=q.PowerXMod(exponent,F);break;
 case 8:x=q.PowerXPlusAMod(scalar,exponent,F);break;
 case 9:x=q.PowerMod(a,exponent,F);break;
 default:throw new Error('unknown NTL word quotient operation');
 }}return JSON.stringify(x,(_,v)=>typeof v==='bigint'||typeof v==='number'?String(v):v);
}
