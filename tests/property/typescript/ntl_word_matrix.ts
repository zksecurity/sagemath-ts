import {mul} from '../../../packages/ntl-ts/src/mat_lzz_p.js';
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
