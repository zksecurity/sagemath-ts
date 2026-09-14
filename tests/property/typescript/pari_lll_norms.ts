import { ZM_lll_norms } from '../../../packages/parigp-ts/src/lll.js';
export function pari_lll_norms(
  op: bigint,
  flag: bigint,
  mode: bigint,
  m: bigint,
  n: bigint,
  p: bigint,
  shift: bigint,
  flat: bigint[]
): string {
  const A = Array.from({ length: Number(n) }, (_, j) =>
    Array.from({ length: Number(m) }, (_, i) => flat[i * Number(n) + j]!)
  );
  const r = ZM_lll_norms(A, op === 75n ? 0.75 : op === 999n ? 0.999 : 0.99, Number(flag));
  const out = (x: any): unknown =>
    typeof x === 'bigint'
      ? String(x)
      : Array.isArray(x)
        ? x.map(out)
        : x && typeof x === 'object'
          ? [x.s, String(x.e), String(x.m), x.p]
          : x;
  return JSON.stringify(out(r));
}

export function pari_lll_norms_resource(...args: Parameters<typeof pari_lll_norms>): string {
  try {
    pari_lll_norms(...args);
  } catch (e) {
    if (
      e instanceof RangeError &&
      e.message === 'fplll_dpe: Gram-Schmidt coefficient requires an unrepresentable shift'
    )
      return 'resource_failure';
    throw e;
  }
  return 'unexpected_success';
}
import { execFileSync } from 'node:child_process';
export function pari_lll_norms_bounded(
  op: bigint,
  flag: bigint,
  mode: bigint,
  m: bigint,
  n: bigint,
  p: bigint,
  shift: bigint,
  flat: bigint[]
): string {
  const A = Array.from({ length: Number(n) }, (_, j) =>
    Array.from({ length: Number(m) }, (_, i) => String(flat[i * Number(n) + j]!))
  );
  const program = `import {ZM_lll_norms}from ${JSON.stringify(new URL('../../../packages/parigp-ts/src/lll.js', import.meta.url).pathname)};
 const out=x=>typeof x==='bigint'?String(x):Array.isArray(x)?x.map(out):x&&typeof x==='object'?[x.s,String(x.e),String(x.m),x.p]:x;
 let value;try{value={kind:'result',value:out(ZM_lll_norms(${JSON.stringify(A)}.map(c=>c.map(BigInt)),${op === 75n ? 0.75 : op === 999n ? 0.999 : 0.99},${Number(flag)}))};}catch(e){value=e instanceof RangeError&&e.message==='fplll_dpe: Gram-Schmidt coefficient requires an unrepresentable shift'?{kind:'resource_failure'}:{kind:'error',type:e.name,message:e.message};}
 console.log(JSON.stringify(value));`;
  try {
    return execFileSync(process.execPath, ['-e', program], {
      encoding: 'utf8',
      timeout: 2000,
    }).trim();
  } catch (e) {
    if ((e as any).code === 'ETIMEDOUT') return JSON.stringify({ kind: 'timeout' });
    throw e;
  }
}
