import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';
const fixtures = await loadLiveNative(import.meta.url, "./lll_norms.bounded.native.json");
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
  const program = `import {ZM_lll_norms}from ${JSON.stringify(new URL('./lll.js', import.meta.url).pathname)};
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

const functions = { pari_lll_norms_bounded };

const arg = (s: string): bigint | bigint[] =>
  s.startsWith('[')
    ? s.slice(1, -1).trim()
      ? s
          .slice(1, -1)
          .split(',')
          .map((v) => BigInt(v.trim()))
      : []
    : BigInt(s);
for (const [i, row] of fixtures.entries())
  test('native adaptive LLL norms ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = (functions[row.function as keyof typeof functions] as any)(...row.args.map(arg));
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
  }, 5000);
