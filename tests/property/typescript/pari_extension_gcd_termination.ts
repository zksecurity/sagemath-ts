/** Isolate the bundled binary GCD cycle; compare termination with native controls. */
import { spawnSync } from 'node:child_process';
import { resolve } from 'node:path';
export function binaryExtensionGcdTermination(degree: number, op: number): string {
  const source = resolve(import.meta.dir, '../../../packages/parigp-ts/src/_extension_gcd.ts');
  const script = `import {extensionGcd} from ${JSON.stringify(source)};
const a = [1n,...new Array(${degree}-1).fill(0n),1n];
const result=extensionGcd(2,${op},2n,7n,a,[0n,0n,1n]);
console.log('OK '+JSON.stringify(result,(_,v)=>typeof v==='bigint'?v.toString():v).replaceAll('"',''));`;
  const result = spawnSync(process.execPath, ['-e', script], {
    encoding: 'utf8',
    timeout: 2000,
  });
  if (result.error && 'code' in result.error && result.error.code === 'ETIMEDOUT') return 'timeout';
  if (result.status !== 0) throw new Error(result.stderr || 'isolated binary GCD probe failed');
  return result.stdout.trim();
}
