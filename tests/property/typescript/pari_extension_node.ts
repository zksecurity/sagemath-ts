/** Node/V8 replay for native extension products beyond the argument-spread limit. */
import { spawnSync } from 'node:child_process';
import { mkdtempSync, writeFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
let executable: string | undefined;
export function nodeWordExtensionProduct(degree: number, zero: boolean): string {
  if (!executable) {
    const folder = mkdtempSync(join(tmpdir(), 'pari-extension-node-'));
    const entry = join(folder, 'entry.ts');
    const source = resolve(import.meta.dir, '../../../packages/parigp-ts/src/FlxX.ts');
    writeFileSync(
      entry,
      `import {FlxqX_mul} from ${JSON.stringify(source)};
import {createHash} from 'node:crypto';
const d=Number(process.argv[2]),zero=process.argv[3]==='1';
const T=[1n,...new Array(d-1).fill(0n),1n];
const value=FlxqX_mul([new Array(d).fill(1n)],zero?[]:[[1n]],T,7n);
const text=JSON.stringify(value,(_,v)=>typeof v==='bigint'?v.toString():v).replaceAll('"','');
process.stdout.write(createHash('sha256').update(text).digest('hex'));
`
    );
    const output = join(folder, 'entry.mjs');
    const build = spawnSync(
      process.execPath,
      ['build', entry, '--target=node', '--outfile=' + output],
      { encoding: 'utf8' }
    );
    if (build.status !== 0) {
      rmSync(folder, { recursive: true, force: true });
      throw new Error(build.stderr || 'Node extension fixture build failed');
    }
    executable = output;
    process.once('exit', () => rmSync(folder, { recursive: true, force: true }));
  }
  const result = spawnSync('node', [executable, String(degree), zero ? '1' : '0'], {
    encoding: 'utf8',
    timeout: 60000,
  });
  if (result.status !== 0) throw new Error(result.stderr || 'Node extension fixture failed');
  return result.stdout;
}
