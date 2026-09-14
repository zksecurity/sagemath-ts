import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { ZZ_pXModulus, build, rem } from './ZZ_pX.js';
import { FFTPrimeContext } from './FFT.js';
import { RandomStream } from './ZZ.js';
import { test, expect } from 'bun:test';

import { fileURLToPath } from 'node:url';

const fixtures = (await loadLiveNative(import.meta.url, "./ZZ_pX_zero_modulus.native.json.gz")) as { args: string[]; result: string | null; error: string | null; errorType: string | null }[];

test('native NTL arbitrary-modulus constant-modulus reduction boundaries terminate with exact results and state', async () => {
  const child = Bun.spawn(
    [process.execPath, fileURLToPath(new URL('./ZZ_pX_zero_modulus.runner.mjs', import.meta.url))],
    { stdin: 'pipe', stdout: 'pipe', stderr: 'pipe' }
  );
  let timedOut = false;
  const timer = setTimeout(() => {
    timedOut = true;
    child.kill('SIGKILL');
  }, 5000);
  try {
    child.stdin.write(JSON.stringify(fixtures.map((row) => row.args)));
    child.stdin.end();
    const [exit, output, errors] = await Promise.all([
      child.exited,
      new Response(child.stdout).text(),
      new Response(child.stderr).text(),
    ]);
    expect(timedOut, 'constant-modulus reduction child exceeded its five-second watchdog').toBe(
      false
    );
    expect(exit, errors).toBe(0);
    const actual = JSON.parse(output);
    expect(actual.length).toBe(fixtures.length);
    for (let i = 0; i < fixtures.length; i++) {
      const { result, error, errorType } = fixtures[i]!;
      expect(actual[i], 'native constant-modulus trace ' + i).toEqual({ result, error, errorType });
    }
    // Only replay in this process after the identical native cohort terminates
    // in the watched child. This also instruments the regression's error branch.
    const state = {
      context: new FFTPrimeContext(),
      stream: new RandomStream(Uint8Array.from({ length: 32 }, (_, i) => i)),
    };
    const F = new ZZ_pXModulus(null, 6n, state);
    build(F, [1n, ...Array<bigint>(22).fill(0n), 1n]);
    expect(() => build(F, [1n])).toThrow('build: deg(f) must be at least 1');
    expect(() => rem([...Array<bigint>(21).fill(0n), 1n, 0n, 6n, 0n], F)).toThrow(
      'negative length in vector::SetLength'
    );
  } finally {
    clearTimeout(timer);
    if (child.exitCode === null) child.kill('SIGKILL');
    await child.exited;
  }
}, 10000);
