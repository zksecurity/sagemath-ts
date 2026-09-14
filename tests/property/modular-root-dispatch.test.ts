import { nativeFixtures as loadLiveNative } from "./native-live.mjs";
import { expect, test } from 'bun:test';
import { execFileSync } from 'node:child_process';
const native = await loadLiveNative(import.meta.url, "./modular-root-dispatch.native.json");

// Recorded from the live Sage runner with the bundled residue-root methods.
// Exercise the real registry: a helper must never overwrite IntegerMod.lift.
test('integer lift and residue-root lifting retain separate runner dispatch', () => {
  const suite = {
    module: 'modular_integers',
    cases: native.map((row) => ({
      function: row.function,
      rows: [[row.seed, ...row.args.map((arg) => (arg.startsWith('[') ? JSON.parse(arg) : arg))]],
    })),
  };
  const output = execFileSync(process.execPath, [`${import.meta.dir}/typescript/runner.ts`], {
    input: JSON.stringify(suite),
    encoding: 'utf8',
    timeout: 10_000,
  });
  expect(JSON.parse(output)).toEqual(native);
}, 15_000);
