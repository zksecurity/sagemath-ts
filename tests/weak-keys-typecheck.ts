/**
 * The existing cheatsheets dependency graph is not tsc-clean. Compare diagnostics live,
 * rather than suppressing dependency checking or checking in a stale diagnostic snapshot.
 * This is a scoped regression gate, NOT a claim that full repository typechecking passes.
 */
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const root = fileURLToPath(new URL('../', import.meta.url));
function diagnostics(project: string): Set<string> {
  const run = spawnSync(
    process.execPath,
    ['x', 'tsc', '--noEmit', '--pretty', 'false', '-p', project],
    { cwd: root, encoding: 'utf8', maxBuffer: 16 * 1024 * 1024 }
  );
  const lines = `${run.stdout ?? ''}\n${run.stderr ?? ''}`.split('\n');
  const errors = new Set(lines.filter((line) => /^.+\(\d+,\d+\): error TS\d+:/.test(line)));
  if (lines.some((line) => /^error TS\d+:/.test(line)))
    throw new Error('Global compiler error; dependency baseline is not usable');
  if (run.error || (run.status !== 0 && run.status !== 2) || (run.status !== 0 && !errors.size))
    throw new Error('Typecheck did not produce a usable result');
  return errors;
}

const baseline = diagnostics('packages/zksecurity-cheatsheets/tsconfig.json');
const current = diagnostics('packages/weak-keys/tsconfig.json');
const introduced = [...current].filter(
  (line) => line.startsWith('packages/weak-keys/') || !baseline.has(line)
);
console.log(
  `Scoped weak-key typecheck: ${current.size} dependency diagnostics; ${introduced.length} introduced (${baseline.size} in existing cheatsheets graph).`
);
if (introduced.length) {
  console.error(introduced.slice(0, 10).join('\n'));
  process.exitCode = 1;
}
