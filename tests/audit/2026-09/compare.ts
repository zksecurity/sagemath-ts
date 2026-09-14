/** Run with `bun tests/audit/2026-09/compare.ts`; exits 1 while discrepancies remain. */
import { join } from 'node:path';

async function run(command: string[]) {
  const process = Bun.spawn(command, { stdout: 'pipe', stderr: 'pipe' });
  const [stdout, stderr, status] = await Promise.all([
    new Response(process.stdout).text(),
    new Response(process.stderr).text(),
    process.exited,
  ]);
  if (status !== 0) throw new Error(`${command[0]} exited ${status}: ${stderr}`);
  return JSON.parse(stdout);
}
const [sage, typescript] = await Promise.all([
  run(['sage', join(import.meta.dir, 'python.py')]),
  run(['bun', join(import.meta.dir, 'typescript.ts')]),
]);
if (sage.length !== typescript.length) throw new Error('Different result counts');
let mismatches = 0;
for (let i = 0; i < sage.length; i++) {
  const expected = sage[i];
  const actual = typescript[i];
  if (expected.id !== actual.id) throw new Error('Different case order');
  // Messages are preserved in the output; classify by value or error type so
  // equivalent exceptions with backend-specific wording do not create findings.
  if (
    expected.error === actual.error &&
    JSON.stringify(expected.value) === JSON.stringify(actual.value)
  )
    continue;
  mismatches++;
  console.log(JSON.stringify({ id: expected.id, sage: expected, typescript: actual }));
}
console.log(`${sage.length} cases: ${sage.length - mismatches} match, ${mismatches} discrepancies`);
process.exitCode = mismatches ? 1 : 0;
