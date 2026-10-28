import { afterEach, expect, spyOn, test } from 'bun:test';
import { spawnSync } from 'node:child_process';
import { createHash, randomInt } from 'node:crypto';
import {
  appendFileSync,
  existsSync,
  linkSync,
  mkdtempSync,
  readFileSync,
  rmSync,
  symlinkSync,
  writeFileSync,
} from 'node:fs';
import { tmpdir } from 'node:os';
import { join, relative } from 'node:path';
import { fileURLToPath } from 'node:url';
import fc from 'fast-check';
import { parse_args, run, runtime, validate_checkpoint } from './cli.js';
import { type Candidate, type Config, attempt_at, candidate, combination_at } from './core.js';
import { json, parse_json } from './json.js';
import { checkpoint_lock, lock_fd } from './storage.js';

const directories: string[] = [];
function temp() {
  const dir = mkdtempSync(join(tmpdir(), 'timestamp-ts-'));
  directories.push(dir);
  return { dir, out: join(dir, 'out'), cp: join(dir, 'cp'), full: join(dir, 'full') };
}
afterEach(() => {
  for (const dir of directories.splice(0)) rmSync(dir, { recursive: true, force: true });
});
const invoke = (...args: string[]) => run(parse_args(args)) as Record<string, unknown>;
const cli = fileURLToPath(new URL('./cli.ts', import.meta.url));
const child = (...args: string[]) =>
  spawnSync(process.execPath, [cli, ...args], { encoding: 'utf8' });
const generate = (out: string, cp: string, ...extra: string[]) =>
  invoke(
    'generate',
    '--start',
    '1',
    '--end',
    '4',
    '--output',
    out,
    '--checkpoint',
    cp,
    '--max-candidates',
    '1',
    '--progress-every',
    '0',
    ...extra
  );

test('estimate has exact ranges, deduplication, and stable attempt ordering', () => {
  const estimate = invoke(
    'estimate',
    '--start',
    '1',
    '--end',
    '6',
    '--step',
    '2',
    '--rng',
    'lrand48',
    'lrand48',
    '--curve',
    'p256',
    'secp256k1'
  );
  expect(estimate.attempts).toBe(6n);
  expect(estimate.timestamps).toBe(3n);
  expect(combination_at(estimate.config as Config, 5n)).toEqual([
    5n,
    'lrand48',
    'concat-be-reject',
    'secp256k1',
  ]);
  const extended = invoke(
    'estimate',
    '--start',
    '1',
    '--end',
    '3',
    '--rng',
    'lrand48',
    'msvc-rand',
    '--skip-draws',
    '0',
    '4',
    '4',
    '--keys-per-seed',
    '3'
  );
  expect(extended.attempts).toBe(24n);
  let index = 0n;
  for (const timestamp of [1n, 2n])
    for (const rng of ['lrand48', 'msvc-rand'] as const)
      for (const skip of [0, 4])
        for (let key = 0; key < 3; key++)
          expect(attempt_at(extended.config as Config, index++)).toEqual([
            timestamp,
            rng,
            'concat-be-reject',
            'p256',
            skip,
            key,
          ]);
  expect(invoke('estimate', '--start', '1', '--end', '18446744073709551658').attempts).toBe(
    18446744073709551657n
  );
});

test('invalid CLI ranges, consumption, options and budgets fail', () => {
  for (const args of [
    ['--start', '-1'],
    ['--step', '0'],
    ['--start', '5', '--end', '4'],
    ['--skip-draws', '-1'],
    ['--skip-draws', '1000001'],
    ['--keys-per-seed', '0'],
    ['--keys-per-seed', '1001'],
    ['--curve', 'unknown'],
    ['--rng', 'unknown'],
    ['--profile', 'python-getrandbits-reject'],
    ['--unit', 'hours'],
    ['--step', '1.5'],
    ['--bogus', '1'],
  ])
    expect(() => invoke('estimate', ...args)).toThrow();
  for (const value of ['nan', 'inf', '-1'])
    expect(() => invoke('generate', '--progress-every', value)).toThrow('Progress interval');
  expect(() => invoke('generate', '--max-candidates', '0')).toThrow('limits');
  expect(() => invoke('generate', '--checkpoint-every', '0')).toThrow('limits');
  expect(() => parse_args(['resume'])).toThrow('checkpoint');
  expect(() => parse_args(['estimate', '--start'])).toThrow('Missing');
  expect(child('--help').status).toBe(0);
});

test('resume verifies durable prefix, discards crash tail, and upgrades schema 1 without changing JSONL bytes', () => {
  const { out, cp, full } = temp();
  generate(out, cp);
  const saved = JSON.parse(readFileSync(cp, 'utf8'));
  saved.schema_version = 1;
  saved.config.skip_draws = undefined;
  saved.config.keys_per_seed = undefined;
  writeFileSync(cp, JSON.stringify(saved));
  appendFileSync(out, 'uncommitted crash tail\n');
  expect(invoke('resume', cp).complete).toBe(true);
  invoke('generate', '--start', '1', '--end', '4', '--output', full);
  expect(readFileSync(out)).toEqual(readFileSync(full));
  expect(JSON.parse(readFileSync(cp, 'utf8')).schema_version).toBe(2);
  const expected =
    [1n, 2n, 3n]
      .map((t) => json(candidate(t, 'seconds', 'lrand48', 'concat-be-reject', 'p256')))
      .join('\n') + '\n';
  expect(readFileSync(out, 'utf8')).toBe(expected);
});

test('corrupted/truncated prefixes and invalid checkpoints leave output untouched', () => {
  const { out, cp } = temp();
  generate(out, cp);
  const saved = JSON.parse(readFileSync(cp, 'utf8'));
  appendFileSync(out, 'uncommitted tail\n');
  const original = readFileSync(out);
  const bad: Record<string, unknown[]> = {
    emitted: ['bad', true, -1, 2],
    next_index: [true, -1, 4],
    output_offset: [true, -1, 9007199254740992],
    output_sha256: ['bad', null],
    output: ['relative', null],
    config: [null, []],
    schema_version: [true, 3],
  };
  for (const [key, values] of Object.entries(bad))
    for (const value of values) {
      writeFileSync(cp, JSON.stringify({ ...saved, [key]: value }));
      expect(() => invoke('resume', cp)).toThrow();
      expect(readFileSync(out)).toEqual(original);
    }
  for (const [key, values] of Object.entries({
    match_fingerprints: ['bad', {}, ['not-a-hash']],
    rngs: [[], ['lrand48', 'lrand48'], ['unknown']],
    start: [true, -1],
    skip_draws: [null, [], [0, 0]],
    keys_per_seed: [null, false, 0],
  }))
    for (const value of values) {
      writeFileSync(cp, JSON.stringify({ ...saved, config: { ...saved.config, [key]: value } }));
      expect(() => invoke('resume', cp)).toThrow();
      expect(readFileSync(out)).toEqual(original);
    }
  for (const value of [[], null, {}, { schema_version: true }]) {
    expect(() => validate_checkpoint(value)).toThrow();
  }
  writeFileSync(cp, JSON.stringify(saved));
  writeFileSync(out, Buffer.concat([Buffer.from('X'), original.subarray(1)]));
  expect(() => invoke('resume', cp)).toThrow('hash mismatch');
  expect(readFileSync(out)[0]).toBe(88);
  writeFileSync(out, 'short');
  expect(() => invoke('resume', cp)).toThrow('shorter');
  expect(readFileSync(out, 'utf8')).toBe('short');
});

test('existing files and checkpoint/output/lock aliases are refused', () => {
  const { dir, out, cp } = temp();
  writeFileSync(out, 'do not overwrite');
  expect(() => generate(out, cp)).toThrow();
  expect(readFileSync(out, 'utf8')).toBe('do not overwrite');
  expect(() => generate(cp, cp)).toThrow('must differ');
  expect(() => generate(`${cp}.lock`, cp)).toThrow('must differ');
  symlinkSync(`${cp}.lock`, join(dir, 'alias'));
  expect(() => generate(join(dir, 'alias'), cp)).toThrow('must differ');
  const fresh = temp();
  generate(fresh.out, fresh.cp);
  const saved = JSON.parse(readFileSync(fresh.cp, 'utf8'));
  linkSync(fresh.cp, join(fresh.dir, 'hardlink'));
  saved.output = join(fresh.dir, 'hardlink');
  writeFileSync(fresh.cp, JSON.stringify(saved));
  const before = readFileSync(fresh.cp);
  expect(() => invoke('resume', fresh.cp)).toThrow('must differ');
  expect(readFileSync(fresh.cp)).toEqual(before);
});

test('separate processes cannot claim a locked checkpoint before first save', () => {
  const { out, cp } = temp();
  const release = checkpoint_lock(cp);
  try {
    const competing = child(
      'generate',
      '--start',
      '1',
      '--end',
      '2',
      '--output',
      out,
      '--checkpoint',
      cp
    );
    expect(competing.status).toBe(2);
    expect(competing.stderr).toContain('Checkpoint is in use');
    expect(existsSync(out)).toBe(false);
  } finally {
    release();
  }
  expect(generate(out, cp).next_index).toBe(1n);
});

test('checkpoint changes between read and output lock never truncate', () => {
  const { out, cp } = temp();
  generate(out, cp);
  appendFileSync(out, 'crash tail');
  const original = readFileSync(out);
  const lock = spyOn(runtime, 'lock_fd').mockImplementation((fd) => {
    const data = JSON.parse(readFileSync(cp, 'utf8'));
    data.emitted++;
    writeFileSync(cp, JSON.stringify(data));
    lock_fd(fd);
  });
  try {
    expect(() => invoke('resume', cp)).toThrow('Checkpoint changed');
    expect(readFileSync(out)).toEqual(original);
  } finally {
    lock.mockRestore();
  }
});

test('interruption after write but before digest keeps last checkpoint', () => {
  const { out, cp, full } = temp();
  const digest = createHash('sha256');
  const update = digest.update.bind(digest);
  let calls = 0;
  digest.update = ((bytes: Buffer) => {
    if (++calls === 2) throw new Error('simulated interruption');
    return update(bytes);
  }) as typeof digest.update;
  const factory = spyOn(runtime, 'digest').mockReturnValue(digest);
  try {
    expect(() => generate(out, cp, '--max-candidates', '3', '--checkpoint-every', '1')).toThrow(
      'simulated interruption'
    );
  } finally {
    factory.mockRestore();
  }
  const saved = JSON.parse(readFileSync(cp, 'utf8'));
  expect(saved.next_index).toBe(1);
  expect(readFileSync(out).length).toBeGreaterThan(saved.output_offset);
  invoke('resume', cp);
  invoke('generate', '--start', '1', '--end', '4', '--output', full);
  expect(readFileSync(out)).toEqual(readFileSync(full));
});

test('matching resume needs no target file; empty matching sets still consume the cap', () => {
  const { dir, out, cp } = temp();
  const target = join(dir, 'match');
  writeFileSync(
    target,
    candidate(2n, 'seconds', 'lrand48', 'concat-be-reject', 'p256')!.spki_sha256
  );
  generate(out, cp, '--match', target);
  rmSync(target);
  expect(invoke('resume', cp).emitted_total).toBe(1n);
  expect(JSON.parse(readFileSync(out, 'utf8')).timestamp).toBe(2);
  const other = temp();
  writeFileSync(target, '# no targets');
  const result = generate(other.out, other.cp, '--match', target, '--max-candidates', '2');
  expect(result.attempted_this_run).toBe(2n);
  expect(result.emitted_total).toBe(0n);
  expect(result.nonmatching_this_run).toBe(2n);
});

test('progress counts rejected/nonmatching attempts and stays off stdout', () => {
  const { dir, out } = temp();
  const target = join(dir, 'match');
  writeFileSync(target, 'a'.repeat(64));
  const records = [
    null,
    { spki_sha256: 'b'.repeat(64) },
    { spki_sha256: 'a'.repeat(64) },
  ] as (Candidate | null)[];
  const generate = spyOn(runtime, 'candidate').mockImplementation(() => records.shift()!);
  let clock = 0;
  const time = spyOn(runtime, 'now').mockImplementation(() => clock++);
  const reports: Record<string, unknown>[] = [];
  const progress = spyOn(runtime, 'progress').mockImplementation((value) => {
    reports.push(value as Record<string, unknown>);
  });
  try {
    const result = invoke(
      'generate',
      '--start',
      '1',
      '--end',
      '10',
      '--max-candidates',
      '3',
      '--match',
      target,
      '--output',
      out,
      '--progress-every',
      '1'
    );
    expect(reports.map((r) => r.attempted_this_run)).toEqual([1n, 2n, 3n]);
    expect(reports.map((r) => r.remaining_this_run)).toEqual([2n, 1n, 0n]);
    expect(reports[0]!.eta_seconds_this_run).toBe(2);
    expect(result).toMatchObject({
      rejected_this_run: 1n,
      nonmatching_this_run: 1n,
      emitted_total: 1n,
      attempts_per_second: 0.75,
      complete: false,
    });
  } finally {
    generate.mockRestore();
    time.mockRestore();
    progress.mockRestore();
  }
  const generated = child('generate', '--start', '1', '--end', '3', '--progress-every', '0.000001');
  expect(generated.status).toBe(0);
  expect(
    generated.stdout
      .trim()
      .split('\n')
      .every((line) => JSON.parse(line).timestamp !== undefined)
  ).toBe(true);
  expect(generated.stderr).toContain('"event":"progress"');
  const quiet = child('generate', '--start', '1', '--end', '2', '--progress-every', '0');
  expect(quiet.status).toBe(0);
  expect(quiet.stderr.trim().split('\n')).toHaveLength(1);
});

test('fresh checkpoints resume mid-sequence with exact >64-bit timestamp JSON', () => {
  const replay = Number(process.env.WEAK_KEYS_SEED ?? randomInt(0, 0x7fffffff));
  console.info(`timestamp resume generator=v1 replay seed=${replay}`);
  fc.assert(
    fc.property(
      fc.bigInt({ min: 1n << 64n, max: (1n << 80n) - 1n }),
      fc.integer({ min: 1, max: 6 }),
      fc.integer({ min: 0, max: 9 }),
      (seed, cutoff, skip) => {
        const { out, cp, full } = temp();
        const options = [
          '--start',
          String(seed),
          '--end',
          String(seed + 1n),
          '--rng',
          'drand48-state',
          '--keys-per-seed',
          '4',
          '--skip-draws',
          '0',
          String(skip + 1),
          '--progress-every',
          '0',
        ];
        invoke(
          'generate',
          ...options,
          '--output',
          out,
          '--checkpoint',
          cp,
          '--max-candidates',
          String(cutoff)
        );
        appendFileSync(out, 'crash tail');
        expect(invoke('resume', cp).complete).toBe(true);
        invoke('generate', ...options, '--output', full);
        expect(readFileSync(out)).toEqual(readFileSync(full));
        const lines = readFileSync(out, 'utf8').trim().split('\n');
        for (const line of lines) expect((parse_json(line) as Candidate).timestamp).toBe(seed);
      }
    ),
    { seed: replay, numRuns: 8 }
  );
});

test('root launcher preserves relative output and checkpoint paths', () => {
  const { dir } = temp();
  const root = fileURLToPath(new URL('../../../../', import.meta.url));
  const result = spawnSync(
    process.execPath,
    ['run', '--cwd', root, 'weak-keys:timestampseeds', 'estimate', '--start', '1', '--end', '3'],
    { cwd: dir, encoding: 'utf8' }
  );
  expect(result.status).toBe(0);
  expect(JSON.parse(result.stdout).attempts).toBe(2);
  const generated = spawnSync(
    process.execPath,
    [
      '--silent',
      'run',
      'weak-keys:timestampseeds',
      'generate',
      '--start',
      '1',
      '--end',
      '2',
      '--output',
      relative(root, join(dir, 'relative')),
      '--checkpoint',
      relative(root, join(dir, 'checkpoint')),
    ],
    { cwd: root, encoding: 'utf8' }
  );
  expect(generated.status).toBe(0);
  expect(existsSync(join(dir, 'relative'))).toBe(true);
  expect(existsSync(join(dir, 'checkpoint'))).toBe(true);
});
