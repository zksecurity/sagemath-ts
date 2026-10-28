#!/usr/bin/env bun
import { createHash } from 'node:crypto';
import {
  closeSync,
  existsSync,
  fstatSync,
  fsyncSync,
  ftruncateSync,
  openSync,
  readFileSync,
  readSync,
  realpathSync,
  statSync,
} from 'node:fs';
import { basename, dirname, isAbsolute, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  CURVES,
  type Config,
  type Curve,
  MAX_KEYS_PER_SEED,
  MAX_SKIP_DRAWS,
  PROFILES,
  type Profile,
  type Unit,
  attempt_at,
  candidate,
  parse_time,
  stream_at,
  timestamp_count,
  total_candidates,
  validate_profile,
} from './core.js';
import { json, parse_json } from './json.js';
import { RNGS, type Rng } from './rngs.js';
import { atomic_json, checkpoint_lock, lock_fd, write_all } from './storage.js';
import { targets, valid_fingerprint } from './targets.js';

export const SCHEMA = 2;
export interface Args {
  command: 'estimate' | 'generate' | 'resume';
  start: string;
  end: string;
  unit: Unit;
  step: bigint;
  rngs: Rng[];
  profiles: Profile[];
  curves: Curve[];
  match?: string[];
  skip_draws: number[];
  keys_per_seed: number;
  max_candidates: bigint;
  output?: string;
  checkpoint?: string;
  checkpoint_every: bigint;
  progress_every: number;
}
const HELP = `Offline public-key candidates for explicit timestamp/RNG/derivation hypotheses.
Usage: timestampseeds estimate|generate [options]
       timestampseeds resume CHECKPOINT [--max-candidates N] [--checkpoint-every N] [--progress-every SECONDS]
--start INTEGER_OR_ISO --end INTEGER_OR_ISO (exclusive) --unit seconds|milliseconds --step N
--rng ${RNGS.join('|')}
--profile ${PROFILES.join('|')}
--curve p256|secp256k1 --skip-draws N [N ...] --keys-per-seed N
--match FILE (repeatable) --output FILE --checkpoint FILE
--max-candidates N (default 10000) --checkpoint-every N (default 100)
--progress-every SECONDS (default 5; 0 disables). Never outputs private keys.`;
function integer(text: string): bigint {
  if (!/^[+-]?\d(?:_?\d)*$/.test(text)) throw new Error(`Expected integer: ${text}`);
  return BigInt(text.replaceAll('_', ''));
}
function small(text: string): number {
  const result = Number(integer(text));
  if (!Number.isSafeInteger(result)) throw new Error('Integer exceeds safe count range');
  return result;
}
export function parse_args(argv: string[]): Args {
  const command = argv[0];
  if (!['estimate', 'generate', 'resume'].includes(command!)) throw new Error(HELP);
  const args: Args = {
    command: command as Args['command'],
    start: '2000-01-01T00:00:00Z',
    end: new Date(Math.floor(Date.now() / 1000) * 1000).toISOString(),
    unit: 'seconds',
    step: 1n,
    rngs: ['lrand48'],
    profiles: ['concat-be-reject'],
    curves: ['p256'],
    skip_draws: [0],
    keys_per_seed: 1,
    max_candidates: 10000n,
    checkpoint_every: 100n,
    progress_every: 5,
  };
  let i = 1;
  if (command === 'resume') {
    if (!argv[i] || argv[i]!.startsWith('--')) throw new Error('resume requires a checkpoint');
    args.checkpoint = argv[i++];
  }
  const common = [
    '--start',
    '--end',
    '--unit',
    '--step',
    '--rng',
    '--profile',
    '--curve',
    '--skip-draws',
    '--keys-per-seed',
    '--match',
  ];
  const generation = ['--max-candidates', '--checkpoint-every', '--progress-every'];
  const allowed =
    command === 'resume'
      ? generation
      : command === 'estimate'
        ? common
        : [...common, ...generation, '--output', '--checkpoint'];
  while (i < argv.length) {
    const token = argv[i++]!;
    const equals = token.indexOf('=');
    const option = equals < 0 ? token : token.slice(0, equals);
    if (!allowed.includes(option)) throw new Error(`Unknown option: ${option}`);
    const values: string[] = equals < 0 ? [] : [token.slice(equals + 1)];
    const multiple = ['--rng', '--profile', '--curve', '--skip-draws'].includes(option);
    if (equals < 0) {
      do {
        if (!argv[i] || argv[i]!.startsWith('--')) break;
        values.push(argv[i++]!);
      } while (multiple);
    }
    if (!values.length) throw new Error(`Missing value for ${option}`);
    const value = values[0]!;
    switch (option) {
      case '--start':
        args.start = value;
        break;
      case '--end':
        args.end = value;
        break;
      case '--unit':
        if (!['seconds', 'milliseconds'].includes(value)) throw new Error('Unknown timestamp unit');
        args.unit = value as Unit;
        break;
      case '--step':
        args.step = integer(value);
        break;
      case '--rng':
        args.rngs = values as Rng[];
        break;
      case '--profile':
        args.profiles = values as Profile[];
        break;
      case '--curve':
        args.curves = values as Curve[];
        break;
      case '--skip-draws':
        args.skip_draws = values.map(small);
        break;
      case '--keys-per-seed':
        args.keys_per_seed = small(value);
        break;
      case '--match':
        args.match ??= [];
        args.match.push(value);
        break;
      case '--output':
        args.output = value;
        break;
      case '--checkpoint':
        args.checkpoint = value;
        break;
      case '--max-candidates':
        args.max_candidates = integer(value);
        break;
      case '--checkpoint-every':
        args.checkpoint_every = integer(value);
        break;
      case '--progress-every':
        args.progress_every = Number(value);
        break;
    }
  }
  return args;
}
export function validate_consumption(config: Pick<Config, 'skip_draws' | 'keys_per_seed'>): void {
  const skips = config.skip_draws ?? [0];
  const keys = config.keys_per_seed ?? 1;
  if (
    !Array.isArray(skips) ||
    !skips.length ||
    !skips.every((n) => Number.isSafeInteger(n) && n >= 0 && n <= MAX_SKIP_DRAWS)
  )
    throw new Error(`Skipped draws must be between 0 and ${MAX_SKIP_DRAWS}`);
  if (new Set(skips).size !== skips.length) throw new Error('Duplicate skipped-draw counts');
  if (!Number.isSafeInteger(keys) || keys < 1 || keys > MAX_KEYS_PER_SEED)
    throw new Error(`Keys per seed must be between 1 and ${MAX_KEYS_PER_SEED}`);
}
export function make_config(args: Args): Config {
  const config: Config = {
    start: parse_time(args.start, args.unit),
    end: parse_time(args.end, args.unit),
    unit: args.unit,
    step: args.step,
    rngs: [...new Set(args.rngs)],
    profiles: [...new Set(args.profiles)],
    curves: [...new Set(args.curves)],
    skip_draws: [...new Set(args.skip_draws)],
    keys_per_seed: args.keys_per_seed,
    match_fingerprints: null,
  };
  if (config.start < 0n || config.end <= config.start || config.step <= 0n)
    throw new Error('Require 0 <= start < end and step > 0');
  if (!config.curves.length || !config.curves.every((curve) => Object.hasOwn(CURVES, curve)))
    throw new Error('Unknown curve');
  for (const rng of config.rngs)
    for (const profile of config.profiles) validate_profile(rng, profile);
  validate_consumption(config);
  const selected = targets(args.match);
  config.match_fingerprints = selected === null ? null : [...selected].sort();
  return config;
}
interface Checkpoint {
  schema_version: number;
  config: Config;
  next_index: bigint;
  emitted: bigint;
  output: string;
  output_offset: number;
  output_sha256: string;
}
function object(value: unknown): Record<string, unknown> {
  if (value === null || typeof value !== 'object' || Array.isArray(value))
    throw new Error('Incomplete checkpoint');
  return value as Record<string, unknown>;
}
function exact(value: unknown, message: string): bigint {
  if (typeof value === 'bigint') return value;
  if (typeof value === 'number' && Number.isSafeInteger(value)) return BigInt(value);
  throw new Error(message);
}
/** Validate completely before opening/truncating output; normalize integer transport types.
 * @see Deviation: Weak-key application package (no SageMath counterpart)
 */
export function validate_checkpoint(input: unknown): Checkpoint {
  const data = object(input);
  const c = object(data.config);
  if (data.schema_version !== 1 && data.schema_version !== SCHEMA)
    throw new Error('Unsupported checkpoint schema');
  if (
    data.schema_version === SCHEMA &&
    (!Object.hasOwn(c, 'skip_draws') || !Object.hasOwn(c, 'keys_per_seed'))
  )
    throw new Error('Missing checkpoint consumption settings');
  if (
    data.schema_version === 1 &&
    (json(c.skip_draws ?? [0]) !== '[0]' || (c.keys_per_seed ?? 1) !== 1)
  )
    throw new Error('Schema 1 cannot contain extended consumption settings');
  // Explicit null is invalid; only absent schema-1 settings receive defaults.
  if (c.skip_draws === null || c.keys_per_seed === null)
    throw new Error('Invalid checkpoint consumption settings');
  validate_consumption(c as Pick<Config, 'skip_draws' | 'keys_per_seed'>);
  const start = exact(c.start, 'Invalid checkpoint range');
  const end = exact(c.end, 'Invalid checkpoint range');
  const step = exact(c.step, 'Invalid checkpoint range');
  if (
    !['seconds', 'milliseconds'].includes(c.unit as string) ||
    start < 0n ||
    end <= start ||
    step <= 0n
  )
    throw new Error('Invalid checkpoint range');
  for (const [key, allowed] of [
    ['rngs', RNGS],
    ['profiles', PROFILES],
    ['curves', Object.keys(CURVES)],
  ] as const) {
    const choices = c[key];
    if (
      !Array.isArray(choices) ||
      !choices.length ||
      !choices.every((x) => typeof x === 'string' && (allowed as readonly string[]).includes(x))
    )
      throw new Error('Unknown checkpoint algorithm');
    if (new Set(choices).size !== choices.length) throw new Error('Duplicate checkpoint algorithm');
  }
  if (
    c.match_fingerprints !== null &&
    (!Array.isArray(c.match_fingerprints) || !c.match_fingerprints.every(valid_fingerprint))
  )
    throw new Error('Invalid checkpoint fingerprints');
  const config = { ...c, start, end, step } as unknown as Config;
  for (const rng of config.rngs)
    for (const profile of config.profiles) validate_profile(rng, profile);
  const next_index = exact(data.next_index, 'Invalid checkpoint index');
  const emitted = exact(data.emitted, 'Invalid checkpoint emitted count');
  if (next_index < 0n || next_index > total_candidates(config))
    throw new Error('Invalid checkpoint index');
  if (emitted < 0n || emitted > next_index) throw new Error('Invalid checkpoint emitted count');
  if (
    typeof data.output_offset !== 'number' ||
    !Number.isSafeInteger(data.output_offset) ||
    data.output_offset < 0
  )
    throw new Error('Invalid checkpoint output offset');
  if (!valid_fingerprint(data.output_sha256)) throw new Error('Invalid checkpoint output hash');
  if (typeof data.output !== 'string' || data.output.includes('\0') || !isAbsolute(data.output))
    throw new Error('Checkpoint output must be an absolute path');
  return {
    schema_version: data.schema_version,
    config,
    next_index,
    emitted,
    output: data.output,
    output_offset: data.output_offset,
    output_sha256: data.output_sha256,
  };
}
function canonical(path: string): string {
  const absolute = resolve(path);
  return existsSync(absolute)
    ? realpathSync(absolute)
    : resolve(realpathSync(dirname(absolute)), basename(absolute));
}
function same_file(left: string, right: string): boolean {
  if (canonical(left) === canonical(right)) return true;
  if (!existsSync(left) || !existsSync(right)) return false;
  const a = statSync(left);
  const b = statSync(right);
  return a.dev === b.dev && a.ino === b.ino;
}
// Explicit hooks permit deterministic crash/progress regressions without replacing globals.
export const runtime = {
  candidate,
  stream_at,
  atomic_json,
  lock_fd,
  now: () => performance.now() / 1000,
  progress: (value: unknown) => {
    write_all(2, Buffer.from(`${json(value)}\n`));
  },
  digest: () => createHash('sha256'),
};
export function run(args: Args): unknown {
  const cp = args.checkpoint ? canonical(args.checkpoint) : undefined;
  if (!cp) return run_locked(args);
  if (
    args.command === 'generate' &&
    args.output &&
    (same_file(args.output, cp) || same_file(args.output, `${cp}.lock`))
  )
    throw new Error('Output, checkpoint and checkpoint lock must differ');
  const release = checkpoint_lock(cp);
  try {
    return run_locked({ ...args, checkpoint: cp });
  } finally {
    release();
  }
}
function run_locked(args: Args): unknown {
  if (args.command !== 'estimate' && (args.max_candidates < 1n || args.checkpoint_every < 1n))
    throw new Error('Candidate and checkpoint limits must be positive');
  if (
    args.command !== 'estimate' &&
    (!Number.isFinite(args.progress_every) || args.progress_every < 0)
  )
    throw new Error('Progress interval must be finite and nonnegative');
  const cp = args.checkpoint;
  let config: Config;
  let index = 0n;
  let emitted = 0n;
  let output: string | undefined;
  let fd = 1;
  let offset = 0;
  const digest = runtime.digest();
  try {
    if (args.command === 'resume') {
      if (!cp) throw new Error('resume requires a checkpoint');
      const original = readFileSync(cp, 'utf8');
      const data = validate_checkpoint(parse_json(original));
      config = {
        ...data.config,
        skip_draws: data.config.skip_draws ?? [0],
        keys_per_seed: data.config.keys_per_seed ?? 1,
      };
      index = data.next_index;
      emitted = data.emitted;
      output = data.output;
      if (same_file(output, cp) || same_file(output, `${cp}.lock`))
        throw new Error('Output, checkpoint and checkpoint lock must differ');
      fd = openSync(output, 'r+');
      runtime.lock_fd(fd);
      if (json(parse_json(readFileSync(cp, 'utf8'))) !== json(parse_json(original)))
        throw new Error('Checkpoint changed before output lock; retry resume');
      let remaining = data.output_offset;
      const block = Buffer.alloc(Math.min(1024 * 1024, remaining));
      while (remaining) {
        const count = readSync(fd, block, 0, Math.min(block.length, remaining), null);
        if (!count) throw new Error('Output is shorter than checkpoint');
        digest.update(block.subarray(0, count));
        remaining -= count;
      }
      if (digest.copy().digest('hex') !== data.output_sha256)
        throw new Error('Checkpoint output prefix hash mismatch');
      ftruncateSync(fd, data.output_offset);
      offset = data.output_offset; // Sequential reads left the descriptor at the durable prefix.
    } else {
      config = make_config(args);
      if (args.command === 'estimate')
        return {
          config,
          attempts: total_candidates(config),
          timestamps: timestamp_count(config),
          note: 'Attempts include rejected candidates; no RNG or EC generation performed. No universal timestamp key mapping.',
        };
      output = args.output ? canonical(args.output) : undefined;
      if (cp && !output) throw new Error('--checkpoint requires --output');
      if (cp && existsSync(cp)) throw new Error('Checkpoint must be new and different from output');
      if (output) {
        fd = openSync(output, 'wx', 0o600);
        runtime.lock_fd(fd);
      }
    }
    const total = total_candidates(config);
    const stop = index + args.max_candidates < total ? index + args.max_candidates : total;
    const start_index = index;
    const selected = config.match_fingerprints === null ? null : new Set(config.match_fingerprints);
    let rejected = 0n;
    let nonmatching = 0n;
    const started = runtime.now();
    let last_report = started;
    const statistics = (now: number) => {
      const elapsed = Math.max(0, now - started);
      const attempted = index - start_index;
      return {
        attempted_this_run: attempted,
        next_index: index,
        total_attempts: total,
        emitted_total: emitted,
        rejected_this_run: rejected,
        nonmatching_this_run: nonmatching,
        elapsed_seconds: elapsed,
        attempts_per_second: elapsed ? Number(attempted) / elapsed : 0,
        complete: index === total,
      };
    };
    const save = () => {
      if (cp) {
        fsyncSync(fd);
        runtime.atomic_json(cp, {
          schema_version: SCHEMA,
          config,
          next_index: index,
          emitted,
          output,
          output_offset: offset,
          output_sha256: digest.copy().digest('hex'),
        });
      }
    };
    let previous_group = '';
    let rng_stream: ReturnType<typeof stream_at> | undefined;
    while (index < stop) {
      const [timestamp, rng, profile, curve, skip, key_index] = attempt_at(config, index);
      const group = json([timestamp, rng, profile, curve, skip]);
      if (group !== previous_group) {
        rng_stream = runtime.stream_at(rng, timestamp, profile, skip, key_index);
        previous_group = group;
      }
      const record = runtime.candidate(
        timestamp,
        config.unit,
        rng,
        profile,
        curve,
        skip,
        key_index,
        rng_stream
      );
      if (record === null) rejected++;
      else if (selected !== null && !selected.has(record.spki_sha256)) nonmatching++;
      else {
        const encoded = Buffer.from(`${json(record)}\n`);
        write_all(fd, encoded);
        digest.update(encoded);
        offset += encoded.length;
        emitted++;
        if (!Number.isSafeInteger(offset))
          throw new Error('Output offset exceeds safe filesystem range');
      }
      index++;
      if ((index - start_index) % args.checkpoint_every === 0n) save();
      if (args.progress_every) {
        const now = runtime.now();
        if (now - last_report >= args.progress_every) {
          const progress = statistics(now);
          runtime.progress({
            ...progress,
            event: 'progress',
            remaining_this_run: stop - index,
            eta_seconds_this_run: progress.attempts_per_second
              ? Number(stop - index) / progress.attempts_per_second
              : null,
          });
          last_report = now;
        }
      }
    }
    save();
    return statistics(runtime.now());
  } finally {
    if (fd !== 1) closeSync(fd);
  }
}
export function main(argv = process.argv.slice(2)): void {
  if (argv.includes('--help') || argv.includes('-h')) {
    console.log(HELP);
    return;
  }
  try {
    const args = parse_args(argv);
    const result = run(args);
    write_all(args.command === 'estimate' ? 1 : 2, Buffer.from(`${json(result)}\n`));
  } catch (error) {
    console.error(`timestampseeds: ${(error as Error).message}`);
    process.exitCode = 2;
  }
}
if (process.argv[1] && resolve(process.argv[1]) === fileURLToPath(import.meta.url)) main();
