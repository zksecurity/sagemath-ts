import { createHash } from 'node:crypto';
import { EcGroup } from '../../../zksecurity-cheatsheets/src/attacks/ec-arith.js';
import { p256, secp256k1 } from '../../../zksecurity-cheatsheets/src/curves/params.js';
import { type Rng, WIDTHS, type WordStream, words } from './rngs.js';

export const PROFILES = [
  'concat-be-reject',
  'concat-le-reject',
  'concat-be-mod',
  'low-byte-reject',
  'buffer-le-reject',
  'python-getrandbits-reject',
  'python-randbytes-reject',
  'java-nextbytes-reject',
] as const;
export type Profile = (typeof PROFILES)[number];
export type Unit = 'seconds' | 'milliseconds';
const PROFILE_RNG: Partial<Record<Profile, Rng>> = {
  'python-getrandbits-reject': 'python-mt19937',
  'python-randbytes-reject': 'python-mt19937',
  'java-nextbytes-reject': 'java-random',
};
export const MAX_SKIP_DRAWS = 1_000_000;
export const MAX_KEYS_PER_SEED = 1_000;
export const CURVES = {
  p256: {
    native: 'prime256v1',
    order: 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551n,
    prefix: '3059301306072a8648ce3d020106082a8648ce3d030107034200',
  },
  secp256k1: {
    native: 'secp256k1',
    order: 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141n,
    prefix: '3056301006072a8648ce3d020106052b8104000a034200',
  },
} as const;
export type Curve = keyof typeof CURVES;
const groups = new Map<Curve, EcGroup>();
function group(curve: Curve): EcGroup {
  let result = groups.get(curve);
  if (!result) {
    const params = curve === 'p256' ? p256 : secp256k1;
    result = new EcGroup({
      p: params.field.p,
      a: params.a,
      b: params.b,
      n: params.order,
      g: params.generator,
    });
    groups.set(curve, result);
  }
  return result;
}
export interface Config {
  start: bigint;
  end: bigint;
  step: bigint;
  unit: Unit;
  rngs: Rng[];
  profiles: Profile[];
  curves: Curve[];
  skip_draws?: number[];
  keys_per_seed?: number;
  match_fingerprints: string[] | null;
}
export interface Candidate {
  timestamp: bigint;
  unit: Unit;
  rng: Rng;
  profile: Profile;
  curve: Curve;
  spki_sha256: string;
  spki_der_base64: string;
  sec1_compressed_hex: string;
  skip_draws?: number;
  key_index?: number;
}

/** Parse datetime.fromisoformat's calendar/basic/week forms without rounding fractional seeds. */
export function parse_time(value: string, unit: Unit): bigint {
  if (unit !== 'seconds' && unit !== 'milliseconds') throw new Error('Unknown timestamp unit');
  if (/^[+-]?\d(?:_?\d)*$/.test(value.trim())) return BigInt(value.trim().replaceAll('_', ''));
  for (const fraction of value.matchAll(/[.,](\d+)/g)) {
    if (/[1-9]/.test(fraction[1]!.slice(6))) throw new Error(`Timestamp is not aligned to ${unit}`);
  }
  if (/[+-]\d{2}:?\d{2}(?::?\d{2})?[.,]/.test(value))
    throw new Error('Fractional timezone offsets are not supported');
  const match =
    /^(\d{4})(?:-?(\d{2})-?(\d{2})|-?W(\d{2})(?:-?([1-7]))?)[\s\S](\d{2})(?:(:?)(\d{2})(?:\7(\d{2}))?)?(?:[.,](\d+))?(Z|[+-]\d{2}(?::?\d{2}(?::?\d{2})?)?)?$/.exec(
      value
    );
  if (!match) throw new Error('Invalid ISO timestamp');
  const [
    ,
    yearText,
    monthText,
    dayText,
    weekText,
    weekDayText,
    hourText,
    ,
    minuteText,
    secondText,
    fraction = '',
    zone,
  ] = match;
  if (!zone) throw new Error('ISO timestamps require Z or an explicit timezone offset');
  const year = Number(yearText);
  const hour = Number(hourText);
  const minute = Number(minuteText ?? 0);
  const second = Number(secondText ?? 0);
  const date = new Date(0);
  date.setUTCFullYear(year, Number(monthText ?? 1) - 1, Number(dayText ?? 4));
  if (year < 1 || hour > 23 || minute > 59 || second > 59) throw new Error('Invalid ISO timestamp');
  if (weekText) {
    const week = Number(weekText);
    date.setUTCDate(
      date.getUTCDate() -
        ((date.getUTCDay() + 6) % 7) +
        (week - 1) * 7 +
        Number(weekDayText ?? 1) -
        1
    );
    const thursday = new Date(date);
    thursday.setUTCDate(date.getUTCDate() + 3 - ((date.getUTCDay() + 6) % 7));
    if (week < 1 || week > 53 || thursday.getUTCFullYear() !== year)
      throw new Error('Invalid ISO timestamp');
  } else if (
    date.getUTCFullYear() !== year ||
    date.getUTCMonth() !== Number(monthText) - 1 ||
    date.getUTCDate() !== Number(dayText)
  )
    throw new Error('Invalid ISO timestamp');
  let offset = 0;
  if (zone !== 'Z') {
    const digits = zone.slice(1).replaceAll(':', '');
    const h = Number(digits.slice(0, 2));
    const m = Number(digits.slice(2, 4));
    const s = Number(digits.slice(4, 6));
    offset = (h * 3600 + m * 60 + s) * (zone[0] === '-' ? -1 : 1);
    if (Math.abs(offset) >= 86400) throw new Error('Invalid timezone offset');
  }
  const micros =
    BigInt(date.getTime()) * 1000n +
    BigInt(hour * 3600 + minute * 60 + second - offset) * 1_000_000n +
    BigInt(fraction.slice(0, 6).padEnd(6, '0'));
  const divisor = unit === 'seconds' ? 1_000_000n : 1000n;
  if (micros % divisor) throw new Error(`Timestamp is not aligned to ${unit}`);
  return micros / divisor;
}

export function validate_profile(rng: Rng, profile: Profile): void {
  if (!Object.hasOwn(WIDTHS, rng)) throw new Error(`Unknown RNG: ${rng}`);
  if (!PROFILES.includes(profile)) throw new Error('Unknown derivation profile');
  if (PROFILE_RNG[profile] && PROFILE_RNG[profile] !== rng)
    throw new Error(`${profile} requires --rng ${PROFILE_RNG[profile]}`);
}
export function draws_per_key(rng: Rng, profile: Profile): number {
  validate_profile(rng, profile);
  return profile === 'low-byte-reject' ? 32 : Math.ceil(32 / WIDTHS[rng]);
}
export function scalar_from_words(
  stream: WordStream,
  profile: Profile,
  order: bigint
): bigint | null {
  if (!PROFILES.includes(profile) || typeof order !== 'bigint' || order < 2n)
    throw new Error('Invalid derivation profile or subgroup order');
  const little = [
    'concat-le-reject',
    'python-getrandbits-reject',
    'python-randbytes-reject',
    'java-nextbytes-reject',
  ].includes(profile);
  const material: number[] = [];
  while (material.length < 32) {
    const next = stream.next();
    if (next.done) throw new Error('RNG stream exhausted');
    const [value, width] = next.value;
    if (profile === 'low-byte-reject') material.push(Number(value & 255n));
    else {
      if (
        !Number.isSafeInteger(width) ||
        width < 1 ||
        width > 32 ||
        value < 0n ||
        value >= 1n << BigInt(width * 8)
      )
        throw new Error('Invalid RNG word');
      for (let i = 0; i < width; i++)
        material.push(Number((value >> BigInt(8 * (little ? i : width - i - 1))) & 255n));
    }
  }
  const bytes = Buffer.from(material.slice(0, 32));
  if (profile === 'buffer-le-reject' || profile === 'python-getrandbits-reject') bytes.reverse();
  const value = BigInt(`0x${bytes.toString('hex')}`);
  if (profile.endsWith('-mod')) return (value % (order - 1n)) + 1n;
  return value >= 1n && value < order ? value : null;
}
export function stream_at(
  rng: Rng,
  seed: bigint,
  profile: Profile,
  skip_draws = 0,
  key_index = 0
): WordStream {
  if (!Number.isSafeInteger(skip_draws) || skip_draws < 0 || skip_draws > MAX_SKIP_DRAWS)
    throw new Error(`Skipped draws must be between 0 and ${MAX_SKIP_DRAWS}`);
  if (!Number.isSafeInteger(key_index) || key_index < 0 || key_index >= MAX_KEYS_PER_SEED)
    throw new Error(`Key index must be between 0 and ${MAX_KEYS_PER_SEED - 1}`);
  const count = skip_draws + key_index * draws_per_key(rng, profile);
  const stream = words(rng, seed);
  for (let i = 0; i < count; i++) stream.next();
  return stream;
}
export function scalar(
  rng: Rng,
  seed: bigint,
  profile: Profile,
  order: bigint,
  skip_draws = 0,
  key_index = 0
): bigint | null {
  return scalar_from_words(stream_at(rng, seed, profile, skip_draws, key_index), profile, order);
}
export function candidate(
  timestamp: bigint,
  unit: Unit,
  rng: Rng,
  profile: Profile,
  curve: Curve,
  skip_draws = 0,
  key_index = 0,
  stream?: WordStream
): Candidate | null {
  const definition = CURVES[curve];
  if (!definition) throw new Error('Unknown curve');
  validate_profile(rng, profile);
  const value = stream
    ? scalar_from_words(stream, profile, definition.order)
    : scalar(rng, timestamp, profile, definition.order, skip_draws, key_index);
  if (value === null) return null;
  const point = group(curve).mul(group(curve).g, value)!;
  const x = point.x.toString(16).padStart(64, '0');
  const y = point.y.toString(16).padStart(64, '0');
  const der = Buffer.from(definition.prefix + '04' + x + y, 'hex');
  const record: Candidate = {
    timestamp,
    unit,
    rng,
    profile,
    curve,
    spki_sha256: createHash('sha256').update(der).digest('hex'),
    spki_der_base64: der.toString('base64'),
    sec1_compressed_hex: (point.y & 1n ? '03' : '02') + x,
  };
  if (skip_draws || key_index) Object.assign(record, { skip_draws, key_index });
  return record;
}
export function combination_at(config: Config, index: bigint): [bigint, Rng, Profile, Curve] {
  const curve = config.curves[Number(index % BigInt(config.curves.length))]!;
  index /= BigInt(config.curves.length);
  const profile = config.profiles[Number(index % BigInt(config.profiles.length))]!;
  index /= BigInt(config.profiles.length);
  const rng = config.rngs[Number(index % BigInt(config.rngs.length))]!;
  index /= BigInt(config.rngs.length);
  return [config.start + index * config.step, rng, profile, curve];
}
export function timestamp_count(config: Config): bigint {
  return (config.end - config.start + config.step - 1n) / config.step;
}
export function total_candidates(config: Config): bigint {
  return (
    timestamp_count(config) *
    BigInt(config.rngs.length) *
    BigInt(config.profiles.length) *
    BigInt(config.curves.length) *
    BigInt((config.skip_draws ?? [0]).length) *
    BigInt(config.keys_per_seed ?? 1)
  );
}
export function attempt_at(
  config: Config,
  index: bigint
): [bigint, Rng, Profile, Curve, number, number] {
  const keys = BigInt(config.keys_per_seed ?? 1);
  const key = Number(index % keys);
  index /= keys;
  const skips = config.skip_draws ?? [0];
  const skip = skips[Number(index % BigInt(skips.length))]!;
  return [...combination_at(config, index / BigInt(skips.length)), skip, key];
}
