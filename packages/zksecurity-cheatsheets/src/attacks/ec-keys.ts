/**
 * Weak elliptic-curve public keys and curve parameters, on any short-Weierstrass curve
 * given as raw parameters.
 *
 * The key checks recover private scalars that were generated badly: tiny, a small value
 * shifted by whole bytes, one 32-bit word repeated, inside a small interval, or close to
 * another key in the batch. They follow the checks of badkeys and Google's
 * paranoid_crypto, reimplemented from the underlying math so they work on curves those
 * tools cannot parse. Every recovered scalar is verified by recomputing the public key.
 *
 * Most searches share one baby-step table of multiples of H = 2^(8J)·G, where 8J is the
 * largest byte shift below the order. Multiplying a key by 2^(8k) costs only doublings,
 * and turns "s = i·2^(8(J−k))" into "2^(8k)·Q = i·H", so one table serves every shift.
 */

import {
  factor,
  inverse_mod,
  is_prime,
  sqrt_mod,
  trial_division,
} from '../../../sagemath-ts/src/arith/misc.js';
import { ValueError } from '../../../sagemath-ts/src/errors.js';
import { type PrimeFieldCurveParams, hasseCofactor } from '../curves.js';
import { EcGroup, EcLanes, type EcPoint } from './ec-arith.js';
import { ecKangaroo } from './ec-kangaroo.js';
import { type EcGiantLane, ecBabyTable, ecGiantSteps } from './ec-table.js';

export type { EcPoint } from './ec-arith.js';

/** An affine point, or a SEC1 encoding (compressed, uncompressed or bare x‖y) as hex or bytes. */
export type EcPublicKey = EcPoint | string | Uint8Array;

export interface EcCurveReport {
  name: string;
  fieldBits: number;
  orderBits: number;
  /** Broken parameters. Key searches refuse to run while any are present. */
  problems: string[];
  /** Known attacks that apply to the curve itself. */
  weaknesses: string[];
  /** From the parameters, or forced by the Hasse bound when the order exceeds 4√p. */
  cofactor: bigint | null;
  curveOrder: bigint | null;
  /** Factorization of the generator order, as far as it could be pushed. */
  orderFactors: [bigint, bigint][];
  /** Composite part of the order left unfactored (null when fully factored). */
  unfactored: bigint | null;
  largestPrimeFactor: bigint | null;
  /** Smallest k ≤ maxEmbeddingDegree with ℓ | p^k − 1, ℓ the largest prime factor. */
  embeddingDegree: number | null;
  /** #E = p, or the order is p: Smart's attack solves discrete logs in linear time. */
  anomalous: boolean;
}

export type EcKeyProblem = 'malformed' | 'not-on-curve' | 'identity' | 'wrong-subgroup';

export interface EcKeyCheck {
  index: number;
  valid: boolean;
  problem?: EcKeyProblem;
  point: EcPoint | null;
}

/** How a recovered scalar s was found. `negated` means s = n − (the described value). */
export type EcKeyPattern =
  | { kind: 'small'; value: bigint; negated: boolean }
  | { kind: 'shifted'; value: bigint; shift: number; negated: boolean }
  | { kind: 'repeated-word'; word: bigint; words: number; negated: boolean }
  | { kind: 'interval'; bits: number }
  | { kind: 'related'; key: number; difference: bigint };

export interface EcRecoveredKey {
  index: number;
  scalar: bigint;
  pattern: EcKeyPattern;
}

/** s_a = s_b + difference (mod n). */
export interface EcKeyRelation {
  a: number;
  b: number;
  difference: bigint;
}

export interface EcTableOptions {
  /**
   * log2 of the baby-step table size, at most 24 (2^22 entries take 64 MB). Default:
   * sized to the work, between 16 and 22.
   */
  tableBits?: number;
}

export interface EcPatternOptions extends EcTableOptions {
  /** Plain small scalars are searched below 2^smallBits. Default 32. */
  smallBits?: number;
}

export interface EcIntervalOptions {
  /** Search s < 2^bits. Default: the largest interval that fits `maxSteps`. */
  bits?: number;
  /** Point additions before giving up. Default 2^22, or enough for `bits` when it is set. */
  maxSteps?: number;
  seed?: bigint;
}

export interface EcIntervalResult {
  bits: number;
  steps: number;
  recovered: EcRecoveredKey[];
  relations: EcKeyRelation[];
}

export interface EcWeakKeyOptions extends EcPatternOptions {
  /** Report pairs whose scalars differ by less than 2^closeBits. Default 24; 0 skips. */
  closeBits?: number;
  /** Also run the kangaroo below 2^intervalBits. Off by default. */
  intervalBits?: number;
  maxSteps?: number;
  seed?: bigint;
}

export interface EcWeakKeyReport {
  curve: EcCurveReport;
  invalid: EcKeyCheck[];
  recovered: EcRecoveredKey[];
  close: EcKeyRelation[];
}

const DEFAULT_KANGAROO_STEPS = 2 ** 22;
const KANGAROO_MARGIN = 6;

function bitLength(v: bigint): number {
  return v === 0n ? 0 : v.toString(2).length;
}

function groupOf(curve: PrimeFieldCurveParams): EcGroup {
  return new EcGroup({
    p: curve.field.p,
    a: curve.a,
    b: curve.b,
    n: curve.order,
    g: curve.generator,
  });
}

function withinHasse(p: bigint, order: bigint): boolean {
  const t = order - p - 1n;
  return t * t <= 4n * p;
}

function factorOrder(n: bigint): { factors: [bigint, bigint][]; rest: bigint | null } {
  const found = new Map<bigint, bigint>();
  let rest = n;
  const push = (q: bigint) => found.set(q, (found.get(q) ?? 0n) + 1n);
  let start = 2n;
  while (rest > 1n && !is_prime(rest) && bitLength(rest) > 160) {
    const q = trial_division(rest, 1n << 20n, start);
    if (q === rest || q > 1n << 20n) break;
    push(q);
    rest /= q;
    start = q;
  }
  if (rest > 1n && is_prime(rest)) {
    push(rest);
    rest = 1n;
  } else if (rest > 1n && bitLength(rest) <= 160) {
    for (const [q, e] of factor(rest)) for (let i = 0n; i < e; i++) push(q);
    rest = 1n;
  }
  const factors = [...found.entries()].sort((x, y) => (x[0] < y[0] ? -1 : 1));
  return { factors, rest: rest > 1n ? rest : null };
}

/**
 * Validate curve parameters and report the attacks that apply to the curve itself:
 * small largest prime factor (Pohlig–Hellman), small embedding degree (MOV/Frey–Rück),
 * anomalous order (Smart), singular curve.
 */
export function ecCurveReport(
  curve: PrimeFieldCurveParams,
  options: { maxEmbeddingDegree?: number } = {}
): EcCurveReport {
  const maxK = options.maxEmbeddingDegree ?? 50;
  const group = groupOf(curve);
  const { p, a, b, n, g } = group;
  const problems: string[] = [];
  const weaknesses: string[] = [];

  if (!is_prime(p)) problems.push('field modulus p is not prime');
  if (!group.isOnCurve(curve.generator)) problems.push('generator is not on the curve');
  else if (n <= 0n || group.mul(g, n) !== null) problems.push('n·G is not the identity');
  if (!is_prime(n)) problems.push('order n is not prime; pass the prime subgroup order');

  if (group.mod(4n * a * a * a + 27n * b * b) === 0n) {
    weaknesses.push(
      'singular curve: discrete logs map to the additive or multiplicative group of a field'
    );
  }

  let cofactor = curve.cofactor ?? hasseCofactor(p, n) ?? null;
  let curveOrder = cofactor === null ? null : cofactor * n;
  if (curveOrder !== null && !withinHasse(p, curveOrder)) {
    problems.push('cofactor·n violates the Hasse bound');
    cofactor = null;
    curveOrder = null;
  }

  const { factors, rest } = n > 1n ? factorOrder(n) : { factors: [], rest: null };
  const largestPrimeFactor =
    rest === null && factors.length > 0
      ? (factors[factors.length - 1] as [bigint, bigint])[0]
      : null;
  if (largestPrimeFactor !== null && bitLength(largestPrimeFactor) < 160) {
    const bits = bitLength(largestPrimeFactor);
    weaknesses.push(
      `largest prime factor of the order has ${bits} bits: Pohlig–Hellman and rho take about 2^${Math.ceil(bits / 2)} steps`
    );
  }

  const anomalous = curveOrder === p || n === p || largestPrimeFactor === p;
  if (anomalous) weaknesses.push("anomalous curve: Smart's p-adic attack solves discrete logs");

  let embeddingDegree: number | null = null;
  const ell = largestPrimeFactor;
  if (ell !== null && p % ell !== 0n) {
    const q = p % ell;
    let t = q;
    for (let k = 1; k <= maxK; k++) {
      if (t === 1n) {
        embeddingDegree = k;
        break;
      }
      t = (t * q) % ell;
    }
  }
  if (embeddingDegree !== null && embeddingDegree * bitLength(p) < 2048) {
    weaknesses.push(
      `embedding degree ${embeddingDegree}: MOV/Frey–Rück moves discrete logs into a ${embeddingDegree * bitLength(p)}-bit field`
    );
  }

  return {
    name: curve.name,
    fieldBits: bitLength(p),
    orderBits: bitLength(n),
    problems,
    weaknesses,
    cofactor,
    curveOrder,
    orderFactors: factors,
    unfactored: rest,
    largestPrimeFactor,
    embeddingDegree,
    anomalous,
  };
}

function decodeKey(
  group: EcGroup,
  key: EcPublicKey
): { point: EcPoint | null } | { problem: EcKeyProblem } {
  if (typeof key === 'object' && !(key instanceof Uint8Array)) return { point: key };
  let bytes: Uint8Array;
  if (typeof key === 'string') {
    const hex = key.replace(/^0x/i, '').replace(/[\s:]/g, '');
    if (!/^[0-9a-f]*$/i.test(hex) || hex.length % 2 !== 0) return { problem: 'malformed' };
    bytes = Uint8Array.from(hex.match(/../g) ?? [], (h) => Number.parseInt(h, 16));
  } else {
    bytes = key;
  }
  const len = Math.ceil(bitLength(group.p) / 8);
  const int = (from: number, to: number) =>
    bytes.subarray(from, to).reduce((acc, byte) => (acc << 8n) | BigInt(byte), 0n);
  if (bytes.length === 1 && bytes[0] === 0) return { point: null };
  if (bytes.length === 2 * len + 1 && bytes[0] === 4) {
    return { point: { x: int(1, len + 1), y: int(len + 1, 2 * len + 1) } };
  }
  if (bytes.length === 2 * len) return { point: { x: int(0, len), y: int(len, 2 * len) } };
  if (bytes.length === len + 1 && (bytes[0] === 2 || bytes[0] === 3)) {
    const x = int(1, len + 1);
    if (x >= group.p) return { problem: 'not-on-curve' };
    const y = sqrt_mod(group.mod(((x * x + group.a) % group.p) * x + group.b), group.p);
    if (y === null) return { problem: 'not-on-curve' };
    return { point: { x, y: (y & 1n) === BigInt(bytes[0] & 1) ? y : group.mod(-y) } };
  }
  return { problem: 'malformed' };
}

/**
 * Check every key: well formed, on the curve, not the identity, and in the subgroup of
 * prime order n. The subgroup test is skipped when the cofactor is known to be 1.
 */
export function validateEcKeys(
  curve: PrimeFieldCurveParams,
  keys: readonly EcPublicKey[]
): EcKeyCheck[] {
  const group = groupOf(curve);
  const cofactor = curve.cofactor ?? hasseCofactor(group.p, group.n);
  const checks: EcKeyCheck[] = keys.map((key, index) => {
    const decoded = decodeKey(group, key);
    if ('problem' in decoded) return { index, valid: false, problem: decoded.problem, point: null };
    const { point } = decoded;
    if (point === null) return { index, valid: false, problem: 'identity', point };
    if (!group.isOnCurve(point)) return { index, valid: false, problem: 'not-on-curve', point };
    return { index, valid: true, point };
  });
  if (cofactor !== 1n) {
    const pending = checks.filter((c) => c.valid);
    const products = group.mulMany(
      pending.map((c) => c.point),
      group.n
    );
    pending.forEach((c, i) => {
      if (products[i] !== null) {
        c.valid = false;
        c.problem = 'wrong-subgroup';
      }
    });
  }
  return checks;
}

interface Prepared {
  group: EcGroup;
  checks: EcKeyCheck[];
  valid: { index: number; point: EcPoint }[];
  /** Largest byte shift: 8·maxShift < bits(n). */
  maxShift: number;
}

function prepare(curve: PrimeFieldCurveParams, keys: readonly EcPublicKey[]): Prepared {
  const group = groupOf(curve);
  if (!group.isOnCurve(curve.generator) || group.mul(group.g, group.n) !== null) {
    throw new ValueError('generator is not on the curve or n·G is not the identity');
  }
  if (!is_prime(group.n))
    throw new ValueError('order n must be prime: pass the prime subgroup order');
  const checks = validateEcKeys(curve, keys);
  const valid = checks
    .filter((c) => c.valid)
    .map((c) => ({ index: c.index, point: c.point as EcPoint }));
  return { group, checks, valid, maxShift: Math.floor((bitLength(group.n) - 1) / 8) };
}

function tableBase(prep: Prepared): EcPoint {
  return prep.group.shl(prep.group.g, 8 * prep.maxShift) as EcPoint;
}

/** log2 of the table size minimizing table + giant steps for the given walk lengths. */
function autoTableBits(totalRange: number): number {
  return Math.min(22, Math.max(16, Math.round(Math.log2(Math.sqrt(totalRange / 2)))));
}

function stepsFor(range: bigint, size: number, width: bigint): number {
  const reach = range - 1n - BigInt(size);
  return reach <= 0n ? 0 : Number((reach + width - 1n) / width);
}

function signedPattern(e: bigint, shift: number): EcKeyPattern {
  const negated = e < 0n;
  let value = negated ? -e : e;
  let s = shift;
  while (value > 0n && (value & 255n) === 0n) {
    value >>= 8n;
    s++;
  }
  return s === 0
    ? { kind: 'small', value, negated }
    : { kind: 'shifted', value, shift: s, negated };
}

interface PatternLane extends EcGiantLane {
  /** Exponents e in [0, range) are reached; steps is derived from it. */
  range: bigint;
  multiplier: bigint;
  describe: (e: bigint) => EcKeyPattern;
}

function runPatternLanes(
  prep: Prepared,
  lanes: PatternLane[],
  bits: number,
  found: EcRecoveredKey[]
): void {
  const { group } = prep;
  const n = group.n;
  const table = ecBabyTable(group, tableBase(prep), bits);
  for (const lane of lanes) lane.steps = stepsFor(lane.range, table.size, table.width);
  ecGiantSteps(table, lanes, (laneIndex, candidates) => {
    const lane = lanes[laneIndex] as PatternLane;
    const key = prep.valid[lane.owner] as { index: number; point: EcPoint };
    for (const e of candidates) {
      const s = (((e * lane.multiplier) % n) + n) % n;
      if (s !== 0n && group.equals(group.mul(group.g, s), key.point)) {
        found.push({ index: key.index, scalar: s, pattern: lane.describe(e) });
        return true;
      }
    }
    return false;
  });
}

/**
 * Instant check: every key whose scalar is ±s with s at most the table size (at least
 * 2^16, more when a larger table is already cached for this curve).
 */
export function findSmallEcKeys(
  curve: PrimeFieldCurveParams,
  keys: readonly EcPublicKey[],
  options: EcTableOptions = {}
): EcRecoveredKey[] {
  const prep = prepare(curve, keys);
  const starts = prep.group.mulMany(
    prep.valid.map((k) => k.point),
    1n << BigInt(8 * prep.maxShift)
  );
  const lanes: PatternLane[] = starts.map((start, owner) => ({
    start,
    steps: 0,
    dir: 1,
    owner,
    range: 0n,
    multiplier: 1n,
    describe: (e) => signedPattern(e, 0),
  }));
  const found: EcRecoveredKey[] = [];
  runPatternLanes(prep, lanes, options.tableBits ?? 16, found);
  return found.sort((x, y) => x.index - y.index);
}

/**
 * Search the patterns of paranoid_crypto's weak EC private key check, for all keys at
 * once: s < 2^smallBits; s = i·2^(8j) for i < 2^32 and every byte shift j; s made of
 * one 32-bit word repeated 2 or more times. Negations of these are found within the
 * table size.
 */
export function findPatternedEcKeys(
  curve: PrimeFieldCurveParams,
  keys: readonly EcPublicKey[],
  options: EcPatternOptions = {}
): EcRecoveredKey[] {
  return patterned(prepare(curve, keys), options);
}

function patterned(prep: Prepared, options: EcPatternOptions): EcRecoveredKey[] {
  const { group, maxShift } = prep;
  const n = group.n;
  const smallBits = options.smallBits ?? 32;
  const shifts = maxShift + 1;
  const chains = group.shlChain(
    prep.valid.map((k) => k.point),
    8,
    shifts
  );
  const lanes: PatternLane[] = [];
  let totalRange = 0;
  prep.valid.forEach((_, owner) => {
    const chain = chains[owner] as (EcPoint | null)[];
    for (let k = shifts - 1; k >= 0; k--) {
      const shift = maxShift - k;
      const rangeBits = shift === 0 ? smallBits : 32;
      totalRange += 2 ** rangeBits;
      lanes.push({
        start: chain[k] as EcPoint | null,
        steps: 0,
        dir: 1,
        owner,
        range: 1n << BigInt(rangeBits),
        multiplier: 1n << BigInt(8 * shift),
        describe: (e) => signedPattern(e, shift),
      });
    }
  });

  const H = 1n << BigInt(8 * maxShift);
  const maxWords = Math.ceil(bitLength(n) / 32);
  const repeats: { words: number; unit: bigint; transform: bigint }[] = [];
  for (let words = 2; words <= maxWords; words++) {
    let unit = 0n;
    for (let w = 0; w < words; w++) unit |= 1n << BigInt(32 * w);
    if (unit % n === 0n) continue;
    repeats.push({ words, unit: unit % n, transform: (H * inverse_mod(unit, n)) % n });
  }
  const repeatStarts = group.mulMany(
    prep.valid.flatMap((k) => repeats.map(() => k.point)),
    prep.valid.flatMap(() => repeats.map((r) => r.transform))
  );
  prep.valid.forEach((_, owner) => {
    repeats.forEach((r, i) => {
      totalRange += 2 ** 32;
      lanes.push({
        start: repeatStarts[owner * repeats.length + i] as EcPoint | null,
        steps: 0,
        dir: 1,
        owner,
        range: 1n << 32n,
        multiplier: r.unit,
        describe: (e) => ({
          kind: 'repeated-word',
          word: e < 0n ? -e : e,
          words: r.words,
          negated: e < 0n,
        }),
      });
    });
  });

  const found: EcRecoveredKey[] = [];
  runPatternLanes(prep, lanes, options.tableBits ?? autoTableBits(totalRange), found);
  return found.sort((x, y) => x.index - y.index);
}

/**
 * Pairs of keys whose scalars differ by less than 2^bits (default 24), such as keys
 * derived by incrementing a counter. A difference of 0 is a repeated key. Quadratic in
 * the number of keys.
 */
export function findCloseEcKeys(
  curve: PrimeFieldCurveParams,
  keys: readonly EcPublicKey[],
  options: EcTableOptions & { bits?: number } = {}
): EcKeyRelation[] {
  const prep = prepare(curve, keys);
  return closePairs(prep, options.bits ?? 24, options.tableBits);
}

function closePairs(prep: Prepared, bits: number, tableBits: number | undefined): EcKeyRelation[] {
  const { group, valid } = prep;
  const scaled = group.mulMany(
    valid.map((k) => k.point),
    1n << BigInt(8 * prep.maxShift)
  );
  const pairs: [number, number][] = [];
  for (let a = 0; a < valid.length; a++)
    for (let b = a + 1; b < valid.length; b++) pairs.push([a, b]);
  if (pairs.length === 0) return [];

  const diffs = new EcLanes(pairs.map(([a]) => scaled[a] as EcPoint));
  diffs.add(
    group,
    pairs.map(([, b]) => (scaled[b] as EcPoint).x),
    pairs.map(([, b]) => group.mod(-(scaled[b] as EcPoint).y))
  );
  const range = 1n << BigInt(bits);
  const totalRange = 2 * pairs.length * Number(range);
  const table = ecBabyTable(group, tableBase(prep), tableBits ?? autoTableBits(totalRange));
  const steps = stepsFor(range, table.size, table.width);
  const lanes: EcGiantLane[] = pairs.flatMap((_, owner) => [
    { start: diffs.point(owner), steps, dir: 1 as const, owner },
    { start: diffs.point(owner), steps, dir: -1 as const, owner },
  ]);

  const out: EcKeyRelation[] = [];
  ecGiantSteps(table, lanes, (laneIndex, candidates) => {
    const [a, b] = pairs[(lanes[laneIndex] as EcGiantLane).owner] as [number, number];
    const qa = (valid[a] as { point: EcPoint }).point;
    const qb = (valid[b] as { point: EcPoint }).point;
    const delta = group.sub(qa, qb);
    for (const d of candidates) {
      if (d < range && -d < range && group.equals(group.mul(group.g, d), delta)) {
        out.push({
          a: (valid[a] as { index: number }).index,
          b: (valid[b] as { index: number }).index,
          difference: d,
        });
        return true;
      }
    }
    return false;
  });
  return out.sort((x, y) => x.a - y.a || x.b - y.b);
}

/**
 * Pollard kangaroo for scalars below 2^bits, all keys sharing one tame herd. Without
 * `bits`, the interval is the largest that `maxSteps` (default 2^22 additions, a few
 * seconds) covers with high probability. A key not found is only probably outside.
 */
export function findEcKeysInInterval(
  curve: PrimeFieldCurveParams,
  keys: readonly EcPublicKey[],
  options: EcIntervalOptions = {}
): EcIntervalResult {
  return interval(prepare(curve, keys), options);
}

function interval(prep: Prepared, options: EcIntervalOptions): EcIntervalResult {
  const K = Math.max(1, prep.valid.length);
  const maxSteps =
    options.maxSteps ??
    (options.bits === undefined
      ? DEFAULT_KANGAROO_STEPS
      : Math.ceil(KANGAROO_MARGIN * Math.sqrt(K * 2 ** options.bits)));
  const bits =
    options.bits ??
    Math.max(1, Math.floor(2 * Math.log2(maxSteps / KANGAROO_MARGIN) - Math.log2(K)));
  const n = prep.group.n;
  const result = ecKangaroo(
    prep.group,
    prep.valid.map((k) => k.point),
    bits,
    maxSteps,
    options.seed ?? 1n
  );
  const recovered = [...result.solved.entries()]
    .map(([k, scalar]) => ({
      index: (prep.valid[k] as { index: number }).index,
      scalar,
      pattern: { kind: 'interval', bits } as EcKeyPattern,
    }))
    .sort((x, y) => x.index - y.index);
  const relations = result.relations.map((r) => ({
    a: (prep.valid[r.a] as { index: number }).index,
    b: (prep.valid[r.b] as { index: number }).index,
    difference: r.difference > n / 2n ? r.difference - n : r.difference,
  }));
  return { bits, steps: result.steps, recovered, relations };
}

/**
 * Every check at once: curve weaknesses, invalid keys, patterned scalars, close pairs
 * (and the scalars they reveal from recovered neighbours), and optionally the kangaroo.
 */
export function findWeakEcKeys(
  curve: PrimeFieldCurveParams,
  keys: readonly EcPublicKey[],
  options: EcWeakKeyOptions = {}
): EcWeakKeyReport {
  const report = ecCurveReport(curve);
  if (report.problems.length > 0) {
    return { curve: report, invalid: [], recovered: [], close: [] };
  }
  const prep = prepare(curve, keys);
  const invalid = prep.checks.filter((c) => !c.valid);
  const recovered = patterned(prep, options);
  const known = new Map(recovered.map((r) => [r.index, r.scalar]));
  if (options.intervalBits !== undefined) {
    const found = interval(
      { ...prep, valid: prep.valid.filter((k) => !known.has(k.index)) },
      { bits: options.intervalBits, maxSteps: options.maxSteps, seed: options.seed }
    );
    for (const r of found.recovered) {
      recovered.push(r);
      known.set(r.index, r.scalar);
    }
  }
  const closeBits = options.closeBits ?? 24;
  const close = closeBits > 0 ? closePairs(prep, closeBits, options.tableBits) : [];

  const { group } = prep;
  const n = group.n;
  const points = new Map(prep.valid.map((k) => [k.index, k.point]));
  let grew = true;
  while (grew) {
    grew = false;
    for (const { a, b, difference } of close) {
      for (const [from, to, d] of [
        [b, a, difference],
        [a, b, -difference],
      ] as const) {
        const base = known.get(from);
        if (base === undefined || known.has(to)) continue;
        const s = (((base + d) % n) + n) % n;
        if (!group.equals(group.mul(group.g, s), points.get(to) as EcPoint)) continue;
        known.set(to, s);
        recovered.push({
          index: to,
          scalar: s,
          pattern: { kind: 'related', key: from, difference: d },
        });
        grew = true;
      }
    }
  }
  recovered.sort((x, y) => x.index - y.index);
  return { curve: report, invalid, recovered, close };
}
