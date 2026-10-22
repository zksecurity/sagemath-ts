/**
 * One baby-step table shared by every key and pattern of a search, plus the giant-step
 * walk that queries it.
 *
 * The table maps x(j·B) to j for 1 ≤ j ≤ size. Since x(−j·B) = x(j·B), a lookup of R
 * answers "is R = ±j·B" and the table spans 2·size + 1 consecutive multiples, so the
 * giant stride is W·B with W = 2·size + 1.
 *
 * Both the table and the walk visit points in chunks around a center E: E ± k·step for
 * k ≤ m come from one shared inversion and need no y-coordinate, so only the centers
 * pay for a full addition.
 */

import { type EcGroup, EcLanes, type EcPoint, ecFingerprint, ecPlusMinusX } from './ec-arith.js';

const TABLE_HALF_CHUNK = 2048;
const WALK_HALF_CHUNK = 32;

/** k·P for k = 1..m, by repeatedly adding the last multiple to the block so far. */
function multiples(group: EcGroup, P: EcPoint, m: number): { xs: bigint[]; ys: bigint[] } {
  const block = new EcLanes([P]);
  while (block.length < m) {
    const last = block.point(block.length - 1);
    if (last === null) throw new RangeError('point order is below the block size');
    const count = Math.min(block.length, m - block.length);
    const next = new EcLanes([]);
    next.xs = block.xs.slice(0, count);
    next.ys = block.ys.slice(0, count);
    next.inf = block.inf.slice(0, count);
    next.add(group, new Array(count).fill(last.x), new Array(count).fill(last.y));
    if (next.inf.includes(1)) throw new RangeError('point order is below the block size');
    block.xs.push(...next.xs);
    block.ys.push(...next.ys);
  }
  return { xs: block.xs, ys: block.ys };
}

export class EcBabyTable {
  readonly group: EcGroup;
  readonly base: EcPoint;
  readonly size: number;
  readonly width: bigint;
  readonly stride: EcPoint;
  private readonly cap: number;
  private readonly mask: number;
  private readonly tags: Uint32Array;
  private readonly values: Uint32Array;

  constructor(group: EcGroup, base: EcPoint, bits: number) {
    if (!Number.isInteger(bits) || bits < 1 || bits > 24) {
      throw new RangeError('table bits must be an integer in [1, 24]');
    }
    this.group = group;
    this.base = base;
    this.size = 2 ** bits;
    this.width = 2n * BigInt(this.size) + 1n;
    this.cap = 2 ** (bits + 1);
    this.mask = this.cap - 1;
    this.tags = new Uint32Array(this.cap);
    this.values = new Uint32Array(this.cap);
    if (2n * this.width >= group.n) throw new RangeError('table is too large for the group order');

    const m = Math.min(TABLE_HALF_CHUNK, this.size);
    const block = multiples(group, base, m);
    const span = 2 * m + 1;
    const step = group.mul(base, BigInt(span)) as EcPoint;
    const center = new EcLanes([group.mul(base, BigInt(m + 1))]);
    const plus: bigint[] = new Array(m);
    const minus: bigint[] = new Array(m);
    for (let c = m + 1; c - m <= this.size; c += span) {
      ecPlusMinusX(group, center, block.xs, block.ys, plus, minus);
      this.insert(center.xs[0] as bigint, c);
      for (let k = 0; k < m; k++) {
        this.insert(plus[k] as bigint, c + k + 1);
        this.insert(minus[k] as bigint, c - k - 1);
      }
      center.set(0, group.add(center.point(0), step));
    }
    this.stride = group.mul(base, this.width) as EcPoint;
  }

  /** j with x(j·B) = x, or 0. A hit is ±j·B up to a 2^-30 false-positive rate. */
  get(x: bigint): number {
    const f = ecFingerprint(x);
    const tag = Math.floor(f / this.cap) >>> 0;
    let i = f % this.cap;
    for (;;) {
      const v = this.values[i] as number;
      if (v === 0 || this.tags[i] === tag) return v;
      i = (i + 1) & this.mask;
    }
  }

  private insert(x: bigint, j: number): void {
    if (j > this.size) return;
    const f = ecFingerprint(x);
    let i = f % this.cap;
    while (this.values[i] !== 0) i = (i + 1) & this.mask;
    this.tags[i] = Math.floor(f / this.cap) >>> 0;
    this.values[i] = j;
  }
}

let cached: { key: string; table: EcBabyTable } | null = null;

/**
 * The table for (group, base) with at least 2^bits entries (fewer for tiny groups),
 * reusing the last one built.
 */
export function ecBabyTable(group: EcGroup, base: EcPoint, bits: number): EcBabyTable {
  const size = Math.min(bits, Math.max(1, group.n.toString(2).length - 4));
  const key = [group.p, group.a, group.b, group.n, base.x, base.y].join(':');
  if (cached !== null && cached.key === key && cached.table.size >= 2 ** size) {
    return cached.table;
  }
  cached = null;
  const table = new EcBabyTable(group, base, size);
  cached = { key, table };
  return table;
}

export interface EcGiantLane {
  start: EcPoint | null;
  /** The lane checks start − dir·u·W·B for u = 0..steps, and possibly a little further. */
  steps: number;
  dir: 1 | -1;
  /** Lanes sharing an owner stop together once `onMatch` accepts one of them. */
  owner: number;
}

/**
 * Walk every lane in lockstep. Whenever a lane's point equals e·B for some e in reach of
 * the table, call `onMatch(lane, candidates)` with the possible e (both signs, since the
 * table stores x only). Returning true retires every lane of that owner. Returns the
 * number of points checked.
 */
export function ecGiantSteps(
  table: EcBabyTable,
  lanes: readonly EcGiantLane[],
  onMatch: (lane: number, candidates: bigint[]) => boolean
): number {
  const group = table.group;
  const W = table.width;
  const maxSteps = lanes.reduce((m, l) => Math.max(m, l.steps), 0);
  const m = Math.min(WALK_HALF_CHUNK, Math.ceil(maxSteps / 2));
  const span = 2 * m + 1;
  const block = m > 0 ? multiples(group, table.stride, m) : { xs: [], ys: [] };
  const jump = group.mul(table.stride, BigInt(span)) as EcPoint;
  const done = new Uint8Array(lanes.reduce((n, l) => Math.max(n, l.owner), 0) + 1);

  let ids = lanes.map((_, i) => i);
  const centers = new EcLanes(lanes.map((l) => l.start));
  if (m > 0) {
    const shift = { x: block.xs[m - 1] as bigint, y: block.ys[m - 1] as bigint };
    centers.add(
      group,
      ids.map(() => shift.x),
      ids.map((id) => ((lanes[id] as EcGiantLane).dir > 0 ? group.mod(-shift.y) : shift.y))
    );
  }
  let plus: bigint[] = [];
  let minus: bigint[] = [];
  let qx: bigint[] = [];
  let qy: bigint[] = [];
  let checked = 0;

  const visit = (t: number, u: number, x: bigint, ownerDone: () => boolean): void => {
    if (ownerDone()) return;
    const id = ids[t] as number;
    const offset = BigInt((lanes[id] as EcGiantLane).dir * u) * W;
    let candidates: bigint[] | null = null;
    if (x < 0n) {
      candidates = [offset];
    } else {
      const j = table.get(x);
      if (j !== 0) candidates = [offset + BigInt(j), offset - BigInt(j)];
    }
    if (candidates !== null && onMatch(id, candidates)) done[(lanes[id] as EcGiantLane).owner] = 1;
  };

  for (let c = m; ids.length > 0; c += span) {
    if (plus.length < ids.length * m) {
      plus = new Array(ids.length * m);
      minus = new Array(ids.length * m);
    }
    ecPlusMinusX(group, centers, block.xs, block.ys, plus, minus);
    for (let t = 0; t < ids.length; t++) {
      const lane = lanes[ids[t] as number] as EcGiantLane;
      const ownerDone = () => done[lane.owner] === 1;
      visit(t, c, centers.inf[t] ? -1n : (centers.xs[t] as bigint), ownerDone);
      for (let k = 0; k < m; k++) {
        const i = t * m + k;
        const ahead = lane.dir > 0 ? minus[i] : plus[i];
        const behind = lane.dir > 0 ? plus[i] : minus[i];
        visit(t, c - k - 1, behind as bigint, ownerDone);
        visit(t, c + k + 1, ahead as bigint, ownerDone);
      }
    }
    checked += ids.length * span;

    const keep = ids.map((id) => {
      const lane = lanes[id] as EcGiantLane;
      return done[lane.owner] === 0 && c + m < lane.steps;
    });
    if (keep.some((k) => !k) || qx.length !== ids.length) {
      centers.filter((t) => keep[t] as boolean);
      ids = ids.filter((_, t) => keep[t]);
      qx = ids.map(() => jump.x);
      qy = ids.map((id) => ((lanes[id] as EcGiantLane).dir > 0 ? group.mod(-jump.y) : jump.y));
    }
    if (ids.length === 0) break;
    centers.add(group, qx, qy);
  }
  return checked;
}
