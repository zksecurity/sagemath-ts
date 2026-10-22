/**
 * Parallel Pollard kangaroo (van Oorschot–Wiener) for many keys at once.
 *
 * Tame kangaroos start at known multiples r·G with r in [0, 2^bits); wild ones start at
 * Q_k + r·G. Every kangaroo jumps by d_i·G with i taken from x, so two kangaroos that
 * ever land on the same point follow the same path and meet at the next distinguished
 * point. A tame–wild meeting reveals a key. The tame herd is shared, so solving K keys
 * costs about √K times one key rather than K times. Wild–wild meetings between two keys
 * give their difference, which solves both once either is known.
 */

import { type EcGroup, EcLanes, type EcPoint, ecFingerprint } from './ec-arith.js';

const JUMPS = 64;

export interface EcKangarooRelation {
  a: number;
  b: number;
  /** s_a − s_b mod n. */
  difference: bigint;
}

export interface EcKangarooResult {
  solved: Map<number, bigint>;
  relations: EcKangarooRelation[];
  steps: number;
}

interface Walker {
  /** −1 for tame; otherwise the target index. */
  key: number;
  dist: bigint;
}

interface Landmark {
  key: number;
  dist: bigint;
  odd: boolean;
}

/** splitmix64, for reproducible start points and jump sizes. */
class Rng {
  private state: bigint;

  constructor(seed: bigint) {
    this.state = BigInt.asUintN(64, seed);
  }

  next(): bigint {
    this.state = BigInt.asUintN(64, this.state + 0x9e3779b97f4a7c15n);
    let z = this.state;
    z = BigInt.asUintN(64, (z ^ (z >> 30n)) * 0xbf58476d1ce4e5b9n);
    z = BigInt.asUintN(64, (z ^ (z >> 27n)) * 0x94d049bb133111ebn);
    return z ^ (z >> 31n);
  }

  below(bound: bigint): bigint {
    let r = 0n;
    let range = 1n;
    while (range < bound << 64n) {
      r = (r << 64n) | this.next();
      range <<= 64n;
    }
    return r % bound;
  }
}

/**
 * Search for s_k in roughly [0, 2^bits) with Q_k = s_k·G, stopping once every target is
 * solved or after about `maxSteps` point additions. Unsolved targets are only probably
 * outside the interval: the walk is randomized.
 */
export function ecKangaroo(
  group: EcGroup,
  targets: readonly EcPoint[],
  bits: number,
  maxSteps: number,
  seed: bigint
): EcKangarooResult {
  const n = group.n;
  const K = targets.length;
  const solved = new Map<number, bigint>();
  const relations: EcKangarooRelation[] = [];
  if (K === 0) return { solved, relations, steps: 0 };

  const rng = new Rng(seed);
  const interval = 1n << BigInt(bits);
  const wildPerKey = Math.max(1, Math.ceil(128 / K));
  const herd = 2 * K * wildPerKey;
  const meanJump = Math.max(1, (herd * Math.sqrt(2 ** bits / K)) / 4);
  const perWalker = (3 * Math.sqrt(K * 2 ** bits)) / herd;
  const dpMask = 2 ** Math.max(0, Math.floor(Math.log2(Math.max(1, perWalker / 16)))) - 1;

  const jumpSize: bigint[] = [];
  for (let i = 0; i < JUMPS; i++) {
    jumpSize.push(1n + rng.below(BigInt(Math.ceil(2 * meanJump))));
  }
  const jumps = group.mulMany(
    jumpSize.map(() => group.g),
    jumpSize
  );
  if (jumps.some((J) => J === null))
    throw new RangeError('generator order is below the jump sizes');
  const jx = jumps.map((J) => (J as EcPoint).x);
  const jy = jumps.map((J) => (J as EcPoint).y);

  const unsolved = (): number[] => targets.map((_, k) => k).filter((k) => !solved.has(k));
  let nextWild = 0;
  const walkers: Walker[] = [];
  for (let i = 0; i < herd; i++) {
    walkers.push(i % 2 === 0 ? { key: -1, dist: 0n } : { key: (i >> 1) % K, dist: 0n });
  }
  const lanes = new EcLanes(walkers.map(() => null));

  const startOf = (w: Walker): EcPoint | null => {
    if (w.key < 0) {
      w.dist = rng.below(interval);
      return group.mul(group.g, w.dist);
    }
    w.dist = rng.below(interval >> 1n);
    return group.add(targets[w.key] as EcPoint, group.mul(group.g, w.dist));
  };
  const reseed = (t: number): void => {
    const w = walkers[t] as Walker;
    if (w.key >= 0 && solved.has(w.key)) {
      const open = unsolved();
      w.key = open.length > 0 ? (open[nextWild++ % open.length] as number) : -1;
    }
    lanes.set(t, startOf(w));
  };
  for (const w of walkers) w.dist = rng.below(w.key < 0 ? interval : interval >> 1n);
  const offsets = group.mulMany(
    walkers.map(() => group.g),
    walkers.map((w) => w.dist)
  );
  walkers.forEach((w, t) => {
    const offset = offsets[t] as EcPoint | null;
    lanes.set(t, w.key < 0 ? offset : group.add(targets[w.key] as EcPoint, offset));
  });

  const known = (key: number, dist: bigint): bigint | null => {
    if (key < 0) return dist;
    const s = solved.get(key);
    return s === undefined ? null : s + dist;
  };
  const accept = (key: number, s: bigint): boolean => {
    if (solved.has(key)) return false;
    const r = ((s % n) + n) % n;
    if (!group.equals(group.mul(group.g, r), targets[key] as EcPoint)) return false;
    solved.set(key, r);
    let grew = true;
    while (grew) {
      grew = false;
      for (const rel of relations) {
        const sa = solved.get(rel.a);
        const sb = solved.get(rel.b);
        if (sa !== undefined && sb === undefined) {
          const sbNew = (((sa - rel.difference) % n) + n) % n;
          if (group.equals(group.mul(group.g, sbNew), targets[rel.b] as EcPoint)) {
            solved.set(rel.b, sbNew);
            grew = true;
          }
        } else if (sb !== undefined && sa === undefined) {
          const saNew = (sb + rel.difference) % n;
          if (group.equals(group.mul(group.g, saNew), targets[rel.a] as EcPoint)) {
            solved.set(rel.a, saNew);
            grew = true;
          }
        }
      }
    }
    for (let t = 0; t < walkers.length; t++) {
      const w = walkers[t] as Walker;
      if (w.key >= 0 && solved.has(w.key)) reseed(t);
    }
    return true;
  };

  const landmarks = new Map<number, Landmark>();
  const qx: bigint[] = new Array(herd).fill(0n);
  const qy: bigint[] = new Array(herd).fill(0n);
  const jumpIndex = new Uint8Array(herd);
  let steps = 0;

  while (steps < maxSteps && solved.size < K) {
    for (let t = 0; t < herd; t++) {
      const w = walkers[t] as Walker;
      if (lanes.inf[t] && !(w.key >= 0 && accept(w.key, -w.dist))) reseed(t);
      let f = ecFingerprint(lanes.xs[t] as bigint);
      if (((f >>> 6) & dpMask) === 0) {
        const odd = ((lanes.ys[t] as bigint) & 1n) === 1n;
        const here = known(w.key, w.dist);
        const seen = landmarks.get(f);
        if (seen === undefined) {
          landmarks.set(f, { key: w.key, dist: w.dist, odd });
        } else {
          const same = seen.odd === odd;
          const there = known(seen.key, seen.dist);
          let restart = false;
          if (here !== null && there !== null) {
            restart = same;
          } else if (here !== null) {
            accept(seen.key, (same ? here : -here) - seen.dist);
            restart = same;
          } else if (there !== null) {
            restart = !accept(w.key, (same ? there : -there) - w.dist);
          } else if (same && seen.key === w.key) {
            restart = true;
          } else if (same) {
            relations.push({
              a: w.key,
              b: seen.key,
              difference: (((seen.dist - w.dist) % n) + n) % n,
            });
            restart = true;
          }
          if (restart) reseed(t);
        }
        if (solved.size === K) break;
        f = ecFingerprint(lanes.xs[t] as bigint);
      }
      const i = f & (JUMPS - 1);
      jumpIndex[t] = i;
      qx[t] = jx[i] as bigint;
      qy[t] = jy[i] as bigint;
    }
    if (solved.size === K) break;
    lanes.add(group, qx, qy);
    for (let t = 0; t < herd; t++) {
      const w = walkers[t] as Walker;
      w.dist += jumpSize[jumpIndex[t] as number] as bigint;
    }
    steps += herd;
  }
  return { solved, relations, steps };
}
