import { gcd, isqrt } from '../../../sagemath-ts/src/arith/misc.js';
import type { RsaSeverity } from './rsa-keys.js';
import { bitLength, exactSqrt, productTree, removeSmallFactors } from './rsa-math.js';

export interface RsaBatchFinding {
  check: 'duplicate-modulus' | 'shared-prime' | 'shared-n-minus-1';
  severity: RsaSeverity;
  /** Indices into the input list. */
  keys: number[];
  message: string;
  /** The shared prime, or the shared factor of n − 1. */
  shared?: bigint;
  /** Keys this finding factored, n = p·q. */
  factored?: Array<{ key: number; p: bigint; q: bigint }>;
}

/**
 * Bernstein's batch GCD: gcd(nᵢ, ∏_{j≠i} nⱼ) for every i, through a product tree and a
 * remainder tree of squares, in quasi-linear time.
 */
export function batchGcd(values: readonly bigint[]): bigint[] {
  if (values.some((v) => v < 1n)) throw new Error('batchGcd needs positive values');
  if (values.length < 2) return values.map(() => 1n);
  const tree = productTree(values);
  let rems = tree.at(-1)!;
  for (let level = tree.length - 2; level >= 0; level--) {
    const parents = rems;
    rems = tree[level]!.map((v, i) => parents[i >> 1]! % (v * v));
  }
  return values.map((v, i) => gcd(v, rems[i]! / v));
}

function groupIndices(moduli: readonly bigint[]): Map<bigint, number[]> {
  const groups = new Map<bigint, number[]>();
  moduli.forEach((n, i) => {
    const list = groups.get(n);
    if (list) list.push(i);
    else groups.set(n, [i]);
  });
  return groups;
}

/** Duplicate moduli, and moduli sharing a prime with another (factored by batch GCD). */
export function findSharedPrimes(moduli: readonly bigint[]): RsaBatchFinding[] {
  const groups = groupIndices(moduli);
  const findings: RsaBatchFinding[] = [];
  for (const keys of groups.values()) {
    if (keys.length > 1) {
      findings.push({
        check: 'duplicate-modulus',
        severity: 'vulnerable',
        keys,
        message: `${keys.length} keys share one modulus: each private key decrypts for the others`,
      });
    }
  }

  const unique = [...groups.keys()];
  const byFactor = new Map<bigint, Set<bigint>>();
  const note = (factor: bigint, n: bigint) => {
    const set = byFactor.get(factor) ?? new Set<bigint>();
    set.add(n);
    byFactor.set(factor, set);
  };
  batchGcd(unique).forEach((g, u) => {
    const n = unique[u]!;
    if (g === 1n) return;
    if (g !== n) return note(g, n);
    // Both primes of n are shared, with different moduli.
    for (const other of unique) {
      const h = other === n ? 1n : gcd(n, other);
      if (h !== 1n && h !== n) note(h, n);
    }
  });

  for (const [factor, sharing] of byFactor) {
    const keys = [...sharing].flatMap((n) => groups.get(n)!).sort((a, b) => a - b);
    findings.push({
      check: 'shared-prime',
      severity: 'broken',
      keys,
      shared: factor,
      message: `a ${bitLength(factor)}-bit prime divides several moduli: batch GCD factors them`,
      factored: keys.map((key) => {
        const [p, q] = [factor, moduli[key]! / factor];
        return p <= q ? { key, p, q } : { key, p: q, q: p };
      }),
    });
  }
  return findings;
}

/**
 * Factors n when p ≡ q ≡ 1 (mod r): (n − 1)/r = (x + y) + r·xy for p = 1 + rx, q = 1 + ry.
 * Succeeds within `maxSteps` when r is above about n^(1/4)/2^8.
 */
export function factorWithNMinus1Divisor(
  n: bigint,
  r: bigint,
  maxSteps = 1 << 16
): [bigint, bigint] | null {
  if (r < 2n || (n - 1n) % r !== 0n) return null;
  const m = (n - 1n) / r;
  const base = m % r;
  // x + y ≥ 2√(xy) ≈ 2√(m/r).
  const lower = 2n * isqrt(m / r);
  const skip = lower > base ? (lower - base) / r - 1n : 0n;
  let s = base + (skip > 0n ? skip : 0n) * r;
  for (let step = 0; step < maxSteps && s <= m; step++, s += r) {
    const xy = (m - s) / r;
    const root = exactSqrt(s * s - 4n * xy);
    if (root === null || ((s + root) & 1n) === 1n) continue;
    const p = 1n + r * ((s - root) / 2n);
    if (p > 1n && p < n && n % p === 0n) return [p, n / p];
  }
  return null;
}

/**
 * Pairs of moduli where n₁ − 1 and n₂ − 1 share a factor of at least `minBits` bits made of
 * primes above `smoothBound`: the primes of both keys are likely 1 modulo a common r.
 * Each such key is then factored with `factorWithNMinus1Divisor` when r is large enough.
 */
export function findSharedNMinus1Factors(
  moduli: readonly bigint[],
  options: { minBits?: number; smoothBound?: number; maxSteps?: number } = {}
): RsaBatchFinding[] {
  const { minBits = 64, smoothBound = 1 << 16, maxSteps } = options;
  const groups = groupIndices(moduli.filter((n) => n > 2n));
  const unique = [...groups.keys()];
  const rough = unique.map((n) => removeSmallFactors(n - 1n, smoothBound));
  const flagged = batchGcd(rough).flatMap((g, u) => (g > 1n ? [u] : []));

  const findings: RsaBatchFinding[] = [];
  for (let a = 0; a < flagged.length; a++) {
    for (let b = a + 1; b < flagged.length; b++) {
      const [u, v] = [flagged[a]!, flagged[b]!];
      const r = gcd(rough[u]!, rough[v]!);
      if (bitLength(r) < minBits) continue;
      const keys = [...groups.get(unique[u]!)!, ...groups.get(unique[v]!)!].sort((x, y) => x - y);
      const factored = keys.flatMap((key) => {
        const pq = factorWithNMinus1Divisor(moduli[key]!, r, maxSteps);
        return pq === null ? [] : [{ key, p: pq[0], q: pq[1] }];
      });
      const what = `n − 1 shares a ${bitLength(r)}-bit factor r across keys: the primes are likely 1 mod r`;
      findings.push({
        check: 'shared-n-minus-1',
        severity: factored.length > 0 ? 'broken' : 'warning',
        keys,
        shared: r,
        message: factored.length > 0 ? `${what}, which factors n` : what,
        ...(factored.length > 0 ? { factored } : {}),
      });
    }
  }
  return findings;
}
