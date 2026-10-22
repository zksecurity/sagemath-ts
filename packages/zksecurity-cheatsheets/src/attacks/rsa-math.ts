import { eratosthenes, gcd, isqrt } from '../../../sagemath-ts/src/arith/misc.js';

export function bitLength(x: bigint): number {
  return x === 0n ? 0 : (x < 0n ? -x : x).toString(2).length;
}

function squareTable(m: number): Uint8Array {
  const table = new Uint8Array(m);
  for (let i = 0; i < m; i++) table[(i * i) % m] = 1;
  return table;
}

const SIEVE_MODULI = [64, 63, 65, 11, 17, 19, 23, 29, 31, 37, 41, 43, 47];
const SIEVE_TABLES = SIEVE_MODULI.map(squareTable);
const SQUARES_64 = SIEVE_TABLES[0]!;
const SQUARES_63 = SIEVE_TABLES[1]!;
const SQUARES_65 = SIEVE_TABLES[2]!;
const SQUARES_11 = SIEVE_TABLES[3]!;

/** The square root of `x` when `x` is a perfect square, else null. */
export function exactSqrt(x: bigint): bigint | null {
  if (x < 0n) return null;
  if (!SQUARES_64[Number(x & 63n)]) return null;
  const r = Number(x % 45045n);
  if (!SQUARES_63[r % 63] || !SQUARES_65[r % 65] || !SQUARES_11[r % 11]) return null;
  const root = isqrt(x);
  return root * root === x ? root : null;
}

/**
 * Fermat's search over a = a0, a0 + step, a0 + 2·step, ...: whenever a² − N is a square b²,
 * `accept(a, b)` decides whether the split a² − b² = N is the one wanted.
 *
 * Squareness is sieved with word-size residues, so a step costs no bigint arithmetic.
 */
export function fermatSearch<T>(
  N: bigint,
  a0: bigint,
  step: bigint,
  iterations: number,
  accept: (a: bigint, b: bigint) => T | null
): T | null {
  const k = SIEVE_MODULI.length;
  const aRes = new Int32Array(k);
  const nRes = new Int32Array(k);
  const sRes = new Int32Array(k);
  for (let j = 0; j < k; j++) {
    const m = BigInt(SIEVE_MODULI[j]!);
    aRes[j] = Number(((a0 % m) + m) % m);
    nRes[j] = Number(N % m);
    sRes[j] = Number(step % m);
  }
  for (let i = 0; i < iterations; i++) {
    let candidate = true;
    for (let j = 0; j < k; j++) {
      const m = SIEVE_MODULI[j]!;
      const r = (aRes[j]! * aRes[j]! - nRes[j]!) % m;
      if (!SIEVE_TABLES[j]![r < 0 ? r + m : r]) {
        candidate = false;
        break;
      }
    }
    if (candidate) {
      const a = a0 + BigInt(i) * step;
      const b = exactSqrt(a * a - N);
      if (b !== null) {
        const result = accept(a, b);
        if (result !== null) return result;
      }
    }
    for (let j = 0; j < k; j++) aRes[j] = (aRes[j]! + sRes[j]!) % SIEVE_MODULI[j]!;
  }
  return null;
}

/** Levels of a product tree, leaves first; the last level holds the product of all values. */
export function productTree(values: readonly bigint[]): bigint[][] {
  const levels: bigint[][] = [values.slice()];
  for (let level = levels[0]!; level.length > 1; ) {
    const next: bigint[] = [];
    for (let i = 0; i < level.length; i += 2) {
      next.push(i + 1 < level.length ? level[i]! * level[i + 1]! : level[i]!);
    }
    levels.push(next);
    level = next;
  }
  return levels;
}

export function product(values: readonly bigint[]): bigint {
  return values.length === 0 ? 1n : productTree(values).at(-1)![0]!;
}

const primeLists = new Map<number, bigint[]>();
const primeProducts = new Map<number, bigint>();

/** Primes up to and including `bound`, cached per bound. */
export function primesUpTo(bound: number): bigint[] {
  let primes = primeLists.get(bound);
  if (primes === undefined) {
    primes = eratosthenes(BigInt(bound));
    primeLists.set(bound, primes);
  }
  return primes;
}

export function primeProduct(bound: number): bigint {
  let p = primeProducts.get(bound);
  if (p === undefined) {
    p = product(primesUpTo(bound));
    primeProducts.set(bound, p);
  }
  return p;
}

/** gcd(a, b) for a small `a` against a possibly huge `b`. */
export function gcdReduced(a: bigint, b: bigint): bigint {
  return a === 0n ? b : gcd(a, b % a);
}

/** `x` with every prime factor up to `bound` divided out. */
export function removeSmallFactors(x: bigint, bound: number): bigint {
  let g = gcdReduced(x, primeProduct(bound));
  while (g > 1n) {
    x /= g;
    g = gcd(x, g);
  }
  return x;
}
