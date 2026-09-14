/**
 * PARI nth-prime lookup and successor search.
 * Reference: reference/pari/src/basemath/prime.c:prime_table_find_n/prime_successor.
 * @see Deviation: PARI prime lookup sizing and cache
 */
import { forprime, isPrime, isqrt } from './ifactor.js';

// The original 64-bit prime_table, stored as [prime, ordinal].
const PRIME_TABLE: readonly (readonly [bigint, bigint])[] = [
  [0n, 0n],
  [7919n, 1000n],
  [17389n, 2000n],
  [27449n, 3000n],
  [37813n, 4000n],
  [48611n, 5000n],
  [59359n, 6000n],
  [70657n, 7000n],
  [81799n, 8000n],
  [93179n, 9000n],
  [104729n, 10000n],
  [224737n, 20000n],
  [350377n, 30000n],
  [479909n, 40000n],
  [611953n, 50000n],
  [746773n, 60000n],
  [882377n, 70000n],
  [1020379n, 80000n],
  [1159523n, 90000n],
  [1299709n, 100000n],
  [2750159n, 200000n],
  [7368787n, 500000n],
  [15485863n, 1000000n],
  [32452843n, 2000000n],
  [86028121n, 5000000n],
  [179424673n, 10000000n],
  [373587883n, 20000000n],
  [982451653n, 50000000n],
  [2038074743n, 100000000n],
  [4000000483n, 189961831n],
  [4222234741n, 200000000n],
  [5336500537n, 250000000n],
  [6461335109n, 300000000n],
  [7594955549n, 350000000n],
  [8736028057n, 400000000n],
  [9883692017n, 450000000n],
  [11037271757n, 500000000n],
  [13359555403n, 600000000n],
  [15699342107n, 700000000n],
  [18054236957n, 800000000n],
  [20422213579n, 900000000n],
  [22801763489n, 1000000000n],
  [47055833459n, 2000000000n],
  [71856445751n, 3000000000n],
  [97011687217n, 4000000000n],
  [122430513841n, 5000000000n],
  [148059109201n, 6000000000n],
  [173862636221n, 7000000000n],
  [200000000507n, 8007105083n],
  [225898512559n, 9000000000n],
  [252097800623n, 10000000000n],
  [384489816343n, 15000000000n],
  [518649879439n, 20000000000n],
  [654124187867n, 25000000000n],
  [790645490053n, 30000000000n],
  [928037044463n, 35000000000n],
  [1066173339601n, 40000000000n],
  [1344326694119n, 50000000000n],
  [1624571841097n, 60000000000n],
  [1906555030411n, 70000000000n],
  [2190026988349n, 80000000000n],
  [2474799787573n, 90000000000n],
  [2760727302517n, 100000000000n],
];

function nativeError(name: string, message: string): never {
  const error = new Error(message);
  error.name = name;
  throw error;
}

/**
 * Iterate a bounded interval using PARI's sieve/next-prime strategy family.
 * Integer sizing keeps every endpoint exact. Short intervals at large offsets
 * use BPSW next-prime search, avoiding a square-root-sized sieve setup.
 */
function* intervalPrimes(a: bigint, b: bigint): Generator<bigint> {
  const root = isqrt(b);
  const width = b - a + 1n;
  if (b <= BigInt(Number.MAX_SAFE_INTEGER) && root <= 4n * width) {
    for (const p of forprime(Number(a), Number(b))) yield BigInt(p);
    return;
  }
  if (a <= 2n && b >= 2n) yield 2n;
  let candidate = a <= 3n ? 3n : a | 1n;
  for (; candidate <= b; candidate += 2n) {
    if (isPrime(candidate)) yield candidate;
  }
}

/** Return the nth prime (one indexed), as PARI's prime(long). */
export function prime(n: bigint): bigint {
  if (n < -(1n << 63n) || n >= 1n << 63n) {
    nativeError('OverflowError', 'Python int too large to convert to C long');
  }
  if (n <= 0n) nativeError('PariError', 'domain error in prime: n <= 0');
  // The source chooses the nearest table row, then backs down if it exceeds n.
  let p = 0n;
  let ordinal = 0n;
  for (const [value, index] of PRIME_TABLE) {
    if (index > n) break;
    p = value;
    ordinal = index;
  }
  if (ordinal === n) return p;
  // With no process-global pari_PRIMES cache, start the first interval at 2.
  let remaining = n - ordinal;
  let left = p < 2n ? 2n : p + 1n;
  while (true) {
    // Source: ceil(4 * remaining * log(Q+1)); bit length is an exact upper bound.
    const extent = 4n * remaining * BigInt((left + 1n).toString(2).length);
    const right = left + extent;
    for (const value of intervalPrimes(left, right)) {
      if (--remaining === 0n) return value;
    }
    left = right + 1n;
  }
}
