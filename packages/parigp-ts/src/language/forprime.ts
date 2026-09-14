/** Word-prime iteration corresponding to PARI language/forprime.c.
 * @see Deviation: PARI Galois integer helper boundaries
 */
import { forprime, nextprime } from '../ifactor.js';

const WORD = 1n << 64n;
// Match the native initialized prime table used by the bundled oracle.
const TABLE_LIMIT = 500000n;
let primeTableCache: { table: number[]; last: bigint; sieveLimit: bigint } | undefined;
/** Initialize on first iteration so bounded factorization can share this iterator. */
function primeTable() {
  if (!primeTableCache) {
    const table = Array.from(forprime(2, Number(TABLE_LIMIT))),
      last = BigInt(table[table.length - 1]!);
    primeTableCache = { table, last, sieveLimit: last * last };
  }
  return primeTableCache;
}
const SEGMENT = 65536n;

/** Private state for the existing number-valued Forprime facade. */
export class WordPrimeIterator {
  private candidate: bigint;
  private tableIndex: number | undefined;
  private segment: Generator<number> | undefined;

  constructor(start: bigint) {
    const { table: TABLE, last: TABLE_LAST } = primeTable();
    this.candidate = start < 2n ? 2n : start;
    if (this.candidate < TABLE_LAST) {
      let lo = 0,
        hi = TABLE.length;
      while (lo < hi) {
        const mid = Math.floor((lo + hi) / 2);
        if (BigInt(TABLE[mid]!) < this.candidate) lo = mid + 1;
        else hi = mid;
      }
      this.tableIndex = lo;
    } else if (this.candidate <= TABLE_LIMIT) {
      // Native a >= maxprime selects the sieve, whose first interval starts
      // after maxprimelim. The final cached prime is skipped in this case.
      this.candidate = TABLE_LIMIT + 1n;
    }
  }

  next(): bigint {
    const { table: TABLE, sieveLimit: SIEVE_LIMIT } = primeTable();
    for (;;) {
      if (this.tableIndex !== undefined) {
        if (this.tableIndex < TABLE.length) return BigInt(TABLE[this.tableIndex++]!);
        this.tableIndex = undefined;
        this.candidate = TABLE_LIMIT + 1n;
      }
      if (this.segment) {
        const value = this.segment.next();
        if (!value.done) return BigInt(value.value);
        this.segment = undefined;
      }
      if (this.candidate >= WORD) return 0n;
      if (this.candidate <= SIEVE_LIMIT) {
        const hi =
          this.candidate + SEGMENT - 1n < SIEVE_LIMIT ? this.candidate + SEGMENT - 1n : SIEVE_LIMIT;
        this.segment = forprime(Number(this.candidate), Number(hi));
        this.candidate = hi + 1n;
      } else {
        const p = nextprime(this.candidate);
        this.candidate = p + 1n;
        if (p >= WORD) {
          this.candidate = WORD;
          return 0n;
        }
        return p;
      }
    }
  }
}
