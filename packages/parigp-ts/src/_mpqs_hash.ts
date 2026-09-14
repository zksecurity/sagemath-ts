import { PariError } from './errors.js';
/** MPQS's relation-key subset of PARI language/hash.c, using 64-bit GEN words.
 * Reference: hash.c:24–189, 340–371 and mpqs.c:frel_add.
 * @see Deviation: MPQS relation hash representation
 */
export interface RelationData {
  Y: bigint;
  relp: number[];
}
const MASK = (1n << 64n) - 1n;
const HASH_PRIMES = [
  53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241, 786433,
  1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319, 201326611, 402653189,
  805306457, 1610612741,
];
const glue = (h: bigint, a: bigint): bigint => (404936533n * h + a) & MASK;
function integerHash(n: bigint): bigint {
  let magnitude = n < 0n ? -n : n;
  const words: bigint[] = [];
  while (magnitude) {
    words.push(magnitude & MASK);
    magnitude >>= 64n;
  }
  const sign = n < 0n ? -1n : n > 0n ? 1n : 0n;
  let h = glue(1n << 57n, ((sign << 62n) & MASK) + BigInt(words.length + 2));
  for (const word of words) h = glue(h, word);
  return h;
}
/** hash_GEN([Y, Vecsmall(relp)]), including native GEN type/length/sign words.
 * @see Deviation: MPQS relation hash representation
 */
export function relationHash(row: RelationData): bigint {
  let h = (22n << 57n) + BigInt(row.relp.length + 1);
  for (const word of row.relp) h = glue(h, BigInt(word) & MASK);
  return glue(glue((17n << 57n) + 3n, integerHash(row.Y)), h);
}
interface Entry<T> {
  row: T;
  hash: bigint;
  next: Entry<T> | undefined;
}
/** Native hash_init/hash_insert2/hash_keys_GEN with MPQS frel_add deduplication.
 * minimum is a nonnegative integer, as for native hash_init's unsigned word.
 * Relation values must remain unchanged while stored, as for native hash keys.
 * @see Deviation: MPQS relation hash representation
 */
export class RelationTable<T extends RelationData = RelationData> {
  private count = 0;
  private index: number;
  private limit: number;
  private buckets: Array<Entry<T> | undefined>;
  constructor(minimum: number) {
    this.index = HASH_PRIMES.findIndex((p) => p > minimum);
    if (this.index < 0) throw new PariError('overflow in hash table [too large].');
    const length = HASH_PRIMES[this.index]!;
    this.limit = Math.ceil(length * 0.65);
    this.buckets = Array(length);
  }
  get size(): number {
    return this.count;
  }
  add(row: T): void {
    const hash = relationHash(row);
    let index = Number(hash % BigInt(this.buckets.length));
    for (let entry = this.buckets[index]; entry; entry = entry.next)
      if (
        entry.hash === hash &&
        entry.row.Y === row.Y &&
        entry.row.relp.length === row.relp.length &&
        entry.row.relp.every((word, i) => word === row.relp[i])
      )
        return;
    this.count++;
    if (this.size > this.limit && this.index < HASH_PRIMES.length - 1) {
      const length = HASH_PRIMES[++this.index]!;
      const next: Array<Entry<T> | undefined> = Array(length);
      for (const head of this.buckets) {
        let entry = head;
        while (entry) {
          const after = entry.next;
          const index = Number(entry.hash % BigInt(length));
          entry.next = next[index];
          next[index] = entry;
          entry = after;
        }
      }
      this.buckets = next;
      this.limit = Math.ceil(length * 0.65);
    }
    index = Number(hash % BigInt(this.buckets.length));
    this.buckets[index] = { row, hash, next: this.buckets[index] };
  }
  *values(): IterableIterator<T> {
    for (const head of this.buckets)
      for (let entry = head; entry; entry = entry.next) yield entry.row;
  }
}
