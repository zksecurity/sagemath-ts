/** CPython 3.11 Objects/setobject.c: iteration order of `set(iterable)`.
 * Source: https://github.com/python/cpython/blob/v3.11.8/Objects/setobject.c
 * Adapted under the Python Software Foundation license; see CPYTHON-LICENSE.txt.
 * Only insertion into a fresh set is modelled (no deletions, so no dummy entries).
 * Iteration walks the hash table in slot order, which depends on the hashes,
 * the insertion order of colliding keys and the resize schedule.
 * @see Deviation: Python set iteration order
 */

const PySet_MINSIZE = 8;
const LINEAR_PROBES = 9;
const PERTURB_SHIFT = 5n;
const PyHASH_MODULUS = (1n << 61n) - 1n;

interface Entry<T> {
  key: T;
  hash: bigint;
}

/** CPython 3.11 Objects/longobject.c long_hash on 64-bit platforms. */
export function python_int_hash(n: bigint): bigint {
  const magnitude = (n < 0n ? -n : n) % PyHASH_MODULUS;
  const h = n < 0n ? -magnitude : magnitude;
  return h === -1n ? -2n : h;
}

/** set_insert_clean: place a key known to be absent into a table without dummies. */
function insertClean<T>(table: (Entry<T> | undefined)[], mask: number, entry: Entry<T>): void {
  let perturb = entry.hash;
  let i = Number(entry.hash & BigInt(mask));
  for (;;) {
    if (table[i] === undefined) {
      table[i] = entry;
      return;
    }
    if (i + LINEAR_PROBES <= mask) {
      for (let j = 1; j <= LINEAR_PROBES; j++) {
        if (table[i + j] === undefined) {
          table[i + j] = entry;
          return;
        }
      }
    }
    perturb >>= PERTURB_SHIFT;
    i = Number((BigInt(i) * 5n + 1n + perturb) & BigInt(mask));
  }
}

/**
 * Return the distinct items of `items` in the order Python iterates
 * `set(items)`, where `hash` is the item's Python hash and `eq` its equality.
 */
export function python_set<T>(
  items: Iterable<T>,
  hash: (item: T) => bigint,
  eq: (a: T, b: T) => boolean
): T[] {
  let mask = PySet_MINSIZE - 1;
  let table: (Entry<T> | undefined)[] = new Array(PySet_MINSIZE);
  let used = 0;
  for (const key of items) {
    // Py_hash_t is converted to size_t for slot selection and perturbation.
    const h = BigInt.asUintN(64, hash(key));
    let perturb = h;
    let i = Number(h & BigInt(mask));
    let inserted = false;
    search: for (;;) {
      let probes = i + LINEAR_PROBES <= mask ? LINEAR_PROBES : 0;
      let j = i;
      do {
        const entry = table[j];
        if (entry === undefined) {
          table[j] = { key, hash: h };
          inserted = true;
          break search;
        }
        if (entry.hash === h && eq(entry.key, key)) break search;
        j++;
      } while (probes-- > 0);
      perturb >>= PERTURB_SHIFT;
      i = Number((BigInt(i) * 5n + 1n + perturb) & BigInt(mask));
    }
    if (!inserted) continue;
    used++;
    if (used * 5 < mask * 3) continue;
    // set_table_resize(so, used > 50000 ? used * 2 : used * 4)
    const minused = used > 50000 ? used * 2 : used * 4;
    let newsize = PySet_MINSIZE;
    while (newsize <= minused) newsize *= 2;
    const old = table;
    table = new Array(newsize);
    mask = newsize - 1;
    for (const entry of old) if (entry !== undefined) insertClean(table, mask, entry);
  }
  const order: T[] = [];
  for (const entry of table) if (entry !== undefined) order.push(entry.key);
  return order;
}
