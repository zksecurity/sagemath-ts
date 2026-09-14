/** Multidimensional enumeration, ported from sage/misc/mrange.py. */

/**
 * Iterate a Cartesian product with the last coordinate varying fastest.
 * Inputs are restarted when a coordinate rolls over; no input is materialized.
 * The default converter copies each row, so callers may reverse/mutate it.
 */
export function* _xmrange_iter<T, U = T[]>(
  iter_list: Iterable<T>[],
  typ: (row: T[]) => U = ((row: T[] = []) => row.slice()) as (row: T[]) => U
): Generator<U> {
  if (iter_list.length === 0) {
    yield (typ as () => U)();
    return;
  }
  // The native infinite-input guard avoids walking an infinite prefix when
  // the final input is empty. Unknown iterables are assumed finite, as in Sage.
  type Sized = Iterable<T> & {
    is_finite?: () => boolean;
    cardinality?: () => unknown;
    length?: number;
  };
  const size = (s: Sized): unknown => {
    try {
      if (s.cardinality) return s.cardinality();
    } catch (error) {
      if (!(error instanceof Error) || error.name !== 'AttributeError') throw error;
    }
    return s.length;
  };
  const infinite = iter_list.some((input) => {
    const s = input as Sized;
    try {
      if (s.is_finite) return !s.is_finite();
    } catch (error) {
      if (error instanceof Error && error.name === 'ValueError') return false;
      if (!(error instanceof Error) || error.name !== 'AttributeError') throw error;
    }
    try {
      const n = size(s);
      return n === Infinity || n === 'Infinity';
    } catch (error) {
      if (
        error instanceof Error &&
        ['TypeError', 'AttributeError', 'NotImplementedError'].includes(error.name)
      )
        return false;
      throw error;
    }
  });
  if (infinite) {
    let n: unknown;
    for (const input of iter_list) {
      try {
        const m = size(input as Sized);
        if (m !== undefined) n = m;
      } catch (error) {
        if (!(error instanceof Error) || error.name !== 'TypeError') throw error;
      }
    }
    if (n === undefined) {
      const error = new Error(
        "cannot access local variable 'n' where it is not associated with a value"
      );
      error.name = 'UnboundLocalError';
      throw error;
    }
    if (n === 0 || n === 0n) return;
  }
  const iters = iter_list.map((input) => input[Symbol.iterator]());
  const row: T[] = [];
  for (let i = 0; i < iters.length - 1; i++) {
    const r = iters[i]!.next();
    if (r.done) return;
    row.push(r.value);
  }
  let place = iters.length - 1;
  while (place >= 0) {
    const r = iters[place]!.next();
    if (r.done) {
      place--;
    } else {
      row[place] = r.value;
      if (place < iters.length - 1) {
        place++;
        iters[place] = iter_list[place]![Symbol.iterator]();
      } else yield typ(row);
    }
  }
}
