import { ValueError } from '../errors.js';
import { Integer } from '../rings/integer_ring.js';
import { Rational } from '../rings/rational.js';

/** The list/scalar MatrixArgs adapter used by the binary matrix constructor.
 * @see Deviation: Binary Constructor and Row Ownership
 */
export function* entries(
  nrows: number,
  ncols: number,
  input: unknown,
  convert: (value: unknown) => number
): Generator<[number, number, number]> {
  // misc_c.sized_iter checks for excess input before yielding the last item.
  function* sized(input: unknown, size: number): Generator<unknown> {
    if (input instanceof Rational) input = [input]; // Rational.__getitem__ exposes its sole coefficient.
    if (input == null || typeof (input as Iterable<unknown>)[Symbol.iterator] !== 'function') {
      const name =
        input == null
          ? 'NoneType'
          : typeof input === 'number'
            ? Number.isInteger(input)
              ? 'int'
              : 'float'
            : typeof input === 'bigint'
              ? 'int'
              : typeof input === 'boolean'
                ? 'bool'
                : input instanceof Integer
                  ? 'sage.rings.integer.Integer'
                  : input instanceof Rational
                    ? 'sage.rings.rational.Rational'
                    : 'sage.rings.finite_rings.integer_mod.IntegerMod_int';
      throw new TypeError(`'${name}' object is not iterable`);
    }
    const it = (input as Iterable<unknown>)[Symbol.iterator]();
    const check = () => {
      if (!it.next().done)
        throw new ValueError(`sequence too long (expected length ${size}, got more)`);
    };
    if (!size) check();
    for (let i = 0; i < size; i++) {
      const item = it.next();
      if (item.done) throw new ValueError(`sequence too short (expected length ${size}, got ${i})`);
      if (i + 1 === size) check();
      yield item.value;
    }
  }
  if (Array.isArray(input)) {
    const nested =
      input.length > 0 && (Array.isArray(input[0]) || (ncols !== 1 && input.length === nrows));
    if (nested) {
      let i = 0;
      for (const row of sized(input, nrows)) {
        let j = 0;
        for (const value of sized(row, ncols)) yield [i, j++, convert(value)];
        i++;
      }
    } else {
      let i = 0;
      for (const value of sized(input, nrows * ncols)) {
        yield [Math.floor(i / ncols), i % ncols, convert(value)];
        i++;
      }
    }
    return;
  }
  const scalarZero =
    input == null ||
    input === false ||
    input === '' ||
    input === 0 ||
    input === 0n ||
    (input instanceof Rational && input.numerator === 0n) ||
    (typeof input === 'object' && 'value' in input && (input.value === 0n || input.value === 0));
  if (scalarZero) return;
  if (nrows !== ncols) throw new TypeError('nonzero scalar matrix must be square');
  const value = convert(input);
  for (let i = 0; i < nrows; i++) yield [i, i, value];
}
