/** matrix0.pyx string-rendering adapter; see DESIGN.md.
 * @see Deviation: Binary Matrix Formatting and Subdivisions
 */
import { AttributeError, IndexError, OverflowError } from '../errors.js';
export interface MatrixStringData {
  nrows: number;
  ncols: number;
  entries: () => readonly number[];
  rowDivisions: readonly bigint[];
  colDivisions: readonly bigint[];
  binaryFast?: boolean;
}
export type MatrixRepresentation = Record<number, string> | ((value: number) => string);
export function str(
  data: MatrixStringData,
  mapping?: MatrixRepresentation | null,
  zero?: string | null,
  one?: string | null
): string {
  if (!data.nrows || !data.ncols) return '[]';
  if (data.binaryFast) {
    const columns = new Set(data.colDivisions.map(Number));
    const entries = data.entries(),
      output: string[] = [];
    let divider = '[';
    for (let c = 0; c < data.ncols; c++) divider += (c ? (columns.has(c) ? '+' : '-') : '') + '-';
    divider += ']';
    for (let r = 0; r < data.nrows; r++) {
      let row = '[';
      for (let c = 0; c < data.ncols; c++)
        row += (c ? (columns.has(c) ? '|' : ' ') : '') + entries[r * data.ncols + c];
      output.push(row + ']');
    }
    // The binary fast path uses Python list.insert, including its clamped negative
    // and oversized row positions. Generic mapped rendering validates count indices.
    for (const line of [...data.rowDivisions].reverse()) {
      if (line < -(1n << 63n) || line >= 1n << 63n)
        throw new OverflowError('Python int too large to convert to C ssize_t');
      const index = Number(line);
      output.splice(
        index < 0 ? Math.max(0, output.length + index) : Math.min(index, output.length),
        0,
        divider
      );
    }
    return output.join('\n');
  }
  const map = mapping ?? {};
  const dictionary = typeof map === 'object' && !Array.isArray(map);
  const length = (value: unknown): number => {
    if (typeof value === 'string') return Array.from(value).length;
    if (Array.isArray(value)) return value.length;
    if (value && Object.getPrototypeOf(value) === Object.prototype)
      return Object.keys(value).length;
    const name =
      value == null
        ? 'NoneType'
        : typeof value === 'boolean'
          ? 'bool'
          : typeof value === 'number'
            ? Number.isInteger(value)
              ? 'int'
              : 'float'
            : typeof value === 'bigint'
              ? 'int'
              : 'object';
    throw new TypeError(`object of type '${name}' has no len()`);
  };
  const lines = (values: readonly bigint[], size: number) => {
    const counts = Array<number>(size + 1).fill(0);
    const indices = values.map((value) => {
      if (value < -(1n << 63n) || value >= 1n << 63n)
        throw new OverflowError('Python int too large to convert to C ssize_t');
      const index = Number(value),
        wrapped = index < 0 ? index + size + 1 : index;
      if (wrapped < 0 || wrapped > size) throw new IndexError('list index out of range');
      counts[wrapped]++;
      return index;
    });
    return { counts, indices };
  };
  const rowDivs = lines(data.rowDivisions, data.nrows),
    colDivs = lines(data.colDivisions, data.ncols);
  if (dictionary) {
    if (zero !== undefined && zero !== null) (map as Record<number, string>)[0] = zero;
    if (one !== undefined && one !== null) (map as Record<number, string>)[1] = one;
  }
  const entries = data.entries().map((value) => {
    if (typeof map === 'function') return map(value);
    if (dictionary)
      return Object.hasOwn(map, value) ? (map as Record<number, string>)[value]! : String(value);
    if (Array.isArray(map)) {
      if (map.includes(value)) throw new AttributeError("'list' object has no attribute 'get'");
      return String(value);
    }
    if (!map) return String(value);
    if (typeof map === 'string')
      throw new TypeError(
        "'in <string>' requires string as left operand, not sage.rings.finite_rings.integer_mod.IntegerMod_int"
      );
    throw new TypeError(
      `argument of type '${typeof map === 'boolean' ? 'bool' : typeof map === 'number' && !Number.isInteger(map) ? 'float' : 'int'}' is not iterable`
    );
  });
  let width = 0;
  for (const entry of entries) width = Math.max(width, length(entry));
  const boundaries = [0, ...colDivs.indices, data.ncols];
  const hline = boundaries
    .slice(0, -1)
    .map((a, i) => '-'.repeat(Math.max(0, (width + 1) * (boundaries[i + 1]! - a) - 1)))
    .join('+');
  const output: string[] = [];
  for (let r = 0; r < data.nrows; r++) {
    for (const div of rowDivs.indices) if (div === r) output.push('[' + hline + ']');
    let row = '';
    for (let c = 0; c < data.ncols; c++) {
      const count = colDivs.counts[c]!;
      const separator = count ? '|'.repeat(count) : c === 0 ? '' : ' ';
      const entry = entries[r * data.ncols + c]!;
      if (typeof entry !== 'string')
        throw new TypeError(
          `can only concatenate str (not "${Array.isArray(entry) ? 'list' : 'dict'}") to str`
        );
      row += separator + ' '.repeat(width - length(entry)) + entry;
    }
    row += '|'.repeat(colDivs.counts[data.ncols]!);
    output.push('[' + row + ']');
  }
  for (const div of rowDivs.indices) if (div === data.nrows) output.push('[' + hline + ']');
  return output.join('\n');
}
