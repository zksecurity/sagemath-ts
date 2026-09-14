/** Python 3.12 float conversion for Sage numeric inputs and text/byte buffers.
 * References: CPython v3.12.5 Objects/floatobject.c, Python/pystrtod.c,
 * Objects/unicodeobject.c; UnicodeData 15.0.0 (Nd zero code points below).
 * @see Deviation: Binary Density Input Conversion
 */
import { OverflowError, ValueError } from '../errors.js';
import { Integer } from '../rings/integer_ring.js';
import { Rational } from '../rings/rational.js';
import { RealNumber } from '../rings/real_mpfr.js';
import { IntegerMod } from '../rings/finite_rings/integer_mod.js';
import { PrimeFieldElement } from '../rings/finite_rings/finite_field_extension.js';
import { GF2Element } from '../rings/finite_rings/gf2.js';
export type FloatInput =
  | number
  | bigint
  | boolean
  | string
  | Uint8Array
  | null
  | Integer
  | Rational
  | RealNumber
  | IntegerMod
  | PrimeFieldElement
  | GF2Element;
import { decimalZeros, printableRanges } from './python_float_data.js';
function printable(cp: number): boolean {
  let lo = 0,
    hi = printableRanges.length;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    const [a, b] = printableRanges[mid]!;
    if (cp < a) hi = mid;
    else if (cp > b) lo = mid + 1;
    else return true;
  }
  return false;
}

function repr(value: string, bytes: boolean): string {
  const quote = value.includes("'") && !value.includes('"') ? '"' : "'";
  let output = bytes ? 'b' + quote : quote;
  for (const c of value) {
    const n = c.codePointAt(0)!;
    if (c === quote || c === '\\') output += '\\' + c;
    else if (c === '\n') output += '\\n';
    else if (c === '\r') output += '\\r';
    else if (c === '\t') output += '\\t';
    else if ((bytes && (n < 32 || n >= 127)) || (!bytes && !printable(n)))
      output +=
        n <= 255
          ? '\\x' + n.toString(16).padStart(2, '0')
          : n <= 65535
            ? '\\u' + n.toString(16).padStart(4, '0')
            : '\\U' + n.toString(16).padStart(8, '0');
    else output += c;
  }
  return output + quote;
}

export function float(value: FloatInput): number {
  if (typeof value === 'number') return value;
  if (typeof value === 'boolean') return value ? 1 : 0;
  if (typeof value === 'bigint') {
    const result = Number(value);
    if (!Number.isFinite(result)) throw new OverflowError('int too large to convert to float');
    return result;
  }
  // Sage integers/rationals permit infinity; Python built-in integers reject it.
  if (value instanceof Integer) return Number(value.value);
  if (value instanceof Rational || value instanceof RealNumber) return value.toNumber();
  if (
    value instanceof IntegerMod ||
    value instanceof PrimeFieldElement ||
    value instanceof GF2Element
  )
    return Number(value.value);
  const bytes = value instanceof Uint8Array;
  if (typeof value === 'string' || bytes) {
    const original = bytes
      ? Array.from(value as Uint8Array, (x) => String.fromCharCode(x)).join('')
      : (value as string);
    let text = original;
    if (!bytes)
      text = Array.from(text, (c) => {
        const n = c.codePointAt(0)!;
        if (n < 128) return c;
        if (
          n === 133 ||
          n === 160 ||
          n === 5760 ||
          (n >= 8192 && n <= 8202) ||
          n === 8232 ||
          n === 8233 ||
          n === 8239 ||
          n === 8287 ||
          n === 12288
        )
          return ' ';
        const zero = decimalZeros.find((z) => n >= z && n < z + 10);
        return zero === undefined ? c : String(n - zero);
      }).join('');
    const space = (c: number) => c === 32 || (c >= 9 && c <= 13);
    let first = 0,
      last = text.length;
    while (first < last && space(text.charCodeAt(first))) first++;
    while (last > first && space(text.charCodeAt(last - 1))) last--;
    text = text.slice(first, last);
    const digits = '[0-9](?:_?[0-9])*';
    const numeric = new RegExp(
      `^[+-]?(?:(?:${digits}(?:\\.(?:${digits})?)?|\\.${digits})(?:[eE][+-]?${digits})?|inf(?:inity)?|nan)$`,
      'i'
    );
    if (!numeric.test(text))
      throw new ValueError('could not convert string to float: ' + repr(original, bytes));
    const clean = text.replaceAll('_', '');
    if (/^[+-]?nan$/i.test(clean))
      return clean.startsWith('-') ? -Number(clean.slice(1)) : Number(clean);
    if (/^[+-]?inf(?:inity)?$/i.test(clean)) return clean.startsWith('-') ? -Infinity : Infinity;
    return Number(clean);
  }
  const name = value == null ? 'NoneType' : Array.isArray(value) ? 'list' : 'dict';
  throw new TypeError(`float() argument must be a string or a real number, not '${name}'`);
}
