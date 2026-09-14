/** Native GEN error display for extension coefficients and polynomials.
 * @see Deviation: PARI extension quotient adapters
 */
import type { ExtensionPolynomial as P } from './_extension_polynomial.js';
function displayNumeric(a: bigint[], variable: string): string {
  const terms: string[] = [];
  for (let i = a.length - 1; i >= 0; i--) {
    const c = a[i]!;
    if (!c) continue;
    const m = c < 0n ? -c : c;
    const v = i ? (m === 1n ? '' : m + '*') + variable + (i === 1 ? '' : '^' + i) : String(m);
    terms.push((terms.length ? (c < 0n ? ' - ' : ' + ') : c < 0n ? '-' : '') + v);
  }
  return terms.join('') || '0';
}
function displayOuter(a: P, mode: 0 | 1 | 2): string {
  const terms: string[] = [];
  for (let i = a.length - 1; i >= 0; i--) {
    let c = a[i]!,
      v: string,
      negative = false;
    if (mode === 0) {
      if (Array.isArray(c) && c.length < 2) c = c[0] ?? 0n;
      if (typeof c === 'bigint') {
        if (!c) continue;
        negative = c < 0n;
        const m = negative ? -c : c;
        v = i ? (m === 1n ? '' : m + '*') + 'x' + (i === 1 ? '' : '^' + i) : String(m);
      } else {
        let text = displayNumeric(c, 'y');
        if (c.filter((x) => x !== 0n).length === 1) {
          negative = text.startsWith('-');
          if (negative) text = text.slice(1);
        } else text = '(' + text + ')';
        v = text + (i ? '*x' + (i === 1 ? '' : '^' + i) : '');
      }
    } else {
      let words: bigint[];
      if (mode === 1) words = c as bigint[];
      else {
        let x = c as bigint;
        words = [];
        while (x) {
          words.push(x & ((1n << 64n) - 1n));
          x >>= 64n;
        }
      }
      v =
        'Vecsmall([' +
        [1n << 46n, ...words].map((x) => BigInt.asIntN(64, x)).join(', ') +
        '])' +
        (i ? '*x' + (i === 1 ? '' : '^' + i) : '');
    }
    terms.push((terms.length ? (negative ? ' - ' : ' + ') : negative ? '-' : '') + v);
  }
  return terms.join('') || '0';
}

import { pariErrorPayload as errorPayload } from './_error_display.js';
export function displayCoefficientError(a: bigint[]): string {
  return errorPayload(displayNumeric(a, 'y'));
}
export function displayExtensionError(a: P, mode: 0 | 1 | 2): string {
  return errorPayload(displayOuter(a, mode));
}
