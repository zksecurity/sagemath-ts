/** PARI Dedekind sums; reference/pari/src/basemath/elltrans.c. */

function normalize(numerator: bigint, denominator: bigint): [bigint, bigint] {
  if (denominator < 0n) {
    numerator = -numerator;
    denominator = -denominator;
  }
  let a = numerator < 0n ? -numerator : numerator;
  let b = denominator;
  while (b !== 0n) [a, b] = [b, a % b];
  return [numerator / a, denominator / a];
}

/** Knuth's quotient recurrence from sumdedekind_coprime/u_sumdedekind_coprime. */
function coprimeSum(h: bigint, k: bigint): [bigint, bigint] {
  const absK = k < 0n ? -k : k;
  h = ((h % absK) + absK) % absK;
  // Keep the native word/generic split: for negative k their conventions differ.
  if (absK <= (2n * ((1n << 63n) - 1n)) / 3n && k < 0n) {
    k = absK;
    h = h === 0n ? 0n : k - h;
  }
  let sign = 1n;
  let s1 = 0n;
  let s2 = h;
  let p = 1n;
  let pp = 0n;
  while (h !== 0n) {
    const nextH = k % h;
    const quotient = k / h;
    if (h === 1n || h === -1n) s2 += sign * p;
    s1 += sign * quotient;
    sign = -sign;
    k = h;
    h = nextH;
    [p, pp] = [quotient * p + pp, p];
  }
  if (sign === -1n) s1 -= 3n;
  return normalize(s2 + p * s1, 12n * p);
}

/**
 * Return the reduced numerator/denominator from PARI's sumdedekind(h,k).
 * @see Deviation: Dedekind sum backend arithmetic
 */
export function sumdedekind(h: bigint, k: bigint): [bigint, bigint] {
  if (k === 0n) {
    const operation = h === 0n ? 'diviiexact' : 'dvmdii';
    const error = new Error(`impossible inverse in ${operation}: 0`);
    error.name = 'PariError';
    throw error;
  }
  let a = h < 0n ? -h : h;
  let b = k < 0n ? -k : k;
  while (b !== 0n) [a, b] = [b, a % b];
  return coprimeSum(h / a, k / a);
}
