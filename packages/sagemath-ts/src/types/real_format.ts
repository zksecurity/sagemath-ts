/** Sage RealNumber(53) repr and its Decimal-backed six-place formatting.
 * Reference: sage/rings/real_mpfr.pyx:str/__format__; MPFR 4.2.1 get_str/set_d.
 * @see Deviation: Gaussian Parameter Representation
 */
import { mpfr_set_d, mpfr_init2, mpfr_get_str } from '@sagemath-ts/mpfr-ts';
function digits53(value: number) {
  const x = mpfr_init2(53);
  mpfr_set_d(x, value);
  const [text, exponent] = mpfr_get_str(10, x.kind === 'zero' ? 16 : 15, x);
  return {
    kind: x.kind,
    sign: x.sign < 0 ? '-' : '',
    text: text.startsWith('-') ? text.slice(1) : text,
    exponent: x.kind === 'zero' ? 1 : exponent,
  };
}
export function repr53(value: number): string {
  const { kind, sign, text, exponent: e } = digits53(value);
  if (kind === 'nan') return 'NaN';
  if (kind === 'inf') return (sign || '+') + 'infinity';
  if (Math.abs(e - 1) >= 6) return sign + text[0] + '.' + text.slice(1) + 'e' + (e - 1);
  if (e <= 0) return sign + '0.' + '0'.repeat(-e) + text;
  if (e >= text.length) return sign + text + '0'.repeat(e - text.length) + '.';
  return sign + text.slice(0, e) + '.' + text.slice(e);
}
/** Sage RealNumber.__format__ -> Decimal(repr(self)) then decimal quantization. */
export function fixed53(value: number): string {
  const { kind, sign, text, exponent } = digits53(value);
  if (kind === 'nan') return 'NaN';
  if (kind === 'inf') return sign + 'Infinity';
  const places = 6,
    shift = exponent - text.length + places;
  let n = BigInt(text),
    q: bigint;
  if (shift >= 0) q = n * 10n ** BigInt(shift);
  else {
    const denominator = 10n ** BigInt(-shift);
    q = n / denominator;
    const r = n % denominator;
    if (2n * r > denominator || (2n * r === denominator && (q & 1n) !== 0n)) q++;
  }
  const result = q.toString().padStart(places + 1, '0');
  return sign + result.slice(0, -places) + '.' + result.slice(-places);
}
