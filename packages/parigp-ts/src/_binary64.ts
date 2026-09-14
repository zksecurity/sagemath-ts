/** Binary64 primitives used by PARI reduction stages. */
const view = new DataView(new ArrayBuffer(8));
function parts(x: number): [bigint, number, bigint] {
  view.setFloat64(0, x, false);
  const raw = view.getBigUint64(0, false),
    sign = raw & (1n << 63n),
    e = Number((raw >> 52n) & 2047n),
    fraction = raw & ((1n << 52n) - 1n);
  return [e ? fraction | (1n << 52n) : fraction, e ? e - 1075 : -1074, sign];
}
/** C ldexp, including a single correctly rounded subnormal result. */
export function ldexp(x: number, e: number): number {
  if (!x || !Number.isFinite(x)) return x;
  const [m, p, sign] = parts(x),
    bits = m.toString(2).length,
    exponent = p + e,
    top = exponent + bits - 1;
  if (top > 1023) return sign ? -Infinity : Infinity;
  let raw: bigint;
  if (top >= -1022)
    raw = (BigInt(top + 1023) << 52n) | ((m << BigInt(53 - bits)) & ((1n << 52n) - 1n));
  else {
    const shift = exponent + 1074;
    if (shift >= 0) raw = m << BigInt(shift);
    else if (-shift > bits) raw = 0n;
    else {
      const k = BigInt(-shift),
        q = m >> k,
        remainder = m - (q << k),
        half = 1n << (k - 1n);
      raw = q + (remainder > half || (remainder === half && (q & 1n) !== 0n) ? 1n : 0n);
    }
  }
  view.setBigUint64(0, raw | sign, false);
  return view.getFloat64(0, false);
}
/** Correctly rounded binary64 a*b+c, including subnormals. */
export function fma(a: number, b: number, c: number): number {
  if (!Number.isFinite(a) || !Number.isFinite(b) || Number.isNaN(c)) return a * b + c;
  if (!Number.isFinite(c)) return c;
  const [am, ae, as] = parts(a),
    [bm, be, bs] = parts(b),
    [cm, ce, cs] = parts(c);
  const e = Math.min(ae + be, ce);
  const exact =
    ((as === bs ? am * bm : -am * bm) << BigInt(ae + be - e)) + ((cs ? -cm : cm) << BigInt(ce - e));
  if (exact === 0n) return a * b + c;
  const negative = exact < 0n,
    m = negative ? -exact : exact,
    bits = m.toString(2).length;
  const unit = Math.max(e + bits - 53, -1074),
    shift = unit - e;
  let rounded: bigint;
  if (shift <= 0) rounded = m << BigInt(-shift);
  else {
    const k = BigInt(shift),
      q = m >> k,
      remainder = m - (q << k),
      half = 1n << (k - 1n);
    rounded = q + (remainder > half || (remainder === half && (q & 1n) !== 0n) ? 1n : 0n);
  }
  return ldexp(negative ? -Number(rounded) : Number(rounded), unit);
}
export function frexp(x: number): [number, number] {
  if (!x || !Number.isFinite(x)) return [x, 0];
  const [m, p, sign] = parts(x),
    bits = m.toString(2).length;
  return [((sign ? -1 : 1) * Number(m)) / 2 ** bits, p + bits];
}
