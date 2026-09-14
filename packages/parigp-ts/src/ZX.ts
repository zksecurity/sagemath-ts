/** PARI ZX.c multiplication and squaring for 64-bit GMP tuning.
 * @see Deviation: PARI polynomial multiplication adapters
 */
import {
  PACKED_PRODUCT_BIT_BUDGET,
  splitProduct,
  splitSquare,
  packSigned,
  unpackSigned,
  trimPolynomial,
} from './_polynomial_packing.js';
const SQUARE_LOW = [
  17, 32, 96, 112, 160, 128, 128, 160, 160, 160, 160, 160, 176, 192, 192, 192, 192, 192, 224, 224,
  224, 240, 240, 240, 272, 288, 288, 240, 288, 304, 304, 304, 304, 304, 304, 352, 352, 368, 352,
  352, 352, 368, 368, 432, 432, 496, 432, 496, 496,
];
const SQUARE_HIGH = [
  102860, 70254, 52783, 27086, 24623, 18500, 15289, 13899, 12635, 11487, 10442, 9493, 8630, 7845,
  7132, 7132, 6484, 6484, 5894, 5894, 4428, 4428, 3660, 4428, 3660, 3660, 2749, 2499, 2272, 2066,
  1282, 1282, 1166, 1166, 1166, 1166, 1166, 1166, 1166, 963, 963, 724, 658, 658, 658, 528, 528, 528,
  528,
];
const exponent = (x: readonly bigint[]): number =>
  x.reduce((m, c) => Math.max(m, (c < 0n ? -c : c).toString(2).length - 1), 0);

/** ZX_mulspec: scalar and three-product linear kernels, then signed Kronecker packing. */
export function ZX_mul(x: readonly bigint[], y: readonly bigint[]): bigint[] {
  if (x === y) return ZX_sqr(x);
  let a = trimPolynomial(x),
    b = trimPolynomial(y);
  if (!a.length || !b.length) return [];
  let va = 0,
    vb = 0;
  while (a[va] === 0n) va++;
  while (b[vb] === 0n) vb++;
  a = a.slice(va);
  b = b.slice(vb);
  const v = va + vb;
  let result: bigint[];
  if (a.length === 1) result = b.map((c) => c * a[0]!);
  else if (b.length === 1) result = a.map((c) => c * b[0]!);
  else if (a.length === 2 && b.length === 2) {
    const lo = a[0]! * b[0]!,
      hi = a[1]! * b[1]!;
    result = [lo, lo + hi - (a[1]! - a[0]!) * (b[1]! - b[0]!), hi];
  } else {
    const e = exponent(a) + exponent(b) + Math.floor(Math.log2(Math.min(a.length, b.length))) + 3;
    const bits = (Math.floor(e / 64) + 1) * 64;
    result =
      (a.length + b.length - 1) * bits > PACKED_PRODUCT_BIT_BUDGET
        ? splitProduct(a, b, ZX_mul)
        : unpackSigned(packSigned(a, bits) * packSigned(b, bits), bits, a.length + b.length - 1);
  }
  return trimPolynomial([...new Array<bigint>(v).fill(0n), ...result]);
}

/** ZX_sqrspec: native exponent/degree dispatch and centered coefficient recovery. */
export function ZX_sqr(x: readonly bigint[]): bigint[] {
  let a = trimPolynomial(x);
  if (!a.length) return [];
  let v = 0;
  while (a[v] === 0n) v++;
  a = a.slice(v);
  let result: bigint[];
  if (a.length === 1) result = [a[0]! * a[0]!];
  else {
    const n = a.length,
      e = exponent(a),
      index = n - 2;
    if (index < SQUARE_LOW.length && SQUARE_LOW[index]! <= e && e <= SQUARE_HIGH[index]!) {
      result = new Array<bigint>(2 * n - 1).fill(0n);
      for (let i = 0; i < n; i++) {
        result[2 * i]! += a[i]! * a[i]!;
        for (let j = 0; j < i; j++) result[i + j]! += 2n * a[i]! * a[j]!;
      }
    } else {
      const bits = (Math.floor((2 * e + Math.floor(Math.log2(n)) + 3) / 64) + 1) * 64;
      if ((2 * n - 1) * bits > PACKED_PRODUCT_BIT_BUDGET) result = splitSquare(a, ZX_sqr);
      else {
        const z = packSigned(a, bits);
        result = unpackSigned(z * z, bits, 2 * n - 1);
      }
    }
  }
  return trimPolynomial([...new Array<bigint>(2 * v).fill(0n), ...result]);
}

export function ZX_deriv(f: bigint[]): bigint[] {
  const r: bigint[] = [];
  for (let i = 1; i < f.length; i++) r.push(f[i]! * BigInt(i));
  return trimPolynomial(r);
}

/** Native ZX_rem quotient-first elimination; divisor must be monic.
 * @see Deviation: PARI rational trace and norm adapters
 */
export function ZX_rem(input: readonly bigint[], divisor: readonly bigint[]): bigint[] {
  const x = trimPolynomial(input),
    y = trimPolynomial(divisor),
    dx = x.length - 1,
    dy = y.length - 1;
  if (y.at(-1) !== 1n) throw new RangeError('ZX_rem requires a monic polynomial');
  if (dx < dy) return x;
  if (dy === 0) return [];
  const dz = dx - dy,
    q = Array<bigint>(dz + 1).fill(0n);
  q[dz] = x[dx]!;
  for (let i = dx - 1; i >= dy; i--) {
    let v = x[i]!;
    for (let j = i - dy + 1; j <= i && j <= dz; j++) v -= q[j]! * y[i - j]!;
    q[i - dy] = v;
  }
  const r = Array<bigint>(dy).fill(0n);
  for (let i = dy - 1; i >= 0; i--) {
    let v = x[i]!;
    for (let j = 0; j <= i && j <= dz; j++) v -= q[j]! * y[i - j]!;
    r[i] = v;
  }
  return trimPolynomial(r);
}
import {
  type RationalPolynomialData,
  primitivePolynomial,
  rationalProduct,
  scaledPolynomial,
} from './_rational_polynomial.js';
/** Native QX_mul: remove rational content, multiply in ZZ[X], restore content.
 * @see Deviation: PARI rational trace and norm adapters
 */
export function QX_mul(
  x: RationalPolynomialData,
  y: RationalPolynomialData
): RationalPolynomialData {
  const [a, ca] = primitivePolynomial(x),
    [b, cb] = primitivePolynomial(y);
  return scaledPolynomial(ZX_mul(a, b), rationalProduct(ca, cb));
}
/** Native QX_ZX_rem, for a monic integral divisor.
 * @see Deviation: PARI rational trace and norm adapters
 */
export function QX_ZX_rem(x: RationalPolynomialData, T: bigint[]): RationalPolynomialData {
  const [a, c] = primitivePolynomial(x);
  return scaledPolynomial(ZX_rem(a, T), c);
}

/** Native ZX_Z_eval skips runs of zero coefficients using integer powers. */
export function ZX_Z_eval(x: readonly bigint[], y: bigint): bigint {
  let i = x.length - 1;
  while (i >= 0 && x[i] === 0n) i--;
  if (i <= 0) return x[0] ?? 0n;
  if (!y) return x[0]!;
  let t = x[i--]!;
  while (i >= 0) {
    let j = i;
    while (!x[j]) {
      if (j === 0) return t * (i === j ? y : y ** BigInt(i - j + 1));
      j--;
    }
    const r = i === j ? y : y ** BigInt(i - j + 1);
    t = t * r + x[j]!;
    i = j - 1;
  }
  return t;
}
