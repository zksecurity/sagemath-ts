/** Linear-time coefficient packing for the native PARI/GMP polynomial schedules. */
export function trimPolynomial(x: readonly bigint[]): bigint[] {
  let n = x.length;
  while (n && x[n - 1] === 0n) n--;
  return x.slice(0, n);
}
export function packUnsigned(x: readonly bigint[], bits: number): bigint {
  if (!x.some((c) => c !== 0n)) return 0n;
  const hex = bits % 4 === 0,
    width = hex ? bits / 4 : bits,
    radix = hex ? 16 : 2;
  return BigInt(
    (hex ? '0x' : '0b') +
      x.slice().reverse().map((c) => c.toString(radix).padStart(width, '0')).join('')
  );
}
export function unpackUnsigned(x: bigint, bits: number, count: number): bigint[] {
  const hex = bits % 4 === 0,
    width = hex ? bits / 4 : bits,
    radix = hex ? 16 : 2;
  const text = x.toString(radix).padStart(width * count, '0');
  return Array.from({ length: count }, (_, i) =>
    BigInt((hex ? '0x' : '0b') + text.slice(text.length - (i + 1) * width, text.length - i * width))
  );
}
/** ZX_eval2BILspec builds positive and negative limb streams separately. */
export function packSigned(x: readonly bigint[], bits: number): bigint {
  return (
    packUnsigned(
      x.map((c) => (c > 0n ? c : 0n)),
      bits
    ) -
    packUnsigned(
      x.map((c) => (c < 0n ? -c : 0n)),
      bits
    )
  );
}
/** Z_mod2BIL_ZX: centered digits with a carry, applied to the product's magnitude. */
export function unpackSigned(x: bigint, bits: number, count: number): bigint[] {
  const negative = x < 0n,
    base = 1n << BigInt(bits),
    half = base >> 1n;
  let carry = 0n;
  return unpackUnsigned(negative ? -x : x, bits, count).map((c) => {
    const digit = c + carry;
    carry = digit >= half ? 1n : 0n;
    const value = carry ? digit - base : digit;
    return negative ? -value : value;
  });
}

/** Bound intermediates below the engine limit, as the portable FLINT/NTL kernels do. */
export const PACKED_PRODUCT_BIT_BUDGET = 1 << 16;
type Multiply = (a: readonly bigint[], b: readonly bigint[]) => bigint[];
type Square = (a: readonly bigint[]) => bigint[];
const coefficient = (x: bigint, p?: bigint): bigint => (p === undefined ? x : ((x % p) + p) % p);
function addCoefficients(a: readonly bigint[], b: readonly bigint[], p?: bigint): bigint[] {
  return trimPolynomial(
    Array.from({ length: Math.max(a.length, b.length) }, (_, i) =>
      coefficient((a[i] ?? 0n) + (b[i] ?? 0n), p)
    )
  );
}

/** Flx_mulspec's balanced/unbalanced Karatsuba split, also valid over Z. */
export function splitProduct(
  a: readonly bigint[],
  b: readonly bigint[],
  multiply: Multiply,
  p?: bigint
): bigint[] {
  if (a.length < b.length) [a, b] = [b, a];
  const cut = Math.ceil(a.length / 2),
    a0 = a.slice(0, cut),
    a1 = a.slice(cut);
  const out = new Array<bigint>(a.length + b.length - 1).fill(0n);
  const insert = (x: readonly bigint[], offset: number, sign = 1n) => {
    for (let i = 0; i < x.length; i++)
      out[i + offset] = coefficient(out[i + offset]! + sign * x[i]!, p);
  };
  if (b.length <= cut) {
    insert(multiply(a0, b), 0);
    insert(multiply(a1, b), cut);
  } else {
    const b0 = b.slice(0, cut),
      b1 = b.slice(cut),
      lo = multiply(a0, b0),
      hi = multiply(a1, b1);
    const middle = multiply(addCoefficients(a0, a1, p), addCoefficients(b0, b1, p));
    insert(lo, 0);
    insert(hi, 2 * cut);
    insert(middle, cut);
    insert(lo, cut, -1n);
    insert(hi, cut, -1n);
  }
  return trimPolynomial(out);
}

/** Flx_sqrspec's Karatsuba square, including its characteristic-two shortcut. */
export function splitSquare(a: readonly bigint[], square: Square, p?: bigint): bigint[] {
  const cut = Math.ceil(a.length / 2),
    a0 = a.slice(0, cut),
    a1 = a.slice(cut);
  const lo = square(a0),
    hi = square(a1),
    out = new Array<bigint>(2 * a.length - 1).fill(0n);
  const insert = (x: readonly bigint[], offset: number, sign = 1n) => {
    for (let i = 0; i < x.length; i++)
      out[i + offset] = coefficient(out[i + offset]! + sign * x[i]!, p);
  };
  insert(lo, 0);
  insert(hi, 2 * cut);
  if (p !== 2n) {
    insert(square(addCoefficients(a0, a1, p)), cut);
    insert(lo, cut, -1n);
    insert(hi, cut, -1n);
  }
  return trimPolynomial(out);
}
