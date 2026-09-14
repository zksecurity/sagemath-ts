import { _fmpz_poly_resultant_kernels } from '../fmpz_poly/gcd.js';
/** @see Deviation: Polynomial Integer Powers and Portable Native Products */
/** Dense-array ports of nmod_poly/mul.c and its Kronecker kernels.
 * @see Deviation: Native Modular Polynomial Products and Fraction Fields
 */
export function _nmod_poly_mul(a: readonly bigint[], b: readonly bigint[], p: bigint): bigint[] {
  const [A, B] = inputs(a, b, p);
  if (!B.length) return [];
  const bits = p.toString(2).length,
    cutoff = Math.min(A.length, 2 * B.length);
  // The portable FLINT dispatch (without architecture-specific fft_small).
  if (B.length <= 5 || 3 * cutoff < 2 * Math.max(bits, 10))
    return _nmod_poly_mul_classical(A, B, p);
  if (cutoff * bits < 800) return _nmod_poly_mul_KS(A, B, p);
  if (cutoff * (bits + 1) * (bits + 1) < 100000) return _nmod_poly_mul_KS2(A, B, p);
  return _nmod_poly_mul_KS4(A, B, p);
}
export function _nmod_poly_mul_classical(
  a: readonly bigint[],
  b: readonly bigint[],
  p: bigint
): bigint[] {
  const [A, B] = inputs(a, b, p);
  if (!B.length) return [];
  const out = Array<bigint>(A.length + B.length - 1).fill(0n);
  if (A === B) {
    for (let i = 0; i < A.length; i++) {
      out[2 * i]! += A[i]! * A[i]!;
      const twice = 2n * A[i]!;
      for (let j = i + 1; j < A.length; j++) out[i + j]! += twice * A[j]!;
    }
  } else {
    for (let i = 0; i < A.length; i++)
      for (let j = 0; j < B.length; j++) out[i + j]! += A[i]! * B[j]!;
  }
  return normalized(out, p);
}
export function _nmod_poly_mul_KS(a: readonly bigint[], b: readonly bigint[], p: bigint): bigint[] {
  const [A, B] = inputs(a, b, p);
  if (!B.length) return [];
  const bits = 2 * p.toString(2).length + B.length.toString(2).length;
  if ((A.length + B.length - 1) * bits > 1 << 19)
    return normalized(_fmpz_poly_resultant_kernels.multiply(A, B), p);
  const x = pack(A, bits),
    y = A === B ? x : pack(B, bits);
  return normalized(digits(x * y, A.length + B.length - 1, bits), p);
}
export function _nmod_poly_mul_KS2(
  a: readonly bigint[],
  b: readonly bigint[],
  p: bigint
): bigint[] {
  const [A, B] = inputs(a, b, p);
  if (!B.length) return [];
  if (B.length === 1)
    return normalized(
      A.map((c) => c * B[0]!),
      p
    );
  const bits = Math.ceil((2 * p.toString(2).length + ceilLog2(B.length)) / 2);
  if ((A.length + B.length - 1) * bits > 1 << 19)
    return normalized(_fmpz_poly_resultant_kernels.multiply(A, B), p);
  const [even, odd] = productParities(A, B, bits),
    n = A.length + B.length - 1;
  return interleave(
    digits(even, Math.ceil(n / 2), 2 * bits),
    digits(odd, Math.floor(n / 2), 2 * bits),
    p
  );
}
export function _nmod_poly_mul_KS4(
  a: readonly bigint[],
  b: readonly bigint[],
  p: bigint
): bigint[] {
  const [A, B] = inputs(a, b, p);
  if (!B.length) return [];
  if (B.length === 1)
    return normalized(
      A.map((c) => c * B[0]!),
      p
    );
  const bits = Math.ceil((2 * p.toString(2).length + ceilLog2(B.length)) / 4),
    n = A.length + B.length - 1;
  if (n * bits > 1 << 19) return normalized(_fmpz_poly_resultant_kernels.multiply(A, B), p);
  const [even, odd] = productParities(A, B, bits);
  const ar = A.slice().reverse(),
    br = A === B ? ar : B.slice().reverse();
  const [reverseEven, reverseOdd] = productParities(ar, br, bits);
  // Reversing an even-length product interchanges its even and odd parts.
  const re = n % 2 ? reverseEven : reverseOdd,
    ro = n % 2 ? reverseOdd : reverseEven;
  return interleave(
    recover(even, re, Math.ceil(n / 2), 2 * bits),
    recover(odd, ro, Math.floor(n / 2), 2 * bits),
    p
  );
}
function normalized(a: readonly bigint[], p: bigint): bigint[] {
  const out = a.map((c) => {
    const r = c % p;
    return r < 0n ? r + p : r;
  });
  while (out.length && out[out.length - 1] === 0n) out.pop();
  return out;
}
function inputs(a: readonly bigint[], b: readonly bigint[], p: bigint): [bigint[], bigint[]] {
  if (p < 1n || p >= 1n << 64n) throw new RangeError('modulus must fit a positive unsigned word');
  const A = normalized(a, p),
    B = a === b ? A : normalized(b, p);
  return A.length >= B.length ? [A, B] : [B, A];
}
function ceilLog2(n: number): number {
  return n <= 1 ? 0 : (n - 1).toString(2).length;
}
/** String packing/unpacking is linear in packed size; no repeated growing BigInt shifts. */
function pack(a: readonly bigint[], bits: number): bigint {
  if (!a.length) return 0n;
  return BigInt(
    '0b' +
      a
        .map((c) => c.toString(2).padStart(bits, '0'))
        .reverse()
        .join('')
  );
}
function digits(a: bigint, n: number, bits: number): bigint[] {
  const s = a.toString(2).padStart(n * bits, '0');
  return Array.from({ length: n }, (_, i) => {
    const end = s.length - i * bits;
    return BigInt('0b' + s.slice(end - bits, end));
  });
}
/** h(B)=f(B)g(B), h(-B)=f(-B)g(-B), followed by even/odd interpolation. */
function productParities(
  a: readonly bigint[],
  b: readonly bigint[],
  bits: number
): [bigint, bigint] {
  const evaluate = (v: readonly bigint[]): [bigint, bigint] => {
    const e = pack(
        v.filter((_, i) => i % 2 === 0),
        2 * bits
      ),
      o =
        pack(
          v.filter((_, i) => i % 2 === 1),
          2 * bits
        ) << BigInt(bits);
    return [e + o, e - o];
  };
  const [ap, am] = evaluate(a),
    [bp, bm] = a === b ? [ap, am] : evaluate(b);
  const plus = ap * bp,
    minus = am * bm;
  return [(plus + minus) >> 1n, (plus - minus) >> BigInt(bits + 1)];
}
/** KS2_reduce.c: two-digit windows reconstruct overlapping coefficients from both ends. */
function recover(normal: bigint, reverse: bigint, n: number, bits: number): bigint[] {
  if (!n) return [];
  const x = digits(normal, n + 1, bits),
    y = digits(reverse, n + 1, bits);
  const shift = BigInt(bits),
    mask = (1n << shift) - 1n,
    out: bigint[] = [];
  let x0 = x[0]!,
    y1 = y[n]!,
    borrow = 0n;
  for (let i = 0; i < n; i++) {
    const y0 = y[n - 1 - i]!;
    let x1 = x[i + 1]!;
    if (y0 < x0) y1--;
    out.push(x0 + (y1 << shift));
    y1 += borrow;
    borrow = x1 < y1 ? 1n : 0n;
    x1 -= y1;
    y1 = (y0 - x0) & mask;
    x0 = x1 & mask;
  }
  return out;
}
function interleave(even: readonly bigint[], odd: readonly bigint[], p: bigint): bigint[] {
  const out = Array<bigint>(even.length + odd.length);
  for (let i = 0; i < even.length; i++) out[2 * i] = even[i]!;
  for (let i = 0; i < odd.length; i++) out[2 * i + 1] = odd[i]!;
  return normalized(out, p);
}
