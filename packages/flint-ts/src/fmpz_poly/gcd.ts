/** @see Deviation: Polynomial Integer Powers and Portable Native Products */
import { _fmpz_poly_divrem } from './divrem.js';
import { _nmod_poly_gcd } from '../nmod_poly/gcd.js';

/** Dense bigint array port of FLINT fmpz_poly/gcd*.c.
 * @see Deviation: Dense FLINT polynomial GCD kernels
 */
export function _fmpz_poly_gcd(a: readonly bigint[], b: readonly bigint[]): bigint[] {
  const [A, B] = ordered(a, b);
  if (!B.length) return A.length && A[A.length - 1]! < 0n ? A.map((c) => -c) : A;
  // gcd.c:24-47 removes powers of x before choosing the integer GCD kernel.
  let va = 0,
    vb = 0;
  while (va < A.length - 1 && A[va] === 0n) va++;
  while (vb < B.length - 1 && B[vb] === 0n) vb++;
  if (va || vb)
    return Array<bigint>(Math.min(va, vb))
      .fill(0n)
      .concat(_fmpz_poly_gcd(A.slice(va), B.slice(vb)));
  if (A.length < 6) return _fmpz_poly_gcd_subresultant(A, B);
  if (maxBits(A) + maxBits(B) < 128) {
    const result = _fmpz_poly_gcd_heuristic(A, B);
    if (result !== null) return result;
  }
  return _fmpz_poly_gcd_modular(A, B);
}
export function _fmpz_poly_gcd_subresultant(a: readonly bigint[], b: readonly bigint[]): bigint[] {
  let [A, B] = ordered(a, b);
  if (!B.length) return A.length && A[A.length - 1]! < 0n ? A.map((c) => -c) : A;
  if (B.length === 1) return [integerGcd(content(A), B[0]!)];
  const ca = content(A),
    cb = content(B),
    d = integerGcd(ca, cb);
  A = A.map((c) => c / ca);
  B = B.map((c) => c / cb);
  let g = 1n,
    h = 1n;
  // gcd_subresultant.c:59-103 is the subresultant PRS, not the primitive
  // pseudo-remainder loop formerly used in the Sage wrapper.
  for (;;) {
    const delta = A.length - B.length;
    const R = pseudoRemainder(A, B);
    if (R.length === 1) return [d];
    if (R.length === 0) return primitive(B).map((c) => c * d);
    A = B;
    const hPower = h ** BigInt(delta);
    B = R.map((c) => c / (g * hPower));
    const lead = A[A.length - 1]!;
    if (delta === 1) h = lead;
    else h = (h * lead ** BigInt(delta)) / hPower;
    g = lead;
  }
}
export function _fmpz_poly_gcd_heuristic(
  a: readonly bigint[],
  b: readonly bigint[]
): bigint[] | null {
  let [A, B] = ordered(a, b);
  if (!B.length) return A.length && A[A.length - 1]! < 0n ? A.map((c) => -c) : A;
  const ca = content(A),
    cb = content(B),
    d = integerGcd(ca, cb);
  if (B.length === 1) return [d];
  A = A.map((c) => c / ca);
  B = B.map((c) => c / cb);
  if (B.length === 2) return divides(A, B) ? primitive(B).map((c) => c * d) : [d];
  const ba = maxBits(A),
    bb = maxBits(B);
  let width = Math.max(Math.min(ba, bb) + 6, Math.max(ba, bb) + 1);
  if (width >= 32) width = 64 * Math.ceil(width / 64);
  const signA = A[A.length - 1]! < 0n ? -1n : 1n;
  const signB = B[B.length - 1]! < 0n ? -1n : 1n;
  const packedA = signA * pack(A, width),
    packedB = signB * pack(B, width);
  let packedG = integerGcd(packedA, packedB);
  const limbs = (x: bigint) => Math.ceil(bits(x) / 64);
  const glen = Math.min(Math.floor((limbs(packedG) * 64) / width) + 1, B.length);
  let G = unpack(packedG, width, glen);
  const gc = content(G);
  packedG /= gc;
  G = G.map((c) => c / gc);
  const verify = (P: bigint[], packedP: bigint): boolean => {
    if (packedP % packedG !== 0n) return false;
    const qlimbs = limbs(packedP) - limbs(packedG) + 1;
    const qlen = Math.min(P.length, Math.floor((qlimbs * 64) / width) + 1);
    const Q = unpack(packedP / packedG, width, qlen);
    if (maxBits(G) + maxBits(Q) + Math.min(bits(BigInt(G.length)), bits(BigInt(Q.length))) < width)
      return true;
    // gcd_heuristic.c uses sign1 for both multiply-out checks.
    const product = multiply(Q, G).map((c) => c * signA);
    return product.length === P.length && product.every((c, i) => c === P[i]);
  };
  return verify(A, packedA) && verify(B, packedB) ? G.map((c) => c * d) : null;
}
export function _fmpz_poly_gcd_modular(a: readonly bigint[], b: readonly bigint[]): bigint[] {
  let [A, B] = ordered(a, b);
  if (!B.length) return A.length && A[A.length - 1]! < 0n ? A.map((c) => -c) : A;
  const ca = content(A),
    cb = content(B),
    d = integerGcd(ca, cb);
  if (B.length === 1) return [d];
  A = A.map((c) => c / ca);
  B = B.map((c) => c / cb);
  const normBits = (P: bigint[]) =>
    A.length < 64 && B.length < 64
      ? bits(P.reduce((s, c) => s + c * c, 0n))
      : 2 * maxBits(P) + bits(BigInt(P.length));
  const g = integerGcd(A[A.length - 1]!, B[B.length - 1]!);
  const leadProduct = A[A.length - 1]! * B[B.length - 1]!;
  const evaluateMinusOne = (P: bigint[]) => P.reduce((s, c, i) => s + (i % 2 ? c : -c), 0n);
  const small = Math.max(bits(integerGcd(evaluateMinusOne(A), evaluateMinusOne(B))), bits(g), 2);
  const bound = (A.length + 2) * Math.max(normBits(A), normBits(B)) + A.length;
  let modulus = 0n,
    result: bigint[] = [],
    n = B.length,
    unlucky = 0,
    currentBits = 0;
  for (const p of modularPrimes()) {
    if (leadProduct % p === 0n) {
      unlucky += 63;
      continue;
    }
    let h = _nmod_poly_gcd(A, B, p);
    if (h.length === 1) return [d];
    if (h.length > n + 1) {
      unlucky += 63;
      continue;
    }
    const scale = mod(g * inverse(h[h.length - 1]!, p), p);
    h = h.map((c) => mod(c * scale, p));
    if (h.length <= n) {
      unlucky += bits(modulus);
      result = h.map((c) => (c > p / 2n ? c - p : c));
      if (g === 1n) {
        if (divides(B, result) && divides(A, result)) return result.map((c) => c * d);
      } else {
        if (63 + unlucky >= bound) return primitive(result).map((c) => c * d);
        if (63 >= small) {
          const candidate = primitive(result);
          if (divides(B, candidate) && divides(A, candidate)) return candidate.map((c) => c * d);
        }
      }
      currentBits = maxBits(result);
      modulus = p;
      n = h.length - 1;
      continue;
    }
    const inv = inverse(modulus, p),
      nextModulus = modulus * p;
    result = result.map((c, i) => {
      let v = c + mod((h[i]! - c) * inv, p) * modulus;
      v = mod(v, nextModulus);
      return v > nextModulus / 2n ? v - nextModulus : v;
    });
    modulus = nextModulus;
    const newBits = maxBits(result);
    if (newBits === currentBits || bits(modulus) >= small) {
      const candidate = g === 1n ? result : primitive(result);
      if (bits(modulus) + unlucky >= bound || (divides(B, candidate) && divides(A, candidate)))
        return candidate.map((c) => c * d);
    }
    currentBits = newBits;
  }
  throw new RangeError('Exception (n_nextprime). No larger single-limb prime exists.');
}

function normalized(a: readonly bigint[]): bigint[] {
  let n = a.length;
  while (n && a[n - 1] === 0n) n--;
  return a.slice(0, n);
}
function integerGcd(a: bigint, b: bigint): bigint {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b) [a, b] = [b, a % b];
  return a;
}
function content(a: readonly bigint[]): bigint {
  let g = 0n;
  for (const c of a) {
    g = integerGcd(g, c);
    if (g === 1n) break;
  }
  return g;
}
function primitive(a: readonly bigint[]): bigint[] {
  const A = normalized(a);
  if (!A.length) return A;
  let c = content(A);
  if (A[A.length - 1]! < 0n) c = -c;
  return A.map((v) => v / c);
}
function ordered(a: readonly bigint[], b: readonly bigint[]): [bigint[], bigint[]] {
  const A = normalized(a),
    B = normalized(b);
  return A.length >= B.length ? [A, B] : [B, A];
}
/** fmpz_poly/pseudo_rem_cohen.c:18-52; the omitted powers matter on degree drops. */
function pseudoRemainder(a: readonly bigint[], b: readonly bigint[]): bigint[] {
  if (a.length < b.length) return a.slice();
  if (b.length === 1) return [];
  let r = a.slice(),
    e = a.length - b.length + 1;
  const leadB = b[b.length - 1]!;
  while (r.length >= b.length) {
    const leadA = r[r.length - 1]!,
      shift = r.length - b.length;
    for (let i = 0; i < r.length - 1; i++) r[i] = r[i]! * leadB;
    for (let i = 0; i < b.length - 1; i++) r[shift + i] = r[shift + i]! - b[i]! * leadA;
    r.pop();
    r = normalized(r);
    e--;
  }
  const scale = leadB ** BigInt(e);
  return r.map((c) => c * scale);
}
function bits(n: bigint): number {
  return n === 0n ? 0 : (n < 0n ? -n : n).toString(2).length;
}
function maxBits(a: readonly bigint[]): number {
  let n = 0;
  for (const c of a) n = Math.max(n, bits(c));
  return n;
}
/** fmpz_poly_bit_pack with native BigInt replacing the GMP limb buffer. */
function pack(a: readonly bigint[], width: number): bigint {
  if (!a.length) return 0n;
  const sign = a[a.length - 1]! < 0n ? -1n : 1n,
    base = 1n << BigInt(width);
  let borrow = 0n;
  const fields: string[] = [];
  for (const c of a) {
    let field = sign * c - borrow;
    borrow = field < 0n ? 1n : 0n;
    if (borrow) field += base;
    fields.push(field.toString(2).padStart(width, '0'));
  }
  return sign * BigInt('0b' + fields.reverse().join(''));
}
/** Signed fields carry a borrow into the next coefficient (bit_unpack.c). */
function unpack(n: bigint, width: number, length: number): bigint[] {
  const base = 1n << BigInt(width),
    half = base >> 1n;
  const source = n.toString(2).padStart(width * length, '0');
  const out: bigint[] = [];
  let borrow = 0n;
  for (let i = 0; i < length; i++) {
    const end = source.length - width * i;
    const field = BigInt('0b' + source.slice(end - width, end)),
      negative = field >= half;
    out.push(field - (negative ? base : 0n) + borrow);
    borrow = negative ? 1n : 0n;
  }
  return normalized(out);
}
/** Exact signed Kronecker substitution; small products use the classical kernel. */
function multiply(a: readonly bigint[], b: readonly bigint[]): bigint[] {
  if (!a.length || !b.length) return [];
  const length = a.length + b.length - 1;
  if (Math.min(a.length, b.length) <= 5) {
    const r = Array<bigint>(length).fill(0n);
    for (let i = 0; i < a.length; i++)
      for (let j = 0; j < b.length; j++) r[i + j] = r[i + j]! + a[i]! * b[j]!;
    return normalized(r);
  }
  const width = maxBits(a) + maxBits(b) + bits(BigInt(Math.min(a.length, b.length))) + 1;
  // Keep temporary packed integers below the engine's BigInt size limit.
  // mul_karatsuba.c splits even/odd coefficients (native buffers use revbin order).
  if (length * width > 1 << 19) {
    const a0 = a.filter((_, i) => i % 2 === 0),
      a1 = a.filter((_, i) => i % 2 === 1);
    const b0 = b.filter((_, i) => i % 2 === 0),
      b1 = b.filter((_, i) => i % 2 === 1);
    const add = (x: readonly bigint[], y: readonly bigint[]) =>
      normalized(
        Array.from({ length: Math.max(x.length, y.length) }, (_, i) => (x[i] ?? 0n) + (y[i] ?? 0n))
      );
    const low = multiply(a0, b0),
      high = multiply(a1, b1),
      middle = multiply(add(a0, a1), add(b0, b1));
    const result = Array<bigint>(length + 2).fill(0n);
    for (let i = 0; i < low.length; i++) {
      result[2 * i] = result[2 * i]! + low[i]!;
      middle[i] = (middle[i] ?? 0n) - low[i]!;
    }
    for (let i = 0; i < high.length; i++) {
      result[2 * i + 2] = result[2 * i + 2]! + high[i]!;
      middle[i] = (middle[i] ?? 0n) - high[i]!;
    }
    for (let i = 0; i < middle.length; i++)
      result[2 * i + 1] = result[2 * i + 1]! + (middle[i] ?? 0n);
    return normalized(result);
  }
  const product = pack(a, width) * pack(b, width);
  const sign = product < 0n ? -1n : 1n;
  return unpack(sign * product, width, length).map((c) => sign * c);
}
/** fmpz_poly/divides.c: reject using constant terms and evaluation at one first. */
function divides(a: readonly bigint[], b: readonly bigint[]): boolean {
  if (a.length < b.length) return false;
  const divisible = (x: bigint, y: bigint) => (y === 0n ? x === 0n : x % y === 0n);
  if (!divisible(a[0] ?? 0n, b[0] ?? 0n)) return false;
  if (
    a.length > 1 &&
    !divisible(
      a.reduce((s, c) => s + c, 0n),
      b.reduce((s, c) => s + c, 0n)
    )
  )
    return false;
  const result = _fmpz_poly_divrem(a, b, true);
  return result !== null && result[1].length === 0;
}

function mod(a: bigint, p: bigint): bigint {
  const r = a % p;
  return r < 0n ? r + p : r;
}
function inverse(a: bigint, p: bigint): bigint {
  let r = p,
    s = mod(a, p),
    x = 0n,
    y = 1n;
  while (s) {
    const q = r / s;
    [r, s] = [s, r - q * s];
    [x, y] = [y, x - q * y];
  }
  if (r !== 1n) throw new RangeError('coefficient is not invertible');
  return mod(x, p);
}
/** ulong_extras/is_prime.c: the 64-bit modular prime table and mod-30 wheel.
 * @see Deviation: Dense FLINT polynomial GCD kernels
 */
function* modularPrimes(): Generator<bigint> {
  const base = 1n << 63n;
  const offsets = [
    29, 99, 123, 131, 155, 255, 269, 359, 435, 449, 453, 485, 491, 543, 585, 599, 753, 849, 879,
    885, 903, 995, 1209, 1251, 1311, 1373, 1403, 1485, 1533, 1535, 1545, 1551, 1575, 1601, 1625,
    1655, 1701, 1709, 1845, 1859, 1913, 1995, 2045, 2219, 2229, 2321, 2363, 2385, 2483, 2499, 2523,
    2543, 2613, 2639, 2679, 2829, 2931, 3089, 3165, 3189, 3245, 3273, 3291, 3341,
  ];
  for (const offset of offsets) yield base + BigInt(offset);
  const wheel = [
    1, 6, 5, 4, 3, 2, 1, 4, 3, 2, 1, 2, 1, 4, 3, 2, 1, 2, 1, 4, 3, 2, 1, 6, 5, 4, 3, 2, 1, 2,
  ];
  let p = base + 3341n;
  while (p < 18446744073709551557n) {
    p += BigInt(wheel[Number(p % 30n)]!);
    if (wordIsPrime(p)) yield p;
  }
}
/** Deterministic Miller-Rabin for unsigned 64-bit inputs. */
function wordIsPrime(n: bigint): boolean {
  let d = n - 1n,
    s = 0;
  while (d % 2n === 0n) {
    d /= 2n;
    s++;
  }
  const power = (a: bigint, e: bigint) => {
    let r = 1n;
    for (; e; e >>= 1n, a = (a * a) % n) if (e & 1n) r = (r * a) % n;
    return r;
  };
  for (const base of [2n, 325n, 9375n, 28178n, 450775n, 9780504n, 1795265022n]) {
    if (base % n === 0n) continue;
    let x = power(base % n, d);
    if (x === 1n || x === n - 1n) continue;
    let passed = false;
    for (let r = 1; r < s; r++) {
      x = (x * x) % n;
      if (x === n - 1n) {
        passed = true;
        break;
      }
    }
    if (!passed) return false;
  }
  return true;
}

/** @internal Shared exact arithmetic and native word-prime stream for resultants. */
export const _fmpz_poly_resultant_kernels = {
  normalized,
  content,
  pseudoRemainder,
  maxBits,
  bits,
  modularPrimes,
  wordIsPrime,
  mod,
  inverse,
  multiply,
};
