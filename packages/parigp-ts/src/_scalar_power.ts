/** Scalar powering schedules and reductions from PARI basemath/arith1.c.
 * @see Deviation: PARI scalar exponentiation boundaries
 */
import { gen_pow_fold, gen_pow_i, gen_powu_i } from './bb_group.js';
import { inverseCoefficient, residue } from './_polynomial_division.js';
import { PariError } from './errors.js';

const WORD = 1n << 64n;
const MASK = WORD - 1n;
const abs = (x: bigint): bigint => (x < 0n ? -x : x);
const bits = (x: bigint): number => (x === 0n ? 0 : abs(x).toString(2).length);
const shift = (x: bigint, n: number): bigint => (n >= 0 ? x << BigInt(n) : x / (1n << BigInt(-n)));
const remainder = (x: bigint, p: bigint): bigint => {
  if (p === 0n) throw new PariError('impossible inverse in dvmdii: 0.');
  return x % p;
};
const double = (x: bigint, p: bigint): bigint => (2n * x >= p ? 2n * x - p : 2n * x);

/** Fl_powu: reduced word input, with right-to-left or fused base-two powering. */
function wordPower(a: bigint, n: bigint, p: bigint): bigint {
  if (n === 0n) return 1n;
  if (n === 1n) return a;
  if (n === 2n) return (a * a) % p;
  if (a <= 1n) return a;
  if (a === 2n)
    return gen_pow_fold(
      a,
      n,
      (x) => (x * x) % p,
      (x) => (2n * ((x * x) % p)) % p
    );
  let y = 1n;
  for (;;) {
    if (n & 1n) y = (y * a) % p;
    n >>= 1n;
    if (n === 0n) return y;
    a = (a * a) % p;
  }
}

type Reduction = {
  base: bigint;
  square: (x: bigint) => bigint;
  multiply: (x: bigint, y: bigint) => bigint;
  squareDouble: (x: bigint) => bigint;
  finish: (x: bigint) => bigint;
};

/** Native 64-bit GMP tuning cutoffs: Barrett 127 GEN words, Montgomery 17. */
function selectReduction(a: bigint, n: bigint, p: bigint): Reduction {
  const size = Math.ceil(bits(p) / 64),
    length = size + 2;
  if (length >= 127 && (n === 0n || Number(n) * (bits(a) - 1) > bits(p) + 1)) {
    const s = 1 + Math.floor((bits(p) - 1) / 2),
      scale = 1n << BigInt(3 * s);
    const P = scale / p,
      Q = scale % p,
      t = bits(P);
    // Modified Barrett, with one fold and at most three conditional subtractions.
    const reduce = (x: bigint): bigint => {
      const A = (x % scale) + Q * shift(x, -3 * s);
      const q = shift(shift(A, t - 3 * s) * P, -t);
      let r = A - q * p;
      for (let i = 0; i < 3; i++) {
        const next = r - p;
        if (next < 0n) return r;
        r = next;
      }
      return r;
    };
    return {
      base: a,
      square: (x) => reduce(x * x),
      multiply: (x, y) => reduce(x * y),
      squareDouble: (x) => double(reduce(x * x), p),
      finish: (x) => x,
    };
  }
  if ((p & 1n) !== 0n && length < 17) {
    const magnitude = abs(p),
      R = 1n << BigInt(64 * size);
    // Newton inverse modulo 2^64, matching invmod2BIL and magnitude-only mod2BIL.
    let inverse = 1n;
    for (let i = 0; i < 6; i++) inverse = (inverse * (2n - magnitude * inverse)) & MASK;
    inverse = -inverse & MASK;
    const reduce = (x: bigint): bigint => {
      // The native GMP kernel copies magnitude limbs, ignoring GEN sign.
      x = abs(x);
      for (let i = 0; i < size; i++) {
        const m = ((x & MASK) * inverse) & MASK;
        x = (x + m * magnitude) >> 64n;
      }
      return x >= R ? x - magnitude : x;
    };
    return {
      base: remainder(a * R, p),
      square: (x) => reduce(x * x),
      multiply: (x, y) => reduce(x * y),
      squareDouble: (x) => {
        let z = 2n * reduce(x * x);
        // For negative N this native subtraction loop cannot decrease z.
        if (p < 0n && z >= R)
          throw new RangeError('Fp_pow: negative modulus makes Montgomery doubling nonterminating');
        while (z >= R) z -= p;
        return z;
      },
      finish: (x) => {
        const y = reduce(x);
        return y >= p ? y - p : y;
      },
    };
  }
  return {
    base: a,
    square: (x) => remainder(x * x, p),
    multiply: (x, y) => remainder(x * y, p),
    squareDouble: (x) => double(remainder(x * x, p), p),
    finish: (x) => x,
  };
}

/** Fp_powu's internal unsigned exponent path; its small powers precede reduction. */
export function unsignedPower(a: bigint, n: bigint, p: bigint): bigint {
  const pp = abs(p);
  if (pp > 0n && pp < WORD) return wordPower(residue(a, pp), n, pp);
  if (n === 2n) return remainder(a * a, p);
  if (n === 1n) return a;
  if (n === 0n) return 1n;
  a = residue(a, p);
  if (abs(a) === 1n) return 1n;
  const baseTwo = abs(a) === 2n,
    D = selectReduction(a, n, p);
  return D.finish(
    baseTwo
      ? gen_pow_fold(D.base, n, D.square, D.squareDouble)
      : gen_powu_i(D.base, n, D.square, D.multiply)
  );
}

/** Fp_pow, including zero-divisibility and exponent-width inverse dispatch. */
export function scalarPower(a: bigint, n: bigint, p: bigint): bigint {
  if (n === 0n) return a === 0n || (p !== 0n && a % p === 0n) ? 0n : 1n;
  const pp = abs(p),
    k = abs(n);
  if (pp > 0n && pp < WORD && k < WORD) {
    let y = residue(a, pp);
    if (n < 0n) y = inverseCoefficient(y, pp, true);
    return y <= 1n ? y : wordPower(y, k, pp);
  }
  let y: bigint;
  if (n < 0n) y = inverseCoefficient(a, p, false);
  else {
    y = residue(a, p);
    if (y === 0n) return 0n;
  }
  if (k < WORD) return unsignedPower(y, k, p);
  const sy = abs(y) > abs(p / 2n),
    negative = sy && (k & 1n) !== 0n;
  if (sy) y = p - y;
  if (abs(y) === 1n) return negative ? p - 1n : 1n;
  const baseTwo = abs(y) === 2n,
    D = selectReduction(y, 0n, p);
  y = D.finish(
    baseTwo
      ? gen_pow_fold(D.base, k, D.square, D.squareDouble)
      : gen_pow_i(D.base, k, D.square, D.multiply)
  );
  return negative ? p - y : y;
}
