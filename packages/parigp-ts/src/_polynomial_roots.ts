/** Native deterministic root splitting and root-count predicates, FpX_factor.c.
 * @see Deviation: PARI polynomial root adapters
 */
import { Flx_divrem, Flx_gcd, Flx_normalize, Flx_rem, Flx_sqr } from './Flx.js';
import { FpX_divrem, FpX_gcd, FpX_normalize, FpX_sqr } from './FpX.js';
import { nonsquare1_Fl } from './_modular_sqrt.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { polynomialQuotient } from './_polynomial_quotient.js';
import { gen_pow_i } from './bb_group.js';
import { PariError } from './errors.js';
import { Fl_sqrt, Fp_sqrt, kronecker } from './ff.js';
import { FpX_add, FpX_red, FpX_sub } from './ffinit.js';
import { isPrime } from './ifactor.js';
type Poly = bigint[];
const mod = (x: bigint, p: bigint) => ((x % p) + p) % p;
const nonsquare = (p: bigint) => ((p & 3n) === 3n ? p - 1n : nonsquare1_Fl(p));
export function polynomialRoots(input: Poly, p: bigint, direct = false): bigint[] {
  if (direct) input = wordRootPolynomial(input, p);
  if (p < 2n) throw new RangeError('polynomial roots require a modulus greater than one');
  const word = direct || p < 1n << 64n;
  let f = direct ? trimPolynomial(input) : FpX_red(input, p);
  if (!f.length) throw new PariError(`zero polynomial in ${direct ? 'Flx' : 'FpX'}_roots.`);
  if (f.length === 1) return [];
  if (p === 2n) {
    const out = [];
    if (!f[0]) out.push(0n);
    if (!f.reduce((a, b) => a ^ b, 0n)) out.push(1n);
    return out;
  }
  if (!direct && word && !(p & 1n)) throw new PariError(`not a prime number in FpX_roots: ${p}.`);
  const norm = (a: Poly) => (word ? Flx_normalize(a, p) : FpX_normalize(a, p));
  const gcd = (a: Poly, b: Poly) => (word ? Flx_gcd(a, b, p) : FpX_gcd(a, b, p));
  const divide = (a: Poly, b: Poly) => (word ? Flx_divrem(a, b, p)[0] : FpX_divrem(a, b, p)[0]);
  const power = (x: Poly, e: bigint, T: Poly): Poly => {
    if (e === 1n && word) return x.slice();
    if (e === 2n && word) return Flx_rem(Flx_sqr(x, p), T, p);
    const ctx = polynomialQuotient(T, p, word);
    return gen_pow_i(
      x,
      e,
      (a) => ctx.reduce(word ? Flx_sqr(a, p) : FpX_sqr(a, p)),
      (a, b) => ctx.reduce(ctx.multiply(a, b))
    );
  };
  const quad = (a: Poly, unknown: boolean): bigint[] => {
    const D = word ? mod(a[1]! * a[1]! - 4n * a[0]!, p) : (a[1]! * a[1]! - 4n * a[0]!) % p;
    if (unknown && kronecker(D, p) === -1) return [];
    const s = word ? Fl_sqrt(D, p) : Fp_sqrt(D, p);
    if (s === null) return [];
    let r = mod(s - a[1]!, p);
    r = (r + (r & 1n ? p : 0n)) >> 1n;
    return [r, mod(-a[1]! - r, p)];
  };
  const done: bigint[] = [];
  const todo: Poly[] = [];
  const strip = (a: Poly) => {
    let v = 0;
    while (a[v] === 0n) v++;
    return a.slice(v);
  };
  if (!word) f = norm(f);
  if (f[0] === 0n) {
    f = strip(f);
    done.push(0n);
  }
  if (word) f = norm(f);
  if (f.length === 1) return done;
  if (f.length === 2) return [...done, p - f[0]!];
  if (f.length === 3) {
    const r = quad(f, true);
    return [...done, ...new Set(r)].sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  }
  const primeError = (name: string): never => {
    throw new PariError(`not a prime number in ${name}: ${p}.`);
  };
  const addRoots = (g: Poly) => {
    const a = power([0n, 1n], p - 1n, g);
    if (!a.length) primeError('rootmod');
    const b = gcd(g, FpX_sub(a, [1n], p));
    if (b.length > 1) todo.push(norm(b));
  };
  const fold = (a: Poly, n: bigint, plus: boolean): Poly => {
    // Flx_mod_Xnm1 / Flx_mod_Xn1 preserve identity when no folding is needed.
    if (BigInt(a.length) <= n) return a;
    const N = Number(n);
    const out = a.slice(0, N);
    for (let i = N; i < a.length; i++)
      out[i % N] = mod(out[i % N]! + (plus && Math.floor(i / N) % 2 ? -a[i]! : a[i]!), p);
    return trimPolynomial(out);
  };
  if (!word) addRoots(f);
  else {
    let g = fold(f, p - 1n, false);
    if (!g.length)
      return Array.from(
        { length: Number(p) - (done.length ? 0 : 1) },
        (_, i) => BigInt(i) + (done.length ? 0n : 1n)
      );
    if (g !== f) g = strip(g);
    if (g.length > 1) {
      // Native small-p cutoff: split squares/nonsquares before the shift loop.
      if (p >> 4n <= BigInt(g.length - 1)) {
        const q = p >> 1n;
        const xt =
          BigInt(g.length - 1) < q ? Flx_rem([...Array(Number(q)).fill(0n), 1n], g, p) : null;
        for (const plus of [false, true]) {
          let a = fold(g, q, plus);
          if (!a.length) {
            const z = plus ? nonsquare(p) : 1n;
            done.push(z);
            for (let i = 2n; i <= q; i++) done.push(mod(z * mod(i * i, p), p));
          } else {
            if (a !== g) a = strip(a);
            if (a.length > 1) {
              const t = xt
                ? FpX_add(xt, [plus ? 1n : p - 1n], p)
                : [plus ? 1n : p - 1n, ...Array(Number(q) - 1).fill(0n), 1n];
              a = gcd(a, t);
              if (a.length > 1) todo.push(norm(a));
            }
          }
        }
      } else addRoots(g);
    }
  }
  for (let c = 1n; todo.length; c++) {
    if (c === 100n && !isPrime(p)) primeError('polrootsmod');
    let l = todo.length;
    for (let j = 0; j < l; j++) {
      const a = todo[j]!;
      if (a.length <= 3) {
        if (a.length === 2) done.push(p - a[0]!);
        else {
          const r = quad(a, false);
          if (!r.length) primeError('polrootsmod');
          done.push(...r);
        }
        todo[j] = todo.at(-1)!;
        todo.pop();
        j--;
        l--;
        continue;
      }
      let b = power([c, 1n], p >> 1n, a);
      if (b.length <= 1) continue;
      b = gcd(a, FpX_sub(b, [1n], p));
      if (b.length <= 1) continue;
      b = norm(b);
      todo[j] = b;
      todo.push(divide(a, b));
    }
  }
  return done.sort((a, b) => {
    // vecsmall_sort compares signed longs; the quadratic shortcut above does not.
    if (word) {
      a = BigInt.asIntN(64, a);
      b = BigInt.asIntN(64, b);
    }
    return a < b ? -1 : a > b ? 1 : 0;
  });
}
export function polynomialTotallySplit(f: Poly, p: bigint, word = false): boolean {
  f = word ? wordRootPolynomial(f, p) : trimPolynomial(f);
  const n = f.length - 1;
  if (n >= 0 && n <= 1) return true;
  // The native degree is cast to ulong, including degree -1 for zero.
  if (BigInt.asUintN(64, BigInt(n)) > (p < 0n ? -p : p)) return false;
  if (p < 2n) throw new RangeError('total splitting requires a modulus greater than one');
  const a = FpX_red(f, p);
  const ctx = polynomialQuotient(a, p, word || p < 1n << 64n);
  const r = ctx.power([0n, 1n]);
  return r.length === 2 && r[0] === 0n && r[1] === 1n;
}

function wordRootPolynomial(f: Poly, p: bigint): Poly {
  if (p < 1n || p >= 1n << 64n) throw new RangeError('modulus must be a positive word integer');
  if (f.some((c) => c < 0n || c >= p))
    throw new RangeError('word polynomial coefficients must be reduced');
  return trimPolynomial(f);
}
export function polynomialWordRootCount(f: Poly, p: bigint): number {
  f = wordRootPolynomial(f, p);
  const n = f.length - 1;
  if (n <= 1) return n;
  if (n === 2) {
    if (p === 2n) return Number(f[0] === 0n) + Number(f[0] !== f[1]);
    return 1 + kronecker(mod(f[1]! * f[1]! - 4n * f[2]! * f[0]!, p), p);
  }
  // Flx_Frobenius -> Flxq_powu_pre has a direct square at p=2: no inverse cache.
  const frobenius =
    p === 2n ? Flx_rem(Flx_sqr([0n, 1n], p), f, p) : polynomialQuotient(f, p, true).power([0n, 1n]);
  return Flx_gcd(FpX_sub(frobenius, [0n, 1n], p), f, p).length - 1;
}
