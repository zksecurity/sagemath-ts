/** Strict bounded integer factorization from PARI ifactor1.c:3369-3507,
 * 4279-4465 and 4493. Returns [known factors, unresolved [base, exponent]].
 * Prime-table and primorial boundaries match the initialized native oracle.
 * @see Deviation: PARI strict partial-factorization state
 */
import { PariError } from './errors.js';
import { forprime, Z_iroot, Z_isanypower, Z_factor, isPrime, isqrt } from './ifactor.js';
import { WordPrimeIterator } from './language/forprime.js';
import { gcd } from './ff.js';
type Factors = [bigint, bigint][];
const cache = new Map<number, bigint>();
function primorial(lim: bigint): [bigint, bigint] {
  const cap = Math.min(500000, 2 ** lim.toString(2).length);
  let g = cache.get(cap);
  const primes = Array.from(forprime(3, cap));
  if (g === undefined) {
    let work = primes.map(BigInt);
    while (work.length > 1) {
      const next = [];
      for (let i = 0; i < work.length; i += 2)
        next.push(i + 1 < work.length ? work[i]! * work[i + 1]! : work[i]!);
      work = next;
    }
    g = work[0] ?? 1n;
    cache.set(cap, g);
  }
  return [g, BigInt(primes.at(-1)!)];
}
function smallPower(n: bigint, all: bigint, divisor?: bigint): [bigint, bigint] {
  let e = 1n;
  const exps =
    divisor === undefined
      ? all > 7129n
        ? [2, 3]
        : all > 563n
          ? [2, 3, 5]
          : [2, 3, 5, 7]
      : [2, 3, 5, 7];
  for (const k of exps)
    for (;;) {
      if (divisor !== undefined && divisor % BigInt(k) !== 0n) break;
      const [r, yes] = Z_iroot(n, k);
      if (!yes) break;
      n = r;
      e *= BigInt(k);
      if (divisor !== undefined) divisor /= BigInt(k);
    }
  return [n, e];
}
export function absZ_factor_limit_strict(
  value: bigint,
  bound: bigint
): [Factors, [bigint, bigint] | null] {
  bound = bound < 0n ? -bound : bound;
  if (bound >= 1n << 64n) throw new PariError('overflow in t_INT-->ulong assignment.');
  let n = value < 0n ? -value : value;
  const all = bound === 0n ? 500001n : bound,
    word = n < 1n << 64n,
    f: Factors = [];
  if (n === 0n) return [[[0n, 1n]], null];
  if (n === 1n) return [f, null];
  let lim = word ? 0n : 1n,
    knownCompositeAt = -1;
  const take = (p: bigint): boolean => {
    let e = 0n;
    while (n % p === 0n) {
      n /= p;
      e++;
    }
    if (e) f.push([p, e]);
    return e !== 0n;
  };
  const stopAfter = (p: bigint) => n <= p * p;
  const primeFinish = (): [Factors, null] => {
    if (n !== 1n) f.push([n, 1n]);
    return [f, null];
  };
  const scan = (lo: bigint, hi: bigint): boolean => {
    const it = new WordPrimeIterator(lo);
    for (let p = it.next(); p && p <= hi; p = it.next()) {
      if (word && p === 673n) {
        if (isPrime(n)) return true;
        knownCompositeAt = f.length;
      }
      take(p);
      let stop = stopAfter(p);
      if (
        !word &&
        !stop &&
        p === 16381n &&
        Math.ceil(n.toString(2).length / 64) * 64 * 0.0578 * Math.LN2 < 10
      ) {
        knownCompositeAt = f.length;
        stop = isPrime(n);
      }
      if (stop) return true;
    }
    return false;
  };
  if (all > 2n) {
    take(2n);
    if (n === 1n) return [f, null];
    if (word && all > 3n) {
      for (const p of [3n, 5n, 7n])
        if (all > p) {
          take(p);
          if (n === 1n) return [f, null];
        }
    }
    if (word && n <= 499979n && isPrime(n)) return primeFinish();
    lim = all - 1n;
    let maxp: bigint;
    if (lim >= 128n && (!word || n >= 691n * 691n)) {
      if (word) {
        const root = isqrt(n);
        if (lim > root) lim = root;
      }
      const [product, cacheLim] = primorial(lim),
        g = gcd(n, product);
      let factors = Z_factor(g).map(([p]) => p);
      if (word) factors = factors.filter((p) => p >= 11n);
      maxp = 500000n;
      if (factors.length) {
        let stop = false;
        for (const p of factors) {
          if (p >= all) {
            if (!word) break;
            let rest = n,
              e = 0n;
            while (rest % p === 0n) {
              rest /= p;
              e++;
            }
            if (rest === 1n) {
              f.push([p, e]);
              return [f, null];
            }
            const [base, k] = smallPower(rest, all, e);
            return [f, [k === 1n ? n : base * p ** (e / k), k]];
          }
          take(p);
          stop = stopAfter(p);
        }
        if (n === 1n) return [f, null];
        if (word && (stop || (n <= 499979n && isPrime(n)))) return primeFinish();
      } else if (word && lim === isqrt(n) && lim <= 500000n) return primeFinish();
    } else {
      maxp = word ? lim : 499979n;
      if (scan(word ? 11n : 3n, lim < maxp ? lim : maxp)) return primeFinish();
    }
    if (lim > maxp && scan(maxp + 1n, lim)) return primeFinish();
  }
  if (word) {
    const [base, e] = smallPower(n, all);
    if (knownCompositeAt === f.length || !isPrime(base)) return [f, [base, e]];
    f.push([base, e]);
    return [f, null];
  }
  const [power, base] = Z_isanypower(n),
    e = BigInt(power || 1),
    known = power > 1 ? -1 : knownCompositeAt;
  if (
    base <= lim * lim ||
    (lim >= 16384n &&
      f.length > known &&
      Math.ceil(base.toString(2).length / 64) * 64 < 2048 &&
      isPrime(base))
  ) {
    f.push([base, e]);
    return [f, null];
  }
  return [f, [base, e]];
}
