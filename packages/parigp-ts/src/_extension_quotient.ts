/** Native extension quotient products, inverses and powering schedules.
 * @see Deviation: PARI extension quotient adapters
 */
import {
  extensionPolynomial as arithmetic,
  validateExtensionInputs,
  type ExtensionPolynomial as P,
} from './_extension_polynomial.js';
import { extensionField, trimExtension as trim } from './_extension_field.js';
import {
  extensionDivision,
  extensionGetRed,
  type ExtensionModulus as Modulus,
} from './_extension_division.js';
import { extensionGcd } from './_extension_gcd.js';
import { polynomialQuotient } from './_polynomial_quotient.js';
import { FpX_red } from './ffinit.js';
import { gen_pow_i, gen_powu_i } from './bb_group.js';
import { displayExtensionError } from './_extension_display.js';
import { PariError } from './errors.js';
export function extensionQuotient(
  mode: 0 | 1 | 2,
  op: 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7,
  p: bigint,
  n: bigint,
  T: bigint[] | bigint,
  a: P,
  b: P,
  modulus: Modulus,
  providedInner?: bigint[]
): P | P[] | null {
  let S = Array.isArray(modulus) ? modulus : modulus.polynomial;
  validateExtensionInputs(mode, p, T, a, b);
  validateExtensionInputs(mode, p, T, S);
  if (Array.isArray(T)) {
    T = T.slice();
    while (T.length && T.at(-1) === 0n) T.pop();
  }
  const canonical = (v: P) =>
    trim(
      v.map((c) => {
        if (typeof c === 'bigint') return c;
        let n = c.length;
        while (n && c[n - 1] === 0n) n--;
        return c.slice(0, n);
      })
    );
  a = canonical(a);
  b = canonical(b);
  S = canonical(S);
  let innerInverse = providedInner;
  const prepareInner = () => {
    if (mode !== 2 && innerInverse === undefined)
      innerInverse = polynomialQuotient(T as bigint[], p, mode === 1).inverse;
  };
  let f = extensionField(mode, T, p, innerInverse);
  const prepare = () => {
    prepareInner();
    f = extensionField(mode, T, p, innerInverse);
  };
  let source: Modulus = Array.isArray(modulus) ? S : modulus;
  const rem = (x: P) => extensionDivision(mode, 1, p, T, x, source, innerInverse) as P;
  const mul = (x: P, y: P) => rem(arithmetic(mode, 0, p, T, x, y, innerInverse));
  const sqr = (x: P) => rem(arithmetic(mode, 1, p, T, x, [], innerInverse));
  const inverse = (x: P, safe = false): P | null => {
    const [g, v] = extensionGcd(mode, 3, p, T, S, x, innerInverse) as [P, P];
    const c = g.length === 1 ? f.invsafe(g[0]!) : null;
    if (c === null) {
      if (safe) return null;
      const message = displayExtensionError(x, mode);
      throw new PariError(
        'impossible inverse in ' +
          (mode === 0 ? 'FpXQXQ_inv' : mode === 1 ? 'FlxqXQ_inv' : 'F2xqXQ_inv') +
          ': ' +
          message +
          '.'
      );
    }
    return f.pscale(v, c);
  };
  if (op === 0) return mul(a, b);
  if (op === 1) return sqr(a);
  if (op === 2) return inverse(a, true);
  if (op === 3) return inverse(a);
  if (op === 4) return mul(a, inverse(b) as P);
  if (op === 6) {
    if (n < 0n || n > BigInt(0xffffffff - 1))
      throw new RangeError('power count must be a nonnegative array length');
    const count = Number(n),
      useSqr = 2 * (a.length - 1) >= S.length - 1;
    prepare();
    source = extensionGetRed(mode, source, T, p, innerInverse);
    const powers: P[] = [[f.u]];
    if (count === 0) return powers;
    powers.push(a.slice());
    if (count === 1) return powers;
    powers.push(sqr(a));
    for (let i = 3; i <= count; i++)
      powers.push(useSqr && i % 2 === 0 ? sqr(powers[i / 2]!) : mul(powers[i - 1]!, a));
    return powers;
  }
  if (op === 7) {
    if (n < 0n || n >= 1n << 64n) throw new RangeError('exponent must be an unsigned word integer');
    if (n === 0n) return [f.u];
    if (n === 1n) return a.slice();
    if (n === 2n) return sqr(a);
    prepare();
    source = extensionGetRed(mode, source, T, p, innerInverse);
    return gen_powu_i(a, n, sqr, mul);
  }
  if (n === 0n) return [f.u];
  if (mode === 0) {
    if (n === 1n) return a.slice();
    if (n === -1n) return inverse(a);
    if (p < 1n << 64n) {
      const convert = (v: P) => trim(v.map((c) => FpX_red(typeof c === 'bigint' ? [c] : c, p)));
      const wordS = Array.isArray(source)
        ? convert(S)
        : (() => {
            if (S.some((c) => typeof c !== 'bigint'))
              throw new RangeError(
                'generic word powering cannot convert cached polynomial coefficients'
              );
            return trim([FpX_red(S as bigint[], p)]);
          })();
      const r = extensionQuotient(
        1,
        op,
        p,
        n,
        FpX_red(T as bigint[], p),
        convert(a),
        convert(b),
        wordS,
        innerInverse === undefined ? undefined : FpX_red(innerInverse, p)
      ) as P;
      return r.map((c) => ((c as bigint[]).length < 2 ? ((c as bigint[])[0] ?? 0n) : c));
    }
    prepare();
    source = extensionGetRed(mode, source, T, p, innerInverse);
    if (n < 0n) a = inverse(a) as P;
  } else {
    if (n < 0n) a = inverse(a) as P;
    if (n === 1n || n === -1n) return a.slice();
    if (a.length >= S.length) a = rem(a);
    prepare();
    source = extensionGetRed(mode, source, T, p, innerInverse);
  }
  return gen_pow_i(a, n, sqr, mul);
}
