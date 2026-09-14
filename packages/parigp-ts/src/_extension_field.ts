/** Native PARI Fq / Flxq / F2xq coefficient operations and tagged outer arithmetic.
 * Word linear operations share the existing reduced-coefficient implementation.
 * @see Deviation: PARI extension-polynomial division adapters
 */
import {
  extensionPolynomial as arithmetic,
  type ExtensionCoefficient as C,
  type ExtensionPolynomial as P,
} from './_extension_polynomial.js';
import { FpX_red, FpX_add, FpX_sub, FpX_neg, FpX_Fp_mul } from './ffinit.js';
import { FpX_mul, FpX_sqr, FpX_rem, FpX_divrem, _FpX_extgcd } from './FpX.js';
import { Flx_mul, Flx_sqr, Flx_rem, Flx_divrem, _Flx_extgcd } from './Flx.js';
import { F2x_mul, F2x_sqr, F2x_rem, F2x_divrem, F2x_degree } from './F2x.js';
import { inverseCoefficient, divisionBarrett } from './_polynomial_division.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { PariError } from './errors.js';
export const extensionZero = (c: C) =>
  typeof c === 'bigint' ? c === 0n : !c.some((x) => x !== 0n);
export const extensionOne = (c: C) =>
  typeof c === 'bigint' ? c === 1n : c[0] === 1n && c.slice(1).every((x) => x === 0n);
export const trimExtension = (a: P) => {
  let n = a.length;
  while (n && zero(a[n - 1]!)) n--;
  return a.slice(0, n);
};
const zero = extensionZero,
  trim = trimExtension;
import {displayCoefficientError as display} from './_extension_display.js';
export function extensionField(mode: 0 | 1 | 2, T: bigint[] | bigint, p: bigint, innerInverse?: bigint[]) {
  const mod = (a: bigint) => ((a % p) + p) % p,
    z: C = mode === 1 ? [] : 0n,
    u: C = mode === 1 ? [1n] : 1n;
  const reducePolynomial = (a: bigint[]): bigint[] =>
    innerInverse === undefined
      ? mode === 1 ? Flx_rem(a, T as bigint[], p) : FpX_rem(a, T as bigint[], p)
      : divisionBarrett(
          a, T as bigint[], innerInverse, p,
          (a,b) => mode === 1 ? Flx_mul(a,b,p) : FpX_mul(a,b,p),
          (a,b) => mode === 1 ? Flx_divrem(a,b,p) : FpX_divrem(a,b,p)
        )[1];
  const red = (a: C): C =>
    mode === 2
      ? F2x_rem(a as bigint, T as bigint)
      : typeof a === 'bigint'
        ? mod(a)
        : reducePolynomial(mode === 1 ? a : FpX_red(a,p));
  const neg = (a: C): C => (mode === 2 ? a : typeof a === 'bigint' ? mod(-a) : FpX_neg(a, p));
  const add = (a: C, b: C): C => {
    if (mode === 2) return (a as bigint) ^ (b as bigint);
    if (typeof a === 'bigint' && typeof b === 'bigint') return mod(a + b);
    if (typeof a !== 'bigint' && typeof b !== 'bigint') return FpX_add(a, b, p);
    if (typeof a === 'bigint') [a, b] = [b, a];
    const x = a as bigint[],
      y = b as bigint;
    return trimPolynomial([mod((x[0] ?? 0n) + y), ...x.slice(1)]);
  };
  const sub = (a: C, b: C): C => {
    if (mode === 2) return (a as bigint) ^ (b as bigint);
    if (typeof a === 'bigint' && typeof b === 'bigint') return mod(a - b);
    if (typeof a !== 'bigint' && typeof b !== 'bigint') return FpX_sub(a, b, p);
    if (typeof b === 'bigint') {
      const x = a as bigint[];
      return trimPolynomial([mod((x[0] ?? 0n) - b), ...x.slice(1)]);
    }
    return trimPolynomial([mod((a as bigint) - (b[0] ?? 0n)), ...b.slice(1).map((c) => mod(-c))]);
  };
  const mul = (a: C, b: C, reduce = true): C => {
    if (mode === 2) {
      const v = F2x_mul(a as bigint, b as bigint);
      return reduce ? F2x_rem(v, T as bigint) : v;
    }
    if (typeof a === 'bigint' && typeof b === 'bigint') return mod(a * b);
    if (typeof a === 'bigint') return FpX_Fp_mul(b as bigint[], a, p);
    if (typeof b === 'bigint') return FpX_Fp_mul(a, b, p);
    const v = mode === 1 ? Flx_mul(a, b, p) : FpX_mul(a, b, p);
    return reduce ? reducePolynomial(v) : v;
  };
  const sqr = (a: C): C => {
    if (mode === 2) return F2x_rem(F2x_sqr(a as bigint), T as bigint);
    if (typeof a === 'bigint') return mod(a * a);
    return reducePolynomial(mode === 1 ? Flx_sqr(a,p) : FpX_sqr(a,p));
  };
  // Native safe inverse returns null for a nonunit; Fl_inv may still throw over composites.
  const scalarSafe = (a:bigint):bigint|null => {
    try { return inverseCoefficient(a,p,false); }
    catch(error) { if(error instanceof PariError)return null;throw error; }
  };
  const invsafe = (a:C):C|null => {
    if(mode===2) {
      let x=T as bigint,y=a as bigint,v=0n,w=1n;
      while(y) {
        const [q,r]=F2x_divrem(x,y);
        [v,w]=[w,v^F2x_mul(q,w)];[x,y]=[y,r];
      }
      return F2x_degree(x)===0?v:null;
    }
    if(typeof a==='bigint')return scalarSafe(a);
    const [g,,v]=mode===1
      ? _Flx_extgcd(T as bigint[],a,p,false)
      : _FpX_extgcd(FpX_red(T as bigint[],p),FpX_red(a,p),p,false);
    if(g.length!==1)return null;
    const inverse=mode===1?inverseCoefficient(g[0]!,p,true):scalarSafe(g[0]!);
    return inverse===null?null:FpX_Fp_mul(v,inverse,p);
  };
  const inv = (a: C): C => {
    const error = () => {
      const f =
        mode === 2
          ? Array.from(
              { length: F2x_degree(a as bigint) + 1 },
              (_, i) => ((a as bigint) >> BigInt(i)) & 1n
            )
          : (a as bigint[]);
      throw new PariError(
        `impossible inverse in ${mode === 0 ? 'FpXQ_inv' : mode === 1 ? 'Flxq_inv' : 'F2xq_inv'}: ${display(f)}.`
      );
    };
    if (mode === 2) {
      let x = T as bigint,
        y = a as bigint,
        v = 0n,
        w = 1n;
      while (y) {
        const [q, r] = F2x_divrem(x, y);
        [v, w] = [w, v ^ F2x_mul(q, w)];
        [x, y] = [y, r];
      }
      if (F2x_degree(x) !== 0) error();
      return v;
    }
    if (typeof a === 'bigint') return inverseCoefficient(a, p, false);
    const [g, , v] =
      mode === 1
        ? _Flx_extgcd(T as bigint[], a, p, false)
        : _FpX_extgcd(FpX_red(T as bigint[], p), FpX_red(a, p), p, false);
    if (g.length !== 1) error();
    return FpX_Fp_mul(v, inverseCoefficient(g[0]!, p, mode === 1), p);
  };
  const pmul = (a: P, b: P) => arithmetic(mode, 0, p, T, a, b, innerInverse);
  const pred = (a: P) => trim(a.map(red));
  const pscale = (a: P, b: C) => trim(a.map((c) => mul(c, b)));
  const psub = (a: P, b: P) =>
    trim(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) =>
        i < a.length ? (i < b.length ? sub(a[i]!, b[i]!) : a[i]!) : neg(b[i]!)
      )
    );
  return { mode, T, p, z, u, red, neg, add, sub, mul, sqr, inv, invsafe, pmul, pred, pscale, psub };
}
export type ExtensionField = ReturnType<typeof extensionField>;
