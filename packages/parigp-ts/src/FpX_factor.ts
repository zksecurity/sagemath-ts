/** Dense polynomial factorization kernels from PARI basemath/FpX_factor.c. */
import { PariError } from './errors.js';
import {
  type FpX, FpX_divrem, FpX_gcd, FpX_normalize, FpX_red, FpX_sub,
} from './ffinit.js';
import { FpX_deriv } from './FpX.js';

/** PARI's Yun algorithm, used when p exceeds a native unsigned word. */
function FpX_factor_Yun(T: FpX, p: bigint): FpX[] {
  let d = FpX_deriv(T, p);
  let a = FpX_gcd(T, d, p);
  if (a.length === 1) return [T];
  let b = FpX_divrem(T, a, p)[0];
  const V: FpX[] = [];
  do {
    const c = FpX_divrem(d, a, p)[0];
    d = FpX_sub(c, FpX_deriv(b, p), p);
    a = FpX_normalize(FpX_gcd(b, d, p), p);
    V.push(a);
    b = FpX_divrem(b, a, p)[0];
  } while (b.length > 1);
  return V;
}

/** PARI Flx_factor_squarefree_pre, with exact BigInt modular arithmetic. */
function Flx_factor_squarefree(f: FpX, p: bigint): FpX[] {
  const n = f.length - 1;
  const u: FpX[] = Array.from({ length: n + 1 }, () => [1n]);
  for (let q = 1; ; q *= Number(p)) {
    let r = Flx_gcd(f, FpX_deriv(f, p), p);
    if (r.length === 1) {
      u[q - 1] = f;
      break;
    }
    let t = Flx_divrem(f, r, p)[0];
    if (t.length > 1) {
      for (let j = 1; ; j++) {
        const v = Flx_gcd(r, t, p);
        const tv = Flx_divrem(t, v, p)[0];
        if (tv.length > 1) u[j * q - 1] = Flx_normalize(tv, p);
        if (v.length <= 1) break;
        r = Flx_divrem(r, v, p)[0];
        t = v;
      }
      if (r.length === 1) break;
    }
    // Here the nonconstant remainder is a p-th power, so p <= deg(r).
    const deflated: FpX = [];
    for (let i = 0; i < r.length; i += Number(p)) deflated.push(r[i]!);
    f = Flx_normalize(deflated, p);
  }
  let i = n;
  while (i > 0 && u[i - 1]!.length === 1) i--;
  return u.slice(0, i);
}

/**
 * PARI FpX_factor_squarefree: component i contains factors of multiplicity i+1.
 * Missing multiplicities contain [1n]; these are not irreducible factors.
 * Coefficients are little endian, p is prime, and large-prime inputs are reduced.
 * The native word branch returns [] for constants; the large-prime branch [T].
 * @see Deviation: PARI squarefree zero-input safety
 */
export function FpX_factor_squarefree(T: FpX, p: bigint): FpX[] {
  const word = p < (1n << 64n);
  const f = word ? FpX_red(T, p) : T;
  if (f.length === 0) {
    if (word) throw new PariError('impossible inverse in Flx_divrem: Vecsmall([0])');
    throw new RangeError('FpX_factor_squarefree requires a nonzero polynomial');
  }
  return word ? Flx_factor_squarefree(f, p) : FpX_factor_Yun(f, p);
}


import { F2x_degree, F2x_add, F2x_deriv, F2x_divrem, F2x_gcd, F2x_sqrt,
  F2x_rem, F2x_valrem, F2x_Frobenius, F2x_matFrobenius, F2xq_sqr } from './F2x.js';
import { F2m_ker_sp } from './F2v.js';
import { random_F2x, random_Fl } from './random.js';

/** PARI's binary squarefree components, indexed by multiplicity minus one.
 * @see Deviation: PARI packed binary-polynomial kernels
 */
export function F2x_factor_squarefree(f:bigint):bigint[]{
  const n=F2x_degree(f),out=Array<bigint>(n+1).fill(1n);
  for(let q=1;;q*=2){
    let r=F2x_gcd(f,F2x_deriv(f));
    if(F2x_degree(r)===0){out[q-1]=f;break;}
    let t=F2x_divrem(f,r)[0];
    if(F2x_degree(t)>0){
      for(let j=1;;j++){
        const v=F2x_gcd(r,t),tv=F2x_divrem(t,v)[0];
        if(F2x_degree(tv)>0)out[j*q-1]=tv;
        if(F2x_degree(v)<=0)break;
        r=F2x_divrem(r,v)[0];t=v;
      }
      if(F2x_degree(r)===0)break;
    }
    f=F2x_sqrt(r);
  }
  let i=n;while(i>0&&F2x_degree(out[i-1]!)===0)i--;
  return out.slice(0,i);
}
function F2x_ddf_simple(T:bigint,XP:bigint):bigint[]{
  const n=F2x_degree(T);
  if(!n)return [];if(n===1)return [T];
  if(n<0)throw new RangeError('polynomial must be nonzero');
  let z=XP,remainder=T;
  const out=Array<bigint>(n).fill(1n);
  for(let j=1;j<=Math.floor(n/2);j++){
    const u=F2x_gcd(remainder,z^2n);
    if(F2x_degree(u)>0){out[j-1]=u;remainder=F2x_divrem(remainder,u)[0];}
    if(F2x_degree(remainder)===0)break;
    z=F2xq_sqr(z,remainder);
  }
  if(F2x_degree(remainder)>0)out[F2x_degree(remainder)-1]=remainder;
  return out;
}
/** DDF pairs contain common irreducible degree, not multiplicity. */
export function F2x_ddf(f:bigint):Array<[bigint,number]>{
  return F2x_ddf_simple(f,F2x_Frobenius(f)).flatMap((g,i)=>F2x_degree(g)>0?[[g,i+1]as [bigint,number]]:[]);
}
function F2x_edf_simple(f:bigint,XP:bigint,d:number):bigint[]{
  const n=F2x_degree(f);
  if(n/d===1)return [f];
  XP=F2x_rem(XP,f);
  let factor:bigint;
  for(;;){
    const g=random_F2x(n);
    let trace=g;
    for(let i=1;i<d;i++)trace=g^F2xq_sqr(trace,f);
    if(!trace)continue;
    factor=F2x_gcd(trace,f);
    const degree=F2x_degree(factor);
    if(degree>0&&degree<n)break;
  }
  const other=F2x_divrem(f,factor)[0];
  return [...F2x_edf_simple(factor,XP,d),...F2x_edf_simple(other,XP,d)];
}
function F2x_factor_Shoup(T:bigint):bigint[]{
  const XP=F2x_Frobenius(T),out:bigint[]=[];
  for(const [i,part]of F2x_ddf_simple(T,XP).entries()){
    const n=F2x_degree(part),d=i+1;
    if(!n)continue;
    if(n===d)out.push(part);else out.push(...F2x_edf_simple(part,XP,d));
  }
  return out;
}
function F2x_split_Berlekamp(u:bigint):bigint[]{
  const degree=F2x_degree(u);
  if(degree===1)return [u];
  if(degree===2)return u===6n?[2n,3n]:[u];
  const Q=F2x_matFrobenius(u).map((col,i)=>col^(1n<<BigInt(i)));
  const basis=F2m_ker_sp(Q,degree)as bigint[],d=basis.length;
  const factors=[u];let ir=0;
  const mark=(i:number)=>{if(i>ir)[factors[i],factors[ir]]=[factors[ir]!,factors[i]!];ir++;};
  while(factors.length<d){
    let polynomial:bigint;
    if(d===2)polynomial=basis[1]!;
    else{
      polynomial=random_Fl(2n);
      for(let i=1;i<d;i++)if(random_Fl(2n))polynomial^=basis[i]!;
    }
    for(let i=ir;i<factors.length&&factors.length<d;i++){
      const a=factors[i]!,n=F2x_degree(a);
      if(n===1)mark(i);
      else if(n===2){if(a===6n){factors[i]=2n;factors.push(3n);}mark(i);}
      else{
        let b=F2x_rem(polynomial,a);
        if(F2x_degree(b)<=0)continue;
        b=F2x_gcd(a,b);
        const m=F2x_degree(b);
        if(m>0&&m<n){factors.push(F2x_divrem(a,b)[0]);factors[i]=b;}
      }
    }
  }
  return factors;
}
/** Native binary factor dispatch: direct degrees <=2, Cantor through 20, Berlekamp above.
 * Returns packed factors with multiplicities, sorted as in PARI.
 * @see Deviation: PARI packed binary-polynomial kernels
 */
export function F2x_factor(f:bigint):Array<[bigint,number]>{
  const degree=F2x_degree(f);
  if(degree<=2){
    if(degree===0)return [];
    if(degree===-1||degree===1||f===7n)return [[f,1]];
    return f===6n?[[2n,1],[3n,1]]:[[2n+(f&1n),2]];
  }
  const out:Array<[bigint,number]>=[];
  if(degree<=20){
    for(const [i,part]of F2x_factor_squarefree(f).entries())
      if(F2x_degree(part)>0)for(const g of F2x_factor_Shoup(part))out.push([g,i+1]);
  }else{
    const [valuation,reduced]=F2x_valrem(f);
    if(valuation)out.push([2n,Number(valuation)]);
    for(const [i,part]of F2x_factor_squarefree(reduced).entries())
      if(F2x_degree(part)>0)for(const g of F2x_split_Berlekamp(part))out.push([g,i+1]);
  }
  out.sort(([a],[b])=>a<b?-1:a>b?1:0);
  return out;
}

import { polynomialQuotient } from './_polynomial_quotient.js';
import { polynomialDistinctDegrees } from './_polynomial_ddf.js';
import { trimPolynomial } from './_polynomial_packing.js';

export type PolynomialDegreeFactor = [FpX, number];
function ddfModulus(p: bigint): void {
  if (p < 2n)
    throw new RangeError('distinct-degree factorization requires a modulus greater than one');
}
function wordDdfPolynomial(T: FpX, p: bigint): FpX {
  if (p < 1n || p >= 1n << 64n) throw new RangeError('modulus must be a positive word integer');
  if (T.some((c) => c < 0n || c >= p))
    throw new RangeError('word polynomial coefficients must be reduced');
  ddfModulus(p);
  return trimPolynomial(T);
}
function distinctDegreeParts(T: FpX, p: bigint, word: boolean): FpX[] {
  const context = polynomialQuotient(T, p, word);
  // Public native wrappers compute Frobenius before Shoup's degree shortcuts.
  return polynomialDistinctDegrees(context, context.power([0n, 1n]));
}
function degreePairs(parts: FpX[]): PolynomialDegreeFactor[] {
  return parts.flatMap((f, i) => (f.length === 1 ? [] : [[f, i + 1] as PolynomialDegreeFactor]));
}
/** PARI Flx_ddf, with unscaled word GCD factors and increasing common degrees.
 * @see Deviation: PARI distinct-degree factor adapters
 */
export function Flx_ddf(T: FpX, p: bigint): PolynomialDegreeFactor[] {
  return degreePairs(distinctDegreeParts(wordDdfPolynomial(T, p), p, true));
}
/** PARI FpX_ddf native dispatch; pairs match the existing F2x_ddf representation.
 * @see Deviation: PARI distinct-degree factor adapters
 */
export function FpX_ddf(T: FpX, p: bigint): PolynomialDegreeFactor[] {
  ddfModulus(p);
  let f = FpX_red(T, p);
  if (p === 2n) {
    let packed = 0n;
    for (let i = f.length - 1; i >= 0; i--) packed = (packed << 1n) | f[i]!;
    return F2x_ddf(packed).map(([a, degree]) => {
      const coefficients: FpX = [];
      while (a) {
        coefficients.push(a & 1n);
        a >>= 1n;
      }
      return [coefficients, degree];
    });
  }
  if (f.length > 1) f = FpX_normalize(f, p);
  return p < 1n << 64n ? Flx_ddf(f, p) : degreePairs(distinctDegreeParts(f, p, false));
}
function degreeCounts(parts: FpX[]): { D: number[]; nb: number } {
  const D = [0, ...parts.map((f, i) => Math.trunc((f.length - 1) / (i + 1)))];
  return { D, nb: D.reduce((a, b) => a + b, 0) };
}
/** PARI Flx_nbfact_by_degree; D[0] is the port's unused zero entry.
 * @see Deviation: PARI distinct-degree factor adapters
 */
export function Flx_nbfact_by_degree(T: FpX, p: bigint): { D: number[]; nb: number } {
  const f = wordDdfPolynomial(T, p);
  if (!f.length) throw new PariError('overflow in lg().');
  return degreeCounts(distinctDegreeParts(f, p, true));
}
/** PARI FpX_nbfact uses the FpX Shoup route directly, including word primes.
 * @see Deviation: PARI distinct-degree factor adapters
 */
export function FpX_nbfact(T: FpX, p: bigint): number {
  ddfModulus(p);
  return degreeCounts(distinctDegreeParts(trimPolynomial(T), p, false)).nb;
}
/** Existing generic count convenience, with native word-library delegation.
 * @see Deviation: PARI distinct-degree factor adapters
 */
export function _FpX_nbfact_by_degree(T: FpX, p: bigint): { D: number[]; nb: number } {
  ddfModulus(p);
  if (p < 1n << 64n) return Flx_nbfact_by_degree(FpX_red(T, p), p);
  return degreeCounts(distinctDegreeParts(trimPolynomial(T), p, false));
}


import { Flx_gcd, Flx_divrem, Flx_normalize } from './Flx.js';
import { Fp_sqrt, Fl_sqrt, kronecker } from './ff.js';
import { polynomialFactorShoup } from './_polynomial_factor.js';

export type PolynomialFactor = [FpX, number];

function factorNormalized(f: FpX, p: bigint, word: boolean): PolynomialFactor[] {
  if (p === 2n) {
    let packed = 0n;
    for (let i = 0; i < f.length; i++) packed |= f[i]! << BigInt(i);
    return F2x_factor(packed).map(([g, e]) => {
      const a: FpX = [];
      while (g) {
        a.push(g & 1n);
        g >>= 1n;
      }
      return [a, e];
    });
  }
  const d = f.length - 1;
  if (d < 0) return [[[], 1]];
  if (!d) return [];
  if (d === 1) return [[f, 1]];
  if (d === 2) {
    const mod = (x: bigint) => ((x % p) + p) % p;
    const b = f[1]!,
      c = f[0]!;
    const discriminant = word ? mod(b * b - 4n * c) : (b * b - 4n * c) % p;
    if (kronecker(discriminant, p) === -1) return [[f, 1]];
    const s = word ? Fl_sqrt(discriminant, p) : Fp_sqrt(discriminant, p);
    if (s === null) return [[f, 1]];
    let root = mod(s - b);
    root = (root + (root & 1n ? p : 0n)) >> 1n;
    const first = mod(-root),
      second = mod(b + root);
    return first === second
      ? [[[first, 1n], 2]]
      : [
          [[first < second ? first : second, 1n], 1],
          [[first < second ? second : first, 1n], 1],
        ];
  }
  const out: PolynomialFactor[] = [];
  for (const [i, part] of FpX_factor_squarefree(f, p).entries())
    if (part.length > 1)
      for (const factor of polynomialFactorShoup(part, p, word)) out.push([factor, i + 1]);
  return out.sort(([a], [b]) => {
    if (a.length !== b.length) return a.length - b.length;
    for (let i = a.length - 1; i >= 0; i--) if (a[i] !== b[i]) return a[i]! < b[i]! ? -1 : 1;
    return 0;
  });
}

/** Native irreducible factor/multiplicity pairs, consuming the global PARI state.
 * @see Deviation: PARI full polynomial factorization adapters
 */
export function FpX_factor(f: FpX, p: bigint): PolynomialFactor[] {
  if (p < 2n) throw new RangeError('polynomial factorization requires a modulus greater than one');
  const word = p < 1n << 64n;
  f = FpX_red(f, p);
  if (p !== 2n && f.length > 1) f = word ? Flx_normalize(f, p) : FpX_normalize(f, p);
  return factorNormalized(f, p, word);
}

/** Native reduced-word factorization, including normalization before degree checks.
 * @see Deviation: PARI full polynomial factorization adapters
 */
export function Flx_factor(f: FpX, p: bigint): PolynomialFactor[] {
  if (p < 1n || p >= 1n << 64n) throw new RangeError('modulus must be a positive word integer');
  if (f.some((c) => c < 0n || c >= p))
    throw new RangeError('word polynomial coefficients must be reduced');
  if (p < 2n) throw new RangeError('polynomial factorization requires a modulus greater than one');
  return factorNormalized(Flx_normalize(trimPolynomial(f), p), p, true);
}

import { polynomialRoots, polynomialTotallySplit, polynomialWordRootCount } from './_polynomial_roots.js';
/** Native reduced roots, including distinct word and large-prime ordering.
 * @see Deviation: PARI polynomial root adapters
 */
export function FpX_roots(f: FpX, p: bigint): bigint[] {
  return polynomialRoots(f, p);
}
/** Native reduced-word roots, with zero-polynomial errors.
 * @see Deviation: PARI polynomial root adapters
 */
export function Flx_roots(f: FpX, p: bigint): bigint[] {
  return polynomialRoots(f, p, true);
}
/** Native Frobenius predicate, preserving raw degree shortcuts.
 * @see Deviation: PARI polynomial root adapters
 */
export function FpX_is_totally_split(f: FpX, p: bigint): boolean {
  return polynomialTotallySplit(f, p);
}
/** Native word Frobenius predicate; repeated factors are not totally split.
 * @see Deviation: PARI polynomial root adapters
 */
export function Flx_is_totally_split(f: FpX, p: bigint): boolean {
  return polynomialTotallySplit(f, p, true);
}
/** Count distinct roots, with PARI's quadratic discriminant shortcut.
 * @see Deviation: PARI polynomial root adapters
 */
export function Flx_nbroots(f: FpX, p: bigint): number {
  return polynomialWordRootCount(f, p);
}

import { polynomialQuotientPower } from './_polynomial_quotient_power.js';
/** PARI FpX_factor.c:738, preserving the small-degree shortcut. */
export function FpX_split_part(f: FpX, p: bigint): FpX {
  f = trimPolynomial(f);
  if (f.length <= 2) return f;
  f = FpX_red(f, p);
  return FpX_gcd(FpX_sub(polynomialQuotientPower([0n, 1n], p, f, p), [0n, 1n], p), f, p);
}
