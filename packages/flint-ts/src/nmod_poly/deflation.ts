/** nmod_poly/deflation.c: gcd of supported exponents, retaining the native skip schedule. */
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
export function nmod_poly_deflation(f: readonly bigint[], p: bigint): number {
  const a = k.normalized(f, p);
  if (a.length <= 1) return a.length;
  const gcd = (a: number, b: number): number => {
    while (b) [a, b] = [b, a % b];
    return a;
  };
  let coeff = 1;
  while (!a[coeff]) coeff++;
  let deflation = gcd(a.length - 1, coeff);
  while (deflation > 1 && coeff + deflation < a.length) {
    let i = 0;
    for (; i < deflation - 1; i++) {
      coeff++;
      if (a[coeff]) deflation = gcd(coeff, deflation);
    }
    if (i === deflation - 1) coeff++;
  }
  return deflation;
}
