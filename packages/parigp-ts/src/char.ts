/** PARI char.c:568-735: znstar's cyclic factors and Conrey generators (flag 0). */
import { Z_factor } from './ifactor.js';
import { Fp_pow, Fp_inv, xgcd } from './ff.js';
import { residue } from './_polynomial_division.js';
import { pgener_Zp } from './arith1.js';
export function znstar(value: bigint): [bigint, bigint[], bigint[]] {
  if (value === 0n) return [2n, [2n], [-1n]];
  const N = value < 0n ? -value : value;
  if (N <= 2n) return [1n, [], []];
  const factors = Z_factor(N),
    cyc: bigint[] = [],
    gen: bigint[] = [],
    mod: bigint[] = [];
  for (const [p, e] of factors) {
    const Q = p ** e;
    if (p === 2n) {
      if (e === 1n) continue;
      cyc.push(2n);
      gen.push(Q - 1n);
      mod.push(Q);
      if (e >= 3n) {
        cyc.push(1n << (e - 2n));
        gen.push(5n);
        mod.push(Q);
      }
    } else {
      cyc.push((p - 1n) * p ** (e - 1n));
      gen.push(pgener_Zp(p));
      mod.push(Q);
    }
  }
  if (factors.length > 1)
    for (let i = 0; i < gen.length; i++) {
      const Q = mod[i]!,
        g = gen[i]!,
        inv = Fp_inv(Q, N / Q);
      gen[i] = residue(g + (1n - g) * inv * Q, N);
    }
  for (let i = cyc.length - 1; i >= 1; i--) {
    let ci = cyc[i]!,
      gi = gen[i]!;
    for (let j = i - 1; j >= 0; j--) {
      const cj = cyc[j]!;
      let [d, , v] = xgcd(ci, cj);
      if (ci === d) continue;
      if (cj === d) {
        [gen[i], gen[j]] = [gen[j]!, gen[i]!];
        gi = gen[i]!;
        [cyc[i], cyc[j]] = [cyc[j]!, cyc[i]!];
        ci = cyc[i]!;
        continue;
      }
      const qj = cj / d;
      cyc[j] = ci * qj;
      cyc[i] = d;
      const gj = residue(gen[j]! * gi, N);
      v = -v;
      if (v < 0n) v = residue(v, ci);
      gen[i] = gi = residue(gi * Fp_pow(gj, qj * v, N), N);
      gen[j] = gj;
      ci = d;
      if (ci === 2n) break;
    }
  }
  return [cyc.reduce((a, b) => a * b, 1n), cyc, gen];
}
