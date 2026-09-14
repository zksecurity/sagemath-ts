/** PARI subcyclo.c:330-374, small znstar/HNF interface. */
import { Fp_pow } from './ff.js';
export interface SmallZnStar {
  modulus: bigint;
  orders: bigint[];
  generators: bigint[];
}
/** Zero-indexed storage adapter for a flag-zero znstar result. */
export function znstar_small(modulus: bigint, Z: [bigint, bigint[], bigint[]]): SmallZnStar {
  return { modulus: modulus < 0n ? -modulus : modulus, orders: Z[1], generators: Z[2] };
}
/** Native HNF generator evaluation, relative-order cosets and sorted elements. */
export function znstar_hnf_elts(Z: SmallZnStar, H: bigint[][]): number[] {
  const N = Z.modulus,
    gen = Z.generators;
  const V = H.map((c) => c.reduce((g, e, i) => (e ? (g * Fp_pow(gen[i]!, e, N)) % N : g), 1n));
  let elts = [1n],
    bits = new Uint8Array(Number(N));
  bits[1] = 1;
  for (const v of V) {
    let g = v,
      o = 0;
    while (!bits[Number(g)]) {
      g = (g * v) % N;
      o++;
    }
    if (!o) continue;
    const old = elts.length,
      count = old * o;
    for (let k = 0; k < count; k++) elts.push((elts[k]! * v) % N);
    bits = new Uint8Array(Number(N));
    for (const x of elts) bits[Number(x)] = 1;
  }
  return elts.map(Number).sort((a, b) => a - b);
}
