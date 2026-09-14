/** The exact ideal-data part of PARI base1.c nfmaxord_to_nf.
 * @see Deviation: Number-field ideal backend adapters
 * @see Deviation: Number-field ideal valuation adapters
 */
export interface NfIdealData {
  polynomial: bigint[];
  basis: bigint[][];
  basisDenominator: bigint;
  discriminant: bigint;
  index?: bigint;
  multiplication: bigint[][][];
  traceInverse: bigint[][];
  traceDenominator: bigint;
  codifferent: bigint[][];
  codifferentTwo: [bigint, bigint[]];
  different: bigint[][];
}

import { integerMatrixInverse } from './_matrix_inverse.js';
import { ZM_hnfmodid } from './hnf_snf.js';
import { idealHNF_inv_Z, mat_ideal_two_elt, zk_multable } from './base4.js';

/** Exact ideal components of nfmaxord_to_nf, supplied the existing maximal basis. */
export function nfmaxord_ideal_data(
  T: bigint[],
  basis: bigint[][],
  den: bigint,
  mul: bigint[][][],
  index: bigint,
  disc: bigint
): NfIdealData {
  const n = mul.length,
    sym = [BigInt(n)];
  // polsym(T,n-1), for the monic integral polynomial produced by nfmaxord.
  for (let k = 1; k < n; k++) {
    let s = BigInt(k) * T[n - k]!;
    for (let i = 1; i < k; i++) s += sym[k - i]! * T[n - i]!;
    sym.push(-s);
  }
  const trace = basis.map((row) => row.reduce((s, c, j) => s + c * sym[j]!, 0n) / den);
  const Tr = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => mul[i]![j]!.reduce((s, c, k) => s + c * trace[k]!, 0n))
  );
  const inverse = integerMatrixInverse(Tr);
  if (!inverse || inverse[1] === 0n) throw new Error('singular number-field trace pairing');
  const [raw, d] = inverse,
    columns = raw[0]!.map((_, i) => raw.map((row) => row[i]!));
  const A = ZM_hnfmodid(columns, Array<bigint>(n).fill(d)),
    two = mat_ideal_two_elt(mul, A);
  const nf: NfIdealData = {
    polynomial: T.slice(),
    basis: basis.map((c) => c.slice()),
    basisDenominator: den,
    discriminant: disc,
    index,
    multiplication: mul,
    traceInverse: columns,
    traceDenominator: d,
    codifferent: A,
    codifferentTwo: two,
    different: [],
  };
  if (index === 1n) {
    const binv = integerMatrixInverse(basis)!;
    const coordinates = Array.from({ length: n }, (_, j) => {
      let v = 0n;
      for (let i = 0; i < n; i++) v += BigInt(i + 1) * T[i + 1]! * binv[0][i]![j]! * den;
      return v / binv[1];
    });
    nf.different = ZM_hnfmodid(
      zk_multable(mul, coordinates),
      Array<bigint>(n).fill(disc < 0n ? -disc : disc)
    );
  } else {
    const c = d / A[0]![0]!;
    nf.different = idealHNF_inv_Z(nf, A).map((col) => col.map((x) => x * c));
  }
  return nf;
}
