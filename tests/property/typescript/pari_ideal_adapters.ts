import { ZM_hnfcenter } from '../../../packages/parigp-ts/src/hnf_snf.js';
import { diviiround } from '../../../packages/parigp-ts/src/gen3.js';
import { idealmul } from '../../../packages/parigp-ts/src/base4.js';
import { setrand, randomi } from '../../../packages/parigp-ts/src/random.js';
import { nfbasis } from '../../../packages/sagemath-ts/src/rings/number_field/pari_nf.js';
import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
const fmt = (value: unknown) =>
  JSON.stringify(value, (_, v) => (typeof v === 'bigint' ? String(v) : v));
const columns = (flat: bigint[], n: number) =>
  Array.from({ length: n }, (_, j) => Array.from({ length: n }, (_, i) => flat[i * n + j]!));
export function pari_hnfcenter(n: bigint, flat: bigint[]): string {
  return fmt(ZM_hnfcenter(columns(flat, Number(n))));
}
export function pari_diviiround(a: bigint, b: bigint): string {
  return fmt(diviiround(a, b));
}
export function pari_centered_basis(cs: bigint[]): string {
  const { basis, den } = nfbasis(cs);
  return JSON.stringify(basis.map((row) => row.map((v) => String(new Rational(v, den)))));
}
const fields = new Map<number, ReturnType<NumberField['_pari_ideal_data']>>();
export function pari_idealmul(
  n: bigint,
  iDen: bigint,
  jDen: bigint,
  ni: bigint,
  nj: bigint,
  a: bigint[],
  b: bigint[]
): string {
  const degree = Number(n);
  let nf = fields.get(degree);
  if (!nf) {
    const cs = Array<bigint>(degree + 1).fill(0n);
    cs[0] = -2n;
    cs[degree] = 1n;
    nf = new NumberField(RationalPolynomial.fromBigInts(cs), 'a')._pari_ideal_data();
    fields.set(degree, nf);
  }
  setrand(42n);
  const [H, den] = idealmul(nf, columns(a, Number(ni)), iDen, columns(b, Number(nj)), jDen);
  return fmt([H, den, randomi(1n << 128n)]);
}
