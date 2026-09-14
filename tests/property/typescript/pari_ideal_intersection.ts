import { idealintersect } from '../../../packages/parigp-ts/src/base4.js';
import { setrand, randomi } from '../../../packages/parigp-ts/src/random.js';
import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
const fields = new Map<number, ReturnType<NumberField['_pari_ideal_data']>>();
export function pari_ideal_intersection(
  n: bigint,
  seed: bigint,
  iDen: bigint,
  jDen: bigint,
  a: bigint[],
  b: bigint[]
): string {
  const degree = Number(n);
  let nf = fields.get(degree);
  if (!nf) {
    const coefficients = Array<bigint>(degree + 1).fill(0n);
    coefficients[0] = -2n;
    coefficients[degree] = 1n;
    nf = new NumberField(RationalPolynomial.fromBigInts(coefficients), 'a')._pari_ideal_data();
    fields.set(degree, nf);
  }
  const columns = (flat: bigint[]) =>
    Array.from({ length: flat.length / degree }, (_, j) =>
      Array.from({ length: degree }, (_, i) => flat[i * (flat.length / degree) + j]!)
    );
  setrand(seed);
  const [H, den] = idealintersect(nf, columns(a), iDen, columns(b), jDen);
  return JSON.stringify([H, den, randomi(1n << 128n)], (_, v) =>
    typeof v === 'bigint' ? String(v) : v
  );
}
