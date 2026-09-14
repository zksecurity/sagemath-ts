import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { getrand, setrand } from '../../../packages/parigp-ts/src/random.js';
export function nf_pari_factor_constructor(degree: bigint, kind: bigint, scale: bigint): string {
  const n = Number(degree),
    c = Array<bigint>(n + 1).fill(0n);
  c[n] = kind === 5n ? 5n : 1n;
  if (kind === 0n || kind === 5n) c[0] = -2n;
  else if (kind === 1n) c[0] = -1n;
  else if (kind === 2n) {
    c[0] = 2n;
    c[1] = -2n;
    c[n - 1] = c[n - 1]! - 1n;
  } else if (kind === 3n) {
    c[n - 2] = 1n;
    c[n - 1] = -2n;
  } else {
    c[0] = 1n;
    c[1] = 1n;
  }
  const f = new RationalPolynomial(c.map((v) => new Rational(v * scale, 7n))),
    result: unknown[] = [f.isIrreducible()];
  try {
    result.push(['degree', String(new NumberField(f, 'a').degree())]);
  } catch (e) {
    result.push([(e as Error).name, (e as Error).message]);
  }
  return JSON.stringify(result);
}
export function nf_irreducible_cache(P: bigint[]): string {
  const f = RationalPolynomial.fromBigInts(P),
    saved = getrand();
  try {
    setrand(1n);
    const first = f.isIrreducible();
    setrand(2n);
    const before = getrand(),
      second = f.isIrreducible();
    return JSON.stringify([first, second, before === getrand()]);
  } finally {
    setrand(saved);
  }
}
