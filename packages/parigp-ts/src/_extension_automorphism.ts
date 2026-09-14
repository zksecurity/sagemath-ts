/** Native extension automorphism powers, additive traces and multiplicative sums.
 * @see Deviation: PARI extension automorphism adapters
 */
import { polynomialQuotient } from './_polynomial_quotient.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { extensionGetRed, type ExtensionModulus } from './_extension_division.js';
import { extensionQuotient } from './_extension_quotient.js';
import { extensionComposition } from './_extension_composition.js';
import { coefficientSubstitution } from './_coefficient_composition.js';
import { extensionField, trimExtension } from './_extension_field.js';
import { validateExtensionInputs, type ExtensionPolynomial as P } from './_extension_polynomial.js';
import { F2x_degree, F2xq_powers, F2x_F2xqV_eval } from './F2x.js';
import { brent_kung_optpow } from './RgX.js';
import { gen_powu_i } from './bb_group.js';

type Inner = bigint | bigint[];
type Aut = [Inner | P, P] | [Inner, P, P];

export function extensionAutomorphism(
  mode: 0 | 1 | 2,
  op: 0 | 1 | 2,
  aut: Aut,
  n: bigint,
  S: ExtensionModulus,
  T: bigint[] | bigint,
  p: bigint,
  providedInner?: bigint[]
): Aut {
  const unsigned = mode === 1 && op === 1;
  if (n === 0n || n < (unsigned ? 1n : -(1n << 63n)) || n >= (unsigned ? 1n << 64n : 1n << 63n))
    throw new RangeError(
      unsigned
        ? 'automorphism exponent must be a nonzero unsigned word'
        : 'automorphism exponent must be a nonzero signed word'
    );
  // The original signed-long entries pass their bit pattern to gen_powu_i.
  n = BigInt.asUintN(64, n);
  const canonical = (a: P) =>
    trimExtension(a.map((c) => (Array.isArray(c) ? trimPolynomial(c) : c)));
  const plainTrace = op === 1 && mode !== 2;
  if (plainTrace) aut = [canonical(aut[0] as P), canonical(aut[1])];
  else {
    const phi = mode === 2 ? (aut[0] as bigint) : trimPolynomial(aut[0] as bigint[]);
    aut = aut.length === 3 ? [phi, canonical(aut[1]), canonical(aut[2])] : [phi, canonical(aut[1])];
  }
  if (Array.isArray(T)) T = trimPolynomial(T);
  if (plainTrace) validateExtensionInputs(mode, p, T, aut[0] as P, aut[1]);
  else {
    validateExtensionInputs(mode, p, T, [aut[0] as Inner]);
    validateExtensionInputs(mode, p, T, aut[1], aut[2] ?? []);
  }
  validateExtensionInputs(mode, p, T, Array.isArray(S) ? S : S.polynomial);
  const inner =
    mode === 2 ? undefined : polynomialQuotient(T as bigint[], p, mode === 1, providedInner);
  const inverse = inner?.inverse;
  S = extensionGetRed(mode, S, T, p, inverse);
  const f = extensionField(mode, T, p, inverse);
  const add = (a: P, b: P) =>
    trimExtension(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) =>
        i < a.length ? (i < b.length ? f.add(a[i]!, b[i]!) : a[i]!) : b[i]!
      )
    );
  const outerPowers = (x: P, count: number) =>
    extensionQuotient(mode, 6, p, BigInt(count), T, x, [], S, inverse) as P[];
  const outerEval = (Q: P, x: P, V: P[], table: boolean) =>
    extensionComposition(mode, table ? 1 : 0, p, T, Q, x, V, S, inverse);
  const innerPowers = (phi: Inner, count: number) =>
    mode === 2
      ? F2xq_powers(phi as bigint, count, T as bigint)
      : inner!.powers(phi as bigint[], count);
  const innerEval = (phi: Inner, V: bigint[] | bigint[][]): Inner =>
    mode === 2
      ? F2x_F2xqV_eval(phi as bigint, V as bigint[], T as bigint)
      : inner!.evaluate(phi as bigint[], V as bigint[][]);
  const substitute = (Q: P, V: bigint[] | bigint[][]) =>
    coefficientSubstitution(mode, 1, p, T, Q, mode === 2 ? 0n : [], V, inverse);
  const degreeT = mode === 2 ? F2x_degree(T as bigint) : (T as bigint[]).length - 1;
  const multiply = (x: Aut, y: Aut): Aut => {
    if (plainTrace) {
      const [S1, a1] = x as [P, P],
        [S2, a2] = y as [P, P];
      const count = brent_kung_optpow(Math.max(S1.length, a1.length) - 1, 2, 1);
      const V = outerPowers(S2, count);
      const S3 = outerEval(S1, [], V, true),
        aS = outerEval(a1, [], V, true);
      return [S3, add(aS, a2)];
    }
    const [phi1, S1, a1] = x,
      [phi2, S2, a2] = y;
    const count = brent_kung_optpow(degreeT - 1, S1.length + (op === 0 ? 0 : a1!.length) + 1, 1);
    const V = innerPowers(phi2 as Inner, count);
    const phi3 = innerEval(phi1 as Inner, V),
      Sphi = substitute(S1, V);
    if (op === 0) return [phi3, outerEval(Sphi, S2, [], false)];
    const aphi = substitute(a1!, V);
    const count2 = brent_kung_optpow(Math.max(Sphi.length, aphi.length) - 1, 2, 1);
    const W = outerPowers(S2, count2);
    const S3 = outerEval(Sphi, [], W, true),
      aS = outerEval(aphi, [], W, true);
    const a3 =
      mode === 2 ? add(aS, a2!) : (extensionQuotient(mode, 0, p, 0n, T, aS, a2!, S, inverse) as P);
    return [phi3, S3, a3];
  };
  const result = gen_powu_i(aut, n, (a) => multiply(a, a), multiply);
  return result.map((c) =>
    Array.isArray(c) ? c.map((x) => (Array.isArray(x) ? x.slice() : x)) : c
  ) as Aut;
}
