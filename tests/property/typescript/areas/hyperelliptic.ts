/**
 * sagemath-ts side of the live hyperelliptic-curve differential area.
 */

import { GF } from '../../../../packages/sagemath-ts/src/rings/finite_rings/index.js';
import type { FiniteFieldElement } from '../../../../packages/sagemath-ts/src/rings/finite_rings/index.js';
import { PolynomialRing } from '../../../../packages/sagemath-ts/src/rings/polynomial/index.js';
import { HyperellipticCurve } from '../../../../packages/sagemath-ts/src/schemes/hyperelliptic_curves/index.js';

type AnyCurve = {
  genus(): number;
  count_points(n: number): bigint[];
  frobenius_polynomial(): bigint[];
  jacobian(): { cardinality(): bigint };
  Cartier_matrix(): FiniteFieldElement[][];
  Hasse_Witt(): FiniteFieldElement[][];
  a_number(): number;
  p_rank(): number;
};

const ints = (xs: Array<bigint | FiniteFieldElement>): string =>
  xs.map((x) => (typeof x === 'bigint' ? x : x.value).toString()).join(',');

function poly(p: bigint, coeffs: bigint[]) {
  const K = GF(p);
  const R = new PolynomialRing<FiniteFieldElement>(K, 'x');
  return R.__call__(coeffs.map((c) => K.__call__(c)));
}

function hyp_summary(p: bigint, fCoeffs: bigint[], hCoeffs: bigint[], extensions: bigint): string {
  const f = poly(p, fCoeffs);
  const h = hCoeffs.length === 0 ? null : poly(p, hCoeffs);
  const H = HyperellipticCurve(f, h) as unknown as AnyCurve;
  return `g=${H.genus()} counts=${ints(H.count_points(Number(extensions)))} frob=${ints(
    H.frobenius_polynomial()
  )} jac=${H.jacobian().cardinality()}`;
}

function hyp_cartier(p: bigint, fCoeffs: bigint[]): string {
  const H = HyperellipticCurve(poly(p, fCoeffs)) as unknown as AnyCurve;
  const C = H.Cartier_matrix().map(ints).join(';');
  const HW = H.Hasse_Witt().map(ints).join(';');
  return `C=${C} HW=${HW} a=${H.a_number()} p=${H.p_rank()}`;
}

export const functions = {
  hyp_summary,
  hyp_cartier,
};
