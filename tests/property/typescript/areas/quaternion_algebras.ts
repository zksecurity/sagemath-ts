import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.ts';
/**
 * sagemath-ts side of the live quaternion-algebra differential area.
 */

import { QuaternionAlgebra } from '../../../../packages/sagemath-ts/src/algebras/quatalg/index.js';

type Printable = { toString(): string };
const vals = (xs: Printable[]): string => xs.map(String).join(',');
const bool = (value: boolean): string => (value ? 'True' : 'False');

function quat_element(a: bigint, b: bigint, coeffs: bigint[]): string {
  const A = QuaternionAlgebra(a, b);
  const x = A.__call__(coeffs);
  const inv = x.is_zero() ? '-' : vals(x.inverse().list());
  const charpoly = vals(x.reduced_characteristic_polynomial().coeffs);
  return `x=${vals(x.list())} conj=${vals(x.conjugate().list())} tr=${x.reduced_trace()} norm=${x.reduced_norm()} inv=${inv} cp=${charpoly}`;
}

function quat_product(a: bigint, b: bigint, left: bigint[], right: bigint[]): string {
  const A = QuaternionAlgebra(a, b);
  const x = A.__call__(left);
  const y = A.__call__(right);
  const z = x.mul(y);
  return `xy=${vals(z.list())} pair=${x.pair(y)} norm=${z.reduced_norm()} anti=${bool(
    z.conjugate().eq(y.conjugate().mul(x.conjugate()))
  )}`;
}

function quat_algebra(a: bigint, b: bigint): string {
  const A = QuaternionAlgebra(a, b);
  return `inv=${vals(A.invariants())} disc=${A.discriminant()} ram=${A.ramified_primes().join(
    ','
  )} definite=${bool(A.is_definite())}`;
}

function quat_lattice_reduction(
  a: bigint,
  b: bigint,
  den: bigint,
  flat: bigint[],
  op: bigint
): string {
  const A = QuaternionAlgebra(a, b);
  const basis = Array.from({ length: 4 }, (_, i) =>
    A.__call__(flat.slice(i * 4, i * 4 + 4)).scalar_mul(new Rational(1n, den))
  );
  const I = A.ideal(basis);
  if (op === 0n)
    return I.reduced_basis()
      .map((x) => vals(x.list()))
      .join(';');
  if (op === 1n) return vals(I.minimal_element().list());
  if (op === 2n)
    return Array.from({ length: 4 }, (_, i) =>
      Array.from({ length: 4 }, (_, j) => String(I.quadratic_form().get(i, j)))
    )
      .flat()
      .join(',');
  return I.theta_series_vector(12).join(',');
}

function quat_order_isomorphism(scale: bigint): string {
  const A = QuaternionAlgebra(-1n, -19n),
    [i, j, k] = A.gens();
  const O0 = A.quaternion_order([
    A.one(),
    i,
    i.add(j).scalar_mul('1/2'),
    A.one().add(k).scalar_mul('1/2'),
  ]);
  let O1 = A.quaternion_order([
    A.one(),
    i.scalar_mul(667n),
    A.__call__('1/2').add(j.scalar_mul('1/2')).add(i.scalar_mul(9n)),
    i.scalar_mul('222075/2').add(j.scalar_mul(333n)).add(k.scalar_mul('1/2')).scalar_mul('1/667'),
  ]);
  if (scale !== 0n) {
    const a = A.__call__([1n, scale, scale + 1n, 2n * scale - 1n]);
    O1 = A.quaternion_order(O1.basis().map((x) => a.inverse().mul(x).mul(a)));
  }
  const gamma = O0.isomorphism_to(O1, { conjugator: true }) as ReturnType<typeof A.__call__>;
  return [gamma, ...A.gens().map((x) => gamma.inverse().mul(x).mul(gamma))]
    .map((x) => vals(x.list()))
    .join(';');
}

export const functions = {
  quat_order_isomorphism,
  quat_lattice_reduction,
  quat_element,
  quat_product,
  quat_algebra,
};
