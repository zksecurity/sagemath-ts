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

export const functions = {
  quat_element,
  quat_product,
  quat_algebra,
};
