import { QQ } from '../../../../packages/sagemath-ts/src/rings/rational_field.js';
import { cardinality_of } from '../../../../packages/sagemath-ts/src/schemes/hyperelliptic_curves/field_ops.js';
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
  hyp_base_cardinality,
  hyp_cartier,
};

function hyp_base_cardinality(p: bigint, f: bigint[], h: bigint[], op: bigint): string {
  const k: any = p === 0n ? QQ : GF(p);
  if (op === 0n) return String(cardinality_of(k));
  const R = new PolynomialRing<any>(k, 'x'),
    H: any = HyperellipticCurve(R.__call__(f), R.__call__(h));
  return String(H) + '|' + H.genus() + '|' + String(H.base_ring());
}

import { GFpn } from '../../../../packages/sagemath-ts/src/rings/finite_rings/index.js';
import { Polynomial } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_element.js';
import {
  sqrt_all_of,
  compare_elements,
} from '../../../../packages/sagemath-ts/src/schemes/hyperelliptic_curves/field_ops.js';
import { cantor_reduction } from '../../../../packages/sagemath-ts/src/schemes/hyperelliptic_curves/jacobian_morphism.js';
import { setrand, getrand } from '@sagemath-ts/parigp-ts';
Object.assign(functions, {
  hyp_root_callers(
    kind: bigint,
    coefficients: bigint[],
    coordinate: bigint,
    operation: bigint,
    seed: bigint
  ) {
    const trace: any[] = [];
    let result: any;
    try {
      const fields: Record<number, () => any> = {
        0: () => QQ,
        2: () => GF(3n),
        3: () => GF(5n),
        4: () => GF(7n),
        5: () => GFpn(3n, 2, [1, 0], 'a'),
        8: () => GFpn(2n, 2, [1, 1], 'a'),
        10: () => GFpn(2n, 3, [1, 1, 0], 'a'),
        11: () => GFpn(3n, 3, [1, 2, 0], 'a'),
      };
      const K = fields[Number(kind)]!(),
        R = new PolynomialRing<any>(K, 'x');
      const decode = (n: bigint) =>
        K.fromInteger ? K.fromInteger(((n % K.order) + K.order) % K.order) : K.__call__(n);
      const f = R.__call__(coefficients.map(decode));
      const proto = Object.getPrototypeOf(K.zero()),
        sqrt = proto.sqrt,
        square = proto.is_square;
      const poly = Polynomial.prototype,
        roots = poly.roots;
      let depth = 0;
      setrand(seed);
      try {
        proto.sqrt = function (options: any) {
          if (!depth)
            trace.push(['sqrt', String(this), options?.all ?? false, options?.extend ?? true]);
          depth++;
          try {
            return sqrt.call(this, options);
          } finally {
            depth--;
          }
        };
        proto.is_square = function () {
          if (!depth) trace.push(['is_square', String(this)]);
          depth++;
          try {
            return square.call(this);
          } finally {
            depth--;
          }
        };
        poly.roots = function (options?: any): any {
          if (!depth)
            trace.push(['roots', this.coeffs.map(String), options?.multiplicities ?? true]);
          depth++;
          try {
            return roots.call(this, options);
          } finally {
            depth--;
          }
        };
        if (operation === 0n) {
          result = { value: sqrt_all_of(K, decode(coordinate)).map(String) };
          if (K.fromInteger) result.state = String(getrand());
        } else if (operation === 3n)
          result = {
            value: [decode(coefficients[0]!), decode(coordinate)]
              .sort(compare_elements)
              .map(String),
          };
        else if (operation === 1n) {
          const H: any = HyperellipticCurve(f),
            changed = H.odd_degree_model();
          result = {
            value: [changed.hyperelliptic_polynomials()[0].coeffs.map(String), changed === H],
          };
        } else if (operation === 4n) {
          const H: any = HyperellipticCurve(f, R.one());
          result = {
            value: H.lift_x(decode(coordinate), { all: true }).map((P: any) => [
              P.coords.map(String),
              P.curve === H,
            ]),
          };
        } else {
          const x = R.gen(),
            t = decode(coordinate),
            a = x.pow(3).add(x).add(R.one());
          const b = R.__call__(t),
            h = R.one(),
            r = decode(coefficients[0]!);
          const q = x.pow(3).scalar_mul(r.mul(r)).add(x).add(R.one());
          const f = a.mul(q).add(h.mul(b)).add(b.mul(b));
          const [aa, bb] = cantor_reduction(a, b, f, h, 2);
          result = { value: [aa.coeffs.map(String), bb.coeffs.map(String)] };
        }
      } finally {
        proto.sqrt = sqrt;
        proto.is_square = square;
        poly.roots = roots;
      }
    } catch (e) {
      result = { error: (e as Error).name, message: (e as Error).message };
    }
    result.calls = trace;
    return JSON.stringify(result);
  },
});
