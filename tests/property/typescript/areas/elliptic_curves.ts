/**
 * sagemath-ts side of the `elliptic_curves` property-test area.
 *
 * Cases: tests/property/cases/elliptic_curves.cases.json
 * SageMath counterpart: tests/property/python/areas/elliptic_curves.py
 */

import {
  EllipticCurve,
  GF,
  embedding_degree as ecEmbeddingDegree,
  is_ordinary as ecIsOrdinary,
  is_supersingular as ecIsSupersingular,
} from '../../../../packages/sagemath-ts/src/index.js';
import type { FiniteFieldElement } from '../../../../packages/sagemath-ts/src/rings/finite_rings/index.js';
import type { EllipticCurveFiniteField } from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/index.js';
import {
  EllipticCurveGeneric,
  EllipticCurveIsogeny,
} from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/index.js';

export const functions = {
  point_add: (p: bigint, a: bigint, b: bigint, x1: bigint, y1: bigint, x2: bigint, y2: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    const P1 = E.point(x1, y1);
    const P2 = E.point(x2, y2);
    const R = P1.add(P2);
    return R.isZero() ? '(0 : 1 : 0)' : `(${R.x!.value} : ${R.y!.value} : 1)`;
  },
  point_add_identity: (p: bigint, a: bigint, b: bigint, x: bigint, y: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    const P = E.point(x, y);
    const R = P.add(E.zero());
    return R.isZero() ? '(0 : 1 : 0)' : `(${R.x!.value} : ${R.y!.value} : 1)`;
  },
  scalar_mul: (p: bigint, a: bigint, b: bigint, x: bigint, y: bigint, n: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    const P = E.point(x, y);
    const R = P.mul(n);
    return R.isZero() ? '(0 : 1 : 0)' : `(${R.x!.value} : ${R.y!.value} : 1)`;
  },
  point_neg: (p: bigint, a: bigint, b: bigint, x: bigint, y: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    const R = E.point(x, y).neg();
    return R.isZero() ? '(0 : 1 : 0)' : `(${R.x!.value} : ${R.y!.value} : 1)`;
  },
  point_order: (p: bigint, a: bigint, b: bigint, x: bigint, y: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    return E.point(x, y).order().toString();
  },
  ellcard: (p: bigint, a: bigint, b: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    return E.cardinality().toString();
  },
  discriminant: (p: bigint, a: bigint, b: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    return E.discriminant().value.toString();
  },
  j_invariant: (p: bigint, a: bigint, b: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    return E.j_invariant().value.toString();
  },
  trace_of_frobenius: (p: bigint, a: bigint, b: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    return E.trace_of_frobenius().toString();
  },
  is_supersingular: (p: bigint, a: bigint, b: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    return ecIsSupersingular(E);
  },
  is_ordinary: (p: bigint, a: bigint, b: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    return ecIsOrdinary(E);
  },
  embedding_degree: (p: bigint, a: bigint, b: bigint, n: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    return ecEmbeddingDegree(E, n).toString();
  },
  division_polynomial: (p: bigint, a: bigint, b: bigint, n: bigint) => {
    const F = GF(p);
    // Use EllipticCurveGeneric which has the division_polynomial method
    const E = new EllipticCurveGeneric(F, [
      F.__call__(0n),
      F.__call__(0n),
      F.__call__(0n),
      F.__call__(a),
      F.__call__(b),
    ]);
    const poly = E.division_polynomial(Number(n));
    // Format polynomial coefficients using .coeffs property
    const coeffs = poly.coeffs.map((c: FiniteFieldElement) => c.value.toString());
    return coeffs.join(',');
  },
  isogeny_degree: (p: bigint, a: bigint, b: bigint, deg: bigint) => {
    try {
      const F = GF(p);
      const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
      const order = E.cardinality();
      if (order % deg !== 0n) return 'no_isogeny';
      const cofactor = order / deg;
      const Egen = new EllipticCurveGeneric(F, [
        F.__call__(0n),
        F.__call__(0n),
        F.__call__(0n),
        F.__call__(a),
        F.__call__(b),
      ]);
      for (let i = 0; i < 100; i++) {
        const P = E.random_point();
        const Q = P.mul(cofactor);
        if (!Q.isZero() && Q.mul(deg).isZero()) {
          const kernel = Egen.point([F.__call__(Q.x!.value), F.__call__(Q.y!.value)], false);
          return new EllipticCurveIsogeny(Egen, kernel).degree().toString();
        }
      }
      return 'no_point_found';
    } catch (e) {
      return `error: ${e instanceof Error ? e.message : String(e)}`;
    }
  },
  isogeny_codomain_j: (p: bigint, a: bigint, b: bigint, deg: bigint) => {
    // Mirrors ec_isogeny_codomain_j in tests/property/python/areas/elliptic_curves.py.
    // EllipticCurveIsogeny takes an EllipticCurveGeneric domain, while
    // cardinality()/random_point() live on EllipticCurveFiniteField, so the
    // kernel point is transferred between the two models of the same curve.
    try {
      const F = GF(p);
      const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
      const order = E.cardinality();
      if (order % deg !== 0n) return 'no_isogeny';
      const cofactor = order / deg;

      const Egen = new EllipticCurveGeneric(F, [
        F.__call__(0n),
        F.__call__(0n),
        F.__call__(0n),
        F.__call__(a),
        F.__call__(b),
      ]);

      for (let i = 0; i < 100; i++) {
        const P = E.random_point();
        const Q = P.mul(cofactor);
        if (!Q.isZero() && Q.mul(deg).isZero()) {
          const kernel = Egen.point([F.__call__(Q.x!.value), F.__call__(Q.y!.value)], false);
          const phi = new EllipticCurveIsogeny(Egen, kernel);
          return phi.codomain().j_invariant().toString();
        }
      }
      return 'no_point_found';
    } catch (e) {
      return `error: ${e instanceof Error ? e.message : String(e)}`;
    }
  },
  torsion_points_count: (p: bigint, a: bigint, b: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    return E.cardinality().toString();
  },
  is_on_curve: (p: bigint, a: bigint, b: bigint, x: bigint, y: bigint) => {
    const F = GF(p);
    const xV = F.__call__(x);
    const yV = F.__call__(y);
    const aV = F.__call__(a);
    const bV = F.__call__(b);
    return yV.mul(yV).eq(xV.mul(xV).mul(xV).add(aV.mul(xV)).add(bV));
  },
  lift_x: (p: bigint, a: bigint, b: bigint, x: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    try {
      const P = E.lift_x(x);
      return P.isZero() ? '(0 : 1 : 0)' : `(${P.x!.value} : ${P.y!.value} : 1)`;
    } catch {
      return 'no_point';
    }
  },
  associativity_check: (
    p: bigint,
    a: bigint,
    b: bigint,
    x1: bigint,
    y1: bigint,
    x2: bigint,
    y2: bigint,
    x3: bigint,
    y3: bigint
  ) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    const P = E.point(x1, y1);
    const Q = E.point(x2, y2);
    const R = E.point(x3, y3);
    return P.add(Q)
      .add(R)
      .eq(P.add(Q.add(R)));
  },
  inverse_check: (p: bigint, a: bigint, b: bigint, x: bigint, y: bigint) => {
    const F = GF(p);
    const E = EllipticCurve(F, [a, b]) as EllipticCurveFiniteField;
    return E.point(x, y).add(E.point(x, y).neg()).isZero();
  },
};
