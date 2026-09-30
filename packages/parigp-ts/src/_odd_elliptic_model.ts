/** PARI ff.c:1145–1197/1363/1466 and elliptic.c:463/798.
 * Copyright (C) The PARI group; adapted under GPL-2.0-or-later.
 * @see Deviation: PARI odd-extension elliptic model adapters
 */
import { extensionField } from './_extension_field.js';
import { EllCurveType } from './elliptic/init.js';
import { type PariFfelt } from './types.js';
import type {
  FFEllipticScalar,
  FFEllipticInvariants,
  FFEllipticInputPoint,
  FFEllipticPoint,
  OddFFEllipticCurve,
} from './_elliptic_finite_field.js';
import { FlxqE_changepoint, FlxqE_changepointinv, FlxqE_mul } from './FlxqE.js';
import { FpXQE_changepoint, FpXQE_changepointinv, FpXQE_mul } from './FpE.js';

function arithmetic(fg: PariFfelt) {
  const p = fg.p,
    T = [...(fg.definingPoly ?? [0n, 1n])],
    word = p < 1n << 64n;
  const F = extensionField(word ? 1 : 0, T, p);
  const coefficient = (x: FFEllipticScalar): bigint[] => {
    const value = typeof x === 'bigint' ? x : x.value;
    // Integer input is a constant; field input uses ascending coefficients.
    const a = typeof value === 'bigint' ? [((value % p) + p) % p] : [...value];
    while (a.length && a[a.length - 1] === 0n) a.pop();
    return F.red(a) as bigint[];
  };
  const element = (value: bigint[]): PariFfelt => ({ ...fg, value });
  const add = (a: bigint[], b: bigint[]) => F.add(a, b) as bigint[];
  const sub = (a: bigint[], b: bigint[]) => F.sub(a, b) as bigint[];
  const mul = (a: bigint[], b: bigint[]) => F.mul(a, b) as bigint[];
  const sqr = (a: bigint[]) => F.sqr(a) as bigint[];
  const scale = (a: bigint[], n: bigint) => F.mul(a, ((n % p) + p) % p) as bigint[];
  const div = (a: bigint[], b: bigint[]) => F.mul(a, F.inv(b)) as bigint[];
  return { p, T, word, coefficient, element, add, sub, mul, sqr, scale, div };
}
export function oddFFInit(E: FFEllipticInvariants, fg: PariFfelt): OddFFEllipticCurve {
  const F = arithmetic(fg),
    { add, sub, mul, sqr, scale, div } = F;
  const a1 = F.coefficient(E.a1),
    a3 = F.coefficient(E.a3),
    b2 = F.coefficient(E.b2);
  let oddModel: OddFFEllipticCurve['oddModel'];
  if (fg.p === 3n) {
    const b4 = F.coefficient(E.b4),
      b6 = F.coefficient(E.b6);
    if (b2.length) {
      const r = div(b4, b2);
      oddModel = [
        [b2],
        sub(b6, mul(r, add(b4, sqr(r)))),
        [[1n], r, scale(a1, -1n), scale(a3, -1n)],
      ];
    } else oddModel = [scale(b4, -1n), b6, [[1n], [], scale(a1, -1n), scale(a3, -1n)]];
  } else {
    oddModel = [
      scale(F.coefficient(E.c4), -27n),
      scale(F.coefficient(E.c6), -54n),
      [[6n % fg.p], scale(b2, 3n), scale(a1, 3n), scale(a3, 108n)],
    ];
  }
  const keys = ['a1', 'a2', 'a3', 'a4', 'a6', 'b2', 'b4', 'b6', 'b8', 'c4', 'c6', 'disc'] as const;
  const values = Object.fromEntries(keys.map((k) => [k, F.element(F.coefficient(E[k]))])) as {
    [K in keyof FFEllipticInvariants]: PariFfelt;
  };
  const D = F.coefficient(E.disc),
    c4 = F.coefficient(E.c4);
  return {
    ...values,
    type: EllCurveType.t_ELL_Fq,
    field: fg,
    oddModel,
    j: F.element(D.length ? div(mul(sqr(c4), c4), D) : []),
  };
}
export function oddFFMul(
  E: OddFFEllipticCurve,
  P: FFEllipticInputPoint,
  n: bigint
): FFEllipticPoint {
  const F = arithmetic(E.field),
    [a, , ch] = E.oddModel;
  const point = P.isInfinity
    ? P
    : { isInfinity: false as const, x: F.coefficient(P.x), y: F.coefficient(P.y) };
  const Q = F.word
    ? FlxqE_changepoint(
        FlxqE_mul(FlxqE_changepointinv(point, ch, F.T, F.p), n, a, F.T, F.p),
        ch,
        F.T,
        F.p
      )
    : FpXQE_changepoint(
        FpXQE_mul(FpXQE_changepointinv(point, ch, F.T, F.p), n, a as bigint[], F.T, F.p),
        ch,
        F.T,
        F.p
      );
  return Q.isInfinity ? Q : { isInfinity: false, x: F.element(Q.x), y: F.element(Q.y) };
}
export function oddInitFq(
  x: readonly FFEllipticScalar[],
  fg: PariFfelt
): OddFFEllipticCurve | null {
  const F = arithmetic(fg),
    { add, sub, mul, sqr, scale } = F;
  if (x.length === 1) {
    // ellinit_Fq coerces j to the supplied field before ellfromj branches.
    const j = F.coefficient(x[0]!);
    if (fg.p === 3n)
      x = j.length ? [0n, F.element(j), 0n, 0n, F.element(scale(sqr(j), -1n))] : [0n, 0n, 0n, 1n, 0n];
    else if (!j.length) x = [0n, 0n, 0n, 0n, 1n];
    else if (!sub(j, F.coefficient(1728n)).length) x = [0n, 0n, 0n, 1n, 0n];
    else {
      const k = sub(F.coefficient(1728n), j),
        kj = mul(k, j);
      x = [0n, 0n, 0n, F.element(scale(kj, 3n)), F.element(scale(mul(kj, k), 2n))];
    }
  }
  const [a1, a2, a3, a4, a6] = x.map(F.coefficient) as [
    bigint[],
    bigint[],
    bigint[],
    bigint[],
    bigint[],
  ];
  const a11 = sqr(a1),
    a13 = mul(a1, a3);
  const b2 = add(a11, scale(a2, 4n)),
    b4 = add(a13, scale(a4, 2n)),
    b6 = add(sqr(a3), scale(a6, 4n));
  const b8 = sub(add(mul(a11, a6), mul(b6, a2)), mul(a4, add(a4, a13))),
    b22 = sqr(b2);
  const c4 = sub(b22, scale(b4, 24n)),
    c6 = sub(mul(b2, sub(scale(b4, 36n), b22)), scale(b6, 216n));
  const disc = sub(
    mul(b4, sub(scale(mul(b2, b6), 9n), scale(sqr(b4), 8n))),
    add(mul(b22, b8), scale(sqr(b6), 27n))
  );
  if (!disc.length) return null;
  const values = { a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, disc };
  const E = Object.fromEntries(
    Object.entries(values).map(([k, v]) => [k, F.element(v)])
  ) as FFEllipticInvariants;
  return oddFFInit(E, fg);
}
