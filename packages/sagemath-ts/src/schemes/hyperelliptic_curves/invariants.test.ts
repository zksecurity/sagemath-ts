/**
 * Tests for the invariants of binary sextics and of genus 2 curves.
 *
 * Expected values are the doctest values of
 * `sage/schemes/hyperelliptic_curves/invariants.py` and
 * `hyperelliptic_g2.py`, re-run against the installed SageMath.
 */

import { describe, expect, it } from 'bun:test';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import type { FiniteFieldElement } from '../../rings/finite_rings/finite_field_prime.js';
import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import { Rational } from '../../rings/rational.js';
import { QQ } from '../../rings/rational_field.js';
import { HyperellipticCurve } from './constructor.js';
import {
  absolute_igusa_invariants_kohel,
  absolute_igusa_invariants_wamelen,
  clebsch_invariants,
  clebsch_to_igusa,
  igusa_clebsch_invariants,
  igusa_to_clebsch,
  sextic_form,
  ubs,
} from './invariants.js';

type FF = FiniteFieldElement;

// biome-ignore lint/suspicious/noExplicitAny: the constructor returns a dynamic subclass
type AnyCurve = any;

const R = new PolynomialRing<Rational>(QQ as never, 'x');
const x = R.gen();
const c = (n: bigint) => R.__call__(new Rational(n));
const Q = (n: bigint) => new Rational(n);
const strs = (a: readonly unknown[]) => a.map(String);

describe('clebsch_invariants (invariants.py:260-290)', () => {
  it('x^6 + 1', () => {
    expect(strs(clebsch_invariants(x.pow(6).add(c(1n))))).toEqual(['2', '2/3', '-2/9', '0']);
  });

  it('x^6 + x^5 + x^4 + x^2 + 2', () => {
    const f = x.pow(6).add(x.pow(5)).add(x.pow(4)).add(x.pow(2)).add(c(2n));
    expect(strs(clebsch_invariants(f))).toEqual([
      '62/15',
      '15434/5625',
      '-236951/140625',
      '229930748/791015625',
    ]);
  });

  it('rejects characteristics 2, 3 and 5', () => {
    for (const p of [2n, 3n, 5n]) {
      const K = GF(p);
      const S = new PolynomialRing<FF>(K, 'x');
      const t = S.gen();
      expect(() => clebsch_invariants(t.pow(5).sub(t))).toThrow(
        'Invariants of binary sextics/genus 2 hyperelliptic curves not implemented in characteristics 2, 3, and 5'
      );
    }
  });
});

describe('igusa_clebsch_invariants (invariants.py:293-331)', () => {
  it('x^6 + 1', () => {
    expect(strs(igusa_clebsch_invariants(x.pow(6).add(c(1n))))).toEqual([
      '-240',
      '1620',
      '-119880',
      '-46656',
    ]);
  });

  it('x^6 + x^5 + x^4 + x^2 + 2', () => {
    const f = x.pow(6).add(x.pow(5)).add(x.pow(4)).add(x.pow(2)).add(c(2n));
    expect(strs(igusa_clebsch_invariants(f))).toEqual(['-496', '6220', '-955932', '-1111784']);
  });
});

describe('clebsch_to_igusa / igusa_to_clebsch (invariants.py:208-257)', () => {
  it('round trips over QQ', () => {
    const I = clebsch_to_igusa(QQ as never, Q(2n), Q(3n), Q(4n), Q(5n));
    expect(strs(I)).toEqual(['-240', '17370', '231120', '-103098906']);
    expect(strs(igusa_to_clebsch(QQ as never, ...I))).toEqual(['2', '3', '4', '5']);
    expect(
      strs(igusa_to_clebsch(QQ as never, Q(-2400n), Q(173700n), Q(23112000n), Q(-10309890600n)))
    ).toEqual(['20', '342/5', '2512/5', '43381012/1125']);
  });

  it('round trips over GF(31)', () => {
    const F = GF(31n);
    const E = (n: bigint) => F.__call__(n);
    const I = clebsch_to_igusa(F, E(2n), E(3n), E(4n), E(5n));
    expect(strs(I)).toEqual(['8', '10', '15', '26']);
    expect(strs(igusa_to_clebsch(F, ...I))).toEqual(['2', '3', '4', '5']);

    const Is: [FF, FF, FF, FF] = [E(-2400n), E(173700n), E(23112000n), E(-10309890600n)];
    expect(strs(Is)).toEqual(['18', '7', '12', '27']);
    const Cs = igusa_to_clebsch(F, ...Is);
    expect(strs(Cs)).toEqual(['20', '25', '25', '12']);
    expect(strs(clebsch_to_igusa(F, ...Cs))).toEqual(['18', '7', '12', '27']);
  });
});

describe('absolute Igusa invariants (invariants.py:334-409)', () => {
  const h = x
    .pow(5)
    .neg()
    .add(x.pow(4).scalar_mul(Q(3n)))
    .add(x.pow(3).scalar_mul(Q(2n)))
    .sub(x.pow(2).scalar_mul(Q(6n)))
    .sub(x.scalar_mul(Q(3n)))
    .add(c(1n));

  it('van Wamelen', () => {
    expect(strs(absolute_igusa_invariants_wamelen(x.pow(5).sub(c(1n))))).toEqual(['0', '0', '0']);
    // 2^7*3^15, 2^5*3^11*5, 2^4*3^9*31
    expect(strs(absolute_igusa_invariants_wamelen(h))).toEqual([
      '1836660096',
      '28343520',
      '9762768',
    ]);
  });

  it('Kohel', () => {
    expect(strs(absolute_igusa_invariants_kohel(x.pow(5).sub(c(1n))))).toEqual(['0', '0', '0']);
    expect(strs(absolute_igusa_invariants_kohel(x.pow(5).sub(x)))).toEqual([
      '100',
      '-20000',
      '-2000',
    ]);
    expect(strs(absolute_igusa_invariants_kohel(h))).toEqual(['150660', '28343520', '9762768']);
  });
});

describe('ubs over GF(31) (invariants.py:176-187)', () => {
  it('reproduces the docstring dictionary', () => {
    const F = GF(31n);
    const S = new PolynomialRing<FF>(F, 't');
    const t = S.gen();
    const E = (n: bigint) => F.__call__(n);
    const p = t
      .pow(6)
      .add(t.pow(5).scalar_mul(E(2n)))
      .add(t.pow(2))
      .add(t.scalar_mul(E(3n)))
      .add(S.one());
    const U = ubs(F, sextic_form(p));
    expect(String(U.A)).toBe('0');
    expect(String(U.B)).toBe(String(E(-12n)));
    expect(String(U.C)).toBe(String(E(-15n)));
    expect(String(U.D)).toBe(String(E(-15n)));
    // Delta = -10*t^4 + 12*t^3*h + 7*t^2*h^2 - 5*t*h^3 + 2*h^4
    expect(strs(U.Delta.c)).toEqual(strs([E(2n), E(-5n), E(7n), E(12n), E(-10n)]));
    // i = -4*t^4 + 10*t^3*h + 2*t^2*h^2 - 9*t*h^3 - 7*h^4
    expect(strs(U.i.c)).toEqual(strs([E(-7n), E(-9n), E(2n), E(10n), E(-4n)]));
    // y1 = 4*t^2 - 10*t*h - 13*h^2
    expect(strs(U.y1.c)).toEqual(strs([E(-13n), E(-10n), E(4n)]));
    // y2 = 6*t^2 - 4*t*h + 2*h^2
    expect(strs(U.y2.c)).toEqual(strs([E(2n), E(-4n), E(6n)]));
    // y3 = 4*t^2 - 4*t*h - 9*h^2
    expect(strs(U.y3.c)).toEqual(strs([E(-9n), E(-4n), E(4n)]));
  });
});

describe('genus 2 curve invariants (hyperelliptic_g2.py:83-191)', () => {
  const f = x.pow(5).sub(x.pow(4)).add(c(3n));
  /** substitute `x -> a*x` into `p` */
  const scale = (p: typeof f, a: bigint) => {
    let out = R.zero();
    for (let i = 0; i <= p.degree(); i++) {
      out = out.add(
        R.__call__(new Rational(a ** BigInt(i)))
          .mul(x.pow(i))
          .scalar_mul(p.getCoeff(i))
      );
    }
    return out;
  };

  it('clebsch_invariants', () => {
    expect(strs((HyperellipticCurve(f) as AnyCurve).clebsch_invariants())).toEqual([
      '0',
      '-2048/375',
      '-4096/25',
      '-4881645568/84375',
    ]);
    expect(strs((HyperellipticCurve(scale(f, 2n)) as AnyCurve).clebsch_invariants())).toEqual([
      '0',
      '-8388608/375',
      '-1073741824/25',
      '-5241627016305836032/84375',
    ]);
    expect(strs((HyperellipticCurve(f, x) as AnyCurve).clebsch_invariants())).toEqual([
      '-8/15',
      '17504/5625',
      '-23162896/140625',
      '-420832861216768/7119140625',
    ]);
    expect(
      strs((HyperellipticCurve(scale(f, 2n), x.scalar_mul(Q(2n))) as AnyCurve).clebsch_invariants())
    ).toEqual([
      '-512/15',
      '71696384/5625',
      '-6072014209024/140625',
      '-451865844002031331704832/7119140625',
    ]);
  });

  it('igusa_clebsch_invariants', () => {
    const f2 = x.pow(5).sub(x).add(c(2n));
    expect(strs((HyperellipticCurve(f2) as AnyCurve).igusa_clebsch_invariants())).toEqual([
      '-640',
      '-20480',
      '1310720',
      '52160364544',
    ]);
    expect(
      strs((HyperellipticCurve(scale(f2, 2n)) as AnyCurve).igusa_clebsch_invariants())
    ).toEqual(['-40960', '-83886080', '343597383680', '56006764965979488256']);
    expect(strs((HyperellipticCurve(f2, x) as AnyCurve).igusa_clebsch_invariants())).toEqual([
      '-640',
      '17920',
      '-1966656',
      '52409511936',
    ]);
    expect(
      strs(
        (
          HyperellipticCurve(scale(f2, 2n), x.scalar_mul(Q(2n))) as AnyCurve
        ).igusa_clebsch_invariants()
      )
    ).toEqual(['-40960', '73400320', '-515547070464', '56274284941110411264']);
  });

  it('absolute invariants are invariant under x -> 3x + 1', () => {
    const sub3x1 = (p: typeof f) => {
      const lin = c(1n).add(x.scalar_mul(Q(3n)));
      let out = R.zero();
      for (let i = 0; i <= p.degree(); i++) {
        out = out.add(lin.pow(i).scalar_mul(p.getCoeff(i)));
      }
      return out;
    };
    const base = x.pow(5).sub(x).add(c(1n));
    const expected = ['-1030567/178769', '259686400/178769', '20806400/178769'];
    expect(
      strs((HyperellipticCurve(base, x.pow(2)) as AnyCurve).absolute_igusa_invariants_kohel())
    ).toEqual(expected);
    expect(
      strs(
        (
          HyperellipticCurve(sub3x1(base), sub3x1(x.pow(2))) as AnyCurve
        ).absolute_igusa_invariants_kohel()
      )
    ).toEqual(expected);
    expect(
      strs(
        (HyperellipticCurve(x.pow(5).sub(c(1n))) as AnyCurve).absolute_igusa_invariants_wamelen()
      )
    ).toEqual(['0', '0', '0']);
  });

  it('is_odd_degree (hyperelliptic_g2.py:16-36)', () => {
    expect((HyperellipticCurve(f) as AnyCurve).is_odd_degree()).toBe(true);
    const sextic = x
      .pow(2)
      .add(c(2n))
      .mul(x.pow(2).add(c(3n)))
      .mul(x.pow(2).add(c(5n)));
    expect((HyperellipticCurve(sextic) as AnyCurve).is_odd_degree()).toBe(false);
  });

  it('kummer_morphism is not ported', () => {
    expect(() => (HyperellipticCurve(f) as AnyCurve).kummer_morphism()).toThrow(
      'SAGE_NOT_IMPLEMENTED: kummer_morphism'
    );
    expect(() => (HyperellipticCurve(f) as AnyCurve).jacobian().kummer_surface()).toThrow(
      'SAGE_NOT_IMPLEMENTED: kummer_surface'
    );
  });
});
