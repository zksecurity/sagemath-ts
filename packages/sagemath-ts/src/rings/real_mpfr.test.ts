/**
 * Unit tests for sage/rings/real_mpfr
 * Tests for real number transcendental functions
 */
import { describe, expect, test } from 'bun:test';
import {
  RR,
  RealField,
  RealNumber,
  RoundingMode,
  mpfr_prec_max,
  mpfr_prec_min,
} from './real_mpfr.js';

const EPSILON = 1e-10; // Tolerance for floating point comparisons

function approxEqual(a: number, b: number, eps: number = EPSILON): boolean {
  if (Number.isNaN(a) && Number.isNaN(b)) return true;
  if (!Number.isFinite(a) || !Number.isFinite(b)) return a === b;
  return Math.abs(a - b) < eps;
}

describe('RealField', () => {
  test('constructor defaults', () => {
    const R = new RealField();
    expect(R.precision()).toBe(53);
    expect(R.characteristic()).toBe(0n);
    expect(R.rounding_mode()).toBe(RoundingMode.RNDN);
  });

  test('constructor with custom precision', () => {
    const R = new RealField(100);
    expect(R.precision()).toBe(100);
  });

  test('invalid precision throws', () => {
    expect(() => new RealField(0)).toThrow();
    expect(() => new RealField(-1)).toThrow();
  });

  test('pi', () => {
    const R = new RealField();
    const pi = R.pi();
    expect(approxEqual(pi.toNumber(), Math.PI)).toBe(true);
  });

  test('euler_constant', () => {
    const R = new RealField();
    const euler = R.euler_constant();
    expect(approxEqual(euler.toNumber(), 0.5772156649015329)).toBe(true);
  });

  test('catalan_constant', () => {
    const R = new RealField();
    const catalan = R.catalan_constant();
    expect(approxEqual(catalan.toNumber(), 0.915965594177219)).toBe(true);
  });

  test('log2', () => {
    const R = new RealField();
    const log2 = R.log2();
    expect(approxEqual(log2.toNumber(), Math.LN2)).toBe(true);
  });

  test('random_element', () => {
    const R = new RealField();
    for (let i = 0; i < 10; i++) {
      const r = R.random_element(-5, 5);
      expect(r.toNumber()).toBeGreaterThanOrEqual(-5);
      expect(r.toNumber()).toBeLessThanOrEqual(5);
    }
  });

  test('factorial', () => {
    const R = new RealField();
    expect(R.factorial(0).toNumber()).toBe(1);
    expect(R.factorial(1).toNumber()).toBe(1);
    expect(R.factorial(5).toNumber()).toBe(120);
    expect(R.factorial(10).toNumber()).toBe(3628800);
    expect(() => R.factorial(-1)).toThrow();
  });

  test('zeta', () => {
    const R = new RealField();
    expect(R.zeta(1).toNumber()).toBe(1);
    expect(R.zeta(2).toNumber()).toBe(-1);
    expect(() => R.zeta(3)).toThrow();
    expect(() => R.zeta(5)).toThrow();
  });

  test('toString', () => {
    const R = new RealField(100);
    expect(R.toString()).toBe('Real Field with 100 bits of precision');

    const Rz = new RealField(53, false, RoundingMode.RNDZ);
    expect(Rz.toString()).toBe('Real Field with 53 bits of precision and rounding RNDZ');
  });
});

describe('RealNumber - Basic operations', () => {
  const R = new RealField();

  test('construction', () => {
    const n1 = R.__call__(3.14);
    expect(n1.toNumber()).toBe(3.14);

    const n2 = R.__call__(42n);
    expect(n2.toNumber()).toBe(42);

    const n3 = R.__call__('2.718');
    expect(approxEqual(n3.toNumber(), 2.718)).toBe(true);
  });

  test('precision', () => {
    const n = R.__call__(1.5);
    expect(n.precision()).toBe(53);
  });

  test('real and imag', () => {
    const n = R.__call__(3.14);
    expect(n.real().toNumber()).toBe(3.14);
    expect(n.imag().toNumber()).toBe(0);
  });

  test('sign', () => {
    expect(R.__call__(5).sign()).toBe(1);
    expect(R.__call__(-5).sign()).toBe(-1);
    expect(R.__call__(0).sign()).toBe(0);
  });

  test('abs', () => {
    expect(R.__call__(-5).abs().toNumber()).toBe(5);
    expect(R.__call__(5).abs().toNumber()).toBe(5);
    expect(R.__call__(0).abs().toNumber()).toBe(0);
  });

  test('floor', () => {
    expect(R.__call__(3.7).floor()).toBe(3n);
    expect(R.__call__(-3.7).floor()).toBe(-4n);
    expect(R.__call__(3).floor()).toBe(3n);
  });

  test('ceil', () => {
    expect(R.__call__(3.2).ceil()).toBe(4n);
    expect(R.__call__(-3.2).ceil()).toBe(-3n);
    expect(R.__call__(3).ceil()).toBe(3n);
  });

  test('round', () => {
    // mpfr_round rounds halfway cases AWAY FROM ZERO, unlike Math.round which
    // rounds them towards +infinity.
    // sage: RR(0.49).round() -> 0 ; RR(0.5).round() -> 1
    // sage: RR(-0.49).round() -> 0 ; RR(-0.5).round() -> -1
    expect(R.__call__(3.4).round()).toBe(3n);
    expect(R.__call__(3.5).round()).toBe(4n);
    expect(R.__call__(-3.5).round()).toBe(-4n);
    expect(R.__call__(-0.5).round()).toBe(-1n);
    expect(R.__call__(0.5).round()).toBe(1n);
    expect(R.__call__(2.5).round()).toBe(3n);
    expect(R.__call__(0.49).round()).toBe(0n);
    expect(R.__call__(-0.49).round()).toBe(0n);
    expect(R.__call__(-1.5).round()).toBe(-2n);
  });

  test('trunc', () => {
    expect(R.__call__(3.7).trunc()).toBe(3n);
    expect(R.__call__(-3.7).trunc()).toBe(-3n);
  });

  test('frac', () => {
    expect(approxEqual(R.__call__(3.7).frac().toNumber(), 0.7)).toBe(true);
    expect(approxEqual(R.__call__(-3.7).frac().toNumber(), -0.7)).toBe(true);
  });

  test('arithmetic operations', () => {
    const a = R.__call__(5);
    const b = R.__call__(3);

    expect(a.add(b).toNumber()).toBe(8);
    expect(a.sub(b).toNumber()).toBe(2);
    expect(a.mul(b).toNumber()).toBe(15);
    expect(approxEqual(a.div(b).toNumber(), 5 / 3)).toBe(true);
    expect(a.neg().toNumber()).toBe(-5);
  });

  test('comparison', () => {
    const a = R.__call__(5);
    const b = R.__call__(3);

    expect(a.cmp(b)).toBe(1);
    expect(b.cmp(a)).toBe(-1);
    expect(a.cmp(5)).toBe(0);
    expect(a.equals(5)).toBe(true);
    expect(a.equals(b)).toBe(false);
  });
});

describe('RealNumber - Power functions', () => {
  const R = new RealField();

  test('sqrt', () => {
    expect(R.__call__(4).sqrt().toNumber()).toBe(2);
    expect(R.__call__(9).sqrt().toNumber()).toBe(3);
    expect(approxEqual(R.__call__(2).sqrt().toNumber(), Math.SQRT2)).toBe(true);
    expect(R.__call__(0).sqrt().toNumber()).toBe(0);
    // Negative numbers return NaN
    expect(Number.isNaN(R.__call__(-1).sqrt().toNumber())).toBe(true);
  });

  test('cube_root', () => {
    expect(R.__call__(8).cube_root().toNumber()).toBe(2);
    expect(R.__call__(-8).cube_root().toNumber()).toBe(-2);
    expect(approxEqual(R.__call__(27).cube_root().toNumber(), 3)).toBe(true);
  });

  test('nth_root', () => {
    expect(R.__call__(16).nth_root(4).toNumber()).toBe(2);
    expect(R.__call__(32).nth_root(5).toNumber()).toBe(2);
    expect(R.__call__(-27).nth_root(3).toNumber()).toBe(-3);
    expect(R.__call__(0).nth_root(5).toNumber()).toBe(0);
    expect(() => R.__call__(16).nth_root(-1)).toThrow();
    expect(() => R.__call__(-16).nth_root(4)).toThrow(); // Even root of negative
  });

  test('pow', () => {
    expect(R.__call__(2).pow(3).toNumber()).toBe(8);
    expect(R.__call__(2).pow(0.5).toNumber()).toBe(Math.SQRT2);
    expect(R.__call__(10).pow(-1).toNumber()).toBe(0.1);
    expect(R.__call__(0).pow(0).toNumber()).toBe(1); // 0^0 = 1 by convention
  });
});

describe('RealNumber - Exponential and logarithmic functions', () => {
  const R = new RealField();

  test('exp', () => {
    expect(R.__call__(0).exp().toNumber()).toBe(1);
    expect(approxEqual(R.__call__(1).exp().toNumber(), Math.E)).toBe(true);
    expect(approxEqual(R.__call__(2).exp().toNumber(), Math.E * Math.E)).toBe(true);
  });

  test('exp2', () => {
    expect(R.__call__(0).exp2().toNumber()).toBe(1);
    expect(R.__call__(3).exp2().toNumber()).toBe(8);
    expect(R.__call__(10).exp2().toNumber()).toBe(1024);
  });

  test('exp10', () => {
    expect(R.__call__(0).exp10().toNumber()).toBe(1);
    expect(R.__call__(2).exp10().toNumber()).toBe(100);
    expect(R.__call__(3).exp10().toNumber()).toBe(1000);
  });

  test('expm1', () => {
    expect(R.__call__(0).expm1().toNumber()).toBe(0);
    expect(approxEqual(R.__call__(1).expm1().toNumber(), Math.E - 1)).toBe(true);
    // For small x, expm1 is more accurate
    const small = R.__call__(1e-15);
    expect(approxEqual(small.expm1().toNumber(), 1e-15, 1e-28)).toBe(true);
  });

  test('log', () => {
    expect(R.__call__(1).log().toNumber()).toBe(0);
    expect(approxEqual(R.__call__(Math.E).log().toNumber(), 1)).toBe(true);
    expect(R.__call__(0).log().toNumber()).toBe(Number.NEGATIVE_INFINITY);
    // Negative numbers return NaN
    expect(Number.isNaN(R.__call__(-1).log().toNumber())).toBe(true);
  });

  test('log with base', () => {
    expect(R.__call__(100).log(10).toNumber()).toBe(2);
    expect(R.__call__(8).log(2).toNumber()).toBe(3);
    expect(approxEqual(R.__call__(27).log(3).toNumber(), 3)).toBe(true);
  });

  test('log2', () => {
    expect(R.__call__(1).log2().toNumber()).toBe(0);
    expect(R.__call__(2).log2().toNumber()).toBe(1);
    expect(R.__call__(8).log2().toNumber()).toBe(3);
    expect(R.__call__(1024).log2().toNumber()).toBe(10);
  });

  test('log10', () => {
    expect(R.__call__(1).log10().toNumber()).toBe(0);
    expect(R.__call__(10).log10().toNumber()).toBe(1);
    expect(R.__call__(100).log10().toNumber()).toBe(2);
  });

  test('log1p', () => {
    expect(R.__call__(0).log1p().toNumber()).toBe(0);
    expect(
      approxEqual(
        R.__call__(Math.E - 1)
          .log1p()
          .toNumber(),
        1
      )
    ).toBe(true);
    expect(R.__call__(-1).log1p().toNumber()).toBe(Number.NEGATIVE_INFINITY);
    // For small x, log1p is more accurate
    const small = R.__call__(1e-15);
    expect(approxEqual(small.log1p().toNumber(), 1e-15, 1e-28)).toBe(true);
  });
});

describe('RealNumber - Trigonometric functions', () => {
  const R = new RealField();
  const pi = Math.PI;

  test('sin', () => {
    expect(R.__call__(0).sin().toNumber()).toBe(0);
    expect(
      approxEqual(
        R.__call__(pi / 2)
          .sin()
          .toNumber(),
        1
      )
    ).toBe(true);
    expect(approxEqual(R.__call__(pi).sin().toNumber(), 0)).toBe(true);
    expect(
      approxEqual(
        R.__call__(pi / 6)
          .sin()
          .toNumber(),
        0.5
      )
    ).toBe(true);
  });

  test('cos', () => {
    expect(R.__call__(0).cos().toNumber()).toBe(1);
    expect(
      approxEqual(
        R.__call__(pi / 2)
          .cos()
          .toNumber(),
        0
      )
    ).toBe(true);
    expect(approxEqual(R.__call__(pi).cos().toNumber(), -1)).toBe(true);
    expect(
      approxEqual(
        R.__call__(pi / 3)
          .cos()
          .toNumber(),
        0.5
      )
    ).toBe(true);
  });

  test('tan', () => {
    expect(R.__call__(0).tan().toNumber()).toBe(0);
    expect(
      approxEqual(
        R.__call__(pi / 4)
          .tan()
          .toNumber(),
        1
      )
    ).toBe(true);
    expect(
      approxEqual(
        R.__call__(pi / 6)
          .tan()
          .toNumber(),
        1 / Math.sqrt(3)
      )
    ).toBe(true);
  });

  test('sincos', () => {
    const [s, c] = R.__call__(pi / 4).sincos();
    expect(approxEqual(s.toNumber(), Math.SQRT2 / 2)).toBe(true);
    expect(approxEqual(c.toNumber(), Math.SQRT2 / 2)).toBe(true);
  });

  test('arcsin', () => {
    expect(R.__call__(0).arcsin().toNumber()).toBe(0);
    expect(approxEqual(R.__call__(1).arcsin().toNumber(), pi / 2)).toBe(true);
    expect(approxEqual(R.__call__(0.5).arcsin().toNumber(), pi / 6)).toBe(true);
  });

  test('arccos', () => {
    expect(approxEqual(R.__call__(1).arccos().toNumber(), 0)).toBe(true);
    expect(approxEqual(R.__call__(0).arccos().toNumber(), pi / 2)).toBe(true);
    expect(approxEqual(R.__call__(0.5).arccos().toNumber(), pi / 3)).toBe(true);
  });

  test('arctan', () => {
    expect(R.__call__(0).arctan().toNumber()).toBe(0);
    expect(approxEqual(R.__call__(1).arctan().toNumber(), pi / 4)).toBe(true);
  });

  test('cot', () => {
    expect(
      approxEqual(
        R.__call__(pi / 4)
          .cot()
          .toNumber(),
        1
      )
    ).toBe(true);
    expect(
      approxEqual(
        R.__call__(pi / 6)
          .cot()
          .toNumber(),
        Math.sqrt(3)
      )
    ).toBe(true);
  });

  test('sec', () => {
    expect(R.__call__(0).sec().toNumber()).toBe(1);
    expect(
      approxEqual(
        R.__call__(pi / 3)
          .sec()
          .toNumber(),
        2
      )
    ).toBe(true);
  });

  test('csc', () => {
    expect(
      approxEqual(
        R.__call__(pi / 2)
          .csc()
          .toNumber(),
        1
      )
    ).toBe(true);
    expect(
      approxEqual(
        R.__call__(pi / 6)
          .csc()
          .toNumber(),
        2
      )
    ).toBe(true);
  });
});

describe('RealNumber - Hyperbolic functions', () => {
  const R = new RealField();

  test('sinh', () => {
    expect(R.__call__(0).sinh().toNumber()).toBe(0);
    expect(approxEqual(R.__call__(1).sinh().toNumber(), (Math.E - 1 / Math.E) / 2)).toBe(true);
  });

  test('cosh', () => {
    expect(R.__call__(0).cosh().toNumber()).toBe(1);
    expect(approxEqual(R.__call__(1).cosh().toNumber(), (Math.E + 1 / Math.E) / 2)).toBe(true);
  });

  test('tanh', () => {
    expect(R.__call__(0).tanh().toNumber()).toBe(0);
    const t = R.__call__(1).tanh().toNumber();
    const expected = (Math.E - 1 / Math.E) / (Math.E + 1 / Math.E);
    expect(approxEqual(t, expected)).toBe(true);
  });

  test('arcsinh', () => {
    expect(R.__call__(0).arcsinh().toNumber()).toBe(0);
    const x = R.__call__(1).arcsinh().toNumber();
    expect(approxEqual(R.__call__(x).sinh().toNumber(), 1)).toBe(true);
  });

  test('arccosh', () => {
    expect(R.__call__(1).arccosh().toNumber()).toBe(0);
    const x = R.__call__(2).arccosh().toNumber();
    expect(approxEqual(R.__call__(x).cosh().toNumber(), 2)).toBe(true);
  });

  test('arctanh', () => {
    expect(R.__call__(0).arctanh().toNumber()).toBe(0);
    const x = R.__call__(0.5).arctanh().toNumber();
    expect(approxEqual(R.__call__(x).tanh().toNumber(), 0.5)).toBe(true);
  });

  test('coth', () => {
    const c = R.__call__(1).coth().toNumber();
    const expected = 1 / R.__call__(1).tanh().toNumber();
    expect(approxEqual(c, expected)).toBe(true);
  });

  test('sech', () => {
    expect(R.__call__(0).sech().toNumber()).toBe(1);
    const s = R.__call__(1).sech().toNumber();
    expect(approxEqual(s, 1 / R.__call__(1).cosh().toNumber())).toBe(true);
  });

  test('csch', () => {
    const c = R.__call__(1).csch().toNumber();
    expect(approxEqual(c, 1 / R.__call__(1).sinh().toNumber())).toBe(true);
  });
});

describe('RealNumber - Special functions', () => {
  const R = new RealField();

  test('agm', () => {
    const a = R.__call__(1);
    const b = R.__call__(2);
    const agm = a.agm(b);
    // AGM(1,2) ≈ 1.4567910310469068...
    expect(approxEqual(agm.toNumber(), 1.4567910310469068, 1e-10)).toBe(true);

    // AGM is symmetric
    expect(approxEqual(b.agm(a).toNumber(), agm.toNumber())).toBe(true);

    // AGM lies between geometric and arithmetic mean
    expect(agm.toNumber()).toBeGreaterThan(Math.sqrt(2));
    expect(agm.toNumber()).toBeLessThan(1.5);
  });

  test('erf', () => {
    expect(approxEqual(R.__call__(0).erf().toNumber(), 0, 1e-6)).toBe(true);
    expect(approxEqual(R.__call__(1).erf().toNumber(), 0.8427007929497149, 1e-6)).toBe(true);
    expect(approxEqual(R.__call__(2).erf().toNumber(), 0.9953222650189527, 1e-6)).toBe(true);
    // erf is odd: erf(-x) = -erf(x)
    expect(approxEqual(R.__call__(-1).erf().toNumber(), -0.8427007929497149, 1e-6)).toBe(true);
  });

  test('erfc', () => {
    expect(approxEqual(R.__call__(0).erfc().toNumber(), 1, 1e-15)).toBe(true);
    expect(approxEqual(R.__call__(1).erfc().toNumber(), 0.15729920705028513, 1e-15)).toBe(true);
    expect(approxEqual(R.__call__(-1).erfc().toNumber(), 1.8427007929497148, 1e-15)).toBe(true);
    expect(approxEqual(R.__call__(2).erfc().toNumber(), 0.004677734981047266, 1e-17)).toBe(true);
  });

  test('erfc keeps relative accuracy in the tail', () => {
    // erfc = 1 - erf underflows to exactly 0 for x >~ 6; mpfr_erfc does not.
    // sage: R(6).erfc()  -> 2.15197367124989e-17
    // sage: R(10).erfc() -> 2.08848758376254e-45
    const e6 = R.__call__(6).erfc().toNumber();
    expect(Math.abs(e6 / 2.1519736712498913e-17 - 1)).toBeLessThan(1e-12);
    const e10 = R.__call__(10).erfc().toNumber();
    expect(Math.abs(e10 / 2.0884875837625446e-45 - 1)).toBeLessThan(1e-12);
    const e25 = R.__call__(25).erfc().toNumber();
    expect(Math.abs(e25 / 8.300172571196522e-274 - 1)).toBeLessThan(1e-12);
    // erf(x) + erfc(x) == 1 for moderate x
    for (const x of [0.25, 0.5, 1, 1.5, 2, 3]) {
      const s = R.__call__(x).erf().toNumber() + R.__call__(x).erfc().toNumber();
      expect(Math.abs(s - 1)).toBeLessThan(1e-15);
    }
  });

  test('gamma', () => {
    // gamma(n) = (n-1)! for positive integers
    expect(approxEqual(R.__call__(1).gamma().toNumber(), 1)).toBe(true);
    expect(approxEqual(R.__call__(2).gamma().toNumber(), 1)).toBe(true);
    expect(approxEqual(R.__call__(3).gamma().toNumber(), 2)).toBe(true);
    expect(approxEqual(R.__call__(6).gamma().toNumber(), 120)).toBe(true);

    // gamma(0.5) = sqrt(pi)
    expect(approxEqual(R.__call__(0.5).gamma().toNumber(), Math.sqrt(Math.PI), 1e-10)).toBe(true);
  });

  test('log_gamma', () => {
    expect(approxEqual(R.__call__(1).log_gamma().toNumber(), 0, 1e-10)).toBe(true);
    expect(approxEqual(R.__call__(6).log_gamma().toNumber(), Math.log(120), 1e-10)).toBe(true);
  });

  test('zeta (Riemann)', () => {
    // zeta(2) = pi^2/6
    expect(approxEqual(R.__call__(2).zeta().toNumber(), Math.PI ** 2 / 6, 1e-15)).toBe(true);
    // zeta(4) = pi^4/90
    expect(approxEqual(R.__call__(4).zeta().toNumber(), Math.PI ** 4 / 90, 1e-15)).toBe(true);
    // zeta(3) = Apery's constant
    expect(approxEqual(R.__call__(3).zeta().toNumber(), 1.2020569031595942, 1e-15)).toBe(true);

    // zeta(1) = infinity
    expect(R.__call__(1).zeta().toNumber()).toBe(Number.POSITIVE_INFINITY);

    // zeta(-2) = 0
    expect(R.__call__(-2).zeta().toNumber()).toBe(0);
    // zeta(-1) = -1/12, zeta(-3) = 1/120
    expect(approxEqual(R.__call__(-1).zeta().toNumber(), -1 / 12, 1e-15)).toBe(true);
    expect(approxEqual(R.__call__(-3).zeta().toNumber(), 1 / 120, 1e-15)).toBe(true);
  });

  test('zeta on the critical strip [0, 1)', () => {
    // mpfr_zeta is defined on the whole real line.  Applying the functional
    // equation for s <= 1 maps (0,1) onto itself and never terminates.
    // sage: RR(0.5).zeta() -> -1.46035450880959
    // sage: RR(0).zeta()   -> -0.5
    expect(approxEqual(R.__call__(0.5).zeta().toNumber(), -1.4603545088095868, 1e-14)).toBe(true);
    expect(approxEqual(R.__call__(0).zeta().toNumber(), -0.5, 1e-15)).toBe(true);
    expect(approxEqual(R.__call__(0.2).zeta().toNumber(), -0.7339209248963406, 1e-14)).toBe(true);
    expect(approxEqual(R.__call__(0.9).zeta().toNumber(), -9.430114019402254, 1e-13)).toBe(true);
    expect(approxEqual(R.__call__(1.5).zeta().toNumber(), 2.612375348685488, 1e-14)).toBe(true);
  });

  test('eint', () => {
    // Ei(1) ≈ 1.8951178163559...
    expect(approxEqual(R.__call__(1).eint().toNumber(), 1.8951178163559, 1e-4)).toBe(true);

    // Ei(0) = -infinity
    expect(R.__call__(0).eint().toNumber()).toBe(Number.NEGATIVE_INFINITY);
  });

  test('j0 (Bessel J_0)', () => {
    // J_0(0) = 1
    expect(approxEqual(R.__call__(0).j0().toNumber(), 1, 1e-10)).toBe(true);
    // J_0(2) ≈ 0.22389077914124
    expect(approxEqual(R.__call__(2).j0().toNumber(), 0.22389077914124, 1e-8)).toBe(true);
  });

  test('j1 (Bessel J_1)', () => {
    // J_1(0) = 0
    expect(R.__call__(0).j1().toNumber()).toBe(0);
    // J_1(2) ≈ 0.5767248078
    expect(approxEqual(R.__call__(2).j1().toNumber(), 0.5767248077568734, 1e-8)).toBe(true);
  });

  test('y0 (Bessel Y_0)', () => {
    // Y_0(2) ≈ 0.51037567...
    expect(approxEqual(R.__call__(2).y0().toNumber(), 0.5103756726497451, 1e-8)).toBe(true);
    // Y_0 at 0 is -infinity
    expect(R.__call__(0).y0().toNumber()).toBe(Number.NEGATIVE_INFINITY);
  });

  test('y1 (Bessel Y_1)', () => {
    // Y_1(2) ≈ -0.10703243...
    expect(approxEqual(R.__call__(2).y1().toNumber(), -0.10703243154093754, 1e-8)).toBe(true);
  });

  test('jn (Bessel J_n) - downward (Miller) recurrence, n > |x|', () => {
    // This branch used to return NaN for every n > |x| because the downward
    // recurrence was seeded from an uninitialised variable.
    // sage: R(2).jn(3)  -> 0.128943249474402
    expect(approxEqual(R.__call__(2).jn(3).toNumber(), 0.12894324947440206, 1e-14)).toBe(true);
    expect(approxEqual(R.__call__(2).jn(4).toNumber(), 0.033995719807568436, 1e-14)).toBe(true);
    expect(approxEqual(R.__call__(2).jn(5).toNumber(), 0.007039629755871686, 1e-15)).toBe(true);
    expect(approxEqual(R.__call__(1).jn(5).toNumber(), 0.00024975773021123444, 1e-16)).toBe(true);

    // Negative order: J_{-n}(x) = (-1)^n J_n(x)
    // sage: R(2).jn(-17) -> -2.65930780516787e-15
    const jm17 = R.__call__(2).jn(-17).toNumber();
    expect(Math.abs(jm17 / -2.6593078051678734e-15 - 1)).toBeLessThan(1e-6);
    expect(approxEqual(R.__call__(2).jn(-3).toNumber(), -0.12894324947440206, 1e-14)).toBe(true);
    expect(approxEqual(R.__call__(2).jn(-4).toNumber(), 0.033995719807568436, 1e-14)).toBe(true);

    // Negative argument: J_n(-x) = (-1)^n J_n(x)
    expect(approxEqual(R.__call__(-2).jn(3).toNumber(), -0.12894324947440206, 1e-14)).toBe(true);
    expect(approxEqual(R.__call__(-2).jn(4).toNumber(), 0.033995719807568436, 1e-14)).toBe(true);
  });

  test('jn (Bessel J_n) - upward recurrence, n <= |x|', () => {
    // This branch is seeded from j0/j1, which use the Numerical Recipes
    // rational approximations (~1e-8 relative accuracy); mpfr_jn is exact.
    // See DEVIATIONS: Bessel functions use approximations.
    expect(approxEqual(R.__call__(5).jn(2).toNumber(), 0.046565116277752214, 1e-7)).toBe(true);
    expect(approxEqual(R.__call__(10).jn(4).toNumber(), -0.21960268610200853, 1e-7)).toBe(true);

    // jn(0) and jn(1) agree with j0 and j1 exactly.
    expect(R.__call__(2).jn(0).toNumber()).toBe(R.__call__(2).j0().toNumber());
    expect(R.__call__(2).jn(1).toNumber()).toBe(R.__call__(2).j1().toNumber());

    // The three-term recurrence J_{n-1}(x) + J_{n+1}(x) = (2n/x) J_n(x) holds.
    const x = 6;
    for (let n = 1; n <= 4; n++) {
      const lhs =
        R.__call__(x)
          .jn(n - 1)
          .toNumber() +
        R.__call__(x)
          .jn(n + 1)
          .toNumber();
      const rhs = ((2 * n) / x) * R.__call__(x).jn(n).toNumber();
      expect(Math.abs(lhs - rhs)).toBeLessThan(1e-7);
    }
  });

  test('yn (Bessel Y_n)', () => {
    // Y_2(2) ≈ -0.617408104, Y_3(2) ≈ -1.127783777
    expect(approxEqual(R.__call__(2).yn(2).toNumber(), -0.6174081041906827, 1e-7)).toBe(true);
    expect(approxEqual(R.__call__(2).yn(3).toNumber(), -1.1277837768404277, 1e-7)).toBe(true);
    expect(R.__call__(2).yn(0).toNumber()).toBe(R.__call__(2).y0().toNumber());
    expect(R.__call__(2).yn(1).toNumber()).toBe(R.__call__(2).y1().toNumber());
    // Y_{-n}(x) = (-1)^n Y_n(x)
    expect(approxEqual(R.__call__(2).yn(-3).toNumber(), 1.1277837768404277, 1e-7)).toBe(true);
    // Y_n is -infinity at 0 and NaN for negative arguments
    expect(R.__call__(0).yn(3).toNumber()).toBe(Number.NEGATIVE_INFINITY);
    expect(Number.isNaN(R.__call__(-1).yn(3).toNumber())).toBe(true);
  });
});

describe('RealNumber - Float representation', () => {
  const R = new RealField();

  test('is_NaN', () => {
    expect(R.__call__(Number.NaN).is_NaN()).toBe(true);
    expect(R.__call__(0).is_NaN()).toBe(false);
    expect(R.__call__(Number.POSITIVE_INFINITY).is_NaN()).toBe(false);
  });

  test('is_positive_infinity', () => {
    expect(R.__call__(Number.POSITIVE_INFINITY).is_positive_infinity()).toBe(true);
    expect(R.__call__(Number.NEGATIVE_INFINITY).is_positive_infinity()).toBe(false);
    expect(R.__call__(0).is_positive_infinity()).toBe(false);
  });

  test('is_negative_infinity', () => {
    expect(R.__call__(Number.NEGATIVE_INFINITY).is_negative_infinity()).toBe(true);
    expect(R.__call__(Number.POSITIVE_INFINITY).is_negative_infinity()).toBe(false);
    expect(R.__call__(0).is_negative_infinity()).toBe(false);
  });

  test('is_infinity', () => {
    expect(R.__call__(Number.POSITIVE_INFINITY).is_infinity()).toBe(true);
    expect(R.__call__(Number.NEGATIVE_INFINITY).is_infinity()).toBe(true);
    expect(R.__call__(0).is_infinity()).toBe(false);
    expect(R.__call__(Number.NaN).is_infinity()).toBe(false);
  });

  test('is_integer', () => {
    expect(R.__call__(5).is_integer()).toBe(true);
    expect(R.__call__(5).is_integer()).toBe(true);
    expect(R.__call__(5.5).is_integer()).toBe(false);
    expect(R.__call__(0).is_integer()).toBe(true);
  });

  test('is_square', () => {
    expect(R.__call__(4).is_square()).toBe(true);
    expect(R.__call__(0).is_square()).toBe(true);
    expect(R.__call__(-1).is_square()).toBe(false);
  });

  test('conjugate', () => {
    const n = R.__call__(3.14);
    expect(n.conjugate().toNumber()).toBe(3.14);
  });

  test('multiplicative_order', () => {
    // sage: RR(1).multiplicative_order() -> 1
    // sage: RR(-1).multiplicative_order() -> 2
    // sage: RR(3).multiplicative_order() -> +Infinity   (it does not raise)
    expect(R.__call__(1).multiplicative_order()).toBe(1);
    expect(R.__call__(-1).multiplicative_order()).toBe(2);
    expect(R.__call__(2).multiplicative_order()).toBe(Number.POSITIVE_INFINITY);
    expect(R.__call__(3).multiplicative_order()).toBe(Number.POSITIVE_INFINITY);
    expect(R.__call__(0).multiplicative_order()).toBe(Number.POSITIVE_INFINITY);
  });

  test('nextabove and nextbelow', () => {
    const n = R.__call__(1);
    const above = n.nextabove();
    const below = n.nextbelow();

    expect(above.toNumber()).toBeGreaterThan(1);
    expect(below.toNumber()).toBeLessThan(1);

    // The difference should be very small
    expect(above.toNumber() - 1).toBeLessThan(1e-15);
    expect(1 - below.toNumber()).toBeLessThan(1e-15);
  });

  test('nexttoward', () => {
    const n = R.__call__(1);
    expect(n.nexttoward(2).toNumber()).toBeGreaterThan(1);
    expect(n.nexttoward(0).toNumber()).toBeLessThan(1);
    expect(n.nexttoward(1).toNumber()).toBe(1);
  });

  test('ulp', () => {
    // ulp(x) = 2^(exponent(x) - prec); it is the FULL gap to the next float
    // with larger magnitude, not the min of the two neighbour gaps.
    // sage: 1.ulp()    -> 2.22044604925031e-16
    // sage: (-1.5).ulp() -> 2.22044604925031e-16
    expect(R.__call__(1).ulp().toNumber()).toBe(2 ** -52);
    expect(R.__call__(-1.5).ulp().toNumber()).toBe(2 ** -52);
    expect(R.__call__(1.5).ulp().toNumber()).toBe(2 ** -52);
    // pi is in [2,4), so its ulp is twice as large as the ulp of 1.
    expect(R.__call__(Math.PI).ulp().toNumber()).toBe(2 ** -51);

    // sage: a = 1; a + a.ulp() == a -> False ; a + a.ulp()/2 == a -> True
    const a = R.__call__(1);
    expect(a.add(a.ulp()).equals(a)).toBe(false);
    expect(a.add(a.ulp().div(2)).equals(a)).toBe(true);

    // sage: RR(infinity).ulp() -> +infinity ; RR('nan').ulp() -> NaN
    expect(R.__call__(Number.POSITIVE_INFINITY).ulp().toNumber()).toBe(Number.POSITIVE_INFINITY);
    expect(R.__call__(Number.NEGATIVE_INFINITY).ulp().toNumber()).toBe(Number.POSITIVE_INFINITY);
    expect(R.__call__(Number.NaN).ulp().is_NaN()).toBe(true);
  });

  test('epsilon', () => {
    // epsilon(x) = |x| / 2^prec -- it is scale dependent, NOT the constant
    // 2^(1-prec).
    // sage: RR(2^53).epsilon() -> 1
    // sage: RR(0).epsilon()    -> 0
    // sage: RR.pi().epsilon()  -> 3.48786849800863e-16
    expect(R.__call__(1).epsilon().toNumber()).toBe(2 ** -53);
    expect(
      R.__call__(2 ** 53)
        .epsilon()
        .toNumber()
    ).toBe(1);
    expect(R.__call__(0).epsilon().toNumber()).toBe(0);
    expect(approxEqual(R.__call__(Math.PI).epsilon().toNumber(), 3.48786849800863e-16, 1e-30)).toBe(
      true
    );
    // epsilon is even
    expect(R.__call__(-Math.PI).epsilon().toNumber()).toBe(
      R.__call__(Math.PI).epsilon().toNumber()
    );
    // sage: a.epsilon() lies in [a.ulp()/2, a.ulp())
    const pi = R.__call__(Math.PI);
    expect(pi.epsilon().toNumber()).toBeGreaterThanOrEqual(pi.ulp().toNumber() / 2);
    expect(pi.epsilon().toNumber()).toBeLessThan(pi.ulp().toNumber());
    // sage: RR('+Inf').epsilon() -> +infinity ; RR('nan').epsilon() -> NaN
    expect(R.__call__(Number.POSITIVE_INFINITY).epsilon().toNumber()).toBe(
      Number.POSITIVE_INFINITY
    );
    expect(R.__call__(Number.NEGATIVE_INFINITY).epsilon().toNumber()).toBe(
      Number.POSITIVE_INFINITY
    );
    expect(R.__call__(Number.NaN).epsilon().is_NaN()).toBe(true);
  });

  test('sign_mantissa_exponent', () => {
    const n = R.__call__(3.5);
    const [sign, mantissa, exp] = n.sign_mantissa_exponent();
    expect(sign).toBe(1);
    // 3.5 = 7 * 2^-1 = (1 + 0.75) * 2^1
    // In IEEE 754: mantissa = 7 * 2^51, exp = 1 - 52 = -51

    // Verify the reconstruction
    const reconstructed = sign * Number(mantissa) * 2 ** Number(exp);
    expect(reconstructed).toBe(3.5);

    // Test negative number
    const neg = R.__call__(-2.5);
    const [negSign, negMant, negExp] = neg.sign_mantissa_exponent();
    expect(negSign).toBe(-1);
    expect(negSign * Number(negMant) * 2 ** Number(negExp)).toBe(-2.5);
  });

  test('fp_rank', () => {
    // sage: RR(0).fp_rank() -> 0
    // sage: RR(1).fp_rank() -> 20769187434139310514121985316880385   (64-bit)
    // sage: RR(-1).fp_rank() -> -20769187434139310514121985316880385
    // sage: RR(-infinity).fp_rank() -> -41538374868278621023740371006390273
    // These are MPFR ranks, not raw IEEE bit patterns (which would give
    // 4607182418800017408 for RR(1)).
    expect(R.__call__(0).fp_rank()).toBe(0n);
    expect(R.__call__(1).fp_rank()).toBe(20769187434139310514121985316880385n);
    expect(R.__call__(-1).fp_rank()).toBe(-20769187434139310514121985316880385n);
    expect(R.__call__(Number.NEGATIVE_INFINITY).fp_rank()).toBe(
      -41538374868278621023740371006390273n
    );
    expect(R.__call__(Number.POSITIVE_INFINITY).fp_rank()).toBe(
      41538374868278621023740371006390273n
    );

    // sage: RR(1).fp_rank() - RR(1).nextbelow().fp_rank() -> 1
    expect(R.__call__(1).fp_rank() - R.__call__(1).nextbelow().fp_rank()).toBe(1n);
    expect(R.__call__(1).nextabove().fp_rank() - R.__call__(1).fp_rank()).toBe(1n);
    // sage: RR(-infinity).fp_rank() - RR(-infinity).nextabove().fp_rank() -> -1
    // (nextabove of -infinity is -MAX_VALUE for us, so we check pi instead)
    const pi = R.__call__(Math.PI);
    expect(pi.nextabove().fp_rank() - pi.fp_rank()).toBe(1n);

    // Ranks are monotone and odd.
    expect(R.__call__(1).fp_rank()).toBeLessThan(R.__call__(2).fp_rank());
    expect(R.__call__(-1).fp_rank()).toBeLessThan(0n);

    // sage: RR('nan').fp_rank() -> ValueError: Cannot compute fp_rank of NaN
    expect(() => R.__call__(Number.NaN).fp_rank()).toThrow('Cannot compute fp_rank of NaN');
  });

  test('exact_rational', () => {
    const n = R.__call__(0.5);
    const [num, den] = n.exact_rational();
    expect(Number(num) / Number(den)).toBe(0.5);

    const zero = R.__call__(0);
    const [zeroNum, zeroDen] = zero.exact_rational();
    expect(zeroNum).toBe(0n);
    expect(zeroDen).toBe(1n);
  });

  test('nearby_rational with max_denominator', () => {
    // sage doctests of RealNumber.nearby_rational
    expect(R.__call__(0.333).nearby_rational(undefined, 100n)).toEqual([1n, 3n]);
    expect(R.__call__(1 / 3 + 1 / 1000000).nearby_rational(undefined, 2999999n)).toEqual([
      777780n,
      2333333n,
    ]);
    expect(R.__call__(1 / 3 + 1 / 1000000).nearby_rational(undefined, 3000000n)).toEqual([
      1000003n,
      3000000n,
    ]);
    expect(R.__call__(-0.333).nearby_rational(undefined, 1000n)).toEqual([-333n, 1000n]);
    expect(R.__call__(3 / 4).nearby_rational(undefined, 2n)).toEqual([1n, 1n]);
    expect(R.__call__(Math.PI).nearby_rational(undefined, 120n)).toEqual([355n, 113n]);
    expect(R.__call__(Math.PI).nearby_rational(undefined, 10000n)).toEqual([355n, 113n]);
    expect(R.__call__(Math.PI).nearby_rational(undefined, 100000n)).toEqual([312689n, 99532n]);
    expect(R.__call__(Math.PI).nearby_rational(undefined, 1n)).toEqual([3n, 1n]);
    expect(R.__call__(-3.5).nearby_rational(undefined, 1n)).toEqual([-3n, 1n]);
  });

  test('nearby_rational returns the CLOSEST rational under the bound', () => {
    // The result must beat every convergent AND semiconvergent, not just the
    // continued fraction convergents.
    const x = 1 / 3 + 1 / 1000000;
    const [num, den] = R.__call__(x).nearby_rational(undefined, 2999999n);
    const err = Math.abs(Number(num) / Number(den) - x);
    // 777780/2333333 is closer than the convergent 333334/1000001
    expect(err).toBeLessThan(2e-13);
    expect(den).toBeLessThanOrEqual(2999999n);
  });

  test('nearby_rational with max_error', () => {
    // sage: (0.333).nearby_rational(max_error=0.001)   -> 1/3
    // sage: (0.333).nearby_rational(max_error=1)       -> 0
    // sage: (-0.333).nearby_rational(max_error=0.0001) -> -257/772
    expect(R.__call__(0.333).nearby_rational(0.001)).toEqual([1n, 3n]);
    expect(R.__call__(0.333).nearby_rational(1)).toEqual([0n, 1n]);
    expect(R.__call__(-0.333).nearby_rational(0.0001)).toEqual([-257n, 772n]);
  });

  test('nearby_rational requires exactly one bound', () => {
    // sage raises ValueError when neither or both are given
    expect(() => R.__call__(1.5).nearby_rational()).toThrow(
      'Must specify exactly one of max_error or max_denominator in nearby_rational()'
    );
    expect(() => R.__call__(1.5).nearby_rational(0.1, 100n)).toThrow(
      'Must specify exactly one of max_error or max_denominator in nearby_rational()'
    );
    expect(() => R.__call__(Number.NaN).nearby_rational(undefined, 1000n)).toThrow(
      'cannot convert NaN or infinity to rational number'
    );
    expect(() => R.__call__(Number.POSITIVE_INFINITY).nearby_rational(0.01)).toThrow(
      'cannot convert NaN or infinity to rational number'
    );
  });

  test('simplest_rational', () => {
    // sage doctests of RealNumber.simplest_rational
    expect(R.__call__(1 / 3).simplest_rational()).toEqual([1n, 3n]);
    expect(R.__call__(-1 / 3).simplest_rational()).toEqual([-1n, 3n]);
    expect(R.__call__(Math.PI).simplest_rational()).toEqual([245850922n, 78256779n]);
    expect(R.__call__(Math.SQRT2).simplest_rational()).toEqual([131836323n, 93222358n]);
    expect(R.__call__(1234).simplest_rational()).toEqual([1234n, 1n]);
    expect(R.__call__(2 ** -210).simplest_rational()).toEqual([
      1n,
      1645504557321205859467264516194506011931735427766374553794641921n,
    ]);
    expect(R.__call__(2 ** 210).simplest_rational()).toEqual([
      1645504557321205950811116849375918117252433820865891134852825088n,
      1n,
    ]);
    // sage: (RR(17).sqrt()).simplest_rational()^2 - 17 -> -1/348729667233025
    const [n17, d17] = R.__call__(Math.sqrt(17)).simplest_rational();
    expect(n17 * n17 - 17n * d17 * d17).toBe(-1n);
    expect(d17 * d17).toBe(348729667233025n);
  });

  test('simplest_rational round-trips back to self', () => {
    // The defining property: the result must be equal to self when coerced
    // back into this field.  A plain continued fraction truncation is not.
    for (const v of [
      Math.PI,
      Math.SQRT2,
      Math.E,
      1 / 3,
      -1 / 3,
      0.1,
      1234,
      1e-30,
      1e30,
      Math.sqrt(17),
      Math.cbrt(23),
    ]) {
      const [num, den] = R.__call__(v).simplest_rational();
      expect(Number(num) / Number(den)).toBe(v);
    }
  });

  test('simplest_rational rejects NaN and infinity', () => {
    expect(() => R.__call__(Number.NaN).simplest_rational()).toThrow(
      'cannot convert NaN or infinity to rational number'
    );
    expect(() => R.__call__(Number.NEGATIVE_INFINITY).simplest_rational()).toThrow(
      'cannot convert NaN or infinity to rational number'
    );
  });
});

describe('RR factory', () => {
  test('creates RealField', () => {
    const R = RR();
    expect(R.precision()).toBe(53);

    const R100 = RR(100);
    expect(R100.precision()).toBe(100);
  });
});

describe('mpfr_prec functions', () => {
  test('mpfr_prec_min', () => {
    expect(mpfr_prec_min()).toBe(1);
  });

  test('mpfr_prec_max', () => {
    expect(mpfr_prec_max()).toBe(9223372036854775551n);
  });
});

describe('RealField - complex_field', () => {
  test('returns ComplexField with same precision', () => {
    const R = new RealField(100);
    const C = R.complex_field();
    expect(C.prec()).toBe(100);
  });

  test('default precision returns 53-bit ComplexField', () => {
    const R = new RealField();
    const C = R.complex_field();
    expect(C.prec()).toBe(53);
  });
});

describe('RealNumber - algebraic_dependency', () => {
  const R = new RealField();

  test('sqrt(2) has algebraic dependency x^2 - 2', () => {
    const sqrt2 = R.__call__(Math.sqrt(2));
    const poly = sqrt2.algebraic_dependency(5);

    // The polynomial should be approximately x^2 - 2
    // Find the ratio to normalize
    let foundQuadratic = false;
    if (poly.length >= 3 && poly[2] !== 0n) {
      // Check if it's x^2 - 2 or a multiple
      const ratio = poly[2]!;
      const expected0 = -2n * ratio;
      const expected1 = 0n;

      foundQuadratic = poly[0] === expected0 && (poly[1] === expected1 || poly[1] === 0n);
    }
    expect(foundQuadratic).toBe(true);
  });

  test('integer values return simple polynomial', () => {
    const three = R.__call__(3);
    const poly = three.algebraic_dependency(2);

    // Should return x - 3 or equivalent
    // poly[0] + poly[1] * x = 0 when x = 3
    // So poly[0] + 3 * poly[1] = 0
    const val = poly[0]! + 3n * (poly[1] ?? 0n);
    expect(val).toBe(0n);
  });

  test('phi (golden ratio) has x^2 - x - 1', () => {
    const phi = R.__call__((1 + Math.sqrt(5)) / 2);
    const poly = phi.algebraic_dependency(5);

    // Check that the polynomial evaluates to approximately 0 at phi
    const phiVal = (1 + Math.sqrt(5)) / 2;
    let evalResult = 0;
    for (let i = 0; i < poly.length; i++) {
      evalResult += Number(poly[i]!) * phiVal ** i;
    }
    expect(Math.abs(evalResult)).toBeLessThan(1e-8);
  });

  test('algdep is alias for algebraic_dependency', () => {
    const sqrt3 = R.__call__(Math.sqrt(3));
    const poly1 = sqrt3.algebraic_dependency(3);
    const poly2 = sqrt3.algdep(3);
    expect(poly1).toEqual(poly2);
  });
});
