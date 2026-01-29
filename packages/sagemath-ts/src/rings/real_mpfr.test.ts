/**
 * Unit tests for sage/rings/real_mpfr
 * Tests for real number transcendental functions
 */
import { describe, expect, test } from 'bun:test';
import { RealField, RealNumber, RR, RoundingMode, mpfr_prec_min, mpfr_prec_max } from './real_mpfr.js';

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
    expect(approxEqual(catalan.toNumber(), 0.9159655941772190)).toBe(true);
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
    expect(R.__call__(3.0).floor()).toBe(3n);
  });

  test('ceil', () => {
    expect(R.__call__(3.2).ceil()).toBe(4n);
    expect(R.__call__(-3.2).ceil()).toBe(-3n);
    expect(R.__call__(3.0).ceil()).toBe(3n);
  });

  test('round', () => {
    expect(R.__call__(3.4).round()).toBe(3n);
    expect(R.__call__(3.5).round()).toBe(4n);
    expect(R.__call__(-3.5).round()).toBe(-3n);
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
    expect(R.__call__(0).log().toNumber()).toBe(-Infinity);
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
    expect(approxEqual(R.__call__(Math.E - 1).log1p().toNumber(), 1)).toBe(true);
    expect(R.__call__(-1).log1p().toNumber()).toBe(-Infinity);
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
    expect(approxEqual(R.__call__(pi / 2).sin().toNumber(), 1)).toBe(true);
    expect(approxEqual(R.__call__(pi).sin().toNumber(), 0)).toBe(true);
    expect(approxEqual(R.__call__(pi / 6).sin().toNumber(), 0.5)).toBe(true);
  });

  test('cos', () => {
    expect(R.__call__(0).cos().toNumber()).toBe(1);
    expect(approxEqual(R.__call__(pi / 2).cos().toNumber(), 0)).toBe(true);
    expect(approxEqual(R.__call__(pi).cos().toNumber(), -1)).toBe(true);
    expect(approxEqual(R.__call__(pi / 3).cos().toNumber(), 0.5)).toBe(true);
  });

  test('tan', () => {
    expect(R.__call__(0).tan().toNumber()).toBe(0);
    expect(approxEqual(R.__call__(pi / 4).tan().toNumber(), 1)).toBe(true);
    expect(approxEqual(R.__call__(pi / 6).tan().toNumber(), 1 / Math.sqrt(3))).toBe(true);
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
    expect(approxEqual(R.__call__(pi / 4).cot().toNumber(), 1)).toBe(true);
    expect(approxEqual(R.__call__(pi / 6).cot().toNumber(), Math.sqrt(3))).toBe(true);
  });

  test('sec', () => {
    expect(R.__call__(0).sec().toNumber()).toBe(1);
    expect(approxEqual(R.__call__(pi / 3).sec().toNumber(), 2)).toBe(true);
  });

  test('csc', () => {
    expect(approxEqual(R.__call__(pi / 2).csc().toNumber(), 1)).toBe(true);
    expect(approxEqual(R.__call__(pi / 6).csc().toNumber(), 2)).toBe(true);
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
    expect(approxEqual(R.__call__(0).erfc().toNumber(), 1, 1e-6)).toBe(true);
    expect(approxEqual(R.__call__(1).erfc().toNumber(), 0.1572992070502851, 1e-6)).toBe(true);
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
    expect(approxEqual(R.__call__(2).zeta().toNumber(), Math.PI ** 2 / 6, 1e-6)).toBe(true);

    // zeta(1) = infinity
    expect(R.__call__(1).zeta().toNumber()).toBe(Infinity);

    // zeta(-2) = 0
    expect(R.__call__(-2).zeta().toNumber()).toBe(0);
  });

  test('eint', () => {
    // Ei(1) ≈ 1.8951178163559...
    expect(approxEqual(R.__call__(1).eint().toNumber(), 1.8951178163559, 1e-4)).toBe(true);

    // Ei(0) = -infinity
    expect(R.__call__(0).eint().toNumber()).toBe(-Infinity);
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
    expect(R.__call__(0).y0().toNumber()).toBe(-Infinity);
  });

  test('y1 (Bessel Y_1)', () => {
    // Y_1(2) ≈ -0.10703243...
    expect(approxEqual(R.__call__(2).y1().toNumber(), -0.10703243154093754, 1e-8)).toBe(true);
  });
});

describe('RealNumber - Float representation', () => {
  const R = new RealField();

  test('is_NaN', () => {
    expect(R.__call__(NaN).is_NaN()).toBe(true);
    expect(R.__call__(0).is_NaN()).toBe(false);
    expect(R.__call__(Infinity).is_NaN()).toBe(false);
  });

  test('is_positive_infinity', () => {
    expect(R.__call__(Infinity).is_positive_infinity()).toBe(true);
    expect(R.__call__(-Infinity).is_positive_infinity()).toBe(false);
    expect(R.__call__(0).is_positive_infinity()).toBe(false);
  });

  test('is_negative_infinity', () => {
    expect(R.__call__(-Infinity).is_negative_infinity()).toBe(true);
    expect(R.__call__(Infinity).is_negative_infinity()).toBe(false);
    expect(R.__call__(0).is_negative_infinity()).toBe(false);
  });

  test('is_infinity', () => {
    expect(R.__call__(Infinity).is_infinity()).toBe(true);
    expect(R.__call__(-Infinity).is_infinity()).toBe(true);
    expect(R.__call__(0).is_infinity()).toBe(false);
    expect(R.__call__(NaN).is_infinity()).toBe(false);
  });

  test('is_integer', () => {
    expect(R.__call__(5).is_integer()).toBe(true);
    expect(R.__call__(5.0).is_integer()).toBe(true);
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
    expect(R.__call__(1).multiplicative_order()).toBe(1);
    expect(R.__call__(-1).multiplicative_order()).toBe(2);
    expect(() => R.__call__(2).multiplicative_order()).toThrow();
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
    const n = R.__call__(1);
    const u = n.ulp();
    expect(u.toNumber()).toBeGreaterThan(0);
    expect(u.toNumber()).toBeLessThan(1e-14);
  });

  test('epsilon', () => {
    const n = R.__call__(1);
    const eps = n.epsilon();
    // For 53-bit precision, epsilon ≈ 2^-52
    expect(approxEqual(eps.toNumber(), Math.pow(2, -52), 1e-20)).toBe(true);
  });

  test('sign_mantissa_exponent', () => {
    const n = R.__call__(3.5);
    const [sign, mantissa, exp] = n.sign_mantissa_exponent();
    expect(sign).toBe(1);
    // 3.5 = 7 * 2^-1 = (1 + 0.75) * 2^1
    // In IEEE 754: mantissa = 7 * 2^51, exp = 1 - 52 = -51

    // Verify the reconstruction
    const reconstructed = sign * Number(mantissa) * Math.pow(2, Number(exp));
    expect(reconstructed).toBe(3.5);

    // Test negative number
    const neg = R.__call__(-2.5);
    const [negSign, negMant, negExp] = neg.sign_mantissa_exponent();
    expect(negSign).toBe(-1);
    expect(negSign * Number(negMant) * Math.pow(2, Number(negExp))).toBe(-2.5);
  });

  test('fp_rank', () => {
    const zero = R.__call__(0);
    expect(zero.fp_rank()).toBe(0n);

    // Positive numbers have increasing rank
    const one = R.__call__(1);
    const two = R.__call__(2);
    expect(one.fp_rank()).toBeLessThan(two.fp_rank());

    // Negative numbers have negative rank
    const negOne = R.__call__(-1);
    expect(negOne.fp_rank()).toBeLessThan(0n);
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

  test('nearby_rational', () => {
    // pi approximation
    const pi = R.__call__(Math.PI);
    const [num, den] = pi.nearby_rational(undefined, 1000n);
    const approx = Number(num) / Number(den);
    expect(Math.abs(approx - Math.PI)).toBeLessThan(0.001);
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

      foundQuadratic =
        poly[0] === expected0 &&
        (poly[1] === expected1 || poly[1] === 0n);
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
      evalResult += Number(poly[i]!) * Math.pow(phiVal, i);
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
