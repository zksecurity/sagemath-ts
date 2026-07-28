/**
 * Tests for p-adic element implementation
 *
 * @see Reference: sage/rings/padics/padic_generic_element.pyx
 */

import { describe, expect, test } from 'bun:test';
import { NotImplementedError, ValueError, ZeroDivisionError } from '../../errors.js';
import { Qp, Zp, pAdicField, pAdicRing } from './padic_generic.js';
import { PrecisionError, pAdicGenericElement } from './padic_generic_element.js';

const INFINITY = Number.POSITIVE_INFINITY;

describe('pAdicRing and pAdicField construction', () => {
  test('Zp construction', () => {
    const R = Zp(5n);
    expect(R.prime()).toBe(5n);
    expect(R.precision_cap()).toBe(20);
    expect(R.is_field()).toBe(false);
  });

  test('Qp construction', () => {
    const K = Qp(7n, 10);
    expect(K.prime()).toBe(7n);
    expect(K.precision_cap()).toBe(10);
    expect(K.is_field()).toBe(true);
  });

  test('invalid prime throws', () => {
    expect(() => new pAdicRing(4n)).toThrow(ValueError);
    expect(() => new pAdicRing(1n)).toThrow(ValueError);
  });
});

describe('pAdicGenericElement basic operations', () => {
  const R = Zp(5n, 10);

  test('element creation', () => {
    const a = R.__call__(17n);
    expect(a.prime()).toBe(5n);
    expect(a.valuation()).toBe(0n);
    expect(a.lift()).toBe(17n);
  });

  test('zero element', () => {
    const zero = R.zero();
    expect(zero.is_zero()).toBe(true);
  });

  test('one element', () => {
    const one = R.one();
    expect(one.is_one()).toBe(true);
    expect(one.is_unit()).toBe(true);
  });

  test('valuation', () => {
    const a = R.__call__(25n); // 5^2
    expect(a.valuation()).toBe(2n);

    const b = R.__call__(125n); // 5^3
    expect(b.valuation()).toBe(3n);

    const c = R.__call__(7n); // coprime to 5
    expect(c.valuation()).toBe(0n);
  });

  test('unit part', () => {
    const a = R.__call__(50n); // 2 * 5^2
    const unit = a.unit_part();
    expect(unit.valuation()).toBe(0n);
    expect(unit.lift() % 5n).toBe(2n);
  });
});

describe('p-adic arithmetic', () => {
  const R = Zp(7n, 10);

  test('addition', () => {
    const a = R.__call__(3n);
    const b = R.__call__(4n);
    const sum = a.add(b);
    expect(sum.lift() % 7n ** 10n).toBe(7n);
  });

  test('subtraction', () => {
    const a = R.__call__(10n);
    const b = R.__call__(3n);
    const diff = a.sub(b);
    expect(diff.lift()).toBe(7n);
  });

  test('multiplication', () => {
    const a = R.__call__(3n);
    const b = R.__call__(4n);
    const prod = a.mul(b);
    expect(prod.lift()).toBe(12n);
  });

  test('division', () => {
    const a = R.__call__(6n);
    const b = R.__call__(2n);
    const quot = a.div(b);
    expect(quot.lift()).toBe(3n);
  });

  test('negation', () => {
    const a = R.__call__(3n);
    const neg = a.neg();
    const sum = a.add(neg);
    expect(sum.is_zero()).toBe(true);
  });

  test('power', () => {
    const a = R.__call__(2n);
    const pow = a.pow(3n);
    expect(pow.lift()).toBe(8n);
  });

  test('division by zero throws', () => {
    const a = R.__call__(5n);
    const zero = R.zero();
    expect(() => a.div(zero)).toThrow();
  });
});

describe('p-adic precision', () => {
  const R = Zp(3n, 5);

  test('absolute precision', () => {
    const a = R.__call__(9n); // 3^2
    expect(a.precision_absolute()).toBe(7); // valuation 2 + relprec 5
  });

  test('relative precision', () => {
    const a = R.__call__(9n);
    expect(a.precision_relative()).toBe(5);
  });

  test('add_bigoh reduces precision', () => {
    const a = R.__call__(1n);
    const b = a.add_bigoh(3);
    expect(b.precision_absolute()).toBe(3);
  });

  test('lift_to_precision increases precision', () => {
    const a = R.__call__(1n);
    const b = a.lift_to_precision(10);
    expect(b.precision_absolute()).toBeGreaterThanOrEqual(10);
  });
});

describe('p-adic expansion', () => {
  const R = Zp(5n, 5);

  test('expansion of integer', () => {
    const a = R.__call__(23n); // 23 = 3 + 4*5
    const exp = a.expansion();
    expect(exp[0]).toBe(3n);
    expect(exp[1]).toBe(4n);
    expect(exp[2]).toBe(0n);
  });

  test('residue', () => {
    const a = R.__call__(23n);
    expect(a.residue(1)).toBe(3n);
    expect(a.residue(2)).toBe(23n);
  });

  test('__getitem__', () => {
    const a = R.__call__(23n);
    expect(a.__getitem__(0)).toBe(3n);
    expect(a.__getitem__(1)).toBe(4n);
  });
});

describe('p-adic predicates', () => {
  const R = Zp(5n, 10);

  test('is_zero', () => {
    expect(R.zero().is_zero()).toBe(true);
    expect(R.one().is_zero()).toBe(false);
    expect(R.__call__(5n).is_zero()).toBe(false);
  });

  test('is_one', () => {
    expect(R.one().is_one()).toBe(true);
    expect(R.zero().is_one()).toBe(false);
    expect(R.__call__(2n).is_one()).toBe(false);
  });

  test('is_unit', () => {
    expect(R.__call__(3n).is_unit()).toBe(true);
    expect(R.__call__(5n).is_unit()).toBe(false);
    expect(R.zero().is_unit()).toBe(false);
  });

  test('is_integral', () => {
    expect(R.__call__(1n).is_integral()).toBe(true);
    expect(R.__call__(5n).is_integral()).toBe(true);
    expect(R.zero().is_integral()).toBe(true);
  });
});

describe('p-adic square root', () => {
  const R = Zp(7n, 10);

  test('is_square for squares', () => {
    const a = R.__call__(4n);
    expect(a.is_square()).toBe(true);

    const b = R.__call__(9n);
    expect(b.is_square()).toBe(true);

    const c = R.__call__(49n); // 7^2
    expect(c.is_square()).toBe(true);
  });

  test('is_square for non-squares', () => {
    const a = R.__call__(3n); // 3 is not a QR mod 7
    expect(a.is_square()).toBe(false);
  });

  test('sqrt of square', () => {
    const a = R.__call__(4n);
    const sqrt = a.sqrt();
    // sqrt^2 should equal a
    expect(sqrt.pow(2n).eq(a)).toBe(true);
  });

  // SageMath's square_root(extend=True) (the default) raises NotImplementedError
  // for a non-square, because it would have to move to an extension field; only
  // extend=False produces the ValueError.
  // Reference: sage/rings/padics/padic_generic_element.pyx:square_root
  //   sage: Zp(3,20)(2).square_root(extend=False)
  //   ValueError: element is not a square
  test('sqrt of non-square throws', () => {
    const a = R.__call__(3n);
    expect(() => a.sqrt()).toThrow(NotImplementedError);
    expect(() => a.sqrt({ extend: false })).toThrow(ValueError);
    expect(a.sqrt({ extend: false, all: true })).toEqual([]);
  });

  test('sqrt of zero', () => {
    const zero = R.zero();
    const sqrt = zero.sqrt();
    expect(sqrt.is_zero()).toBe(true);
  });
});

describe('p-adic Teichmuller lift', () => {
  const R = Zp(5n, 10);

  test('Teichmuller of residue', () => {
    const t = R.teichmuller(2n);
    // Teichmuller lift satisfies t^(p-1) = 1
    const pow = t.pow(4n); // 5-1 = 4
    expect(pow.is_one()).toBe(true);
  });

  test('Teichmuller system', () => {
    const system = R.teichmuller_system();
    expect(system.length).toBe(4); // p-1 = 4 elements

    // Each should satisfy t^(p-1) = 1
    for (const t of system) {
      expect(t.pow(4n).is_one()).toBe(true);
    }
  });

  test('Teichmuller of zero is zero', () => {
    const t = R.teichmuller(0n);
    expect(t.is_zero()).toBe(true);
  });
});

describe('p-adic equality', () => {
  const R = Zp(5n, 10);

  test('equal elements', () => {
    const a = R.__call__(17n);
    const b = R.__call__(17n);
    expect(a.eq(b)).toBe(true);
  });

  test('different elements', () => {
    const a = R.__call__(17n);
    const b = R.__call__(18n);
    expect(a.eq(b)).toBe(false);
  });

  test('equality up to precision', () => {
    const a = R.__call__(1n);
    const b = R.__call__(1n + 5n ** 10n); // Differ only at precision 10
    // They should be equal since default precision is 10
    // Actually for our implementation they might differ
    // Let's use add_bigoh to make precision explicit
    const a5 = a.add_bigoh(5);
    const b5 = b.add_bigoh(5);
    expect(a5.eq(b5)).toBe(true);
  });
});

describe('p-adic string representation', () => {
  const R = Zp(5n, 5);

  test('toString of simple element', () => {
    const a = R.__call__(3n);
    const str = a.toString();
    expect(str).toContain('3');
    expect(str).toContain('O(5^');
  });

  test('toString of zero', () => {
    const zero = R.zero();
    const str = zero.toString();
    expect(str).toBe('0');
  });
});

describe('p-adic field operations', () => {
  const K = Qp(5n, 10);

  test('field is_field', () => {
    expect(K.is_field()).toBe(true);
  });

  test('fraction_field returns self', () => {
    expect(K.fraction_field()).toBe(K);
  });

  test('integer_ring returns Zp', () => {
    const R = K.integer_ring();
    expect(R.is_field()).toBe(false);
    expect(R.prime()).toBe(5n);
  });
});

describe('p-adic exponential and logarithm', () => {
  const R = Zp(5n, 10);

  test('exp converges for p*x', () => {
    const x = R.__call__(5n); // v(x) = 1 > e/(p-1) = 1/4
    // exp should converge
    const expX = x.exp();
    expect(expX.valuation()).toBe(0n);
  });

  test('exp does not converge for unit', () => {
    const x = R.__call__(1n); // v(x) = 0
    expect(() => x.exp()).toThrow(ValueError);
  });

  test('log of 1-unit', () => {
    const x = R.__call__(1n + 5n); // 1 + 5 is a 1-unit
    const logX = x.log();
    // log(1+5) should have valuation 1
    expect(logX.valuation()).toBeGreaterThan(0n);
  });

  test('log of zero throws', () => {
    expect(() => R.zero().log()).toThrow(ValueError);
  });
});

describe('p-adic norm and trace', () => {
  const R = Zp(7n, 10);

  test('norm returns self for base field', () => {
    const a = R.__call__(5n);
    expect(a.norm().eq(a)).toBe(true);
  });

  test('trace returns self for base field', () => {
    const a = R.__call__(5n);
    expect(a.trace().eq(a)).toBe(true);
  });
});

describe('p-adic abs', () => {
  const R = Zp(5n, 10);

  test('abs of unit', () => {
    const a = R.__call__(3n);
    expect(a.abs()).toBe(1);
  });

  test('abs of p', () => {
    const a = R.__call__(5n);
    expect(a.abs()).toBe(0.2); // 5^(-1) = 0.2
  });

  test('abs of p^2', () => {
    const a = R.__call__(25n);
    expect(a.abs()).toBe(0.04); // 5^(-2) = 0.04
  });

  test('abs of zero', () => {
    expect(R.zero().abs()).toBe(0);
  });
});

// ---------------------------------------------------------------------------
// Value assertions ported from SageMath doctests.
//
// The tests above only exercise valuation-0 operands and check valuations
// rather than values, which is how the p^v double-multiplication in add() and
// the wrong nth_root/exp/log values went unnoticed.
// ---------------------------------------------------------------------------

describe('p-adic arithmetic with nonzero valuation', () => {
  test('addition aligns at the common valuation', () => {
    const R = Zp(5n, 10);
    expect(R.__call__(10n).add(R.__call__(15n)).lift()).toBe(25n);
    expect(R.__call__(25n).add(R.__call__(1n)).toString()).toBe('1 + 5^2 + O(5^10)');
    expect(R.__call__(25n).add(R.__call__(1n)).lift()).toBe(26n);

    const S = Zp(7n, 10);
    expect(S.__call__(7n).add(S.__call__(14n)).lift()).toBe(21n);
    expect(S.__call__(49n).add(S.__call__(343n)).lift()).toBe(392n);
  });

  test('subtraction with nonzero valuation', () => {
    const R = Zp(5n, 10);
    expect(R.__call__(125n).sub(R.__call__(25n)).lift()).toBe(100n);
    expect(R.__call__(25n).sub(R.__call__(25n)).is_zero()).toBe(true);
  });

  test('multiplication and division with nonzero valuation', () => {
    const R = Zp(5n, 10);
    expect(R.__call__(25n).mul(R.__call__(15n)).lift()).toBe(375n);
    expect(R.__call__(375n).div(R.__call__(25n)).lift()).toBe(15n);
    // 375 = 3*5^3 and 25 = 5^2, so the quotient is 15 = 3*5 with valuation 1
    expect(R.__call__(375n).div(R.__call__(25n)).valuation()).toBe(1n);
  });

  test('powers with nonzero valuation', () => {
    const R = Zp(5n, 10);
    expect(R.__call__(10n).pow(3n).lift()).toBe(1000n);
    expect(R.__call__(10n).pow(3n).valuation()).toBe(3n);
  });

  test('negative valuations in Qp', () => {
    const K = Qp(5n, 20);
    const inv25 = K.one().div(K.__call__(25n));
    expect(inv25.valuation()).toBe(-2n);
    // SageMath: Qp(5)(1/25) prints 5^-2 + O(5^18)
    expect(inv25.toString()).toBe('5^-2 + O(5^18)');
    expect(inv25.add(K.__call__(5n)).toString()).toBe('5^-2 + 5 + O(5^18)');
    expect(inv25.mul(K.__call__(25n)).is_one()).toBe(true);
  });
});

describe('p-adic expansion (SageMath doctests)', () => {
  test('ring expansion starts at p^0', () => {
    // sage: R = Zp(7,6); a = R(12837162817); a
    // 3 + 4*7 + 4*7^2 + 4*7^4 + O(7^6)
    const R = Zp(7n, 6);
    const a = R.__call__(12837162817n);
    expect(a.expansion()).toEqual([3n, 4n, 4n, 0n, 4n, 0n]);
    expect(a.toString()).toBe('3 + 4*7 + 4*7^2 + 4*7^4 + O(7^6)');

    const b = Zp(7n, 4).__call__(6n * 7n + 49n);
    expect(b.expansion()).toEqual([0n, 6n, 1n, 0n, 0n]);
  });

  test('field expansion starts at p^valuation', () => {
    const K = Qp(7n, 4);
    const a = K.__call__(6n * 7n + 49n);
    expect(a.expansion()).toEqual([6n, 1n, 0n, 0n]);
  });

  test('expansion of zero is empty', () => {
    const R = Zp(5n, 10);
    expect(R.zero().expansion()).toEqual([]);
    expect(R.zero().add_bigoh(4).expansion()).toEqual([]);
  });
});

describe('p-adic slice (SageMath doctests)', () => {
  // sage: R = Zp(5, 6, 'capped-rel'); a = R(1/2)
  const R = Zp(5n, 6);
  const a = R.__call__((5n ** 6n + 1n) / 2n); // 1/2 in Zp(5,6)

  test('R(1/2) expansion', () => {
    expect(a.toString()).toBe('3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + O(5^6)');
  });

  test('slice keeps the p^i weighting', () => {
    // sage: a.slice(2, 4) -> 2*5^2 + 2*5^3 + O(5^4)
    expect(a.slice(2, 4).toString()).toBe('2*5^2 + 2*5^3 + O(5^4)');
    expect(a.slice(2, 4).lift()).toBe(300n);
  });

  test('slice with a step', () => {
    expect(a.slice(1, 6, 2).toString()).toBe('2*5 + 2*5^3 + 2*5^5 + O(5^6)');
    expect(a.slice(0, 5, 2).toString()).toBe('3 + 2*5^2 + 2*5^4 + O(5^5)');
    expect(a.slice(0, 6, 2).toString()).toBe('3 + 2*5^2 + 2*5^4 + O(5^6)');
    expect(a.slice(0, 7, 2).toString()).toBe('3 + 2*5^2 + 2*5^4 + O(5^6)');
  });

  test('slice step must be positive', () => {
    expect(() => a.slice(0, 3, 0)).toThrow(ValueError);
    expect(() => a.slice(0, 3, -1)).toThrow(ValueError);
  });

  test('empty slices carry the precision given by j', () => {
    expect(a.slice(5, 4).toString()).toBe('O(5^4)');
    expect(a.slice(6, 5).toString()).toBe('O(5^5)');
    expect(a.slice(101, 100).toString()).toBe('O(5^6)');
  });

  test('slices over fields', () => {
    const K = Qp(5n, 6);
    const x = K.one().div(K.__call__(25n)); // 1/25
    const b = K.__call__(25n);
    expect(x.toString()).toBe('5^-2 + O(5^4)');
    expect(b.toString()).toBe('5^2 + O(5^8)');
    expect(x.slice(2, 4).toString()).toBe('O(5^4)');
    expect(b.slice(2, 4).toString()).toBe('5^2 + O(5^4)');
    expect(x.slice(-3, -1).toString()).toBe('5^-2 + O(5^-1)');
    expect(b.slice(-1, 1).toString()).toBe('O(5)');
    expect(b.slice(-3, -1).toString()).toBe('O(5^-1)');
    expect(b.slice(101, 100).toString()).toBe('O(5^8)');
    expect(b.slice(0, 7, 2).toString()).toBe('5^2 + O(5^7)');
    expect(b.slice(0, 9, 2).toString()).toBe('5^2 + O(5^8)');
  });

  test('slice with i = null starts at the valuation', () => {
    const K = Qp(5n, 6);
    const x = K.one().div(K.__call__(25n)).add(K.__call__(5n));
    expect(x.toString()).toBe('5^-2 + 5 + O(5^4)');
    expect(x.slice(null, 3).toString()).toBe('5^-2 + 5 + O(5^3)');

    const S = Zp(5n, 7);
    expect(S.__call__(300n).slice(null, 5).toString()).toBe('2*5^2 + 2*5^3 + O(5^5)');
  });

  test('slice of an exact zero is an exact zero', () => {
    expect(Qp(3n).zero().slice(0, null).toString()).toBe('0');
  });
});

describe('p-adic valuation and orders of zero', () => {
  const R = Zp(5n, 10);

  test('valuation of an exact zero is +Infinity', () => {
    expect(R.zero().valuation()).toBe(INFINITY);
    expect(R.zero().ordp()).toBe(INFINITY);
  });

  test('residue and exp of an exact zero', () => {
    expect(R.zero().residue()).toBe(0n);
    expect(R.zero().exp().toString()).toBe('1 + O(5^10)');
  });

  test('additive_order (SageMath doctests)', () => {
    // sage: R = Zp(7, 4); a = R(7^3); a.additive_order(3) -> 1 ; a.additive_order(4) -> +Infinity
    const S = Zp(7n, 4);
    expect(S.__call__(343n).additive_order(3)).toBe(1n);
    expect(S.__call__(343n).additive_order(4)).toBe(INFINITY);
    expect(R.zero().additive_order()).toBe(1n);
    // An inexact zero is indistinguishable from zero, so its additive order is 1
    expect(R.zero().add_bigoh(5).additive_order()).toBe(1n);
    expect(R.__call__(3n).additive_order()).toBe(INFINITY);
  });

  test('multiplicative_order (SageMath doctests)', () => {
    // sage: K = Qp(5,20,'capped-rel')
    const K = Qp(5n, 20);
    expect(K.__call__(-1n).multiplicative_order(20)).toBe(2n);
    expect(K.__call__(1n).multiplicative_order(20)).toBe(1n);
    expect(K.__call__(2n).multiplicative_order(20)).toBe(INFINITY);
    expect(K.__call__(5n).multiplicative_order(20)).toBe(INFINITY);
    expect(K.one().div(K.__call__(5n)).multiplicative_order(20)).toBe(INFINITY);
    // K.zeta() is the primitive 4th root of unity
    expect(K.teichmuller(2n).multiplicative_order(20)).toBe(4n);

    const S = Zp(5n, 20);
    for (const v of [2n, 3n, 4n, 5n, 25n]) {
      expect(S.__call__(v).multiplicative_order(20)).toBe(INFINITY);
    }
    expect(S.__call__(-1n).multiplicative_order(20)).toBe(2n);
    // p = 2 is the case where p-adic roots of unity of p-power order exist
    expect(Zp(2n, 20).__call__(-1n).multiplicative_order()).toBe(2n);
    expect(Zp(2n, 20).__call__(3n).multiplicative_order()).toBe(INFINITY);
  });

  test('is_square distinguishes exact from inexact zero', () => {
    expect(R.zero().is_square()).toBe(true);
    expect(() => R.zero().add_bigoh(5).is_square()).toThrow(PrecisionError);
  });
});

describe('p-adic square roots (SageMath doctests)', () => {
  test('odd p', () => {
    // sage: R = Zp(3, 20)
    const R = Zp(3n, 20);
    expect(R.zero().square_root().toString()).toBe('0');
    expect(R.__call__(1n).square_root().toString()).toBe('1 + O(3^20)');
    expect(R.__call__(4n).square_root().neg().toString()).toBe('2 + O(3^20)');
    expect(R.__call__(9n).square_root().toString()).toBe('3 + O(3^21)');
    expect(() => R.__call__(2n).square_root({ extend: false })).toThrow(ValueError);
  });

  test('p = 2 loses one digit of relative precision', () => {
    // sage: R2 = Zp(2, 20); R2(1).square_root() -> 1 + O(2^19); R2(4).square_root() -> 2 + O(2^20)
    const R2 = Zp(2n, 20);
    expect(R2.__call__(1n).square_root().toString()).toBe('1 + O(2^19)');
    expect(R2.__call__(4n).square_root().toString()).toBe('2 + O(2^20)');
    // sage: Z2(17).square_root()
    expect(R2.__call__(17n).square_root().toString()).toBe(
      '1 + 2^3 + 2^5 + 2^6 + 2^7 + 2^9 + 2^10 + 2^13 + 2^16 + 2^17 + O(2^19)'
    );
    expect(R2.__call__(9n).square_root().pow(2n).eq(R2.__call__(9n).add_bigoh(19))).toBe(true);
  });

  test('is_square for p = 2 (SageMath doctests)', () => {
    const R2 = Zp(2n, 20);
    const expected: Record<string, boolean> = {
      '0': true,
      '1': true,
      '2': false,
      '3': false,
      '4': true,
      '5': false,
      '6': false,
      '7': false,
      '8': false,
      '9': true,
    };
    for (const [k, want] of Object.entries(expected)) {
      expect(R2.__call__(BigInt(k)).is_square()).toBe(want);
    }
  });

  test('all square roots', () => {
    const R = Zp(3n, 20);
    const roots = R.__call__(4n).square_root_all();
    expect(roots.length).toBe(2);
    expect(roots[0]!.pow(2n).eq(R.__call__(4n))).toBe(true);
    expect(roots[1]!.pow(2n).eq(R.__call__(4n))).toBe(true);
    expect(roots[0]!.eq(roots[1]!)).toBe(false);
  });
});

describe('p-adic nth roots', () => {
  const R = Zp(5n, 10);

  test('cube roots actually cube back', () => {
    for (const v of [8n, 27n]) {
      const r = R.__call__(v).nth_root(3n);
      expect(r.pow(3n).eq(R.__call__(v))).toBe(true);
      expect(r.lift() ** 3n % 5n ** 10n).toBe(v);
    }
  });

  test('SageMath doctest: A = Zp(5,10); A(61376).nth_root(4)', () => {
    const x = R.__call__(61376n);
    expect(x.toString()).toBe('1 + 5^3 + 3*5^4 + 4*5^5 + 3*5^6 + O(5^10)');
    const y = x.nth_root(4n);
    expect(y.pow(4n).eq(x)).toBe(true);
    expect(y.precision_absolute()).toBe(10);

    const all = x.nth_root_all(4n);
    expect(all.length).toBe(4);
    for (const r of all) {
      expect(r.pow(4n).eq(x)).toBe(true);
    }
  });

  test('nth roots when p divides n lose one digit per p-th root', () => {
    const u = R.__call__(7n);
    const u5 = u.pow(5n);
    const r = u5.nth_root(5n);
    expect(r.precision_relative()).toBe(9);
    expect(r.pow(5n).eq(u5.add_bigoh(r.precision_absolute()))).toBe(true);
    expect(() => R.__call__(6n).nth_root(5n)).toThrow(ValueError);
  });

  test('nth root with nonzero valuation', () => {
    const v = R.__call__(8n * 125n); // 8 * 5^3
    const r = v.nth_root(3n);
    expect(r.valuation()).toBe(1n);
    expect(r.pow(3n).eq(v)).toBe(true);
    expect(() => R.__call__(8n * 25n).nth_root(3n)).toThrow(ValueError);
  });

  test('nth root of a non-power throws', () => {
    // Cubing is a bijection on (Z/5)^*, so every 5-adic unit is a cube; but the
    // 4th powers are exactly the units congruent to 1 mod 5.
    expect(R.__call__(2n).is_nth_power(3n)).toBe(true);
    expect(() => R.__call__(2n).nth_root(4n)).toThrow(ValueError);
    expect(R.__call__(2n).is_nth_power(4n)).toBe(false);
    expect(R.__call__(16n).is_nth_power(4n)).toBe(true);
    expect(R.__call__(16n).nth_root(4n).pow(4n).eq(R.__call__(16n))).toBe(true);
  });

  test('n = 0 and exact zero', () => {
    expect(() => R.__call__(2n).nth_root(0n)).toThrow(ValueError);
    expect(R.zero().nth_root(3n).is_zero()).toBe(true);
  });
});

describe('p-adic exp and log values', () => {
  const R = Zp(5n, 10);

  test('exp(5) has the exact expected lift', () => {
    // sum_{n>=0} 5^n/n! mod 5^10 computed with exact rational arithmetic
    expect(R.__call__(5n).exp().lift()).toBe(3474831n);
  });

  test('log(1+5) has the exact expected lift', () => {
    // sum_{n>=1} (-1)^(n+1) 5^n/n mod 5^10
    expect(R.__call__(6n).log().lift()).toBe(6970555n);
  });

  test('exp and log are inverse (SageMath: Z13(14).log().exp() == 14)', () => {
    expect(R.__call__(6n).log().exp().eq(R.__call__(6n))).toBe(true);
    expect(R.__call__(5n).exp().log().eq(R.__call__(5n))).toBe(true);

    const S = Zp(13n, 10);
    expect(S.__call__(14n).log().exp().eq(S.__call__(14n))).toBe(true);
  });

  test('exp does not converge for valuation 0', () => {
    expect(() => Zp(2n, 5).__call__(2n).exp()).toThrow(ValueError);
  });
});

describe('p-adic uniformizer powers', () => {
  test('negative powers work in a field', () => {
    const K = Qp(5n, 20);
    expect(K.uniformizer_pow(-2n).toString()).toBe('5^-2 + O(5^18)');
    expect(K.uniformizer_pow(3n).lift()).toBe(125n);
  });

  test('negative powers throw in a ring', () => {
    expect(() => Zp(5n, 20).uniformizer_pow(-2n)).toThrow(ValueError);
  });

  test('roots of unity', () => {
    const K = Qp(5n, 20);
    const mu4 = K.roots_of_unity(4n);
    expect(mu4.length).toBe(4);
    for (const z of mu4) {
      expect(z.pow(4n).is_one()).toBe(true);
    }
    // SageMath: Zp(5,10).roots_of_unity() is [1, w(2), w(4), w(3)] -- the powers
    // of the Teichmuller lift of a primitive root mod p.
    const R = Zp(5n, 10);
    expect(R.roots_of_unity().map((e) => e.residue(1))).toEqual([1n, 2n, 4n, 3n]);
    expect(R.roots_of_unity(10n).map((e) => e.residue(1))).toEqual([1n, 4n]);
    expect(K.roots_of_unity(3n).length).toBe(1);
    expect(Qp(2n, 20).roots_of_unity(2n).length).toBe(2);
  });
});
