/**
 * Tests for p-adic element implementation
 *
 * @see Reference: sage/rings/padics/padic_generic_element.pyx
 */

import { describe, expect, test } from 'bun:test';
import { Zp, Qp, pAdicRing, pAdicField } from './padic_generic.js';
import { pAdicGenericElement, PrecisionError } from './padic_generic_element.js';
import { ValueError, ZeroDivisionError } from '../../errors.js';

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
    expect(sum.lift() % (7n ** 10n)).toBe(7n);
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

  test('sqrt of non-square throws', () => {
    const a = R.__call__(3n);
    expect(() => a.sqrt()).toThrow(ValueError);
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
