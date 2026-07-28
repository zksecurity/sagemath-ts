/**
 * Tests for rational_field.ts (the field QQ).
 *
 * Port of the doctests in reference/sage/src/sage/rings/rational_field.py.
 */

import { describe, expect, test } from 'bun:test';
import { legendre_symbol } from '../arith/misc.js';
import { TypeError, ValueError } from '../errors.js';
import { Integer } from './integer_ring.js';
import { Rational } from './rational.js';
import { QQ } from './rational_field.js';

const R = (n: bigint, d: bigint = 1n) => new Rational(n, d);

describe('RationalField basics', () => {
  test('properties', () => {
    expect(QQ.is_field()).toBe(true);
    expect(QQ.characteristic()).toBe(0n);
    expect(QQ.is_prime_field()).toBe(true);
    expect(QQ.is_finite()).toBe(false);
    expect(QQ.degree()).toBe(1n);
    expect(QQ.discriminant()).toBe(1n);
    expect(QQ.class_number()).toBe(1n);
    expect(QQ.signature()).toEqual([1n, 0n]);
  });

  test('coercion', () => {
    expect(QQ.__call__('3/4').toString()).toBe('3/4');
    expect(QQ.__call__(7n).toString()).toBe('7');
    expect(QQ.__call__(3n, 4n).toString()).toBe('3/4');
    expect(QQ.zero().toString()).toBe('0');
    expect(QQ.one().toString()).toBe('1');
  });
});

/**
 * `QQ(x)` is `Rational(x)` -- `_element_constructor_ = Rational`
 * (reference/sage/src/sage/rings/rational_field.py:243) -- so it must accept
 * every shape `Rational.__set_value`
 * (reference/sage/src/sage/rings/rational.pyx:590-704) accepts.
 */
describe('QQ.__call__ accepts every form Sage accepts', () => {
  test('Sage doctest values', () => {
    // rational_field.py:226 -- QQ._element_constructor_((2, 3)) == 2/3
    expect(QQ.__call__([2n, 3n]).toString()).toBe('2/3');
    // rational.pyx:548-552 (issue 19835) -- QQ((0r,-1r)) == 0, QQ((-1r,-1r)) == 1
    expect(QQ.__call__([0n, -1n]).toString()).toBe('0');
    expect(QQ.__call__([-1n, -1n]).toString()).toBe('1');
    // rational.pyx:531-533 -- a.__init__(7); a == 7
    expect(QQ.__call__(7n).toString()).toBe('7');
    // rational.pyx:673 -- a length-1 list is unwrapped
    expect(QQ.__call__([5n]).toString()).toBe('5');
  });

  test('integers, Integer wrappers, floats and strings', () => {
    expect(QQ.__call__(new Integer(5n)).toString()).toBe('5');
    expect(QQ.__call__(new Integer(-12n)).toString()).toBe('-12');
    expect(QQ.__call__(5).toString()).toBe('5');
    expect(QQ.__call__(1.5).toString()).toBe('3/2');
    expect(QQ.__call__('-17/34').toString()).toBe('-1/2');
    expect(QQ.__call__('0.125').toString()).toBe('1/8');
  });

  test('pairs are reduced and sign-normalized', () => {
    expect(QQ.__call__(6n, 8n).toString()).toBe('3/4');
    expect(QQ.__call__(6n, -8n).toString()).toBe('-3/4');
    expect(QQ.__call__([new Integer(6n), new Integer(-8n)]).toString()).toBe('-3/4');
    // a Rational entry is read as a quotient
    expect(QQ.__call__(R(1n, 2n), 3n).toString()).toBe('1/6');
  });

  test('objects: numerator/denominator, numer/denom, _rational_', () => {
    // rational.pyx:692 -- fractions.Fraction
    expect(QQ.__call__({ numerator: 6n, denominator: -8n }).toString()).toBe('-3/4');
    // the property names Rational itself exposes (used by the polynomial code)
    expect(QQ.__call__({ numer: 6n, denom: 8n }).toString()).toBe('3/4');
    // rational.pyx:637 -- hasattr(x, "_rational_")
    expect(QQ.__call__({ _rational_: () => R(7n, 9n) }).toString()).toBe('7/9');
    // a Rational is returned unchanged
    const r = R(2n, 3n);
    expect(QQ.__call__(r)).toBe(r);
  });

  test('errors match Sage', () => {
    // rational.pyx:704 -- TypeError("unable to convert {!r} to a rational")
    expect(() => QQ.__call__(null as never)).toThrow(TypeError);
    expect(() => QQ.__call__(null as never)).toThrow('unable to convert None to a rational');
    expect(() => QQ.__call__({ foo: 1 } as never)).toThrow('to a rational');
    // rational.pyx:426-429 -- QQ('1/0') is a TypeError, not a division error
    expect(() => QQ.__call__('1/0')).toThrow("unable to convert '1/0' to a rational");
    expect(() => QQ.__call__('not a number')).toThrow('to a rational');
    // rational.pyx:663 -- ValueError("denominator must not be 0")
    expect(() => QQ.__call__(3n, 0n)).toThrow(ValueError);
    expect(() => QQ.__call__([3n, 0n])).toThrow('denominator must not be 0');
  });

  test('round trip through QQ preserves arithmetic', () => {
    for (let num = -12n; num <= 12n; num++) {
      for (let den = 1n; den <= 7n; den++) {
        const direct = R(num, den);
        expect(QQ.__call__(num, den).eq(direct)).toBe(true);
        expect(QQ.__call__([num, den]).eq(direct)).toBe(true);
        expect(QQ.__call__(direct.toString()).eq(direct)).toBe(true);
        expect(QQ.__call__({ numer: num, denom: den }).eq(direct)).toBe(true);
      }
    }
  });
});

/**
 * `quadratic_defect` -- reference/sage/src/sage/rings/rational_field.py:1457.
 *
 * Sage's algorithm:
 *   v, u = self(a).val_unit(p)
 *   if v % 2 == 1: return v                       # Python's non-negative %
 *   if p != 2: return Infinity if legendre_symbol(u, p) == 1 else v
 *   # p == 2, u % 8 computed with Rational.__mod__
 *   1 -> Infinity, 5 -> v+2, 3 or 7 -> v+1
 *
 * `legendre_symbol(x, p)` evaluates the Kronecker symbol of
 * `numerator(x) * denominator(x)` (sage/arith/misc.py:4365), so the *rational*
 * unit must be passed, not just its numerator.
 */
describe('QQ.quadratic_defect', () => {
  test('Sage doctests', () => {
    expect(QQ.quadratic_defect(0n, 7n)).toBe('Infinity');
    expect(QQ.quadratic_defect(5n, 7n)).toBe(0n);
    expect(QQ.quadratic_defect(5n, 2n)).toBe(2n);
    expect(QQ.quadratic_defect(5n, 5n)).toBe(1n);
  });

  test('non-integral input (audit counterexamples)', () => {
    // v = 0, u = 1/3, legendre(1*3, 7) = legendre(3,7) = -1 ... but
    // legendre(1/3, 7) = legendre(3, 7)? kronecker(1*3,7) = -1 -> defect 0.
    expect(QQ.quadratic_defect(R(1n, 3n), 7n)).toBe(0n);
    // p = 2: u % 8 = 1 * inverse(3, 8) = 3 -> v + 1 = 1
    expect(QQ.quadratic_defect(R(1n, 3n), 2n)).toBe(1n);
    // negative odd valuation: v(1/5) at 5 is -1, and -1 % 2 == 1 in Python
    expect(QQ.quadratic_defect(R(1n, 5n), 5n)).toBe(-1n);
  });

  test('negative odd valuations are returned verbatim', () => {
    expect(QQ.quadratic_defect(R(1n, 7n), 7n)).toBe(-1n);
    expect(QQ.quadratic_defect(R(1n, 8n), 2n)).toBe(-3n);
    expect(QQ.quadratic_defect(R(3n, 2n), 2n)).toBe(-1n);
  });

  test('rejects non-primes when check is set', () => {
    expect(() => QQ.quadratic_defect(5n, 9n)).toThrow('9 must be prime');
  });

  /**
   * Independent structural oracle: the quadratic defect of `a` is `+Infinity`
   * exactly when `a` is a square in Q_p.
   */
  test('defect is Infinity iff a is a square in Q_p', () => {
    const isSquareInQp = (num: bigint, den: bigint, p: bigint): boolean => {
      if (num === 0n) return true;
      // a and a*den^2 = num*den lie in the same square class
      let x = num * den;
      let v = 0n;
      while (x % p === 0n) {
        x /= p;
        v++;
      }
      if (v % 2n !== 0n) return false;
      if (p === 2n) return ((x % 8n) + 8n) % 8n === 1n;
      return legendre_symbol(x, p) === 1n;
    };

    for (const p of [2n, 3n, 5n, 7n, 11n, 13n]) {
      for (let num = -40n; num <= 40n; num++) {
        for (let den = 1n; den <= 12n; den++) {
          const r = R(num, den);
          const isInf = QQ.quadratic_defect(r, p, false) === 'Infinity';
          expect(isInf).toBe(isSquareInQp(r.numerator, r.denominator, p));
        }
      }
    }
  });

  /**
   * Second oracle: transcription of Sage's branch structure, computed
   * independently of the implementation under test.
   */
  test('agrees with a direct transcription of the Sage algorithm', () => {
    const oracle = (num: bigint, den: bigint, p: bigint): bigint | 'Infinity' => {
      if (num === 0n) return 'Infinity';
      let v = 0n;
      let n = num;
      let d = den;
      while (n % p === 0n) {
        n /= p;
        v++;
      }
      while (d % p === 0n) {
        d /= p;
        v--;
      }
      if (((v % 2n) + 2n) % 2n === 1n) return v;
      if (p !== 2n) return legendre_symbol(n * d, p) === 1n ? 'Infinity' : v;
      const nm = ((n % 8n) + 8n) % 8n;
      const dm = ((d % 8n) + 8n) % 8n;
      let dinv = 1n;
      for (let k = 1n; k < 8n; k += 2n) if ((dm * k) % 8n === 1n) dinv = k;
      const u = (nm * dinv) % 8n;
      if (u === 1n) return 'Infinity';
      if (u === 5n) return v + 2n;
      return v + 1n;
    };

    for (const p of [2n, 3n, 5n, 7n, 11n, 13n]) {
      for (let num = -60n; num <= 60n; num++) {
        for (let den = 1n; den <= 20n; den++) {
          const r = R(num, den);
          expect(QQ.quadratic_defect(r, p, false)).toBe(oracle(r.numerator, r.denominator, p));
        }
      }
    }
  });
});
