import { expect, test } from 'bun:test';
import { GF } from '../finite_rings/index.js';
import { QQ } from '../rational_field.js';
import { AssertionError, IndexError, NotImplementedError, ZeroDivisionError } from '../../errors.js';
import { FunctionField } from './constructor.js';
import { FunctionFieldElement_rational } from './element_rational.js';
import { IdealMonoid } from './ideal.js';
import type { ConstantField, ConstantFieldElement } from './constant_field.js';

const field = (p: bigint) => FunctionField((p ? GF(p) : QQ) as unknown as ConstantField<ConstantFieldElement>, 'x');

test('rational function parents reuse factory keys and return both renaming maps', () => {
  const K = field(0n), k = K.constant_base_field();
  expect(FunctionField(k, ['x'])).toBe(K);
  expect(K.degree(FunctionField(k, 'x'))).toBe(1n);
  const [L, fromL, toL] = K.change_variable_name('y');
  expect(L.change_variable_name('x')[0]).toBe(K);
  expect(toL(K.gen()).toString()).toBe('y');
  expect(fromL(L.gen()).toString()).toBe('x');
  const f = K.gen().add(K.one()).div(K.gen().sub(K.one()));
  expect(fromL(toL(f)).eq(f)).toBe(true);
  const [same, from, to] = K.change_variable_name('x');
  expect(same).toBe(K);
  expect(from).toBe(to);
  expect(from(f)).toBe(f);
});

test('function field and order generators retain Sage IndexError', () => {
  const K = field(7n);
  expect(() => K.gen(1)).toThrow(new IndexError('Only one generator.'));
  expect(() => K.maximal_order().gen(-1)).toThrow(new IndexError('there is only one generator'));
});

test('fraction comparison selects FpT ordering and compares rational coefficients numerically', () => {
  const Q = field(0n);
  expect(Q.one().cmp(Q.one().div(Q.__call__(2n)))).toBe(1);
  for (const p of [5n, 46337n]) {
    const K = field(p), x = K.gen();
    const f = K.__call__(-2n).div(K.__call__(3n));
    const g = x.add(K.one()).pow(2).div(x.sub(K.one()).pow(2));
    expect(f.cmp(g)).toBe(-1);
  }
});

test('square roots and zero errors retain the selected fraction backend behavior', () => {
  const K = field(2n);
  expect((K.gen().pow(2).sqrt(true) as unknown[]).length).toBe(1);
  expect(() => K.gen().nth_root(2)).toThrow('element is not an n-th power');
  expect(() => field(5n).gen().sqrt()).toThrow(new NotImplementedError('function fields not yet implemented'));
  expect(() => field(0n).gen().sqrt()).toThrow(TypeError);
  for (const p of [0n,2n,5n,46349n]) {
    const F = field(p);
    expect(() => F.zero().inv()).toThrow(new ZeroDivisionError(p===5n ? '' : 'fraction field element division by zero'));
    expect(F.zero().is_nth_power(-1)).toBe(false);
  }
});

test('valuation rejects constants promptly and validates zero inputs at the divisor', () => {
  for (const p of [0n,2n,5n,65537n]) {
    const K = field(p);
    for (const f of [K.zero(),K.one(),K.gen()]) {
      expect(() => f.valuation(K.one())).toThrow('The polynomial, p, must be non-constant.');
      expect(() => f.valuation(K.zero())).toThrow(ZeroDivisionError);
    }
    expect(K.zero().valuation(K.gen())).toBe(Number.POSITIVE_INFINITY);
    expect(K.gen().pow(3).valuation(K.gen())).toBe(3n);
  }
});

test('valuation converts the complete fraction argument through the polynomial ring', () => {
  for (const p of [0n,2n,5n,65537n]) {
    const K = field(p);
    for (const f of [K.zero(),K.gen()]) {
      expect(() => f.valuation(K.gen().inv())).toThrow(p===5n?'not integral':'fraction must have unit denominator');
    }
  }
});

test('zero factorization reaches the bundled polynomial guard, including finite ideals', () => {
  for (const p of [0n,2n,4n,5n,9n,65537n]) {
    const K = field(p);
    expect(() => K.zero().factor()).toThrow('factorization of 0 is not defined');
    expect(() => K.maximal_order().ideal(K.zero()).factor()).toThrow('factorization of 0 is not defined');
    expect(() => K.maximal_order().ideal(K.zero()).divisor()).toThrow('factorization of 0 is not defined');
    expect(K.maximal_order_infinite().ideal(K.zero()).divisor().toString()).toBe('Place (1/x)');
  }
});


test('function fields expose and coerce their actual fraction representation', () => {
  const K=field(0n),F=K.field(),f=K.__call__('1/x');
  expect(F).toBe(K._ring.fraction_field());
  expect(f.element().parent).toBe(F);
  expect(f.element()).toBe(f.element());
  expect(K.__call__(f.element()).element()).toBe(f.element());
  expect(F.__call__('x==x','1/x').toString()).toBe('x');
  expect(field(2n).field().__call__('1/x','1/2').toString()).toBe('0');
  expect(K.__call__('(x*x-1)/(x-1)').toString()).toBe('x + 1');
  expect(() => K.__call__('y')).toThrow("unable to evaluate 'y' in "+F);
});

test('fraction string fallback compiles before looking up unknown names', () => {
  const F=field(0n).field();
  expect(() => F.__call__("y'")).toThrow('unterminated string literal');
  expect(() => F.__call__('y\nx')).toThrow('invalid syntax (<string>, line 2)');
  expect(() => F.__call__('y\tx')).toThrow('invalid syntax (<string>, line 1)');
  expect(() => F.__call__('y\\x')).toThrow('unexpected character after line continuation character');
});

test('wrapped fractions preserve native normalization, predicates and power identity', () => {
  const K=field(5n),F=K.field(),R=K._ring;
  const raw=new F._element_class(F,R.zero(),R.gen(),{reduce:false});
  const f=new FunctionFieldElement_rational(K,raw);
  expect(f.is_zero()).toBe(false);
  expect(f.toString()).toBe('0/x');
  expect(f.pow(1n)).toBe(f);
  expect(f.neg().toString()).toBe('0');
  const unit=new FunctionFieldElement_rational(K,new F._element_class(F,R.gen(),R.gen(),{reduce:false}));
  expect(unit.is_one()).toBe(false);
  R.__call__(unit.element());
  expect(unit.numerator().toString()).toBe('1');
  expect(unit.denominator().toString()).toBe('1');
  expect(unit.is_one()).toBe(true);
});


test('modular inversion keeps native assertions and unit gcd scaling at zero modulus', () => {
  for(const p of [0n,3n,5n,46349n]){
    const K=field(p),O=K.maximal_order(),two=K.__call__(2n);
    expect(two.inverse_mod(O.ideal(K.zero())).eq(two.inv())).toBe(true);
    expect(() => K.gen().inv().inverse_mod(O.ideal(K.gen()))).toThrow(new AssertionError());
    expect(() => K.one().inverse_mod(O.ideal(K.gen().inv()))).toThrow(new AssertionError());
  }
});

test('orders propagate syntax and division errors from element construction', () => {
  const K=field(5n);
  for(const O of [K.maximal_order(),K.maximal_order_infinite()]){
    expect(() => O.__call__('1/0')).toThrow(ZeroDivisionError);
    expect(() => O.__call__('x+')).toThrow(SyntaxError);
    expect(() => O.__call__('y')).toThrow('unable to convert to an element of '+K);
    expect(() => O.is_subring(O)).toThrow(new NotImplementedError(''));
  }
});

test('orders convert ideals across finite and infinite parents and retain monoid identity', () => {
  const K=field(5n),O=K.maximal_order(),P=K.maximal_order_infinite(),I=O.ideal(K.gen());
  expect(P.ideal(I).gens()[0]!.toString()).toBe('x');
  const M=O.ideal_monoid();
  expect(M).toBe(O.ideal_monoid());
  expect(M).toBe(new IdealMonoid(O));
  expect(M.__call__(I)).toBe(I);
  expect(M.__call__({gens:()=>[K.gen()]}).eq(I)).toBe(true);
});

test('the bundled infinite order basis is the integer singleton tuple', () => {
  expect(field(5n).maximal_order_infinite().basis()).toEqual([1n]);
});
