/**
 * Tests for `sage/rings/function_field/` (rational function fields).
 *
 * Every expected value below was produced by running SageMath 10.3 itself
 * (`/usr/local/bin/sage`) on the corresponding expression, or is copied
 * verbatim from an upstream doctest (cited where that is the case).  A
 * 34 800-line transcript covering GF(2), GF(3), GF(5), GF(7), GF(11),
 * GF(65537), GF(131101) and QQ was diffed against SageMath while developing
 * this module; the cases below are the distilled regression set.
 */
import { describe, expect, test } from 'bun:test';
import { NotImplementedError, ValueError } from '../../errors.js';
import { GF } from '../finite_rings/finite_field_constructor.js';
import { QQ } from '../rational_field.js';
import type { ConstantField, ConstantFieldElement } from './constant_field.js';
import { FunctionField } from './constructor.js';
import type { FunctionFieldElement_rational } from './element_rational.js';
import { is_FunctionField } from './function_field.js';
import type {
  RationalFunctionField,
  RationalFunctionField_global,
} from './function_field_rational.js';
import { RationalFunctionField_char_zero } from './function_field_rational.js';

type CE = ConstantFieldElement;
type FF = FunctionFieldElement_rational<CE>;

function gf(q: bigint): RationalFunctionField_global<CE> {
  return FunctionField(
    GF(q) as unknown as ConstantField<CE>,
    'x'
  ) as RationalFunctionField_global<CE>;
}

function qq(name = 'x'): RationalFunctionField<CE> {
  return FunctionField(QQ as unknown as ConstantField<CE>, name) as RationalFunctionField<CE>;
}

/** Divisor repr with the newlines upstream inserts flattened. */
function d(x: { toString(): string }): string {
  return x.toString().replace(/\n/g, ' ~ ');
}

describe('RationalFunctionField', () => {
  test('repr, genus, characteristic (function_field_rational.py:54, :743)', () => {
    const K = gf(3n);
    // sage: K.<t> = FunctionField(GF(3)); K
    // Rational function field in t over Finite Field of size 3
    expect(gf(3n).toString()).toBe('Rational function field in x over Finite Field of size 3');
    expect(K.genus()).toBe(0n);
    expect(K.characteristic()).toBe(3n);
    expect(K.degree()).toBe(1n);
    expect(K.ngens()).toBe(1);
    expect(K.is_global()).toBe(true);
    expect(K.is_finite()).toBe(false);
    expect(K.is_perfect()).toBe(false);
    expect(K.base_field()).toBe(K);
    expect(K.rational_function_field()).toBe(K);
    expect(is_FunctionField(K)).toBe(true);
    expect(is_FunctionField(GF(3n))).toBe(false);
  });

  test('QQ constant field is char zero and not global', () => {
    const K = qq();
    expect(K.toString()).toBe('Rational function field in x over Rational Field');
    expect(K.characteristic()).toBe(0n);
    expect(K.is_global()).toBe(false);
    expect(K.is_perfect()).toBe(true);
    expect(K instanceof RationalFunctionField_char_zero).toBe(true);
  });

  test('constant_field must be a field (function_field_rational.py:146)', () => {
    // sage: FunctionField(ZZ, 't')
    // TypeError: constant_field must be a field
    const notAField = {
      zero: () => 0,
      one: () => 1,
      __call__: (x: unknown) => x,
      is_field: () => false,
      toString: () => 'Integer Ring',
    };
    expect(() => FunctionField(notAField as unknown as ConstantField<CE>, 't')).toThrow(
      'constant_field must be a field'
    );
  });

  test('gen(1) raises (function_field_rational.py:567)', () => {
    // sage: K.gen(1)
    // IndexError: Only one generator.
    expect(() => gf(5n).gen(1)).toThrow('Only one generator.');
  });

  test('orders (order.py:15, order_rational.py:41)', () => {
    const K = gf(19n);
    expect(K.maximal_order().toString()).toBe(
      'Maximal order of Rational function field in x over Finite Field of size 19'
    );
    expect(K.maximal_order_infinite().toString()).toBe(
      'Maximal infinite order of Rational function field in x over Finite Field of size 19'
    );
    expect(K.equation_order()).toBe(K.maximal_order());
    expect(K.maximal_order().basis().map(String)).toEqual(['1']);
    expect(K.maximal_order().ngens()).toBe(1);
    expect(K.maximal_order().gen().toString()).toBe('x');
    // sage: 1/x in O -> False ; 1/x in Oinf -> True   (order.py:21)
    const x = K.gen();
    expect(K.maximal_order().contains(x.inv())).toBe(false);
    expect(K.maximal_order_infinite().contains(x.inv())).toBe(true);
    expect(K.maximal_order_infinite().prime_ideal().toString()).toBe(
      'Ideal (1/x) of Maximal infinite order of Rational function field in x over Finite Field of size 19'
    );
  });

  test('divisor group / place set reprs', () => {
    const K = gf(5n);
    expect(K.divisor_group().toString()).toBe(
      'Divisor group of Rational function field in x over Finite Field of size 5'
    );
    expect(K.place_set().toString()).toBe(
      'Set of places of Rational function field in x over Finite Field of size 5'
    );
    expect(d(K.different())).toBe('0');
  });
});

describe('elements', () => {
  test('repr follows Frac(k[x]) (fraction_field_element.pyx:523)', () => {
    const K = gf(5n);
    const x = K.gen();
    const one = K.one();
    const c = (n: number) => K.__call__(n);
    // Values produced by SageMath: K.<x> = FunctionField(GF(5))
    expect(x.toString()).toBe('x');
    expect(x.mul(c(2)).toString()).toBe('2*x');
    expect(one.div(x).toString()).toBe('1/x');
    expect(one.div(x.mul(c(2))).toString()).toBe('3/x'); //  1/(2x) = 3/x over GF(5)
    expect(one.div(x.pow(2)).toString()).toBe('1/x^2');
    expect(x.pow(2).add(c(2)).div(x.pow(2)).toString()).toBe('(x^2 + 2)/x^2');
    expect(x.mul(c(2)).div(x.pow(2).add(one)).toString()).toBe('2*x/(x^2 + 1)');
    expect(c(3).div(x.add(one)).toString()).toBe('3/(x + 1)');
    expect(x.add(one).div(x.sub(one)).toString()).toBe('(x + 1)/(x + 4)');
    expect(one.div(x.pow(3).add(x.pow(2)).add(x)).toString()).toBe('1/(x^3 + x^2 + x)');
    expect(x.div(x.pow(2).add(x)).toString()).toBe('1/(x + 1)');
    expect(K.zero().toString()).toBe('0');
  });

  test('upstream doctests verbatim', () => {
    // sage: K.<t> = FunctionField(GF(3)); 1/t + t^3 + 5
    // (t^4 + 2*t + 1)/t                     (function_field_rational.py:58)
    const K3 = FunctionField(
      GF(3n) as unknown as ConstantField<CE>,
      't'
    ) as RationalFunctionField_global<CE>;
    const t = K3.gen();
    expect(t.inv().add(t.pow(3)).add(K3.__call__(5)).toString()).toBe('(t^4 + 2*t + 1)/t');
    // sage: K.<t> = FunctionField(QQ); 1/t + t^3 + 5
    // (t^4 + 5*t + 1)/t                     (function_field_rational.py:66)
    const KQ = qq('t');
    const u = KQ.gen();
    expect(u.inv().add(u.pow(3)).add(KQ.__call__(5)).toString()).toBe('(t^4 + 5*t + 1)/t');
    // sage: K.<x> = FunctionField(GF(5^2,'a')); f = (x^2 + x + 1)/(x^3 + 1); f^3
    // (x^6 + 3*x^5 + x^4 + 2*x^3 + x^2 + 3*x + 1)/(x^9 + 3*x^6 + 3*x^3 + 1)
    //                                       (function_field.py:16)
    // All coefficients lie in the prime field, so GF(5) reproduces it exactly.
    const K5 = gf(5n);
    const x = K5.gen();
    const one = K5.one();
    const f = x.pow(2).add(x).add(one).div(x.pow(3).add(one));
    expect(f.toString()).toBe('(x^2 + x + 1)/(x^3 + 1)');
    expect(f.pow(3).toString()).toBe(
      '(x^6 + 3*x^5 + x^4 + 2*x^3 + x^2 + 3*x + 1)/(x^9 + 3*x^6 + 3*x^3 + 1)'
    );
  });

  test('degree, numerator, denominator (element.pyx:561)', () => {
    // sage: FF.<t> = FunctionField(QQ); f = (t^2 + 3)/(t^3 - 1/3); f.degree() -> 3
    const QQK = qq('t');
    const u = QQK.gen();
    expect(
      u
        .pow(2)
        .add(QQK.__call__(3))
        .div(u.pow(3).sub(QQK.__call__(1).div(QQK.__call__(3))))
        .degree()
    ).toBe(3n);
    // Over GF(7) the same expression reduces: x^2 + 3 and x^3 - 1 share x - 2,
    // so SageMath reports (x + 2)/(x^2 + 2*x + 4) with degree 2.
    const K = gf(7n);
    const t = K.gen();
    const one = K.one();
    const f = t.pow(2).add(K.__call__(3)).div(t.pow(3).sub(one));
    expect(f.toString()).toBe('(x + 2)/(x^2 + 2*x + 4)');
    expect(f.degree()).toBe(2n);
    expect(t.add(K.__call__(1)).degree()).toBe(1n);
    expect(K.zero().degree()).toBe(0n);
    const g = t.add(one).div(t.pow(2).add(K.__call__(3)));
    expect(g.numerator().toString()).toBe('x + 1');
    expect(g.denominator().toString()).toBe('x^2 + 3');
  });

  test('arithmetic normalises like Sage', () => {
    const K = gf(11n);
    const x = K.gen();
    const one = K.one();
    // sage: (t+1)/(t^2-1) -> 1/(t - 1)
    expect(x.add(one).div(x.pow(2).sub(one)).toString()).toBe('1/(x + 10)');
    expect(x.pow(-2).toString()).toBe('1/x^2');
    expect(x.add(one).pow(3).toString()).toBe('x^3 + 3*x^2 + 3*x + 1');
    expect(x.mul(x.inv()).is_one()).toBe(true);
  });

  test('valuation at an irreducible polynomial (element_rational.pyx:278)', () => {
    // sage: K.<t> = FunctionField(QQ)
    // sage: f = (t - 1)^2*(t + 1)/(t^2 - 1/3)^3
    // sage: f.valuation(t - 1) -> 2 ; f.valuation(t) -> 0 ; f.valuation(t^2-1/3) -> -3
    const K = gf(7n);
    const t = K.gen();
    const one = K.one();
    const f = t
      .sub(one)
      .pow(2)
      .mul(t.add(one))
      .div(t.pow(2).add(K.__call__(3)).pow(3));
    expect(f.valuation(t.sub(one))).toBe(2n);
    expect(f.valuation(t)).toBe(0n);
    expect(f.valuation(t.pow(2).add(K.__call__(3)))).toBe(-3n);
  });

  test('valuation at places (element_rational.pyx:300)', () => {
    // sage: K.<x> = FunctionField(GF(2)); p = K.places_finite()[0]
    // sage: (1/x^2).valuation(p) -> -2
    const K = gf(2n);
    const p = K.places_finite()[0]!;
    expect(K.gen().pow(-2).valuation(p)).toBe(-2n);
  });

  test('inverse_mod (element_rational.pyx:494)', () => {
    // sage: K.<x> = FunctionField(QQ)
    // sage: O = K.maximal_order(); I = O.ideal(x^2 + 1)
    // sage: t = O(x + 1).inverse_mod(I); t -> -1/2*x + 1/2
    const K = qq();
    const x = K.gen();
    const O = K.maximal_order();
    const I = O.ideal(x.pow(2).add(K.one()));
    const t = O.__call__(x.add(K.one())).inverse_mod(I);
    expect(t.numerator().getCoeff(1).toString()).toBe('-1/2');
    expect(t.numerator().getCoeff(0).toString()).toBe('1/2');
    expect(I.contains(t.mul(x.add(K.one())).sub(K.one()) as FF)).toBe(true);
  });

  test('inverse_mod without an inverse raises Sage message', () => {
    // sage: K.<x> = FunctionField(GF(2)); O = K.maximal_order()
    // sage: O(x+1).inverse_mod(O.ideal(x^2+1))
    // ValueError: Impossible inverse modulo
    const K = gf(2n);
    const x = K.gen();
    const O = K.maximal_order();
    expect(() => O.__call__(x.add(K.one())).inverse_mod(O.ideal(x.pow(2).add(K.one())))).toThrow(
      'Impossible inverse modulo'
    );
  });

  test('is_square / sqrt (element_rational.pyx:316)', () => {
    // sage: K.<t> = FunctionField(GF(5)); (-t^2).is_square() -> True; (-t^2).sqrt() -> 2*t
    const K = gf(5n);
    const t = K.gen();
    const f = t.pow(2).neg();
    expect(f.is_square()).toBe(true);
    expect((f.sqrt() as FF).pow(2).eq(f)).toBe(true);
    expect(t.is_square()).toBe(false);
  });

  test('factor sorts by (degree, exponent, prime) like Factorization.sort', () => {
    // sage: K.<x> = FunctionField(GF(2)); (x^2/(x+1)).factor()
    // (x + 1)^-1 * x^2      [factors in this order]
    const K = gf(2n);
    const x = K.gen();
    const one = K.one();
    const f = x.pow(2).div(x.add(one));
    const fa = f.factor();
    expect(fa.factors.map(([p, e]) => `${p} ^ ${e}`)).toEqual(['x + 1 ^ -1', 'x ^ 2']);
    expect(fa.unit.toString()).toBe('1');
    // sage: ((x^3+x+1)/(x^2+1)).factor() over GF(2)
    const g = x.pow(3).add(x).add(one).div(x.pow(2).add(one));
    expect(g.factor().factors.map(([p, e]) => `${p} ^ ${e}`)).toEqual([
      'x + 1 ^ -2',
      'x^3 + x + 1 ^ 1',
    ]);
  });

  test('trace and norm over the field itself (element.pyx:525, :538)', () => {
    const K = gf(5n);
    const x = K.gen();
    expect(x.trace().toString()).toBe('x');
    expect(x.norm().toString()).toBe('x');
    expect(x.matrix().map((r) => r.map(String))).toEqual([['x']]);
  });

  test('divisor of zero raises (element.pyx:756)', () => {
    const K = gf(5n);
    expect(() => K.zero().divisor()).toThrow('divisor not defined for zero');
    expect(() => K.zero().divisor_of_zeros()).toThrow('divisor of zeros not defined for zero');
    expect(() => K.zero().divisor_of_poles()).toThrow('divisor of poles not defined for zero');
  });

  test('nth_root/is_nth_power stubs name what is missing', () => {
    const K = gf(3n);
    const x = K.gen();
    expect(x.is_nth_power(1)).toBe(true);
    expect(() => x.is_nth_power(3)).toThrow(NotImplementedError);
    expect(() => x.nth_root(5)).toThrow(NotImplementedError);
  });
});

describe('places', () => {
  test('places of GF(5)(x) (function_field_rational.py:863)', () => {
    // sage: F.<x> = FunctionField(GF(5)); F.places()
    // [Place (1/x), Place (x), Place (x + 1), Place (x + 2), Place (x + 3), Place (x + 4)]
    expect(gf(5n).places().map(String)).toEqual([
      'Place (1/x)',
      'Place (x)',
      'Place (x + 1)',
      'Place (x + 2)',
      'Place (x + 3)',
      'Place (x + 4)',
    ]);
    // sage: F.places_finite()
    expect(gf(5n).places_finite().map(String)).toEqual([
      'Place (x)',
      'Place (x + 1)',
      'Place (x + 2)',
      'Place (x + 3)',
      'Place (x + 4)',
    ]);
  });

  test('places of GF(2)(x) (place_rational.py:48)', () => {
    // sage: F.<x> = FunctionField(GF(2)); F.places()
    // [Place (1/x), Place (x), Place (x + 1)]
    // sage: [p.is_infinite_place() for p in F.places()] -> [True, False, False]
    // sage: [p.local_uniformizer() for p in F.places()] -> [1/x, x, x + 1]
    const F = gf(2n);
    expect(F.places().map(String)).toEqual(['Place (1/x)', 'Place (x)', 'Place (x + 1)']);
    expect(F.places().map((p) => p.is_infinite_place())).toEqual([true, false, false]);
    expect(F.places().map((p) => p.local_uniformizer().toString())).toEqual(['1/x', 'x', 'x + 1']);
    expect(F.places().map((p) => p.degree())).toEqual([1n, 1n, 1n]);
  });

  test('degree-2 places of GF(5)(x) keep Sage enumeration order', () => {
    // sage: F.<x> = FunctionField(GF(5)); F.places_finite(2)
    expect(gf(5n).places_finite(2).map(String)).toEqual([
      'Place (x^2 + 2)',
      'Place (x^2 + 3)',
      'Place (x^2 + x + 1)',
      'Place (x^2 + x + 2)',
      'Place (x^2 + 2*x + 3)',
      'Place (x^2 + 2*x + 4)',
      'Place (x^2 + 3*x + 3)',
      'Place (x^2 + 3*x + 4)',
      'Place (x^2 + 4*x + 1)',
      'Place (x^2 + 4*x + 2)',
    ]);
  });

  test('get_place over GF(2) (function_field_rational.py:938)', () => {
    // sage: F.<a> = GF(2); K.<x> = FunctionField(F)
    // sage: K.get_place(1) -> Place (x)
    // sage: K.get_place(2) -> Place (x^2 + x + 1)
    // sage: K.get_place(3) -> Place (x^3 + x + 1)
    // sage: K.get_place(4) -> Place (x^4 + x + 1)
    // sage: K.get_place(5) -> Place (x^5 + x^2 + 1)
    const K = gf(2n);
    expect([1, 2, 3, 4, 5].map((n) => K.get_place(n).toString())).toEqual([
      'Place (x)',
      'Place (x^2 + x + 1)',
      'Place (x^3 + x + 1)',
      'Place (x^4 + x + 1)',
      'Place (x^5 + x^2 + 1)',
    ]);
  });

  test('degree of a degree-2 place (place_rational.py:30)', () => {
    // sage: F.<x> = FunctionField(GF(2)); O = F.maximal_order()
    // sage: O.ideal(x^2 + x + 1).place().degree() -> 2
    const F = gf(2n);
    const x = F.gen();
    const p = F.maximal_order().ideal(x.pow(2).add(x).add(F.one())).place();
    expect(p.degree()).toBe(2n);
    expect(p.toString()).toBe('Place (x^2 + x + 1)');
    expect(p.valuation_ring().toString()).toBe('Valuation ring at Place (x^2 + x + 1)');
  });

  test('places sort with infinity first', () => {
    const K = gf(3n);
    const shuffled = [...K.places()].reverse();
    shuffled.sort((a, b) => a.cmp(b));
    expect(shuffled.map(String)).toEqual(K.places().map(String));
  });

  test('non-prime ideals have no place (ideal.py:313)', () => {
    // sage: K.<x> = FunctionField(GF(4)); O.ideal(x^2 + x + 1).place()
    // TypeError: not a prime ideal        (x^2+x+1 splits over GF(4))
    const K = gf(2n);
    const x = K.gen();
    const O = K.maximal_order();
    expect(() => O.ideal(x.mul(x.add(K.one()))).place()).toThrow('not a prime ideal');
  });

  test('residue field of a degree-1 place (place_rational.py:105)', () => {
    const F = gf(5n);
    const x = F.gen();
    // finite place (x + 3): residue field is GF(5), and x |-> -3 = 2
    const p = F.maximal_order()
      .ideal(x.add(F.__call__(3)))
      .place();
    const [R, from_R, to_R] = p._residue_field();
    expect(R.toString()).toBe('Finite Field of size 5');
    expect(to_R(x).toString()).toBe('2');
    expect(from_R(R.__call__(4)).toString()).toBe('4');
    // infinite place: x has a pole, 1/x |-> 0, (x+1)/x |-> 1
    const pinf = F.place_infinite();
    const [Rinf, , to_inf] = pinf._residue_field();
    expect(Rinf.toString()).toBe('Finite Field of size 5');
    expect(to_inf(x.inv()).toString()).toBe('0');
    expect(to_inf(x.add(F.one()).div(x)).toString()).toBe('1');
  });

  test('residue fields of higher-degree places are an honest stub', () => {
    const F = gf(2n);
    const p = F.get_place(3);
    expect(() => p._residue_field()).toThrow(NotImplementedError);
    expect(() => p._residue_field()).toThrow(/_residue_field for places of degree 3 > 1/);
  });

  test('evaluate (element.pyx:909)', () => {
    // sage: K.<t> = FunctionField(GF(5)); p = K.place_infinite()
    // sage: f = 1/t^2 + 3; f.evaluate(p) -> 3
    const K = gf(5n);
    const t = K.gen();
    const f = t.pow(-2).add(K.__call__(3));
    expect(f.evaluate(K.place_infinite()).toString()).toBe('3');
    // a pole gives ValueError
    expect(() => t.evaluate(K.place_infinite())).toThrow('has a pole at the place');
  });
});

describe('ideals', () => {
  test('reprs and arithmetic (ideal.py:13, order_rational.py:407)', () => {
    // sage: K.<x> = FunctionField(QQ); O = K.maximal_order()
    // sage: I = O.ideal(x^3 + 1); I    -> Ideal (x^3 + 1) of Maximal order ...
    // sage: I^2                        -> Ideal (x^6 + 2*x^3 + 1) ...
    // sage: ~I                         -> Ideal (1/(x^3 + 1)) ...
    // sage: ~I * I                     -> Ideal (1) ...
    const K = qq();
    const x = K.gen();
    const O = K.maximal_order();
    const I = O.ideal(x.pow(3).add(K.one()));
    const suffix = ' of Maximal order of Rational function field in x over Rational Field';
    expect(I.toString()).toBe(`Ideal (x^3 + 1)${suffix}`);
    expect(I.pow(2).toString()).toBe(`Ideal (x^6 + 2*x^3 + 1)${suffix}`);
    expect(I.inv().toString()).toBe(`Ideal (1/(x^3 + 1))${suffix}`);
    expect(I.inv().mul(I).toString()).toBe(`Ideal (1)${suffix}`);
    expect(I.div(I).toString()).toBe(`Ideal (1)${suffix}`);
  });

  test('generated ideal is the gcd/lcm normal form (order_rational.py:407)', () => {
    // sage: K.<x> = FunctionField(QQ); O = K.maximal_order()
    // sage: O.ideal(x^3 + 1, x^3 + 6)     -> Ideal (1)
    // sage: O.ideal((x^2+1)*(x^3+1), (x^3+6)*(x^2+1)) -> Ideal (x^2 + 1)
    const K = qq();
    const x = K.gen();
    const one = K.one();
    const O = K.maximal_order();
    const c = (n: number) => K.__call__(n);
    expect(
      O.ideal(x.pow(3).add(one), x.pow(3).add(c(6)))
        .gen()
        .toString()
    ).toBe('1');
    const a = x.pow(2).add(one).mul(x.pow(3).add(one));
    const b = x.pow(3).add(c(6)).mul(x.pow(2).add(one));
    expect(O.ideal(a, b).gen().toString()).toBe('x^2 + 1');
    expect(O.ideal([x, x.inv()]).eq(O.ideal(x, x.inv()))).toBe(true);
    // an ideal argument re-uses its generators
    expect(O.ideal(O.ideal(a, b)).gen().toString()).toBe('x^2 + 1');
  });

  test('infinite ideals are powers of x (order_rational.py:542)', () => {
    // sage: K.<x> = FunctionField(QQ); O = K.maximal_order_infinite()
    // sage: O.ideal(x^3 + 1, x^3 + 6)                  -> Ideal (x^3)
    // sage: O.ideal((x^2+1)*(x^3+1), (x^3+6)*(x^2+1))  -> Ideal (x^5)
    const K = qq();
    const x = K.gen();
    const one = K.one();
    const c = (n: number) => K.__call__(n);
    const O = K.maximal_order_infinite();
    expect(
      O.ideal(x.pow(3).add(one), x.pow(3).add(c(6)))
        .gen()
        .toString()
    ).toBe('x^3');
    const a = x.pow(2).add(one).mul(x.pow(3).add(one));
    const b = x.pow(3).add(c(6)).mul(x.pow(2).add(one));
    expect(O.ideal(a, b).gen().toString()).toBe('x^5');
  });

  test('infinite ideal reductions over GF(2) (ideal_rational.py:545)', () => {
    // sage: K.<x> = FunctionField(GF(2)); Oinf = K.maximal_order_infinite()
    // sage: Oinf.ideal((x+1)/(x^3+x), (x^2+1)/x^4).gen()  -> 1/x^2
    // sage: Oinf.ideal(x/(x^2+1)) + Oinf.ideal(1/(x+1))   -> Ideal (1/x)
    // sage: Oinf.ideal(x/(x^2+1)) * Oinf.ideal(1/(x+1))   -> Ideal (1/x^2)
    // sage: ~Oinf.ideal(x/(x^2 + 1))                      -> Ideal (x)
    const K = gf(2n);
    const x = K.gen();
    const one = K.one();
    const Oinf = K.maximal_order_infinite();
    const g1 = x.add(one).div(x.pow(3).add(x));
    const g2 = x.pow(2).add(one).div(x.pow(4));
    expect(Oinf.ideal(g1, g2).gen().toString()).toBe('1/x^2');
    const I = Oinf.ideal(x.div(x.pow(2).add(one)));
    const J = Oinf.ideal(one.div(x.add(one)));
    expect(I.add(J).gens()[0]!.toString()).toBe('1/x');
    expect(I.mul(J).gens()[0]!.toString()).toBe('1/x^2');
    expect(I.inv().gens()[0]!.toString()).toBe('x');
    expect(I.is_prime()).toBe(true);
  });

  test('is_prime and factorisation (ideal_rational.py:201, :333)', () => {
    // sage: K.<x> = FunctionField(GF(2)); O = K.maximal_order()
    // sage: I = O.ideal(x^3*(x+1)^2); I.factor()
    // (Ideal (x))^3 * (Ideal (x + 1))^2
    const K = gf(2n);
    const x = K.gen();
    const one = K.one();
    const O = K.maximal_order();
    const I = O.ideal(x.pow(3).mul(x.add(one).pow(2)));
    expect(I.factor().map(([p, m]) => `${p.gens()[0]} ^ ${m}`)).toEqual(['x ^ 3', 'x + 1 ^ 2']);
    expect(I.factor().every(([p]) => p.is_prime())).toBe(true);
    expect(O.ideal(x).is_prime()).toBe(true);
    expect(O.ideal(x.pow(2)).is_prime()).toBe(false);
    expect(O.ideal(x.inv()).is_prime()).toBe(false);
  });

  test('valuation of ideals (ideal_rational.py:286, :582)', () => {
    // sage: F.<x> = FunctionField(QQ); O = F.maximal_order()
    // sage: I = O.ideal(x^2*(x^2+x+1)^3); [f.valuation(I) for f,_ in I.factor()] -> [2, 3]
    const F = qq();
    const x = F.gen();
    const one = F.one();
    const O = F.maximal_order();
    const I = O.ideal(x.pow(2).mul(x.pow(2).add(x).add(one).pow(3)));
    expect(I.factor().map(([p]) => p.valuation(I))).toEqual([2n, 3n]);
    // sage: p = O.ideal(x); p.valuation(O.ideal(x + 1)) -> 0, (x^2) -> 2, (1/x^3) -> -3, (0) -> +Infinity
    const p = O.ideal(x);
    expect(p.valuation(O.ideal(x.add(one)))).toBe(0n);
    expect(p.valuation(O.ideal(x.pow(2)))).toBe(2n);
    expect(p.valuation(O.ideal(x.pow(-3)))).toBe(-3n);
    expect(p.valuation(O.ideal(F.zero()))).toBe(Number.POSITIVE_INFINITY);
    // sage: Oinf = F.maximal_order_infinite(); q = Oinf.ideal(1/x)
    // sage: q.valuation(Oinf.ideal(x/(x+1))) -> 0 ; q.valuation(Oinf.ideal(0)) -> +Infinity
    const Oinf = F.maximal_order_infinite();
    const q = Oinf.ideal(x.inv());
    expect(q.valuation(Oinf.ideal(x.div(x.add(one))))).toBe(0n);
    expect(q.valuation(Oinf.ideal(F.zero()))).toBe(Number.POSITIVE_INFINITY);
    expect(() => O.ideal(x.pow(2)).valuation(I)).toThrow('not a prime ideal');
  });

  test('the zero ideal keeps SageMath quirks', () => {
    // Verified with SageMath 10.3: `O.ideal(K(0))` prints as `Ideal (0) of ...`
    // and `is_zero()` is False (Element.__bool__ has no zero to compare to).
    const K = gf(5n);
    const O = K.maximal_order();
    const I = O.ideal(K.zero());
    expect(I.is_zero()).toBe(false);
    expect(I.toString()).toBe(
      'Ideal (0) of Maximal order of Rational function field in x over Finite Field of size 5'
    );
    // sage: I.divisor() -> 0   (empty factorization over odd prime fields)
    expect(d(I.divisor())).toBe('0');
    // sage: Oinf.ideal(K(0)).divisor() -> Place (1/x)
    expect(d(K.maximal_order_infinite().ideal(K.zero()).divisor())).toBe('Place (1/x)');
  });

  test('ideal monoid (ideal.py:1017)', () => {
    const K = gf(2n);
    const O = K.maximal_order();
    const M = O.ideal_monoid();
    expect(M.toString()).toBe(
      'Monoid of ideals of Maximal order of Rational function field in x over Finite Field of size 2'
    );
    expect(M.ring()).toBe(O);
    expect(M.__call__(K.gen()).toString()).toBe(
      'Ideal (x) of Maximal order of Rational function field in x over Finite Field of size 2'
    );
  });
});

describe('divisors', () => {
  test('divisor of an element (element.pyx:756)', () => {
    // sage: K.<x> = FunctionField(GF(2)); f = 1/(x^3 + x^2 + x); f.divisor()
    // 3*Place (1/x) - Place (x) - Place (x^2 + x + 1)
    const K = gf(2n);
    const x = K.gen();
    const f = K.one().div(x.pow(3).add(x.pow(2)).add(x));
    expect(d(f.divisor())).toBe('3*Place (1/x) ~  - Place (x) ~  - Place (x^2 + x + 1)');
    expect(d(f.divisor_of_zeros())).toBe('3*Place (1/x)');
    expect(d(f.divisor_of_poles())).toBe('Place (x) ~  + Place (x^2 + x + 1)');
    expect(f.zeros().map(String)).toEqual(['Place (1/x)']);
    expect(f.poles().map(String)).toEqual(['Place (x)', 'Place (x^2 + x + 1)']);
    expect(f.divisor().degree()).toBe(0n);
    // sage: x.divisor() -> - Place (1/x) + Place (x)          (divisor.py:198)
    expect(d(x.divisor())).toBe('- Place (1/x) ~  + Place (x)');
  });

  test('divisor arithmetic and formatting (divisor.py:182)', () => {
    const K = gf(5n);
    const O = K.maximal_order();
    const x = K.gen();
    const p1 = O.ideal(x).place();
    const p2 = O.ideal(x.add(K.one())).place();
    const pinf = K.place_infinite();
    const D = p1.divisor(2).sub(p2.divisor(3));
    expect(d(D)).toBe('2*Place (x) ~  - 3*Place (x + 1)');
    expect(D.degree()).toBe(-1n);
    expect(D.multiplicity(p1)).toBe(2n);
    expect(D.multiplicity(p2)).toBe(-3n);
    expect(D.multiplicity(pinf)).toBe(0n);
    expect(D.is_effective()).toBe(false);
    expect(p1.divisor(2).add(p2.divisor(3)).is_effective()).toBe(true);
    expect(d(D.numerator())).toBe('2*Place (x)');
    expect(d(D.denominator())).toBe('3*Place (x + 1)');
    expect(D.numerator().sub(D.denominator()).eq(D)).toBe(true);
    expect(D.neg().neg().eq(D)).toBe(true);
    expect(D.scalar_mul(3).sub(D).eq(D.scalar_mul(2))).toBe(true);
    expect(D.add(D.scalar_mul(2)).eq(D.scalar_mul(3))).toBe(true);
    expect(D.support().map(String)).toEqual(['Place (x)', 'Place (x + 1)']);
    expect(d(K.divisor_group().zero())).toBe('0');
    // upstream keeps zero multiplicities produced by `place.divisor(0)`
    expect(d(pinf.divisor(0))).toBe('0*Place (1/x)');
    expect(pinf.divisor(0).support().map(String)).toEqual(['Place (1/x)']);
  });

  test('Riemann-Roch on P^1 over GF(5) (divisor.py:554)', () => {
    // All values from SageMath: K.<x> = FunctionField(GF(5))
    const K = gf(5n);
    const O = K.maximal_order();
    const x = K.gen();
    const pinf = K.place_infinite();
    const P0 = O.ideal(x).place();
    const Pd = K.get_place(2); // Place (x^2 + 2)

    const D1 = pinf.divisor(3);
    expect(D1.degree()).toBe(3n);
    expect(D1.dimension()).toBe(4n);
    expect(D1.basis_function_space().map(String)).toEqual(['x^3', 'x^2', 'x', '1']);

    const D2 = P0.divisor(2);
    expect(D2.dimension()).toBe(3n);
    expect(D2.basis_function_space().map(String)).toEqual(['1', '1/x', '1/x^2']);

    const D3 = pinf.divisor(-1).add(P0.divisor());
    expect(D3.degree()).toBe(0n);
    expect(D3.dimension()).toBe(1n);
    expect(D3.basis_function_space().map(String)).toEqual(['1/x']);

    const D4 = P0.divisor(2).add(pinf.divisor(3)).sub(Pd.divisor());
    expect(d(D4)).toBe('3*Place (1/x) ~  + 2*Place (x) ~  - Place (x^2 + 2)');
    expect(D4.degree()).toBe(3n);
    expect(D4.dimension()).toBe(4n);
    expect(D4.basis_function_space().map(String)).toEqual([
      'x^3 + 2*x',
      'x^2 + 2',
      '(x^2 + 2)/x',
      '(x^2 + 2)/x^2',
    ]);

    // negative degree divisors have empty Riemann-Roch spaces
    expect(pinf.divisor(-1).dimension()).toBe(0n);
    expect(P0.divisor(-3).dimension()).toBe(0n);
    // deg D = -1 still gives l(D) = 0 = deg D + 1 - g on P^1
    expect(K.divisor_group().zero().dimension()).toBe(1n);
  });

  test('Riemann-Roch theorem l(D) - i(D) = deg D + 1 - g on P^1', () => {
    for (const q of [2n, 3n, 5n, 7n]) {
      const K = gf(q);
      const O = K.maximal_order();
      const x = K.gen();
      const pinf = K.place_infinite();
      const p1 = O.ideal(x).place();
      const pd = K.get_place(2);
      const cands = [
        K.divisor_group().zero(),
        pinf.divisor(3),
        pinf.divisor(-3),
        p1.divisor(4),
        p1.divisor(-2),
        pd.divisor(2),
        pd.divisor(-1),
        p1.divisor(2).sub(pinf.divisor(3)),
        pinf.divisor(3).add(p1.divisor(2)).sub(pd.divisor()),
        pd.divisor(4).sub(p1.divisor()),
      ];
      for (const D of cands) {
        const l = D.dimension();
        const i = BigInt(D._differential_space()[0].length);
        expect(l - i).toBe(D.degree() - K.genus() + 1n);
      }
    }
  });

  test('Riemann-Roch coordinates round-trip (divisor.py:632)', () => {
    const K = gf(7n);
    const O = K.maximal_order();
    const x = K.gen();
    const D = K.place_infinite().divisor(2).add(O.ideal(x).place().divisor(1));
    const [basis, coordinates] = D._function_space();
    expect(basis.length).toBe(Number(D.dimension()));
    for (let i = 0; i < basis.length; i++) {
      const coords = coordinates(basis[i]!);
      expect(coords.map((c) => c.toString())).toEqual(basis.map((_, j) => (j === i ? '1' : '0')));
    }
    const v = basis[0]!.add(basis[1]!) as FF;
    expect(coordinates(v).map(String).slice(0, 2)).toEqual(['1', '1']);
    const [n, from_V, to_V] = D.function_space();
    expect(n).toBe(basis.length);
    const k = K.constant_base_field();
    const vec = basis.map((_, j) => k.__call__(j + 1));
    expect(to_V(from_V(vec)).map(String)).toEqual(vec.map(String));
  });

  test('divisors compare and the group rejects nonzero construction', () => {
    const K = gf(3n);
    const O = K.maximal_order();
    const x = K.gen();
    const A = O.ideal(x).place().divisor(1);
    const B = K.place_infinite().divisor(1);
    expect(A.cmp(B)).toBe(1); // finite places sort after the infinite place
    expect(B.cmp(A)).toBe(-1);
    expect(A.cmp(A)).toBe(0);
    expect(K.divisor_group().__call__(0).eq(K.divisor_group().zero())).toBe(true);
    expect(() => K.divisor_group().__call__(1)).toThrow(ValueError);
  });

  test('differential space basis is an honest stub', () => {
    const K = gf(5n);
    expect(() => K.divisor_group().zero().basis_differential_space()).toThrow(NotImplementedError);
    expect(() => K.space_of_differentials()).toThrow(/differential\.py/);
  });
});

describe('not-yet-ported entry points throw', () => {
  test('extensions and companions name their upstream file', () => {
    const K = gf(2n);
    expect(() => K.extension(null)).toThrow(/function_field_polymod\.py/);
    expect(() => K.completion(null)).toThrow(NotImplementedError);
    expect(() => K.jacobian()).toThrow(NotImplementedError);
    expect(() => K.valuation(null)).toThrow(/valuation\.py/);
    expect(() => K.higher_derivation()).toThrow(NotImplementedError);
    expect(() => K.polynomial_ring()).toThrow(NotImplementedError);
    expect(() => K.gen().differential()).toThrow(/differential\.py/);
  });
});
