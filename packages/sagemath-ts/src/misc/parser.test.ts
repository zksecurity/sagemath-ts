import { describe, expect, test } from 'bun:test';
import { GFpn, PrimeField } from '../rings/finite_rings/finite_field_extension.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import { LookupNameMaker, Parser, Tokenizer } from './parser.js';

describe('Sage expression grammar', () => {
  const parser = () =>
    new Parser<string>(
      (s) => 'I' + s,
      (s) => 'F' + s,
      (s) => 'V' + s,
      {},
      true,
      {
        binary: (op, a, b) => '(' + op + ' ' + a + ' ' + b + ')',
        unary: (op, a) => '(' + op + ' ' + a + ')',
      }
    );
  test('exponentiation is right associative and precedes unary negation', () => {
    expect(parser().parse('-2^3^2')).toBe('(- (^ I2 (^ I3 I2)))');
    expect(parser().parse('2^-3')).toBe('(^ I2 (- I3))');
    expect(parser().parse('2a+a*b')).toBe('(+ (* I2 Va) (* Va Vb))');
  });
  test('implicit multiplication applies only before names', () => {
    expect(() => parser().parse('2(a+1)')).toThrow('Malformed expression');
    expect(() => parser().parse('a(a+1)')).toThrow("Unknown variable: 'a'");
    expect(() => parser().parse('(a')).toThrow('Mismatched parentheses');
  });
  test('Unicode digits and whitespace follow Python', () => {
    expect(new Tokenizer('²+½+⑨').test()).toEqual(['INT(²)', '+', 'NAME(½)', '+', 'INT(⑨)']);
    expect(new Tokenizer('a\u0085+\u001cb').test()).toEqual(['NAME(a)', '+', 'NAME(b)']);
    expect(new Tokenizer('a\ufeffb').test()).toEqual(['NAME(a)', 'ERROR', 'NAME(b)']);
  });
  test('lookahead, reset and one-step backtracking preserve token state', () => {
    const t = new Tokenizer('12+a');
    expect(t.peek()).toBe(128);
    expect(t.next()).toBe(128);
    expect(t.last_token_string()).toBe('12');
    expect(t.backtrack()).toBe(false);
    expect(() => t.backtrack()).toThrow('Can only backtrack once.');
    t.reset(1);
    expect(t.test()).toEqual(['INT(2)', '+', 'NAME(a)']);
    expect(t.test()).toEqual([]);
  });
  test('lookup ignores inherited JavaScript properties', () => {
    const names = new LookupNameMaker<string>({ a: 'old' }, (x) => 'fallback:' + x);
    expect(names.__call__('constructor')).toBe('fallback:constructor');
    names.set_names({ ['__proto__']: 'value' });
    expect(names.__call__('__proto__')).toBe('value');
  });
});

describe('finite-field polynomial string construction', () => {
  test('fractions cancel before coercion back into the polynomial ring', () => {
    const R = new PolynomialRing(new PrimeField(7n), 'a');
    expect(R.__call__('(a^2-1)/(a-1)').toString()).toBe('a + 1');
    expect(R.__call__('1/a*a').toString()).toBe('1');
    expect(() => R.__call__('1/a')).toThrow('not integral');
    expect(() => R.__call__('a+1/2')).toThrow("and 'Rational Field'");
  });
  test('string evaluation takes place in the polynomial ring before extension reduction', () => {
    const F = GFpn(3n, 2);
    expect(F.__call__('a^2+1').toString()).toBe('a + 2');
    expect(() => F.__call__('a^-1')).toThrow('not integral');
    expect(F.gen().pow(-1n).mul(F.gen()).toString()).toBe('1');
  });
});
