import { expect, test } from 'bun:test';
import { ValueError, ZeroDivisionError } from '../errors.js';
import { matrix_gf2_from_entries, zero_matrix_gf2 } from '../matrix/matrix_mod2.js';
import { current_randstate, set_random_seed } from '../misc/randstate.js';
import { Mod } from '../rings/finite_rings/integer_mod.js';
import { discrete_log_lambda, discrete_log_rho } from './generic.js';

test('singleton lambda bounds hash once and preserve the untouched random stream', () => {
  set_random_seed(0);
  const hashes: bigint[] = [];
  expect(() =>
    discrete_log_lambda(
      Mod(13n, 509n),
      Mod(1n, 509n),
      [13n, 13n],
      '+',
      undefined,
      undefined,
      undefined,
      (x) => {
        hashes.push(x.value);
        return x.value;
      }
    )
  ).toThrow(new ZeroDivisionError('integer modulo by zero'));
  expect(hashes).toEqual([13n]);
  expect(current_randstate().python_random().getrandbits(64)).toBe(15049713646458849547n);
});

test('singleton lambda invokes the hash before the zero-modulus error', () => {
  expect(() =>
    discrete_log_lambda(
      Mod(1n, 11n),
      Mod(1n, 11n),
      [1n, 1n],
      '+',
      undefined,
      undefined,
      undefined,
      () => {
        throw new TypeError('hash sentinel');
      }
    )
  ).toThrow('hash sentinel');
});

test('lambda freezes actual binary matrices before hashing and freezes the wild target', () => {
  set_random_seed(0);
  const base = matrix_gf2_from_entries([
    [0, 1],
    [1, 1],
  ]);
  const identity = zero_matrix_gf2(2, 2);
  const target = identity.add(base);
  const result = discrete_log_lambda(
    target,
    base,
    [0n, 16n],
    'other',
    identity,
    (x) => x.neg(),
    (x, y) => x.add(y),
    (x) => {
      expect(x.is_immutable()).toBe(true);
      let value = 0n;
      for (let i = 0; i < 2; i++)
        for (let j = 0; j < 2; j++) value = 2n * value + BigInt(x.get(i, j));
      return value;
    }
  );
  expect(result).toBe(15n);
  expect(base.is_immutable()).toBe(false);
  expect(target.is_immutable()).toBe(true);
});

test('singleton lambda freezes a reused binary-matrix base before failing', () => {
  const base = matrix_gf2_from_entries([
    [0, 1],
    [1, 1],
  ]);
  const identity = zero_matrix_gf2(2, 2);
  expect(() =>
    discrete_log_lambda(
      identity,
      base,
      [1n, 1n],
      'other',
      identity,
      (x) => x.neg(),
      (x, y) => x.add(y),
      () => 0n
    )
  ).toThrow('integer modulo by zero');
  expect(base.is_immutable()).toBe(true);
  expect(identity.is_immutable()).toBe(false);
});

function mutableGroup(
  collisions = false,
  freezeThrows = false,
  parentThrows = false,
  modulus = 509n
) {
  let freezes = 0;
  const parent = {
    zero: () => {
      if (parentThrows) throw new ValueError('identity sentinel');
      return new Element(0n);
    },
  };
  const norm = (v: bigint) => ((v % modulus) + modulus) % modulus;
  class Element {
    frozen = false;
    parent = parent;
    constructor(readonly v: bigint) {}
    add(other: Element) {
      return new Element(norm(this.v + other.v));
    }
    mul(n: bigint) {
      return new Element(norm(this.v * n));
    }
    neg() {
      return new Element(norm(-this.v));
    }
    eq(other: Element) {
      return this.v === other.v;
    }
    set_immutable() {
      freezes++;
      this.frozen = true;
      if (freezeThrows) throw new ValueError('freeze sentinel');
    }
    toString() {
      if (!this.frozen) throw new TypeError('mutable element');
      return String(collisions ? this.v % 4n : this.v);
    }
  }
  return { element: (n: bigint) => new Element(norm(n)), freezes: () => freezes };
}

test('rho freezes walk points before hashing without freezing the input target', () => {
  set_random_seed(0);
  const group = mutableGroup();
  const base = group.element(1n),
    target = group.element(17n);
  expect(discrete_log_rho(target, base, 509n, '+')).toBe(17n);
  expect(group.freezes()).toBeGreaterThan(0);
  expect(base.frozen).toBe(false);
  expect(target.frozen).toBe(false);
});

test('rho distinguishes unequal stored points with colliding string keys', () => {
  set_random_seed(0);
  const group = mutableGroup(true);
  expect(discrete_log_rho(group.element(17n), group.element(1n), 509n, '+')).toBe(17n);
  expect(current_randstate().python_random().getrandbits(64)).toBe(16564833271613235632n);
});

test('rho propagates mutation errors after native walk initialization', () => {
  set_random_seed(0);
  const group = mutableGroup(false, true);
  expect(() => discrete_log_rho(group.element(17n), group.element(1n), 509n, '+')).toThrow(
    'freeze sentinel'
  );
  expect(group.freezes()).toBe(1);
  expect(current_randstate().python_random().getrandbits(64)).toBe(16564833271613235632n);
});

test('rho parses the parent before prime-order validation', () => {
  const group = mutableGroup(false, false, true);
  expect(() => discrete_log_rho(group.element(17n), group.element(1n), 510n, '+')).toThrow(
    'identity sentinel'
  );
});

test('rho small-order BSGS fallback does not introduce mutation calls', () => {
  const group = mutableGroup(false, false, false, 11n);
  expect(discrete_log_rho(group.element(6n), group.element(1n), 11n, '+')).toBe(6n);
  expect(group.freezes()).toBe(0);
});
