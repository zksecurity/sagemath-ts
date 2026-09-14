import { nativeFixtures as loadLiveNative } from "../../../../../tests/property/native-live.mjs";


import { expect, test } from 'bun:test';
import { execFileSync } from 'node:child_process';
import { GF } from '../../rings/finite_rings/finite_field_constructor.js';
import { GFpn, PrimeField } from '../../rings/finite_rings/finite_field_extension.js';
import { GF2 } from '../../rings/finite_rings/gf2.js';
import { QQ } from '../../rings/rational_field.js';
import { Polynomial as IsomorphismTracePolynomial } from '../../rings/polynomial/polynomial_element.js';
import { EllipticCurve } from './constructor.js';
import { _isomorphisms, baseWI, WeierstrassIsomorphism } from './weierstrass_morphism.js';
const native = (await loadLiveNative(import.meta.url, "./ell_generic_isomorphism_order.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];

// Runtime adapters cover the port's separate QQ and finite-field interfaces.
type Any = any;
const functions: Record<string, (...args: Any[]) => string> = {};
const field = (p: bigint): Any => (p === 0n ? QQ : GF(p));
const curve = (p: bigint, c: bigint[]): Any => EllipticCurve(field(p), c as Any);
functions.ec_generic_isomorphism_order = (
  p: bigint,
  degree: bigint,
  modulus: bigint[],
  left: bigint[],
  right: bigint[],
  transform: bigint[]
) => {
  try {
    const K: Any =
      degree > 1n ? GFpn(p, Number(degree), modulus.slice(0, -1).map(Number), 'a') : field(p);
    const decode = (v: bigint) => (degree > 1n ? K.fromInteger(v) : K.__call__(v));
    const encode = (v: Any) => (degree > 1n ? String(v.integer_representation()) : String(v));
    const E = EllipticCurve(K, left.map(decode) as Any);
    const F = EllipticCurve(
      K,
      transform.length
        ? (new baseWI(...(transform.map(decode) as [Any, Any, Any, Any])).call(
            E.a_invariants()
          ) as Any)
        : (right.map(decode) as Any)
    );
    let first;
    try {
      first = { value: E.isomorphism_to(F).map(encode) };
    } catch (e) {
      first = { error: (e as Error).name, message: (e as Error).message };
    }
    const ordered = E.isomorphisms(F).map((t) => [
      t.map(encode),
      t.every((v) => p === 0n || v.parent === K),
    ]);
    const same = E.a_invariants().every((v, i) => v.eq(F.a_invariants()[i]!));
    const autos = same ? E.automorphisms().map((t) => t.map(encode)) : null;
    return JSON.stringify({ value: [first, ordered, autos, E.is_isomorphic(F)] });
  } catch (e) {
    return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
  }
};
functions.wm_isomorphism_comparisons = (
  p: bigint,
  degree: bigint,
  modulus: bigint[],
  left: bigint[],
  transform: bigint[]
) => {
  try {
    const K: Any =
      degree > 1n ? GFpn(p, Number(degree), modulus.slice(0, -1).map(Number), 'a') : field(p);
    const decode = (v: bigint) => (degree > 1n ? K.fromInteger(v) : K.__call__(v));
    const E = EllipticCurve(K, left.map(decode) as Any);
    const F = transform.length
      ? EllipticCurve(
          K,
          new baseWI(...(transform.map(decode) as [Any, Any, Any, Any])).call(
            E.a_invariants()
          ) as Any
        )
      : E;
    const morphisms = [..._isomorphisms(E, F)].map((t) => new WeierstrassIsomorphism(E, t, F));
    const compare = (a: Any, b: Any) =>
      ['lt', 'le', 'eq', 'ne', 'gt', 'ge'].map((op) => {
        try {
          return WeierstrassIsomorphism._comparison_impl(a, b, op);
        } catch (e) {
          return { error: (e as Error).name, message: (e as Error).message };
        }
      });
    const matrix = morphisms.map((a) => morphisms.map((b) => compare(a, b)));
    const identityE = new WeierstrassIsomorphism(E, [K.one(), K.zero(), K.zero(), K.zero()], E);
    const identityF = new WeierstrassIsomorphism(F, [K.one(), K.zero(), K.zero(), K.zero()], F);
    const domains = compare(identityE, identityF),
      codomains = morphisms.length ? compare(identityE, morphisms[0]) : null;
    const invalid = [compare(null, identityE), compare(identityE, null), compare(null, 7)];
    return JSON.stringify({ value: [matrix, domains, codomains, invalid] });
  } catch (e) {
    return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
  }
};
functions.ec_is_isomorphic_arguments = (kind: bigint) => {
  try {
    const E = curve(5n, [0n, 1n]),
      other = [null, 7, [], {}, E][Number(kind)];
    return JSON.stringify({ value: E.is_isomorphic(other) });
  } catch (e) {
    return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
  }
};
functions.ec_isomorphism_parent_guards = (left: bigint, right: bigint) => {
  try {
    const field = (kind: number): Any => {
      if (kind === 0) return QQ;
      if (kind >= 16) {
        const [p, g]: [bigint, bigint] = (
          {
            16: [7n, 2n],
            17: [7n, 3n],
            18: [7n, 2n],
            19: [7n, 1n],
            20: [2n, 0n],
            21: [2n, 1n],
          } as Any
        )[kind];
        const R = new IsomorphismParentPolynomialRing(new PrimeField(p), 'x');
        return new PrimeField(p, { modulus: R.__call__([-g, 1n]) });
      }
      if (kind === 12) return new PrimeField(5n);
      if (kind === 13) return GF2;
      const primes: Record<number, bigint> = { 1: 2n, 2: 3n, 3: 5n, 4: 7n };
      if (primes[kind]) return GF(primes[kind]);
      const [p, d, m, name]: Any = (
        {
          5: [2n, 2, [1, 1], 'a'],
          6: [2n, 2, [1, 1], 'b'],
          7: [3n, 2, [1, 0], 'a'],
          8: [3n, 2, [2, 1], 'a'],
          9: [3n, 2, [1, 0], 'b'],
          10: [5n, 2, [2, 0], 'a'],
          11: [3n, 2, [1, 0], 'a'],
          14: [2n, 3, [1, 1, 0], 'a'],
          15: [3n, 3, [1, 2, 0], 'a'],
        } as Any
      )[kind];
      return GFpn(p, d, m, name);
    };
    const curve = (K: Any) =>
      EllipticCurve(
        K,
        K.characteristic === 2n
          ? [0n, 0n, 1n, 0n, 0n]
          : K.characteristic === 3n
            ? [0n, 0n, 0n, 1n, 0n]
            : [0n, 0n, 0n, 0n, 1n]
      );
    const E = curve(field(Number(left))),
      F = curve(field(Number(right)));
    const identity = (E: Any) =>
      new WeierstrassIsomorphism(
        E,
        [E.base_ring.one(), E.base_ring.zero(), E.base_ring.zero(), E.base_ring.zero()],
        E
      );
    const a = identity(E),
      b = identity(F);
    const comparison = ['lt', 'le', 'eq', 'ne', 'gt', 'ge'].map((op) => {
      try {
        return WeierstrassIsomorphism._comparison_impl(a, b, op);
      } catch (e) {
        return { error: (e as Error).name, message: (e as Error).message };
      }
    });
    let isomorphic;
    try {
      isomorphic = { value: E.is_isomorphic(F) };
    } catch (e) {
      isomorphic = { error: (e as Error).name, message: (e as Error).message };
    }
    return JSON.stringify({ value: [isomorphic, comparison] });
  } catch (e) {
    return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
  }
};
functions.ec_isomorphism_root_trace = (
  p: bigint,
  degree: bigint,
  modulus: bigint[],
  coefficients: bigint[],
  operation: bigint
) => {
  try {
    const K: Any =
      degree > 1n ? GFpn(p, Number(degree), modulus.slice(0, -1).map(Number), 'a') : field(p);
    const decode = (v: bigint) => (degree > 1n ? K.fromInteger(v) : K.__call__(v)),
      encode = (v: Any) => (degree > 1n ? String(v.integer_representation()) : String(v));
    const E = EllipticCurve(K, coefficients.map(decode) as Any),
      trace: Any[] = [];
    const prototype: Any = IsomorphismTracePolynomial.prototype,
      old = prototype.roots;
    try {
      prototype.roots = function (this: Any, options: Any) {
        trace.push([this.coeffs.map(encode), options?.multiplicities ?? true]);
        return old.call(this, options);
      };
      const value = operation === 0n ? E.is_isomorphic(E) : E.isomorphism_to(E).map(encode);
      return JSON.stringify({ value: [value, trace] });
    } finally {
      prototype.roots = old;
    }
  } catch (e) {
    return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
  }
};

for (const row of native) {
  test(`native generic isomorphism ordering ${row.function} ${row.seed}`, () => {
    // Preserve native integer-list literals before parsing 61/89-bit inputs.
    const raw = row.args.map((arg) =>
      arg.startsWith('[') ? JSON.parse(arg.replace(/-?\d+/g, '"$&"')) : arg
    );
    const decode = (v: Any): Any => (Array.isArray(v) ? v.map(decode) : BigInt(v));
    if (
      (row.function === 'ec_generic_isomorphism_order' ||
        row.function === 'wm_isomorphism_comparisons') &&
      (BigInt(raw[0]) > 1_000_000n || Number(raw[1]) >= 8)
    ) {
      const script = `
        import { GF } from '${import.meta.dir}/../../rings/finite_rings/finite_field_constructor.ts';
        import { GFpn, PrimeField } from '${import.meta.dir}/../../rings/finite_rings/finite_field_extension.ts';
        import { GF2 } from '${import.meta.dir}/../../rings/finite_rings/gf2.ts';
        import { QQ } from '${import.meta.dir}/../../rings/rational_field.ts';
        import { EllipticCurve } from '${import.meta.dir}/constructor.ts';
        import { _isomorphisms, baseWI, WeierstrassIsomorphism } from '${import.meta.dir}/weierstrass_morphism.ts';
        const field = p => p === 0n ? QQ : GF(p);
        const decode = v => Array.isArray(v) ? v.map(decode) : BigInt(v);
        console.log((${functions[row.function]!.toString()})(...${JSON.stringify(raw)}.map(decode)));`;
      const output = execFileSync(process.execPath, ['--eval', script], {
        encoding: 'utf8',
        timeout: 10_000,
      });
      expect(output.trim()).toBe(row.result);
    } else {
      expect(functions[row.function]!(...raw.map(decode))).toBe(row.result);
    }
  }, 15_000);
}

test('native rich-comparison values also apply to punctuation aliases', () => {
  const row = native.find((r) => r.seed === 28800000)!;
  const expected = JSON.parse(row.result).value;
  const K = GFpn(2n, 2, [1, 1], 'a');
  const E = EllipticCurve(K, [0n, 0n, 1n, 0n, 0n]);
  const morphisms = [..._isomorphisms(E, E)].map((t) => new WeierstrassIsomorphism(E, t, E));
  for (let i = 0; i < 2; i++)
    for (let j = 0; j < 2; j++) {
      ['<', '<=', '==', '!=', '>', '>='].forEach((op, k) => {
        expect(WeierstrassIsomorphism._comparison_impl(morphisms[i]!, morphisms[j]!, op)).toBe(
          expected[0][i][j][k]
        );
      });
    }
  expect(WeierstrassIsomorphism._comparison_impl(null, morphisms[0], 'eq')).toBe(expected[3][0][2]);
  for (const op of ['toString', 'constructor', '__proto__', 'invalid']) {
    expect(() => WeierstrassIsomorphism._comparison_impl(morphisms[0]!, morphisms[1]!, op)).toThrow(
      `unsupported comparison operator: ${op}`
    );
  }
});

import { PolynomialRing as IsomorphismParentPolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
