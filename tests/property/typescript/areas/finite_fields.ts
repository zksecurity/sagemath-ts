import { FpXQ_pow, FpXQ_powBig, FpXQ_inv } from '../../../../packages/parigp-ts/src/index.js';
/**
 * sagemath-ts side of the `finite_fields` property-test area.
 *
 * Cases: tests/property/cases/finite_fields.cases.json
 * SageMath counterpart: tests/property/python/areas/finite_fields.py
 */

import { Fp_order } from '../../../../packages/parigp-ts/src/ff.js';
import { Z_isanypower } from '../../../../packages/parigp-ts/src/index.js';
import { primitive_root } from '../../../../packages/sagemath-ts/src/arith/misc.js';
import { order_from_multiple } from '../../../../packages/sagemath-ts/src/groups/generic.js';
import { GF, generic_discrete_log, sqrt_mod } from '../../../../packages/sagemath-ts/src/index.js';
import { set_random_seed } from '../../../../packages/sagemath-ts/src/misc/randstate.js';
import {
  conway_polynomial,
  has_conway_polynomial,
} from '../../../../packages/sagemath-ts/src/rings/finite_rings/conway_polynomials.js';
import {
  analyzeFiniteFieldOrder,
  isValidFiniteFieldOrder,
} from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_constructor.js';
import {
  FiniteFieldElement as ExtensionElement,
  GFpn,
  PrimeField,
  PrimeFieldElement,
} from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
import {
  FiniteFieldPrime,
  FiniteFieldElement as LegacyElement,
} from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_prime.js';
import { Zmod } from '../../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod_ring.js';
import { Integer, ZZ } from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { Polynomial } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_ring.js';
import { QQ } from '../../../../packages/sagemath-ts/src/rings/rational_field.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';

export const functions: Record<string, (...args: never[]) => unknown> = {
  // Prime field arithmetic
  ff_add: (p: bigint, a: bigint, b: bigint) => {
    const F = GF(p);
    return F.__call__(a).add(F.__call__(b)).value;
  },
  ff_mul: (p: bigint, a: bigint, b: bigint) => {
    const F = GF(p);
    return F.__call__(a).mul(F.__call__(b)).value;
  },
  ff_inv: (p: bigint, a: bigint) => {
    const F = GF(p);
    return F.__call__(a).inv().value;
  },
  ff_pow: (p: bigint, a: bigint, n: bigint) => {
    const F = GF(p);
    return F.__call__(a).pow(n).value;
  },
  ff_pow_neg: (p: bigint, a: bigint, n: bigint) => {
    const F = GF(p);
    return F.__call__(a).pow(n).value;
  },
  // Square root
  sqrt_mod: (a: bigint, p: bigint) => {
    const r = sqrt_mod(a, p);
    if (r === null) return null;
    // Return canonical (smaller) root
    return r > (p - 1n) / 2n ? p - r : r;
  },
  sqrt_mod_p3mod4: (a: bigint, p: bigint) => {
    const r = sqrt_mod(a, p);
    if (r === null) return null;
    return r > (p - 1n) / 2n ? p - r : r;
  },
  sqrt_mod_p5mod8: (a: bigint, p: bigint) => {
    const r = sqrt_mod(a, p);
    if (r === null) return null;
    return r > (p - 1n) / 2n ? p - r : r;
  },
  sqrt_mod_general: (a: bigint, p: bigint) => {
    const r = sqrt_mod(a, p);
    if (r === null) return null;
    return r > (p - 1n) / 2n ? p - r : r;
  },
  sqrt_mod_nonresidue: (a: bigint, p: bigint) => {
    const r = sqrt_mod(a, p);
    if (r === null) return null;
    return r > (p - 1n) / 2n ? p - r : r;
  },
  // Primitive roots
  primitive_root: (n: bigint) => primitive_root(n),
  ff_multiplicative_generator: (p: bigint) => {
    const F = GF(p);
    return F.multiplicative_generator().value;
  },
  // Discrete logarithm
  discrete_log: (p: bigint, base: bigint, target: bigint) => {
    const F = GF(p);
    const g = F.__call__(base);
    const h = F.__call__(target);
    return generic_discrete_log(h, g, p - 1n, '*');
  },
  // Extension field operations
  ff_ext_add: (p: bigint, n: bigint, a: bigint, b: bigint) => {
    const F = GFpn(p, Number(n));
    const x = F.fromInteger(a);
    const y = F.fromInteger(b);
    return x.add(y).integer_representation();
  },
  ff_ext_mul: (p: bigint, n: bigint, a: bigint, b: bigint) => {
    const F = GFpn(p, Number(n));
    const x = F.fromInteger(a);
    const y = F.fromInteger(b);
    return x.mul(y).integer_representation();
  },
  ff_ext_inv: (p: bigint, n: bigint, a: bigint) => {
    const F = GFpn(p, Number(n));
    const x = F.fromInteger(a);
    return x.inv().integer_representation();
  },
  ff_ext_pow: (p: bigint, n: bigint, a: bigint, e: bigint) => {
    const F = GFpn(p, Number(n));
    const x = F.fromInteger(a);
    return x.pow(e).integer_representation();
  },
  ff_ext_order: (p: bigint, n: bigint) => GFpn(p, Number(n)).cardinality(),
  ff_ext_frobenius: (p: bigint, n: bigint, a: bigint) => {
    const F = GFpn(p, Number(n));
    const x = F.fromInteger(a);
    return x.frobenius().integer_representation();
  },
  ff_generator_order: (p: bigint, n: bigint) => {
    const F = GFpn(p, Number(n));
    return order_from_multiple(F.multiplicative_generator(), F.cardinality() - 1n, undefined, '*');
  },
};

function invoke(x: object, name: string, ...args: unknown[]): unknown {
  return (x as Record<string, (...values: unknown[]) => unknown>)[name]!.apply(x, args);
}

// This port uses the PARI polynomial-basis backend. The paired oracle selects
// that backend and the same Conway modulus, including backend-specific errors.
function comparison(run: () => unknown): string {
  function normalized(x: unknown): unknown {
    if (Array.isArray(x)) return x.map(normalized);
    if (typeof x === 'boolean' || x === null) return x;
    return String(x);
  }
  try {
    return JSON.stringify({ value: normalized(run()) });
  } catch (e) {
    return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
  }
}
const extensionUnary = [
  '_integer_',
  'neg',
  'isZero',
  'isOne',
  'trace',
  'norm',
  'coefficients',
  'integer_representation',
  'toString',
  'repr',
  'minpoly',
  'minimal_polynomial',
  'minimalPolynomial',
];
for (const method of extensionUnary) {
  functions[`ff_ext_${method}`] = (p: bigint, n: bigint, a: bigint) =>
    comparison(() => {
      const x = GFpn(p, Number(n)).fromInteger(a);
      const value = invoke(x, method);
      if (value instanceof Polynomial) return value.coeffs.map((c) => String(c));
      return value;
    });
}
functions.ff_ext_fromInteger = (p: bigint, n: bigint, a: bigint) =>
  comparison(() => GFpn(p, Number(n)).fromInteger(a).integer_representation());
functions.ff_ext_iteration = (p: bigint, n: bigint) =>
  comparison(() => [...GFpn(p, Number(n))].map((x) => x.integer_representation()));
functions.ff_ext_parent = (p: bigint, n: bigint) =>
  comparison(() => {
    const F = GFpn(p, Number(n));
    return [
      F.zero(),
      F.one(),
      F.gen(),
      F.cardinality(),
      F.is_field(),
      F.toString(),
      F.characteristic,
      F.degree,
      F.order,
      [...F.elements()],
      F.primitiveElement(),
      F.primitive_element(),
      F.multiplicative_generator(),
    ];
  });
functions.ff_scalar_constructor = (p: bigint, n: bigint, mode: bigint, a: bigint, b: bigint) =>
  comparison(() => {
    const F = n === 0n ? new FiniteFieldPrime(p) : GF(p ** n);
    let x: unknown;
    switch (Number(mode)) {
      case 0:
        x = new Integer(a);
        break;
      case 1:
        x = new Rational(a, b);
        break;
      case 2:
        x = a !== 0n;
        break;
      case 3:
        x = null;
        break;
      case 4:
        x = Number(a) / Number(b);
        break;
      case 5:
        x = a.toString();
        break;
      case 6:
        x = [new Rational(a, b), 1n];
        break;
      case 7:
        x = [];
        break;
      default:
        throw new Error('Unknown constructor mode');
    }
    return F.__call__(x).toString();
  });

const primeUnary = [
  '_integer_',
  'neg',
  'inv',
  'isZero',
  'isOne',
  'is_square',
  'sqrt',
  'toString',
  'repr',
  'toBigInt',
];
for (const method of primeUnary) {
  functions[`ff_prime_${method}`] = (p: bigint, a: bigint, legacy: bigint) =>
    comparison(() => {
      const F = legacy ? new FiniteFieldPrime(p) : new PrimeField(p);
      return method === 'sqrt'
        ? invoke(F.__call__(a), method, { extend: false })
        : invoke(F.__call__(a), method);
    });
}
for (const method of ['isUnit', 'lift', 'multiplicative_order']) {
  functions[`ff_prime_${method}`] = (p: bigint, a: bigint) =>
    comparison(() => invoke(new FiniteFieldPrime(p).__call__(a), method));
}
for (const method of ['add', 'sub', 'mul', 'div', 'eq', 'pow']) {
  functions[`ff_prime_${method}`] = (p: bigint, a: bigint, b: bigint, legacy: bigint) =>
    comparison(() => {
      const F = legacy ? new FiniteFieldPrime(p) : new PrimeField(p);
      return invoke(F.__call__(a), method, b);
    });
}
functions.ff_ext_conversion = (p: bigint, n: bigint, a: bigint, mode: bigint) =>
  comparison(() => {
    const F = GFpn(p, Number(n));
    let x: unknown;
    switch (Number(mode)) {
      case 0:
        x = GFpn(p, Number(n)).fromInteger(a);
        break;
      case 1:
        x = GFpn(p, Number(n), undefined, 'b').fromInteger(a);
        break;
      case 2:
        x = GFpn(p, Number(n) + 1, undefined, 'b').fromInteger(a);
        break;
      case 3:
        x = new PrimeField(p === 2n ? 3n : 2n).__call__(a);
        break;
      case 4:
        x = new Polynomial(
          [new PrimeField(p).__call__(a), new PrimeField(p).__call__(1n)],
          new PolynomialRing(new PrimeField(p), 'z')
        );
        break;
      case 5: {
        const P = new PrimeField(p === 2n ? 3n : 2n);
        x = new Polynomial([P.__call__(a), P.one()], new PolynomialRing(P, 'z'));
        break;
      }
      case 6:
        return new ExtensionElement(
          new Polynomial(
            [F.baseField.__call__(a), F.baseField.one()],
            new PolynomialRing(F.baseField, 'z')
          ),
          F
        ).toString();
      default:
        throw new Error('Unknown conversion mode');
    }
    return F.__call__(x).toString();
  });

functions.ff_prime_constructor = (p: bigint, legacy: bigint) =>
  comparison(() => (legacy ? new FiniteFieldPrime(p) : new PrimeField(p)).toString());
functions.ff_prime_parent = (p: bigint, legacy: bigint) =>
  comparison(() => {
    const F = legacy ? new FiniteFieldPrime(p) : new PrimeField(p);
    return [
      F.zero(),
      F.one(),
      F.gen(),
      F.cardinality(),
      F.is_field(),
      F.toString(),
      F.characteristic,
      F.degree,
      F.order,
      [...F],
      [...F.elements()],
      F.multiplicative_generator(),
      legacy ? (F as FiniteFieldPrime).list() : (F as PrimeField).primitive_element(),
    ];
  });
functions.ff_direct_element = (p: bigint, a: bigint, b: bigint, mode: bigint, legacy: bigint) =>
  comparison(() => {
    const value =
      mode === 0n
        ? new Integer(a)
        : mode === 1n
          ? new Rational(a, b)
          : mode === 2n
            ? null
            : a !== 0n;
    return legacy
      ? new LegacyElement(value, new FiniteFieldPrime(p)).toString()
      : new PrimeFieldElement(value, new PrimeField(p)).toString();
  });
functions.ff_random = (p: bigint, n: bigint, seed: bigint) =>
  comparison(() => {
    const F =
      n === 0n ? new FiniteFieldPrime(p) : n === 1n ? new PrimeField(p) : GFpn(p, Number(n));
    set_random_seed(seed);
    return Array.from({ length: 20 }, () => F.random_element().toString());
  });
functions.ff_pari_order = (p: bigint, a: bigint, multiple: bigint) =>
  Fp_order(a, (p - 1n) * multiple, p);
for (const method of ['add', 'sub', 'mul', 'div', 'eq']) {
  functions[`ff_scalar_${method}`] = (
    p: bigint,
    n: bigint,
    a: bigint,
    b: bigint,
    d: bigint,
    mode: bigint
  ) =>
    comparison(() => {
      const F =
        n === 0n ? new FiniteFieldPrime(p) : n === 1n ? new PrimeField(p) : GFpn(p, Number(n));
      const value = mode === 0n ? new Integer(b) : mode === 1n ? new Rational(b, d) : b !== 0n;
      return invoke(F.__call__(a), method, value);
    });
}

functions.ff_prime_p = (p: bigint, a: bigint) => new FiniteFieldPrime(p).__call__(a).p;
functions.ff_prime_quadratic_non_residue = (p: bigint) =>
  new FiniteFieldPrime(p).quadratic_non_residue().value;
for (const method of ['sub', 'div', 'eq']) {
  functions[`ff_ext_${method}`] = (p: bigint, n: bigint, a: bigint, b: bigint) =>
    comparison(() => {
      const F = GFpn(p, Number(n));
      return invoke(F.fromInteger(a), method, F.fromInteger(b));
    });
}

const factoryAlgorithms = [
  'conway',
  'first_lexicographic',
  'minimal_weight',
  'adleman-lenstra',
  'primitive',
  'unknown',
];
functions.ff_factory = (q: bigint, mode: bigint, codes: bigint[], coeffs: bigint[]) =>
  comparison(() => {
    const name = String.fromCodePoint(...codes.map(Number));
    const factory = GF as unknown as (...args: unknown[]) => PrimeField | ReturnType<typeof GFpn>;
    let F: PrimeField | ReturnType<typeof GFpn>;
    if (mode === 0n) F = factory(q, name);
    else if (mode === 1n) F = factory(q, { name, modulus: coeffs });
    else if (mode === 2n) F = factory(q, name, { modulus: coeffs });
    else if (mode === 3n) F = factory(q, { name, modulus: factoryAlgorithms[Number(coeffs[0])] });
    else if (mode === 4n) {
      const base = new PrimeField(BigInt(coeffs[0]!));
      const R = new PolynomialRing(base, 'z');
      F = factory(q, name, {
        modulus: new Polynomial(
          coeffs.slice(1).map((c) => base.__call__(c)),
          R
        ),
      });
    } else if (mode === 5n) F = factory(q, { name });
    else if (mode === 6n) F = factory(new Integer(q), name);
    else if (mode === 7n) F = factory(q, { name, check: false });
    else throw new Error('Unknown factory mode');
    return [
      F.toString(),
      F.gen().toString(),
      F.gen().pow(2n).toString(),
      F.cardinality(),
      F instanceof PrimeField ? [] : F.modulus.coeffs.map((c) => c.value),
    ];
  });

functions.ff_conway = (p: bigint, n: bigint) =>
  comparison(() => [
    has_conway_polynomial(Number(p), Number(n)),
    conway_polynomial(Number(p), Number(n)),
  ]);
functions.ff_perfect_power_backend = (n: bigint) =>
  comparison(() => {
    // This factory adapter compares magnitude decomposition, like pari(abs(n)).
    const [k, p] = Z_isanypower(n < 0n ? -n : n);
    return [p, k || 1];
  });

functions.ff_order_helpers = (q: bigint) =>
  comparison(() => {
    const info = analyzeFiniteFieldOrder(q);
    return [
      isValidFiniteFieldOrder(q),
      info && [info.order, info.characteristic, info.degree, info.isPrimeField],
    ];
  });

functions.ff_sqrt_options = (p: bigint, a: bigint, mode: bigint, legacy: bigint) =>
  comparison(() => {
    const F = legacy ? new FiniteFieldPrime(p) : new PrimeField(p);
    const x = F.__call__(a);
    const options =
      mode === 0n
        ? undefined
        : mode === 1n
          ? { extend: false }
          : mode === 2n
            ? { extend: true }
            : mode === 3n
              ? { all: true }
              : mode === 4n
                ? { extend: false, all: true }
                : mode === 5n
                  ? { extend: true, all: true }
                  : { unknown: true };
    const result = options === undefined ? invoke(x, 'sqrt') : invoke(x, 'sqrt', options);
    if (Array.isArray(result)) return result.map(String);
    const root = result as PrimeFieldElement | LegacyElement | ExtensionElement;
    return [root.toString(), root.pow(2n).toString(), root.parent.toString()];
  });

functions.ff_string_expression = (p: bigint, n: bigint, codes: bigint[], polynomial: bigint) =>
  comparison(() => {
    const F = GFpn(p, Number(n));
    const input = String.fromCodePoint(...codes.map(Number));
    if (polynomial === 2n) return F.polynomialRing.toString();
    if (polynomial) {
      const value = (F.polynomialRing.__call__ as (x: unknown) => Polynomial<PrimeFieldElement>)(
        input
      );
      return value.coeffs.map(String);
    }
    return F.__call__(input).toString();
  });

functions.ff_prime_power_scalar = (
  p: bigint,
  a: bigint,
  mode: bigint,
  b: bigint,
  d: bigint,
  legacy: bigint
) =>
  comparison(() => {
    const exponent =
      mode === 0n
        ? new Integer(b)
        : mode === 1n
          ? new Rational(b, d)
          : mode === 2n
            ? Boolean(b)
            : mode === 3n
              ? null
              : mode === 4n
                ? String(b)
                : mode === 6n
                  ? []
                  : Number(b) / Number(d);
    return invoke(
      (legacy ? new FiniteFieldPrime(p) : new PrimeField(p)).__call__(a),
      'pow',
      exponent
    );
  });

for (const method of ['add', 'sub', 'mul', 'div', 'eq']) {
  functions['ff_mixed_' + method] = (
    kind: bigint,
    p: bigint,
    n: bigint,
    a: bigint,
    otherKind: bigint,
    q: bigint,
    m: bigint,
    b: bigint
  ) =>
    comparison(() => {
      const x =
        kind === 2n
          ? GFpn(p, Number(n)).fromInteger(a)
          : kind === 3n
            ? Zmod(p).__call__(a)
            : kind === 1n
              ? new FiniteFieldPrime(p).__call__(a)
              : new PrimeField(p).__call__(a);
      const y =
        otherKind === 0n
          ? new Integer(b)
          : otherKind === 1n
            ? new Rational(b, m)
            : otherKind === 2n
              ? Boolean(b)
              : otherKind === 3n
                ? null
                : otherKind === 4n
                  ? String(b)
                  : otherKind === 5n
                    ? Number(b) / Number(m)
                    : otherKind === 6n
                      ? Zmod(q).__call__(b)
                      : otherKind === 7n
                        ? new PrimeField(q).__call__(b)
                        : otherKind === 8n
                          ? GFpn(q, Number(m)).fromInteger(b)
                          : otherKind === 10n
                            ? new FiniteFieldPrime(q).__call__(b)
                            : b === 0n
                              ? []
                              : [b];
      const result = invoke(x as object, method, y);
      return typeof result === 'boolean'
        ? result
        : [
            Array.isArray(result) ? '[' + result.map(String).join(', ') + ']' : String(result),
            typeof result === 'string'
              ? "<class 'str'>"
              : Array.isArray(result)
                ? "<class 'list'>"
                : String((result as { parent: unknown }).parent),
          ];
    });
}

functions.ff_sequence_multiply = (kind: bigint, p: bigint, a: bigint, mode: bigint) =>
  comparison(() => {
    const R =
      kind === 3n
        ? Zmod(p)
        : kind === 2n
          ? GFpn(p, 2)
          : kind === 1n
            ? new FiniteFieldPrime(p)
            : new PrimeField(p);
    const sequence = ['', 'ab', '𝄞é', [], [2n, 7n]][Number(mode)];
    const result = invoke(R.__call__(a) as object, 'mul', sequence);
    return [result, typeof result === 'string' ? "<class 'str'>" : "<class 'list'>"];
  });

function parentScalar(mode: bigint, a: bigint, b: bigint): unknown {
  if (mode === 0n) return new Integer(a);
  if (mode === 1n) return new Rational(a, b);
  if (mode === 2n) return Boolean(a);
  if (mode === 3n) return null;
  if (mode === 4n) return '';
  if (mode === 5n) return String(a);
  if (mode === 6n) return Number(a) / Number(b);
  if (mode === 7n) return a === 0n ? [] : [a];
  if (mode === 8n) return Zmod(b).__call__(a);
  if (mode === 9n) return new PrimeField(b).__call__(a);
  return GFpn(b, 2).fromInteger(a);
}
functions.ff_generator_index = (kind: bigint, p: bigint, mode: bigint, a: bigint, b: bigint) =>
  comparison(() => {
    const R =
      kind === 4n
        ? QQ
        : kind === 3n
          ? Zmod(p)
          : kind === 2n
            ? GFpn(p, 2)
            : kind === 1n
              ? new FiniteFieldPrime(p)
              : new PrimeField(p);
    return invoke(R, 'gen', parentScalar(mode, a, b));
  });
functions.ff_parent_constructor = (
  kind: bigint,
  mode: bigint,
  a: bigint,
  b: bigint,
  check: bigint
) =>
  comparison(() => {
    const p = parentScalar(mode, a, b);
    const R = kind
      ? new FiniteFieldPrime(p as bigint, Boolean(check))
      : new PrimeField(p as bigint);
    return [String(R), R.characteristic, R.order, R.zero(), R.one(), R.is_field()];
  });

functions.ff_integer_conversion = (kind: bigint, p: bigint, n: bigint, a: bigint) =>
  comparison(() => {
    const x = kind === 2n ? GFpn(p, Number(n)).fromInteger(a) : new PrimeField(p).__call__(a);
    return [ZZ.__call__(x), new Integer(x)];
  });

for (const method of ['pow', 'powBig', 'inv']) {
  functions['ff_pari_quotient_' + method] = (
    p: bigint,
    coeffs: bigint[],
    modulus: bigint[],
    exponent: bigint
  ) =>
    comparison(() => {
      return method === 'inv'
        ? FpXQ_inv(coeffs, modulus, p)
        : method === 'powBig'
          ? FpXQ_powBig(coeffs, exponent, modulus, p)
          : FpXQ_pow(coeffs, exponent, modulus, p);
    });
}
functions.ff_extension_power_scalar = (
  p: bigint,
  n: bigint,
  a: bigint,
  mode: bigint,
  b: bigint,
  d: bigint
) => comparison(() => invoke(GFpn(p, Number(n)).fromInteger(a), 'pow', parentScalar(mode, b, d)));
functions.ff_pari_inverse_lift = (p: bigint, degree: bigint, coeffs: bigint[], modulus: bigint[]) =>
  comparison(() => FpXQ_inv(coeffs, modulus, p ** degree, p));

import { FF_issquareall, PariType } from '../../../../packages/parigp-ts/src/index.js';
import { setrand as rootSetrand, getrand as rootGetrand } from '../../../../packages/parigp-ts/src/random.js';
import { FiniteFieldExtension as RootField } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
functions.pari_ff_square_root = (p: bigint, T: bigint[], a: bigint[], seed: bigint) => {
  rootSetrand(seed);
  const result = FF_issquareall({type:PariType.t_FFELT,p,degree:T.length-1,value:a,definingPoly:T});
  return JSON.stringify({value:result === null? null : (typeof result.value === 'bigint'?[result.value]:result.value).map(String),state:String(rootGetrand())});
};
functions.ff_extension_sqrt = (p: bigint,T:bigint[],a:bigint[],seed:bigint,mode:bigint) => {
  const K = new PrimeField(p), R = new PolynomialRing(K,'x');
  const F = new RootField(p,T.length-1,R.__call__(T.map(c=>K.__call__(c))),'a');
  const x = F.__call__(R.__call__(a.map(c=>K.__call__(c))));
  const options=mode===0n?undefined:mode===1n?{extend:false}:mode===2n?{extend:true}:mode===3n?{all:true}:mode===4n?{extend:false,all:true}:mode===5n?{extend:true,all:true}:{unknown:true};
  rootSetrand(seed);
  try {
    const result=invoke(x,'sqrt',...(options===undefined?[]:[options]));
    const encode=(v:ExtensionElement)=>v.lift.coeffs.map(c=>String(c.value));
    return JSON.stringify({value:Array.isArray(result)?result.map(encode):encode(result as ExtensionElement),state:String(rootGetrand())});
  } catch(e) {
    return JSON.stringify({error:(e as Error).name,message:(e as Error).message,state:String(rootGetrand())});
  }
};

import { FpX_resultant, Flx_resultant, FpXQ_norm, Flxq_norm, FpXQ_issquare, Flxq_issquare, FF_norm, FF_issquare } from '../../../../packages/parigp-ts/src/index.js';
functions.pari_field_predicates=(op:bigint,p:bigint,a:bigint[],T:bigint[],seed:bigint)=>{
  rootSetrand(seed);
  const x={type:PariType.t_FFELT as PariType.t_FFELT,p,degree:T.length-1,value:a,definingPoly:T};
  const value=op===0n?FpX_resultant(a,T,p):op===1n?Flx_resultant(a,T,p):op===2n?FpXQ_norm(a,T,p):op===3n?Flxq_norm(a,T,p):op===4n?FpXQ_issquare(a,T,p):op===5n?Flxq_issquare(a,T,p):op===6n?FF_norm(x):FF_issquare(x);
  return JSON.stringify({value:typeof value==='boolean'?value:String(value),state:String(rootGetrand())});
};
functions.ff_extension_is_square=(p:bigint,T:bigint[],a:bigint[],seed:bigint)=>{
  const K=new PrimeField(p),R=new PolynomialRing(K,'x');
  const F=new RootField(p,T.length-1,R.__call__(T.map(c=>K.__call__(c))),'a');
  const x=F.__call__(R.__call__(a.map(c=>K.__call__(c))));
  rootSetrand(seed);
  return JSON.stringify({value:x.is_square(),state:String(rootGetrand())});
};
functions.pari_field_record_scalar=(op:bigint,p:bigint,a:bigint,seed:bigint)=>{
  rootSetrand(seed);
  const x={type:PariType.t_FFELT as PariType.t_FFELT,p,degree:1,value:a};
  if(op===0n){
    const root=FF_issquareall(x);
    return JSON.stringify({value:root===null?null:(typeof root.value==='bigint'?[root.value]:root.value).map(String),state:String(rootGetrand())});
  }
  const value=op===1n?FF_norm(x):FF_issquare(x);
  return JSON.stringify({value:typeof value==='boolean'?value:String(value),state:String(rootGetrand())});
};
