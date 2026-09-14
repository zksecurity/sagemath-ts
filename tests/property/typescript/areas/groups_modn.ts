import { pollard, pollard_matrix } from '../group_pollard.js';
import { group_parent } from '../group_parent.js';
import { group_schedule, rdf_schedule } from '../group_schedule.js';
import { rdf_multiple } from '../real_double_group.js';
import { integer_group } from '../integer_rational.js';
import { nf_scalar, nf_group } from '../group_number_field.js';
import { iterator_state, iterator_standard, parent_field } from '../group_iterator.js';
import { current_randstate, set_random_seed } from '../../../../packages/sagemath-ts/src/misc/randstate.js';
/**
 * sagemath-ts side of the `groups_modn` property-test area.
 *
 * Differential oracle for
 *   - packages/sagemath-ts/src/groups/generic.ts
 *   - packages/sagemath-ts/src/rings/finite_rings/integer_mod.ts
 *   - packages/sagemath-ts/src/rings/finite_rings/integer_mod_ring.ts
 *   - packages/sagemath-ts/src/rings/finite_rings/roots_of_unity.ts
 *   - packages/sagemath-ts/src/rings/finite_rings/tower_field.ts
 *
 * Cases: tests/property/cases/groups_modn.cases.json
 * SageMath counterpart: tests/property/python/areas/groups_modn.py
 *
 * Every function returns an already-formatted string, and every function
 * funnels through `run()`, which renders a thrown error as
 * `"<ErrorClass>: <message>"`.  That is deliberate: `compare.ts` scores
 * "both sides raised" as a pass regardless of *what* was raised, so letting
 * exceptions escape would make every error case vacuous.  Here a wrong error
 * class or a drifted message is a failure.
 */

import {
  bsgs,
  parseGroupOps,
  discrete_log,
  discrete_log_lambda,
  discrete_log_rho,
  has_order,
  multiple,
  multiple_of_order,
  multiples,
  order_from_bounds,
  order_from_multiple,
} from '../../../../packages/sagemath-ts/src/groups/generic.js';
import { GF } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_constructor.js';
import {
  GF2,
  type GF2Element,
} from '../../../../packages/sagemath-ts/src/rings/finite_rings/gf2.js';
import { Mod } from '../../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod.js';
import { Integers } from '../../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod_ring.js';
import {
  type CosetDomain,
  FFTDomain,
  type FiniteFieldElementLike,
  cyclotomic_polynomial,
  elementOfOrder,
  findMultiplicativeGenerator,
  has_primitive_root,
  maxFFTSize,
  multiplicative_order,
  primitive_nth_root,
  roots_of_unity,
  twoAdicity,
  validFFTSizes,
} from '../../../../packages/sagemath-ts/src/rings/finite_rings/roots_of_unity.js';
import {
  Ti,
  type TowerField,
  type TowerFieldElement,
} from '../../../../packages/sagemath-ts/src/rings/finite_rings/tower_field.js';
import { Polynomial } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_element.js';
import type {
  QuotientRing,
  QuotientRingElement,
} from '../../../../packages/sagemath-ts/src/rings/polynomial/quotient_ring.js';

// ------------------------------------------------------------- formatting

/** Call `f`, formatting either its result or the error it threw. */
function run(f: () => unknown): string {
  try {
    return String(f());
  } catch (e) {
    if (e instanceof Error) {
      return `${e.name}: ${e.message}`;
    }
    return `Error: ${String(e)}`;
  }
}

function fmtList(xs: readonly unknown[]): string {
  return `[${xs.map((x) => String(x)).join(', ')}]`;
}

function fmtBool(b: boolean): string {
  return b ? 'True' : 'False';
}

// --------------------------------------------------------- groups/generic

export const functions = {
  gg_pollard: pollard,
  gg_pollard_matrix: pollard_matrix,
  gg_parent: group_parent,
  gg_schedule: group_schedule,
  gg_rdf_schedule: rdf_schedule,
  gg_rdf_multiple: rdf_multiple,
  gg_integer: integer_group,
  gg_nf_scalar: nf_scalar,
  gg_nf_group: nf_group,
  gg_iterator_state: iterator_state,
  gg_iterator_standard: iterator_standard,
  gg_parent_field: parent_field,

  gg_order_list_trace: (modulus: bigint, value: bigint, multipleValue: bigint, factorMode: bigint, packedFactors: bigint[], plistMode: bigint, plist: bigint[]) => {
    const trace: number[][] = [];
    const factors: [bigint, bigint][] | undefined = factorMode === 0n ? undefined : [];
    if (factors) for (let i=0; i<packedFactors.length; i+=2) factors.push([packedFactors[i]!, packedFactors[i+1]!]);
    let result: string | null = null, error: string | null = null;
    try {
      result = String(order_from_multiple(Mod(value, modulus), multipleValue, factors, 'other', Mod(0n, modulus), x => x.neg(),
        (x,y) => { trace.push([Number(x.value), Number(y.value)]); return x.add(y); },
        { plist: plistMode === 0n ? undefined : plist }));
    } catch (e) { error = e instanceof Error ? `${e.name}: ${e.message}` : String(e); }
    return JSON.stringify([result, error, trace]);
  },

  gg_group_parse: (fn: bigint, mode: bigint, mask: bigint, variant: bigint) => {
    set_random_seed(0);
    const a = Mod(2n, 11n), b = Mod(4n, 11n), v = Number(variant);
    const operation = ['+', '*', 'other', 'custom'][Number(mode)] as 'other' | '+' | '*';
    const identity = mask & 1n ? Mod(0n, 11n) : undefined;
    const inverse = mask & 2n ? (x: typeof a) => x.neg() : undefined;
    const op = mask & 4n ? (x: typeof a, y: typeof a) => x.add(y) : undefined;
    const kw = [operation, identity, inverse, op] as const;
    const bounds = [[0n, 8n], [-1n, 8n], [8n, 0n]][v] as [bigint, bigint];
    let result: string | null = null, error: string | null = null;
    try {
      let r: unknown;
      if (fn === 0n) r = multiple(a, [-1n, 0n, 3n][v]!, ...kw);
      else if (fn === 1n) r = bsgs(a, b, bounds, ...kw);
      else if (fn === 2n) r = discrete_log(b, a, [operation === '*' ? 10n : 11n, 0n, -1n][v]!, ...kw);
      else if (fn === 3n) r = order_from_multiple(a, [110n, 3n, 0n][v]!, undefined, ...kw);
      else if (fn === 4n) r = order_from_bounds(a, [[1n,22n],[-1n,22n],[22n,1n]][v] as [bigint,bigint], undefined, ...kw);
      else if (fn === 5n) r = discrete_log_lambda(b, a, bounds, ...kw, x => x.value ** 2n);
      else if (fn === 6n) r = discrete_log_rho(b, a, [11n, 0n, -1n][v]!, ...kw);
      else {
        const o = parseGroupOps(...kw, a);
        r = JSON.stringify([Number(o.identity.value), Number(o.inverse(a).value), Number(o.op(a, b).value)]);
      }
      result = String(r);
    } catch (e) { error = e instanceof Error ? `${e.name}: ${e.message}` : String(e); }
    return JSON.stringify([result, error]);
  },
  gg_bounds_trace: (value: bigint, lb: bigint, ub: bigint, d: bigint) => {
    const trace: number[][] = [];
    let result: string | null = null, error: string | null = null;
    try {
      result = String(order_from_bounds(Mod(value, 11n), [lb, ub], d, 'other', Mod(0n, 11n), x => x.neg(),
        (x,y) => { trace.push([Number(x.value), Number(y.value)]); return x.add(y); }));
    } catch (e) { error = e instanceof Error ? `${e.name}: ${e.message}` : String(e); }
    return JSON.stringify([result, error, trace]);
  },
  gg_multiple_trace: (p: bigint, value: bigint, mode: bigint, n: bigint) => {
    const trace: number[][] = [];
    let result: string | null = null, error: string | null = null;
    try {
      result = String(multiple(Mod(value, p), n, 'other', Mod(mode === 1n ? 1n : mode === 2n ? p-1n : 0n, p),
        x => mode === 1n ? x.inv() : x.neg(), (x,y) => {
          trace.push([Number(x.value), Number(y.value)]);
          return mode === 0n ? x.add(y) : mode === 1n ? x.mul(y) : x.value < y.value ? x : y;
        }));
    } catch (e) { error = e instanceof Error ? `${e.name}: ${e.message}` : String(e); }
    return JSON.stringify([result, error, trace]);
  },

  gg_rho_state: (seed: bigint, mode: bigint, p: bigint, b: bigint, order: bigint, target: bigint) => {
    set_random_seed(seed);
    const base = Mod(b, p);
    const operation = (['+', '*', 'other'] as const)[Number(mode)]!;
    let result: string | null = null;
    let error: string | null = null;
    try {
      result = String(discrete_log_rho(
        Mod(target, p), base, order, operation,
        operation === 'other' ? Mod(0n, p) : undefined,
        operation === 'other' ? v => v.neg() : undefined,
        operation === 'other' ? (u, v) => u.add(v) : undefined
      ));
    } catch (e) {
      error = e instanceof Error ? `${e.name}: ${e.message}` : String(e);
    }
    return JSON.stringify([result, error, String(current_randstate().python_random().getrandbits(64))]);
  },

  gg_lambda_trace: (seed: bigint, mode: bigint, p: bigint, b: bigint, lb: bigint, ub: bigint, x: bigint, hashMode: bigint) => {
    set_random_seed(seed);
    const trace: number[] = [];
    const base = Mod(b, p);
    const operation = (['+', '*', 'other'] as const)[Number(mode)]!;
    let result: string | null = null;
    let error: string | null = null;
    try {
      result = String(discrete_log_lambda(
        operation === '*' ? base.pow(x) : base.mul(x), base, [lb, ub], operation,
        operation === 'other' ? Mod(0n, p) : undefined,
        operation === 'other' ? v => v.neg() : undefined,
        operation === 'other' ? (u, v) => u.add(v) : undefined,
        v => {
          trace.push(Number(v.value));
          const z = v.value ** 2n + 1n;
          return [z, -z, -1n, -(1n << 100n) + z, 0n][Number(hashMode)]!;
        }
      ));
    } catch (e) {
      error = e instanceof Error ? `${e.name}: ${e.message}` : String(e);
    }
    return JSON.stringify([result, error, trace]);
  },

  gg_bsgs_mul: (p: bigint, a: bigint, b: bigint, lb: bigint, ub: bigint) =>
    run(() => bsgs(Mod(a, p), Mod(b, p), [lb, ub], '*')),

  gg_bsgs_add: (n: bigint, a: bigint, b: bigint, lb: bigint, ub: bigint) =>
    run(() => bsgs(Mod(a, n), Mod(b, n), [lb, ub], '+')),

  gg_multiple_mul: (p: bigint, a: bigint, k: bigint) => run(() => multiple(Mod(a, p), k, '*')),

  gg_multiple_add: (n: bigint, a: bigint, k: bigint) => run(() => multiple(Mod(a, n), k, '+')),

  gg_dlog_mul: (p: bigint, a: bigint, b: bigint, ord: bigint) =>
    run(() => discrete_log(Mod(a, p), Mod(b, p), ord, '*')),

  gg_dlog_mul_noord: (p: bigint, a: bigint, b: bigint) =>
    run(() => discrete_log(Mod(a, p), Mod(b, p), undefined, '*')),

  gg_dlog_add: (n: bigint, a: bigint, b: bigint, ord: bigint) =>
    run(() => discrete_log(Mod(a, n), Mod(b, n), ord, '+')),

  gg_order_from_multiple_mul: (n: bigint, a: bigint, m: bigint) =>
    run(() => order_from_multiple(Mod(a, n), m, undefined, '*')),

  // No `operation`: exercises the port's default, which must be '+' like Sage.
  gg_order_from_multiple_default: (n: bigint, a: bigint, m: bigint) =>
    run(() => order_from_multiple(Mod(a, n), m)),

  gg_order_from_bounds_mul: (n: bigint, a: bigint, lb: bigint, ub: bigint) =>
    run(() => order_from_bounds(Mod(a, n), [lb, ub], undefined, '*')),

  gg_order_from_bounds_d: (n: bigint, a: bigint, lb: bigint, ub: bigint, d: bigint) =>
    run(() => order_from_bounds(Mod(a, n), [lb, ub], d, '*')),

  gg_order_from_bounds_add: (n: bigint, a: bigint, lb: bigint, ub: bigint) =>
    run(() => order_from_bounds(Mod(a, n), [lb, ub], undefined, '+')),

  gg_order_from_bounds_nobounds: (n: bigint, a: bigint) =>
    run(() => order_from_bounds(Mod(a, n), undefined, undefined, '*')),

  gg_multiple_of_order: (n: bigint, a: bigint) => run(() => multiple_of_order(Mod(a, n), '*')),

  gg_has_order_mul: (n: bigint, a: bigint, m: bigint) =>
    run(() => fmtBool(has_order(Mod(a, n), m, '*'))),

  // No `operation`: exercises the port's default, which must be '+' like Sage.
  gg_has_order_default: (n: bigint, a: bigint, m: bigint) =>
    run(() => fmtBool(has_order(Mod(a, n), m))),

  // Sage raises ValueError('unknown group operation') (generic.py:1584).
  gg_has_order_other: (n: bigint, a: bigint, m: bigint) =>
    run(() => fmtBool(has_order(Mod(a, n), m, 'other'))),

  // No `P0`/`indexed`/`operation`: Sage yields 0, a, 2a, ... (additive).
  gg_multiples_default: (n: bigint, a: bigint, k: bigint) =>
    run(() => fmtList([...multiples(Mod(a, n), k)])),

  gg_multiples_mul: (p: bigint, a: bigint, k: bigint) =>
    run(() => fmtList([...multiples(Mod(a, p), k, undefined, false, '*')])),

  gg_multiples_indexed: (n: bigint, a: bigint, k: bigint) =>
    run(() =>
      fmtList([...multiples(Mod(a, n), k, undefined, true)].map(([i, e]) => `(${i}, ${e})`))
    ),

  gg_dlog_lambda: (p: bigint, base: bigint, x: bigint, lb: bigint, ub: bigint) =>
    run(() => discrete_log_lambda(Mod(base, p).pow(x), Mod(base, p), [lb, ub], '*')),

  gg_dlog_rho: (p: bigint, base: bigint, x: bigint, ord: bigint) =>
    run(() => discrete_log_rho(Mod(base, p).pow(x), Mod(base, p), ord, '*')),

  // ------------------------------------------------------------ integer_mod

  im_log: (n: bigint, a: bigint, b: bigint) => run(() => Mod(a, n).log(Mod(b, n))),

  im_log_nobase: (n: bigint, a: bigint) => run(() => Mod(a, n).log()),

  im_log_check: (n: bigint, a: bigint, b: bigint, order: bigint) =>
    run(() => Mod(a, n).log(Mod(b, n), order, { check: true })),

  im_mult_order: (n: bigint, a: bigint) => run(() => Mod(a, n).multiplicative_order()),

  im_inv: (n: bigint, a: bigint) => run(() => Mod(a, n).inv().value),

  im_div: (n: bigint, a: bigint, b: bigint) => run(() => Mod(a, n).div(Mod(b, n)).value),

  im_pow: (n: bigint, a: bigint, e: bigint) => run(() => Mod(a, n).pow(e).value),

  im_add: (n: bigint, a: bigint, b: bigint) => run(() => Mod(a, n).add(Mod(b, n)).value),

  im_sub: (n: bigint, a: bigint, b: bigint) => run(() => Mod(a, n).sub(Mod(b, n)).value),

  im_mul: (n: bigint, a: bigint, b: bigint) => run(() => Mod(a, n).mul(Mod(b, n)).value),

  im_is_unit: (n: bigint, a: bigint) => run(() => fmtBool(Mod(a, n).isUnit())),

  im_lift: (n: bigint, a: bigint) => run(() => Mod(a, n).lift()),

  // ------------------------------------------------------- integer_mod_ring

  imr_unit_gens: (n: bigint) =>
    run(() =>
      fmtList(
        Integers(n)
          .unit_gens()
          .map((g) => g.value)
      )
    ),

  imr_mult_gen: (n: bigint) => run(() => Integers(n).multiplicative_generator().value),

  imr_is_cyclic: (n: bigint) => run(() => fmtBool(Integers(n).multiplicative_group_is_cyclic())),

  imr_is_field: (n: bigint) => run(() => fmtBool(Integers(n).is_field())),

  imr_cardinality: (n: bigint) => run(() => Integers(n).cardinality()),

  imr_units: (n: bigint) =>
    run(() =>
      fmtList(
        Integers(n)
          .units()
          .map((u) => u.value)
          .sort((x, y) => (x < y ? -1 : x > y ? 1 : 0))
      )
    ),

  imr_list: (n: bigint) =>
    run(() =>
      fmtList(
        Integers(n)
          .list()
          .map((x) => x.value)
      )
    ),

  // --------------------------------------------------------- roots_of_unity

  ru_has_primitive_root: (p: bigint, n: bigint) => run(() => fmtBool(has_primitive_root(GF(p), n))),

  ru_primitive_nth_root: (p: bigint, n: bigint) => run(() => valueOf(primitive_nth_root(GF(p), n))),

  ru_roots_of_unity: (p: bigint, n: bigint) =>
    run(() => fmtList(roots_of_unity(GF(p), n).map(valueOf))),

  ru_mult_order: (p: bigint, a: bigint) =>
    run(() => multiplicative_order(GF(p).__call__(a) as unknown as FiniteFieldElementLike)),

  ru_find_mult_gen: (p: bigint) => run(() => valueOf(findMultiplicativeGenerator(GF(p)))),

  ru_element_of_order: (p: bigint, n: bigint) => run(() => valueOf(elementOfOrder(GF(p), n))),

  ru_cyclotomic: (n: bigint) => run(() => fmtList(cyclotomic_polynomial(n) as number[])),

  ru_max_fft_size: (p: bigint) => run(() => maxFFTSize(GF(p))),

  ru_two_adicity: (p: bigint) => run(() => twoAdicity(GF(p))),

  ru_valid_fft_sizes: (p: bigint) => run(() => fmtList(validFFTSizes(GF(p)))),

  ru_fft_domain: (p: bigint, n: bigint) =>
    run(() => fmtList(new FFTDomain(GF(p), n).elements().map(valueOf))),

  ru_coset: (p: bigint, n: bigint, offset: bigint) =>
    run(() => {
      const F = GF(p);
      const domain = new FFTDomain(F, n);
      const coset = domain.coset(F.__call__(offset) as unknown as FiniteFieldElementLike);
      return fmtList(coset.elements().map(valueOf));
    }),

  ru_coset_fold: (p: bigint, n: bigint, offset: bigint) =>
    run(() => {
      const F = GF(p);
      const domain = new FFTDomain(F, n);
      const coset = domain.coset(F.__call__(offset) as unknown as FiniteFieldElementLike);
      // The FRI fold ignores the challenge for the *domain*; it squares the
      // offset and the generator and halves the size (roots_of_unity.ts:723).
      const folded: CosetDomain = coset.fold(F.one() as unknown as FiniteFieldElementLike);
      return fmtList(folded.elements().map(valueOf));
    }),

  // ----------------------------------------------------------- tower_field

  tf_cardinality: (i: bigint) => run(() => towerRing(i).cardinality()),

  tf_add: (i: bigint, a: bigint, b: bigint) =>
    run(() => {
      const T = towerRing(i);
      return towerEncode(i, towerAdd(towerDecode(T, i, a), towerDecode(T, i, b)));
    }),

  tf_mul: (i: bigint, a: bigint, b: bigint) =>
    run(() => {
      const T = towerRing(i);
      return towerEncode(i, towerMul(towerDecode(T, i, a), towerDecode(T, i, b)));
    }),

  tf_inv: (i: bigint, a: bigint) =>
    run(() => {
      const T = towerRing(i);
      const x = towerDecode(T, i, a);
      return towerEncode(i, towerInv(x));
    }),

  tf_pow: (i: bigint, a: bigint, e: bigint) =>
    run(() => {
      const T = towerRing(i);
      return towerEncode(i, towerPow(towerDecode(T, i, a), e));
    }),

  tf_gen: (i: bigint) => run(() => towerEncode(i, towerRing(i).gen() as TowerFieldElement)),

  tf_elements: (i: bigint) =>
    run(() =>
      fmtList([...(towerRing(i) as Iterable<TowerFieldElement>)].map((e) => towerEncode(i, e)))
    ),

  tf_mult_order: (i: bigint, a: bigint) =>
    run(() => {
      const T = towerRing(i);
      const x = towerDecode(T, i, a);
      if (towerIsZero(x)) {
        throw new ArithmeticErrorShim('Multiplicative order of 0 not defined.');
      }
      let k = 1n;
      let y = x;
      while (!towerIsOne(y)) {
        y = towerMul(y, x);
        k += 1n;
      }
      return k;
    }),
};

// ------------------------------------------------------------------ helpers

/** `Error` whose `name` matches Python's `ArithmeticError`. */
class ArithmeticErrorShim extends Error {
  override name = 'ArithmeticError';
}

/** Integer representative of a prime-field element. */
function valueOf(x: FiniteFieldElementLike | { value?: bigint }): bigint {
  const v = (x as { value?: bigint }).value;
  if (typeof v !== 'bigint') {
    throw new Error(`element has no integer representative: ${String(x)}`);
  }
  return v;
}

const towerCache = new Map<string, TowerField>();

function towerRing(i: bigint): TowerField {
  const key = i.toString();
  const cached = towerCache.get(key);
  if (cached !== undefined) {
    return cached;
  }
  const T = Ti(Number(i));
  towerCache.set(key, T);
  return T;
}

function towerAdd(x: TowerFieldElement, y: TowerFieldElement): TowerFieldElement {
  return (x as { add(o: unknown): TowerFieldElement }).add(y);
}

function towerMul(x: TowerFieldElement, y: TowerFieldElement): TowerFieldElement {
  return (x as { mul(o: unknown): TowerFieldElement }).mul(y);
}

function towerInv(x: TowerFieldElement): TowerFieldElement {
  return (x as { inv(): TowerFieldElement }).inv();
}

function towerPow(x: TowerFieldElement, e: bigint): TowerFieldElement {
  return (x as { pow(n: bigint): TowerFieldElement }).pow(e);
}

function towerIsZero(x: TowerFieldElement): boolean {
  return (x as { isZero(): boolean }).isZero();
}

function towerIsOne(x: TowerFieldElement): boolean {
  return (x as { isOne(): boolean }).isOne();
}

/** Integer -> element of T_i, little-endian in the tower basis. */
function towerDecode(T: TowerField, i: bigint, v: bigint): TowerFieldElement {
  if (i === 0n) {
    return GF2.__call__(((v % 2n) + 2n) % 2n) as TowerFieldElement;
  }
  const Q = T as QuotientRing<TowerFieldElement>;
  const S = Q.polynomial_ring.base_ring as unknown as TowerField;
  const half = 2n ** (2n ** (i - 1n));
  const lo = towerDecode(S, i - 1n, ((v % half) + half) % half);
  const hi = towerDecode(S, i - 1n, v / half);
  return Q.__call__(new Polynomial<TowerFieldElement>([lo, hi], Q.polynomial_ring));
}

/** Element of T_i -> integer (inverse of `towerDecode`). */
function towerEncode(i: bigint, e: TowerFieldElement): bigint {
  if (i === 0n) {
    return BigInt((e as GF2Element).value);
  }
  const half = 2n ** (2n ** (i - 1n));
  const elem = e as QuotientRingElement<TowerFieldElement>;
  const coeffs = elem.lift.coeffs;
  const zero = (elem.parent as QuotientRing<TowerFieldElement>).polynomial_ring.base_ring.zero();
  const c0 = (coeffs[0] ?? zero) as TowerFieldElement;
  const c1 = (coeffs[1] ?? zero) as TowerFieldElement;
  return towerEncode(i - 1n, c0) + half * towerEncode(i - 1n, c1);
}
