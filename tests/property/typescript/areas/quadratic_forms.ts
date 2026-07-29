/**
 * sagemath-ts side of the `quadratic_forms` property-test area.
 *
 * Cases: tests/property/cases/quadratic_forms.cases.json
 * SageMath counterpart: tests/property/python/areas/quadratic_forms.py
 *
 * Covers `sage.quadratic_forms.binary_qf`: Gauss composition (`__mul__` ->
 * PARI `qfbcompraw`), reduction (`reduced_form` -> PARI `qfbred` / `qfbredsl2`
 * / Sage's own `_reduce_indef`), the class-group enumerators
 * (`BinaryQF_reduced_representatives`), cycles, proper/improper equivalence
 * and `solve_integer` (PARI `qfbsolve` / `qfbcornacchia` plus Sage's
 * elementary algorithm for square discriminants).
 *
 * Every function takes flat integer lists (the runners can only generate
 * `bigint` and `bigint[]`) and returns an already-formatted newline-separated
 * string, so the two sides are compared byte-for-byte.  The formatting helpers
 * below are mirrored verbatim from the Python area module; see its docstring
 * for the three SageMath-10.3-vs-vendored-10.9 / PARI-2.15-vs-2.18 divergences
 * that the case list deliberately routes around.
 */

import {
  BinaryQF,
  BinaryQF_reduced_representatives,
  class_number,
} from '../../../../packages/sagemath-ts/src/quadratic_forms/index.js';

type Triple = [bigint, bigint, bigint];
type Matrix2 = [[bigint, bigint], [bigint, bigint]];

// ---------------------------------------------------------------------------
// formatting helpers (mirrored verbatim in the Python area module)
// ---------------------------------------------------------------------------

const F = (f: BinaryQF | Triple): string =>
  f instanceof BinaryQF ? `(${f.a},${f.b},${f.c})` : `(${f[0]},${f[1]},${f[2]})`;

const T = (x: boolean): string => (x ? 'True' : 'False');

/**
 * Run `fn`, returning its string, or `!<message>` if it threw.
 *
 * Error messages are compared byte-for-byte; `compare.ts` would otherwise
 * score "both sides raised" as a pass no matter how different the reasons.
 */
function attempt(fn: () => string): string {
  try {
    return fn();
  } catch (e) {
    return `!${e instanceof Error ? e.message : String(e)}`;
  }
}

function chunks(flat: bigint[], k: number): bigint[][] {
  const out: bigint[][] = [];
  for (let i = 0; i + k <= flat.length; i += k) out.push(flat.slice(i, i + k));
  return out;
}

function forms(flat: bigint[]): BinaryQF[] {
  return chunks(flat, 3).map((t) => new BinaryQF(t[0]!, t[1]!, t[2]!));
}

const lines = (items: string[]): string => items.join('\n');

function abs(n: bigint): bigint {
  return n < 0n ? -n : n;
}

/** Exact integer square root (floor). */
function isqrt(n: bigint): bigint {
  if (n < 0n) throw new RangeError('isqrt of negative');
  if (n < 2n) return n;
  let x = n;
  let y = (x + 1n) / 2n;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }
  return x;
}

function isSquare(n: bigint): boolean {
  if (n < 0n) return false;
  const s = isqrt(n);
  return s * s === n;
}

/** Python `//`: floor division rather than BigInt's truncating division. */
function floorDiv(a: bigint, b: bigint): bigint {
  const q = a / b;
  return a % b !== 0n && a < 0n !== b < 0n ? q - 1n : q;
}

/** Python `%`: remainder with the sign of the divisor. */
function floorMod(a: bigint, b: bigint): bigint {
  return a - floorDiv(a, b) * b;
}

/**
 * Whether PARI's `Qfb` constructor accepts this form.
 *
 * `reference/pari/src/basemath/Qfb.c:174-176`: a negative definite form raises
 * `pari_err_IMPL` and a square (hence also zero) discriminant raises
 * `pari_err_DOMAIN`.  Sage reaches that constructor through
 * `_pari_init_`/`__pari__` on both operands of `__mul__`.
 */
function pariQfbOk(f: BinaryQF): boolean {
  const D = f.discriminant();
  if (D < 0n) return f.a > 0n;
  return !isSquare(D);
}

/**
 * A version-independent name for the *class* of `f`.
 *
 * For an indefinite non-square form the reduced representative returned by
 * PARI depends on the PARI version (`reference/pari/CHANGES-2.16:142`), but
 * the proper cycle it lives in does not, so we sort the cycle.  `D == 0`
 * short-circuits to `singular` (SageMath 10.3 and the vendored 10.9 disagree
 * on the message raised there).
 */
function canonical(f: BinaryQF): string {
  const D = f.discriminant();
  if (D === 0n) return 'singular';
  if (D > 0n && !isSquare(D)) {
    return f.reduced_form().cycle({ proper: true }).map(F).sort().join(',');
  }
  return F(f.reduced_form());
}

// ---------------------------------------------------------------------------
// predicates / accessors
// ---------------------------------------------------------------------------

function qf_predicates(flat: bigint[]): string {
  return lines(
    forms(flat).map((f) => {
      const flags = [
        T(f.is_primitive()),
        T(f.is_zero()),
        T(f.is_reducible()),
        T(f.is_positive_definite()),
        T(f.is_negative_definite()),
        T(f.is_indefinite()),
        T(f.is_singular()),
        T(f.is_nonsingular()),
      ];
      return `${F(f)} D=${f.discriminant()} content=${f.content()} ${flags.join(' ')}`;
    })
  );
}

function qf_is_reduced(flat: bigint[]): string {
  return lines(forms(flat).map((f) => `${F(f)} ${attempt(() => T(f.is_reduced()))}`));
}

function qf_evaluate(flat: bigint[]): string {
  return lines(
    chunks(flat, 5).map((t) => {
      const f = new BinaryQF(t[0]!, t[1]!, t[2]!);
      return `${F(f)}(${t[3]},${t[4]})=${f.evaluate(t[3]!, t[4]!)}`;
    })
  );
}

function qf_matrix_action(flat: bigint[]): string {
  return lines(
    chunks(flat, 7).map((t) => {
      const f = new BinaryQF(t[0]!, t[1]!, t[2]!);
      const M: Matrix2 = [
        [t[3]!, t[4]!],
        [t[5]!, t[6]!],
      ];
      return (
        `${F(f)} [${t[3]},${t[4]};${t[5]},${t[6]}] ` +
        `right=${attempt(() => F(f.matrix_action_right(M)))} ` +
        `left=${attempt(() => F(f.matrix_action_left(M)))}`
      );
    })
  );
}

function qf_principal(Ds: bigint[]): string {
  return lines(Ds.map((D) => `${D} ${attempt(() => F(BinaryQF.principal(D)))}`));
}

// ---------------------------------------------------------------------------
// reduction
// ---------------------------------------------------------------------------

const ALGORITHMS = ['default', 'pari', 'sage'] as const;

function qf_reduced_form(flat: bigint[]): string {
  const out: string[] = [];
  for (const f of forms(flat)) {
    for (const alg of ALGORITHMS) {
      out.push(`${F(f)} ${alg} -> ${attempt(() => F(f.reduced_form({ algorithm: alg })))}`);
    }
  }
  return lines(out);
}

function qf_reduced_transformation(flat: bigint[]): string {
  const out: string[] = [];
  for (const f of forms(flat)) {
    for (const alg of ALGORITHMS) {
      out.push(
        `${F(f)} ${alg} -> ` +
          attempt(() => {
            const [g, M] = f.reduced_form({ transformation: true, algorithm: alg });
            const ok = T(F(f.matrix_action_right(M)) === F(g));
            const det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
            return `${F(g)} [${M[0][0]},${M[0][1]};${M[1][0]},${M[1][1]}] det=${det} action=${ok}`;
          })
      );
    }
  }
  return lines(out);
}

/**
 * `reduced_form(transformation=True)` under the default algorithm only.
 *
 * Used for large indefinite discriminants: `algorithm='sage'` runs Sage's
 * `_reduce_indef` loop, whose termination test is the 53-bit floating-point
 * `is_reduced` in SageMath 10.3 (fixed in the vendored 10.9, Sage issue
 * #37635), so it does not terminate there.  The default algorithm is PARI's
 * `qfbredsl2`, which is exact in both.
 */
function qf_reduced_transformation_default(flat: bigint[]): string {
  return lines(
    forms(flat).map(
      (f) =>
        `${F(f)} default -> ` +
        attempt(() => {
          const [g, M] = f.reduced_form({ transformation: true });
          const ok = T(F(f.matrix_action_right(M)) === F(g));
          const det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
          return `${F(g)} [${M[0][0]},${M[0][1]};${M[1][0]},${M[1][1]}] det=${det} action=${ok}`;
        })
    )
  );
}

function qf_reduced_class(flat: bigint[]): string {
  return lines(forms(flat).map((f) => `${F(f)} -> ${attempt(() => canonical(f))}`));
}

function qf_cycle(flat: bigint[]): string {
  const out: string[] = [];
  for (const f of forms(flat)) {
    for (const proper of [false, true]) {
      out.push(
        `${F(f)} proper=${T(proper)} -> ` + attempt(() => f.cycle({ proper }).map(F).join(','))
      );
    }
  }
  return lines(out);
}

// ---------------------------------------------------------------------------
// composition
// ---------------------------------------------------------------------------

function qf_compose(flat: bigint[]): string {
  return lines(
    chunks(flat, 6).map((t) => {
      const f = new BinaryQF(t[0]!, t[1]!, t[2]!);
      const g = new BinaryQF(t[3]!, t[4]!, t[5]!);
      return `${F(f)} * ${F(g)} = ${attempt(() => F(f.compose(g)))}`;
    })
  );
}

/**
 * Composition on inputs PARI's `Qfb` constructor rejects.
 *
 * Sage builds both operands with `Qfb(a, b, c)` before calling `qfbcompraw`,
 * so a negative definite operand or a square/zero discriminant raises a
 * `PariError` rather than returning a form.
 */
function qf_compose_domain(flat: bigint[]): string {
  return qf_compose(flat);
}

function qf_compose_table(D: bigint): string {
  const R = BinaryQF_reduced_representatives(D, { primitive_only: false, proper: true });
  const out = [
    `D=${D} n=${R.length} reps=${R.map(F).join(',')}`,
    `contents=${R.map((f) => f.content()).join(',')}`,
  ];
  for (const f of R) {
    for (const g of R) {
      out.push(
        `${F(f)} * ${F(g)} = ` +
          attempt(() => {
            const h = f.compose(g);
            return `${F(h)} ~ ${canonical(h)}`;
          })
      );
    }
  }
  return lines(out);
}

function qf_compose_powers(flat: bigint[]): string {
  return lines(
    chunks(flat, 4).map((t) => {
      const f = new BinaryQF(t[0]!, t[1]!, t[2]!);
      const k = Number(t[3]!);
      return (
        `${F(f)} -> ` +
        attempt(() => {
          let acc = f;
          const parts: string[] = [];
          for (let i = 1; i <= k; i++) {
            if (i > 1) acc = acc.compose(f);
            parts.push(`${i}:${canonical(acc)}`);
          }
          return parts.join(' ');
        })
      );
    })
  );
}

// ---------------------------------------------------------------------------
// class groups
// ---------------------------------------------------------------------------

function qf_reduced_representatives(Ds: bigint[]): string {
  const out: string[] = [];
  for (const D of Ds) {
    for (const primitive_only of [false, true]) {
      for (const proper of [true, false]) {
        out.push(
          `D=${D} primitive_only=${T(primitive_only)} proper=${T(proper)} -> ` +
            attempt(() => {
              const R = BinaryQF_reduced_representatives(D, { primitive_only, proper });
              return `${R.length} ${R.map(F).join(',')}`;
            })
        );
      }
    }
  }
  return lines(out);
}

/**
 * `class_number(D)`.
 *
 * Upstream has no such helper, so the oracle is the port's definition:
 * `len(BinaryQF_reduced_representatives(D, primitive_only=True, proper=True))`.
 */
function qf_class_number(Ds: bigint[]): string {
  return lines(Ds.map((D) => `${D} ` + attempt(() => String(class_number(D)))));
}

function qf_is_equivalent(flat: bigint[]): string {
  const out: string[] = [];
  for (const t of chunks(flat, 6)) {
    const f = new BinaryQF(t[0]!, t[1]!, t[2]!);
    const g = new BinaryQF(t[3]!, t[4]!, t[5]!);
    for (const proper of [true, false]) {
      out.push(
        `${F(f)} ~ ${F(g)} proper=${T(proper)} -> ` +
          attempt(() => T(f.is_equivalent(g, { proper })))
      );
    }
  }
  return lines(out);
}

function qf_class_group_closure(D: bigint): string {
  const R = BinaryQF_reduced_representatives(D, { primitive_only: false, proper: true });
  const P = BinaryQF.principal(D);
  const out = [`D=${D} principal=${F(P)} reps=${R.length}`];
  for (const f of R) {
    out.push(
      `${F(f)} ` +
        attempt(
          () =>
            `reduced=${T(f.is_reduced())} canon=${canonical(f)} ` +
            `id=${canonical(f.compose(P))} ` +
            `inv=${canonical(f.compose(new BinaryQF(f.a, -f.b, f.c)))}`
        )
    );
  }
  return lines(out);
}

// ---------------------------------------------------------------------------
// solve_integer
// ---------------------------------------------------------------------------

/**
 * All `(x, y)` with `a x^2 + b x y + c y^2 == n`, for `a > 0` and `D < 0`.
 *
 * `4 a n = (2 a x + b y)^2 - D y^2` bounds `|y|`, so the search is finite and
 * exhaustive.  Pure integer arithmetic, identical on both sides.
 */
function bruteForceSolutions(a: bigint, b: bigint, c: bigint, n: bigint): Array<[bigint, bigint]> {
  const D = b * b - 4n * a * c;
  const sols: Array<[bigint, bigint]> = [];
  if (n < 0n) return sols;
  const ymax = isqrt(floorDiv(4n * a * n, -D));
  for (let y = -ymax; y <= ymax; y++) {
    const disc = D * y * y + 4n * a * n;
    if (disc < 0n) continue;
    const s = isqrt(disc);
    if (s * s !== disc) continue;
    const roots = s === 0n ? [s] : [s, -s];
    for (const r of roots) {
      const num = -b * y + r;
      if (floorMod(num, 2n * a) === 0n) sols.push([floorDiv(num, 2n * a), y]);
    }
  }
  const seen = new Set<string>();
  const uniq = sols.filter((p) => {
    const k = `${p[0]},${p[1]}`;
    if (seen.has(k)) return false;
    seen.add(k);
    return true;
  });
  uniq.sort((p, q) =>
    p[0] !== q[0] ? (p[0] < q[0] ? -1 : 1) : p[1] < q[1] ? -1 : p[1] > q[1] ? 1 : 0
  );
  return uniq;
}

function qf_solve_integer_definite(flat: bigint[]): string {
  return lines(
    chunks(flat, 4).map((t) => {
      const [a, b, c, n] = [t[0]!, t[1]!, t[2]!, t[3]!];
      const f = new BinaryQF(a, b, c);
      return (
        `${F(f)} n=${n} -> ` +
        attempt(() => {
          const sols = bruteForceSolutions(a, b, c, n);
          const xy = f.solve_integer(n);
          const found = xy !== null;
          const valid = found && f.evaluate(xy[0], xy[1]) === n;
          const member = found && sols.some((p) => p[0] === xy[0] && p[1] === xy[1]);
          const agree = found === sols.length > 0 && (!found || member);
          const first = sols.length ? `(${sols[0]![0]},${sols[0]![1]})` : 'None';
          return `nsols=${sols.length} min=${first} found=${T(found)} valid=${T(valid)} agree=${T(agree)}`;
        })
      );
    })
  );
}

function qf_solve_integer_existence(flat: bigint[]): string {
  return lines(
    chunks(flat, 4).map((t) => {
      const f = new BinaryQF(t[0]!, t[1]!, t[2]!);
      const n = t[3]!;
      return (
        `${F(f)} n=${n} -> ` +
        attempt(() => {
          const xy = f.solve_integer(n);
          if (xy === null) return 'found=False';
          return `found=True valid=${T(f.evaluate(xy[0], xy[1]) === n)}`;
        })
      );
    })
  );
}

function qf_solve_integer_square_disc(flat: bigint[]): string {
  return lines(
    chunks(flat, 4).map((t) => {
      const f = new BinaryQF(t[0]!, t[1]!, t[2]!);
      const n = t[3]!;
      return (
        `${F(f)} n=${n} -> ` +
        attempt(() => {
          const xy = f.solve_integer(n);
          if (xy === null) return 'None';
          return `(${xy[0]},${xy[1]}) check=${T(f.evaluate(xy[0], xy[1]) === n)}`;
        })
      );
    })
  );
}

function qf_solve_integer_cornacchia(flat: bigint[]): string {
  return lines(
    chunks(flat, 4).map((t) => {
      const f = new BinaryQF(t[0]!, t[1]!, t[2]!);
      const n = t[3]!;
      return (
        `${F(f)} n=${n} -> ` +
        attempt(() => {
          const xy = f.solve_integer(n, { algorithm: 'cornacchia' });
          if (xy === null) return 'None';
          return `(${xy[0]},${xy[1]}) check=${T(f.evaluate(xy[0], xy[1]) === n)}`;
        })
      );
    })
  );
}

function qf_solve_integer_domain(flat: bigint[]): string {
  return lines(
    chunks(flat, 4).map((t) => {
      const f = new BinaryQF(t[0]!, t[1]!, t[2]!);
      const n = t[3]!;
      return (
        `${F(f)} n=${n} -> ` +
        attempt(() => {
          const xy = f.solve_integer(n);
          return xy === null ? 'None' : `(${xy[0]},${xy[1]})`;
        })
      );
    })
  );
}

// ---------------------------------------------------------------------------
// structured / random batteries
// ---------------------------------------------------------------------------

/**
 * Every `(a, b, c)` of discriminant `D` with `|a| <= bound` and
 * `|b| <= 2*bound` (including `a = 0` when `D` is a square).
 *
 * Deliberately *structured*: it produces negative-`a`, `a = 0`, `|a| > |c|`
 * and imprimitive shapes that uniform random sampling of small triples almost
 * never hits.
 */
function formsOfDiscriminant(D: bigint, bound: bigint): Triple[] {
  const out: Triple[] = [];
  for (let a = -bound; a <= bound; a++) {
    for (let b = -2n * bound; b <= 2n * bound; b++) {
      const num = b * b - D;
      if (a === 0n) {
        if (num === 0n) {
          for (let c = -bound; c <= bound; c++) out.push([0n, b, c]);
        }
        continue;
      }
      const den = 4n * a;
      if (floorMod(num, den) === 0n) out.push([a, b, floorDiv(num, den)]);
    }
  }
  return out;
}

function qf_sweep_discriminant(args: bigint[]): string {
  const D = args[0]!;
  const bound = args[1]!;
  const out: string[] = [];
  for (const t of formsOfDiscriminant(D, bound)) {
    const f = new BinaryQF(t[0], t[1], t[2]);
    const parts = [
      `content=${f.content()}`,
      `reduced=${attempt(() => T(f.is_reduced()))}`,
      `canon=${attempt(() => canonical(f))}`,
      `trans=${attempt(() => {
        const [g, M] = f.reduced_form({ transformation: true });
        return `${F(g)}|${M[0][0]},${M[0][1]},${M[1][0]},${M[1][1]}`;
      })}`,
      `sq=${pariQfbOk(f) ? attempt(() => F(f.compose(f))) : 'domain'}`,
    ];
    out.push(`${F(f)} ${parts.join(' ')}`);
  }
  return lines(out);
}

function qf_random_battery(values: bigint[]): string {
  return lines(
    chunks(values, 3).map((t) => {
      const f = new BinaryQF(t[0]!, t[1]!, t[2]!);
      const D = f.discriminant();
      const ok = pariQfbOk(f);
      const parts = [
        `D=${D}`,
        `content=${f.content()}`,
        `primitive=${T(f.is_primitive())}`,
        `reducible=${T(f.is_reducible())}`,
        `posdef=${T(f.is_positive_definite())}`,
        `canon=${attempt(() => canonical(f))}`,
        `sq=${ok ? attempt(() => F(f.compose(f))) : 'domain'}`,
        `sqcanon=${ok ? attempt(() => canonical(f.compose(f))) : 'domain'}`,
        `equiv=${D !== 0n ? attempt(() => T(f.is_equivalent(f))) : 'singular'}`,
      ];
      return `${F(f)} ${parts.join(' ')}`;
    })
  );
}

/**
 * Random walks in the class group of a random negative discriminant.
 *
 * For each random value `v` the discriminant `D = -(|v| + 3)` is rounded down
 * to `0` or `1` mod 4, its reduced representatives are enumerated, and a
 * deterministic product of representatives is reduced.  Exercises composition
 * on non-fundamental and imprimitive inputs.
 */
function qf_random_class_group(values: bigint[]): string {
  return lines(
    values.map((v) => {
      let D = -(abs(v) + 3n);
      while (floorMod(D, 4n) !== 0n && floorMod(D, 4n) !== 1n) D -= 1n;
      return (
        `D=${D} ` +
        attempt(() => {
          const R = BinaryQF_reduced_representatives(D, {
            primitive_only: false,
            proper: true,
          });
          if (R.length === 0) return 'empty';
          let acc = R[0]!;
          const names: string[] = [];
          for (let i = 0; i < R.length; i++) {
            const idx = Number(floorMod(abs(v) + BigInt(i), BigInt(R.length)));
            acc = acc.compose(R[idx]!);
            names.push(canonical(acc));
          }
          return `h=${R.length} ${names.join(' ')}`;
        })
      );
    })
  );
}

export const functions = {
  qf_predicates,
  qf_is_reduced,
  qf_evaluate,
  qf_matrix_action,
  qf_principal,
  qf_reduced_form,
  qf_reduced_transformation,
  qf_reduced_transformation_default,
  qf_reduced_class,
  qf_cycle,
  qf_compose,
  qf_compose_domain,
  qf_compose_table,
  qf_compose_powers,
  qf_reduced_representatives,
  qf_class_number,
  qf_is_equivalent,
  qf_class_group_closure,
  qf_solve_integer_definite,
  qf_solve_integer_existence,
  qf_solve_integer_square_disc,
  qf_solve_integer_cornacchia,
  qf_solve_integer_domain,
  qf_sweep_discriminant,
  qf_random_battery,
  qf_random_class_group,
};
