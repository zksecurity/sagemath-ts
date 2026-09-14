import { GF2 } from '../../../../packages/sagemath-ts/src/rings/finite_rings/gf2.js';
import { create_RealNumber } from '../../../../packages/sagemath-ts/src/rings/real_mpfr.js';
/**
 * sagemath-ts side of the `matrix_ops` property-test area.
 *
 * Cases: tests/property/cases/matrix_ops.cases.json
 * SageMath counterpart: tests/property/python/areas/matrix_ops.py
 */

import { IntegerMatrix } from '../../../../packages/sagemath-ts/src/index.js';
import {
  LLL,
  elementary_divisors_integer,
  hermite_normal_form,
  rank_integer,
  smith_form_integer,
} from '../../../../packages/sagemath-ts/src/matrix/matrix_integer.js';

export const functions = {
  determinant_2x2: (a: bigint, b: bigint, c: bigint, d: bigint) => {
    const M = new IntegerMatrix(2, 2, [
      [a, b],
      [c, d],
    ]);
    return M.determinant().value;
  },
  determinant_3x3: (
    a11: bigint,
    a12: bigint,
    a13: bigint,
    a21: bigint,
    a22: bigint,
    a23: bigint,
    a31: bigint,
    a32: bigint,
    a33: bigint
  ) => {
    const M = new IntegerMatrix(3, 3, [
      [a11, a12, a13],
      [a21, a22, a23],
      [a31, a32, a33],
    ]);
    return M.determinant().value;
  },
  determinant_4x4: (...args: bigint[]) => {
    const entries: bigint[][] = [];
    for (let i = 0; i < 4; i++) {
      entries.push(args.slice(i * 4, (i + 1) * 4));
    }
    const M = new IntegerMatrix(4, 4, entries);
    return M.determinant().value;
  },
  rank_2x3: (a11: bigint, a12: bigint, a13: bigint, a21: bigint, a22: bigint, a23: bigint) => {
    const M = new IntegerMatrix(2, 3, [
      [a11, a12, a13],
      [a21, a22, a23],
    ]);
    return BigInt(rank_integer(M));
  },
  rank_3x3: (
    a11: bigint,
    a12: bigint,
    a13: bigint,
    a21: bigint,
    a22: bigint,
    a23: bigint,
    a31: bigint,
    a32: bigint,
    a33: bigint
  ) => {
    const M = new IntegerMatrix(3, 3, [
      [a11, a12, a13],
      [a21, a22, a23],
      [a31, a32, a33],
    ]);
    return BigInt(rank_integer(M));
  },
  hnf_2x2: (a: bigint, b: bigint, c: bigint, d: bigint) => {
    const M = new IntegerMatrix(2, 2, [
      [a, b],
      [c, d],
    ]);
    const H = hermite_normal_form(M);
    return [
      [H.get(0, 0).value, H.get(0, 1).value],
      [H.get(1, 0).value, H.get(1, 1).value],
    ];
  },
  hnf_3x3: (...args: bigint[]) => {
    const entries: bigint[][] = [];
    for (let i = 0; i < 3; i++) {
      entries.push(args.slice(i * 3, (i + 1) * 3));
    }
    const M = new IntegerMatrix(3, 3, entries);
    const H = hermite_normal_form(M);
    const result: bigint[][] = [];
    for (let i = 0; i < 3; i++) {
      result.push([H.get(i, 0).value, H.get(i, 1).value, H.get(i, 2).value]);
    }
    return result;
  },
  snf_2x2: (a: bigint, b: bigint, c: bigint, d: bigint) => {
    const M = new IntegerMatrix(2, 2, [
      [a, b],
      [c, d],
    ]);
    const result = smith_form_integer(M, true);
    const D = Array.isArray(result) ? result[0] : result;
    return [D.get(0, 0).value, D.get(1, 1).value];
  },
  snf_3x3: (...args: bigint[]) => {
    const entries: bigint[][] = [];
    for (let i = 0; i < 3; i++) {
      entries.push(args.slice(i * 3, (i + 1) * 3));
    }
    const M = new IntegerMatrix(3, 3, entries);
    const result = smith_form_integer(M, true);
    const D = Array.isArray(result) ? result[0] : result;
    return [D.get(0, 0).value, D.get(1, 1).value, D.get(2, 2).value];
  },
  lll_2x2: (a: bigint, b: bigint, c: bigint, d: bigint) => {
    const M = new IntegerMatrix(2, 2, [
      [a, b],
      [c, d],
    ]);
    const L = LLL(M);
    return [
      [L.get(0, 0).value, L.get(0, 1).value],
      [L.get(1, 0).value, L.get(1, 1).value],
    ];
  },
  lll_3x3: (...args: bigint[]) => {
    const entries: bigint[][] = [];
    for (let i = 0; i < 3; i++) {
      entries.push(args.slice(i * 3, (i + 1) * 3));
    }
    const M = new IntegerMatrix(3, 3, entries);
    const L = LLL(M);
    const result: bigint[][] = [];
    for (let i = 0; i < 3; i++) {
      result.push([L.get(i, 0).value, L.get(i, 1).value, L.get(i, 2).value]);
    }
    return result;
  },
  elementary_divisors_2x2: (a: bigint, b: bigint, c: bigint, d: bigint) => {
    const M = new IntegerMatrix(2, 2, [
      [a, b],
      [c, d],
    ]);
    const divs = elementary_divisors_integer(M);
    return divs.map((d) => d.value);
  },
  elementary_divisors_3x3: (...args: bigint[]) => {
    const entries: bigint[][] = [];
    for (let i = 0; i < 3; i++) {
      entries.push(args.slice(i * 3, (i + 1) * 3));
    }
    const M = new IntegerMatrix(3, 3, entries);
    const divs = elementary_divisors_integer(M);
    return divs.map((d) => d.value);
  },
};

import {
  mzd_init,
  mzd_add,
  mzd_init_window,
  mzd_mul,
  mzd_mul_m4rm,
  mzd_mul_naive,
  mzd_make_table,
  type mzd_t,
} from '../../../../packages/m4ri-ts/src/index.js';
function m4riPackedFixture(m: number, n: number, seed: bigint): bigint[] {
  let state = seed;
  const rows: bigint[] = [];
  for (let i = 0; i < m; i++) {
    let row = 0n;
    for (let j = 0; j < n; j += 64) {
      state = BigInt.asUintN(64, state * 6364136223846793005n + 1442695040888963407n);
      row |= state << BigInt(j);
    }
    rows.push(BigInt.asUintN(n, row));
  }
  return rows;
}
Object.assign(functions, {
  m4ri_product: (
    method: bigint,
    m: bigint,
    k: bigint,
    n: bigint,
    seed: bigint,
    parameter: bigint,
    square: bigint
  ) => {
    const M = Number(m),
      K = Number(k),
      N = Number(n),
      p = Number(parameter),
      a = mzd_init(M, K, m4riPackedFixture(M, K, seed)),
      b = square
        ? a
        : mzd_init(
            method === 3n ? M : K,
            method === 3n ? K : N,
            m4riPackedFixture(method === 3n ? M : K, method === 3n ? K : N, seed + 1n)
          );
    const frame = (h: mzd_t) => [String(h.nrows), String(h.ncols), h.rows.map(String)];
    if (method === 5n) {
      const table = mzd_make_table(a, 0, p);
      return JSON.stringify([frame(mzd_init(table.rows.length, K, table.rows)), table.lookup]);
    }
    const h =
      method === 0n
        ? mzd_mul(a, b, p)
        : method === 1n
          ? mzd_mul_m4rm(a, b, p)
          : method === 2n
            ? mzd_mul_naive(a, b)
            : method === 3n
              ? mzd_add(a, b)
              : mzd_init_window(a, p, 64, M, K);
    return JSON.stringify(frame(h));
  },
});
Object.assign(functions, {
  m4ri_strassen_boundary: (
    m: bigint,
    k: bigint,
    n: bigint,
    seed: bigint,
    parameter: bigint,
    square: bigint
  ) => {
    const M = Number(m),
      K = Number(k),
      N = Number(n),
      a = mzd_init(M, K, m4riPackedFixture(M, K, seed)),
      b = square ? a : mzd_init(K, N, m4riPackedFixture(K, N, seed + 1n)),
      h = mzd_mul(a, b, Number(parameter));
    return JSON.stringify([String(h.nrows), String(h.ncols), h.rows.map(String)]);
  },
});

import { Matrix_mod2_dense } from '../../../../packages/sagemath-ts/src/matrix/matrix_mod2.js';
Object.assign(functions, {
  binary_matrix_product: (
    method: bigint,
    m: bigint,
    k: bigint,
    rightRows: bigint,
    n: bigint,
    seed: bigint,
    parameter: bigint,
    square: bigint
  ) => {
    const M = Number(m),
      K = Number(k),
      R = Number(rightRows),
      N = Number(n),
      aa = m4riPackedFixture(M, K, seed),
      bb = square ? aa : m4riPackedFixture(R, N, seed + 1n);
    const entries = (a: bigint[], width: number) =>
      a.map((row) => Array.from({ length: width }, (_, j) => Number((row >> BigInt(j)) & 1n)));
    const A = new Matrix_mod2_dense(M, K, entries(aa, K)),
      B = square ? A : new Matrix_mod2_dense(R, N, entries(bb, N));
    if (method === 1n && M && K && N && K !== R)
      return JSON.stringify('invalid-native-classical-dimensions');
    const C =
      method === 0n
        ? A.mul(B)
        : method === 1n
          ? A._multiply_classical(B)
          : method === 2n
            ? A._multiply_m4rm(B, parameter)
            : A._multiply_strassen(B, parameter);
    const values = Array.from({ length: C.nrows }, (_, i) =>
      String(
        Array.from({ length: C.ncols }, (_, j) => BigInt(C.get(i, j)) << BigInt(j)).reduce(
          (s, c) => s + c,
          0n
        )
      )
    );
    return JSON.stringify([String(C.nrows), String(C.ncols), values, C === A, C === B]);
  },
});

Object.assign(functions, {
  binary_matrix_access: (action: bigint, m: bigint, n: bigint, seed: bigint, args: bigint[]) => {
    const M = Number(m),
      N = Number(n),
      a = args.map(Number);
    const A = new Matrix_mod2_dense(
      M,
      N,
      Array.from({ length: M }, (_, i) =>
        Array.from({ length: N }, (_, j) => Number((BigInt(i * N + j) + seed) % 2n))
      )
    );
    let result: unknown;
    if (action === 0n) result = A.row(a[0]!, Boolean(a[1]));
    else if (action === 1n) {
      const B = A.submatrix(...a);
      result = [String(B.nrows), String(B.ncols), B.list().map(String)];
    } else if (action === 2n) result = String(A.get(a[0]!, a[1]!));
    else {
      A.set(a[0]!, a[1]!, a[2]!);
      result = A.list().map(String);
    }
    return JSON.stringify(result);
  },
});

import { Integer } from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';
Object.assign(functions, {
  binary_matrix_index_conversion: (
    action: bigint,
    m: bigint,
    n: bigint,
    kind: bigint,
    position: bigint,
    value: bigint
  ) => {
    const x: any =
      kind === 0n
        ? Number(value)
        : kind === 1n
          ? value
          : kind === 2n
            ? new Integer(value)
            : kind === 3n
              ? Boolean(value)
              : kind === 4n
                ? Number(value) + 0.5
                : kind === 5n
                  ? NaN
                  : kind === 6n
                    ? Infinity
                    : kind === 7n
                      ? String(value)
                      : kind === 8n
                        ? null
                        : new Rational(value, 2n);
    const M = Number(m),
      N = Number(n),
      p = Number(position),
      A = new Matrix_mod2_dense(
        M,
        N,
        Array.from({ length: M }, (_, i) => Array.from({ length: N }, (_, j) => (i * N + j) % 2))
      );
    let result: unknown;
    if (action === 0n) result = A.row(x);
    else if (action === 1n) {
      const args: any[] = [0, 0];
      args[p] = x;
      result = String(A.get(args[0], args[1]));
    } else if (action === 2n) {
      const args: any[] = [0, 0];
      args[p] = x;
      A.set(args[0], args[1], 1);
      result = A.list().map(String);
    } else {
      const args: any[] = [0, 0, 0, 0];
      args[p] = x;
      const B = A.submatrix(...args);
      result = [String(B.nrows), String(B.ncols), B.list().map(String)];
    }
    return JSON.stringify(result);
  },
});

import { mzd_submatrix } from '../../../../packages/m4ri-ts/src/index.js';
Object.assign(functions, {
  binary_matrix_slice_words: (
    m: bigint,
    n: bigint,
    seed: bigint,
    r0: bigint,
    c0: bigint,
    r1: bigint,
    c1: bigint
  ) => {
    const M = Number(m),
      N = Number(n),
      rows = m4riPackedFixture(M, N, seed),
      A = new Matrix_mod2_dense(
        M,
        N,
        rows.map((row) => Array.from({ length: N }, (_, j) => Number((row >> BigInt(j)) & 1n)))
      );
    const B = A.submatrix(r0, c0, r1 - r0, c1 - c0),
      p = mzd_submatrix(mzd_init(M, N, rows), Number(r0), Number(c0), Number(r1), Number(c1));
    const values = Array.from({ length: B.nrows }, (_, i) =>
      String(B.row(i).reduce((v, b, j) => v + (BigInt(b) << BigInt(j)), 0n))
    );
    return JSON.stringify([
      [String(B.nrows), String(B.ncols), values],
      [String(p.nrows), String(p.ncols), p.rows.map(String)],
      B === A,
    ]);
  },
});

import { mzd_density } from '../../../../packages/m4ri-ts/src/index.js';
import { realLiteralFrame } from './real_literals.js';
Object.assign(functions, {
  binary_matrix_density: (
    method: bigint,
    m: bigint,
    n: bigint,
    pattern: bigint,
    seed: bigint,
    resolution: bigint
  ) => {
    const M = Number(m),
      N = Number(n),
      rows = m4riPackedFixture(M, N, seed);
    const entries = Array.from({ length: M }, (_, i) =>
      Array.from({ length: N }, (_, j) =>
        pattern === 5n
          ? Number((rows[i]! >> BigInt(j)) & 1n)
          : pattern === 0n
            ? 0
            : pattern === 1n
              ? 1
              : pattern === 2n
                ? Number(j >= Math.floor(N / 2))
                : pattern === 3n
                  ? Number(j < Math.floor(N / 2))
                  : pattern === 4n
                    ? (i * N + j) % 2
                    : Number(BigInt(i * N + j) < seed)
      )
    );
    if (method === 2n) {
      const packed = entries.map((row) => (N ? BigInt('0b' + row.slice().reverse().join('')) : 0n)),
        value = mzd_density(mzd_init(M, N, packed), Number(resolution));
      if (Number.isNaN(value)) return JSON.stringify('NaN');
      const buffer = new ArrayBuffer(8),
        view = new DataView(buffer);
      view.setFloat64(0, value, false);
      return JSON.stringify(view.getBigUint64(0, false).toString(16).padStart(16, '0'));
    }
    const A = new Matrix_mod2_dense(M, N, entries);
    if (method === 0n) {
      const v = A.density();
      return JSON.stringify([
        typeof v === 'bigint' ? 'int' : typeof v === 'number' ? 'number' : v.constructor.name,
        String(v),
      ]);
    }
    return JSON.stringify(realLiteralFrame(A.density(true) as any));
  },
});

const BINARY_ALGORITHMS = [
  'heuristic',
  'm4ri',
  'pluq',
  'classical',
  'ple',
  'linbox',
  'default',
  'bogus',
];
functions.binary_matrix_elimination = (
  method: bigint,
  m: bigint,
  n: bigint,
  seed: bigint,
  algorithm: bigint,
  reduced: bigint,
  sequence: bigint
) => {
  const M = Number(m),
    N = Number(n),
    rows = m4riPackedFixture(M, N, seed);
  const A = new Matrix_mod2_dense(
    M,
    N,
    rows.map((row) => Array.from({ length: N }, (_, j) => Number((row >> BigInt(j)) & 1n)))
  );
  const algo = BINARY_ALGORITHMS[Number(algorithm)]!;
  const frame = (B: Matrix_mod2_dense) => B.list().map(String);
  const observe = (call: () => unknown) => {
    try {
      return ['ok', call()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const call = () =>
    method === 0n
      ? A.rank(algo)
      : method === 1n
        ? (A.echelonize(algo, 0, Boolean(reduced)), frame(A))
        : frame(A.echelon_form(algo, 0, Boolean(reduced)));
  if (sequence === 1n) A.rank();
  else if (sequence === 2n) A.echelonize('m4ri', 0, false);
  else if (sequence === 3n) A.echelon_form('m4ri', 0, false);
  const first = observe(call);
  return JSON.stringify([
    first,
    frame(A),
    observe(() => A.rank('bogus')),
    observe(() => {
      A.echelonize('bogus');
      return null;
    }),
    frame(A),
  ]);
};

import {
  mzd_echelonize_m4ri,
  mzd_echelonize,
  mzd_echelonize_pluq,
  mzd_ple,
} from '../../../../packages/m4ri-ts/src/index.js';
functions.m4ri_echelon = (
  m: bigint,
  n: bigint,
  seed: bigint,
  algorithm: bigint,
  full: bigint,
  k: bigint
) => {
  const M = Number(m),
    N = Number(n),
    A = mzd_init(M, N, m4riPackedFixture(M, N, seed));
  const result =
    algorithm === 0n
      ? mzd_echelonize_m4ri(A, Boolean(full), Number(k))
      : algorithm === 1n
        ? mzd_echelonize(A, Boolean(full))
        : mzd_echelonize_pluq(A, Boolean(full));
  return JSON.stringify([result.rank, result.matrix.rows.map(String)]);
};

import { _mzd_ple_russian } from '../../../../packages/m4ri-ts/src/ple_russian.js';
functions.m4ri_ple = (m: bigint, n: bigint, seed: bigint, kind: bigint, k: bigint) => {
  const M = Number(m),
    N = Number(n),
    A = mzd_init(M, N, m4riPackedFixture(M, N, seed));
  const result = kind === 0n ? _mzd_ple_russian(A, Number(k)) : mzd_ple(A, Number(k));
  return JSON.stringify([result.rank, result.P, result.Q, result.matrix.rows.map(String)]);
};

import {
  mzd_trsm_upper_left,
  mzd_trsm_lower_left,
} from '../../../../packages/m4ri-ts/src/index.js';
functions.m4ri_trsm = (m: bigint, n: bigint, seed: bigint, upper: bigint, cutoff: bigint) => {
  const M = Number(m),
    N = Number(n),
    rows = m4riPackedFixture(M, M, seed).map((row, i) => {
      const diagonal = 1n << BigInt(i);
      return (row & (upper ? -1n << BigInt(i) : diagonal - 1n)) | diagonal;
    });
  const A = mzd_init(M, M, rows),
    B = mzd_init(M, N, m4riPackedFixture(M, N, seed + 1n));
  return JSON.stringify(
    (upper ? mzd_trsm_upper_left : mzd_trsm_lower_left)(A, B, Number(cutoff)).rows.map(String)
  );
};

functions.binary_matrix_cache = (
  m: bigint,
  n: bigint,
  seed: bigint,
  prepare: bigint,
  operation: bigint
) => {
  const M = Number(m),
    N = Number(n),
    packed = m4riPackedFixture(M, N, seed);
  let A = new Matrix_mod2_dense(
    M,
    N,
    packed.map((row) => Array.from({ length: N }, (_, j) => Number((row >> BigInt(j)) & 1n)))
  );
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.echelonize('m4ri', 0, false);
  else if (prepare === 3n) A.echelon_form('m4ri', 0, false);
  else if (prepare === 4n) A = A.echelon_form('m4ri', 0, false);
  else if (prepare === 5n) A.set_immutable();
  else if (prepare === 6n) A.echelonize('classical');
  const before = [A.is_immutable(), A.is_mutable()];
  const observe = (call: () => unknown) => {
    try {
      return ['ok', call()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const output = observe(() => {
    if (operation === 0n) A.set(0, 0, A.get(0, 0) ^ 1);
    else if (operation === 1n) A.set(M, N, 1);
    else if (operation === 2n) A.swap_rows(0, Math.max(0, M - 1));
    else if (operation === 3n) A.swap_rows(-1, M);
    else if (operation === 4n) A.swap_columns(0, Math.max(0, N - 1));
    else if (operation === 5n) A.swap_columns(-1, N);
    else if (operation === 6n)
      A.permute_rows(Array.from({ length: M }, (_, i) => (M ? 1 + ((i + 1) % M) : 0)));
    else if (operation === 7n)
      A.permute_columns(Array.from({ length: N }, (_, i) => (N ? 1 + ((i + 1) % N) : 0)));
    else if (operation === 8n) A.randomize(0);
    else if (operation === 9n) A.randomize(-1);
    else if (operation === 10n) A.echelonize('bogus');
    else if (operation === 11n) A.echelonize('m4ri');
    else if (operation === 12n) return A.echelon_form('bogus').list().map(String);
    else if (operation === 13n) A._clear_cache();
    else if (operation === 14n) A.set_immutable();
    else if (operation === 15n) A = A.copy();
    else if (operation === 16n) A.doubly_lexical_ordering(true);
    else if (operation === 17n) A.doubly_lexical_ordering(false);
    return null;
  });
  return JSON.stringify([
    before,
    output,
    A.list().map(String),
    [A.is_immutable(), A.is_mutable()],
    observe(() => A.rank('bogus')),
    observe(() => {
      A.echelonize('bogus');
      return null;
    }),
  ]);
};
functions.binary_matrix_echelon_options = (
  m: bigint,
  n: bigint,
  seed: bigint,
  algorithm: bigint,
  reduced: bigint,
  k: bigint,
  prepare: bigint
) => {
  const M = Number(m),
    N = Number(n),
    packed = m4riPackedFixture(M, N, seed),
    A = new Matrix_mod2_dense(
      M,
      N,
      packed.map((row) => Array.from({ length: N }, (_, j) => Number((row >> BigInt(j)) & 1n)))
    );
  if (prepare === 1n) A.echelonize();
  else if (prepare === 2n) A.set_immutable();
  else if (prepare === 3n) A.rank();
  let returned: unknown;
  const observe = (call: () => unknown) => {
    try {
      return ['ok', call()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const result = observe(() => {
    returned = A.echelonize(BINARY_ALGORITHMS[Number(algorithm)], 0, Boolean(reduced), { k });
    return returned === A ? 'self' : 'None';
  });
  return JSON.stringify([result, A.list().map(String), observe(() => A.rank('bogus'))]);
};

functions.m4ri_echelon_pattern = (
  m: bigint,
  n: bigint,
  seed: bigint,
  pattern: bigint,
  algorithm: bigint,
  full: bigint,
  k: bigint
) => {
  const M = Number(m),
    N = Number(n),
    rows = m4riPackedFixture(M, N, seed).map((row, i) => {
      if (pattern === 0n) return row;
      if (pattern === 1n) return 0n;
      if (pattern === 2n) return N ? 1n << BigInt(N - 1) : 0n;
      if (pattern === 3n) return i < 400 ? (i < N ? 1n << BigInt(i) : 0n) : row;
      if (pattern === 4n) {
        const start = Math.max(0, N - 128);
        return (row & (-1n << BigInt(start))) | (i < start ? 1n << BigInt(i) : 0n);
      }
      if (pattern === 5n) return row & (-1n << 128n);
      return i >= Math.floor(M / 2) ? 0n : row;
    });
  const A = mzd_init(M, N, rows),
    result =
      algorithm === 0n
        ? mzd_echelonize_m4ri(A, Boolean(full), Number(k))
        : algorithm === 1n
          ? mzd_echelonize(A, Boolean(full))
          : mzd_echelonize_pluq(A, Boolean(full));
  return JSON.stringify([result.rank, result.matrix.rows.map(String)]);
};

import {
  ple as sagePle,
  pluq as sagePluq,
} from '../../../../packages/sagemath-ts/src/matrix/matrix_mod2.js';
functions.binary_factorization = (
  method: bigint,
  m: bigint,
  n: bigint,
  seed: bigint,
  algorithm: bigint,
  param: bigint,
  kind: bigint
) => {
  const M = Number(m),
    N = Number(n),
    rows = m4riPackedFixture(M, N, seed),
    A = new Matrix_mod2_dense(
      M,
      N,
      rows.map((row) => Array.from({ length: N }, (_, j) => Number((row >> BigInt(j)) & 1n)))
    );
  const name = ['standard', method ? 'mmpf' : 'russian', 'naive', 'bogus'][Number(algorithm)]!;
  const value = kind === 0n ? Number(param) : kind === 1n ? param : new Integer(param);
  const [B, P, Q] = (method ? sagePluq : sagePle)(A, name, value as any);
  return JSON.stringify([
    B.list().map(String),
    P,
    Q,
    A.list().map(String),
    B === A,
    B.is_mutable(),
  ]);
};

import {
  _mzd_ple_naive,
  _mzd_pluq_naive,
  _mzd_pluq_russian,
  mzd_pluq,
} from '../../../../packages/m4ri-ts/src/index.js';
functions.m4ri_factorization = (m: bigint, n: bigint, seed: bigint, kind: bigint, k: bigint) => {
  const M = Number(m),
    N = Number(n),
    A = mzd_init(M, N, m4riPackedFixture(M, N, seed));
  const result =
    kind === 0n
      ? mzd_ple(A, Number(k))
      : kind === 1n
        ? _mzd_ple_russian(A, Number(k))
        : kind === 2n
          ? _mzd_ple_naive(A)
          : kind === 3n
            ? mzd_pluq(A, Number(k))
            : kind === 4n
              ? _mzd_pluq_russian(A, Number(k))
              : _mzd_pluq_naive(A);
  return JSON.stringify([result.rank, result.P, result.Q, result.matrix.rows.map(String)]);
};

functions.binary_factorization_conversion = (
  method: bigint,
  algorithm: bigint,
  kind: bigint,
  value: bigint
) => {
  const param =
    kind === 0n
      ? Number(value) / 2
      : kind === 1n
        ? NaN
        : kind === 2n
          ? Infinity
          : kind === 3n
            ? -Infinity
            : kind === 4n
              ? String(value)
              : kind === 5n
                ? null
                : kind === 6n
                  ? Boolean(value)
                  : kind === 7n
                    ? new Rational(value, 2n)
                    : kind === 8n
                      ? new Integer(value)
                      : value;
  const A = new Matrix_mod2_dense(3, 5, [
    [0, 1, 0, 1, 0],
    [0, 1, 1, 0, 1],
    [0, 0, 1, 1, 1],
  ]);
  const [B, P, Q] = (method ? sagePluq : sagePle)(
    A,
    ['standard', method ? 'mmpf' : 'russian', 'naive', 'bogus'][Number(algorithm)],
    param as any
  );
  return JSON.stringify([B.list().map(String), P, Q]);
};

functions.binary_swap_conversion = (
  axis: bigint,
  m: bigint,
  n: bigint,
  prepare: bigint,
  kind: bigint,
  value: bigint,
  position: bigint
) => {
  const M = Number(m),
    N = Number(n),
    rows = m4riPackedFixture(M, N, 42n);
  let A = new Matrix_mod2_dense(
    M,
    N,
    rows.map((r) => Array.from({ length: N }, (_, j) => Number((r >> BigInt(j)) & 1n)))
  );
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.echelonize('m4ri');
  else if (prepare === 3n) A.set_immutable();
  else if (prepare === 4n) A = A.echelon_form('m4ri');
  const v =
    kind === 0n
      ? value
      : kind === 1n
        ? new Integer(value)
        : kind === 2n
          ? new Rational(value, 2n)
          : kind === 3n
            ? Number(value) / 2
            : kind === 4n
              ? NaN
              : kind === 5n
                ? Infinity
                : kind === 6n
                  ? String(value)
                  : kind === 7n
                    ? null
                    : Boolean(value);
  const args = position === 0n ? [v, 0] : position === 1n ? [0, v] : [v, v];
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const result = observe(() => {
    if (axis) A.swap_columns(args[0] as any, args[1] as any);
    else A.swap_rows(args[0] as any, args[1] as any);
    return null;
  });
  return JSON.stringify([
    result,
    A.list().map(String),
    observe(() => A.rank('bogus')),
    A.is_immutable(),
  ]);
};
functions.binary_inverse = (
  m: bigint,
  n: bigint,
  seed: bigint,
  pattern: bigint,
  prepare: bigint
) => {
  const M = Number(m),
    N = Number(n),
    rows = m4riPackedFixture(M, N, seed).map((r, i) =>
      pattern === 1n ? (r & (-1n << BigInt(i))) | (1n << BigInt(i)) : pattern === 2n ? 0n : r
    );
  const A = new Matrix_mod2_dense(
    M,
    N,
    rows.map((r) => Array.from({ length: N }, (_, j) => Number((r >> BigInt(j)) & 1n)))
  );
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.set_immutable();
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const result = observe(() => {
    const B = A.inverse();
    return [B.nrows, B.ncols, B.list().map(String), B === A, B.is_mutable()];
  });
  return JSON.stringify([result, A.list().map(String), observe(() => A.rank('bogus'))]);
};

import { mzd_inv_m4ri } from '../../../../packages/m4ri-ts/src/index.js';
functions.m4ri_inverse = (n: bigint, seed: bigint, pattern: bigint, k: bigint) => {
  const N = Number(n),
    rows = m4riPackedFixture(N, N, seed).map((r, i) =>
      pattern === 1n ? (r & (-1n << BigInt(i))) | (1n << BigInt(i)) : pattern === 2n ? 0n : r
    );
  return JSON.stringify(mzd_inv_m4ri(mzd_init(N, N, rows), Number(k)).rows.map(String));
};

functions.binary_permutation = (
  axis: bigint,
  m: bigint,
  n: bigint,
  prepare: bigint,
  perm: bigint[]
) => {
  const M = Number(m),
    N = Number(n),
    rows = m4riPackedFixture(M, N, 42n);
  let A = new Matrix_mod2_dense(
    M,
    N,
    rows.map((r) => Array.from({ length: N }, (_, j) => Number((r >> BigInt(j)) & 1n)))
  );
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.echelonize('m4ri');
  else if (prepare === 3n) A.set_immutable();
  else if (prepare === 4n) A = A.echelon_form('m4ri');
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const result = observe(() => {
    if (axis) A.permute_columns(perm.map(Number));
    else A.permute_rows(perm.map(Number));
    return null;
  });
  return JSON.stringify([
    result,
    A.nrows,
    A.ncols,
    A.list().map(String),
    observe(() => A.rank('bogus')),
  ]);
};

functions.binary_solve = (
  m: bigint,
  n: bigint,
  brows: bigint,
  bcols: bigint,
  seed: bigint,
  pattern: bigint,
  check: bigint,
  prepare: bigint
) => {
  const M = Number(m),
    N = Number(n),
    BR = Number(brows),
    BC = Number(bcols);
  const make = (r: number, c: number, rows: readonly bigint[]) =>
    new Matrix_mod2_dense(
      r,
      c,
      rows.map((row) => Array.from({ length: c }, (_, j) => Number((row >> BigInt(j)) & 1n)))
    );
  const rows = m4riPackedFixture(M, N, seed).map((r, i) =>
    pattern === 1n ? 0n : pattern === 2n ? (i < N ? 1n << BigInt(i) : 0n) : r
  );
  const A = make(M, N, rows);
  let B = make(BR, BC, m4riPackedFixture(BR, BC, seed + 1n));
  if (pattern === 3n && BR === M) B = A.mul(make(N, BC, m4riPackedFixture(N, BC, seed + 2n)));
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) {
    A.set_immutable();
    B.set_immutable();
  } else if (prepare === 3n) A.echelon_form();
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const result = observe(() => {
    const X = check === 2n ? A.solve_right(B) : A.solve_right(B, Boolean(check));
    return [X.nrows, X.ncols, X.list().map(String), X.is_mutable()];
  });
  return JSON.stringify([
    result,
    A.list().map(String),
    B.list().map(String),
    observe(() => A.rank('bogus')),
    observe(() => A.echelon_form('bogus').list().map(String)),
  ]);
};
functions.binary_kernel = (
  m: bigint,
  n: bigint,
  seed: bigint,
  pattern: bigint,
  basis: bigint,
  algorithm: bigint,
  prepare: bigint
) => {
  const M = Number(m),
    N = Number(n),
    rows = m4riPackedFixture(M, N, seed).map((r) =>
      pattern === 1n ? 0n : pattern === 2n ? (N ? 1n << BigInt(N - 1) : 0n) : r
    );
  const A = new Matrix_mod2_dense(
    M,
    N,
    rows.map((r) => Array.from({ length: N }, (_, j) => Number((r >> BigInt(j)) & 1n)))
  );
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.set_immutable();
  else if (prepare === 3n) A.echelon_form();
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const result = observe(() => {
    const B = A.right_kernel_matrix({
      basis: ['echelon', 'pivot', 'computed', 'default', 'bogus', 'LLL'][Number(basis)],
      algorithm: ['default', 'generic', 'pluq', 'bogus', 'flint', 'linbox', 'pari', 'padic'][
        Number(algorithm)
      ],
    } as any);
    return [B.nrows, B.ncols, B.list().map(String), B.is_mutable()];
  });
  return JSON.stringify([
    result,
    A.list().map(String),
    observe(() => A.rank('bogus')),
    observe(() => A.echelon_form('bogus').list().map(String)),
  ]);
};

import { mzd_solve_left, mzd_kernel_left_pluq } from '../../../../packages/m4ri-ts/src/index.js';
functions.m4ri_solve = (
  m: bigint,
  n: bigint,
  bc: bigint,
  seed: bigint,
  pattern: bigint,
  cutoff: bigint,
  check: bigint
) => {
  const M = Number(m),
    N = Number(n),
    C = Number(bc),
    R = Math.max(M, N);
  const a = mzd_init(
    M,
    N,
    m4riPackedFixture(M, N, seed).map((r, i) =>
      pattern === 1n ? 0n : pattern === 2n ? (i < N ? 1n << BigInt(i) : 0n) : r
    )
  );
  let b = mzd_init(R, C, m4riPackedFixture(R, C, seed + 1n));
  if (pattern === 3n) {
    const product = mzd_mul(a, mzd_init(N, C, m4riPackedFixture(N, C, seed + 2n)));
    b = mzd_init(R, C, [...product.rows, ...Array<bigint>(R - M).fill(0n)]);
  }
  const result = mzd_solve_left(a, b, Number(cutoff), Boolean(check));
  return JSON.stringify([
    result.status,
    result.matrix.rows.map(String),
    result.rhs.rows.map(String),
  ]);
};
functions.m4ri_kernel = (m: bigint, n: bigint, seed: bigint, pattern: bigint, cutoff: bigint) => {
  const M = Number(m),
    N = Number(n),
    a = mzd_init(
      M,
      N,
      m4riPackedFixture(M, N, seed).map((r, i) =>
        pattern === 1n ? 0n : pattern === 2n ? (i < N ? 1n << BigInt(i) : 0n) : r
      )
    );
  const result = mzd_kernel_left_pluq(a, Number(cutoff));
  return JSON.stringify([
    result.matrix.rows.map(String),
    result.kernel === null ? null : [result.kernel.ncols, result.kernel.rows.map(String)],
  ]);
};

functions.binary_kernel_proof = (
  m: bigint,
  n: bigint,
  basis: bigint,
  algorithm: bigint,
  proof: bigint
) => {
  const A = new Matrix_mod2_dense(Number(m), Number(n));
  const B = A.right_kernel_matrix({
    basis: ['default', 'computed', 'bogus', 'LLL'][Number(basis)],
    algorithm: ['default', 'bogus', 'linbox'][Number(algorithm)],
    proof: [null, true, false, 0n, 1n, 2n, 'bad'][Number(proof)],
  } as any);
  return JSON.stringify([B.nrows, B.ncols, B.list().map(String), B.is_mutable()]);
};

functions.binary_basic = (
  m: bigint,
  n: bigint,
  r: bigint,
  c: bigint,
  seed: bigint,
  operation: bigint,
  prepare: bigint
) => {
  const make = (m: bigint, n: bigint, seed: bigint) =>
    new Matrix_mod2_dense(
      Number(m),
      Number(n),
      m4riPackedFixture(Number(m), Number(n), seed).map((v) =>
        Array.from({ length: Number(n) }, (_, j) => Number((v >> BigInt(j)) & 1n))
      )
    );
  const A = make(m, n, seed),
    B = make(r, c, seed + 1n);
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.set_immutable();
  else if (prepare === 3n) A.echelon_form();
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const result = observe(() => {
    let X: Matrix_mod2_dense;
    if (operation === 0n) X = A.add(B);
    else if (operation === 1n) X = A.sub(B);
    else if (operation === 2n) X = A.transpose();
    else if (operation === 3n) X = A.augment(B);
    else if (operation === 4n) return A.determinant();
    else if (operation === 5n) X = A.neg();
    else if (operation === 6n) X = A.copy();
    else if (operation === 7n) return A.columns().map((row) => row.map(String));
    else if (operation === 8n) return A.list().map(String);
    else if (operation === 9n) return A.eq(B);
    else return A.toString();
    return [X.nrows, X.ncols, X.list().map(String), X.is_mutable()];
  });
  return JSON.stringify([
    result,
    A.list().map(String),
    observe(() => A.rank('bogus')),
    observe(() => A.echelon_form('bogus').list().map(String)),
  ]);
};

functions.binary_columns_state = (
  m: bigint,
  n: bigint,
  copy: bigint,
  operation: bigint,
  prepare: bigint
) => {
  const M = Number(m),
    N = Number(n),
    A = new Matrix_mod2_dense(
      M,
      N,
      m4riPackedFixture(M, N, 42n).map((r) =>
        Array.from({ length: N }, (_, j) => Number((r >> BigInt(j)) & 1n))
      )
    );
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.echelon_form();
  const C = copy === 2n ? A.columns() : A.columns(Boolean(copy)),
    D = A.columns(false);
  const before = [C === D, Boolean(C.length && D.length && C[0] === D[0])];
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const result = observe(() => {
    if (operation === 1n && C.length && M) C[0]![0] = 1 - C[0]![0]!;
    else if (operation === 2n && C.length) C[0] = Array<number>(M).fill(0);
    else if (operation === 3n) C.push(Array<number>(M).fill(1));
    else if (operation === 4n && C.length) C.pop();
    else if (operation === 5n && M && N) A.set(0, 0, 1 - A.get(0, 0));
    else if (operation === 6n) A.set(M, N, 1);
    else if (operation === 7n) A._clear_cache();
    return null;
  });
  const E = A.columns(false),
    F = A.columns();
  return JSON.stringify([
    before,
    result,
    [E === D, F === E, Boolean(E.length && F.length && E[0] === F[0])],
    C.map((row) => row.map(String)),
    E.map((row) => row.map(String)),
    A.list().map(String),
    observe(() => A.rank('bogus')),
  ]);
};

import { mzd_transpose, mzd_concat, mzd_add } from '../../../../packages/m4ri-ts/src/index.js';
functions.m4ri_basic = (m: bigint, n: bigint, right: bigint, seed: bigint, op: bigint) => {
  const M = Number(m),
    N = Number(n),
    R = Number(op === 1n ? n : right),
    A = mzd_init(M, N, m4riPackedFixture(M, N, seed)),
    B = mzd_init(M, R, m4riPackedFixture(M, R, seed + 1n));
  const C = op === 0n ? mzd_transpose(A) : op === 1n ? mzd_add(A, B) : mzd_concat(A, B);
  return JSON.stringify([C.nrows, C.ncols, C.rows.map(String)]);
};

functions.binary_format = (
  m: bigint,
  n: bigint,
  seed: bigint,
  mapping: bigint,
  zero: bigint,
  one: bigint,
  minus: bigint
) => {
  const M = Number(m),
    N = Number(n),
    A = new Matrix_mod2_dense(
      M,
      N,
      m4riPackedFixture(M, N, seed).map((r) =>
        Array.from({ length: N }, (_, j) => Number((r >> BigInt(j)) & 1n))
      )
    );
  const texts = [undefined, '', '.', 'zero', '1', 'one', '🙂', 'a\nb'];
  let counter = 0;
  const mp =
    mapping === 0n
      ? undefined
      : mapping === 1n
        ? {}
        : mapping === 2n
          ? { 0: 'xx' }
          : mapping === 3n
            ? { 1: 'y' }
            : mapping === 4n
              ? { 0: 'xx', 1: 'y' }
              : (x: number) => String((++counter + x) % 3);
  let result;
  try {
    result = [
      'ok',
      A.str(mp as any, texts[Number(zero)], texts[Number(one)], texts[Number(minus)]),
    ];
  } catch (e) {
    result = ['error', (e as Error).name, (e as Error).message];
  }
  const state = typeof mp === 'object' ? Object.entries(mp).sort() : counter;
  return JSON.stringify([result, state]);
};
functions.binary_subdivision = (
  m: bigint,
  n: bigint,
  rows: bigint[],
  cols: bigint[],
  operation: bigint,
  prepare: bigint
) => {
  const M = Number(m),
    N = Number(n);
  let A = new Matrix_mod2_dense(
    M,
    N,
    m4riPackedFixture(M, N, 42n).map((r) =>
      Array.from({ length: N }, (_, j) => Number((r >> BigInt(j)) & 1n))
    )
  );
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.echelon_form();
  else if (prepare === 3n) A.set_immutable();
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  let result = observe(() => {
    (A as any).subdivide(rows, cols);
    return null;
  });
  if (operation === 1n) A = A.copy();
  else if (operation === 2n) A = A.transpose();
  else if (operation === 3n) A = A.neg();
  else if (operation === 4n) A = A.add(A);
  else if (operation === 5n) A = A.augment(A, true);
  else if (operation === 6n) A = A.augment(A, false);
  else if (operation === 7n) A = A.echelon_form();
  else if (operation === 8n) {
    const [r, c] = (A as any).get_subdivisions();
    r.push(0);
    c.push(0);
  } else if (operation === 9n)
    result = observe(() => {
      (A as any).subdivide();
      return null;
    });
  else if (operation === 10n)
    result = observe(() => {
      (A as any).subdivide([1], undefined);
      return null;
    });
  return JSON.stringify([
    result,
    (A as any).subdivisions().map((r: bigint[]) => r.map(String)),
    A.is_mutable(),
    observe(() => A.str()),
    observe(() => A.str(undefined, '.')),
    observe(() => A.rank('bogus')),
  ]);
};

functions.binary_subdivision_conversion = (
  kind: bigint,
  value: bigint,
  side: bigint,
  prepare: bigint
) => {
  const A = new Matrix_mod2_dense(2, 3, [
    [1, 0, 1],
    [0, 1, 0],
  ]);
  A.subdivide(1, 1);
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.set_immutable();
  else if (prepare === 3n) A.echelon_form();
  const convert = (x: bigint): any =>
    kind === 0n
      ? Number(x) / 2
      : kind === 1n
        ? new Rational(x, 2n)
        : kind === 2n
          ? new Integer(x)
          : kind === 3n
            ? x
            : kind === 4n
              ? String(x)
              : kind === 5n
                ? null
                : kind === 6n
                  ? Boolean(x)
                  : kind === 7n
                    ? NaN
                    : Infinity;
  const v = convert(value);
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const result = observe(() => {
    if (side === 0n) A.subdivide(v, null);
    else if (side === 1n) A.subdivide(null, v);
    else if (side === 2n) A.subdivide([v, convert(0n)], []);
    else if (side === 3n) A.subdivide([], [v, convert(2n)]);
    else if (side === 4n) A.subdivide([v, convert(2n)], ['bad', 0]);
    else A.subdivide([v, 0], []);
    return null;
  });
  return JSON.stringify([
    result,
    A.subdivisions().map((r) => r.map(String)),
    observe(() => A.rank('bogus')),
  ]);
};

functions.binary_format_invalid = (m: bigint, n: bigint, kind: bigint) => {
  const M = Number(m),
    N = Number(n),
    A = new Matrix_mod2_dense(
      M,
      N,
      Array.from({ length: M }, (_, i) => Array.from({ length: N }, (_, j) => (i + j) % 2))
    );
  const mp = [
    { 0: null },
    { 1: 1 },
    () => null,
    () => 1,
    () => ['x'],
    undefined,
    undefined,
    1,
    '0',
    [],
    [0],
    [1],
    { 0: { x: 1 } },
  ][Number(kind)];
  return A.str(
    mp as any,
    (kind === 5n ? 1 : undefined) as any,
    (kind === 6n ? false : undefined) as any
  );
};
functions.binary_augmentation_subdivisions = (
  m: bigint,
  n: bigint,
  c: bigint,
  lr: bigint[],
  lc: bigint[],
  rr: bigint[],
  rc: bigint[],
  flag: bigint,
  operation: bigint
) => {
  const make = (m: bigint, n: bigint, seed: bigint) =>
    new Matrix_mod2_dense(
      Number(m),
      Number(n),
      m4riPackedFixture(Number(m), Number(n), seed).map((r) =>
        Array.from({ length: Number(n) }, (_, j) => Number((r >> BigInt(j)) & 1n))
      )
    );
  const A = make(m, n, 42n),
    B = make(m, c, 43n);
  A.subdivide(lr, lc);
  B.subdivide(rr, rc);
  A.set_immutable();
  B.set_immutable();
  let C = A.augment(B, Boolean(flag));
  if (operation === 1n) C = C.copy();
  else if (operation === 2n) C = C.transpose();
  else if (operation === 3n) C = C.neg();
  else if (operation === 4n) C = C.submatrix(0, 0, C.nrows, C.ncols);
  return JSON.stringify([
    C.nrows,
    C.ncols,
    C.subdivisions().map((r) => r.map(String)),
    C.str(),
    C.str(undefined, '.'),
    C.is_mutable(),
    A.subdivisions().map((r) => r.map(String)),
    B.subdivisions().map((r) => r.map(String)),
  ]);
};

functions.binary_row_state = (
  m: bigint,
  n: bigint,
  index: bigint,
  from_list: bigint,
  operation: bigint,
  prepare: bigint
) => {
  const M = Number(m),
    N = Number(n);
  const A = new Matrix_mod2_dense(
    M,
    N,
    m4riPackedFixture(M, N, 42n).map((r) =>
      Array.from({ length: N }, (_, j) => Number((r >> BigInt(j)) & 1n))
    )
  );
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.echelon_form();
  else if (prepare === 3n) A.set_immutable();
  const C = A.row(index, Boolean(from_list)),
    D = A.row(index, true),
    before = C === D;
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const result = observe(() => {
    if (operation === 1n && N) C[0] = 1 - C[0]!;
    else if (operation === 2n && M && N) A.set(0, 0, 1 - A.get(0, 0));
    else if (operation === 3n) A.set(M, N, 1);
    else if (operation === 4n) A._clear_cache();
    else if (operation === 5n) A.swap_rows(0, M - 1);
    else if (operation === 6n) A.subdivide([0], [0]);
    else if (operation === 7n) A.echelonize();
    return null;
  });
  const E = A.row(index, Boolean(from_list)),
    F = A.row(index, true);
  return JSON.stringify([
    before,
    result,
    [C === E, D === F, E === F],
    C.map(String),
    E.map(String),
    A.list().map(String),
    observe(() => A.rank('bogus')),
  ]);
};

import { PrimeField } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
import { Mod } from '../../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod.js';
function binaryEntryValue(kind: bigint, value: bigint): any {
  return kind === 0n
    ? Number(value)
    : kind === 1n
      ? new Integer(value)
      : kind === 2n
        ? new Rational(value, 3n)
        : kind === 3n
          ? new Rational(value, 2n)
          : kind === 4n
            ? Number(value) / 2
            : kind === 5n
              ? Boolean(value)
              : kind === 6n
                ? String(value)
                : kind === 7n
                  ? null
                  : kind === 8n
                    ? new PrimeField(2n).__call__(value)
                    : kind === 9n
                      ? new PrimeField(3n).__call__(value)
                      : kind === 10n
                        ? Mod(value, 4n)
                        : kind === 11n
                          ? NaN
                          : Infinity;
}
functions.binary_constructor = (
  m: bigint,
  n: bigint,
  layout: bigint,
  kind: bigint,
  value: bigint
) => {
  const M = Number(m),
    N = Number(n),
    v = binaryEntryValue(kind, value),
    flat = Array(M * N).fill(v),
    nested = Array.from({ length: M }, () => Array(N).fill(v));
  const entries =
    layout === 0n
      ? undefined
      : layout === 1n
        ? null
        : layout === 2n
          ? v
          : layout === 3n
            ? flat
            : layout === 4n
              ? nested
              : layout === 5n
                ? flat.slice(0, -1)
                : layout === 6n
                  ? [...flat, v]
                  : layout === 7n
                    ? nested.slice(0, -1)
                    : layout === 8n
                      ? [...nested, Array(N).fill(v)]
                      : layout === 9n
                        ? nested.map((row) => row.slice(0, -1))
                        : layout === 10n
                          ? nested.map((row) => [...row, v])
                          : [nested[0] ?? [], v];
  const A = new Matrix_mod2_dense(M, N, entries);
  return JSON.stringify([A.nrows, A.ncols, A.list().map(String), A.is_mutable()]);
};
functions.binary_set_entry = (
  m: bigint,
  n: bigint,
  index: bigint,
  kind: bigint,
  value: bigint,
  prepare: bigint
) => {
  const M = Number(m),
    N = Number(n),
    A = new Matrix_mod2_dense(
      M,
      N,
      Array.from({ length: M }, (_, i) => Array.from({ length: N }, (_, j) => (i * N + j + 1) % 2))
    );
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.echelon_form();
  else if (prepare === 3n) A.set_immutable();
  const v = binaryEntryValue(kind, value),
    observe = (f: () => unknown) => {
      try {
        return ['ok', f()];
      } catch (e) {
        return ['error', (e as Error).name, (e as Error).message];
      }
    };
  return JSON.stringify([
    observe(() => {
      A.set(index, index, v);
      return null;
    }),
    A.list().map(String),
    observe(() => A.rank('bogus')),
  ]);
};
functions.binary_dimensions = (kind: bigint, m: bigint, n: bigint) => {
  const cv = (v: bigint): any =>
    kind === 0n
      ? Number(v)
      : kind === 1n
        ? v
        : kind === 2n
          ? new Integer(v)
          : kind === 3n
            ? new Rational(v, 2n)
            : kind === 4n
              ? Number(v) / 2
              : kind === 5n
                ? NaN
                : Infinity;
  const A = new Matrix_mod2_dense(cv(m), cv(n));
  return JSON.stringify([A.nrows, A.ncols, A.list().map(String)]);
};

import {
  set_random_seed,
  current_randstate,
} from '../../../../packages/sagemath-ts/src/misc/randstate.js';
functions.binary_random_state = (
  m: bigint,
  n: bigint,
  seed: bigint,
  density: bigint,
  nonzero: bigint,
  prepare: bigint,
  repeats: bigint
) => {
  const M = Number(m),
    N = Number(n),
    A = new Matrix_mod2_dense(
      M,
      N,
      Array.from({ length: M }, (_, i) => Array.from({ length: N }, (_, j) => (i * N + j + 1) % 2))
    );
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.echelon_form();
  else if (prepare === 3n) A.set_immutable();
  const densities = [undefined, -1, 0, 1 / 100, 1 / 3, 0.5, 0.999, 1, 2, NaN, Infinity, -Infinity];
  set_random_seed(seed);
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const states = [];
  for (let i = 0n; i < repeats; i++)
    states.push([
      observe(() => {
        A.randomize(densities[Number(density)], Boolean(nonzero));
        return null;
      }),
      A.list().map(String),
    ]);
  return JSON.stringify([
    states,
    Array.from({ length: 3 }, () => String(current_randstate().c_random())),
    observe(() => A.rank('bogus')),
    A.is_mutable(),
  ]);
};

import {
  zero_matrix_gf2,
  identity_matrix_gf2,
  random_matrix_gf2,
} from '../../../../packages/sagemath-ts/src/matrix/matrix_mod2.js';
functions.binary_factory = (
  op: bigint,
  m: bigint,
  n: bigint,
  kind: bigint,
  seed: bigint,
  density: bigint
) => {
  const cv = (v: bigint): any =>
    kind === 0n
      ? Number(v)
      : kind === 1n
        ? v
        : kind === 2n
          ? new Integer(v)
          : kind === 3n
            ? new Rational(v, 2n)
            : Number(v) / 2;
  const rows = cv(m),
    cols = n === -999n ? undefined : cv(n),
    densities = [undefined, -1, 0, 1 / 100, 1 / 3, 0.5, 0.999, 1, 2, NaN, Infinity, -Infinity];
  set_random_seed(seed);
  const A =
    op === 0n
      ? zero_matrix_gf2(rows, cols)
      : op === 1n
        ? identity_matrix_gf2(rows)
        : random_matrix_gf2(rows, cols, densities[Number(density)]);
  return JSON.stringify([
    A.nrows,
    A.ncols,
    A.list().map(String),
    A.is_mutable(),
    Array.from({ length: 3 }, () => String(current_randstate().c_random())),
  ]);
};

import { matrix_gf2_from_entries } from '../../../../packages/sagemath-ts/src/matrix/matrix_mod2.js';
functions.binary_from_entries = (
  m: bigint,
  n: bigint,
  layout: bigint,
  kind: bigint,
  value: bigint
) => {
  const M = Number(m),
    N = Number(n),
    v = binaryEntryValue(kind, value);
  let entries = Array.from({ length: M }, () => Array(N).fill(v));
  if (layout === 1n) entries = entries.slice(0, -1);
  else if (layout === 2n) entries.push(Array(N).fill(v));
  else if (layout === 3n && entries.length)
    entries[entries.length - 1] = entries.at(-1)!.slice(0, -1);
  else if (layout === 4n && entries.length) entries.at(-1)!.push(v);
  else if (layout === 5n && entries.length) entries[0] = [];
  const A = matrix_gf2_from_entries(entries);
  return JSON.stringify([A.nrows, A.ncols, A.list().map(String), A.is_mutable()]);
};

const binaryDensityStrings = [
  'bad',
  '',
  '0.5',
  ' 1_0e-2 ',
  'inf',
  '-Infinity',
  'nan',
  'NaN',
  '0x1',
  '1e9999',
  '-1e9999',
  '-0.0',
  '1__0',
  '1_.0',
  '\u0661.\u0665',
  '\uff11.\uff10',
  '\u00a00.5\u2003',
  '0.5\u001c',
  '\u001c0.5',
  '0.5\ufeff',
  '\ufeff0.5',
  '1\u0000',
  "a'b",
  'a"b',
  'a\'"b',
  '\ud83d\ude42',
  '\ud800',
];
function binaryDensityValue(kind: bigint): any {
  const values: any[] = [
    null,
    false,
    true,
    2n ** 2048n,
    new Integer(2n ** 2048n),
    -(2n ** 2048n),
    new Integer(-(2n ** 2048n)),
    new Rational(1n, 3n),
    new Rational(2n ** 2048n),
    new Rational(-1n, 2n ** 2048n),
    create_RealNumber('0.333333333333333333333333333'),
    create_RealNumber('0.5'),
    new PrimeField(2n).__call__(1n),
    new PrimeField(3n).__call__(2n),
    Mod(0n, 4n),
    [],
    {},
    new TextEncoder().encode('0.5'),
    new TextEncoder().encode('bad'),
    GF2.one(),
  ];
  return kind < 20n ? values[Number(kind)] : binaryDensityStrings[Number(kind) - 20];
}
functions.binary_density_input = (
  m: bigint,
  n: bigint,
  kind: bigint,
  prepare: bigint,
  nonzero: bigint
) => {
  const density = binaryDensityValue(kind);
  const M = Number(m),
    N = Number(n),
    A = new Matrix_mod2_dense(
      M,
      N,
      Array.from({ length: M }, (_, i) => Array.from({ length: N }, (_, j) => (i * N + j + 1) % 2))
    );
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.echelon_form();
  else if (prepare === 3n) A.set_immutable();
  set_random_seed(42n);
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  const result = observe(() => {
    A.randomize(density as any, Boolean(nonzero));
    return null;
  });
  return JSON.stringify([
    result,
    A.list().map(String),
    Array.from({ length: 3 }, () => String(current_randstate().c_random())),
    observe(() => A.rank('bogus')),
  ]);
};

import { float as python_float } from '../../../../packages/sagemath-ts/src/types/python_float.js';
functions.python_float_text = (codes: bigint[], as_bytes: bigint) => {
  const input = as_bytes
    ? Uint8Array.from(codes.map(Number))
    : codes.map((c) => String.fromCodePoint(Number(c))).join('');
  try {
    const value = python_float(input),
      view = new DataView(new ArrayBuffer(8));
    view.setFloat64(0, value);
    return JSON.stringify(['ok', view.getBigUint64(0).toString(16).padStart(16, '0')]);
  } catch (e) {
    return JSON.stringify(['error', (e as Error).name, (e as Error).message]);
  }
};

functions.binary_predicates = (
  m: bigint,
  n: bigint,
  mask: bigint,
  certificate: bigint,
  prepare: bigint
) => {
  const M = Number(m),
    N = Number(n),
    A = new Matrix_mod2_dense(
      M,
      N,
      Array.from({ length: M }, (_, i) =>
        Array.from({ length: N }, (_, j) => Number((mask >> BigInt(i * N + j)) & 1n))
      )
    );
  if (prepare === 1n) A.rank();
  else if (prepare === 2n) A.set_immutable();
  const observe = (f: () => unknown) => {
    try {
      return ['ok', f()];
    } catch (e) {
      return ['error', (e as Error).name, (e as Error).message];
    }
  };
  return JSON.stringify([
    A.is_Gamma_free(Boolean(certificate)),
    A.eq(A.copy()),
    A.eq(A.transpose()),
    A.eq(new Matrix_mod2_dense(M, N)),
    A.list().map(String),
    A.is_mutable(),
    observe(() => A.rank('bogus')),
  ]);
};

import {
  from_png_data,
  to_png_data,
} from '../../../../packages/sagemath-ts/src/matrix/matrix_mod2.js';
functions.binary_png_data = (m: bigint, n: bigint, mask: bigint, mode: bigint) => {
  const M = Number(m),
    N = Number(n),
    A = new Matrix_mod2_dense(
      M,
      N,
      Array.from({ length: M }, (_, i) =>
        Array.from({ length: N }, (_, j) => Number((mask >> BigInt(i * N + j)) & 1n))
      )
    );
  const { width, height, pixels } =
    mode === 0n
      ? to_png_data(A)
      : { width: N, height: M, pixels: Uint8Array.from(A.list(), (x) => (x ? 0 : 255)) };
  const B = from_png_data(width, height, pixels);
  return JSON.stringify([width, height, [...pixels], B.list().map(String), B.is_mutable()]);
};

functions.binary_density_factory = (m: bigint, n: bigint, kind: bigint) => {
  const density = binaryDensityValue(kind);
  set_random_seed(42n);
  let result;
  try {
    const A = random_matrix_gf2(m, n, density);
    result = ['ok', A.list().map(String), A.is_mutable()];
  } catch (e) {
    result = ['error', (e as Error).name, (e as Error).message];
  }
  return JSON.stringify([
    result,
    Array.from({ length: 3 }, () => String(current_randstate().c_random())),
  ]);
};
