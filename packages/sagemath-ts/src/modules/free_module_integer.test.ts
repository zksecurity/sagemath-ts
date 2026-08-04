/**
 * Tests for sage/modules/free_module_integer.py
 *
 * The expected values are SageMath's own doctest outputs (marked "Sage:")
 * or exact mathematical invariants (lattice volume, LLL-reducedness,
 * Voronoi relevance) verified against brute force.
 */

import { describe, expect, it } from 'bun:test';
import { type IntegerMatrix, IntegerMatrixFromEntries } from '../matrix/index.js';
import {
  IntegerLattice,
  bkzRootHermiteFactor,
  estimateBKZBlockSize,
  genLattice,
  isLLLReduced,
  lllReduce,
  qaryDualLattice,
  qaryLattice,
} from './free_module_integer.js';

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

function rows(M: IntegerMatrix): bigint[][] {
  const out: bigint[][] = [];
  for (let i = 0; i < M.nrows; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < M.ncols; j++) {
      row.push(M.get(i, j).value);
    }
    out.push(row);
  }
  return out;
}

/** Exact determinant by fraction-free (Bareiss) elimination. */
function det(M: IntegerMatrix): bigint {
  const n = M.nrows;
  if (n !== M.ncols) throw new Error('not square');
  const a = rows(M);
  let sign = 1n;
  let prev = 1n;
  for (let k = 0; k < n - 1; k++) {
    if (a[k]![k] === 0n) {
      let s = -1;
      for (let i = k + 1; i < n; i++) {
        if (a[i]![k] !== 0n) {
          s = i;
          break;
        }
      }
      if (s < 0) return 0n;
      const t = a[k]!;
      a[k] = a[s]!;
      a[s] = t;
      sign = -sign;
    }
    for (let i = k + 1; i < n; i++) {
      for (let j = k + 1; j < n; j++) {
        a[i]![j] = (a[i]![j]! * a[k]![k]! - a[i]![k]! * a[k]![j]!) / prev;
      }
    }
    prev = a[k]![k]!;
  }
  return sign * a[n - 1]![n - 1]!;
}

function absBig(x: bigint): bigint {
  return x < 0n ? -x : x;
}

function normSq(v: bigint[]): bigint {
  return v.reduce((s, x) => s + x * x, 0n);
}

/** A q-ary lattice basis [[q I_n], [A | I_{m-n}]] with deterministic entries. */
function qaryBasis(n: number, m: number, q: bigint, seed: number): IntegerMatrix {
  let state = seed >>> 0;
  const rng = () => {
    state ^= state << 13;
    state ^= state >>> 17;
    state ^= state << 5;
    return (state >>> 0) / 0x100000000;
  };
  const B: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const row: bigint[] = new Array(m).fill(0n);
    row[i] = q;
    B.push(row);
  }
  for (let i = 0; i < m - n; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < n; j++) {
      let v = BigInt(Math.floor(rng() * Number(q)));
      if (v > q / 2n) v -= q;
      row.push(v);
    }
    for (let j = 0; j < m - n; j++) row.push(i === j ? 1n : 0n);
    B.push(row);
  }
  return IntegerMatrixFromEntries(B);
}

// ---------------------------------------------------------------------------
// LLL (C13, H68)
// ---------------------------------------------------------------------------

describe('lllReduce', () => {
  it('reproduces the IntegerLattice doctest basis exactly', () => {
    // Sage: IntegerLattice([[1,0,3], [0,2,1], [0,2,7]])
    // User basis matrix:
    // [-2  0  0]
    // [ 0  2  1]
    // [ 1 -2  2]
    const L = IntegerLattice([
      [1, 0, 3],
      [0, 2, 1],
      [0, 2, 7],
    ]);
    expect(rows(L.reducedBasis)).toEqual([
      [-2n, 0n, 0n],
      [0n, 2n, 1n],
      [1n, -2n, 2n],
    ]);
  });

  it('reproduces the second IntegerLattice doctest basis exactly', () => {
    // Sage: IntegerLattice([[1,0,-2], [0,2,5], [0,0,7]])
    // [ 1  0 -2]
    // [ 1 -2  0]
    // [ 2  2  1]
    const L = IntegerLattice([
      [1, 0, -2],
      [0, 2, 5],
      [0, 0, 7],
    ]);
    expect(rows(L.reducedBasis)).toEqual([
      [1n, 0n, -2n],
      [1n, -2n, 0n],
      [2n, 2n, 1n],
    ]);
  });

  it('reproduces the voronoi_cell doctest basis exactly', () => {
    // Sage: IntegerLattice(Matrix(ZZ, 4, 4, [[0,0,1,-1],[1,-1,2,1],[-6,0,3,3],[-6,-24,-6,-5]]))
    // prints back the same user basis matrix
    const L = IntegerLattice([
      [0, 0, 1, -1],
      [1, -1, 2, 1],
      [-6, 0, 3, 3],
      [-6, -24, -6, -5],
    ]);
    expect(rows(L.reducedBasis)).toEqual([
      [0n, 0n, 1n, -1n],
      [1n, -1n, 2n, 1n],
      [-6n, 0n, 3n, 3n],
      [-6n, -24n, -6n, -5n],
    ]);
    // sqrt(L.discriminant()) == V.volume() == 678
    expect(L.discriminant()).toBe(678n * 678n);
    expect(L.volume()).toBe(678n);
  });

  it('defaults to delta = 0.99 (SageMath default), not 0.75', () => {
    const A = IntegerMatrixFromEntries([
      [1, 0, 3],
      [0, 2, 1],
      [0, 2, 7],
    ]);
    expect(rows(lllReduce(A))).toEqual(rows(lllReduce(A, { delta: 0.99 })));
    // The default output is 0.99-reduced, which 0.75-reduced output need not be.
    expect(isLLLReduced(lllReduce(A), 0.99)).toBe(true);
  });

  it('puts zero rows first for a rank-deficient generating set', () => {
    // Sage: matrix(4,3,[1,2,3,2,4,6,7,0,1,-1,-2,-3]).LLL()[0:2] == zero rows
    const R = lllReduce(
      IntegerMatrixFromEntries([
        [1, 2, 3],
        [2, 4, 6],
        [7, 0, 1],
        [-1, -2, -3],
      ])
    );
    expect(R.nrows).toBe(4);
    expect(rows(R)[0]).toEqual([0n, 0n, 0n]);
    expect(rows(R)[1]).toEqual([0n, 0n, 0n]);
    expect(rows(R)[2]!.some((x) => x !== 0n)).toBe(true);

    // Sage: matrix(ZZ, [[1,2,3],[31,41,51],[101,201,301]]).LLL() has one zero row
    const R2 = lllReduce(
      IntegerMatrixFromEntries([
        [1, 2, 3],
        [31, 41, 51],
        [101, 201, 301],
      ])
    );
    expect(rows(R2)[0]).toEqual([0n, 0n, 0n]);
    // Sage's remaining rows are (-1,0,1) and (1,1,1); ours agree up to sign
    expect(normSq(rows(R2)[1]!)).toBe(2n);
    expect(normSq(rows(R2)[2]!)).toBe(3n);
  });

  it('does not lose precision on entries far above 2^53 (C13)', () => {
    // The old floating point Gram-Schmidt either crashed ("Not an integer")
    // or returned a basis of a *different* lattice for these.
    for (const [m, e] of [
      [10, 30],
      [16, 30],
      [30, 40],
      [20, 60],
      [40, 60],
    ] as [number, number][]) {
      const q = 2n ** BigInt(e);
      const A = qaryBasis(1, m, q, 42);
      const R = lllReduce(A);
      expect(isLLLReduced(R, 0.99)).toBe(true);
      expect(absBig(det(R))).toBe(absBig(det(A)));
    }
  });

  it('constructs the lattice from the SageMath BKZ/HKZ doctest generator', () => {
    // Sage: A = sage.crypto.gen_lattice(type='random', n=1, m=30, q=2^40, seed=42)
    //       L = IntegerLattice(A)
    const A = genLattice({ type: 'random', n: 1, m: 30, q: 2n ** 40n, seed: 42 }) as IntegerMatrix;
    const L = IntegerLattice(A);
    expect(L.rank()).toBe(30);
    expect(isLLLReduced(L.reducedBasis, 0.99)).toBe(true);
    expect(absBig(det(L.reducedBasis))).toBe(absBig(det(A)));
  });

  it('rejects out-of-range parameters', () => {
    const I = IntegerMatrixFromEntries([
      [1, 0],
      [0, 1],
    ]);
    expect(() => lllReduce(I, { delta: 2 })).toThrow();
    expect(() => lllReduce(I, { delta: 0.25 })).toThrow();
    expect(() => lllReduce(I, { eta: 0.4 })).toThrow();
  });
});

describe('isLLLReduced', () => {
  it('rejects an unreduced basis and accepts its reduction', () => {
    const A = IntegerMatrixFromEntries([
      [1, 1000],
      [0, 1],
    ]);
    expect(isLLLReduced(A)).toBe(false);
    expect(isLLLReduced(lllReduce(A))).toBe(true);
  });

  it('raises on linearly dependent nonzero rows, like SageMath', () => {
    // Sage: matrix(ZZ, [[1,2,3],[2,4,6]]).is_LLL_reduced(algorithm='sage')
    // ValueError: linearly dependent input for module version of Gram-Schmidt
    expect(() =>
      isLLLReduced(
        IntegerMatrixFromEntries([
          [1, 2, 3],
          [2, 4, 6],
        ])
      )
    ).toThrow(/linearly dependent input/);
  });
});

// ---------------------------------------------------------------------------
// rank-deficient lattices (H70)
// ---------------------------------------------------------------------------

describe('rank-deficient lattices', () => {
  it('drops the zero rows LLL produces (Sage: rank 1)', () => {
    const L = IntegerLattice([
      [3, 0],
      [4, 0],
    ]);
    expect(L.rank()).toBe(1);
    expect(L.degree()).toBe(2);
    expect(rows(L.reducedBasis)).toEqual([[1n, 0n]]);
    expect(L.volume()).toBe(1n);
    expect(L.discriminant()).toBe(1n);
    expect(L.isUnimodular()).toBe(true);
    expect(Number.isFinite(L.hadamardRatio())).toBe(true);
  });

  it('computes the Voronoi relevant vectors of a rank 1 lattice', () => {
    // Sage: IntegerLattice([[3, 0], [4, 0]]).voronoi_relevant_vectors()
    //       [(-1, 0), (1, 0)]
    const L = IntegerLattice([
      [3, 0],
      [4, 0],
    ]);
    expect(L.voronoiRelevantVectors()).toEqual([
      [-1n, 0n],
      [1n, 0n],
    ]);
  });

  it('rejects a linearly dependent basis when lllReduce is off', () => {
    // `FreeModule_submodule_with_basis_pid.__init__` is called with
    // `check=True` (`free_module_integer.py:305-311`) and raises
    // (`free_module.py:6737-6738`).  Verified: SageMath 10.3 raises
    // `ValueError: The given basis vectors must be linearly independent.` for
    // `IntegerLattice([[1,2,3],[2,4,6],[7,0,1]], lll_reduce=False)`.  (This
    // construction previously succeeded here, which is why the zero-row test
    // below used it.)
    expect(() =>
      IntegerLattice(
        [
          [1, 2, 3],
          [2, 4, 6],
          [7, 0, 1],
        ],
        { lllReduce: false }
      )
    ).toThrow('the given basis vectors must be linearly independent');
  });

  it('keeps LLL/BKZ/HKZ output free of zero rows', () => {
    // With `lllReduce` on (the default) SageMath drops the zero rows LLL
    // produces for a rank-deficient generating set, leaving a genuine basis.
    const L = IntegerLattice([
      [1, 2, 3],
      [2, 4, 6],
      [7, 0, 1],
    ]);
    for (const M of [L.LLL(), L.BKZ({ blockSize: 2 }), L.HKZ()]) {
      expect(M.nrows).toBe(2);
      for (const row of rows(M)) {
        expect(row.some((x) => x !== 0n)).toBe(true);
      }
    }
  });
});

// ---------------------------------------------------------------------------
// Voronoi cell (H71)
// ---------------------------------------------------------------------------

describe('voronoiRelevantVectors / voronoiCell', () => {
  it('gives exactly the four unit vectors and the unit square for Z^2', () => {
    const L = IntegerLattice([
      [1, 0],
      [0, 1],
    ]);
    expect(L.voronoiRelevantVectors()).toEqual([
      [-1n, 0n],
      [0n, -1n],
      [0n, 1n],
      [1n, 0n],
    ]);

    // Voronoi cell of Z^2 is [-1/2, 1/2]^2: 2 x_i <= 1 and -2 x_i <= 1
    const { normals, offsets } = L.voronoiCell();
    expect(normals.length).toBe(4);
    const facets = normals.map((v, i) => `${v.join(',')}|${offsets[i]}`).sort();
    expect(facets).toEqual(['-2,0|1', '0,-2|1', '0,2|1', '2,0|1'].sort());
  });

  it('matches the SageMath H-representation of a non full rank lattice', () => {
    // Sage: IntegerLattice([[2,0,0],[0,2,0]]).voronoi_cell().Hrepresentation()
    // (An inequality (-1, 0, 0) x + 1 >= 0, (0, -1, 0) x + 1 >= 0,
    //  (1, 0, 0) x + 1 >= 0, (0, 1, 0) x + 1 >= 0)
    const L = IntegerLattice([
      [2, 0, 0],
      [0, 2, 0],
    ]);
    const { normals, offsets } = L.voronoiCell();
    const facets = normals.map((v, i) => `${v.join(',')}|${offsets[i]}`).sort();
    expect(facets).toEqual(['-1,0,0|1', '0,-1,0|1', '0,1,0|1', '1,0,0|1'].sort());
  });

  it('returns a set closed under negation whose members really are relevant', () => {
    const L = IntegerLattice([
      [2, 1, 0],
      [0, 3, 1],
      [1, 0, 4],
    ]);
    const relevant = L.voronoiRelevantVectors();
    expect(relevant.length).toBe(14);

    const key = (v: bigint[]) => v.join(',');
    const set = new Set(relevant.map(key));
    for (const v of relevant) {
      expect(set.has(key(v.map((x) => -x)))).toBe(true);
    }

    // v is relevant iff 0 and v are the only lattice points w with
    // |v/2 - w| <= |v/2|.  Check that by brute force over a coefficient box.
    const basis = rows(L.reducedBasis);
    for (const v of relevant) {
      let strict = true;
      for (let a = -3n; a <= 3n; a++) {
        for (let b = -3n; b <= 3n; b++) {
          for (let c = -3n; c <= 3n; c++) {
            const w = [0n, 1n, 2n].map(
              (j) =>
                a * basis[0]![Number(j)]! + b * basis[1]![Number(j)]! + c * basis[2]![Number(j)]!
            );
            if (w.every((x) => x === 0n) || w.every((x, j) => x === v[j])) continue;
            // 4|v/2 - w|^2 vs 4|v/2|^2, kept integral
            const lhs = normSq(v.map((x, j) => x - 2n * w[j]!));
            const rhs = normSq(v);
            if (lhs <= rhs) strict = false;
          }
        }
      }
      expect(strict).toBe(true);
    }
  });
});

// ---------------------------------------------------------------------------
// CVP (M80, M81, M82)
// ---------------------------------------------------------------------------

describe('approximateClosestVector', () => {
  const mk = () =>
    IntegerLattice(
      [
        [101, 0, 0, 0],
        [0, 101, 0, 0],
        [0, 0, 101, 0],
        [-28, 39, 45, 1],
      ],
      { lllReduce: false }
    );
  const t = [1337, 1337, 1337, 1337];

  it('honours delta (SageMath doctest)', () => {
    // Sage: L.approximate_closest_vector(t, delta=0.26) -> (1331, 1324, 1349, 1334)
    //       L.approximate_closest_vector(t, delta=0.99) -> (1326, 1349, 1339, 1345)
    expect(mk().approximateClosestVector(t, { delta: 0.26 })).toEqual([1331n, 1324n, 1349n, 1334n]);
    expect(mk().approximateClosestVector(t, { delta: 0.99 })).toEqual([1326n, 1349n, 1339n, 1345n]);
  });

  it("defaults to the 'embedding' algorithm", () => {
    expect(mk().approximateClosestVector(t)).toEqual(
      mk().approximateClosestVector(t, { algorithm: 'embedding' })
    );
    expect(mk().approximateClosestVector(t)).toEqual([1326n, 1349n, 1339n, 1345n]);
  });

  it('implements nearest_plane and rounding_off as in SageMath', () => {
    // Sage: nearest_plane -> (1326, 1349, 1339, 1345)
    //       rounding_off  -> (1331, 1324, 1349, 1334)
    expect(mk().approximateClosestVector(t, { algorithm: 'nearest_plane' })).toEqual([
      1326n,
      1349n,
      1339n,
      1345n,
    ]);
    expect(mk().approximateClosestVector(t, { algorithm: 'rounding_off' })).toEqual([
      1331n,
      1324n,
      1349n,
      1334n,
    ]);
  });

  it('rounds halves to even, like QQ.round("even")', () => {
    const Z2 = IntegerLattice([
      [1, 0],
      [0, 1],
    ]);
    // Sage: (2.5, 0.5) rounds to (2, 0), not (3, 1)
    expect(Z2.approximateClosestVector([2.5, 0.5], { algorithm: 'nearest_plane' })).toEqual([
      2n,
      0n,
    ]);
    expect(Z2.approximateClosestVector([2.5, 0.5], { algorithm: 'rounding_off' })).toEqual([
      2n,
      0n,
    ]);
  });

  it('takes exact rational targets', () => {
    // Sage: IntegerLattice([[1,0],[0,1]]).approximate_closest_vector((-6, 5/3)) -> (-6, 2)
    const Z2 = IntegerLattice([
      [1, 0],
      [0, 1],
    ]);
    expect(Z2.approximateClosestVector([-6, '5/3'])).toEqual([-6n, 2n]);
  });

  it('rejects unknown algorithms', () => {
    const Z2 = IntegerLattice([
      [1, 0],
      [0, 1],
    ]);
    expect(() =>
      Z2.approximateClosestVector([1, 1], {
        algorithm: 'nonsense' as unknown as 'embedding',
      })
    ).toThrow(/algorithm must be one of/);
  });
});

describe('closestVector', () => {
  it('solves the SageMath doctest', () => {
    // Sage: IntegerLattice([[1,0],[0,1]]).closest_vector((-6, 5/3)) -> (-6, 2)
    const Z2 = IntegerLattice([
      [1, 0],
      [0, 1],
    ]);
    expect(Z2.closestVector([-6, '5/3'])).toEqual([-6n, 2n]);
  });

  it('is exact, not merely Babai', () => {
    // Brute force the closest vector over a coefficient box and compare.
    const basisList: number[][][] = [
      [
        [1, 0, 3],
        [0, 2, 1],
        [0, 2, 7],
      ],
      [
        [4, 1],
        [3, 2],
      ],
      [
        [2, 1, 0],
        [0, 3, 1],
        [1, 0, 4],
      ],
    ];
    const targets: number[][] = [
      [5.5, 3.3, 10.1],
      [10.5, -7.25],
      [11.5, -4.25, 6.75],
    ];

    for (let k = 0; k < basisList.length; k++) {
      const L = IntegerLattice(basisList[k]!);
      const t = targets[k]!;
      const got = L.closestVector(t).map(Number);
      let gotDist = 0;
      for (let j = 0; j < t.length; j++) gotDist += (got[j]! - t[j]!) ** 2;

      const basis = rows(L.reducedBasis).map((r) => r.map(Number));
      let best = Number.POSITIVE_INFINITY;
      const recur = (i: number, cur: number[]) => {
        if (i === basis.length) {
          let d = 0;
          for (let j = 0; j < t.length; j++) d += (cur[j]! - t[j]!) ** 2;
          if (d < best) best = d;
          return;
        }
        for (let c = -8; c <= 8; c++) {
          recur(
            i + 1,
            cur.map((x, j) => x + c * basis[i]![j]!)
          );
        }
      };
      recur(0, new Array(t.length).fill(0));

      expect(gotDist).toBeCloseTo(best, 9);
    }
  });

  it('rejects targets of the wrong dimension', () => {
    const Z2 = IntegerLattice([
      [1, 0],
      [0, 1],
    ]);
    expect(() => Z2.closestVector([1, 2, 3])).toThrow();
  });
});

// ---------------------------------------------------------------------------
// SVP
// ---------------------------------------------------------------------------

describe('shortestVector', () => {
  it('returns an exact shortest vector', () => {
    // lambda_1(L)^2 = 5 for the lattice generated by (4,1) and (3,2):
    // (1,-1) = (4,1)-(3,2) has norm^2 2 ... check by brute force instead
    const L = IntegerLattice([
      [4, 1],
      [3, 2],
    ]);
    const sv = L.shortestVector({ updateReducedBasis: false });
    const basis = rows(L.reducedBasis);
    let best = normSq(sv);
    for (let a = -6n; a <= 6n; a++) {
      for (let b = -6n; b <= 6n; b++) {
        if (a === 0n && b === 0n) continue;
        const v = [0, 1].map((j) => a * basis[0]![j]! + b * basis[1]![j]!);
        const n = normSq(v);
        if (n > 0n && n < best) best = n;
      }
    }
    expect(normSq(sv)).toBe(best);
  });

  it('updates the reduced basis with the vector found', () => {
    const L = IntegerLattice(
      [
        [1, 1000],
        [0, 1],
      ],
      { lllReduce: false }
    );
    const before = normSq(rows(L.reducedBasis)[0]!);
    L.shortestVector();
    const after = normSq(rows(L.reducedBasis)[0]!);
    expect(after < before).toBe(true);
    expect(L.reducedBasis.nrows).toBe(2);
    expect(absBig(det(L.reducedBasis))).toBe(1n);
  });

  it('rejects unknown algorithms', () => {
    const L = IntegerLattice([
      [1, 0],
      [0, 1],
    ]);
    expect(() => L.shortestVector({ algorithm: 'nonsense' as unknown as 'pari' })).toThrow(
      /unknown/
    );
  });
});

// ---------------------------------------------------------------------------
// gen_lattice (H69)
// ---------------------------------------------------------------------------

describe('genLattice', () => {
  it('builds the primal basis [[q I_n], [A | I]] with |det| = q^n', () => {
    const B = genLattice({ type: 'modular', n: 4, m: 10, q: 11n, seed: 42 }) as IntegerMatrix;
    expect(B.nrows).toBe(10);
    expect(B.ncols).toBe(10);
    expect(absBig(det(B))).toBe(11n ** 4n);
    // First n rows are q*e_i
    expect(rows(B)[0]).toEqual([11n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n]);
  });

  it('builds the dual basis [[I_n, -A^T], [0, q I_{m-n}]] with |det| = q^(m-n)', () => {
    const B = genLattice({
      type: 'modular',
      n: 4,
      m: 10,
      q: 11n,
      seed: 42,
      dual: true,
    }) as IntegerMatrix;
    expect(B.nrows).toBe(10);
    expect(absBig(det(B))).toBe(11n ** 6n);

    // Sage's doctest for the dual has m - n = 6 rows containing q
    const withQ = rows(B).filter((row) => row.some((x) => x === 11n)).length;
    expect(withQ).toBe(6);
  });

  it('satisfies the primal/dual relation of the SageMath doctest', () => {
    // Sage: transpose(q*B_primal.inverse()).hermite_form() == B_dual.hermite_form()
    // Equivalently: B_dual * B_primal^T = 0 (mod q) and |det B_dual| = q^(m-n).
    const q = 11n;
    const P = genLattice({ type: 'modular', n: 4, m: 10, q, seed: 42 }) as IntegerMatrix;
    const D = genLattice({
      type: 'modular',
      n: 4,
      m: 10,
      q,
      seed: 42,
      dual: true,
    }) as IntegerMatrix;

    const pr = rows(P);
    const dr = rows(D);
    for (const d of dr) {
      for (const p of pr) {
        let s = 0n;
        for (let k = 0; k < 10; k++) s += d[k]! * p[k]!;
        expect(((s % q) + q) % q).toBe(0n);
      }
    }
    expect(absBig(det(D))).toBe(q ** 6n);
  });

  it('requires n = 1 for random bases', () => {
    expect(() => genLattice({ type: 'random', n: 2, m: 8, q: 11n })).toThrow(
      /random bases require n = 1/
    );
  });
});

// ---------------------------------------------------------------------------
// q-ary lattices (M83)
// ---------------------------------------------------------------------------

describe('qaryLattice', () => {
  it('uses A: Lambda_q(A) = {x : A x = 0 mod q}', () => {
    // A = [[1,2],[2,4]], q = 7: index 7 in Z^2, not 49
    const L = qaryLattice(
      [
        [1, 2],
        [2, 4],
      ],
      7n
    );
    expect(L.rank()).toBe(2);
    expect(L.volume()).toBe(7n);
    for (const row of rows(L.reducedBasis)) {
      expect((row[0]! + 2n * row[1]!) % 7n === 0n).toBe(true);
      expect((2n * row[0]! + 4n * row[1]!) % 7n === 0n).toBe(true);
    }
  });

  it('handles a composite modulus', () => {
    // {x : 3 x_0 = 0 mod 15, 5 x_1 = 0 mod 15} = 5Z x 3Z, volume 15
    const L = qaryLattice(
      [
        [3, 0],
        [0, 5],
      ],
      15n
    );
    expect(L.volume()).toBe(15n);
  });

  it('contains q Z^n', () => {
    const q = 7n;
    const L = qaryLattice(
      [
        [1, 2],
        [2, 4],
      ],
      q
    );
    // q*e_0 must be in the lattice: it is the closest vector to itself
    expect(L.closestVector([Number(q), 0])).toEqual([q, 0n]);
  });

  it('qaryDualLattice uses A^T', () => {
    const L = qaryDualLattice(
      [
        [1, 2],
        [2, 4],
      ],
      7n
    );
    expect(L.volume()).toBe(7n);
    for (const row of rows(L.reducedBasis)) {
      expect((row[0]! + 2n * row[1]!) % 7n === 0n).toBe(true);
    }
  });
});

// ---------------------------------------------------------------------------
// BKZ block size estimate (L51)
// ---------------------------------------------------------------------------

describe('estimateBKZBlockSize', () => {
  it('has a strictly decreasing root Hermite factor', () => {
    let prev = Number.POSITIVE_INFINITY;
    for (let beta = 2; beta <= 500; beta++) {
      const d = bkzRootHermiteFactor(beta);
      expect(d).toBeLessThanOrEqual(prev);
      prev = d;
    }
    expect(bkzRootHermiteFactor(2)).toBeCloseTo(1.0219, 4);
  });

  it('returns 2 when LLL already reaches the target', () => {
    // docstring example: dimension 100, volume 2^100, target 100
    expect(estimateBKZBlockSize(100, 2n ** 100n, 100)).toBe(2);
  });

  it('returns intermediate block sizes, not just 2 or n', () => {
    // vol^(1/n) = 2, target = 2 * delta^(n-1) for a prescribed delta
    const n = 200;
    const beta = (delta: number) => estimateBKZBlockSize(n, 2n ** 200n, 2 * delta ** (n - 1));
    const b1 = beta(1.02);
    const b2 = beta(1.014);
    const b3 = beta(1.01);
    expect(b1).toBeGreaterThan(2);
    expect(b1).toBeLessThan(b2);
    expect(b2).toBeLessThan(b3);
    expect(b3).toBeLessThan(n);
  });

  it('returns n when even full reduction cannot reach the target', () => {
    expect(estimateBKZBlockSize(20, 2n ** 200n, 1)).toBe(20);
  });

  it('does not overflow for cryptographic volumes', () => {
    // Number(2^2000) is Infinity; the estimate must still be finite.
    const beta = estimateBKZBlockSize(1000, 2n ** 2000n, 2 ** 10);
    expect(Number.isFinite(beta)).toBe(true);
    expect(beta).toBeGreaterThanOrEqual(2);
  });
});
