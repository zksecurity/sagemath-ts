/** PARI kernel/none/halfgcd.c; matrices are zero-indexed rows.
 * Native Lehmer scalar batches are expanded into their exact Euclidean steps;
 * the native 66-word recursive split and fixups remain in place.
 * @see Deviation: Number-field ideal backend adapters
 */
type Matrix = [[bigint, bigint], [bigint, bigint]];
type State = [Matrix, bigint, bigint];
const abs = (a: bigint) => (a < 0n ? -a : a);
const exp = (a: bigint): number => (a === 0n ? -1000000000 : abs(a).toString(2).length - 1);
const ceilExp = (a: bigint): number => exp(a) + ((a & (a - 1n)) === 0n ? 0 : 1);
const threshold = (a: bigint) => Math.floor((3 + ceilExp(a)) / 2);
const identity = (): Matrix => [
  [1n, 0n],
  [0n, 1n],
];
const determinant = (M: Matrix) => M[0][0] * M[1][1] - M[0][1] * M[1][0];
function mulq(M: Matrix, q: bigint): Matrix {
  return [
    [M[0][0] * q + M[0][1], M[0][0]],
    [M[1][0] * q + M[1][1], M[1][0]],
  ];
}
function undo([M, a, b]: State, q: bigint): State {
  return [
    [
      [M[0][1], M[0][0] - M[0][1] * q],
      [M[1][1], M[1][0] - M[1][1] * q],
    ],
    a * q + b,
    a,
  ];
}
function step([M, a, b]: State, q: bigint): State {
  return [mulq(M, q), b, a - b * q];
}
function lastq(M: Matrix): bigint {
  const p = M[0][0],
    q = M[0][1],
    s = M[1][1];
  if (!q) throw new Error('bug in halfgcd');
  return s === 0n ? p : q === 1n ? p - 1n : p / q;
}
function multiply(A: Matrix, B: Matrix): Matrix {
  if ([...A.flat(), ...B.flat()].some((v) => Math.ceil((exp(v) + 1) / 64) < 14))
    return [
      [A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
      [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]],
    ];
  const [[a, b], [c, d]] = A,
    [[e, f], [g, h]] = B,
    M1 = (a + d) * (e + h),
    M2 = (c + d) * e,
    M3 = a * (f - h),
    M4 = d * (g - e),
    M5 = (a + b) * h,
    M6 = (c - a) * (e + f),
    M7 = (b - d) * (g + h);
  return [
    [M1 + M4 + M7 - M5, M3 + M5],
    [M2 + M4, M1 - M2 + M3 + M6],
  ];
}
function fix0(M: Matrix, a: bigint, b: bigint, m: number): State {
  let count = 0;
  while (exp(b) >= m) {
    [M, a, b] = step([M, a, b], a / b);
    count++;
  }
  if (count > 6) throw new Error('bug in FIXUP0');
  return [M, a, b];
}
function fix1([Q, a0, b0]: State, a: bigint, b: bigint, m: number, t: number): State {
  const shift = BigInt(m),
    mask = (1n << shift) - 1n,
    am = a & mask,
    bm = b & mask;
  let ap: bigint, bp: bigint;
  if (determinant(Q) === -1n) {
    ap = bm * Q[0][1] - am * Q[1][1] + ((a0 + Q[1][1]) << shift);
    bp = am * Q[1][0] - bm * Q[0][0] + ((b0 - Q[1][0]) << shift);
    if (bp >= 0n) return [Q, ap, bp];
    if (exp(ap + bp) >= m + t)
      return [
        [
          [Q[0][0], Q[0][0] - Q[0][1]],
          [Q[1][0], Q[1][0] - Q[1][1]],
        ],
        ap + bp,
        -bp,
      ];
    const q = lastq(Q),
      state = undo([Q, ap, bp], q);
    return q >= 2n ? step(state, q - 1n) : undo(state, lastq(state[0]));
  }
  ap = am * Q[1][1] - bm * Q[0][1] + ((a0 - Q[1][1]) << shift);
  bp = bm * Q[0][0] - am * Q[1][0] + ((b0 + Q[1][0]) << shift);
  return exp(ap) >= m + t
    ? fix0(Q, ap, bp, m + t)
    : Q[0][1] === 0n
      ? [Q, ap, bp]
      : undo([Q, ap, bp], lastq(Q));
}
function hgcd(a: bigint, b: bigint): State {
  const m = threshold(a);
  if (Math.ceil((exp(b) + 1) / 64) < 66) {
    let M: Matrix =
      a === b
        ? [
            [0n, 1n],
            [1n, 0n],
          ]
        : identity();
    while (exp(b) >= m) [M, a, b] = step([M, a, b], a / b);
    return [M, a, b];
  }
  if (exp(b) < m) return [identity(), a, b];
  const a0 = (a >> BigInt(m)) + 1n;
  if (a0 <= 7n) return fix0(identity(), a, b, m);
  const b0 = b >> BigInt(m),
    t = threshold(a0);
  let [R, ap, bp] = fix1(hgcd(a0, b0), a, b, m, t);
  if (exp(bp) < m) return [R, ap, bp];
  const q = ap / bp,
    c = bp,
    d = ap % bp;
  if (c >> BigInt(m) <= 6n) return fix0(mulq(R, q), c, d, m);
  const k = 2 * m - ceilExp(c) - 1;
  if (k < 0) throw new Error('bug in halfgcd');
  const c0 = (c >> BigInt(k)) + 1n,
    d0 = d >> BigInt(k),
    tp = threshold(c0);
  if (c0 < 8n) throw new Error('bug in halfgcd');
  const [S, cp, dp] = fix1(hgcd(c0, d0), c, d, k, tp);
  if (!(exp(cp) >= m + 1 && m + 1 > exp(dp))) throw new Error('bug in halfgcd');
  return fix0(multiply(mulq(R, q), S), cp, dp, m);
}
function hgcd0(x: bigint, y: bigint): State {
  if (y >= 0n && x >= y) return hgcd(x, y);
  if (x < y) {
    const [Q, a, b] = hgcd0(y, x);
    return [[Q[1], Q[0]], a, b];
  }
  if (x <= 0n) {
    const [Q, a, b] = hgcd(-y, -x);
    return [
      [
        [-Q[1][0], -Q[1][1]],
        [-Q[0][0], -Q[0][1]],
      ],
      a,
      b,
    ];
  }
  const [Q, a, b] = hgcd0(x, -y);
  return [[Q[0], [-Q[1][0], -Q[1][1]]], a, b];
}
/** Return M,V with V=M*[A,B], matching the native stopping point. */
export function halfgcdii(A: bigint, B: bigint): [Matrix, [bigint, bigint]] {
  const m = abs(A) > abs(B) ? abs(A) : abs(B);
  let [Q, a, b] = hgcd0(A, B);
  while (b !== 0n && b * b >= m) [Q, a, b] = step([Q, a, b], a / b);
  const s = determinant(Q);
  return [
    [
      [s * Q[1][1], -s * Q[0][1]],
      [-s * Q[1][0], s * Q[0][0]],
    ],
    [a, b],
  ];
}
