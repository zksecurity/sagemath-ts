/** Portable single-thread nmod_mat/mul.c dispatch and mul_strassen.c schedule.
 * Dense rows replace native matrix windows; ncols retains empty matrix shapes.
 * @see Deviation: Native FLINT factor array kernels
 */
export function nmod_mat_mul(
  a: readonly (readonly bigint[])[],
  b: readonly (readonly bigint[])[],
  p: bigint,
  ncols = b[0]?.length ?? 0
): bigint[][] {
  const rows = a.length,
    inner = b.length,
    cols = ncols;
  if (p < 2n) throw new RangeError('modulus must be at least 2');
  if (
    !Number.isSafeInteger(cols) ||
    cols < 0 ||
    a.some((r) => r.length !== inner) ||
    b.some((r) => r.length !== cols)
  )
    throw new RangeError('incompatible matrix dimensions');
  const mod = (x: bigint) => {
    const r = x % p;
    return r < 0n ? r + p : r;
  };
  const out = Array.from({ length: rows }, () => Array<bigint>(cols).fill(0n));
  if (Math.min(rows, inner, cols) < (p < 2048n ? 400 : 200)) {
    for (let i = 0; i < rows; i++)
      for (let j = 0; j < cols; j++) {
        let sum = 0n;
        for (let k = 0; k < inner; k++) sum += a[i]![k]! * b[k]![j]!;
        out[i]![j] = mod(sum);
      }
    return out;
  }
  const ar = Math.floor(rows / 2),
    ac = Math.floor(inner / 2),
    bc = Math.floor(cols / 2);
  const block = (v: readonly (readonly bigint[])[], r: number, c: number, nr: number, nc: number) =>
    v.slice(r, r + nr).map((row) => row.slice(c, c + nc));
  const add = (x: bigint[][], y: bigint[][], sign = 1n) =>
    x.map((r, i) => r.map((v, j) => mod(v + sign * y[i]![j]!)));
  const mul = (x: bigint[][], y: bigint[][]) => nmod_mat_mul(x, y, p);
  const A11 = block(a, 0, 0, ar, ac),
    A12 = block(a, 0, ac, ar, ac),
    A21 = block(a, ar, 0, ar, ac),
    A22 = block(a, ar, ac, ar, ac);
  const B11 = block(b, 0, 0, ac, bc),
    B12 = block(b, 0, bc, ac, bc),
    B21 = block(b, ac, 0, ac, bc),
    B22 = block(b, ac, bc, ac, bc);
  let X1 = add(A22, A12),
    X2 = add(B22, B12),
    C21 = mul(X1, X2);
  X1 = add(A22, A21, -1n);
  X2 = add(B22, B21, -1n);
  let C22 = mul(X1, X2);
  X1 = add(X1, A12);
  X2 = add(X2, B12);
  let C11 = mul(X1, X2);
  X1 = add(X1, A11, -1n);
  let C12 = mul(X1, B12);
  X1 = mul(A12, B21);
  C11 = add(C11, X1);
  C12 = add(C12, C22);
  C12 = add(C11, C12, -1n);
  C11 = add(C21, C11, -1n);
  X2 = add(X2, B11, -1n);
  C21 = mul(A21, X2);
  C21 = add(C11, C21, -1n);
  C22 = add(C22, C11);
  C11 = add(mul(A11, B11), X1);
  for (let i = 0; i < ar; i++)
    for (let j = 0; j < bc; j++) {
      out[i]![j] = C11[i]![j]!;
      out[i]![j + bc] = C12[i]![j]!;
      out[i + ar]![j] = C21[i]![j]!;
      out[i + ar]![j + bc] = C22[i]![j]!;
    }
  if (cols % 2)
    for (let i = 0; i < rows; i++) {
      let s = 0n;
      for (let k = 0; k < inner; k++) s += a[i]![k]! * b[k]![cols - 1]!;
      out[i]![cols - 1] = mod(s);
    }
  if (rows % 2)
    for (let j = 0; j < cols; j++) {
      let s = 0n;
      for (let k = 0; k < inner; k++) s += a[rows - 1]![k]! * b[k]![j]!;
      out[rows - 1]![j] = mod(s);
    }
  if (inner % 2)
    for (let i = 0; i < 2 * ar; i++)
      for (let j = 0; j < 2 * bc; j++)
        out[i]![j] = mod(out[i]![j]! + a[i]![inner - 1]! * b[inner - 1]![j]!);
  return out;
}
