/** Internal row-oriented transcription of PARI ZV.c/FpV.c matrix products. */
export type Rows = bigint[][];
export type Product = (a: Rows, b: Rows, m: number, n: number, p: number) => Rows;
export const reduceMod = (x: bigint, p: bigint): bigint => ((x % p) + p) % p;
export function fromColumns(a: bigint[][]): Rows {
  const rows = a.length > 1 ? a[1]!.length - 1 : 0;
  return Array.from({ length: rows }, (_, i) => a.slice(1).map((c) => c[i + 1]!));
}
export function toColumns(a: Rows, columns: number): bigint[][] {
  return [[], ...Array.from({ length: columns }, (_, j) => [0n, ...a.map((r) => r[j]!)])];
}
export function classical(
  a: Rows,
  b: Rows,
  m: number,
  n: number,
  p: number,
  modulus?: bigint
): Rows {
  return Array.from({ length: m }, (_, i) =>
    Array.from({ length: p }, (_, j) => {
      let s = 0n;
      for (let k = 0; k < n; k++) s += a[i]![k]! * b[k]![j]!;
      return modulus === undefined ? s : reduceMod(s, modulus);
    })
  );
}
/** Seven-product Winograd schedule, including native padding for odd dimensions. */
export function winograd(
  a: Rows,
  b: Rows,
  m: number,
  n: number,
  p: number,
  mul: Product,
  modulus?: bigint
): Rows {
  const m1 = Math.ceil(m / 2),
    m2 = Math.floor(m / 2),
    n1 = Math.ceil(n / 2),
    n2 = Math.floor(n / 2),
    p1 = Math.ceil(p / 2),
    p2 = Math.floor(p / 2);
  const slice = (v: Rows, i: number, j: number, h: number, w: number): Rows =>
    Array.from({ length: h }, (_, r) =>
      Array.from({ length: w }, (_, c) => v[i + r]?.[j + c] ?? 0n)
    );
  const add = (v: Rows, w: Rows, h: number, k: number, sign = 1n): Rows =>
    Array.from({ length: h }, (_, i) =>
      Array.from({ length: k }, (_, j) => {
        const z = (v[i]?.[j] ?? 0n) + sign * (w[i]?.[j] ?? 0n);
        return modulus === undefined ? z : reduceMod(z, modulus);
      })
    );
  const A11 = slice(a, 0, 0, m1, n1),
    A12 = slice(a, 0, n1, m1, n2),
    A21 = slice(a, m1, 0, m2, n1),
    A22 = slice(a, m1, n1, m2, n2);
  const B11 = slice(b, 0, 0, n1, p1),
    B12 = slice(b, 0, p1, n1, p2),
    B21 = slice(b, n1, 0, n2, p1),
    B22 = slice(b, n1, p1, n2, p2);
  const T2 = add(B12, B22, n1, p2, -1n),
    S1 = add(A21, A11, m2, n1, -1n),
    M2 = mul(S1, T2, m2, n1, p2);
  const T3 = add(T2, B11, n1, p1, -1n),
    S2 = add(A21, A22, m2, n1),
    T1 = add(B12, B11, n1, p1, -1n),
    M3 = mul(S2, T1, m2, n1, p2);
  const S3 = add(S2, A11, m1, n1, -1n),
    M1 = mul(A11, B11, m1, n1, p1),
    M4 = mul(A12, B21, m1, n2, p1);
  const C11 = add(M1, M4, m1, p1),
    M5 = mul(S3, T3, m1, n1, p1),
    S4 = add(A12, S3, m1, n2, -1n),
    T4 = add(B21, T3, n2, p1);
  const V1 = add(M1, M5, m1, p1, -1n),
    M6 = mul(S4, B22, m1, n2, p2),
    M7 = mul(A22, T4, m2, n2, p1);
  const C12 = add(add(V1, M3, m1, p2), M6, m1, p2),
    V2 = add(V1, M2, m2, p1),
    C21 = add(V2, M7, m2, p1),
    C22 = add(V2, M3, m2, p2);
  return [...C11.map((r, i) => r.concat(C12[i]!)), ...C21.map((r, i) => r.concat(C22[i]!))];
}
