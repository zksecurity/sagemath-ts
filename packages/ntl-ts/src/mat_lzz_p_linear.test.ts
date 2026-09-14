import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';


const fixtures = (await loadLiveNative(import.meta.url, "./mat_lzz_p_linear.native.json.gz")) as { args: string[]; result: string | null; error: string | null; errorType: string | null }[];
import {
  gauss,
  image,
  kernel,
  relaxed_inv,
  inv,
  relaxed_determinant,
  determinant,
  relaxed_solve,
  solve,
} from './mat_lzz_p.js';
export function ntl_word_linear(
  p: bigint,
  dims: bigint[],
  av: bigint[],
  b: bigint[],
  prev: bigint[],
  xprev: bigint[],
  cmd: bigint[]
): string {
  if (p <= 1n) throw new Error('zz_pContext: p must be > 1');
  if (p >= 1n << 60n) throw new Error('zz_pContext: modulus too big');
  const [n, m, r, c] = dims.map(Number);
  if (n! < 0 || m! < 0 || r! < 0 || c! < 0) throw new Error('SetDims: bad args');
  const norm = (x: bigint) => ((x % p) + p) % p;
  const A = Array.from({ length: n! }, (_, i) => av.slice(i * m!, (i + 1) * m!).map(norm)),
    previous = Array.from({ length: r! }, (_, i) => prev.slice(i * c!, (i + 1) * c!).map(norm)),
    out: unknown[] = [];
  for (let i = 0; i < cmd.length; i += 4) {
    const op = Number(cmd[i]),
      w = Number(cmd[i + 1]),
      relax = !!cmd[i + 2],
      left = !!cmd[i + 3];
    let result: unknown = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      if (op === 0) result = gauss(A, p, { w, columns: m });
      else if (op === 1) result = gauss(A, p, { columns: m });
      else if (op === 2) result = image(A, p, { columns: m });
      else if (op === 3) result = kernel(A, p, { columns: m });
      else if (op === 4) result = relaxed_inv(A, p, { columns: m, relax, previous });
      else if (op === 5) result = inv(A, p, { columns: m });
      else if (op === 6) result = relaxed_determinant(A, p, { columns: m, relax });
      else if (op === 7) result = determinant(A, p, { columns: m });
      else if (op === 8)
        result = relaxed_solve(A, b, p, { columns: m, relax, left, previous: xprev });
      else if (op === 9) result = solve(A, b, p, { columns: m, left, previous: xprev });
      else if (op === 10) result = relaxed_inv(A, p, { columns: m, previous });
      else if (op === 11) result = relaxed_determinant(A, p, { columns: m });
      else if (op === 12) result = relaxed_solve(A, b, p, { columns: m, previous: xprev });
      else if (op === 13 || op === 14) {
        const bad = op === 13 ? [[0n], [0n, 0n]] : [],
          options = op === 14 ? { columns: -1 } : {};
        if (w === 0) result = gauss(bad, p, options);
        else if (w === 1) result = image(bad, p, options);
        else if (w === 2) result = kernel(bad, p, options);
        else if (w === 3) result = relaxed_inv(bad, p, options);
        else if (w === 4) result = inv(bad, p, options);
        else if (w === 5) result = relaxed_determinant(bad, p, options);
        else if (w === 6) result = determinant(bad, p, options);
        else if (w === 7) result = relaxed_solve(bad, [], p, options);
        else result = solve(bad, [], p, options);
      } else throw new Error('unknown word linear-algebra operation');
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    out.push([errorType, error, result]);
  }
  out.push([null, null, A]);
  return JSON.stringify(out, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}

const parse = (s: string): bigint[] =>
  s
    .slice(1, -1)
    .split(',')
    .filter((x) => x.trim())
    .map((x) => BigInt(x.trim()));
for (const [i, row] of fixtures.entries())
  test('native NTL word linear algebra ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_word_linear(
        BigInt(row.args[0]!),
        parse(row.args[1]!),
        parse(row.args[2]!),
        parse(row.args[3]!),
        parse(row.args[4]!),
        parse(row.args[5]!),
        parse(row.args[6]!)
      );
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
  }, 30000);

test('NTL elimination outputs own rows independently of input and other results', () => {
  const A = [
      [1n, 2n],
      [2n, 4n],
      [0n, 1n],
    ],
    before = structuredClone(A);
  const G = gauss(A, 5n)[1],
    I = image(A, 5n),
    K = kernel(A, 5n);
  G[0]![0] = 999n;
  expect(I[0]![0]).toBe(1n);
  I[0]![0] = 777n;
  K[0]![0] = 888n;
  expect(A).toEqual(before);
});
test('NTL singular inverse and solve retain independent normalized previous outputs', () => {
  const A = [[3n]],
    previous = [[7n]],
    xprevious = [11n, 12n],
    b = [1n];
  const [d, X] = relaxed_inv(A, 9n, { previous });
  expect(d).toBe(0n);
  expect(X).toEqual([[7n]]);
  X[0]![0] = 999n;
  const [e, x] = relaxed_solve(A, b, 9n, { previous: xprevious });
  expect(e).toBe(0n);
  expect(x).toEqual([2n, 3n]);
  x[0] = 999n;
  expect(previous).toEqual([[7n]]);
  expect(xprevious).toEqual([11n, 12n]);
  expect(A).toEqual([[3n]]);
  expect(b).toEqual([1n]);
});
