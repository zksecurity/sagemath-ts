/** Native delayed-reduction determinant over an arbitrary-size NTL modulus.
 * @see Deviation: NTL integer CRT and modular determinant adapters
 */
export function determinant(
  A: readonly (readonly bigint[])[],
  p: bigint,
  options: { columns?: number } = {}
): bigint {
  if (p <= 1n) throw new Error('ZZ_pContext: p must be > 1');
  const columns = options.columns ?? A[0]?.length ?? 0;
  if (columns < 0) throw new Error('SetDims: bad args');
  if (A.some((row) => row.length !== columns)) throw new Error('nonrectangular matrix');
  const n = A.length;
  if (columns !== n) throw new Error('determinant: nonsquare matrix');
  if (!n) return 1n;
  const norm = (x: bigint) => ((x % p) + p) % p;
  const M = A.map((row) => row.map(norm));
  let det = 1n;
  for (let k = 0; k < n; k++) {
    let pos = -1;
    for (let i = k; i < n; i++) {
      M[i]![k] = norm(M[i]![k]!);
      if (pos === -1 && M[i]![k] !== 0n) pos = i;
    }
    if (pos === -1) return 0n;
    if (k !== pos) {
      [M[k], M[pos]] = [M[pos]!, M[k]!];
      det = norm(-det);
    }
    det = (det * M[k]![k]!) % p;
    let u = M[k]![k]!,
      v = p,
      s = 1n,
      t = 0n;
    while (v) {
      const q = u / v;
      [u, v] = [v, u - q * v];
      [s, t] = [t, s - q * t];
    }
    if (u !== 1n) throw new Error('InvMod: inverse undefined');
    const inverse = norm(-s);
    for (let j = k + 1; j < n; j++) M[k]![j] = (norm(M[k]![j]!) * inverse) % p;
    for (let i = k + 1; i < n; i++) {
      const factor = M[i]![k]!;
      for (let j = k + 1; j < n; j++) M[i]![j] = M[i]![j]! + M[k]![j]! * factor;
    }
  }
  return det;
}
