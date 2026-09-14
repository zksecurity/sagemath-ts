/** PARI extension-polynomial composition, retaining native evaluation schedules. */
import { wordExtensionMatrix } from './_extension_matrix.js';
import { extensionField, trimExtension as trim } from './_extension_field.js';
import {
  extensionPolynomial as arithmetic,
  validateExtensionInputs,
  type ExtensionPolynomial as P,
} from './_extension_polynomial.js';
import {
  extensionGetRed,
  type ExtensionModulus as M,
} from './_extension_division.js';
import { extensionQuotient } from './_extension_quotient.js';
import { polynomialQuotient } from './_polynomial_quotient.js';
import { PariError } from './errors.js';
/** Native Brent–Kung blocks or packed word matrix evaluation.
 * @see Deviation: PARI extension composition adapters
 */
export function extensionComposition(
  mode: 0 | 1 | 2,
  op: 0 | 1,
  p: bigint,
  T: bigint[] | bigint,
  Q: P,
  x: P,
  V: P[],
  S: M,
  providedInverse?: bigint[]
): P {
  validateExtensionInputs(mode, p, T, Q, x);
  validateExtensionInputs(mode, p, T, Array.isArray(S) ? S : S.polynomial);
  for (const v of V) validateExtensionInputs(mode, p, T, v);
  const canonical = (x: P) =>
    trim(
      x.map((c) =>
        Array.isArray(c)
          ? (() => {
              let n = c.length;
              while (n && c[n - 1] === 0n) n--;
              return c.slice(0, n);
            })()
          : c
      )
    );
  Q = canonical(Q);
  x = canonical(x);
  V = V.map(canonical);
  S = Array.isArray(S) ? canonical(S) : S;
  if (Array.isArray(T)) {
    T = T.slice();
    while (T.length && T.at(-1) === 0n) T.pop();
  }
  let inverse = providedInverse;
  const prepareT = () => {
    if (mode !== 2 && inverse === undefined)
      inverse = polynomialQuotient(T as bigint[], p, mode === 1).inverse;
  };
  if (mode === 1 && op === 0 && !Q.length) return [];
  const ready = () => {
    prepareT();
    S = extensionGetRed(mode, S, T, p, inverse);
  };
  if (mode !== 1 || op === 0) ready();
  if (op === 0) {
    if (!Q.length) return [];
    V = extensionQuotient(
      mode,
      6,
      p,
      BigInt(Math.floor(Math.sqrt(Q.length - 1))),
      T,
      x,
      [],
      S,
      inverse
    ) as P[];
  }
  let f = extensionField(mode, T, p, inverse);
  const mul = (a: P, b: P) => extensionQuotient(mode, 0, p, 0n, T, a, b, S, inverse) as P;
  const add = (a: P, b: P) =>
    trim(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) =>
        i < a.length ? (i < b.length ? f.add(a[i]!, b[i]!) : a[i]!) : b[i]!
      )
    );
  if (mode === 1) {
    if (!Q.length) return [];
    const raw = Array.isArray(S) ? S : S.polynomial;
    const m = raw.length - 1,
      l = V.length;
    if (!l || (l === 1 && Q.length > 1)) throw new RangeError('power table is too short');
    if (m < 0) throw new RangeError('word composition modulus must be nonzero');
    const n = Q.length <= l ? l : l - 1,
      d = Q.length <= l ? 1 : Math.ceil(Q.length / n);
    const A = V.slice(0, n).map((a) =>
      Array.from({ length: m }, (_, i) => (a[i] ?? []) as bigint[])
    );
    const B = Array.from({ length: d }, (_, i) =>
      Array.from({ length: n }, (_, j) => (Q[i * n + j] ?? []) as bigint[])
    );
    // matrix retains caller-provided inner cache; ordinary matrix n>1 does not prepare one.
    const C = wordExtensionMatrix(A, B, T as bigint[], p, inverse);
    ready();
    f = extensionField(mode, T, p, inverse);
    let R = trim(C[d - 1]!);
    for (let i = d - 2; i >= 0; i--) R = add(mul(R, V[l - 1]!), trim(C[i]!));
    return R;
  }
  const l = V.length,
    d0 = Q.length - 1;
  if (d0 < 0) return [];
  const block = (a: number, n: number) => {
    let z = trim([f.mul(f.u, Q[a]!, false)]);
    for (let i = 1; i <= n; i++) z = add(z, trim(V[i]!.map((c) => f.mul(c, Q[a + i]!, false))));
    return arithmetic(mode, 2, p, T, z, [], inverse);
  };
  if (d0 < l) return block(0, d0);
  if (l < 2) throw new PariError('domain error in gen_RgX_bkeval_powers: #powers < 2');
  let d = d0 - l,
    z = block(d + 1, l - 1);
  while (d >= l - 1) {
    d -= l - 1;
    z = add(block(d + 1, l - 2), mul(z, V[l - 1]!));
  }
  z = add(block(0, d), mul(z, V[d + 1]!));
  return arithmetic(mode, 2, p, T, z, [], inverse);
}
