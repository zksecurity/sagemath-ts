/** Word modular square root from FLINT ulong_extras/sqrtmod.c.
 * The input residue must lie in [0,p), with p prime. The original's small
 * search, congruence shortcuts and bounded Tonelli iteration are preserved.
 * @see Deviation: FLINT word square roots and arithmetic
 */
import { n_jacobi_unsigned } from './jacobi.js';
import { n_is_square } from './is_square.js';
import { n_preinvert_limb } from './preinvert_limb.js';
import { n_powmod2_ui_preinv } from './powmod2_ui_preinv.js';
export function n_sqrtmod(a: bigint, p: bigint): bigint {
  if (p < 2n || p >= 1n << 64n || a < 0n || a >= p)
    throw new RangeError('square-root arguments require an unsigned word modulus and a reduced residue');
  if (a <= 1n) return a;
  if (p < 600n) {
    if (p > 50n && n_jacobi_unsigned(a, p) === -1) return 0n;
    let t = 0n, t2 = 0n;
    while (t < (p - 1n) / 2n) {
      t2 = (t2 + 2n * t + 1n) % p;
      t++;
      if (t2 === a) return t;
    }
    return 0n;
  }
  if (n_is_square(p) || !(p & 1n)) return 0n;
  const pinv = n_preinvert_limb(p);
  if (n_jacobi_unsigned(a, p) === -1) return 0n;
  if ((p & 3n) === 3n) return n_powmod2_ui_preinv(a, (p + 1n) / 4n, p, pinv);
  if ((p & 7n) === 5n) {
    const b = n_powmod2_ui_preinv(a, (p + 3n) / 8n, p, pinv);
    if (b * b % p === a) return b;
    return n_powmod2_ui_preinv(2n, (p - 1n) / 4n, p, pinv) * b % p;
  }
  let r = 0, p1 = p - 1n;
  do { p1 >>= 1n; r++; } while (!(p1 & 1n));
  let b = n_powmod2_ui_preinv(a, p1, p, pinv), k = 3n;
  while (n_jacobi_unsigned(k, p) !== -1) k += 2n;
  let g = n_powmod2_ui_preinv(k, p1, p, pinv);
  let res = n_powmod2_ui_preinv(a, (p1 + 1n) / 2n, p, pinv), iter = r - 1;
  while (b !== 1n) {
    let bpow = b, m = 0;
    do { bpow = bpow * bpow % p; m++; } while (m < r && bpow !== 1n);
    let gpow = g;
    for (let i = 1; i < r - m; i++) gpow = gpow * gpow % p;
    res = res * gpow % p;
    g = gpow * gpow % p;
    b = b * g % p;
    r = m;
    if (iter-- === 0) return 0n;
  }
  return res;
}
