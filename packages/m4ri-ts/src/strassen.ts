import { type mzd_t, mzd_init, mzd_init_window as window, mzd_add as add } from './mzd.js';
import { mzd_mul_m4rm } from './brilliantrussian.js';
/** M4RI Strassen-Winograd/Bodrato multiplication, with a distinct squaring schedule.
 * Source: strassen.c _mzd_mul_even, _mzd_sqr_even and mzd_mul (20251207).
 * Row values are immutable; owned submatrices replace native writable windows.
 */
export function mzd_mul(a: mzd_t, b: mzd_t, cutoff = 0): mzd_t {
  if (a.ncols !== b.nrows) throw new RangeError('incompatible M4RI multiplication dimensions');
  if (!Number.isInteger(cutoff)) throw new RangeError('M4RI cutoff must be an integer');
  if (cutoff < 0) throw Object.assign(new Error('Aborted'), { name: 'RuntimeError' });
  cutoff = Math.max(64, Math.floor((cutoff === 0 ? 4096 : cutoff) / 64) * 64);
  const m = a.nrows,
    k = a.ncols,
    n = b.ncols;
  if (m === 0 || n === 0) return mzd_init(m, n);
  if (3 * m < 4 * cutoff || 3 * k < 4 * cutoff || 3 * n < 4 * cutoff) return mzd_mul_m4rm(a, b);
  let multiple = 64,
    width = Math.floor(Math.min(m, k, n) / 2);
  while (width > cutoff) {
    width = Math.floor(width / 2);
    multiple *= 2;
  }
  const mh = Math.floor((m - (m % multiple)) / 128) * 64,
    kh = Math.floor((k - (k % multiple)) / 128) * 64,
    nh = Math.floor((n - (n % multiple)) / 128) * 64;
  // The pinned native code reaches invalid empty windows at very small
  // cutoffs. Preserve the signal-guarded error instead of terminating JS:
  // empty right-column transpose -> mzd_copy abort; empty inner window
  // with nonempty output -> mzd_copy dereferences a zero-width row.
  if (kh === 0 && mh > 0 && nh > 0)
    throw Object.assign(new Error('Segmentation fault'), { name: 'RuntimeError' });
  if (nh === 0 && m > 2 * mh) throw Object.assign(new Error('Aborted'), { name: 'RuntimeError' });
  if (mh === 0 || kh === 0 || nh === 0) return mzd_mul_m4rm(a, b);
  const A11 = window(a, 0, 0, mh, kh),
    A12 = window(a, 0, kh, mh, 2 * kh),
    A21 = window(a, mh, 0, 2 * mh, kh),
    A22 = window(a, mh, kh, 2 * mh, 2 * kh);
  let C11: mzd_t, C12: mzd_t, C21: mzd_t, C22: mzd_t, Wmk: mzd_t, Wkn: mzd_t;
  if (a === b) {
    Wkn = add(A22, A12);
    C21 = mzd_mul(Wkn, Wkn, cutoff);
    Wkn = add(A22, A21);
    C22 = mzd_mul(Wkn, Wkn, cutoff);
    Wkn = add(Wkn, A12);
    C11 = mzd_mul(Wkn, Wkn, cutoff);
    Wkn = add(Wkn, A11);
    C12 = add(mzd_mul(Wkn, A12, cutoff), C22);
    Wmk = mzd_mul(A12, A21, cutoff);
    C11 = add(C11, Wmk);
    C12 = add(C11, C12);
    C11 = add(C21, C11);
    C21 = add(C11, mzd_mul(A21, Wkn, cutoff));
    C22 = add(C22, C11);
    C11 = add(mzd_mul(A11, A11, cutoff), Wmk);
  } else {
    const B11 = window(b, 0, 0, kh, nh),
      B12 = window(b, 0, nh, kh, 2 * nh),
      B21 = window(b, kh, 0, 2 * kh, nh),
      B22 = window(b, kh, nh, 2 * kh, 2 * nh);
    Wkn = add(B22, B12);
    Wmk = add(A22, A12);
    C21 = mzd_mul(Wmk, Wkn, cutoff);
    Wmk = add(A22, A21);
    Wkn = add(B22, B21);
    C22 = mzd_mul(Wmk, Wkn, cutoff);
    Wkn = add(Wkn, B12);
    Wmk = add(Wmk, A12);
    C11 = mzd_mul(Wmk, Wkn, cutoff);
    Wmk = add(Wmk, A11);
    C12 = add(mzd_mul(Wmk, B12, cutoff), C22);
    Wmk = mzd_mul(A12, B21, cutoff);
    C11 = add(C11, Wmk);
    C12 = add(C11, C12);
    C11 = add(C21, C11);
    Wkn = add(Wkn, B11);
    C21 = add(C11, mzd_mul(A21, Wkn, cutoff));
    C22 = add(C22, C11);
    C11 = add(mzd_mul(A11, B11, cutoff), Wmk);
  }
  const rows = Array<bigint>(m).fill(0n);
  for (let i = 0; i < mh; i++) {
    rows[i] = C11.rows[i]! | (C12.rows[i]! << BigInt(nh));
    rows[mh + i] = C21.rows[i]! | (C22.rows[i]! << BigInt(nh));
  }
  const mm = 2 * mh,
    kk = 2 * kh,
    nn = 2 * nh;
  if (n > nn) {
    const tail = mzd_mul_m4rm(a, window(b, 0, nn, k, n));
    for (let i = 0; i < m; i++) rows[i] |= tail.rows[i]! << BigInt(nn);
  }
  if (m > mm) {
    const tail = mzd_mul_m4rm(window(a, mm, 0, m, k), window(b, 0, 0, k, nn));
    for (let i = 0; i < m - mm; i++) rows[mm + i] |= tail.rows[i]!;
  }
  if (k > kk) {
    const tail = mzd_mul_m4rm(window(a, 0, kk, mm, k), window(b, kk, 0, k, nn));
    for (let i = 0; i < mm; i++) rows[i] ^= tail.rows[i]!;
  }
  return { nrows: m, ncols: n, rows };
}
