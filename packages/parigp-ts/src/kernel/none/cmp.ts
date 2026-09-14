import type { MpReal } from '../../qfb.js';

/** PARI `cmprr` (`kernel/none/cmp.c`), including finite-accuracy zero comparisons. */
export function cmprr(x: MpReal, y: MpReal): number {
  if (!x.s) return !y.s || x.e >= y.e ? 0 : -y.s;
  if (!y.s) return y.e >= x.e ? 0 : x.s;
  if (x.s !== y.s) return x.s < y.s ? -1 : 1;
  if (x.e !== y.e) return x.e > y.e ? x.s : -x.s;
  // Native compares common words, then any nonzero trailing words. Packed
  // alignment performs the same comparison with work bounded by precision.
  const p = Math.max(x.p, y.p);
  const a = x.m << BigInt(p - x.p), b = y.m << BigInt(p - y.p);
  return a === b ? 0 : a > b ? x.s : -x.s;
}

/** PARI abscmprr, whose zero comparisons ignore the accuracy exponent. */
export function abscmprr(x: MpReal, y: MpReal): number {
  if (!x.s) return y.s ? -1 : 0;
  if (!y.s) return 1;
  return cmprr({ ...x, s: 1 }, { ...y, s: 1 });
}
