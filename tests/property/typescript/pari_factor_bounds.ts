import { rationalPair } from '../../../packages/parigp-ts/src/_rational_polynomial.js';
import {
  Mignotte_bound,
  Beauzamy_bound,
  factor_bound,
  root_bound,
} from '../../../packages/parigp-ts/src/QX_factor.js';
import { vecbinomial } from '../../../packages/parigp-ts/src/bibli2.js';
import { ZX_Z_eval } from '../../../packages/parigp-ts/src/ZX.js';
import { ceil_safe } from '../../../packages/parigp-ts/src/gen3.js';
import { powruhalf } from '../../../packages/parigp-ts/src/trans1.js';
import { cmpir, cmpri } from '../../../packages/parigp-ts/src/kernel/none/level1.js';
import { itor, shiftr, type MpReal, type MpComplex } from '../../../packages/parigp-ts/src/qfb.js';
export function pari_factor_bounds(
  op: bigint,
  p: bigint,
  shift: bigint,
  n: bigint,
  a: bigint | bigint[]
): string {
  let result: unknown;
  if (op === 0n)
    result = ceil_safe(
      n === 0n
        ? shiftr(itor(a as bigint, Number(p)), Number(shift))
        : (a as bigint | [bigint, bigint])
    );
  else if (op === 1n) result = powruhalf(shiftr(itor(a as bigint, Number(p)), Number(shift)), n);
  else if (op === 6n) result = vecbinomial(Number(n));
  else if (op === 10n) {
    const [a0, b0] = a as bigint[];
    result = rationalPair(a0!, b0!);
  } else if (op === 7n || op === 8n) {
    const [x, y] = a as bigint[],
      r = shiftr(itor(y!, Number(p)), Number(shift));
    result = op === 7n ? cmpir(x!, r) : cmpri(r, x!);
  } else if (op === 9n) {
    const [y, ...coefficients] = a as bigint[];
    result = ZX_Z_eval(coefficients, y!);
  } else
    result =
      op === 2n
        ? Mignotte_bound(a as bigint[])
        : op === 3n
          ? Beauzamy_bound(a as bigint[])
          : op === 4n
            ? factor_bound(a as bigint[])
            : root_bound(a as bigint[]);
  const encode = (x: unknown): unknown =>
    typeof x === 'bigint' || typeof x === 'number'
      ? String(x)
      : Array.isArray(x)
        ? x.map(encode)
        : x && typeof x === 'object'
          ? 're' in x
            ? [encode((x as MpComplex).re), encode((x as MpComplex).im)]
            : [(x as MpReal).s, String((x as MpReal).e), String((x as MpReal).m), (x as MpReal).p]
          : x;
  return JSON.stringify(encode(result));
}
