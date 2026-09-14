import { _fmpz_poly_pseudo_divrem_divconquer as dc } from '../../../packages/flint-ts/src/fmpz_poly/pseudo_divrem_divconquer.js';
import { _fmpz_poly_pseudo_divrem_basecase as bc } from '../../../packages/flint-ts/src/fmpz_poly/pseudo_divrem_basecase.js';
import { _fmpz_poly_pseudo_rem as rem } from '../../../packages/flint-ts/src/fmpz_poly/pseudo_rem.js';
import { _fmpq_poly_rem as qr } from '../../../packages/flint-ts/src/fmpq_poly/rem.js';
const trim = (x: bigint[]) => {
  while (x.length && x[x.length - 1] === 0n) x.pop();
  return x;
};

export function flint_remainder(
  op: bigint,
  A: bigint[],
  B: bigint[],
  da: bigint,
  db: bigint
): string {
  const r = op === 0n ? dc(A, B) : op === 1n ? bc(A, B) : op === 2n ? rem(A, B) : qr(A, da, B, db);
  return JSON.stringify(r.map((x) => (Array.isArray(x) ? trim(x).map(String) : String(x))));
}
