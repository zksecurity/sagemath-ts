import { zz_pXModulus, PowerXMod } from '../../../packages/ntl-ts/src/lzz_pX.js';
import { NewDDF, SFCanZass1 } from '../../../packages/ntl-ts/src/lzz_pXFactoring.js';
export function ntl_word_ddf(op: bigint, p: bigint, f: bigint[], h: bigint[]): string {
  let x: unknown;
  if (op === 0n) x = NewDDF(f, h, p);
  else if (op === 1n) x = SFCanZass1(f, p);
  else if (op === 2n) x = NewDDF(f, PowerXMod(p, new zz_pXModulus(f, p)), p);
  else throw new Error('unknown NTL distinct-degree operation');
  return JSON.stringify(x, (_, v) =>
    typeof v === 'bigint' || typeof v === 'number' ? String(v) : v
  );
}
