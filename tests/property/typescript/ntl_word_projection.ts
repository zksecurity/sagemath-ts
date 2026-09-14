import {
  zz_pXModulus,
  zz_pXMultiplier,
  build,
  MulMod,
} from '../../../packages/ntl-ts/src/lzz_pX.js';
import {
  UpdateMap,
  ProjectPowers,
  zz_pXNewArgument,
  build as buildArgument,
} from '../../../packages/ntl-ts/src/lzz_pX1.js';
export function ntl_word_projection(
  op: bigint,
  p: bigint,
  param: bigint[],
  packed: bigint[]
): string {
  const w: bigint[][] = [];
  for (let i = 0, k = 1; i < Number(packed[0]); i++) {
    const n = Number(packed[k++]!);
    w.push(packed.slice(k, k + n));
    k += n;
  }
  const F = new zz_pXModulus(param[2] ? w[0]! : null, p),
    B = new zz_pXMultiplier(),
    H = new zz_pXNewArgument();
  const m = Number(param[0]),
    count = Number(param[1]);
  let x: unknown;
  if (op <= 2n || op === 5n) {
    if (param[3]) build(B, w[1]!, F);
  }
  const state = () => [B.b, B.UseFFT, B.val()];
  if (op === 0n) x = state();
  else if (op === 1n) x = MulMod(w[2]!, B, F);
  else if (op === 2n) x = UpdateMap(w[2]!, B, F);
  else if (op === 3n) x = ProjectPowers(w[2]!, count, w[3]!, F);
  else if (op === 4n) {
    if (param[3]) buildArgument(H, w[3]!, F, m);
    x = ProjectPowers(w[2]!, count, H, F);
  } else if (op === 5n) {
    let error: string | null = null;
    try {
      build(B, w[4]!, F);
    } catch (e) {
      error = (e as Error).message;
    }
    const result = (run: () => bigint[]) => {
      try {
        return [null, run()];
      } catch (e) {
        return [(e as Error).message, null];
      }
    };
    x = [error, state(), result(() => MulMod(w[2]!, B, F)), result(() => UpdateMap(w[2]!, B, F))];
  } else throw new Error('unknown NTL projection operation');
  return JSON.stringify(x, (_, v) =>
    typeof v === 'bigint' || typeof v === 'number' ? String(v) : v
  );
}
