import { numberofconjugates } from '../../../packages/parigp-ts/src/galconj.js';
import { Flx_deriv, Flx_is_squarefree } from '../../../packages/parigp-ts/src/Flx.js';
import { FpX_red } from '../../../packages/parigp-ts/src/ffinit.js';
import { getrand, setrand } from '../../../packages/parigp-ts/src/random.js';
import { numberofconjugates as legacyCount } from '../../../packages/sagemath-ts/src/rings/number_field/pari_nf.js';

export function nf_conjugate_bound(
  op: bigint,
  p: bigint,
  start: bigint,
  coefficients: bigint[]
): string {
  const saved = getrand();
  try {
    setrand(1n);
    const value =
      op === 0n
        ? numberofconjugates(coefficients, start)
        : op === 1n
          ? legacyCount(coefficients, start)
          : op === 2n
            ? Flx_deriv(FpX_red(coefficients, p), p)
            : Flx_is_squarefree(FpX_red(coefficients, p), p)
              ? 1
              : 0;
    return JSON.stringify([value, getrand()], (_, v) =>
      typeof v === 'bigint' || typeof v === 'number' ? String(v) : v
    );
  } finally {
    setrand(saved);
  }
}
