import {
  zpFactorSquarefree,
  zpIsIrreducibleOverQ,
} from '../../../packages/sagemath-ts/src/rings/number_field/pari_nf.js';
import { getrand, setrand } from '../../../packages/parigp-ts/src/random.js';

export function nf_legacy_polynomial_factor(op: bigint, coefficients: bigint[]): string {
  const saved = getrand();
  try {
    setrand(1n);
    const result =
      op === 1n ? zpIsIrreducibleOverQ(coefficients) : zpFactorSquarefree(coefficients);
    return JSON.stringify([result, getrand()], (_, v) => (typeof v === 'bigint' ? String(v) : v));
  } finally {
    setrand(saved);
  }
}
