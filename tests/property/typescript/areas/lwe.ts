/**
 * sagemath-ts side of the `lwe` property-test area.
 *
 * Cases: tests/property/cases/lwe.cases.json
 * SageMath counterpart: tests/property/python/areas/lwe.py
 */

import {
  LindnerPeikert,
  Regev,
  UniformNoiseLWE,
} from '../../../../packages/sagemath-ts/src/crypto/lwe.js';

export const functions = {
  regev_q: (n: bigint): bigint => new Regev(n).K.order,
  regev_dimension: (n: bigint): number => new Regev(n).n,
  lindner_peikert_m: (n: bigint): number | null => new LindnerPeikert(n).m,
  sample_vector_length: (n: bigint): number => new Regev(n).call()[0].length,
  uniform_noise_min_n: (n: bigint): string => {
    new UniformNoiseLWE(n);
    return 'valid';
  },
};
