/**
 * @module sage/coding
 * @description Coding theory - error-correcting codes
 *
 * This module provides implementations of error-correcting codes,
 * including Reed-Solomon codes for ZK applications (FRI, STARKs),
 * BCH (Bose-Chaudhuri-Hocquenghem) codes, and Goppa codes for
 * post-quantum cryptography (McEliece cryptosystem).
 *
 * @example
 * ```typescript
 * import { ReedSolomonCode, BCHCode, GoppaCode, createFRIReedSolomonCode } from 'sagemath-ts/coding';
 * import { GF, GFExtended } from 'sagemath-ts';
 *
 * // Create a Reed-Solomon code over GF(17) with length 16 and dimension 8
 * const F = GF(17n);
 * const rs = new ReedSolomonCode(F, 16, 8);
 *
 * // Create a BCH(15, 7) code over GF(2)
 * const F2 = GF(2n);
 * const bch = new BCHCode(15, 7, F2);
 *
 * // Create a binary Goppa code for McEliece cryptosystem
 * const F64 = GFExtended(64);  // GF(2^6)
 * const goppa = createBinaryGoppaCode(g, F64);
 *
 * // Encode and decode with error correction
 * const message = bch.encode([...]);
 * const decoded = bch.decode(received);
 * ```
 */

export {
  ReedSolomonCode,
  DecodingError as RSDecodingError,
  createClassicalReedSolomonCode,
  createFRIReedSolomonCode,
  type FieldElement,
  type FiniteField,
} from './reed_solomon.js';

export {
  BCHCode,
  DecodingError,
  createBCHCode,
  createPrimitiveBCHCode,
  isReedSolomonCode,
} from './bch_code.js';

export {
  ReedMullerCode,
  PuncturedReedMullerCode,
  BinaryReedMullerCode,
  FirstOrderReedMullerCode,
  RepetitionCode,
} from './reed_muller_code.js';

export {
  GoppaCode,
  BinaryGoppaCode,
  DecodingError as GoppaDecodingError,
  createGoppaCode,
  createBinaryGoppaCode,
  type FieldElement as GoppaFieldElement,
  type FiniteField as GoppaFiniteField,
} from './goppa_code.js';
