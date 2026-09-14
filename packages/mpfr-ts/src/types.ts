/** MPFR 4.2.1 decimal conversion adapter. Mutable values preserve native setters.
 * Supported: nearest-even rounding, decimal inputs of at most 4096 digits,
 * precisions 1..4096 and decimal exponent magnitude at most 4096.
 * Other rounding modes/bases and larger working domains fail explicitly.
 * @see Deviation: Native Real Literal Conversion
 */
export type mpfr_rnd_t = 'RNDN' | 'RNDZ' | 'RNDU' | 'RNDD' | 'RNDA' | 'RNDF';

export interface mpfr_t {
  precision: number;
  sign: 1 | -1;
  kind: 'finite' | 'zero' | 'inf' | 'nan';
  mantissa: bigint;
  /** Native exponent: value = sign * mantissa * 2^(exponent - precision). */
  exponent: number;
}
