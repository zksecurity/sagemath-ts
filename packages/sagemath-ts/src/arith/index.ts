/**
 * @module sage/arith
 * @description Arithmetic functions
 *
 * Port of: sage/arith/
 */

// Re-export IntegerLike and RationalLike types for use in arithmetic functions
export { type IntegerLike, type RationalLike, toBigInt, toRational } from '../types/coercion.js';

export * from './misc.js';
