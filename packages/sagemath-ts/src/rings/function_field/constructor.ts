/**
 * @module sage/rings/function_field/constructor
 * @description Factories for function fields
 *
 * Port of: sage/rings/function_field/constructor.py
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import type { ConstantField, ConstantFieldElement } from './constant_field.js';
import { makeRationalFunctionField } from './function_field_rational.js';
import type { RationalFunctionField } from './function_field_rational.js';

/**
 * Create a function field.
 *
 * INPUT:
 *
 * - ``F`` -- field
 * - ``names`` -- name of the generator
 *
 * ```typescript
 * const K = FunctionField(GF(7n), 'x');
 * K.toString(); // 'Rational function field in x over Finite Field of size 7'
 * ```
 *
 * @see Reference: sage/rings/function_field/constructor.py:115 (FunctionField)
 */
export function FunctionField<C extends ConstantFieldElement>(
  F: ConstantField<C>,
  names: string | [string]
): RationalFunctionField<C> {
  if (names === null || names === undefined) {
    throw new ValueError('variable name must be specified');
  }
  return makeRationalFunctionField(F, names);
}

/**
 * Create a function field defined as an extension of another function field by
 * adjoining a root of a univariate polynomial.
 *
 * @see Reference: sage/rings/function_field/constructor.py:216 (FunctionFieldExtension)
 */
export function FunctionFieldExtension(_polynomial: unknown, _names?: string | [string]): never {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: FunctionFieldExtension ' +
      '(sage/rings/function_field/function_field_polymod.py)'
  );
}
