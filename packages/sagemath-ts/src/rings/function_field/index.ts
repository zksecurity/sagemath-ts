/**
 * @module sage/rings/function_field
 * @description Function fields
 *
 * Port of: sage/rings/function_field/
 *
 * Implemented: rational function fields `k(x)` over an arbitrary constant
 * field, their elements, both maximal orders, the ideals of those orders,
 * places, divisors and Riemann-Roch spaces (Hess' algorithm 6.1 specialised to
 * degree one).  Finite extensions (`function_field_polymod.py` and everything
 * that depends on it) are not ported; see `function_field_polymod.ts`.
 */

export type {
  ConstantField,
  ConstantFieldElement,
} from './constant_field.js';
export {
  compare_constants,
  constant_field_cardinality,
  constant_field_characteristic,
  constant_field_element_list,
  constant_field_is_finite,
  divide_constants,
} from './constant_field.js';

export { FunctionField as FunctionFieldBase, is_FunctionField } from './function_field.js';
export {
  makeRationalFunctionField,
  RationalFunctionField,
  RationalFunctionField_char_zero,
  RationalFunctionField_global,
} from './function_field_rational.js';
export { FunctionField, FunctionFieldExtension } from './constructor.js';

export { FunctionFieldElement, is_FunctionFieldElement } from './element.js';
export {
  compare_polynomials,
  FunctionFieldElement_rational,
} from './element_rational.js';

export {
  FunctionFieldIdeal,
  FunctionFieldIdealInfinite,
  IdealMonoid,
} from './ideal.js';
export {
  FunctionFieldIdealInfinite_rational,
  FunctionFieldIdeal_rational,
} from './ideal_rational.js';

export {
  FunctionFieldMaximalOrder,
  FunctionFieldMaximalOrderInfinite,
  FunctionFieldOrder,
  FunctionFieldOrderInfinite,
  FunctionFieldOrder_base,
} from './order.js';
export {
  FunctionFieldMaximalOrderInfinite_rational,
  FunctionFieldMaximalOrder_rational,
} from './order_rational.js';

export { FunctionFieldPlace, PlaceSet } from './place.js';
export { FunctionFieldPlace_rational } from './place_rational.js';

export {
  divisor,
  DivisorGroup,
  FunctionFieldDivisor,
  prime_divisor,
} from './divisor.js';

export { FunctionFieldValuationRing } from './valuation_ring.js';

export { FunctionField_global, FunctionField_polymod } from './function_field_polymod.js';
