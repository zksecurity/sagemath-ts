/**
 * @module sage/rings/function_field/function_field_polymod
 * @description Function fields: extension fields (NOT YET PORTED)
 *
 * Port of: sage/rings/function_field/function_field_polymod.py
 *
 * Finite extensions `L = K[y]/(f)` of a rational function field are not ported
 * yet.  Everything they need — `order_polymod.py` (maximal order via round-two /
 * Hermite forms), `ideal_polymod.py`, `place_polymod.py`, `differential.py`,
 * `derivations_polymod.py` and the general (degree `n > 1`) weak-Popov step of
 * Hess' algorithm 6.1 in `divisor.py` — is likewise missing.  Every entry point
 * below throws `NotImplementedError` naming the upstream file that has to be
 * ported.
 */

import { NotImplementedError } from '../../errors.js';

/**
 * @see Reference: sage/rings/function_field/function_field_polymod.py:71 (FunctionField_polymod)
 */
export function FunctionField_polymod(_polynomial: unknown, _names?: string | [string]): never {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: FunctionField_polymod ' +
      '(sage/rings/function_field/function_field_polymod.py)'
  );
}

/**
 * @see Reference: sage/rings/function_field/function_field_polymod.py (FunctionField_global)
 */
export function FunctionField_global(_polynomial: unknown, _names?: string | [string]): never {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: FunctionField_global ' +
      '(sage/rings/function_field/function_field_polymod.py)'
  );
}
