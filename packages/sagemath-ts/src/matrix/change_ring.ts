import { fmpz_get_d } from '@sagemath-ts/flint-ts';
import { RDF, type RealDoubleElement } from '../rings/real_double.js';
import { Matrix } from './matrix_generic.js';
import type { IntegerMatrix } from './matrix_integer.js';

/** sage/matrix/change_ring.pyx: integer_to_real_double_dense.
 * @see Deviation: Polynomial Matrix Evaluation Actions
 */
export function integer_to_real_double_dense(a: IntegerMatrix): Matrix<RealDoubleElement> {
  return new Matrix(
    RDF,
    a.nrows,
    a.ncols,
    Array.from({ length: a.nrows }, (_, i) =>
      Array.from({ length: a.ncols }, (_, j) => RDF.__call__(fmpz_get_d(a.get(i, j).value)))
    )
  );
}
