/**
 * @module sage/schemes/hyperelliptic_curves/hyperelliptic_rational_field
 * @description Hyperelliptic curves over the rationals
 *
 * Port of: `sage/schemes/hyperelliptic_curves/hyperelliptic_rational_field.py`
 *
 * Both methods upstream delegate to machinery that is not part of this port:
 * `matrix_of_frobenius` to Monsky-Washnitzer cohomology, and `lseries` to
 * PARI's `lfun_genus2`.
 */

import { NotImplementedError } from '../../errors.js';
import type { RingElement } from '../../rings/polynomial/polynomial_element.js';
import { HyperellipticCurve_generic } from './hyperelliptic_generic.js';

export class HyperellipticCurve_rational_field<
  C extends RingElement,
> extends HyperellipticCurve_generic<C> {
  /** `hyperelliptic_rational_field.py:23-66` */
  matrix_of_frobenius(_p: unknown, _prec = 20): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: matrix_of_frobenius ' +
        '(sage.schemes.hyperelliptic_curves.monsky_washnitzer is not ported)'
    );
  }

  /** `hyperelliptic_rational_field.py:68-83` */
  lseries(_prec = 53): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: lseries (PARI lfun_genus2 is not available in parigp-ts)'
    );
  }
}
