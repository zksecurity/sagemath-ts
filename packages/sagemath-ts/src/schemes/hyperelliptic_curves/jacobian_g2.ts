/**
 * @module sage/schemes/hyperelliptic_curves/jacobian_g2
 * @description Jacobian of a hyperelliptic curve of genus 2
 *
 * Port of: `sage/schemes/hyperelliptic_curves/jacobian_g2.py`
 */

import { NotImplementedError } from '../../errors.js';
import type { RingElement } from '../../rings/polynomial/polynomial_element.js';
import { HyperellipticJacobian_generic } from './jacobian_generic.js';

export class HyperellipticJacobian_g2<
  C extends RingElement,
> extends HyperellipticJacobian_generic<C> {
  /** `jacobian_g2.py:26-31` — needs `kummer_surface`, which is not ported. */
  kummer_surface(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: kummer_surface (sage.schemes.hyperelliptic_curves.kummer_surface is not ported)'
    );
  }
}
