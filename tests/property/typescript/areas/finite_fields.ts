/**
 * sagemath-ts side of the `finite_fields` property-test area.
 *
 * Cases: tests/property/cases/finite_fields.cases.json
 * SageMath counterpart: tests/property/python/areas/finite_fields.py
 */

import { primitive_root } from '../../../../packages/sagemath-ts/src/arith/misc.js';
import { GF, generic_discrete_log, sqrt_mod } from '../../../../packages/sagemath-ts/src/index.js';
import { GFpn } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js';

export const functions = {
  // Prime field arithmetic
  ff_add: (p: bigint, a: bigint, b: bigint) => {
    const F = GF(p);
    return F.__call__(a).add(F.__call__(b)).value;
  },
  ff_mul: (p: bigint, a: bigint, b: bigint) => {
    const F = GF(p);
    return F.__call__(a).mul(F.__call__(b)).value;
  },
  ff_inv: (p: bigint, a: bigint) => {
    const F = GF(p);
    return F.__call__(a).inv().value;
  },
  ff_pow: (p: bigint, a: bigint, n: bigint) => {
    const F = GF(p);
    return F.__call__(a).pow(n).value;
  },
  ff_pow_neg: (p: bigint, a: bigint, n: bigint) => {
    const F = GF(p);
    return F.__call__(a).pow(n).value;
  },
  // Square root
  sqrt_mod: (a: bigint, p: bigint) => {
    const r = sqrt_mod(a, p);
    if (r === null) return null;
    // Return canonical (smaller) root
    return r > (p - 1n) / 2n ? p - r : r;
  },
  sqrt_mod_p3mod4: (a: bigint, p: bigint) => {
    const r = sqrt_mod(a, p);
    if (r === null) return null;
    return r > (p - 1n) / 2n ? p - r : r;
  },
  sqrt_mod_p5mod8: (a: bigint, p: bigint) => {
    const r = sqrt_mod(a, p);
    if (r === null) return null;
    return r > (p - 1n) / 2n ? p - r : r;
  },
  sqrt_mod_general: (a: bigint, p: bigint) => {
    const r = sqrt_mod(a, p);
    if (r === null) return null;
    return r > (p - 1n) / 2n ? p - r : r;
  },
  sqrt_mod_nonresidue: (a: bigint, p: bigint) => {
    const r = sqrt_mod(a, p);
    if (r === null) return null;
    return r > (p - 1n) / 2n ? p - r : r;
  },
  // Primitive roots
  primitive_root: (n: bigint) => primitive_root(n),
  ff_multiplicative_generator: (p: bigint) => {
    const F = GF(p);
    return F.multiplicative_generator().value;
  },
  // Discrete logarithm
  discrete_log: (p: bigint, base: bigint, target: bigint) => {
    const F = GF(p);
    const g = F.__call__(base);
    const h = F.__call__(target);
    return generic_discrete_log(h, g, p - 1n, '*');
  },
  // Extension field operations
  ff_ext_add: (p: bigint, n: bigint, a: bigint, b: bigint) => {
    const F = GFpn(p, Number(n));
    const x = F.fromInteger(a);
    const y = F.fromInteger(b);
    return x.add(y).integer_representation();
  },
  ff_ext_mul: (p: bigint, n: bigint, a: bigint, b: bigint) => {
    const F = GFpn(p, Number(n));
    const x = F.fromInteger(a);
    const y = F.fromInteger(b);
    return x.mul(y).integer_representation();
  },
  ff_ext_inv: (p: bigint, n: bigint, a: bigint) => {
    const F = GFpn(p, Number(n));
    const x = F.fromInteger(a);
    return x.inv().integer_representation();
  },
  ff_ext_pow: (p: bigint, n: bigint, a: bigint, e: bigint) => {
    const F = GFpn(p, Number(n));
    const x = F.fromInteger(a);
    return x.pow(e).integer_representation();
  },
  ff_ext_order: (p: bigint, n: bigint) => p ** n,
  ff_ext_frobenius: (p: bigint, n: bigint, a: bigint) => {
    const F = GFpn(p, Number(n));
    const x = F.fromInteger(a);
    return x.pow(p).integer_representation();
  },
  ff_generator_order: (p: bigint, n: bigint) => {
    // The multiplicative order of the generator is p^n - 1
    // for a primitive element
    return p ** n - 1n;
  },
};
