import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
export function nf_trace_norm(op: bigint, cs: bigint[], flat: bigint[], d: bigint): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a', undefined, false),
    a = K.__call__(flat.map((c) => new Rational(c, d)));
  return String(
    op === 0n
      ? a.trace()
      : op === 1n
        ? a.norm()
        : op === 2n
          ? a.absolute_norm()
          : op === 3n
            ? a.absolute_trace()
            : op === 4n
              ? a.relative_norm()
              : a.relative_trace()
  );
}
