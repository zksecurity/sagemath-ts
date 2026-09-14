import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
import {
  NumberField,
  NumberFieldElement,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { QQ } from '../../../packages/sagemath-ts/src/rings/rational_field.js';
export function nf_base_trace_norm(
  op: bigint,
  nn: bigint,
  cs: bigint[],
  d: bigint,
  mode: bigint
): string {
  const n = Number(nn),
    K = new NumberField(
      RationalPolynomial.fromBigInts([-2n, ...Array(n - 1).fill(0n), 1n]),
      'a',
      undefined,
      false
    ),
    B = new NumberField(RationalPolynomial.fromBigInts([-3n, 1n]), 'b', undefined, false),
    a = K.__call__(cs.map((c) => new Rational(c, d)));
  const args =
    mode === 0n
      ? []
      : mode === 1n
        ? [null]
        : mode === 2n
          ? [QQ]
          : mode === 3n
            ? [K]
            : mode === 4n
              ? [B]
              : mode === 5n
                ? [
                    new NumberField(
                      RationalPolynomial.fromBigInts([-3n, 0n, 0n, 0n, 0n, 0n, 0n, 1n]),
                      'b',
                      undefined,
                      false
                    ),
                  ]
                : mode === 6n
                  ? [QQ, QQ]
                  : mode === 7n
                    ? [4n]
                    : mode === 8n
                      ? [0n]
                      : mode === 9n
                        ? ['bad']
                        : mode === 10n
                          ? [true]
                          : mode === 11n
                            ? [1.5]
                            : mode === 12n
                              ? [[]]
                              : mode === 13n
                                ? [{}]
                                : mode === 14n
                                  ? [new Integer(4n)]
                                  : [new Rational(1n, 2n)];
  const z = (op === 0n ? a.trace : a.norm).apply(a, args as any) as Rational | NumberFieldElement;
  return JSON.stringify({
    parent: z instanceof Rational ? 'QQ' : z.parent() === K ? 'K' : 'B',
    coefficients: z instanceof Rational ? [String(z)] : z.list().map(String),
  });
}
