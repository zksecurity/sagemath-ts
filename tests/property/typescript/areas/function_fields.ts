/**
 * sagemath-ts side of the live rational-function-field differential area.
 */

import { GF } from '../../../../packages/sagemath-ts/src/rings/finite_rings/index.js';
import type {
  ConstantField,
  ConstantFieldElement,
} from '../../../../packages/sagemath-ts/src/rings/function_field/constant_field.js';
import type { FunctionFieldElement_rational } from '../../../../packages/sagemath-ts/src/rings/function_field/element_rational.js';
import type { RationalFunctionField_global } from '../../../../packages/sagemath-ts/src/rings/function_field/function_field_rational.js';
import { FunctionField } from '../../../../packages/sagemath-ts/src/rings/function_field/index.js';

type CE = ConstantFieldElement;
type FF = FunctionFieldElement_rational<CE>;
const bool = (value: boolean): string => (value ? 'True' : 'False');

function field(p: bigint): RationalFunctionField_global<CE> {
  return FunctionField(
    GF(p) as unknown as ConstantField<CE>,
    'x'
  ) as RationalFunctionField_global<CE>;
}

function element(K: RationalFunctionField_global<CE>, coeffs: bigint[]): FF {
  const k = K.constant_base_field();
  return K.__call__(coeffs.map((c) => k.__call__(c))) as FF;
}

const ints = (xs: Array<{ value?: bigint; toString(): string }>): string =>
  xs.map((x) => (typeof x.value === 'bigint' ? x.value : x.toString())).join(',');

function parts(f: FF): string {
  return `${ints(f.numerator().coeffs)}/${ints(f.denominator().coeffs)}`;
}

function ff_arithmetic(
  p: bigint,
  numerator: bigint[],
  denominator: bigint[],
  exponent: bigint
): string {
  const K = field(p);
  const f = element(K, numerator).div(element(K, denominator));
  const g = f.pow(exponent);
  return `f=${parts(f)} pow=${parts(g)} degree=${f.degree()} vx=${f.valuation(
    K.gen()
  )} square=${bool(f.is_square())}`;
}

function ff_factor(p: bigint, numerator: bigint[], denominator: bigint[]): string {
  const K = field(p);
  const f = element(K, numerator).div(element(K, denominator));
  const F = f.factor();
  const factors = F.factors.map(([q, e]) => `${parts(q)}:${e}`).join(';');
  return `unit=${parts(F.unit)} factors=${factors}`;
}

function ff_places(p: bigint, degree: bigint): string {
  return field(p).places(Number(degree)).map(String).join(';');
}

export const functions = {
  ff_arithmetic,
  ff_factor,
  ff_places,
};
