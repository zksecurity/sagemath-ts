/** Permanent source/oracle regressions for core arithmetic and coercion boundaries. */
import { crt } from '../../../../packages/sagemath-ts/src/arith/misc.js';
import { discrete_log } from '../../../../packages/sagemath-ts/src/groups/generic.js';
import {
  GF,
  GFpn,
} from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
import { FiniteFieldPrime } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_prime.js';
import { Integer } from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { PolynomialRing } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_ring.js';

function run(operation: string, args: number[]): unknown {
  const [a, b, c, d] = args;
  switch (operation) {
    case 'modulus': {
      const F = GFpn(BigInt(a), 2, [a - 1, 0]);
      const u = F.gen().sub(F.one());
      const v = F.gen().add(F.one());
      return [F.is_field(), u.isZero(), v.isZero(), u.mul(v).isZero()];
    }
    case 'field_add':
      return String(
        GF(BigInt(a))
          .one()
          .add(GF(BigInt(b)).one())
      );
    case 'field_div':
      return String(
        GF(BigInt(a))
          .one()
          .div(GF(BigInt(b)).__call__(2))
      );
    case 'legacy_field_eq':
      return new FiniteFieldPrime(BigInt(a)).one().eq(new FiniteFieldPrime(BigInt(b)).one());
    case 'legacy_field_add':
      return String(
        new FiniteFieldPrime(BigInt(a)).one().add(new FiniteFieldPrime(BigInt(b)).one())
      );
    case 'extension_eq':
      return GF(BigInt(a), 'a')
        .gen()
        .eq(GF(BigInt(b), 'b').gen());
    case 'extension_scalar_eq':
      return GF(BigInt(a)).__call__(b).eq(BigInt(b));
    case 'field_eq':
      return GF(BigInt(a))
        .one()
        .eq(GF(BigInt(b)).one());
    case 'frobenius':
      return String(GF(BigInt(a)).gen().frobenius(b).integer_representation());
    case 'constant_polynomial_eq':
      return new PolynomialRing(GF(BigInt(a)), 'x')
        .__call__(c)
        .eq(new PolynomialRing(GF(BigInt(b)), d === 0 ? 'x' : 'y').__call__(c));
    case 'polynomial_eq':
    case 'polynomial_sub': {
      const F = GF(BigInt(a));
      const x = new PolynomialRing(F, 'x').gen();
      const y = new PolynomialRing(F, 'y').gen();
      return operation === 'polynomial_eq' ? x.eq(y) : x.sub(y).isZero();
    }
    case 'squarefree':
    case 'factor_product':
    case 'squarefree_inseparable': {
      const F = GF(BigInt(a));
      const R = new PolynomialRing(F, 'x');
      const f = R.gen()
        .add(R.one())
        .pow(operation === 'squarefree_inseparable' ? a : 2)
        .scalar_mul(F.__call__(b));
      const factors = operation === 'factor_product' ? f.factor() : f.squarefree_decomposition();
      const product = factors.reduce((v, [g, e]) => v.mul(g.pow(e)), R.one());
      return product.eq(f);
    }
    case 'crt':
      return String(crt(BigInt(a), BigInt(b), BigInt(c), BigInt(d)));
    case 'additive_log': {
      const x = GF(BigInt(a)).gen();
      return String(discrete_log(x.add(x), x, BigInt(b), '+'));
    }
    case 'pow_error':
    case 'integer_pow':
      return String(new Integer(BigInt(a)).pow(BigInt(b)));
    case 'invert_error':
      return String(new Integer(BigInt(a)).__invert__());
    case 'exact_log':
      return String(new Integer(BigInt(a)).exact_log(BigInt(b)));
    default:
      throw new Error(`Unknown operation: ${operation}`);
  }
}
// Encode exceptions as ordinary results; the main harness otherwise treats any two errors alike.
const operations = [
  'invert_error',
  'pow_error',
  'legacy_field_eq',
  'legacy_field_add',
  'extension_eq',
  'extension_scalar_eq',
  'constant_polynomial_eq',
  'factor_product',
  'squarefree_inseparable',
  'modulus',
  'field_add',
  'field_div',
  'field_eq',
  'frobenius',
  'polynomial_eq',
  'polynomial_sub',
  'squarefree',
  'crt',
  'additive_log',
  'integer_pow',
  'exact_log',
];
export const functions = Object.fromEntries(
  operations.map((operation) => [
    operation,
    (...args: bigint[]) => {
      try {
        return JSON.stringify({ value: run(operation, args.map(Number)) });
      } catch (error) {
        return JSON.stringify({
          error: error instanceof Error ? error.name : String(error),
          ...(operation.endsWith('_error') && error instanceof Error
            ? { message: error.message }
            : {}),
        });
      }
    },
  ])
);
