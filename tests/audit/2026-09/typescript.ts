/** Standalone differential audit; intentionally exercises unresolved fidelity gaps. */
import { readFileSync } from 'node:fs';
import { crt } from '../../../packages/sagemath-ts/src/arith/misc.js';
import { discrete_log } from '../../../packages/sagemath-ts/src/groups/generic.js';
import {
  GF,
  GFpn,
} from '../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { PolynomialRing } from '../../../packages/sagemath-ts/src/rings/polynomial/polynomial_ring.js';

const cases = JSON.parse(readFileSync(new URL('./cases.json', import.meta.url), 'utf8'));
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
    case 'field_eq':
      return GF(BigInt(a))
        .one()
        .eq(GF(BigInt(b)).one());
    case 'frobenius':
      return String(GF(BigInt(a)).gen().frobenius(b).integer_representation());
    case 'polynomial_eq':
    case 'polynomial_sub': {
      const F = GF(BigInt(a));
      const x = new PolynomialRing(F, 'x').gen();
      const y = new PolynomialRing(F, 'y').gen();
      return operation === 'polynomial_eq' ? x.eq(y) : x.sub(y).isZero();
    }
    case 'squarefree': {
      const F = GF(BigInt(a));
      const R = new PolynomialRing(F, 'x');
      const f = R.gen().add(R.one()).pow(2).scalar_mul(F.__call__(b));
      const product = f.squarefree_decomposition().reduce((v, [g, e]) => v.mul(g.pow(e)), R.one());
      return product.eq(f);
    }
    case 'crt':
      return String(crt(BigInt(a), BigInt(b), BigInt(c), BigInt(d)));
    case 'additive_log': {
      const x = GF(BigInt(a)).gen();
      return String(discrete_log(x.add(x), x, BigInt(b), '+'));
    }
    case 'integer_pow':
      return String(new Integer(BigInt(a)).pow(BigInt(b)));
    case 'exact_log':
      return String(new Integer(BigInt(a)).exact_log(BigInt(b)));
    default:
      throw new Error(`Unknown operation: ${operation}`);
  }
}
console.log(
  JSON.stringify(
    cases.map(({ id, operation, args }) => {
      try {
        return { id, value: run(operation, args) };
      } catch (error) {
        return { id, error: error.name, message: error.message };
      }
    })
  )
);
