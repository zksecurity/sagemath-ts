/** @see Deviation: Polynomial Roots and Truncated Series */
/** Shared finite-ring coefficient root algorithms from element_base.pyx/integer_mod.pyx. */
import { crt, factor, gcd, inverse_mod, xgcd } from '../../arith/misc.js';
import { AssertionError, ValueError } from '../../errors.js';
import { discrete_log, type GroupElement } from '../../groups/generic.js';
import type { RingElement } from '../polynomial/polynomial_element.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import { FiniteFieldElement, PrimeFieldElement } from './finite_field_extension.js';
import { FiniteFieldElement as LegacyPrimeElement } from './finite_field_prime.js';
import { GF2Element } from './gf2.js';
import { IntegerMod } from './integer_mod.js';
import { IntegerModRing } from './integer_mod_ring.js';
import { pAdicGeneric } from '../padics/padic_generic.js';

type RootElement = {
  pow(n: bigint): RootElement;
  mul(b: RootElement): RootElement;
  add(b: RootElement): RootElement;
  isZero(): boolean;
  isOne(): boolean;
};

/** Positive coefficient roots used by Polynomial.nth_root and its Newton series. */
export function _coefficient_nth_root(c: RingElement, n: bigint): RingElement {
  if (n <= 0n) throw new ValueError('root exponent must be positive');
  if (c instanceof GF2Element) return c;
  if (c instanceof IntegerMod && !c.parent.is_field?.()) {
    const factors = factor(c.modulus);
    if (!factors.length) return c;
    let value = 0n,
      modulus = 1n;
    for (const [p, k] of factors) {
      const pk = p ** BigInt(k);
      const component = new IntegerModRing(pk).__call__(c.value);
      let root: IntegerMod;
      if (component.isZero()) root = component;
      else if (k === 1n) root = _coefficient_nth_root(component, n) as IntegerMod;
      else {
        let valuation = 0n,
          unit = component.value;
        while (unit % p === 0n) {
          unit /= p;
          valuation++;
        }
        if (valuation % n) throw new ValueError('no nth root');
        if (valuation > 0n) {
          const smaller = new IntegerModRing(p ** (k - valuation)).__call__(unit);
          const lifted = _coefficient_nth_root(smaller, n) as IntegerMod;
          root = component.parent.__call__(p ** (valuation / n) * lifted.value);
        } else {
          // integer_mod.pyx uses fixed-modulus p-adic log/exp for prime powers.
          // The quotient of two fixed-modulus lifts follows FM_template._quo_rem:
          // shift both representatives, then invert the divisor's unit modulo p^k.
          const R = new pAdicGeneric(p, Number(k));
          let sign = 1n,
            modp: bigint;
          if (p === 2n) {
            if (unit % 4n === 3n) {
              if (n % 2n === 0n) throw new ValueError('no nth root');
              sign = -1n;
              unit = pk - unit;
            } else if (n % 2n === 0n && k > 2n && unit % 8n === 5n)
              throw new ValueError('no nth root');
            if (k === 2n) {
              root = component.parent.__call__(sign);
              value = crt(value, root.value, modulus, pk);
              modulus *= pk;
              continue;
            }
            modp = unit % 8n;
          } else {
            modp = (_coefficient_nth_root(new IntegerModRing(p).__call__(unit), n) as IntegerMod)
              .value;
            const teich = R.teichmuller(unit % p).lift() as bigint;
            unit = (unit * inverse_mod(teich, pk)) % pk;
          }
          const logarithm = R.__call__(unit).log({ aprec: Number(k) }),
            logLift = logarithm.lift() as bigint;
          let nval = 0n,
            nu = n;
          while (nu % p === 0n) {
            nu /= p;
            nval++;
          }
          const logval = Number(logarithm.valuation());
          if (Number(nval) >= logval + (p === 2n ? -1 : 0)) {
            if (unit !== 1n) throw new ValueError('no nth root');
            root = component.parent.__call__(modp);
          } else {
            const ppart = p ** nval,
              right = (n % pk) / ppart;
            const quotient = ((logLift / ppart) * inverse_mod(right, pk)) % pk;
            const exponential = R.__call__(quotient).exp().lift() as bigint;
            const teich = R.teichmuller(modp).lift() as bigint;
            root = component.parent.__call__(sign * teich * exponential);
          }
        }
      }
      value = crt(value, root.value, modulus, pk);
      modulus *= pk;
    }
    return c.parent.__call__(value);
  }
  if (
    !(
      c instanceof IntegerMod ||
      c instanceof PrimeFieldElement ||
      c instanceof LegacyPrimeElement ||
      c instanceof FiniteFieldElement
    )
  )
    return (c as unknown as { nth_root(n: bigint): RingElement }).nth_root(n);
  if (c.isZero()) return c;
  const extension = c instanceof FiniteFieldElement;
  const q = extension
    ? c.parent.order
    : c instanceof IntegerMod
      ? c.modulus
      : c instanceof PrimeFieldElement
        ? c.parent.characteristic
        : c.p;
  const parent = c.parent as unknown as {
    __call__(x: unknown): RootElement;
    gen(): RootElement;
    multiplicative_generator(): RootElement;
    [Symbol.iterator](): Iterator<RootElement>;
  };
  const zeta = (order: bigint): RootElement => {
    // IntegerModRing retains the general Rings category's zeta protocol:
    // its first linear cyclotomic factor can select a different primitive root
    // from FiniteField.zeta even when the modulus is prime.
    if (c instanceof IntegerMod) {
      if (order === 1n) return parent.__call__(1n);
      if (order === 2n) return parent.__call__(-1n);
      const R = new PolynomialRing(c.parent, 'x');
      for (const [f] of R.cyclotomic_polynomial(Number(order)).factor())
        if (f.degree() === 1) return f.getCoeff(0).neg() as unknown as RootElement;
      throw new ValueError(`no ${order}th root of unity in ${c.parent}`);
    }
    const coorder = (q - 1n) / order;
    if (coorder <= 500000n) return parent.multiplicative_generator().pow(coorder);
    const primes = factor(order).map(([p]) => p);
    const gen = c instanceof IntegerMod ? parent.multiplicative_generator() : parent.gen();
    if (extension) {
      for (const x of parent) {
        const a = gen.add(x).pow(coorder);
        if (!a.isZero() && primes.every((p) => !a.pow(order / p).isOne())) return a;
      }
    } else {
      for (let i = 0n; i < q; i++) {
        const a = gen.add(parent.__call__(i)).pow(coorder);
        if (!a.isZero() && primes.every((p) => !a.pow(order / p).isOne())) return a;
      }
    }
    throw new AssertionError('no element found');
  };
  let self = c as unknown as RootElement,
    d = gcd(n, q - 1n);
  if (self.isOne()) return (d === 1n ? self : zeta(d)) as unknown as RingElement;
  if (d === q - 1n) throw new ValueError('no nth root');
  const [, alpha] = xgcd(n, q - 1n);
  if (d === 1n) return self.pow(alpha) as unknown as RingElement;
  if (!self.pow((q - 1n) / d).isOne()) throw new ValueError('no nth root');
  self = self.pow(alpha);
  for (const [r, v0] of factor(d)) {
    const v = BigInt(v0);
    let k = 0n,
      h = q - 1n;
    while (h % r === 0n) {
      h /= r;
      k++;
    }
    const rv = r ** v,
      hinv = inverse_mod(-h, rv),
      x = (1n + h * hinv) / rv;
    if (k === v) self = self.pow(x);
    else {
      const gh = zeta(r ** k);
      const t = discrete_log(
        self.pow(h) as unknown as GroupElement,
        gh.pow(rv) as unknown as GroupElement,
        r ** (k - v),
        '*'
      );
      self = self.pow(x).mul(gh.pow(-hinv * t));
    }
  }
  return self as unknown as RingElement;
}
