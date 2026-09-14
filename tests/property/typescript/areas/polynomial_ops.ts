import { PariError } from '../../../../packages/parigp-ts/src/errors.js';
/**
 * sagemath-ts side of the `polynomial_ops` property-test area (ZZ[x] and GF(p)[x]).
 *
 * Cases: tests/property/cases/polynomial_ops.cases.json
 * SageMath counterpart: tests/property/python/areas/polynomial_ops.py
 */

import { Integer, ZZ } from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { PolynomialRing } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_ring.js';
import {
  polyDerivative,
  polyEval,
  polyFactor,
  polyGcd,
  polyIsIrreducible,
  polyMul,
  polyPow,
  polyQuoRem,
  polyRoots,
} from './_helpers.js';

// Use the public polynomial implementation; test-local GCD/Horner algorithms
// previously bypassed both production methods and concealed missing content.
const integerCoefficients = {
  zero: () => new Integer(0n),
  one: () => new Integer(1n),
  __call__: (x: unknown) => new Integer(ZZ.__call__(x as Parameters<typeof ZZ.__call__>[0])),
  is_field: () => false,
  toString: () => 'Integer Ring',
};
const ZZx = new PolynomialRing(integerCoefficients, 'x');
export const functions = {
  poly_gcd_zz: (a: bigint[], b: bigint[]) =>
    ZZx.__call__(a)
      .gcd(ZZx.__call__(b))
      .coeffs.map((c) => c.value),
  poly_eval_zz: (a: bigint[], x: bigint) => ZZx.__call__(a).evaluate(new Integer(x)).value,
  poly_gcd_ff: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) => polyGcd(p, coeffs1, coeffs2),
  poly_eval_ff: (p: bigint, coeffs: bigint[], x: bigint) => polyEval(p, coeffs, x),
  poly_factor_ff: (p: bigint, coeffs: bigint[]) => polyFactor(p, coeffs),
  poly_roots_ff: (p: bigint, coeffs: bigint[]) => polyRoots(p, coeffs),
  poly_derivative_ff: (p: bigint, coeffs: bigint[]) => polyDerivative(p, coeffs),
  poly_is_irreducible_ff: (p: bigint, coeffs: bigint[]) => polyIsIrreducible(p, coeffs),
  poly_mul_ff: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) => polyMul(p, coeffs1, coeffs2),
  poly_quo_rem_ff: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) =>
    polyQuoRem(p, coeffs1, coeffs2),
  poly_pow_ff: (p: bigint, coeffs: bigint[], n: bigint) => polyPow(p, coeffs, n),
};

import {
  _fmpz_poly_gcd,
  _fmpz_poly_gcd_subresultant,
  _fmpz_poly_gcd_heuristic,
  _fmpz_poly_gcd_modular,
} from '../../../../packages/flint-ts/src/fmpz_poly/gcd.js';
Object.assign(functions, {
  flint_poly_gcd: _fmpz_poly_gcd,
  flint_poly_gcd_subresultant: _fmpz_poly_gcd_subresultant,
  flint_poly_gcd_heuristic: _fmpz_poly_gcd_heuristic,
  flint_poly_gcd_modular: _fmpz_poly_gcd_modular,
});

import { _nmod_poly_gcd } from '../../../../packages/flint-ts/src/nmod_poly/gcd.js';
Object.assign(functions, {
  flint_nmod_poly_gcd: (p: bigint, a: bigint[], b: bigint[]) => {
    const g = _nmod_poly_gcd(a, b, p);
    if (!g.length) return g;
    // Normalize the raw remainder to compare with FLINT's public monic wrapper.
    let r = p,
      s = g[g.length - 1]!,
      x = 0n,
      y = 1n;
    while (s) {
      const q = r / s;
      [r, s] = [s, r - q * s];
      [x, y] = [y, x - q * y];
    }
    return g.map((c) => (((c * x) % p) + p) % p);
  },
});

import { _fmpz_poly_divrem } from '../../../../packages/flint-ts/src/fmpz_poly/divrem.js';
Object.assign(functions, {
  flint_poly_divrem: (a: bigint[], b: bigint[], exact: bigint) => {
    const result = exact ? _fmpz_poly_divrem(a, b, true) : _fmpz_poly_divrem(a, b);
    return JSON.stringify(result?.map((poly) => poly.map(String)) ?? null);
  },
});

import { _nmod_poly_divrem } from '../../../../packages/flint-ts/src/nmod_poly/divrem.js';
Object.assign(functions, {
  flint_nmod_poly_divrem: (p: bigint, a: bigint[], b: bigint[]) =>
    JSON.stringify(_nmod_poly_divrem(a, b, p).map((poly) => poly.map(String))),
});

import { _fmpq_poly_gcd } from '../../../../packages/flint-ts/src/fmpq_poly/gcd.js';
Object.assign(functions, {
  flint_fmpq_poly_gcd: (a: bigint[], b: bigint[]) => {
    const [g, den] = _fmpq_poly_gcd(a, b);
    return JSON.stringify([g.map(String), String(den)]);
  },
});

import { _nmod_poly_xgcd } from '../../../../packages/flint-ts/src/nmod_poly/xgcd.js';
Object.assign(functions, {
  flint_nmod_poly_xgcd: (n: bigint, a: bigint[], b: bigint[]) => {
    try {
      return JSON.stringify({ value: _nmod_poly_xgcd(a, b, n).map((poly) => poly.map(String)) });
    } catch (e) {
      if (e instanceof RangeError) return JSON.stringify({ error: 'noninvertible' });
      throw e;
    }
  },
});

import { _fmpz_poly_resultant } from '../../../../packages/flint-ts/src/fmpz_poly/resultant.js';
import { _nmod_poly_resultant } from '../../../../packages/flint-ts/src/nmod_poly/resultant.js';
Object.assign(functions, {
  flint_poly_resultant: (a: bigint[], b: bigint[]) =>
    JSON.stringify({ value: String(_fmpz_poly_resultant(a, b)) }),
  flint_nmod_poly_resultant: (n: bigint, a: bigint[], b: bigint[], same: bigint) => {
    try {
      return JSON.stringify({ value: String(_nmod_poly_resultant(a, same ? a : b, n)) });
    } catch (e) {
      if (e instanceof RangeError) return JSON.stringify({ error: 'noninvertible' });
      throw e;
    }
  },
});

import { _fmpz_poly_xgcd } from '../../../../packages/flint-ts/src/fmpz_poly/xgcd.js';
import { _fmpq_poly_xgcd } from '../../../../packages/flint-ts/src/fmpq_poly/xgcd.js';
Object.assign(functions, {
  flint_poly_xgcd: (a: bigint[], b: bigint[]) => {
    const [r, s, t] = _fmpz_poly_xgcd(a, b);
    return JSON.stringify([String(r), s.map(String), t.map(String)]);
  },
  flint_fmpq_poly_xgcd: (a: bigint[], denA: bigint, b: bigint[], denB: bigint) =>
    JSON.stringify(
      _fmpq_poly_xgcd(a, denA, b, denB).map(([poly, den]) => [poly.map(String), String(den)])
    ),
});

import { _fmpq_poly_resultant } from '../../../../packages/flint-ts/src/fmpq_poly/resultant.js';
Object.assign(functions, {
  flint_fmpq_poly_resultant: (a: bigint[], denA: bigint, b: bigint[], denB: bigint) =>
    JSON.stringify(_fmpq_poly_resultant(a, denA, b, denB).map(String)),
});

import { dbltor, rtodbl } from '../../../../packages/parigp-ts/src/kernel/none/mp_indep.js';
import { resultant as pariRealResultant } from '../../../../packages/parigp-ts/src/polarit2.js';
function realFromBits(bits: bigint): number {
  const v = new DataView(new ArrayBuffer(8));
  v.setBigUint64(0, bits, false);
  return v.getFloat64(0, false);
}
function realToBits(x: number): string {
  const v = new DataView(new ArrayBuffer(8));
  v.setFloat64(0, x, false);
  return String(v.getBigUint64(0, false));
}
function nativeRealResult(run: () => number): string {
  try {
    return JSON.stringify({ value: realToBits(run()) });
  } catch (e) {
    if (e instanceof PariError && e.message.includes('overflow'))
      return JSON.stringify({ error: 'overflow' });
    throw e;
  }
}
Object.assign(functions, {
  pari_double_roundtrip: (a: bigint) => nativeRealResult(() => rtodbl(dbltor(realFromBits(a)))),
  pari_real_resultant: (a: bigint[], b: bigint[]) =>
    nativeRealResult(() => pariRealResultant(a.map(realFromBits), b.map(realFromBits))),
});

import { addrr as nativeAddrr } from '../../../../packages/parigp-ts/src/kernel/none/add.js';
import { divrr as nativeDivrr } from '../../../../packages/parigp-ts/src/kernel/gmp/mp.js';
import { mulrr, negr, real_0_bit, type MpReal } from '../../../../packages/parigp-ts/src/qfb.js';
function nativeRealFrame(z: MpReal): string[] {
  let a = z.m,
    b = 1n;
  const k = z.e + 1 - z.p;
  if (z.s) {
    if (k >= 0) a <<= BigInt(k);
    else b <<= BigInt(-k);
    while (!(a & 1n) && b > 1n) {
      a >>= 1n;
      b >>= 1n;
    }
    if (z.s < 0) a = -a;
  }
  return [String(z.s), String(z.e), String(a), String(b)];
}
Object.assign(functions, {
  pari_real_arithmetic: (operation: bigint, a: bigint, b: bigint) => {
    const cv = (raw: bigint) => {
      const x = realFromBits(raw);
      return x === 0 ? real_0_bit(-53) : dbltor(x);
    };
    const x = cv(a),
      y = cv(b);
    try {
      const z =
        operation === 0n
          ? nativeAddrr(x, y)
          : operation === 1n
            ? nativeAddrr(x, negr(y))
            : operation === 2n
              ? mulrr(x, y)
              : operation === 3n
                ? nativeDivrr(x, y)
                : operation === 4n
                  ? nativeAddrr(nativeDivrr(x, y), mulrr(x, y))
                  : nativeDivrr(nativeAddrr(x, y), nativeAddrr(x, negr(y)));
      return JSON.stringify(nativeRealFrame(z));
    } catch (e) {
      if (e instanceof RangeError && e.message.includes('inverse'))
        return JSON.stringify({ error: 'inverse' });
      throw e;
    }
  },
});

import { _fmpz_poly_derivative } from '../../../../packages/flint-ts/src/fmpz_poly/derivative.js';
import { _fmpq_poly_derivative } from '../../../../packages/flint-ts/src/fmpq_poly/derivative.js';
Object.assign(functions, {
  flint_poly_derivative: (a: bigint[]) => JSON.stringify(_fmpz_poly_derivative(a).map(String)),
  flint_fmpq_poly_derivative: (a: bigint[], den: bigint) => {
    const [b, d] = _fmpq_poly_derivative(a, den);
    return JSON.stringify([b.map(String), String(d)]);
  },
});

import {
  _nmod_poly_mul,
  _nmod_poly_mul_classical,
  _nmod_poly_mul_KS,
  _nmod_poly_mul_KS2,
  _nmod_poly_mul_KS4,
} from '../../../../packages/flint-ts/src/nmod_poly/mul.js';
import { _nmod_poly_pow } from '../../../../packages/flint-ts/src/nmod_poly/pow.js';
Object.assign(functions, {
  flint_nmod_poly_product: (p: bigint, a: bigint[], b: bigint[], method: bigint, same: bigint) =>
    JSON.stringify(
      [
        _nmod_poly_mul,
        _nmod_poly_mul_classical,
        _nmod_poly_mul_KS,
        _nmod_poly_mul_KS2,
        _nmod_poly_mul_KS4,
        _nmod_poly_add,
        _nmod_poly_sub,
        (a: readonly bigint[], _b: readonly bigint[], p: bigint) => _nmod_poly_sub([], a, p),
      ][Number(method)]!(a, same ? a : b, p).map(String)
    ),
  flint_nmod_poly_power: (p: bigint, a: bigint[], e: bigint) =>
    JSON.stringify(_nmod_poly_pow(a, e, p).map(String)),
});

import { _nmod_poly_add } from '../../../../packages/flint-ts/src/nmod_poly/add.js';
import { _nmod_poly_sub } from '../../../../packages/flint-ts/src/nmod_poly/sub.js';

import {
  _fmpz_poly_pow,
  _fmpz_poly_pow_small,
  _fmpz_poly_pow_binomial,
  _fmpz_poly_pow_multinomial,
  _fmpz_poly_pow_binexp,
  _fmpq_poly_pow,
} from '../../../../packages/flint-ts/src/index.js';
function nativePowerFrame(a: bigint[]): [number, string] {
  return [
    a.length,
    new Bun.CryptoHasher('sha256').update(a.map((c) => c.toString(16)).join(',')).digest('hex'),
  ];
}
Object.assign(functions, {
  flint_fmpz_poly_power: (a: bigint[], e: bigint, method: bigint) =>
    JSON.stringify(
      nativePowerFrame(
        [
          _fmpz_poly_pow,
          _fmpz_poly_pow_binomial,
          _fmpz_poly_pow_multinomial,
          _fmpz_poly_pow_binexp,
          _fmpz_poly_pow_small,
        ][Number(method)]!(a, e)
      )
    ),
  flint_fmpq_poly_power: (a: bigint[], d: bigint, e: bigint) => {
    const [b, den] = _fmpq_poly_pow(a, d, e);
    return JSON.stringify([nativePowerFrame(b), String(den)]);
  },
});

import {
  GF2X as NativeGF2X,
  GF2X_power,
  ZZ_pX_power,
  ZZ_pEX_power,
} from '../../../../packages/ntl-ts/src/index.js';
Object.assign(functions, {
  ntl_gf2x_power: (a: bigint, e: bigint) => GF2X_power(new NativeGF2X(a), e).rep().toString(16),
  ntl_gf2x_product: (a: bigint, b: bigint, same: bigint) => {
    const x = new NativeGF2X(a);
    return x
      .mul(same ? x : new NativeGF2X(b))
      .rep()
      .toString(16);
  },
  ntl_zzpx_power: (p: bigint, a: bigint[], e: bigint) =>
    JSON.stringify(nativePowerFrame(ZZ_pX_power(a, e, p))),
  ntl_zzpex_power: (p: bigint, f: bigint[], a: bigint[], e: bigint) => {
    const coefficients = a.map((c) => {
      const r: bigint[] = [];
      while (c) {
        r.push(c % p);
        c /= p;
      }
      return r;
    });
    return JSON.stringify(ZZ_pEX_power(coefficients, e, f, p).map((row) => row.map(String)));
  },
});

import {
  _fmpz_poly_mullow,
  _fmpq_poly_mullow,
  _nmod_poly_mullow,
  _fmpz_poly_pow_trunc,
  _nmod_poly_pow_trunc,
  _fmpz_poly_inv_series,
  _fmpq_poly_inv_series_newton,
} from '../../../../packages/flint-ts/src/index.js';
functions.flint_polynomial_series = (
  kind: bigint,
  p: bigint,
  a: bigint[],
  da: bigint,
  b: bigint[],
  db: bigint,
  e: bigint,
  n: bigint
) => {
  const result: [bigint[], bigint] =
    kind === 0n
      ? [_fmpz_poly_mullow(a, b, Number(n)), 1n]
      : kind === 1n
        ? _fmpq_poly_mullow(a, da, b, db, Number(n))
        : kind === 2n
          ? [_nmod_poly_mullow(a, b, Number(n), p), 1n]
          : kind === 3n
            ? [_fmpz_poly_pow_trunc(a, e, Number(n)), 1n]
            : kind === 4n
              ? [_nmod_poly_pow_trunc(a, e, Number(n), p), 1n]
              : kind === 5n
                ? [_fmpz_poly_inv_series(a, Number(n)), 1n]
                : _fmpq_poly_inv_series_newton(a, da, Number(n));
  return JSON.stringify([nativePowerFrame(result[0]), String(result[1])]);
};

import { XGCD as denseXGCD } from '../../../../packages/ntl-ts/src/ZZ_pX1.js';
import { InvTrunc as extensionInvTrunc } from '../../../../packages/ntl-ts/src/ZZ_pEX.js';
functions.ntl_zzpx_dense_xgcd = (p: bigint, a: bigint[], b: bigint[]) =>
  JSON.stringify(denseXGCD(a, b, p).map((v) => v.map(String)));
functions.ntl_zzpex_inverse_series = (p: bigint, f: bigint[], a: bigint[], n: bigint) => {
  const coefficients = a.map((c) => {
    const r: bigint[] = [];
    while (c) {
      r.push(c % p);
      c /= p;
    }
    return r;
  });
  return JSON.stringify(extensionInvTrunc(coefficients, Number(n), f, p).map((v) => v.map(String)));
};

import { _fmpz_poly_mul, _fmpq_poly_mul } from '../../../../packages/flint-ts/src/index.js';
import { ZZ_pX_mul, ZZ_pEX_mul } from '../../../../packages/ntl-ts/src/index.js';
Object.assign(functions, {
  flint_full_product: (
    a: bigint[],
    da: bigint,
    b: bigint[],
    db: bigint,
    rational: bigint,
    same: bigint
  ) => {
    const [out, d] = rational
      ? _fmpq_poly_mul(a, da, same ? a : b, same ? da : db)
      : ([_fmpz_poly_mul(a, same ? a : b), 1n] as const);
    return JSON.stringify([nativePowerFrame(out), String(d)]);
  },
  ntl_zzpx_full_product: (p: bigint, a: bigint[], b: bigint[], same: bigint) =>
    JSON.stringify(nativePowerFrame(ZZ_pX_mul(a, same ? a : b, p))),
  ntl_zzpex_full_product: (p: bigint, f: bigint[], a: bigint[], b: bigint[], same: bigint) => {
    const pack = (v: bigint[]) =>
      v.map((c) => {
        const r: bigint[] = [];
        while (c) {
          r.push(c % p);
          c /= p;
        }
        return r;
      });
    const A = pack(a),
      B = same ? A : pack(b);
    return JSON.stringify(ZZ_pEX_mul(A, B, f, p).map((row) => row.map(String)));
  },
});

import {
  _nmod_poly_inv_series_newton,
  _nmod_poly_powmod_ui_binexp,
  _nmod_poly_powmod_fmpz_binexp_preinv,
  _nmod_poly_powmod_x_fmpz_preinv,
} from '../../../../packages/flint-ts/src/index.js';
import {
  ZZ_pEX_PowerMod,
  ZZ_pEX_PowerXMod,
  ZZ_pEX_XGCD,
} from '../../../../packages/ntl-ts/src/index.js';
Object.assign(functions, {
  flint_modular_power: (
    p: bigint,
    a: bigint[],
    m: bigint[],
    e: bigint,
    method: bigint,
    precision: bigint
  ) => {
    while (m.length && m[m.length - 1] === 0n) m = m.slice(0, -1);
    const inverse =
      method === 1n || method === 2n
        ? _nmod_poly_inv_series_newton(m.slice().reverse(), m.length, p)
        : [];
    if (method === 3n && precision === 0n) {
      try {
        _nmod_poly_inv_series_newton(a, 0, p);
      } catch (err) {
        if (err instanceof RangeError && err.message === 'invalid truncation length')
          return JSON.stringify({ invalid_precision: true });
        throw err;
      }
    }
    const h =
      method === 0n
        ? _nmod_poly_powmod_ui_binexp(a, e, m, p)
        : method === 1n
          ? _nmod_poly_powmod_fmpz_binexp_preinv(a, e, m, inverse, p)
          : method === 2n
            ? _nmod_poly_powmod_x_fmpz_preinv(e, m, inverse, p)
            : _nmod_poly_inv_series_newton(a, Number(precision), p);
    return JSON.stringify(nativePowerFrame(h));
  },
  ntl_extension_modular: (
    p: bigint,
    f: bigint[],
    a: bigint[],
    b: bigint[],
    e: bigint,
    method: bigint
  ) => {
    const pack = (v: bigint[]) =>
      v.map((c) => {
        const out: bigint[] = [];
        while (c) {
          out.push(c % p);
          c /= p;
        }
        return out;
      });
    const A = pack(a),
      B = pack(b),
      frame = (h: bigint[][]) => h.map((row) => row.map(String));
    try {
      return JSON.stringify({
        value:
          method === 2n
            ? ZZ_pEX_XGCD(A, B, f, p).map(frame)
            : frame(method === 0n ? ZZ_pEX_PowerMod(A, e, B, f, p) : ZZ_pEX_PowerXMod(e, B, f, p)),
      });
    } catch (err) {
      return JSON.stringify({ error: 'NTLError', message: (err as Error).message });
    }
  },
});

import { _nmod_poly_xgcd_kernels as nmodGcdKernels } from '../../../../packages/flint-ts/src/nmod_poly/gcd.js';
Object.assign(functions, {
  flint_nmod_gcd_product: (p: bigint, a: bigint[], b: bigint[], same: bigint) =>
    JSON.stringify(nmodGcdKernels.mul(a, same ? a : b, p).map(String)),
});

import { _fmpz_poly_evaluate_fmpz } from '../../../../packages/flint-ts/src/fmpz_poly/evaluate_fmpz.js';
import { _fmpz_poly_evaluate_fmpq } from '../../../../packages/flint-ts/src/fmpz_poly/evaluate_fmpq.js';
import {
  _fmpq_poly_evaluate_fmpz,
  _fmpq_poly_evaluate_fmpq,
} from '../../../../packages/flint-ts/src/fmpq_poly/evaluate.js';
import { _nmod_poly_evaluate_nmod } from '../../../../packages/flint-ts/src/nmod_poly/evaluate_nmod.js';
import { eval as extensionEvaluate } from '../../../../packages/ntl-ts/src/ZZ_pEX.js';
Object.assign(functions, {
  flint_polynomial_evaluation: (
    kind: bigint,
    a: bigint[],
    den: bigint,
    num: bigint,
    xden: bigint
  ) =>
    JSON.stringify(
      (kind === 0n
        ? [_fmpz_poly_evaluate_fmpz(a, num), 1n]
        : kind === 1n
          ? _fmpz_poly_evaluate_fmpq(a, num, xden)
          : kind === 2n
            ? _fmpq_poly_evaluate_fmpq(a, den, num, xden)
            : _fmpq_poly_evaluate_fmpz(a, den, num)
      ).map((v) => v.toString(16))
    ),
  flint_nmod_evaluation: (p: bigint, a: bigint[], x: bigint) =>
    String(_nmod_poly_evaluate_nmod(a, x, p)),
  ntl_extension_evaluation: (p: bigint, f: bigint[], a: bigint[], x: bigint) => {
    const digits = (n: bigint) => {
      const out: bigint[] = [];
      while (n) {
        out.push(n % p);
        n /= p;
      }
      return out;
    };
    return JSON.stringify(extensionEvaluate(a.map(digits), digits(x), f, p).map(String));
  },
});

import { _fmpz_poly_compose } from '../../../../packages/flint-ts/src/fmpz_poly/compose.js';
import { _fmpq_poly_compose } from '../../../../packages/flint-ts/src/fmpq_poly/compose.js';
import { _nmod_poly_compose } from '../../../../packages/flint-ts/src/nmod_poly/compose.js';
import { _fmpz_poly_taylor_shift } from '../../../../packages/flint-ts/src/fmpz_poly/taylor_shift.js';
Object.assign(functions, {
  flint_polynomial_composition: (
    kind: bigint,
    a: bigint[],
    da: bigint,
    b: bigint[],
    db: bigint,
    p: bigint
  ) => {
    const [coefficients, denominator] =
      kind === 0n
        ? ([_fmpz_poly_compose(a, b), 1n] as const)
        : kind === 1n
          ? _fmpq_poly_compose(a, da, b, db)
          : kind === 2n
            ? ([_nmod_poly_compose(a, b, p), 1n] as const)
            : ([_fmpz_poly_taylor_shift(a, b[0]!), 1n] as const);
    return JSON.stringify([coefficients.map((c) => c.toString(16)), denominator.toString(16)]);
  },
});

import { ZZ_pX_evaluate } from '../../../../packages/ntl-ts/src/index.js';
Object.assign(functions, {
  ntl_prime_evaluation: (p: bigint, a: bigint[], x: bigint) => String(ZZ_pX_evaluate(a, x, p)),
});

import { fmpz_get_d } from '../../../../packages/flint-ts/src/index.js';
functions.flint_integer_double = (value: bigint) => {
  const view = new DataView(new ArrayBuffer(8));
  view.setFloat64(0, fmpz_get_d(value));
  return String(view.getBigUint64(0));
};
