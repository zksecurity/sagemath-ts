import {exp1r_abs} from '../../../../packages/parigp-ts/src/trans1.js';
import { abpq_init, abpq_sum } from '../../../../packages/parigp-ts/src/trans1.js';
import { atanhuu } from '../../../../packages/parigp-ts/src/trans2.js';
import { mplog2 as sharedLog2 } from '../../../../packages/parigp-ts/src/qfb.js';
import { trunc2nr as nativeTruncateShift } from '../../../../packages/parigp-ts/src/kernel/none/mp_indep.js';
import { truncr as qfbTruncate } from '../../../../packages/parigp-ts/src/qfb.js';
import * as buchReal from '../../../../packages/parigp-ts/src/buch.js';
import { cmprr as buchRealCompare } from '../../../../packages/parigp-ts/src/buch.js';
import { dbltor as nativeDoubleToReal, rtodbl as nativeRealToDouble } from '../../../../packages/parigp-ts/src/kernel/none/mp_indep.js';
import { dbltor as buchDoubleToReal, rtodbl as buchRealToDouble } from '../../../../packages/parigp-ts/src/buch.js';
import { gcvtoi as qfbIntegerError } from '../../../../packages/parigp-ts/src/qfb.js';
import { gcvtoi as buchIntegerError } from '../../../../packages/parigp-ts/src/buch.js';
import { sqrti as nativeIntegerSqrt } from '../../../../packages/parigp-ts/src/qfb.js';
import { sqrti as buchIntegerSqrt } from '../../../../packages/parigp-ts/src/buch.js';
import { sqrtremi as nativeIntegerSqrtRem } from '../../../../packages/parigp-ts/src/kernel/gmp/mp.js';
import { sqrtr as nativeSqrt, sqrtr_abs as nativeAbsSqrt, type MpComplex } from '../../../../packages/parigp-ts/src/qfb.js';
import { mulri as nativeRealIntegerMul, mulsr as nativeSignedRealMul, mulrs as nativeRealSignedMul } from '../../../../packages/parigp-ts/src/qfb.js';
import { rtor as nativeRealConvert, real_0 as nativeZero, real_0_bit as nativeZeroBit, real_1 as nativeOne, mulir as nativeIntegerMul, mulrr as nativeRealMul, sqrr as nativeRealSquare } from '../../../../packages/parigp-ts/src/qfb.js';
import { PariInvError } from '../../../../packages/parigp-ts/src/matkermod.js';
import { PariError } from '../../../../packages/parigp-ts/src/errors.js';
import { divrr as realDivReal, divir as integerDivReal, divri as realDivInteger, divru as realDivWord } from '../../../../packages/parigp-ts/src/qfb.js';
import { gen_product as nativeProductSchedule, gen_powu_i as nativePowerSchedule } from '../../../../packages/parigp-ts/src/bb_group.js';
import { ZV_prod as nativeVectorProduct } from '../../../../packages/parigp-ts/src/ZV.js';
import { invr as nativeReciprocal } from '../../../../packages/parigp-ts/src/kernel/none/mp_indep.js';
import { mppi as nativePi, mpexp as nativeExp, exp1r_abs as nativeExpm1Abs, agm1r_abs as nativeAgm, powru as nativeRealPower } from '../../../../packages/parigp-ts/src/trans1.js';
import { mpfactr as nativeRealFactorial } from '../../../../packages/parigp-ts/src/trans2.js';
import { bernfrac as nativeBernoulli } from '../../../../packages/parigp-ts/src/bern.js';
import { mpfact as nativeIntegerFactorial, mulu_interval_step as nativeIntervalProduct } from '../../../../packages/parigp-ts/src/arith1.js';
import { rdivii as nativeRealQuotient } from '../../../../packages/parigp-ts/src/kernel/none/level1.js';
import { quadratic_prec_mask as nativePrecisionMask } from '../../../../packages/parigp-ts/src/Zp.js';
import { itor as nativeIntegerReal, shiftr as nativeShift, logr_abs as nativeLog } from '../../../../packages/parigp-ts/src/qfb.js';
import { createHash } from 'node:crypto';
import { gmp_primesieve, mpz_fac_ui, mpz_oddfac_1 } from '../../../../packages/sagemath-ts/src/types/gmp_factorial.js';
import { mpz_remove } from '../../../../packages/sagemath-ts/src/types/gmp.js';
import { sorted as pythonSorted } from '../../../../packages/sagemath-ts/src/types/python_sort.js';
import { two_squares_pyx, three_squares_pyx, four_squares_pyx, is_sum_of_two_squares_pyx } from '../../../../packages/sagemath-ts/src/rings/sum_of_squares.js';
import {
  eulerphi as nativeEulerphi,
  numdiv as nativeNumdiv,
} from '../../../../packages/parigp-ts/src/arith2.js';
import { Z_factor as nativeFactor } from '../../../../packages/parigp-ts/src/ifactor.js';
import { addrr, subrr, addir, subir, addrs, subrs, type MpReal, Zn_quad_roots } from '../../../../packages/parigp-ts/src/qfb.js';
import {
  ProductTree,
  prod_with_derivative,
} from '../../../../packages/sagemath-ts/src/rings/generic.js';
import { PrimeField } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
import { FiniteFieldPrime } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_prime.js';
import { Mod } from '../../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod.js';
import { AttributeError, ValueError, OverflowError } from '../../../../packages/sagemath-ts/src/errors.js';
/**
 * sagemath-ts side of the `arith_special` property-test area.
 *
 * Cases: tests/property/cases/arith_special.cases.json
 * SageMath counterpart: tests/property/python/areas/arith_special.py
 */

import {
  two_squares, three_squares, four_squares, sum_of_k_squares,
  mqrr_rational_reconstruction,
  sort_complex_numbers_for_display,
  valuation,
  algebraic_dependency,
  get_gcd,
  xlcm,
  euler_phi,
  radical,
  number_of_divisors,
  prime_factors,
  smooth_part,
  coprime_part,
  crt,
  CRT,
  CRT_list,
  CRT_basis,
  CRT_vectors,
  get_inverse_mod,
  quadratic_residues,
  fundamental_discriminant,
  odd_part,
  prime_to_m_part,
  continuant,
  squarefree_divisors,
  hilbert_symbol,
  hilbert_conductor,
  hilbert_conductor_inverse,
  next_prime_power,
  previous_prime_power,
  next_probable_prime,
  is_power_of_two,
  is_pseudoprime_power,
  is_prime_power,
  next_prime,
  previous_prime,
  prime_powers,
  eratosthenes,
  bernoulli,
  binomial,
  carmichael_lambda,
  dedekind_psi,
  dedekind_sum,
  differences,
  integer_floor,
  integer_ceil,
  integer_trunc,
  primes,
  primes_first_n,
  factorial,
  fibonacci,
  lucas_number,
  multinomial,
  nth_prime,
  primitive_root,
  subfactorial,
} from '../../../../packages/sagemath-ts/src/arith/misc.js';

import { Z_pvalrem } from '../../../../packages/parigp-ts/src/gen2.js';
import { hilbert as pariHilbert } from '../../../../packages/parigp-ts/src/arith1.js';
import { sumdedekind } from '../../../../packages/parigp-ts/src/elltrans.js';
import { fmpq_dedekind_sum } from '../../../../packages/flint-ts/src/fmpq/dedekind_sum.js';
import {
  nextprime as pariNextprime,
  precprime as pariPrecprime,
  isprimepower as pariPrimePower,
} from '../../../../packages/parigp-ts/src/ifactor.js';
import { algdep as pariAlgdep } from '../../../../packages/parigp-ts/src/bibli1.js';
import { prime as pariPrime } from '../../../../packages/parigp-ts/src/prime.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';
import { RealField } from '../../../../packages/sagemath-ts/src/rings/real_mpfr.js';
import { GF, Mod } from '../../../../packages/sagemath-ts/src/index.js';
import { Integer } from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';

const roundingTexts = [
  '',
  'foo',
  '-2.5',
  '+2.5',
  '1_000.5',
  ' nan ',
  '-Infinity',
  '0e1000',
  '1e1000',
  '-1e-1000',
  '\uff11\uff12.\uff15',
  '\u0661.\u0665',
];
const algdepTexts = [
  '0',
  '-0',
  '1',
  '-1',
  '1.5',
  '-1.5',
  '1.4142135623730950488016887242096980785696718753769',
  '3.1415926535897932384626433832795028841971693993751',
  '1e-30',
  '1e30',
  'NaN',
  'inf',
  '-inf',
  '1.6180339887498948482045868343656381177203091798058',
];
const crtBases = [
  [],
  [3n],
  [3n, 5n],
  [60n, 90n, 150n],
  [7n, 6n, 10n],
  [2n, 4n, 8n, 16n],
  [-3n, 5n],
  [-3n, -5n],
  [0n],
  [0n, 3n],
  [1n, -1n, 3n],
  [2n, 3n, 5n, 7n, 11n],
  [4n, 6n, 9n, 10n],
];
const productBases = [
  [],
  [2, 3],
  [3, 2],
  [6],
  [6, 14],
  [-6, 14],
  [6, -14],
  [2, 2],
  [0],
  [1],
  [-1],
  [-2, 3],
  [5, 7, 11, 13],
  [4, 9, 25],
  [2, 3, 5, 7, 11],
  [2, 3, 5, 7, 11, 13, 17],
  [2, 3, 5, 7, 11, 13, 17, 19],
  [2, 3, 5, 7, 11, 13, 17, 19, 23],
  [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47],
  [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53],
  [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59],
  [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127,
  ],
  [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131,
  ],
  [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137,
  ],
  [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
  ],
  [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311,
  ],
  [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313,
  ],
  [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
    661, 673, 677, 683, 691, 701, 709,
  ],
  [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
    661, 673, 677, 683, 691, 701, 709, 719,
  ],
  [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
    661, 673, 677, 683, 691, 701, 709, 719, 727,
  ],
].map((xs) => xs.map(BigInt));
const productFrame = (value: unknown) =>
  JSON.stringify(value, (_, v) => (typeof v === 'bigint' || typeof v === 'number' ? String(v) : v));
const nativeFactorArithmetic = (n: bigint, op: bigint) =>
  op === 0n ? nativeEulerphi(n) : op === 1n ? nativeNumdiv(n) : productFrame(nativeFactor(n));
export const functions = {
  pari_real_double,
  pari_real_integer_error,
  pari_real_sqrt,
  pari_real_mul,
  pari_real_zero,
  pari_real_div,
  pari_product_trace,
  pari_real_trans,
  pari_real_add,
  gmp_factorial: gmpFactorial,
  valuation_dispatch: valuationDispatch,
  valuation_power: (op: bigint, p: bigint, e: bigint, unit: bigint, d: bigint) => op === 10n ? JSON.stringify(mpz_remove(unit * p ** e, p).map(String)) : valuationDispatch(op, unit * p ** e, d, p, 3n),
  valuation_protocol: valuationProtocol,
  gmp_remove: (n: bigint, p: bigint) => JSON.stringify(mpz_remove(n,p).map(String)),
  integer_bit_size: (bits: bigint, delta: bigint, sign: bigint, op: bigint) => new Integer(sign * ((1n << bits) + delta))[op === 0n ? 'nbits' : 'bit_length'](),
  complex_display_order: complexDisplayOrder,
  python_sort_dispatch: pythonSortDispatch,
  mqrr_dispatch: (u: bigint, m: bigint, t: bigint, wrapped: bigint) => {
    const cv = (x: bigint) => wrapped ? new Integer(x) : x;
    return mqrr_rational_reconstruction(cv(u), cv(m), cv(t));
  },
  square_decomposition: squareDecomposition,
  square_decomposition_power: (op: bigint, n: bigint, e: bigint, mode: bigint) => squareDecomposition(op, n * 4n ** e, 0n, mode),
  native_factor_arithmetic: nativeFactorArithmetic,
  native_factor_power: (p: bigint, e: bigint, unit: bigint, op: bigint) =>
    nativeFactorArithmetic(unit * p ** e, op),
  native_quadratic_roots: (p: bigint, e: bigint, v: bigint, shift: bigint, sign: bigint) => {
    const result = Zn_quad_roots(nativeFactor(sign * p ** e), 2n * shift, shift ** 2n - p ** v);
    return productFrame(result === null ? null : [result.Np, result.roots]);
  },
  xlcm_dispatch: (m: bigint, n: bigint, wrapped: bigint) =>
    xlcm((wrapped & 1n ? new Integer(m) : m) as any, (wrapped & 2n ? new Integer(n) : n) as any),
  factor_scalar_dispatch: (n: bigint, operation: bigint, wrapped: bigint) =>
    productFrame(
      (
        [dedekind_psi, carmichael_lambda, euler_phi, radical, number_of_divisors, prime_factors][
          Number(operation)
        ]! as any
      )(wrapped ? new Integer(n) : n)
    ),
  product_tree_dispatch: (
    code: bigint,
    x: bigint,
    operation: bigint,
    delta: bigint,
    wrapped: bigint
  ) => {
    const leaves = productBases[Number(code)]!.map((v) => (wrapped & 1n ? new Integer(v) : v));
    const tree = new ProductTree(wrapped & 2n ? leaves.values() : leaves);
    if (operation === 0n)
      return productFrame([
        tree.__len__(),
        [...tree.__iter__()],
        tree.leaves(),
        tree.layers,
        tree.leaves() === tree.layers[0],
      ]);
    if (operation === 1n) return tree.root();
    if (operation === 2n) return tree.remainders(wrapped & 1n ? new Integer(x) : x);
    const xs = Array.from({ length: Math.max(0, leaves.length + Number(delta)) }, (_, i) =>
      wrapped & 1n ? new Integer(x + BigInt(i * i)) : x + BigInt(i * i)
    );
    if (operation === 3n) {
      const first = tree.interpolation(xs),
        cache = tree._crt_bases;
      return productFrame([first, tree.interpolation(xs), cache === tree._crt_bases]);
    }
    let error = null,
      answer;
    try {
      tree.interpolation(xs);
    } catch (e) {
      error = [(e as Error).name, (e as Error).message];
    }
    try {
      answer = [String(tree.interpolation(leaves.map(() => x)))];
    } catch (e) {
      answer = [(e as Error).name, (e as Error).message];
    }
    return productFrame([error, tree._crt_bases, answer]);
  },
  smooth_part_dispatch: (x: bigint, code: bigint, operation: bigint, mode: bigint) => {
    const leaves = productBases[Number(code)]!.map((v) => (mode & 2n ? new Integer(v) : v)),
      base = mode & 1n ? new ProductTree(leaves) : mode & 4n ? leaves.values() : leaves;
    const result = ([smooth_part, coprime_part][Number(operation)]! as any)(
      mode & 2n ? new Integer(x) : x,
      base
    );
    return operation === 0n ? productFrame(result) : result;
  },
  smooth_nontermination: (x: bigint, code: bigint, operation: bigint) =>
    ([smooth_part, coprime_part][Number(operation)]! as any)(x, productBases[Number(code)]!),
  product_derivative_dispatch: (length: bigint, x: bigint, wrapped: bigint) => {
    const pairs = Array.from({ length: Number(length) }, (_, i) => {
      const a = x + BigInt(i),
        b = (-1n) ** BigInt(i) * BigInt(i);
      return [wrapped & 1n ? new Integer(a) : a, wrapped & 1n ? new Integer(b) : b] as const;
    });
    return prod_with_derivative(wrapped & 2n ? pairs.values() : pairs);
  },
  crt_scalar_dispatch: (a: bigint, b: bigint, m: bigint, n: bigint, wrapped: bigint) =>
    crt(
      ...([a, b, m, n].map((v, i) => (wrapped & (1n << BigInt(i)) ? new Integer(v) : v)) as [
        any,
        any,
        any,
        any,
      ])
    ),
  crt_list_dispatch: (
    code: bigint,
    offset: bigint,
    delta: bigint,
    mode: bigint,
    wrapped: bigint
  ) => {
    const moduli = crtBases[Number(code)]!.map((v) => (wrapped ? new Integer(v) : v));
    const values = Array.from({ length: Math.max(0, moduli.length + Number(delta)) }, (_, i) => {
      const v = offset + BigInt(i * (i + 1));
      return wrapped ? new Integer(v) : v;
    });
    if (mode === 3n || mode === 7n || mode === 8n) {
      const input = values
        .slice(0, moduli.length)
        .map((v, i) =>
          mode === 3n
            ? Mod(v, moduli[i]!)
            : mode === 7n
              ? new PrimeField(moduli[i]!).__call__(v)
              : new FiniteFieldPrime(moduli[i]!).__call__(v)
        );
      const result = (CRT_list as any)(input);
      return [
        result.value,
        result.modulus ?? result.parent.characteristic,
        input.length === 1 && result === input[0],
      ];
    }
    if (mode === 4n) return (CRT_list as any)('bad', moduli);
    if (mode === 5n) return (CRT_list as any)(values, 'bad');
    if (mode === 6n) return (CRT_list as any)(values);
    return ([CRT_list, crt, CRT][Number(mode)]! as any)(values, moduli);
  },
  crt_basis_dispatch: (length: bigint, encoding: bigint, required: bigint, wrapped: bigint) => {
    const moduli = [];
    for (let i = 0n; i < length; i++) {
      const v = (encoding % 7n) - 3n;
      encoding /= 7n;
      moduli.push(wrapped ? new Integer(v) : v);
    }
    if (length < 0n)
      moduli.push(...crtBases[Number(encoding)]!.map((v) => (wrapped ? new Integer(v) : v)));
    return JSON.stringify(CRT_basis(moduli as any, !!required), (_, v) =>
      typeof v === 'bigint' ? String(v) : v
    );
  },
  crt_vector_dispatch: (
    code: bigint,
    offset: bigint,
    columns: bigint,
    shape: bigint,
    wrapped: bigint
  ) => {
    const moduli = crtBases[Number(code)]!.map((v) => (wrapped ? new Integer(v) : v));
    const values = moduli.map((_, i) =>
      Array.from({ length: Number(columns) }, (_, j) => {
        const v = offset + BigInt(i * j + j);
        return wrapped ? new Integer(v) : v;
      })
    );
    if (shape === 1n)
      values.push(
        Array.from({ length: Number(columns) }, () => (wrapped ? new Integer(offset) : offset))
      );
    if (shape === 2n && values.length) values[values.length - 1] = values.at(-1)!.slice(0, -1);
    if (shape === 3n && values.length) values.at(-1)!.push(wrapped ? new Integer(offset) : offset);
    return CRT_vectors(values as any, moduli as any);
  },
  arithmetic_factory: (order: bigint, a: bigint, b: bigint, operation: bigint, wrapped: bigint) => {
    const factory = [get_gcd, get_inverse_mod][Number(operation)]!;
    const f = factory(wrapped & 1n ? new Integer(order) : order);
    return f(wrapped & 2n ? new Integer(a) : a, wrapped & 4n ? new Integer(b) : b);
  },
  arithmetic_factory_selection: (order: bigint, operation: bigint, wrapped: bigint) => {
    const maker = [get_gcd, get_inverse_mod][Number(operation)]!,
      size = wrapped ? new Integer(order) : order;
    const f = maker(size),
      g = maker(size);
    return f.name.replace(/^bound /, '') + '|' + Number(f === g);
  },
  arith_scalar_edge: (n: bigint, m: bigint, operation: bigint, wrapped: bigint) => {
    const x = wrapped & 1n ? new Integer(n) : n,
      y = wrapped & 2n ? new Integer(m) : m;
    if (operation === 0n) return JSON.stringify(quadratic_residues(x as any).map(String));
    const result =
      operation === 1n
        ? fundamental_discriminant(x as any)
        : operation === 2n
          ? odd_part(x as any)
          : prime_to_m_part(x as any, y as any);
    return JSON.stringify([
      String(result),
      typeof result === 'bigint' ? 'integer' : 'wrapped_integer',
    ]);
  },
  arith_continuant: (
    length: bigint,
    power: bigint,
    offset: bigint,
    order: bigint,
    mode: bigint,
    wrapped: bigint
  ) => {
    const values = Array.from({ length: Number(length) }, (_, i) => {
      const x = BigInt(i) ** power + offset;
      return wrapped & 1n ? new Integer(x) : x;
    });
    const result =
      mode === 0n
        ? continuant(values as any)
        : continuant(values as any, (wrapped & 2n ? new Integer(order) : order) as any);
    return JSON.stringify([
      String(result),
      typeof result === 'bigint' ? 'integer' : 'wrapped_integer',
    ]);
  },
  arith_squarefree_prefix: (n: bigint, take: bigint, wrapped: bigint) => {
    const iterator = squarefree_divisors((wrapped ? new Integer(n) : n) as any),
      values: string[] = [];
    let done = false;
    for (let i = 0n; i < take; i++) {
      const value = iterator.next();
      if (value.done) {
        done = true;
        break;
      }
      values.push(String(value.value));
    }
    return JSON.stringify([values, done]);
  },
  native_valuation_unit: (n: bigint, p: bigint) => {
    const [v, u] = Z_pvalrem(n, p);
    return [BigInt(v), u];
  },
  native_hilbert: (a: bigint, b: bigint, p: bigint) => pariHilbert(a, b, p),
  hilbert_dispatch: (
    an: bigint,
    ad: bigint,
    bn: bigint,
    bd: bigint,
    p: bigint,
    mode: bigint,
    wrapped: bigint
  ) => {
    const a = ad !== 1n ? new Rational(an, ad) : wrapped & 1n ? new Integer(an) : an;
    const b = bd !== 1n ? new Rational(bn, bd) : wrapped & 2n ? new Integer(bn) : bn;
    const prime = wrapped & 4n ? new Integer(p) : p;
    return mode === 0n
      ? hilbert_symbol(a as any, b as any, prime as any)
      : hilbert_symbol(
          a as any,
          b as any,
          prime as any,
          ['pari', 'direct', 'all', 'invalid', null, ''][Number(mode) - 1] as any
        );
  },
  hilbert_conductor_dispatch: (a: bigint, b: bigint, operation: bigint, wrapped: bigint) => {
    const x = wrapped & 1n ? new Integer(a) : a,
      y = wrapped & 2n ? new Integer(b) : b;
    return operation ? hilbert_conductor_inverse(x as any) : hilbert_conductor(x as any, y as any);
  },
  native_prime_traversal: (value: bigint, operation: bigint) => {
    if (operation === 0n) return pariNextprime(value);
    if (operation === 1n) return pariPrecprime(value);
    const data = pariPrimePower(value);
    return data === null ? null : [data[0], BigInt(data[1])];
  },
  prime_traversal: (value: bigint, operation: bigint, wrapped: bigint, data: bigint) => {
    const n = wrapped ? new Integer(value) : value;
    const f = [
      next_prime_power,
      previous_prime_power,
      next_probable_prime,
      is_power_of_two,
      is_pseudoprime_power,
      is_prime_power,
      next_prime,
      previous_prime,
    ][Number(operation)]!;
    return operation === 4n || operation === 5n
      ? data
        ? (f as typeof is_prime_power)(n, true)
        : (f as typeof is_prime_power)(n, false)
      : (f as (n: any) => unknown)(n);
  },
  prime_power_range: (start: bigint, stop: bigint, mode: bigint, wrapped: bigint) => {
    const a = wrapped & 1n ? new Integer(start) : start,
      b = wrapped & 2n ? new Integer(stop) : stop;
    return mode === 0n
      ? prime_powers(a as any)
      : mode === 1n
        ? prime_powers(a as any, b as any)
        : eratosthenes(a as any);
  },
  native_algdep: (bits: bigint, degree: bigint) => {
    const view = new DataView(new ArrayBuffer(8));
    view.setBigUint64(0, bits);
    return pariAlgdep(view.getFloat64(0), degree);
  },
  algdep_real: (
    kind: bigint,
    value: bigint,
    denominator: bigint,
    precision: bigint,
    degree: bigint,
    mode: bigint,
    bound: bigint,
    first: bigint,
    second: bigint
  ) => {
    const z =
      kind === 0n
        ? value
        : kind === 1n
          ? new Integer(value)
          : kind === 2n
            ? new Rational(value, denominator)
            : new RealField(Number(precision)).__call__(algdepTexts[Number(value)]!);
    const options =
      mode === 0n
        ? {}
        : mode === 1n
          ? { height_bound: bound }
          : mode === 2n
            ? { height_bound: bound, proof: true }
            : mode === 3n
              ? { proof: true }
              : mode === 4n
                ? { known_bits: first }
                : mode === 5n
                  ? { use_bits: first }
                  : mode === 6n
                    ? { known_digits: first }
                    : mode === 7n
                      ? { use_digits: first }
                      : {
                          known_bits: first,
                          use_bits: second,
                          known_digits: first,
                          use_digits: second,
                        };
    if (mode === 9n || mode === 10n) {
      const real = z as ReturnType<RealField['__call__']>;
      return mode === 9n ? real.algebraic_dependency(Number(degree)) : real.algdep(Number(degree));
    }
    return algebraic_dependency(z, degree, options);
  },
  algdep_dispatch: (
    bits: bigint,
    degree: bigint,
    mode: bigint,
    bound: bigint,
    first: bigint,
    second: bigint
  ) => {
    const view = new DataView(new ArrayBuffer(8));
    view.setBigUint64(0, bits);
    const z = view.getFloat64(0);
    const options =
      mode === 0n
        ? {}
        : mode === 1n
          ? { height_bound: bound }
          : mode === 2n
            ? { height_bound: bound, proof: true }
            : mode === 3n
              ? { proof: true }
              : mode === 4n
                ? { known_bits: first }
                : mode === 5n
                  ? { use_bits: first }
                  : mode === 6n
                    ? { known_digits: first }
                    : mode === 7n
                      ? { use_digits: first }
                      : {
                          known_bits: first,
                          use_bits: second,
                          known_digits: first,
                          use_digits: second,
                        };
    return algebraic_dependency(z, degree, options);
  },
  integer_rounding: (kind: bigint, a: bigint, b: bigint, precision: bigint, operation: bigint) => {
    let value: any;
    if (kind === 0n) value = a;
    else if (kind === 1n) value = new Integer(a);
    else if (kind === 2n) value = new Rational(a, b);
    else if (kind === 3n) value = Number(a) / Number(b);
    else if (kind === 4n) value = new RealField(Number(precision)).__call__(`${a}e${b}`);
    else if (kind === 5n)
      value = new RealField(Number(precision)).__call__(['NaN', 'inf', '-inf', '-0'][Number(a)]!);
    else if (kind === 6n) {
      const data = new DataView(new ArrayBuffer(8));
      data.setBigUint64(0, a);
      value = data.getFloat64(0);
    } else if (kind === 7n) value = roundingTexts[Number(a)];
    else if (kind === 8n)
      value = [
        true,
        false,
        null,
        new TextEncoder().encode('2.5'),
        new TextEncoder().encode('x'),
        new TextEncoder().encode('nan'),
      ][Number(a)];
    else if (kind === 12n || kind === 13n) {
      value = new Integer(a);
      const method = () => {
        throw kind === 12n
          ? new AttributeError('rounding sentinel')
          : new ValueError('rounding sentinel');
      };
      Object.defineProperty(value, 'floor', { value: method });
      Object.defineProperty(value, 'ceil', { value: method });
    } else if (kind === 9n) value = Mod(a, 7n);
    else value = GF(kind === 10n ? 7n : 2n).__call__(a);
    return [integer_floor, integer_ceil, integer_trunc][Number(operation)]!(value);
  },
  dedekind_native: (h: bigint, k: bigint, backend: bigint) =>
    backend === 0n ? fmpq_dedekind_sum(h, k) : sumdedekind(h, k),
  dedekind_sum: (h: bigint, k: bigint, mode: bigint, wrapped: bigint) => {
    const algorithm = [undefined, 'default', 'flint', 'pari', 'invalid', null][Number(mode)];
    const r = dedekind_sum(
      wrapped ? new Integer(h) : h,
      wrapped ? new Integer(k) : k,
      algorithm as any
    );
    return [r.numerator, r.denominator];
  },
  prime_native: (n: bigint) => pariPrime(n),
  first_primes_nonfinite: (mode: bigint) =>
    primes_first_n(mode === 0n ? NaN : mode === 1n ? Infinity : -Infinity),
  primes_iterator: (
    mode: bigint,
    start: bigint,
    stop: bigint,
    count: bigint,
    proof: bigint,
    wrapped: bigint
  ) => {
    const a = wrapped ? new Integer(start) : start;
    const b = wrapped ? new Integer(stop) : stop;
    const p = proof === 0n ? undefined : proof === 1n;
    const iterator =
      mode === 0n ? primes() : mode === 1n ? primes(a) : primes(a, mode === 3n ? Infinity : b, p);
    const values: bigint[] = [];
    let done = false;
    for (let i = 0n; i < count; ++i) {
      const item = iterator.next();
      if (item.done) {
        done = true;
        break;
      }
      values.push(item.value);
    }
    return `${values.join(',')}|${done ? 1 : 0}`;
  },
  differences: (length: bigint, degree: bigint, offset: bigint, order: bigint, wrapped: bigint) => {
    const values = Array.from({ length: Number(length) }, (_, i) => (BigInt(i) + offset) ** degree);
    return differences(
      wrapped ? values.map((x) => new Integer(x)) : values,
      wrapped ? new Integer(order) : order
    );
  },
  binomial: (n: bigint, k: bigint) => binomial(n, k),
  fibonacci: (n: bigint) => fibonacci(n),
  lucas_number: (n: bigint) => lucas_number(n),
  factorial: (n: bigint) => factorial(n),
  bernoulli_numerator: (n: bigint) => {
    const b = bernoulli(n);
    return b.numerator;
  },
  bernoulli_denominator: (n: bigint) => {
    const b = bernoulli(n);
    return b.denominator;
  },
  multinomial: (...args: bigint[]) => multinomial(...args),
  primitive_root: (n: bigint) => primitive_root(n),
  nth_prime: (n: bigint) => nth_prime(n),
  nth_prime_wrapped: (n: bigint, wrapped: bigint) => nth_prime(wrapped ? new Integer(n) : n),
  primes_first_n: (n: bigint, mode: bigint) =>
    primes_first_n(
      mode === 0n ? Number(n) : mode === 1n ? n : mode === 2n ? new Integer(n) : Number(n) / 1000
    ),
  subfactorial: (n: bigint) => subfactorial(n),
  carmichael_lambda: (n: bigint) => carmichael_lambda(n),
  dedekind_psi: (n: bigint) => dedekind_psi(n),
};
function squareDecomposition(op: bigint, n: bigint, k: bigint, mode: bigint): unknown {
  const input = mode % 2n ? new Integer(n) : n;
  if (op === 3n) {
    const count = mode / 2n === 0n ? Number(k)
      : mode / 2n === 1n ? k
      : mode / 2n === 2n ? new Integer(k)
      : mode / 2n === 3n ? Number(k) + (k >= 0n ? 0.5 : -0.5)
      : [NaN, Infinity, -Infinity][Number(k)]!;
    return sum_of_k_squares(count, input);
  }
  return [two_squares, three_squares, four_squares, undefined,
    two_squares_pyx, three_squares_pyx, four_squares_pyx,
    is_sum_of_two_squares_pyx][Number(op)]!(input);
}

function complexDisplayOrder(base: bigint, step: bigint, imag: bigint, mode: bigint, order: bigint, size: bigint){
 const ims=[[0,-0,1,-1,1e-12,-1e-12,2],[1,-1,2,-2,3,-3,4],[1e-12,-1e-12,0,1,-1,2,-2],[1e-10,-1e-10,1e-11,-1e-11,1,0,-0],[1e-300,-1e-300,0,1e300,-1e300,1,-1]][Number(imag)]!;
 let vals:any[]=[];const view=new DataView(new ArrayBuffer(8));
 for(let i=0;i<Number(size);i++){view.setBigUint64(0,BigInt.asUintN(64,base+(BigInt(i%7)-3n)*step),false);const z={re:view.getFloat64(0,false),im:ims[i%7]!};vals.push(mode?[z,i]:z);}
 if(order===1n)vals.reverse();else if(order===2n)vals=[...vals.filter((_,i)=>i%2===0),...vals.filter((_,i)=>i%2===1)];
 const before=vals.slice(),out=sort_complex_numbers_for_display(vals);
 return JSON.stringify([out.map(x=>String(before.indexOf(x))),out===vals,vals.every((x,i)=>x===before[i])]);
}

function pythonSortDispatch(seed: bigint, len: bigint, pat: bigint){const n=Number(len),pattern=Number(pat);let x=seed;const values:number[]=[];
 for(let i=0;i<n;i++){
  x=(1664525n*x+1013904223n)&0xffffffffn;const r=Number(x);
  let v:number;
  if(pattern===0)v=r%101-50;
  else if(pattern===1)v=i+r%3;
  else if(pattern===2)v=n-i-r%3;
  else if(pattern===3)v=r%11===0?NaN:r%101-50;
  else if(pattern===4)v=i%4===0?NaN:r%101-50;
  else if(pattern===5)v=Math.floor(i/64)*128+63-i%64;
  else v=(i%4)*n+Math.floor(i/4);
  values.push(v);
 }
 let trace=14695981039346656037n;
 const out=pythonSorted(values.map((_,i)=>i),(a,b)=>{
  trace=BigInt.asUintN(64,((trace^BigInt(a))*1099511628211n)^BigInt(b));
  return values[a]!<values[b]!;
 });
 return JSON.stringify([out.map(String),trace.toString()]);
}

function valuationDispatch(op: bigint, n: bigint, d: bigint, p: bigint, mode: bigint){const base=mode>=2n?new Integer(p):p;
 let r:any;
 if(op<3n){const z=mode%2n?new Integer(n):n;r=op===0n?valuation(z,base):new Integer(n)[op===1n?'valuation':'ord'](base);}
 else{const z=new Rational(n,d);r=op===3n?valuation(z as any,base):z[['valuation','ord','val_unit','numerator_valuation','denominator_valuation','padic_valuation'][Number(op)-4] as any](base);}
 return JSON.stringify({value:Array.isArray(r)?r.map(String):String(r)});
}

function valuationProtocol(n: bigint, p: bigint, mode: bigint): string {
  const calls: string[] = [], saved = [n,p];
  const obj = {
    _integer_() { calls.push('integer'); return n; },
    valuation(_p: bigint): bigint | bigint[] {
      calls.push('valuation');
      if (mode === 1n) throw new AttributeError('inside');
      if (mode === 2n) throw new ValueError('valuation failed');
      return mode === 4n ? saved : n + _p;
    }
  };
  if (mode === 3n) Object.defineProperty(obj, 'valuation', {get(){calls.push('get');throw new AttributeError('missing');}});
  try {
    const value = mode >= 5n ? valuation({_integer_(){
      calls.push('integer');
      if (mode === 6n) throw new AttributeError('integer failed');
      if (mode === 7n) throw new ValueError('integer failed');
      return n;
    }}, p) : valuation(obj, p);
    return JSON.stringify({value: Array.isArray(value) ? value.map(String) : String(value), same: value === saved, calls});
  } catch(e) {
    return JSON.stringify({error: e instanceof Error ? e.name : String(e), message: e instanceof Error ? e.message : '', calls});
  }
}

function gmpFactorial(op: bigint, n: bigint): string {
  if (op === 3n) {
    const mask = gmp_primesieve(n);
    let compositeCount = 0;
    const bytes = Buffer.alloc(mask.byteLength);
    for (let i = 0; i < mask.length; i++) {
      let word = mask[i]!;
      bytes.writeUInt32LE(word, i * 4);
      while (word) { word &= word - 1; compositeCount++; }
    }
    return JSON.stringify([String(mask.length * 32 - compositeCount), createHash('sha256').update(bytes).digest('hex')]);
  }
  const value = op >= 5n ? factorial(op >= 10n ? new Integer(n) : n, [undefined, 'gmp', 'unknown', 'pari', null][Number(op >= 10n ? op-10n : op-5n)] as any) : op === 0n ? new Integer(n).factorial().value : op === 1n ? mpz_fac_ui(n) : mpz_oddfac_1(n, op === 4n ? 1 : 0);
  let hex = value.toString(16);
  if (hex.length % 2) hex = '0' + hex;
  return JSON.stringify([String(value.toString(2).length), createHash('sha256').update(Buffer.from(hex,'hex')).digest('hex')]);
}

/** Exact native real storage, including precision and zero accuracy. */
function pari_real_add(op: bigint, px0: bigint, py0: bigint, ex0: bigint, ey0: bigint,
  mx0: bigint, my0: bigint, sx0: bigint, sy0: bigint): string {
  const [px, py, ex, ey, mx, my, sx, sy] = [px0, py0, ex0, ey0, mx0, my0, sx0, sy0].map(Number);
  const real = (p: number, e: number, mode: number, sign: number): MpReal => {
    const low = 1n << BigInt(p - 1), mask = low - 1n;
    const m = low + [0n, mask, 1n, mask - 1n, mask / 3n, mask / 3n * 2n,
      mask - ((1n << BigInt(Math.min(64, p - 1))) - 1n), mask / 7n, low / 2n - 1n][mode]!;
    return { s: sign as -1 | 0 | 1, e, p: sign ? p : 0, m: sign ? m : 0n };
  };
  const x = real(px!, ex!, mx!, sx!), y = real(py!, ey!, my!, sy!);
  const shift = ex! - px! + 1;
  const n = BigInt(sx!) * (shift >= 0 ? x.m << BigInt(shift) : x.m >> BigInt(-shift));
  const z = op === 0n ? addrr(x, y) : op === 1n ? subrr(x, y)
    : op === 2n ? addir(n, y) : op === 3n ? subir(n, y)
    : op === 4n ? addrs(y, Number(n)) : subrs(y, Number(n));
  return JSON.stringify([z.s, String(z.e), String(z.m), z.p]);
}

/** Match complete original real storage or hash the entire exact integer value. */
function pari_real_trans(op: bigint, p0: bigint, n: bigint, shift0: bigint, k: bigint): string {
  const p = Number(p0), shift = Number(shift0);
  const frame = (x: MpReal | MpReal<bigint>) => [x.s, String(x.e), String(x.m), x.p];
  const hash = (value: bigint) => {
    const mag = value < 0n ? -value : value;
    const bits = mag ? mag.toString(2).length : 0;
    let hex = mag ? mag.toString(16) : '';
    if (hex.length % 2) hex = '0' + hex;
    return [String(value < 0n ? -1 : value > 0n ? 1 : 0), String(bits), createHash('sha256').update(Buffer.from(hex, 'hex')).digest('hex')];
  };
  let result: unknown;
  if (op === 0n) result = frame(nativePi(p));
  else if (op === 3n) {
    // cypari2 performs this signed-C-long conversion before invoking the kernel.
    if (n < -(1n << 63n) || n >= 1n << 63n) throw new OverflowError('Python int too large to convert to C long');
    result = frame(nativeRealFactorial(n, p));
  } else if (op === 7n || op === 8n) result = frame(factorial(op === 8n ? new Integer(n) : n, 'pari'));
  else if (op === 9n) result = hash(nativeIntegerFactorial(n));
  else if (op === 10n) result = nativeBernoulli(Number(n)).map(String);
  else if (op === 11n) result = frame(nativeRealQuotient(shift >= 0 ? n << BigInt(shift) : n >> BigInt(-shift), k, p));
  else if (op === 12n) result = String(nativePrecisionMask(Number(n)));
  else if (op === 13n) result = hash(nativeIntervalProduct(n, shift0, k));
  else {
    const x = nativeShift(nativeIntegerReal(n, p), shift);
    result = frame(op === 1n ? nativeExp(x) : op === 2n ? nativeExpm1Abs(x)
      : op === 4n ? nativeLog(x) : op === 5n ? nativeAgm(x) : op === 14n ? nativeReciprocal(x) : nativeRealPower(x, k));
  }
  return JSON.stringify(result);
}

function pari_product_trace(op: bigint, size: bigint, seed: bigint): string {
  const trace: string[] = [], mask = (1n << 64n) - 1n;
  const multiply = (a: bigint, b: bigint) => {
    trace.push(`M:${a},${b}`);
    return op === 2n ? a * b % ((1n << 61n) - 1n) : ((a * 65599n) ^ b) & mask;
  };
  let result: bigint;
  if (op === 2n) result = nativePowerSchedule(seed, size, (a) => {
    trace.push(`S:${a}`);
    return a * a % ((1n << 61n) - 1n);
  }, multiply);
  else {
    let state = seed;
    const values: bigint[] = [];
    for (let i = 0n; i < size; i++) {
      state = (1664525n * state + 1013904223n) & ((1n << 32n) - 1n);
      values.push(state);
    }
    result = op === 0n ? nativeVectorProduct(values) : nativeProductSchedule(values, multiply);
  }
  const bits = result ? result.toString(2).length : 0;
  let hex = result ? result.toString(16) : '';
  if (hex.length % 2) hex = '0' + hex;
  return JSON.stringify([String(bits), createHash('sha256').update(Buffer.from(hex, 'hex')).digest('hex'), createHash('sha256').update(trace.join('|')).digest('hex')]);
}

function pari_real_div(op: bigint, px0: bigint, py0: bigint, ex0: bigint, ey0: bigint,
  mx0: bigint, my0: bigint, sx0: bigint, sy0: bigint): string {
  const real = (p: number, e: number, mode: bigint, sign: number): MpReal => {
    const low = 1n << BigInt(p - 1), mask = low - 1n;
    const m = mode < 0n ? -mode : low + [0n, mask, 1n, mask - 1n, mask / 3n, mask / 3n * 2n,
      mask - ((1n << BigInt(Math.min(64, p - 1))) - 1n), mask / 7n][Number(mode)]!;
    return { s: sign as -1 | 0 | 1, e, p, m };
  };
  const integer = (x: MpReal) => {
    const shift = x.e - x.p + 1;
    return BigInt(x.s) * (shift >= 0 ? x.m << BigInt(shift) : x.m >> BigInt(-shift));
  };
  const x = real(Number(px0), Number(ex0), mx0, Number(sx0));
  const y = real(Number(py0), Number(ey0), my0, Number(sy0));
  const z = op === 0n ? realDivReal(x, y) : op === 1n ? integerDivReal(integer(x), y)
    : op === 2n ? realDivInteger(x, integer(y)) : realDivWord(x, integer(y));
  return JSON.stringify([z.s, String(z.e), String(z.m), z.p]);
}


function pari_real_zero(op: bigint, p0: bigint, q0: bigint, e0: bigint, n: bigint): string {
  const p = Number(p0), q = Number(q0), e = Number(e0);
  const x = p ? { ...nativeIntegerReal(0n, p), e } : nativeZeroBit(e);
  try {
    const z = op === 0n ? nativeIntegerReal(0n, q) : op === 1n ? nativeRealConvert(x, q)
      : op === 2n ? nativeIntegerMul(n, x) : op === 3n ? integerDivReal(n, x)
      : op === 4n ? nativeReciprocal(x) : op === 5n ? realDivReal(nativeOne(q), x)
      : op === 6n ? nativeRealMul(x, nativeOne(q)) : op === 7n ? nativeRealSquare(x)
      : op === 8n ? addrr(x, nativeOne(q)) : op === 9n ? nativeRealQuotient(0n, n, q)
      : nativeZero(q);
    return JSON.stringify([z.s, String(z.e), String(z.m), z.p]);
  } catch (error) {
    if (!(error instanceof PariInvError || error instanceof PariError) || !error.message.startsWith('impossible inverse in ')) throw error;
    return JSON.stringify(['inverse', error.message.split(':', 1)[0]]);
  }
}


function pari_real_mul(op: bigint, px0: bigint, py0: bigint, ex0: bigint, ey0: bigint,
  mx0: bigint, my0: bigint, sx0: bigint, sy0: bigint): string {
  const real = (p: number, e: number, mode: bigint, sign: number): MpReal => {
    const low = 1n << BigInt(p - 1), mask = low - 1n;
    const m = mode < 0n ? -mode : low + [0n, mask, 1n, mask - 1n, mask / 3n, mask / 3n * 2n,
      mask - ((1n << BigInt(Math.min(64, p - 1))) - 1n), mask / 7n][Number(mode)]!;
    return { s: sign as -1 | 0 | 1, e, p, m };
  };
  const integer = (x: MpReal) => {
    const shift = x.e - x.p + 1;
    return BigInt(x.s) * (shift >= 0 ? x.m << BigInt(shift) : x.m >> BigInt(-shift));
  };
  const x = real(Number(px0), Number(ex0), mx0, Number(sx0));
  const y = real(Number(py0), Number(ey0), my0, Number(sy0));
  const z = op === 0n ? nativeRealMul(x, y) : op === 1n ? nativeIntegerMul(integer(x), y)
    : op === 2n ? nativeRealIntegerMul(x, integer(y)) : op === 3n ? nativeRealSquare(x)
    : op === 4n ? nativeRealMul(x, x) : op === 5n ? nativeSignedRealMul(Number(integer(x)), y)
    : nativeRealSignedMul(x, Number(integer(y)));
  return JSON.stringify([z.s, String(z.e), String(z.m), z.p]);
}


function pari_real_sqrt(op: bigint, p: bigint, e: bigint, m: bigint, sign: bigint): string {
  if (op === 4n || op === 5n) return JSON.stringify(['integer-rootonly', String((op === 4n ? nativeIntegerSqrt : buchIntegerSqrt)((m << e) + sign))]);
  if (op === 2n) return JSON.stringify(['integer-root', ...nativeIntegerSqrtRem((m << e) + sign).map(String)]);
  if (op === 3n) {
    const digest = createHash('sha256');
    for (let i = 0n; i < m; i++) digest.update(nativeIntegerSqrtRem(e + i * sign).map(String).join(',') + ';');
    return digest.digest('hex');
  }
  const x: MpReal = { s: Number(sign) as -1 | 0 | 1, p: Number(p), e: Number(e), m };
  const frame = (z: bigint | MpReal | MpComplex): unknown[] => typeof z === 'bigint' ? ['integer', String(z)]
    : 're' in z ? ['complex', frame(z.re), frame(z.im)] : ['real', z.s, String(z.e), String(z.m), z.p];
  return JSON.stringify(frame(op ? nativeSqrt(x) : nativeAbsSqrt(x)));
}


function pari_real_integer_error(op: bigint, p0: bigint, e: bigint, mode: bigint, s: bigint): string {
  const p = Number(p0), low = 1n << BigInt(Math.max(0, p - 1)), tail = low - 1n;
  const m = s ? low + [0n, tail, 1n, tail ? tail - 1n : 0n, tail / 3n, tail / 3n * 2n, tail / 7n, tail / 5n][Number(mode)]! : 0n;
  const x: MpReal = { s: Number(s) as -1 | 0 | 1, p, e: Number(e), m };
  if (!op) return JSON.stringify(qfbIntegerError(x).map(String));
  const z = buchIntegerError(x);
  return JSON.stringify([String(z.z), String(z.e)]);
}


function pari_real_double(op: bigint, p: bigint, e: bigint, m: bigint, sign: bigint): string {
  const view = new DataView(new ArrayBuffer(8));
  if (op === 0n || op === 2n) {
    view.setBigUint64(0, m, false);
    const d = view.getFloat64(0, false);
    const z = op === 0n ? nativeDoubleToReal(d) : buchDoubleToReal(d, p ? Number(p) : undefined);
    return JSON.stringify([z.s, String(z.e), String(z.m), z.p]);
  }
  const x: MpReal<number | bigint> = { s: Number(sign) as -1 | 0 | 1, p: Number(p), e: op >= 4n ? e : Number(e), m };
  view.setFloat64(0, (op === 1n || op === 4n ? nativeRealToDouble : buchRealToDouble)(x as MpReal | MpReal<bigint>), false);
  return JSON.stringify(view.getBigUint64(0, false).toString(16).padStart(16, '0'));
}

functions.pari_real_compare = (p: bigint, e: bigint, m: bigint, s: bigint, q: bigint, f: bigint, n: bigint, t: bigint) =>
  buchRealCompare({p:Number(p),e:Number(e),m,s:Number(s)}, {p:Number(q),e:Number(f),m:n,s:Number(t)});

function pari_buch_real_add(op: bigint, px0: bigint, py0: bigint, ex0: bigint, ey0: bigint,
  mx0: bigint, my0: bigint, sx0: bigint, sy0: bigint): string {
  const [px, py, ex, ey, mx, my, sx, sy] = [px0, py0, ex0, ey0, mx0, my0, sx0, sy0].map(Number);
  const real = (p: number, e: number, mode: number, sign: number): MpReal => {
    const low = 1n << BigInt(p - 1), mask = low - 1n;
    const m = low + [0n, mask, 1n, mask - 1n, mask / 3n, mask / 3n * 2n,
      mask - ((1n << BigInt(Math.min(64, p - 1))) - 1n), mask / 7n, low / 2n - 1n][mode]!;
    return { s: sign as -1 | 0 | 1, e, p: sign ? p : 0, m: sign ? m : 0n };
  };
  const x = real(px!, ex!, mx!, sx!), y = real(py!, ey!, my!, sy!);
  const shift = ex! - px! + 1;
  const n = BigInt(sx!) * (shift >= 0 ? x.m << BigInt(shift) : x.m >> BigInt(-shift));
  const z = op === 0n ? buchReal.addrr(x, y) : op === 1n ? buchReal.subrr(x, y)
    : op === 2n ? addir(n, y) : op === 3n ? subir(n, y)
    : op === 4n ? addrs(y, Number(n)) : subrs(y, Number(n));
  return JSON.stringify([z.s, String(z.e), String(z.m), z.p]);
}

functions.pari_buch_real_add = pari_buch_real_add;
function pari_buch_real_div(op: bigint, px0: bigint, py0: bigint, ex0: bigint, ey0: bigint,
  mx0: bigint, my0: bigint, sx0: bigint, sy0: bigint): string {
  const real = (p: number, e: number, mode: bigint, sign: number): MpReal => {
    const low = 1n << BigInt(p - 1), mask = low - 1n;
    const m = mode < 0n ? -mode : low + [0n, mask, 1n, mask - 1n, mask / 3n, mask / 3n * 2n,
      mask - ((1n << BigInt(Math.min(64, p - 1))) - 1n), mask / 7n][Number(mode)]!;
    return { s: sign as -1 | 0 | 1, e, p, m };
  };
  const integer = (x: MpReal) => {
    const shift = x.e - x.p + 1;
    return BigInt(x.s) * (shift >= 0 ? x.m << BigInt(shift) : x.m >> BigInt(-shift));
  };
  const x = real(Number(px0), Number(ex0), mx0, Number(sx0));
  const y = real(Number(py0), Number(ey0), my0, Number(sy0));
  const z = op === 0n ? buchReal.divrr(x, y) : op === 1n ? buchReal.divir(integer(x), y)
    : op === 2n ? buchReal.divri(x, integer(y)) : buchReal.divru(x, integer(y));
  return JSON.stringify([z.s, String(z.e), String(z.m), z.p]);
}

functions.pari_buch_real_div = pari_buch_real_div;
function pari_buch_real_mul(op: bigint, px0: bigint, py0: bigint, ex0: bigint, ey0: bigint,
  mx0: bigint, my0: bigint, sx0: bigint, sy0: bigint): string {
  const real = (p: number, e: number, mode: bigint, sign: number): MpReal => {
    const low = 1n << BigInt(p - 1), mask = low - 1n;
    const m = mode < 0n ? -mode : low + [0n, mask, 1n, mask - 1n, mask / 3n, mask / 3n * 2n,
      mask - ((1n << BigInt(Math.min(64, p - 1))) - 1n), mask / 7n][Number(mode)]!;
    return { s: sign as -1 | 0 | 1, e, p, m };
  };
  const integer = (x: MpReal) => {
    const shift = x.e - x.p + 1;
    return BigInt(x.s) * (shift >= 0 ? x.m << BigInt(shift) : x.m >> BigInt(-shift));
  };
  const x = real(Number(px0), Number(ex0), mx0, Number(sx0));
  const y = real(Number(py0), Number(ey0), my0, Number(sy0));
  const z = op === 0n ? buchReal.mulrr(x, y) : op === 1n ? buchReal.mulir(integer(x), y)
    : op === 2n ? nativeRealIntegerMul(x, integer(y)) : op === 3n ? buchReal.sqrr(x)
    : op === 4n ? buchReal.mulrr(x, x) : op === 5n ? nativeSignedRealMul(Number(integer(x)), y)
    : nativeRealSignedMul(x, Number(integer(y)));
  return JSON.stringify([z.s, String(z.e), String(z.m), z.p]);
}

functions.pari_buch_real_mul = pari_buch_real_mul;
function pari_buch_real_zero(op: bigint, p0: bigint, q0: bigint, e0: bigint, n: bigint): string {
  const p = Number(p0), q = Number(q0), e = Number(e0);
  const x = p ? { ...nativeIntegerReal(0n, p), e } : nativeZeroBit(e);
  try {
    const z = op === 0n ? buchReal.itor(0n, q) : op === 1n ? buchReal.setprec(x, q)
      : op === 2n ? buchReal.mulir(n, x) : op === 3n ? buchReal.divir(n, x)
      : op === 4n ? nativeReciprocal(x) : op === 5n ? buchReal.divrr(nativeOne(q), x)
      : op === 6n ? buchReal.mulrr(x, nativeOne(q)) : op === 7n ? buchReal.sqrr(x)
      : op === 8n ? buchReal.addrr(x, nativeOne(q)) : op === 9n ? nativeRealQuotient(0n, n, q)
      : buchReal.real_0(q);
    return JSON.stringify([z.s, String(z.e), String(z.m), z.p]);
  } catch (error) {
    if (!(error instanceof PariInvError || error instanceof PariError) || !error.message.startsWith('impossible inverse in ')) throw error;
    return JSON.stringify(['inverse', error.message.split(':', 1)[0]]);
  }
}

functions.pari_buch_real_zero = pari_buch_real_zero;

functions.pari_buch_convert = (op: bigint, p: bigint, e: bigint, m: bigint, s: bigint, q: bigint, n: bigint) => {
  const x = {p:Number(p),e:Number(e),m,s:Number(s)};
  if (op === 12n) return JSON.stringify(['integer',String(nativeTruncateShift(x as MpReal, Number(n)))]);
  if (op >= 8n && op <= 10n) return JSON.stringify(['integer',String(op === 8n ? buchReal.truncr(x) : op === 9n ? buchReal.real_sign(x) : buchReal.real_expo(x))]);
  const z = op === 0n ? buchReal.itor(n, Number(q)) : op === 1n ? buchReal.real_0(Number(q))
    : op === 2n ? buchReal.real_1(Number(q)) : op === 3n ? buchReal.real_0_bit(Number(e))
    : op === 4n ? buchReal.setprec(x, Number(q)) : op === 5n ? buchReal.real_neg(x)
    : op === 6n ? buchReal.real_abs(x) : op === 7n ? buchReal.shiftr(x, Number(n)) : buchReal.mulur(n, x);
  return JSON.stringify(['real',z.s,String(z.e),String(z.m),z.p]);
};

functions.pari_real_truncate = (op: bigint,p: bigint,e: bigint,m: bigint,s: bigint,q: bigint,n: bigint) => JSON.stringify(['integer',String(qfbTruncate({p:Number(p),e:Number(e),m,s:Number(s) as -1|0|1}))]);

functions.pari_buch_real_sqrt = (op: bigint,p: bigint,e: bigint,m: bigint,s: bigint) => {
  const x = {p:Number(p),e:Number(e),m,s:Number(s)};
  const frame = (z: any): unknown[] => typeof z === 'bigint' ? ['integer',String(z)]
    : 're' in z ? ['complex',frame(z.re),frame(z.im)] : ['real',z.s,String(z.e),String(z.m),z.p];
  return JSON.stringify(frame(buchReal.sqrtr(x)));
};

functions.pari_real_logexp = (op: bigint,p: bigint,e: bigint,m: bigint,s: bigint,q: bigint) => {
  const x: MpReal = {p:Number(p),e:Number(e),m,s:Number(s) as -1|0|1};
  const z = op === 0n ? nativeLog(x) : op === 1n ? buchReal.logr_abs(x)
    : op === 2n ? nativeExp(x) : op === 3n ? buchReal.expr(x)
    : op === 4n ? sharedLog2(Number(q)) : buchReal.mplog2(Number(q));
  return JSON.stringify([z.s,String(z.e),String(z.m),z.p]);
};

functions.pari_binary_split = (op: bigint,n: bigint,a: bigint,b: bigint) => {
  if (!op) {
    const z = atanhuu(n,a,Number(b));
    return JSON.stringify([z.s,String(z.e),String(z.m),z.p]);
  }
  const size = Number(n), offset = Number(b), A = abpq_init(size + offset);
  let state = a;
  const next = (mod: bigint, shift: bigint) => {state = (1664525n * state + 1013904223n) & 0xffffffffn;return state % mod - shift;};
  for (let i = 0; i <= size + offset; i++) {A.a[i] = next(17n,8n);A.b[i] = next(19n,9n);A.p[i] = next(23n,11n);A.q[i] = next(29n,14n);}
  const R = abpq_sum(offset,offset + size,A);
  return JSON.stringify([R.P,R.Q,R.B,R.T].map(String));
};


functions.pari_real_exp1 = (p: bigint,e: bigint,m: bigint,s: bigint) => {
  const z = exp1r_abs({p:Number(p),e:Number(e),m,s:Number(s) as -1|0|1});
  return JSON.stringify([z.s,String(z.e),String(z.m),z.p]);
};
