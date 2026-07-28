/**
 * Unit tests for sage/rings/complex_mpfr
 * Tests for complex number transcendental functions
 */
import { describe, expect, test } from 'bun:test';
import { CC, ComplexField, type ComplexNumber } from './complex_mpfr.js';

// Helper to check if two complex numbers are approximately equal
function approxEqual(a: ComplexNumber, b: ComplexNumber, tol: number = 1e-10): boolean {
  return Math.abs(a.real() - b.real()) < tol && Math.abs(a.imag() - b.imag()) < tol;
}

// Helper to check if a number is approximately equal to a value
function numApproxEqual(a: number, b: number, tol: number = 1e-10): boolean {
  return Math.abs(a - b) < tol;
}

describe('ComplexField', () => {
  test('creation and basic properties', () => {
    const C = CC();
    expect(C.prec()).toBe(53);
    expect(C.precision()).toBe(53);
    expect(C.is_exact()).toBe(false);
    expect(C.characteristic()).toBe(0n);
    expect(C.ngens()).toBe(1);
  });

  test('pi', () => {
    const C = CC();
    const pi = C.pi();
    expect(numApproxEqual(pi.real(), Math.PI)).toBe(true);
    expect(pi.imag()).toBe(0);
  });

  test('zeta (roots of unity)', () => {
    const C = CC();

    // zeta(1) = 1
    const z1 = C.zeta(1);
    expect(numApproxEqual(z1.real(), 1)).toBe(true);
    expect(numApproxEqual(z1.imag(), 0)).toBe(true);

    // zeta(2) = -1
    const z2 = C.zeta(2);
    expect(numApproxEqual(z2.real(), -1)).toBe(true);
    expect(numApproxEqual(z2.imag(), 0)).toBe(true);

    // zeta(4) = i
    const z4 = C.zeta(4);
    expect(numApproxEqual(z4.real(), 0)).toBe(true);
    expect(numApproxEqual(z4.imag(), 1)).toBe(true);

    // zeta(3)^3 = 1
    const z3 = C.zeta(3);
    const z3cubed = z3.mul(z3).mul(z3);
    expect(numApproxEqual(z3cubed.real(), 1)).toBe(true);
    expect(numApproxEqual(z3cubed.imag(), 0)).toBe(true);
  });

  test('random_element', () => {
    const C = CC();
    const r = C.random_element(10);
    expect(Math.abs(r.real()) <= 10).toBe(true);
    expect(Math.abs(r.imag()) <= 10).toBe(true);
  });
});

describe('ComplexNumber basic operations', () => {
  test('creation', () => {
    const C = CC();
    const z = C.__call__(3, 4);
    expect(z.real()).toBe(3);
    expect(z.imag()).toBe(4);
  });

  test('arithmetic', () => {
    const C = CC();
    const z1 = C.__call__(3, 4);
    const z2 = C.__call__(1, 2);

    // Addition
    const sum = z1.add(z2);
    expect(sum.real()).toBe(4);
    expect(sum.imag()).toBe(6);

    // Subtraction
    const diff = z1.sub(z2);
    expect(diff.real()).toBe(2);
    expect(diff.imag()).toBe(2);

    // Multiplication: (3+4i)(1+2i) = 3 + 6i + 4i + 8i^2 = 3 + 10i - 8 = -5 + 10i
    const prod = z1.mul(z2);
    expect(prod.real()).toBe(-5);
    expect(prod.imag()).toBe(10);
  });

  test('norm and abs', () => {
    const C = CC();
    const z = C.__call__(3, 4);
    expect(z.norm()).toBe(25);
    expect(z.abs()).toBe(5);
  });

  test('conjugate', () => {
    const C = CC();
    const z = C.__call__(3, 4);
    const conj = z.conjugate();
    expect(conj.real()).toBe(3);
    expect(conj.imag()).toBe(-4);
  });
});

describe('Exponential and Logarithmic functions', () => {
  test('exp', () => {
    const C = CC();

    // exp(0) = 1
    const z0 = C.__call__(0, 0);
    const exp0 = z0.exp();
    expect(numApproxEqual(exp0.real(), 1)).toBe(true);
    expect(numApproxEqual(exp0.imag(), 0)).toBe(true);

    // exp(1) = e
    const z1 = C.__call__(1, 0);
    const exp1 = z1.exp();
    expect(numApproxEqual(exp1.real(), Math.E)).toBe(true);
    expect(numApproxEqual(exp1.imag(), 0)).toBe(true);

    // exp(i*pi) = -1 (Euler's identity)
    const iPi = C.__call__(0, Math.PI);
    const expIPi = iPi.exp();
    expect(numApproxEqual(expIPi.real(), -1)).toBe(true);
    expect(numApproxEqual(expIPi.imag(), 0)).toBe(true);

    // exp(1+i) test from SageMath docs
    // sage: (1+I).exp()
    // 1.46869393991589 + 2.28735528717884*I
    const z = C.__call__(1, 1);
    const expZ = z.exp();
    expect(numApproxEqual(expZ.real(), 1.46869393991589, 1e-10)).toBe(true);
    expect(numApproxEqual(expZ.imag(), 2.28735528717884, 1e-10)).toBe(true);
  });

  test('log', () => {
    const C = CC();

    // log(1) = 0
    const z1 = C.__call__(1, 0);
    const log1 = z1.log();
    expect(numApproxEqual(log1.real(), 0)).toBe(true);
    expect(numApproxEqual(log1.imag(), 0)).toBe(true);

    // log(e) = 1
    const ze = C.__call__(Math.E, 0);
    const loge = ze.log();
    expect(numApproxEqual(loge.real(), 1)).toBe(true);
    expect(numApproxEqual(loge.imag(), 0)).toBe(true);

    // log(-1) = i*pi
    const zm1 = C.__call__(-1, 0);
    const logm1 = zm1.log();
    expect(numApproxEqual(logm1.real(), 0)).toBe(true);
    expect(numApproxEqual(logm1.imag(), Math.PI)).toBe(true);

    // log(2+i)
    // sage: (2+I).log()
    // 0.80471895621705 + 0.463647609000806*I
    const z = C.__call__(2, 1);
    const logZ = z.log();
    expect(numApproxEqual(logZ.real(), 0.80471895621705, 1e-10)).toBe(true);
    expect(numApproxEqual(logZ.imag(), 0.463647609000806, 1e-10)).toBe(true);
  });

  test('log with base', () => {
    const C = CC();

    // log_2(4) = 2
    const z4 = C.__call__(4, 0);
    const log2_4 = z4.log(2);
    expect(numApproxEqual(log2_4.real(), 2)).toBe(true);
    expect(numApproxEqual(log2_4.imag(), 0)).toBe(true);
  });

  test('sqrt', () => {
    const C = CC();

    // sqrt(4) = 2
    const z4 = C.__call__(4, 0);
    const sqrt4 = z4.sqrt() as ComplexNumber;
    expect(numApproxEqual(sqrt4.real(), 2)).toBe(true);
    expect(numApproxEqual(sqrt4.imag(), 0)).toBe(true);

    // sqrt(-1) = i
    const zm1 = C.__call__(-1, 0);
    const sqrtm1 = zm1.sqrt() as ComplexNumber;
    expect(numApproxEqual(sqrtm1.real(), 0)).toBe(true);
    expect(numApproxEqual(sqrtm1.imag(), 1)).toBe(true);

    // sqrt(i) = (1+i)/sqrt(2)
    const i = C.__call__(0, 1);
    const sqrtI = i.sqrt() as ComplexNumber;
    const expected = 1 / Math.sqrt(2);
    expect(numApproxEqual(sqrtI.real(), expected)).toBe(true);
    expect(numApproxEqual(sqrtI.imag(), expected)).toBe(true);

    // sqrt(1+i) test from SageMath docs
    // sage: (1+I).sqrt()
    // 1.0986841 + 0.45508986*I
    const z = C.__call__(1, 1);
    const sqrtZ = z.sqrt() as ComplexNumber;
    expect(numApproxEqual(sqrtZ.real(), 1.0986841134678, 1e-6)).toBe(true);
    expect(numApproxEqual(sqrtZ.imag(), 0.4550898605622, 1e-6)).toBe(true);
  });

  test('sqrt all roots', () => {
    const C = CC();
    const z = C.__call__(4, 0);
    const roots = z.sqrt(true) as ComplexNumber[];
    expect(roots.length).toBe(2);
    expect(numApproxEqual(roots[0].real(), 2)).toBe(true);
    expect(numApproxEqual(roots[1].real(), -2)).toBe(true);
  });

  test('nth_root', () => {
    const C = CC();

    // 27^(1/3) = 3
    const z27 = C.__call__(27, 0);
    const cbrt27 = z27.nth_root(3) as ComplexNumber;
    expect(numApproxEqual(cbrt27.real(), 3)).toBe(true);
    expect(numApproxEqual(cbrt27.imag(), 0)).toBe(true);

    // All cube roots of 27
    const allRoots = z27.nth_root(3, true) as ComplexNumber[];
    expect(allRoots.length).toBe(3);

    // Verify each root cubed equals 27
    for (const root of allRoots) {
      const cubed = root.mul(root).mul(root);
      expect(numApproxEqual(cubed.real(), 27, 1e-8)).toBe(true);
      expect(numApproxEqual(cubed.imag(), 0, 1e-8)).toBe(true);
    }
  });
});

describe('Trigonometric functions', () => {
  test('sin', () => {
    const C = CC();

    // sin(0) = 0
    const z0 = C.__call__(0, 0);
    expect(numApproxEqual(z0.sin().real(), 0)).toBe(true);

    // sin(pi/2) = 1
    const zPi2 = C.__call__(Math.PI / 2, 0);
    expect(numApproxEqual(zPi2.sin().real(), 1)).toBe(true);

    // sin(1+i) test from SageMath docs
    // sage: (1+I).sin()
    // 1.29845758141598 + 0.634963914784736*I
    const z = C.__call__(1, 1);
    const sinZ = z.sin();
    expect(numApproxEqual(sinZ.real(), 1.29845758141598, 1e-10)).toBe(true);
    expect(numApproxEqual(sinZ.imag(), 0.634963914784736, 1e-10)).toBe(true);
  });

  test('cos', () => {
    const C = CC();

    // cos(0) = 1
    const z0 = C.__call__(0, 0);
    expect(numApproxEqual(z0.cos().real(), 1)).toBe(true);

    // cos(pi) = -1
    const zPi = C.__call__(Math.PI, 0);
    expect(numApproxEqual(zPi.cos().real(), -1)).toBe(true);

    // cos(1+i) test from SageMath docs
    // sage: (1+I).cos()
    // 0.833730025131149 - 0.988897705762865*I
    const z = C.__call__(1, 1);
    const cosZ = z.cos();
    expect(numApproxEqual(cosZ.real(), 0.833730025131149, 1e-10)).toBe(true);
    expect(numApproxEqual(cosZ.imag(), -0.988897705762865, 1e-10)).toBe(true);
  });

  test('tan', () => {
    const C = CC();

    // tan(0) = 0
    const z0 = C.__call__(0, 0);
    expect(numApproxEqual(z0.tan().real(), 0)).toBe(true);

    // tan(1+i) test from SageMath docs
    // sage: (1+I).tan()
    // 0.271752585319512 + 1.08392332733869*I
    const z = C.__call__(1, 1);
    const tanZ = z.tan();
    expect(numApproxEqual(tanZ.real(), 0.271752585319512, 1e-10)).toBe(true);
    expect(numApproxEqual(tanZ.imag(), 1.08392332733869, 1e-10)).toBe(true);
  });

  test('cot', () => {
    const C = CC();
    // cot(1+i) = 1/tan(1+i)
    // sage: (1+I).cot()
    // 0.217621561854403 - 0.868014142895925*I
    const z = C.__call__(1, 1);
    const cotZ = z.cot();
    expect(numApproxEqual(cotZ.real(), 0.217621561854403, 1e-10)).toBe(true);
    expect(numApproxEqual(cotZ.imag(), -0.868014142895925, 1e-10)).toBe(true);
  });

  test('sec', () => {
    const C = CC();
    // sec(1+i) = 1/cos(1+i)
    // sage: ComplexField(100)(1,1).sec()
    // 0.4983370305551868 + 0.591083841721045*I
    const z = C.__call__(1, 1);
    const secZ = z.sec();
    expect(numApproxEqual(secZ.real(), 0.498337030555187, 1e-10)).toBe(true);
    expect(numApproxEqual(secZ.imag(), 0.591083841721045, 1e-10)).toBe(true);
  });

  test('csc', () => {
    const C = CC();
    // csc(1+i) = 1/sin(1+i)
    // sage: ComplexField(100)(1,1).csc()
    // 0.6215180171704284 - 0.30393100162842646*I
    const z = C.__call__(1, 1);
    const cscZ = z.csc();
    expect(numApproxEqual(cscZ.real(), 0.621518017170428, 1e-10)).toBe(true);
    expect(numApproxEqual(cscZ.imag(), -0.303931001628426, 1e-10)).toBe(true);
  });
});

describe('Hyperbolic functions', () => {
  test('sinh', () => {
    const C = CC();

    // sinh(0) = 0
    const z0 = C.__call__(0, 0);
    expect(numApproxEqual(z0.sinh().real(), 0)).toBe(true);

    // sinh(1+i) test from SageMath docs
    // sage: (1+I).sinh()
    // 0.634963914784736 + 1.29845758141598*I
    const z = C.__call__(1, 1);
    const sinhZ = z.sinh();
    expect(numApproxEqual(sinhZ.real(), 0.634963914784736, 1e-10)).toBe(true);
    expect(numApproxEqual(sinhZ.imag(), 1.29845758141598, 1e-10)).toBe(true);
  });

  test('cosh', () => {
    const C = CC();

    // cosh(0) = 1
    const z0 = C.__call__(0, 0);
    expect(numApproxEqual(z0.cosh().real(), 1)).toBe(true);

    // cosh(1+i) test from SageMath docs
    // sage: (1+I).cosh()
    // 0.833730025131149 + 0.988897705762865*I
    const z = C.__call__(1, 1);
    const coshZ = z.cosh();
    expect(numApproxEqual(coshZ.real(), 0.833730025131149, 1e-10)).toBe(true);
    expect(numApproxEqual(coshZ.imag(), 0.988897705762865, 1e-10)).toBe(true);
  });

  test('tanh', () => {
    const C = CC();

    // tanh(0) = 0
    const z0 = C.__call__(0, 0);
    expect(numApproxEqual(z0.tanh().real(), 0)).toBe(true);

    // tanh(1+i) test from SageMath docs
    // sage: (1+I).tanh()
    // 1.08392332733869 + 0.271752585319512*I
    const z = C.__call__(1, 1);
    const tanhZ = z.tanh();
    expect(numApproxEqual(tanhZ.real(), 1.08392332733869, 1e-10)).toBe(true);
    expect(numApproxEqual(tanhZ.imag(), 0.271752585319512, 1e-10)).toBe(true);
  });

  test('coth', () => {
    const C = CC();
    // coth(1+i) = 1/tanh(1+i)
    // sage: ComplexField(100)(1,1).coth()
    // 0.8680141428959249 - 0.21762156185440268*I
    const z = C.__call__(1, 1);
    const cothZ = z.coth();
    expect(numApproxEqual(cothZ.real(), 0.868014142895925, 1e-10)).toBe(true);
    expect(numApproxEqual(cothZ.imag(), -0.217621561854403, 1e-10)).toBe(true);
  });

  test('sech', () => {
    const C = CC();
    // sech(1+i) = 1/cosh(1+i)
    // sage: ComplexField(100)(1,1).sech()
    // 0.4983370305551868 - 0.591083841721045*I
    const z = C.__call__(1, 1);
    const sechZ = z.sech();
    expect(numApproxEqual(sechZ.real(), 0.498337030555187, 1e-10)).toBe(true);
    expect(numApproxEqual(sechZ.imag(), -0.591083841721045, 1e-10)).toBe(true);
  });

  test('csch', () => {
    const C = CC();
    // csch(1+i) = 1/sinh(1+i)
    // sage: ComplexField(100)(1,1).csch()
    // 0.30393100162842646 - 0.6215180171704284*I
    const z = C.__call__(1, 1);
    const cschZ = z.csch();
    expect(numApproxEqual(cschZ.real(), 0.303931001628426, 1e-10)).toBe(true);
    expect(numApproxEqual(cschZ.imag(), -0.621518017170428, 1e-10)).toBe(true);
  });
});

describe('Inverse trigonometric functions', () => {
  test('arcsin', () => {
    const C = CC();

    // arcsin(0) = 0
    const z0 = C.__call__(0, 0);
    expect(numApproxEqual(z0.arcsin().real(), 0)).toBe(true);

    // arcsin(1+i) test from SageMath docs
    // sage: (1+I).arcsin()
    // 0.666239432492515 + 1.06127506190504*I
    const z = C.__call__(1, 1);
    const arcsinZ = z.arcsin();
    expect(numApproxEqual(arcsinZ.real(), 0.666239432492515, 1e-10)).toBe(true);
    expect(numApproxEqual(arcsinZ.imag(), 1.06127506190504, 1e-10)).toBe(true);
  });

  test('arccos', () => {
    const C = CC();

    // arccos(1) = 0
    const z1 = C.__call__(1, 0);
    expect(numApproxEqual(z1.arccos().real(), 0)).toBe(true);

    // arccos(1+i) test from SageMath docs
    // sage: (1+I).arccos()
    // 0.904556894302381 - 1.06127506190504*I
    const z = C.__call__(1, 1);
    const arccosZ = z.arccos();
    expect(numApproxEqual(arccosZ.real(), 0.904556894302381, 1e-10)).toBe(true);
    expect(numApproxEqual(arccosZ.imag(), -1.06127506190504, 1e-10)).toBe(true);
  });

  test('arctan', () => {
    const C = CC();

    // arctan(0) = 0
    const z0 = C.__call__(0, 0);
    expect(numApproxEqual(z0.arctan().real(), 0)).toBe(true);

    // arctan(1+i) test from SageMath docs
    // sage: (1+I).arctan()
    // 1.01722196789785 + 0.402359478108525*I
    const z = C.__call__(1, 1);
    const arctanZ = z.arctan();
    expect(numApproxEqual(arctanZ.real(), 1.01722196789785, 1e-10)).toBe(true);
    expect(numApproxEqual(arctanZ.imag(), 0.402359478108525, 1e-10)).toBe(true);
  });
});

describe('Inverse hyperbolic functions', () => {
  test('arcsinh', () => {
    const C = CC();

    // arcsinh(0) = 0
    const z0 = C.__call__(0, 0);
    expect(numApproxEqual(z0.arcsinh().real(), 0)).toBe(true);

    // arcsinh(1+i) test from SageMath docs
    // sage: (1+I).arcsinh()
    // 1.06127506190504 + 0.666239432492515*I
    const z = C.__call__(1, 1);
    const arcsinhZ = z.arcsinh();
    expect(numApproxEqual(arcsinhZ.real(), 1.06127506190504, 1e-10)).toBe(true);
    expect(numApproxEqual(arcsinhZ.imag(), 0.666239432492515, 1e-10)).toBe(true);
  });

  test('arccosh', () => {
    const C = CC();

    // arccosh(1) = 0
    const z1 = C.__call__(1, 0);
    expect(numApproxEqual(z1.arccosh().real(), 0)).toBe(true);

    // arccosh(1+i) test from SageMath docs
    // sage: (1+I).arccosh()
    // 1.06127506190504 + 0.904556894302381*I
    const z = C.__call__(1, 1);
    const arccoshZ = z.arccosh();
    expect(numApproxEqual(arccoshZ.real(), 1.06127506190504, 1e-10)).toBe(true);
    expect(numApproxEqual(arccoshZ.imag(), 0.904556894302381, 1e-10)).toBe(true);
  });

  test('arctanh', () => {
    const C = CC();

    // arctanh(0) = 0
    const z0 = C.__call__(0, 0);
    expect(numApproxEqual(z0.arctanh().real(), 0)).toBe(true);

    // arctanh(1+i) test from SageMath docs
    // sage: (1+I).arctanh()
    // 0.402359478108525 + 1.01722196789785*I
    const z = C.__call__(1, 1);
    const arctanhZ = z.arctanh();
    expect(numApproxEqual(arctanhZ.real(), 0.402359478108525, 1e-10)).toBe(true);
    expect(numApproxEqual(arctanhZ.imag(), 1.01722196789785, 1e-10)).toBe(true);
  });
});

describe('Special functions', () => {
  test('gamma basic cases', () => {
    const C = CC();

    // gamma(1) = 1
    const z1 = C.__call__(1, 0);
    const gamma1 = z1.gamma();
    expect(numApproxEqual(gamma1.real(), 1, 1e-6)).toBe(true);
    expect(numApproxEqual(gamma1.imag(), 0, 1e-6)).toBe(true);

    // gamma(2) = 1! = 1
    const z2 = C.__call__(2, 0);
    const gamma2 = z2.gamma();
    expect(numApproxEqual(gamma2.real(), 1, 1e-6)).toBe(true);

    // gamma(5) = 4! = 24
    const z5 = C.__call__(5, 0);
    const gamma5 = z5.gamma();
    expect(numApproxEqual(gamma5.real(), 24, 1e-6)).toBe(true);

    // gamma(1+i) test from SageMath docs (approximate)
    // sage: (1+I).gamma()
    // 0.498015668... - 0.154949828...*I
    const z = C.__call__(1, 1);
    const gammaZ = z.gamma();
    expect(numApproxEqual(gammaZ.real(), 0.498015668, 1e-6)).toBe(true);
    expect(numApproxEqual(gammaZ.imag(), -0.154949828, 1e-6)).toBe(true);
  });

  test('zeta basic cases', () => {
    const C = CC();

    // zeta(2) = pi^2/6
    const z2 = C.__call__(2, 0);
    const zeta2 = z2.zeta();
    expect(numApproxEqual(zeta2.real(), (Math.PI * Math.PI) / 6, 1e-6)).toBe(true);

    // zeta(1+i) test from SageMath docs (approximate)
    // sage: (1+I).zeta()
    // 0.582158... - 0.926848...*I
    const z = C.__call__(1, 1);
    const zetaZ = z.zeta();
    expect(numApproxEqual(zetaZ.real(), 0.582158, 1e-4)).toBe(true);
    expect(numApproxEqual(zetaZ.imag(), -0.926848, 1e-4)).toBe(true);
  });

  test('dilog', () => {
    const C = CC();

    // dilog(0) = 0
    const z0 = C.__call__(0, 0);
    const dilog0 = z0.dilog();
    expect(numApproxEqual(dilog0.real(), 0)).toBe(true);
    expect(numApproxEqual(dilog0.imag(), 0)).toBe(true);

    // dilog(1) = pi^2/6
    const z1 = C.__call__(1, 0);
    const dilog1 = z1.dilog();
    expect(numApproxEqual(dilog1.real(), (Math.PI * Math.PI) / 6, 1e-6)).toBe(true);
  });

  test('agm', () => {
    const C = CC();

    // AGM(a, a) = a
    const z = C.__call__(2, 0);
    const agmSame = z.agm(z);
    expect(numApproxEqual(agmSame.real(), 2, 1e-10)).toBe(true);

    // AGM(a, 0) = 0
    const z0 = C.__call__(0, 0);
    const agmZero = z.agm(z0);
    expect(numApproxEqual(agmZero.real(), 0)).toBe(true);

    // AGM(1, sqrt(2)) is known value
    const sqrt2 = C.__call__(Math.SQRT2, 0);
    const z1 = C.__call__(1, 0);
    const agmVal = z1.agm(sqrt2);
    // AGM(1, sqrt(2)) = pi / (2 * K(1/2)) approximately 1.1981
    expect(numApproxEqual(agmVal.real(), 1.19814023473559, 1e-6)).toBe(true);
  });

  test('eta', () => {
    const C = CC();

    // eta(i) should be defined (i is in upper half plane)
    const i = C.__call__(0, 1);
    const etaI = i.eta();
    expect(Number.isFinite(etaI.real())).toBe(true);
    expect(Number.isFinite(etaI.imag())).toBe(true);

    // eta(1+i) test from SageMath docs
    // sage: (1+I).eta()
    // 0.742048775836565 + 0.198831370229911*I
    const z = C.__call__(1, 1);
    const etaZ = z.eta();
    expect(numApproxEqual(etaZ.real(), 0.742048775836565, 1e-6)).toBe(true);
    expect(numApproxEqual(etaZ.imag(), 0.198831370229911, 1e-6)).toBe(true);
  });

  test('eta throws for lower half plane', () => {
    const C = CC();
    const z = C.__call__(1, -1);
    expect(() => z.eta()).toThrow('value must be in the upper half plane');
  });
});

describe('Properties and predicates', () => {
  test('is_real', () => {
    const C = CC();
    expect(C.__call__(3, 0).is_real()).toBe(true);
    expect(C.__call__(3, 1).is_real()).toBe(false);
  });

  test('is_imaginary', () => {
    const C = CC();
    expect(C.__call__(0, 3).is_imaginary()).toBe(true);
    expect(C.__call__(3, 0).is_imaginary()).toBe(false);
    expect(C.__call__(0, 0).is_imaginary()).toBe(false);
  });

  test('is_integer', () => {
    const C = CC();
    expect(C.__call__(3, 0).is_integer()).toBe(true);
    expect(C.__call__(3.5, 0).is_integer()).toBe(false);
    expect(C.__call__(3, 1).is_integer()).toBe(false);
  });

  test('is_square', () => {
    const C = CC();
    // All complex numbers are squares
    expect(C.__call__(3, 4).is_square()).toBe(true);
    expect(C.__call__(-1, 0).is_square()).toBe(true);
  });

  test('additive_order', () => {
    // sage: CC(0).additive_order()    -> 1
    // sage: CC.gen().additive_order() -> +Infinity
    const C = CC();
    expect(C.__call__(0, 0).additive_order()).toBe(1);
    expect(C.__call__(1, 0).additive_order()).toBe(Number.POSITIVE_INFINITY);
    expect(C.__call__(0, 1).additive_order()).toBe(Number.POSITIVE_INFINITY);
  });

  test('multiplicative_order', () => {
    // sage: C.<i> = ComplexField()
    // sage: i.multiplicative_order()      -> 4
    // sage: C(1).multiplicative_order()   -> 1
    // sage: C(-1).multiplicative_order()  -> 2
    // sage: C(-i).multiplicative_order()  -> 4
    // sage: C(2).multiplicative_order()   -> +Infinity
    const C = CC();
    expect(C.__call__(1, 0).multiplicative_order()).toBe(1);
    expect(C.__call__(-1, 0).multiplicative_order()).toBe(2);
    expect(C.__call__(0, 1).multiplicative_order()).toBe(4);
    expect(C.__call__(0, -1).multiplicative_order()).toBe(4);
    expect(C.__call__(2, 0).multiplicative_order()).toBe(Number.POSITIVE_INFINITY);
    expect(C.gen().multiplicative_order()).toBe(4);
  });

  test('multiplicative_order raises on the unit circle', () => {
    // sage: w = (1 + sqrt(-3))/2
    // sage: w.multiplicative_order()
    // NotImplementedError: order of element not known
    //
    // A floating point number on the unit circle carries no proof that it is a
    // root of unity, so SageMath refuses rather than brute-forcing n <= 1000.
    const C = CC();
    const w = C.__call__(0.5, Math.sqrt(3) / 2);
    expect(Math.abs(w.abs() - 1)).toBeLessThan(1e-15);
    expect(() => w.multiplicative_order()).toThrow('order of element not known');
    // e^{2 pi i / 5} is a genuine 5th root of unity, but still not known
    expect(() => C.zeta(5).multiplicative_order()).toThrow('order of element not known');
  });
});

describe('Identities', () => {
  test('Euler identity: e^(i*pi) + 1 = 0', () => {
    const C = CC();
    const iPi = C.__call__(0, Math.PI);
    const result = iPi.exp().add(C.__call__(1, 0));
    expect(numApproxEqual(result.real(), 0)).toBe(true);
    expect(numApproxEqual(result.imag(), 0)).toBe(true);
  });

  test('sin^2 + cos^2 = 1', () => {
    const C = CC();
    const z = C.__call__(1, 1);
    const sinZ = z.sin();
    const cosZ = z.cos();
    const sum = sinZ.mul(sinZ).add(cosZ.mul(cosZ));
    expect(numApproxEqual(sum.real(), 1, 1e-10)).toBe(true);
    expect(numApproxEqual(sum.imag(), 0, 1e-10)).toBe(true);
  });

  test('cosh^2 - sinh^2 = 1', () => {
    const C = CC();
    const z = C.__call__(1, 1);
    const sinhZ = z.sinh();
    const coshZ = z.cosh();
    const diff = coshZ.mul(coshZ).sub(sinhZ.mul(sinhZ));
    expect(numApproxEqual(diff.real(), 1, 1e-10)).toBe(true);
    expect(numApproxEqual(diff.imag(), 0, 1e-10)).toBe(true);
  });

  test('exp(log(z)) = z', () => {
    const C = CC();
    const z = C.__call__(2, 3);
    const result = z.log().exp();
    expect(numApproxEqual(result.real(), z.real(), 1e-10)).toBe(true);
    expect(numApproxEqual(result.imag(), z.imag(), 1e-10)).toBe(true);
  });

  test('sqrt(z)^2 = z', () => {
    const C = CC();
    const z = C.__call__(2, 3);
    const sqrtZ = z.sqrt() as ComplexNumber;
    const result = sqrtZ.mul(sqrtZ);
    expect(numApproxEqual(result.real(), z.real(), 1e-10)).toBe(true);
    expect(numApproxEqual(result.imag(), z.imag(), 1e-10)).toBe(true);
  });
});

describe('ComplexNumber - sqrt near the negative real axis (avoid_branch path)', () => {
  // Re(z) < 0 and |Im(z)| < |Re(z)| takes the "avoid_branch" path, which
  // computes sqrt(-z) and then divides by i.  Negating the real part instead of
  // just swapping the two parts returns sqrt(conj(z)).
  test('specific values', () => {
    const C = CC();
    const cases: Array<[number, number, number, number]> = [
      [-3, 1, 0.28484878459314106, 1.755317301824428],
      [-3, -1, 0.28484878459314106, -1.755317301824428],
      [-1, 0.5, 0.24293413587832283, 1.0290855136357462],
      [-5, 2, 0.4388421169022545, 2.27872385417085],
    ];
    for (const [re, im, er, ei] of cases) {
      const s = C.__call__(re, im).sqrt() as ComplexNumber;
      expect(numApproxEqual(s.real(), er, 1e-14)).toBe(true);
      expect(numApproxEqual(s.imag(), ei, 1e-14)).toBe(true);
    }
  });

  test('round trip over a grid covering both signs of Im', () => {
    const C = CC();
    for (let re = -5; re <= 5; re += 0.5) {
      for (let im = -5; im <= 5; im += 0.5) {
        if (re === 0 && im === 0) continue;
        const z = C.__call__(re, im);
        const s = z.sqrt() as ComplexNumber;
        const back = s.mul(s);
        const scale = Math.max(1, z.abs());
        expect(Math.abs(back.real() - re)).toBeLessThan(1e-13 * scale);
        expect(Math.abs(back.imag() - im)).toBeLessThan(1e-13 * scale);
        // Principal branch: Re(sqrt(z)) >= 0.
        expect(s.real()).toBeGreaterThanOrEqual(-1e-15);
      }
    }
  });

  test('sqrt(all=True) returns both roots', () => {
    const C = CC();
    const roots = C.__call__(-3, 1).sqrt(true) as ComplexNumber[];
    expect(roots.length).toBe(2);
    for (const r of roots) {
      const sq = r.mul(r);
      expect(numApproxEqual(sq.real(), -3, 1e-13)).toBe(true);
      expect(numApproxEqual(sq.imag(), 1, 1e-13)).toBe(true);
    }
  });

  test('arccos stays in the principal strip', () => {
    // arccos calls sqrt(1 - z^2); a wrong branch pushes Re(arccos) out of [0, pi].
    const C = CC();
    for (const [re, im] of [
      [-2, 0.5],
      [-2, -0.5],
      [-3, 1],
      [2, 0.5],
      [0.5, -3],
    ] as Array<[number, number]>) {
      const z = C.__call__(re, im);
      const w = z.arccos();
      expect(w.real()).toBeGreaterThanOrEqual(0);
      expect(w.real()).toBeLessThanOrEqual(Math.PI);
      // cos(arccos(z)) == z
      const back = w.cos();
      expect(numApproxEqual(back.real(), re, 1e-12)).toBe(true);
      expect(numApproxEqual(back.imag(), im, 1e-12)).toBe(true);
    }
  });
});

describe('ComplexNumber - dilog', () => {
  // Reference values from PARI/mpmath (SageMath delegates to PARI's dilog).
  const cases: Array<[number, number, number, number]> = [
    [0, 1, -0.2056167583560283, 0.915965594177219],
    [-1, 0, -0.8224670334241132, 0],
    [0.5, 0, 0.5822405264650125, 0],
    [0.1, 0, 0.10261779109939113, 0],
    // |z| > 1: these used to overflow to ~1e13 / ~1e26 through a divergent series
    [-1.5, 0, -1.1473806603755707, 0],
    [2, 0, 2.4674011002723395, -2.177586090303602],
    [0, 1.5, -0.3927071122175517, 1.27496944849438],
    [-3, 2, -2.0713071652315143, 0.8922731679007034],
    [5, -2, 0.5705771689906514, -4.724338407765386],
    [1, 1, 0.6168502750680849, 1.4603621167531196],
    // 0.5 < |z| <= 1
    [0.9, 0.1, 1.264186732338754, 0.24373567998101406],
    [0.99, 0, 1.5886254480763753, 0],
    [0.3, 0.4, 0.2665968667427404, 0.461362891819109],
  ];

  test('special values', () => {
    const C = CC();
    // sage: ComplexNumber(1,0).dilog() -> 1.64493406684823 (= pi^2/6)
    const d1 = C.__call__(1, 0).dilog();
    expect(numApproxEqual(d1.real(), Math.PI ** 2 / 6, 1e-14)).toBe(true);
    expect(d1.imag()).toBe(0);
    // sage: ComplexNumber(0,0).dilog() -> 0
    const d0 = C.__call__(0, 0).dilog();
    expect(d0.real()).toBe(0);
    expect(d0.imag()).toBe(0);
    // sage: ComplexNumber(0,1).dilog() -> -0.205616758356028 + 0.915965594177219*I
    const di = C.__call__(0, 1).dilog();
    expect(numApproxEqual(di.real(), -0.205616758356028, 1e-14)).toBe(true);
    expect(numApproxEqual(di.imag(), 0.915965594177219, 1e-14)).toBe(true);
  });

  test('matches PARI across all branches', () => {
    const C = CC();
    for (const [re, im, er, ei] of cases) {
      const g = C.__call__(re, im).dilog();
      const scale = Math.max(1, Math.hypot(er, ei));
      expect(Math.abs(g.real() - er)).toBeLessThan(1e-13 * scale);
      expect(Math.abs(g.imag() - ei)).toBeLessThan(1e-13 * scale);
    }
  });

  test('satisfies the inversion identity', () => {
    // Li2(z) + Li2(1/z) = -pi^2/6 - (1/2) log(-z)^2
    const C = CC();
    for (const [re, im] of [
      [3, 1],
      [-4, 2],
      [2.5, -3.5],
      [-1.5, -0.25],
    ] as Array<[number, number]>) {
      const z = C.__call__(re, im);
      const lhs = z.dilog().add(z.inv().dilog());
      const logNegZ = C.__call__(-re, im === 0 ? 0 : -im).log();
      const rhs = C.__call__(-(Math.PI ** 2) / 6, 0).sub(
        logNegZ.mul(logNegZ).mul(C.__call__(0.5, 0))
      );
      expect(numApproxEqual(lhs.real(), rhs.real(), 1e-12)).toBe(true);
      expect(numApproxEqual(lhs.imag(), rhs.imag(), 1e-12)).toBe(true);
    }
  });
});

describe('ComplexNumber - gamma_inc', () => {
  // Reference values from PARI's incgam (which SageMath delegates to).
  test('SageMath doctests', () => {
    const C = CC();
    const cases: Array<[number, number, number, number, number, number]> = [
      // sage: (1+i).gamma_inc(2 + 3*i) -> 0.0020969149 - 0.059981914*I
      [1, 1, 2, 3, 0.0020969148636468277, -0.059981913655449706],
      // sage: (1+i).gamma_inc(5) -> -0.0013781309 + 0.00651982*I
      [1, 1, 5, 0, -0.001378130936215849, 0.006519820023119819],
      // sage: C(2).gamma_inc(1 + i) -> 0.7070921 - 0.42035364*I
      [2, 0, 1, 1, 0.7070920963459381, -0.4203536409598115],
      // sage: CC(2).gamma_inc(5) -> 0.0404276819945128
      [2, 0, 5, 0, 0.040427681994512805, 0],
      // sage: C(2 + I).gamma_inc(C(3 + I))
      [2, 1, 3, 1, 0.12151564466450869, 0.10153390907982604],
      [3, -1, 0.5, 0.5, 1.1078103895745073, -1.4111479439739232],
      [0.5, 0, 0.3, 0, 0.7773593112498081, 0],
    ];
    for (const [ar, ai, tr, ti, er, ei] of cases) {
      const g = C.__call__(ar, ai).gamma_inc(C.__call__(tr, ti));
      const scale = Math.hypot(er, ei);
      expect(Math.abs(g.real() - er)).toBeLessThan(1e-12 * scale);
      expect(Math.abs(g.imag() - ei)).toBeLessThan(1e-12 * scale);
    }
  });

  test('keeps relative accuracy for large |t|', () => {
    // Gamma(2, 50) = 51 * e^-50 = 9.8366242246e-21.  Subtracting a truncated
    // lower-gamma series from Gamma(2) = 1 loses 10 orders of magnitude.
    const C = CC();
    const g = C.__call__(2, 0).gamma_inc(C.__call__(50, 0));
    expect(Math.abs(g.real() / 9.83662422461598e-21 - 1)).toBeLessThan(1e-12);
    expect(Math.abs(g.imag())).toBeLessThan(1e-30);
    // Gamma(1, t) = e^-t exactly
    for (const t of [3, 10, 30, 80]) {
      const e = C.__call__(1, 0).gamma_inc(C.__call__(t, 0));
      expect(Math.abs(e.real() / Math.exp(-t) - 1)).toBeLessThan(1e-12);
    }
  });

  test('gamma_inc(a, 0) = gamma(a)', () => {
    const C = CC();
    const a = C.__call__(2.5, 0);
    const g = a.gamma_inc(C.__call__(0, 0));
    expect(numApproxEqual(g.real(), a.gamma().real(), 1e-12)).toBe(true);
  });
});

describe('ComplexNumber - algebraic_dependency', () => {
  test('primitive 6th root of unity', () => {
    // sage: C = ComplexField()
    // sage: z = (1/2)*(1 + sqrt(3)*C.0)
    // sage: p = z.algebraic_dependency(5); p
    // x^2 - x + 1
    // sage: p(z)
    // 1.11022302462516e-16
    const C = CC();
    const z = C.__call__(0.5, Math.sqrt(3) / 2);
    for (const degree of [2, 3, 5, 6]) {
      expect(z.algebraic_dependency(degree)).toEqual([1n, -1n, 1n]);
    }

    // p(z) must actually vanish
    const p = z.algebraic_dependency(5);
    let acc = C.__call__(0, 0);
    for (let i = p.length - 1; i >= 0; i--) {
      acc = acc.mul(z).add(C.__call__(Number(p[i]!), 0));
    }
    expect(acc.abs()).toBeLessThan(1e-14);
  });

  test('other algebraic numbers', () => {
    const C = CC();
    // sage: algdep(complex("1+2j"), 4) -> x^2 - 2*x + 5
    expect(C.__call__(1, 2).algebraic_dependency(4)).toEqual([5n, -2n, 1n]);
    // i is a root of x^2 + 1
    expect(C.__call__(0, 1).algebraic_dependency(4)).toEqual([1n, 0n, 1n]);
    // sqrt(2) is a root of x^2 - 2
    expect(C.__call__(Math.SQRT2, 0).algebraic_dependency(4)).toEqual([-2n, 0n, 1n]);
    // Integers give x - n
    expect(C.__call__(3, 0).algebraic_dependency(2)).toEqual([-3n, 1n]);
    expect(C.__call__(-7, 0).algebraic_dependency(5)).toEqual([7n, 1n]);
  });

  test('the returned polynomial vanishes at z', () => {
    const C = CC();
    const values: Array<[number, number]> = [
      [Math.SQRT2 / 2, Math.SQRT2 / 2],
      [(1 + Math.sqrt(5)) / 2, 0],
      [0, Math.cbrt(2)],
    ];
    for (const [re, im] of values) {
      const z = C.__call__(re, im);
      const p = z.algebraic_dependency(6);
      let acc = C.__call__(0, 0);
      for (let i = p.length - 1; i >= 0; i--) {
        acc = acc.mul(z).add(C.__call__(Number(p[i]!), 0));
      }
      expect(acc.abs()).toBeLessThan(1e-8);
      // Non-constant
      expect(p.slice(1).some((c) => c !== 0n)).toBe(true);
    }
  });
});
