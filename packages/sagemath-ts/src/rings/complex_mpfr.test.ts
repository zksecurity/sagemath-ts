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
    // 0.804718956217050 + 0.463647609000806*I
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
    // 0.49833703055518678521380589177 + 0.59108384172104504805039169297*I
    const z = C.__call__(1, 1);
    const secZ = z.sec();
    expect(numApproxEqual(secZ.real(), 0.498337030555187, 1e-10)).toBe(true);
    expect(numApproxEqual(secZ.imag(), 0.591083841721045, 1e-10)).toBe(true);
  });

  test('csc', () => {
    const C = CC();
    // csc(1+i) = 1/sin(1+i)
    // sage: ComplexField(100)(1,1).csc()
    // 0.62151801717042842123490780586 - 0.30393100162842645033448560451*I
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
    // 0.86801414289592494863584920892 - 0.21762156185440268136513424361*I
    const z = C.__call__(1, 1);
    const cothZ = z.coth();
    expect(numApproxEqual(cothZ.real(), 0.868014142895925, 1e-10)).toBe(true);
    expect(numApproxEqual(cothZ.imag(), -0.217621561854403, 1e-10)).toBe(true);
  });

  test('sech', () => {
    const C = CC();
    // sech(1+i) = 1/cosh(1+i)
    // sage: ComplexField(100)(1,1).sech()
    // 0.49833703055518678521380589177 - 0.59108384172104504805039169297*I
    const z = C.__call__(1, 1);
    const sechZ = z.sech();
    expect(numApproxEqual(sechZ.real(), 0.498337030555187, 1e-10)).toBe(true);
    expect(numApproxEqual(sechZ.imag(), -0.591083841721045, 1e-10)).toBe(true);
  });

  test('csch', () => {
    const C = CC();
    // csch(1+i) = 1/sinh(1+i)
    // sage: ComplexField(100)(1,1).csch()
    // 0.30393100162842645033448560451 - 0.62151801717042842123490780586*I
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
    const C = CC();
    expect(C.__call__(0, 0).additive_order()).toBe(1);
    expect(C.__call__(1, 0).additive_order()).toBe('Infinity');
    expect(C.__call__(0, 1).additive_order()).toBe('Infinity');
  });

  test('multiplicative_order', () => {
    const C = CC();
    expect(C.__call__(1, 0).multiplicative_order()).toBe(1);
    expect(C.__call__(-1, 0).multiplicative_order()).toBe(2);
    expect(C.__call__(2, 0).multiplicative_order()).toBe('Infinity');
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
