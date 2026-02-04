/**
 * # Gate Identity Check in HyperPlonk
 *
 * The gate identity check verifies that all gates in the circuit
 * are correctly computed. This is done via ZeroCheck on the
 * gate constraint polynomial.
 *
 * ## The Gate Constraint
 *
 * At each gate b in the boolean hypercube, we have:
 *   q_L(b)*w_1(b) + q_R(b)*w_2(b) + q_O(b)*w_3(b)
 *   + q_M(b)*w_1(b)*w_2(b) + q_C(b) = 0
 *
 * This constraint must hold for ALL gates.
 *
 * ## Using ZeroCheck
 *
 * We reduce the gate check to ZeroCheck:
 * - Define G(x) = gate constraint as a polynomial
 * - Prove G(b) = 0 for all b in {0,1}^mu
 *
 * @module tutorial/part7-hyperplonk/23-hyperplonk/gate-identity
 */

import { NotImplementedError } from 'sagemath-ts';

console.log("=== Gate Identity Check ===\n");

// ============================================================================
// Field Arithmetic (simple modular arithmetic)
// ============================================================================

const p = 97n;

const mod = (x: bigint): bigint => ((x % p) + p) % p;
const add = (a: bigint, b: bigint): bigint => mod(a + b);
const sub = (a: bigint, b: bigint): bigint => mod(a - b);
const mul = (a: bigint, b: bigint): bigint => mod(a * b);

interface HypercubeFunction {
  numVars: number;
  values: bigint[];
}

function booleanHypercube(n: number): number[][] {
  if (n === 0) return [[]];
  const smaller = booleanHypercube(n - 1);
  return [...smaller.map(v => [0, ...v]), ...smaller.map(v => [1, ...v])];
}

function eq(x: bigint[], b: bigint[]): bigint {
  let result = 1n;
  for (let i = 0; i < x.length; i++) {
    result = mul(result, add(mul(x[i]!, b[i]!), mul(sub(1n, x[i]!), sub(1n, b[i]!))));
  }
  return result;
}

function evaluateMLE(f: HypercubeFunction, x: bigint[]): bigint {
  const hypercube = booleanHypercube(f.numVars);
  let result = 0n;
  for (let i = 0; i < hypercube.length; i++) {
    const b = hypercube[i]!.map(v => BigInt(v));
    result = add(result, mul(f.values[i]!, eq(x, b)));
  }
  return result;
}

// ============================================================================
// Section 1: The Gate Constraint Polynomial
// ============================================================================

console.log("--- Section 1: Gate Constraint Polynomial ---\n");

/**
 * GATE CONSTRAINT POLYNOMIAL
 *
 * For a PLONK-style gate, the constraint at position b is:
 *
 *   G(b) = q_L(b)*w_1(b) + q_R(b)*w_2(b) + q_O(b)*w_3(b)
 *        + q_M(b)*w_1(b)*w_2(b) + q_C(b)
 *
 * where:
 * - w_1, w_2, w_3: wire polynomials (left, right, output)
 * - q_L, q_R, q_O: selector polynomials for linear terms
 * - q_M: selector for multiplication term
 * - q_C: constant term
 *
 * DEGREE ANALYSIS:
 *
 * - Each w_i is multilinear (degree 1 in each variable)
 * - Each q_* is multilinear
 * - Product q_M * w_1 * w_2 has degree 3 in each variable
 * - Overall G has degree at most 3 in each variable
 *
 * For ZeroCheck on G:
 * - h(x) = eq(x, r) * G(x) has degree 4
 * - Each round sends 5 evaluations (degree 4 polynomial)
 */

console.log("Gate constraint polynomial G(b):\n");
console.log("  G(b) = q_L(b)*w_1(b) + q_R(b)*w_2(b) + q_O(b)*w_3(b)");
console.log("       + q_M(b)*w_1(b)*w_2(b) + q_C(b)\n");

console.log("Degree analysis:");
console.log("  - Linear terms: degree 2 (q * w)");
console.log("  - Multiplication term: degree 3 (q_M * w_1 * w_2)");
console.log("  - Overall: degree 3 per variable\n");

console.log("ZeroCheck on G:");
console.log("  - h(x) = eq(x, r) * G(x) has degree 4");
console.log("  - Each round: 5 field elements\n");

// ============================================================================
// Section 2: Example Circuit
// ============================================================================

console.log("--- Section 2: Example Circuit ---\n");

/**
 * EXAMPLE CIRCUIT
 *
 * Let's create a simple circuit with 4 gates (mu = 2):
 *
 * Gate 0: w_1 + w_2 - w_3 = 0  (addition)
 *         q_L=1, q_R=1, q_O=-1, q_M=0, q_C=0
 *
 * Gate 1: w_1 * w_2 - w_3 = 0  (multiplication)
 *         q_L=0, q_R=0, q_O=-1, q_M=1, q_C=0
 *
 * Gate 2: w_1 - 5 = 0  (constant check)
 *         q_L=1, q_R=0, q_O=0, q_M=0, q_C=-5
 *
 * Gate 3: w_3 - 40 = 0  (output check)
 *         q_L=0, q_R=0, q_O=1, q_M=0, q_C=-40
 *
 * Valid assignment:
 *   Gate 0: w_1=5, w_2=3, w_3=8 (5 + 3 = 8)
 *   Gate 1: w_1=8, w_2=5, w_3=40 (8 * 5 = 40)
 *   Gate 2: w_1=5 (check w_1 = 5)
 *   Gate 3: w_3=40 (check output = 40)
 */

const mu = 2;
const n = 2 ** mu;  // 4 gates
const cube = booleanHypercube(mu);

// Selector polynomials
const q_L: HypercubeFunction = { numVars: mu, values: [1n, 0n, 1n, 0n] };
const q_R: HypercubeFunction = { numVars: mu, values: [1n, 0n, 0n, 0n] };
const q_O: HypercubeFunction = { numVars: mu, values: [mod(-1n), mod(-1n), 0n, 1n] };  // -1 mod p
const q_M: HypercubeFunction = { numVars: mu, values: [0n, 1n, 0n, 0n] };
const q_C: HypercubeFunction = { numVars: mu, values: [0n, 0n, mod(-5n), mod(-40n)] };  // -5, -40 mod p

// Wire polynomials (valid assignment)
const w_1: HypercubeFunction = { numVars: mu, values: [5n, 8n, 5n, 0n] };
const w_2: HypercubeFunction = { numVars: mu, values: [3n, 5n, 0n, 0n] };
const w_3: HypercubeFunction = { numVars: mu, values: [8n, 40n, 0n, 40n] };

console.log("Example circuit (4 gates):\n");

console.log("Gate 0 (add): w_1 + w_2 - w_3 = 0");
console.log("  w_1=5, w_2=3, w_3=8 => 5 + 3 - 8 = 0 OK\n");

console.log("Gate 1 (mul): w_1 * w_2 - w_3 = 0");
console.log("  w_1=8, w_2=5, w_3=40 => 8*5 - 40 = 0 OK\n");

console.log("Gate 2 (const): w_1 - 5 = 0");
console.log("  w_1=5 => 5 - 5 = 0 OK\n");

console.log("Gate 3 (output): w_3 - 40 = 0");
console.log("  w_3=40 => 40 - 40 = 0 OK\n");

// ============================================================================
// Section 3: Evaluating the Gate Constraint
// ============================================================================

console.log("--- Section 3: Evaluating Gate Constraint ---\n");

/**
 * Evaluates the gate constraint at a point.
 * G(x) = q_L(x)*w_1(x) + q_R(x)*w_2(x) + q_O(x)*w_3(x) + q_M(x)*w_1(x)*w_2(x) + q_C(x)
 */
function evaluateGateConstraint(x: bigint[]): bigint {
  const qL = evaluateMLE(q_L, x);
  const qR = evaluateMLE(q_R, x);
  const qO = evaluateMLE(q_O, x);
  const qM = evaluateMLE(q_M, x);
  const qC = evaluateMLE(q_C, x);

  const w1 = evaluateMLE(w_1, x);
  const w2 = evaluateMLE(w_2, x);
  const w3 = evaluateMLE(w_3, x);

  // G = q_L*w_1 + q_R*w_2 + q_O*w_3 + q_M*w_1*w_2 + q_C
  return add(
    add(
      add(mul(qL, w1), mul(qR, w2)),
      add(mul(qO, w3), mul(mul(qM, w1), w2))
    ),
    qC
  );
}

// Verify constraint on boolean points
console.log("Gate constraint G(b) on boolean points:");
let allZero = true;
for (let i = 0; i < cube.length; i++) {
  const b = cube[i]!.map(v => BigInt(v));
  const G = evaluateGateConstraint(b);
  const status = G === 0n ? 'OK' : 'FAIL';
  console.log(`  G(${cube[i]!.join(',')}) = ${G} ${status}`);
  if (G !== 0n) allZero = false;
}
console.log(`\nAll gates satisfied? ${allZero ? 'YES' : 'NO'}\n`);

// ============================================================================
// Section 4: ZeroCheck for Gate Identity
// ============================================================================

console.log("--- Section 4: ZeroCheck for Gate Identity ---\n");

/**
 * GATE IDENTITY CHECK VIA ZEROCHECK
 *
 * To prove all gates are satisfied:
 * 1. Verifier sends random r in F^mu
 * 2. Define h(x) = eq(x, r) * G(x)
 * 3. Run SumCheck to prove sum_{b in B_mu} h(b) = 0
 *
 * Since G(b) = 0 for all b, we have h(b) = 0 for all b,
 * so the sum is 0.
 *
 * If any G(b) != 0, then with high probability
 * sum h(b) != 0, and the prover is caught.
 */

console.log("ZeroCheck for gate identity:\n");

console.log("  1. Verifier sends random r");
console.log("  2. Define h(x) = eq(x, r) * G(x)");
console.log("  3. SumCheck: prove sum_{b} h(b) = 0");
console.log("  4. Final check: h(s) = eq(s, r) * G(s)\n");

// Demonstrate with a random r
const r = [3n, 7n];
console.log(`Random challenge: r = (${r.join(', ')})\n`);

// Compute sum_{b} eq(b, r) * G(b)
let hSum = 0n;
console.log("Computing sum_{b} eq(b, r) * G(b):");
for (let i = 0; i < cube.length; i++) {
  const b = cube[i]!.map(v => BigInt(v));
  const eqVal = eq(b, r);
  const G = evaluateGateConstraint(b);
  const h = mul(eqVal, G);
  console.log(`  b=${cube[i]!.join(',')}: eq(b,r)=${eqVal}, G(b)=${G}, h(b)=${h}`);
  hSum = add(hSum, h);
}
console.log(`\nSum = ${hSum} (should be 0 for valid circuit)\n`);

// ============================================================================
// Section 5: Detecting Invalid Gate
// ============================================================================

console.log("--- Section 5: Detecting Invalid Gate ---\n");

// Create bad wire assignment
const w_3_bad: HypercubeFunction = {
  numVars: mu,
  values: [8n, 41n, 0n, 40n]  // Gate 1 wrong: 8*5=40, not 41
};

console.log("Bad assignment: w_3(gate 1) = 41 instead of 40\n");

function evaluateGateConstraintBad(x: bigint[]): bigint {
  const qL = evaluateMLE(q_L, x);
  const qR = evaluateMLE(q_R, x);
  const qO = evaluateMLE(q_O, x);
  const qM = evaluateMLE(q_M, x);
  const qC = evaluateMLE(q_C, x);

  const w1 = evaluateMLE(w_1, x);
  const w2 = evaluateMLE(w_2, x);
  const w3 = evaluateMLE(w_3_bad, x);  // Using bad w_3

  return add(
    add(
      add(mul(qL, w1), mul(qR, w2)),
      add(mul(qO, w3), mul(mul(qM, w1), w2))
    ),
    qC
  );
}

// Verify on boolean points
console.log("Gate constraint with bad assignment:");
for (let i = 0; i < cube.length; i++) {
  const b = cube[i]!.map(v => BigInt(v));
  const G = evaluateGateConstraintBad(b);
  const status = G === 0n ? 'OK' : 'FAIL';
  console.log(`  G(${cube[i]!.join(',')}) = ${G} ${status}`);
}

// ZeroCheck would catch this
console.log("\nZeroCheck sum with bad assignment:");
let hSumBad = 0n;
for (let i = 0; i < cube.length; i++) {
  const b = cube[i]!.map(v => BigInt(v));
  const eqVal = eq(b, r);
  const G = evaluateGateConstraintBad(b);
  hSumBad = add(hSumBad, mul(eqVal, G));
}
console.log(`  Sum = ${hSumBad} (NOT zero => caught!)\n`);

// ============================================================================
// Section 6: Custom Gates
// ============================================================================

console.log("--- Section 6: Custom Gates ---\n");

/**
 * CUSTOM GATES IN HYPERPLONK
 *
 * The standard PLONK gate is:
 *   q_L*w_1 + q_R*w_2 + q_O*w_3 + q_M*w_1*w_2 + q_C = 0
 *
 * But we can use ANY polynomial constraint!
 *
 * EXAMPLE: Boolean gate
 *   w_1 * (1 - w_1) = 0
 *   Enforces w_1 in {0, 1}
 *
 * EXAMPLE: Range gate (w in {0, 1, 2, 3})
 *   w * (w-1) * (w-2) * (w-3) = 0
 *   Degree 4
 *
 * EXAMPLE: Poseidon S-box
 *   w_out = w_in^5
 *   Or: w_out - w_in^5 = 0
 *
 * IMPLEMENTATION:
 *
 * Use a selector q_custom that activates the custom constraint:
 *   q_custom(b) * (custom constraint) = 0
 *
 * When q_custom = 0, the constraint is inactive.
 * When q_custom = 1, the constraint must hold.
 */

console.log("Custom gate examples:\n");

console.log("Boolean: w * (1 - w) = 0");
console.log("  Degree 2, enforces w in {0, 1}\n");

console.log("Range [0,4): w * (w-1) * (w-2) * (w-3) = 0");
console.log("  Degree 4, enforces w in {0, 1, 2, 3}\n");

console.log("Poseidon S-box: w_out - w_in^5 = 0");
console.log("  Degree 5, implements x^5 permutation\n");

// Boolean constraint example
const w_bool: HypercubeFunction = {
  numVars: mu,
  values: [0n, 1n, 1n, 0n]  // All 0 or 1
};

console.log("Boolean constraint check:");
for (let i = 0; i < cube.length; i++) {
  const w = w_bool.values[i]!;
  const constraint = mul(w, sub(1n, w));  // w * (1 - w)
  console.log(`  w(${cube[i]!.join(',')}) = ${w}: constraint = ${constraint}`);
}
console.log();

// ============================================================================
// Section 7: Implementation Stubs
// ============================================================================

console.log("--- Section 7: Implementation Stubs ---\n");

/**
 * Complete gate identity check using ZeroCheck.
 *
 * @throws NotImplementedError - Full implementation pending
 */
export function gateIdentityCheck(
  selectors: { qL: HypercubeFunction; qR: HypercubeFunction; qO: HypercubeFunction; qM: HypercubeFunction; qC: HypercubeFunction },
  wires: { w1: HypercubeFunction; w2: HypercubeFunction; w3: HypercubeFunction },
  r: bigint[],
  getChallenge: () => bigint
): { accepted: boolean } {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: gateIdentityCheck - ' +
    'Full gate identity check with ZeroCheck'
  );
}

/**
 * Custom gate check with arbitrary polynomial constraint.
 *
 * @throws NotImplementedError - Custom gate implementation pending
 */
export function customGateCheck(
  constraint: (wires: bigint[]) => bigint,
  degree: number,
  wires: HypercubeFunction[],
  r: bigint[],
  getChallenge: () => bigint
): { accepted: boolean } {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: customGateCheck - ' +
    'Custom gate constraint verification'
  );
}

console.log("Stub functions:");
console.log("  - gateIdentityCheck: full gate check");
console.log("  - customGateCheck: arbitrary constraints\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("Gate Identity Check:");
console.log("  Constraint: G(b) = sum of q*w terms = 0 at each gate");
console.log("  Method: ZeroCheck on G\n");

console.log("Key points:");
console.log("  - G has degree 3 (from q_M*w_1*w_2)");
console.log("  - ZeroCheck h(x) = eq(x,r)*G(x) has degree 4");
console.log("  - Custom gates: just change the constraint polynomial\n");

console.log("Next: Wiring identity check");

// Exports
export {
  evaluateGateConstraint,
  type HypercubeFunction
};
