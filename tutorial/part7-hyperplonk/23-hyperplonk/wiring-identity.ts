/**
 * # Wiring Identity Check in HyperPlonk
 *
 * The wiring identity check ensures that wire values are correctly
 * "copied" between different positions in the circuit. In PLONK-style
 * systems, this is called the copy constraint or permutation argument.
 *
 * ## The Problem
 *
 * In a circuit, the output of one gate often becomes the input of another.
 * For example:
 *   - Gate i has output w_3(i)
 *   - Gate j uses this value as input w_1(j)
 *   - We need w_3(i) = w_1(j)
 *
 * ## The Solution
 *
 * Encode all such equalities as a permutation and use ProductCheck
 * to verify the permutation is respected.
 *
 * @module tutorial/part7-hyperplonk/23-hyperplonk/wiring-identity
 */

import { NotImplementedError } from 'sagemath-ts';

console.log("=== Wiring Identity Check ===\n");

// ============================================================================
// Field Arithmetic (simple modular arithmetic)
// ============================================================================

const p = 97n;

const mod = (x: bigint): bigint => ((x % p) + p) % p;
const add = (a: bigint, b: bigint): bigint => mod(a + b);
const sub = (a: bigint, b: bigint): bigint => mod(a - b);
const mul = (a: bigint, b: bigint): bigint => mod(a * b);

function powMod(base: bigint, exp: bigint, modulus: bigint): bigint {
  let result = 1n;
  base = ((base % modulus) + modulus) % modulus;
  while (exp > 0n) {
    if (exp % 2n === 1n) {
      result = (result * base) % modulus;
    }
    exp = exp / 2n;
    base = (base * base) % modulus;
  }
  return result;
}

const inv = (a: bigint): bigint => powMod(a, p - 2n, p);

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
    // eq_i = x_i * b_i + (1 - x_i) * (1 - b_i)
    const term = add(mul(x[i]!, b[i]!), mul(sub(1n, x[i]!), sub(1n, b[i]!)));
    result = mul(result, term);
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
// Section 1: Copy Constraints
// ============================================================================

console.log("--- Section 1: Copy Constraints ---\n");

/**
 * COPY CONSTRAINTS IN CIRCUITS
 *
 * A circuit has wires connecting gates. When we represent the circuit
 * as polynomials, we have:
 *
 *   - Wire polynomials: w_1, w_2, w_3 (left input, right input, output)
 *   - Each w_j: B_mu -> F gives the value of wire j at each gate
 *
 * COPY CONSTRAINT:
 * When the output of gate i connects to the input of gate j:
 *   w_3(i) = w_1(j)  or  w_3(i) = w_2(j)
 *
 * PERMUTATION ENCODING:
 * We encode all copy constraints as a single permutation sigma.
 *
 * sigma acts on the set of all wire positions:
 *   {(b, j) : b in B_mu, j in {1, 2, 3}}
 *
 * If positions (b, j) and (b', j') must have equal values,
 * sigma puts them in the same cycle.
 *
 * EXAMPLE:
 * If w_3(00) = w_1(01) = w_2(10), then sigma has cycle:
 *   (00, 3) -> (01, 1) -> (10, 2) -> (00, 3)
 *
 * All positions in a cycle must have the same wire value.
 */

console.log("Copy constraints:");
console.log("  When output of gate i -> input of gate j:");
console.log("  w_3(i) = w_1(j) or w_2(j)\n");

console.log("Permutation encoding:");
console.log("  sigma acts on {(gate, wire_index)}");
console.log("  Positions in same cycle must have equal values\n");

console.log("Example cycle:");
console.log("  (00, 3) -> (01, 1) -> (10, 2) -> (00, 3)");
console.log("  Means: w_3(00) = w_1(01) = w_2(10)\n");

// ============================================================================
// Section 2: The Permutation Argument
// ============================================================================

console.log("--- Section 2: Permutation Argument ---\n");

/**
 * PERMUTATION ARGUMENT
 *
 * To verify copy constraints, we check that wire values respect sigma.
 *
 * CLAIM: For all positions (b, j):
 *   w_j(b) = w_{sigma_j(b)}(sigma_b(b))
 *
 * where sigma maps (b, j) to (sigma_b(b), sigma_j(b)).
 *
 * GRAND PRODUCT ARGUMENT:
 *
 * Create tagged values:
 *   p(b, j) = w_j(b) + gamma * ID(b, j)
 *   q(b, j) = w_j(b) + gamma * sigma(b, j)
 *
 * where ID(b, j) is a unique identifier for position (b, j).
 *
 * If copy constraints are satisfied:
 *   {p(b, j)} = {q(b, j)} as multisets
 *
 * Because for each position, the value is the same, just tagged differently.
 *
 * Using ProductCheck:
 *   prod_{b,j} (beta + p(b, j)) = prod_{b,j} (beta + q(b, j))
 */

console.log("Permutation argument:\n");

console.log("Tagged values:");
console.log("  p(b, j) = w_j(b) + gamma * ID(b, j)");
console.log("  q(b, j) = w_j(b) + gamma * sigma(b, j)\n");

console.log("If copy constraints hold:");
console.log("  {p} = {q} as multisets");
console.log("  => prod (beta + p) = prod (beta + q)\n");

// ============================================================================
// Section 3: Example Circuit Wiring
// ============================================================================

console.log("--- Section 3: Example Circuit Wiring ---\n");

/**
 * EXAMPLE CIRCUIT
 *
 * 4 gates (mu = 2):
 *
 * Gate 0: w_1=5, w_2=3, w_3=8     (5 + 3 = 8)
 * Gate 1: w_1=8, w_2=5, w_3=40    (8 * 5 = 40)
 * Gate 2: w_1=5, w_2=0, w_3=0     (constant 5)
 * Gate 3: w_1=0, w_2=0, w_3=40    (output 40)
 *
 * Copy constraints:
 * - w_3(0) = w_1(1)    (output of add -> first input of mul)
 * - w_1(0) = w_2(1)    (both use 5)
 * - w_1(0) = w_1(2)    (both use 5)
 * - w_3(1) = w_3(3)    (multiplication result is final output)
 *
 * Permutation cycles:
 * - Cycle 1: (0,3) -> (1,1) -> (0,3)  [value 8]
 * - Cycle 2: (0,1) -> (1,2) -> (2,1) -> (0,1)  [value 5]
 * - Cycle 3: (1,3) -> (3,3) -> (1,3)  [value 40]
 * - All other positions are fixed points (self-cycles)
 */

const mu = 2;
const n = 2 ** mu;  // 4 gates
const numWires = 3;  // w_1, w_2, w_3
const cube = booleanHypercube(mu);

// Wire polynomials
const wires: HypercubeFunction[] = [
  { numVars: mu, values: [5n, 8n, 5n, 0n] },  // w_1
  { numVars: mu, values: [3n, 5n, 0n, 0n] },  // w_2
  { numVars: mu, values: [8n, 40n, 0n, 40n] } // w_3
];

console.log("Wire values:");
for (let j = 0; j < numWires; j++) {
  console.log(`  w_${j + 1}: [${wires[j]!.values.join(', ')}]`);
}
console.log();

// Define permutation
// Position encoding: pos = gate * 3 + wire_index
// For gate b (0-3) and wire j (0-2): pos = b * 3 + j
const sigma: number[] = new Array(n * numWires);

// Initialize as identity
for (let i = 0; i < n * numWires; i++) {
  sigma[i] = i;
}

// Cycle 1: (0,3) -> (1,1) -> (0,3)
// In our encoding: (0,2) -> (1,0) -> (0,2)
// pos(0,2) = 0*3+2 = 2, pos(1,0) = 1*3+0 = 3
sigma[2] = 3;  // (0, w_3) -> (1, w_1)
sigma[3] = 2;  // (1, w_1) -> (0, w_3)

// Cycle 2: (0,1) -> (1,2) -> (2,1) -> (0,1)
// pos(0,0) = 0, pos(1,1) = 4, pos(2,0) = 6
sigma[0] = 4;  // (0, w_1) -> (1, w_2)
sigma[4] = 6;  // (1, w_2) -> (2, w_1)
sigma[6] = 0;  // (2, w_1) -> (0, w_1)

// Cycle 3: (1,3) -> (3,3) -> (1,3)
// pos(1,2) = 5, pos(3,2) = 11
sigma[5] = 11;  // (1, w_3) -> (3, w_3)
sigma[11] = 5;  // (3, w_3) -> (1, w_3)

console.log("Copy constraint cycles:");
console.log("  Cycle 1 (value 8):  w_3(00) <-> w_1(01)");
console.log("  Cycle 2 (value 5):  w_1(00) -> w_2(01) -> w_1(10) -> w_1(00)");
console.log("  Cycle 3 (value 40): w_3(01) <-> w_3(11)\n");

// ============================================================================
// Section 4: Verifying Copy Constraints
// ============================================================================

console.log("--- Section 4: Verifying Copy Constraints ---\n");

// Helper to get wire value at position
function getWireValue(pos: number): bigint {
  const gate = Math.floor(pos / numWires);
  const wire = pos % numWires;
  return wires[wire]!.values[gate]!;
}

// Verify each cycle
console.log("Verifying cycles:");

// Check that sigma(pos) has same wire value as pos
let wiringValid = true;
for (let pos = 0; pos < n * numWires; pos++) {
  const sigmaPos = sigma[pos]!;
  const val1 = getWireValue(pos);
  const val2 = getWireValue(sigmaPos);

  if (val1 !== val2) {
    const gate1 = Math.floor(pos / numWires);
    const wire1 = pos % numWires;
    const gate2 = Math.floor(sigmaPos / numWires);
    const wire2 = sigmaPos % numWires;
    console.log(`  FAIL: w_${wire1 + 1}(${cube[gate1]!.join('')})=${val1} != w_${wire2 + 1}(${cube[gate2]!.join('')})=${val2}`);
    wiringValid = false;
  }
}

if (wiringValid) {
  console.log("  All copy constraints satisfied!\n");
}

// ============================================================================
// Section 5: ProductCheck for Wiring
// ============================================================================

console.log("--- Section 5: ProductCheck for Wiring ---\n");

/**
 * PRODUCTCHECK FOR WIRING
 *
 * Random challenges: beta, gamma
 *
 * For each position pos = (gate, wire):
 *   p[pos] = w[pos] + gamma * ID[pos]
 *   q[pos] = w[pos] + gamma * sigma[pos]
 *
 * where ID[pos] = pos (unique identifier).
 *
 * ProductCheck:
 *   prod_pos (beta + p[pos]) = prod_pos (beta + q[pos])
 */

const beta = 11n;
const gamma = 17n;
console.log(`Random challenges: beta=${beta}, gamma=${gamma}\n`);

let prodP = 1n;
let prodQ = 1n;

console.log("Computing products:");
for (let pos = 0; pos < n * numWires; pos++) {
  const gate = Math.floor(pos / numWires);
  const wire = pos % numWires;
  const w = wires[wire]!.values[gate]!;

  const id = BigInt(pos);
  const sigmaId = BigInt(sigma[pos]!);

  const pVal = add(add(beta, w), mul(gamma, id));
  const qVal = add(add(beta, w), mul(gamma, sigmaId));

  prodP = mul(prodP, pVal);
  prodQ = mul(prodQ, qVal);
}

console.log(`  prod (beta + w + gamma*ID) = ${prodP}`);
console.log(`  prod (beta + w + gamma*sigma) = ${prodQ}`);
console.log(`  Equal? ${prodP === prodQ ? 'YES' : 'NO'}\n`);

// ============================================================================
// Section 6: Detecting Wiring Violation
// ============================================================================

console.log("--- Section 6: Detecting Wiring Violation ---\n");

// Create bad wires where copy constraint is violated
const wiresBad: HypercubeFunction[] = [
  { numVars: mu, values: [5n, 9n, 5n, 0n] },  // w_1(1) = 9 instead of 8
  { numVars: mu, values: [3n, 5n, 0n, 0n] },
  { numVars: mu, values: [8n, 40n, 0n, 40n] }
];

console.log("Bad wires: w_1(01) = 9 instead of 8");
console.log("  Violates: w_3(00) = 8 != w_1(01) = 9\n");

function getWireValueBad(pos: number): bigint {
  const gate = Math.floor(pos / numWires);
  const wire = pos % numWires;
  return wiresBad[wire]!.values[gate]!;
}

let prodPBad = 1n;
let prodQBad = 1n;

for (let pos = 0; pos < n * numWires; pos++) {
  const gate = Math.floor(pos / numWires);
  const wire = pos % numWires;
  const w = wiresBad[wire]!.values[gate]!;

  const id = BigInt(pos);
  const sigmaId = BigInt(sigma[pos]!);

  const pVal = add(add(beta, w), mul(gamma, id));
  const qVal = add(add(beta, w), mul(gamma, sigmaId));

  prodPBad = mul(prodPBad, pVal);
  prodQBad = mul(prodQBad, qVal);
}

console.log("ProductCheck with bad wires:");
console.log(`  prod P = ${prodPBad}`);
console.log(`  prod Q = ${prodQBad}`);
console.log(`  Equal? ${prodPBad === prodQBad ? 'YES (false positive!)' : 'NO (caught!)'}\n`);

// ============================================================================
// Section 7: Full Wiring Protocol
// ============================================================================

console.log("--- Section 7: Full Wiring Protocol ---\n");

/**
 * FULL WIRING PROTOCOL
 *
 * SETUP:
 * - Permutation sigma encoded as polynomial
 * - ID polynomials for each wire column
 *
 * PROTOCOL:
 * 1. V -> P: beta, gamma (random challenges)
 *
 * 2. P computes tagged values:
 *    f(b) = prod_j (beta + w_j(b) + gamma * ID_j(b))
 *    g(b) = prod_j (beta + w_j(b) + gamma * sigma_j(b))
 *
 * 3. P runs ProductCheck on (f, g):
 *    - Computes running product h
 *    - ZeroCheck on constraint
 *
 * 4. V verifies:
 *    - ProductCheck verification
 *    - Polynomial openings
 *
 * OPTIMIZATION:
 *
 * Instead of separate products, combine all wires:
 *   f(b) = prod_{j=1}^k (beta + w_j(b) + gamma * ID_j(b))
 *
 * This gives a single polynomial with degree k in each variable.
 */

console.log("Wiring protocol:\n");

console.log("  1. V sends beta, gamma");
console.log("  2. P computes tagged products:");
console.log("     f(b) = prod_j (beta + w_j(b) + gamma * ID_j(b))");
console.log("     g(b) = prod_j (beta + w_j(b) + gamma * sigma_j(b))");
console.log("  3. P runs ProductCheck on (f, g)");
console.log("  4. V verifies ProductCheck\n");

// Compute f and g polynomials
const fValues: bigint[] = [];
const gValues: bigint[] = [];

for (let gate = 0; gate < n; gate++) {
  let fProd = 1n;
  let gProd = 1n;

  for (let wire = 0; wire < numWires; wire++) {
    const pos = gate * numWires + wire;
    const w = wires[wire]!.values[gate]!;

    const id = BigInt(pos);
    const sigmaId = BigInt(sigma[pos]!);

    fProd = mul(fProd, add(add(beta, w), mul(gamma, id)));
    gProd = mul(gProd, add(add(beta, w), mul(gamma, sigmaId)));
  }

  fValues.push(fProd);
  gValues.push(gProd);
}

console.log("Combined polynomials f and g:");
for (let gate = 0; gate < n; gate++) {
  console.log(`  Gate ${cube[gate]!.join('')}: f=${fValues[gate]!}, g=${gValues[gate]!}`);
}

const finalFProd = fValues.reduce((a, b) => mul(a, b), 1n);
const finalGProd = gValues.reduce((a, b) => mul(a, b), 1n);
console.log(`\nprod f = ${finalFProd}`);
console.log(`prod g = ${finalGProd}`);
console.log(`Equal? ${finalFProd === finalGProd ? 'YES' : 'NO'}\n`);

// ============================================================================
// Section 8: Implementation Stubs
// ============================================================================

console.log("--- Section 8: Implementation Stubs ---\n");

/**
 * Full wiring identity check using ProductCheck.
 *
 * @throws NotImplementedError - Full implementation pending
 */
export function wiringIdentityCheck(
  wires: HypercubeFunction[],
  sigma: number[],
  beta: bigint,
  gamma: bigint,
  getChallenge: () => bigint
): { accepted: boolean } {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: wiringIdentityCheck - ' +
    'Full wiring check with ProductCheck'
  );
}

/**
 * Compute combined wiring polynomials f and g.
 *
 * @throws NotImplementedError - Combination logic pending
 */
export function computeWiringPolynomials(
  wires: HypercubeFunction[],
  sigma: number[],
  beta: bigint,
  gamma: bigint
): { f: HypercubeFunction; g: HypercubeFunction } {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: computeWiringPolynomials - ' +
    'Compute tagged wiring polynomials'
  );
}

console.log("Stub functions:");
console.log("  - wiringIdentityCheck: full wiring verification");
console.log("  - computeWiringPolynomials: tagged polynomial computation\n");

// ============================================================================
// Summary
// ============================================================================

console.log("=== Summary ===\n");

console.log("Wiring Identity Check:");
console.log("  Verifies copy constraints between wire positions");
console.log("  Encoded as permutation sigma\n");

console.log("Method:");
console.log("  1. Tag values with positions: w + gamma * pos");
console.log("  2. ProductCheck: prod(tagged) = prod(sigma-tagged)");
console.log("  3. Reduces to ZeroCheck on constraint polynomial\n");

console.log("Efficiency:");
console.log("  - Combines all wires into one ProductCheck");
console.log("  - Communication: O(mu) field elements");
console.log("  - Prover: O(2^mu) with optimization\n");

console.log("Next: Full HyperPlonk demonstration");

// Exports
export { type HypercubeFunction };
