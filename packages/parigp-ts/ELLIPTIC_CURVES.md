# Elliptic Curves Design Document for parigp-ts

This document describes the elliptic curve functionality to be ported from PARI/GP to TypeScript, based on analysis of the PARI/GP source code in `reference/pari/src/basemath/elliptic.c`, `FpE.c`, `Fle.c`, and related files.

## Table of Contents

1. [PARI/GP Functions to Port](#1-parigp-functions-to-port)
2. [Data Structures](#2-data-structures)
3. [Dependencies](#3-dependencies)
4. [Implementation Order](#4-implementation-order)
5. [Test Vectors](#5-test-vectors)
6. [References](#6-references)

---

## 1. PARI/GP Functions to Port

### 1.1 Core Curve Initialization

#### `ellinit(x, D, prec)` - Initialize Elliptic Curve
**Source:** `elliptic.c:887-893`

Creates an elliptic curve structure from coefficients or j-invariant.

**Input formats:**
- `[a1, a2, a3, a4, a6]` - General Weierstrass coefficients
- `[a4, a6]` - Short Weierstrass form (y^2 = x^3 + a4*x + a6)
- `[j]` - j-invariant (curve from j-invariant)

**Domain `D` types:**
- `t_INT` p: Prime, creates curve over Fp
- `t_FFELT`: Finite field element, creates curve over Fq
- `t_PADIC`: p-adic field
- `NULL`: Curve over Q or generic ring

**Returns:** 16-component vector:
```
[a1, a2, a3, a4, a6, b2, b4, b6, b8, c4, c6, disc, j, type, domain_data, dynamic_data]
```

Where the derived quantities are:
- `b2 = a1^2 + 4*a2`
- `b4 = a1*a3 + 2*a4`
- `b6 = a3^2 + 4*a6`
- `b8 = a1^2*a6 + 4*a2*a6 + a2*a3^2 - a4*(a4 + a1*a3)`
- `c4 = b2^2 - 24*b4`
- `c6 = 36*b2*b4 - b2^3 - 216*b6`
- `disc = (c4^3 - c6^2) / 1728`
- `j = c4^3 / disc`

**TypeScript signature:**
```typescript
function ellinit(
  x: [bigint, bigint, bigint, bigint, bigint] | [bigint, bigint] | [bigint],
  D?: bigint | FiniteField,
  prec?: number
): EllipticCurve;
```

---

### 1.2 Point Operations

#### `ellisoncurve(E, P)` - Check if Point is on Curve
**Source:** `elliptic.c:2003-2037`

Tests if point P lies on curve E by checking:
- For point at infinity: always true
- Otherwise: `y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6`

**TypeScript signature:**
```typescript
function ellisoncurve(E: EllipticCurve, P: EllipticPoint): boolean;
```

---

#### `elladd(E, P, Q)` - Point Addition
**Source:** `elliptic.c:2061-2094`

Adds two points on an elliptic curve using the standard chord-and-tangent method.

**Algorithm:**
1. If P = O (infinity), return Q
2. If Q = O, return P
3. If P.x == Q.x:
   - If P.y != Q.y, return O (point at infinity)
   - Otherwise use tangent (doubling formula)
4. Otherwise compute slope `s = (Q.y - P.y) / (Q.x - P.x)`
5. `R.x = s*(s + a1) - P.x - Q.x - a2`
6. `R.y = -(P.y + a1*R.x + a3 + s*(R.x - P.x))`

**TypeScript signature:**
```typescript
function elladd(E: EllipticCurve, P: EllipticPoint, Q: EllipticPoint): EllipticPoint;
```

---

#### `ellneg(E, P)` - Point Negation
**Source:** `elliptic.c:2113-2126`

Computes the negation of a point.

**Formula:**
- `-P = [P.x, -(P.y + a1*P.x + a3)]`

**TypeScript signature:**
```typescript
function ellneg(E: EllipticCurve, P: EllipticPoint): EllipticPoint;
```

---

#### `ellsub(E, P, Q)` - Point Subtraction
**Source:** `elliptic.c:2129-2136`

Computes P - Q = P + (-Q).

**TypeScript signature:**
```typescript
function ellsub(E: EllipticCurve, P: EllipticPoint, Q: EllipticPoint): EllipticPoint;
```

---

#### `ellmul(E, P, n)` - Scalar Multiplication
**Source:** `elliptic.c:2306-2316, 2283-2316`

Computes [n]P using double-and-add algorithm.

For finite fields, uses Jacobian coordinates internally for efficiency:
- `FpE_mul` in FpE.c (line 360-365)
- Uses `gen_pow_i` with Jacobian `FpJ_dbl` and `FpJ_add`

**TypeScript signature:**
```typescript
function ellmul(E: EllipticCurve, P: EllipticPoint, n: bigint): EllipticPoint;
```

---

### 1.3 Group Structure

#### `ellcard(E, p?)` - Cardinality (Number of Points)
**Source:** `elliptic.c:6332-6388`

Computes the number of points on E(Fq), including the point at infinity.

**Algorithms used:**
- Small primes: Naive counting (`Fl_elltrace_naive` in FpE.c:812-832)
- Medium primes: Baby-step giant-step (Shanks/Mestre) (`Fp_ellcard_Shanks` in FpE.c:921-1135)
- Large primes: Schoof-Elkies-Atkin (SEA) algorithm (`Fp_ellcard_SEA` in ellsea.c)

**TypeScript signature:**
```typescript
function ellcard(E: EllipticCurve, p?: bigint): bigint;
```

---

#### `ellgroup(E, p?)` - Group Structure
**Source:** `elliptic.c:6462-6495`

Returns the structure of E(Fq) as a product of cyclic groups [d1] or [d1, d2] where d2 | d1.

Uses `Fp_ellgroup` from FpE.c:
```c
GEN Fp_ellgroup(GEN a4, GEN a6, GEN N, GEN p, GEN *pt_m);
```

**TypeScript signature:**
```typescript
function ellgroup(E: EllipticCurve, p?: bigint): bigint[];
```

---

#### `ellgenerators(E)` - Find Generators
**Source:** `elliptic.c:6444-6459`

Returns generators for the group E(Fq).

**TypeScript signature:**
```typescript
function ellgenerators(E: EllipticCurve): EllipticPoint[];
```

---

#### `ellorder(E, P, o?)` - Order of a Point
**Source:** `elliptic.c` (via `FpE_order` in FpE.c:406-423)

Computes the order of point P in E(Fq).

**Algorithm:**
Uses generic group order algorithm with factorization of the known group order.

**TypeScript signature:**
```typescript
function ellorder(E: EllipticCurve, P: EllipticPoint, o?: bigint | [bigint, number[][]]): bigint;
```

---

### 1.4 Other Essential Functions

#### `ellordinate(E, x)` - Find y-coordinates
**Source:** `elliptic.c:2140-2211`

Given x-coordinate, finds all valid y-coordinates (0, 1, or 2 values).

Solves: `y^2 + (a1*x + a3)*y - (x^3 + a2*x^2 + a4*x + a6) = 0`

**TypeScript signature:**
```typescript
function ellordinate(E: EllipticCurve, x: FieldElement): FieldElement[];
```

---

#### `ellrandom(E)` - Random Point
**Source:** `elliptic.c:2214-2229`, `FpE.c:369-385`

Returns a random point on the curve (excluding point at infinity).

**Algorithm:**
1. Pick random x
2. Compute RHS = x^3 + a4*x + a6
3. Check if RHS is a quadratic residue
4. If yes, compute y = sqrt(RHS)
5. Otherwise repeat

**TypeScript signature:**
```typescript
function ellrandom(E: EllipticCurve): EllipticPoint;
```

---

#### `elllog(E, P, G, o?)` - Discrete Logarithm
**Source:** `FpE.c:425-443`

Computes n such that P = [n]G (if it exists).

Uses Pohlig-Hellman algorithm via `gen_PH_log` in bb_group.c.

**TypeScript signature:**
```typescript
function elllog(E: EllipticCurve, P: EllipticPoint, G: EllipticPoint, o?: bigint): bigint;
```

---

#### `elltatepairing(E, P, Q, m)` - Tate Pairing
**Source:** `FpE.c:609-621`

Computes the Tate pairing of P and Q for m-torsion points.

**TypeScript signature:**
```typescript
function elltatepairing(E: EllipticCurve, P: EllipticPoint, Q: EllipticPoint, m: bigint): FieldElement;
```

---

#### `ellweilpairing(E, P, Q, m)` - Weil Pairing
**Source:** `FpE.c:588-606`

Computes the Weil pairing of P and Q for m-torsion points.

Uses Miller's algorithm (`FpE_Miller` in FpE.c:574-586).

**TypeScript signature:**
```typescript
function ellweilpairing(E: EllipticCurve, P: EllipticPoint, Q: EllipticPoint, m: bigint): FieldElement;
```

---

### 1.5 Curve Transformations

#### `ellchangecurve(E, v)` - Change of Variables
**Source:** `elliptic.c:1468-1498`

Applies the change of variables [u, r, s, t] to the curve.

The transformation (x, y) -> (x', y') is:
- `x = u^2 * x' + r`
- `y = u^3 * y' + s*u^2*x' + t`

**TypeScript signature:**
```typescript
function ellchangecurve(E: EllipticCurve, v: [bigint, bigint, bigint, bigint]): EllipticCurve;
```

---

#### `ellchangepoint(P, v)` - Transform Point
**Source:** `FpE.c:190-208`

Applies change of variables to a point.

**TypeScript signature:**
```typescript
function ellchangepoint(P: EllipticPoint, v: ChangeOfVariables): EllipticPoint;
```

---

## 2. Data Structures

### 2.1 Elliptic Curve Structure

PARI represents an elliptic curve as a 16-element vector (from `elliptic.c:443-495`):

```typescript
interface EllipticCurve {
  // Weierstrass coefficients (y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6)
  a1: FieldElement;  // index 1
  a2: FieldElement;  // index 2
  a3: FieldElement;  // index 3
  a4: FieldElement;  // index 4
  a6: FieldElement;  // index 5

  // Derived quantities
  b2: FieldElement;  // index 6: a1^2 + 4*a2
  b4: FieldElement;  // index 7: a1*a3 + 2*a4
  b6: FieldElement;  // index 8: a3^2 + 4*a6
  b8: FieldElement;  // index 9: a1^2*a6 + 4*a2*a6 + a2*a3^2 - a4*(a4 + a1*a3)
  c4: FieldElement;  // index 10: b2^2 - 24*b4
  c6: FieldElement;  // index 11: 36*b2*b4 - b2^3 - 216*b6
  disc: FieldElement; // index 12: discriminant
  j: FieldElement;   // index 13: j-invariant

  // Type information
  type: EllCurveType; // index 14: t_ELL_Rg, t_ELL_Q, t_ELL_Qp, t_ELL_Fp, t_ELL_Fq, t_ELL_NF

  // Domain-specific data
  domainData: DomainData; // index 15

  // Dynamic (cached) data
  dynamicData: DynamicData; // index 16
}

// Curve types from paridecl.h:3340
enum EllCurveType {
  t_ELL_Rg = 0,  // Generic ring
  t_ELL_Q = 1,   // Rationals
  t_ELL_Qp = 2,  // p-adics
  t_ELL_Fp = 3,  // Prime field
  t_ELL_Fq = 4,  // Extension field
  t_ELL_NF = 5,  // Number field
}
```

### 2.2 Point Representation

From `pariinl.h:2753-2754`:

```typescript
// Point at infinity represented as single-element vector [0]
// Finite point represented as [x, y]

type EllipticPoint =
  | { isInfinity: true }
  | { isInfinity: false; x: FieldElement; y: FieldElement };

// Helper functions
function ellinf(): EllipticPoint {
  return { isInfinity: true };
}

function ell_is_inf(P: EllipticPoint): boolean {
  return P.isInfinity;
}
```

### 2.3 Jacobian Coordinates (Internal)

From `FpE.c:26-31`:

```typescript
// Jacobian coordinates: (X : Y : Z) represents affine (X/Z^2, Y/Z^3)
// Point at infinity: Z = 0, (1 : 1 : 0)
interface JacobianPoint {
  X: bigint;
  Y: bigint;
  Z: bigint;
}

// Conversion functions
function FpE_to_FpJ(P: AffinePoint): JacobianPoint;
function FpJ_to_FpE(P: JacobianPoint, p: bigint): AffinePoint;
```

### 2.4 Short Weierstrass Form

For finite fields, PARI internally converts to short Weierstrass form (y^2 = x^3 + a4*x + a6).

From `elliptic.c:63-68`:
```typescript
// Convert from general form using c4, c6
// a4' = -27 * c4
// a6' = -54 * c6
function ell_to_a4a6(E: EllipticCurve, p: bigint): [bigint, bigint];
```

### 2.5 Change of Variables

```typescript
// [u, r, s, t] where u != 0
// Transformation: x = u^2*x' + r, y = u^3*y' + s*u^2*x' + t
type ChangeOfVariables = [bigint, bigint, bigint, bigint];
```

---

## 3. Dependencies

### 3.1 Required Before Implementing EC

1. **Finite Field Arithmetic** (Critical)
   - `Fp_add`, `Fp_sub`, `Fp_mul`, `Fp_div`, `Fp_neg`
   - `Fp_sqr`, `Fp_sqrt`, `Fp_inv`
   - `Fp_pow` (modular exponentiation)
   - Quadratic residue testing (`kronecker` symbol)

2. **Big Integer Arithmetic**
   - Basic operations (+, -, *, /, mod)
   - GCD, extended GCD
   - Integer square root

3. **Extension Field Arithmetic** (for Fq support)
   - Polynomial arithmetic over Fp
   - `FpXQ_*` operations (quotient ring Fp[x]/(T))

### 3.2 Algorithms Needed

1. **Square Root in Finite Fields**
   - Tonelli-Shanks algorithm
   - For p = 3 mod 4: sqrt(a) = a^((p+1)/4)

2. **Point Counting** (for `ellcard`)
   - Baby-step giant-step
   - Optional: SEA algorithm for large primes

3. **Generic Group Operations** (from bb_group.c)
   - `gen_pow` - Generic exponentiation
   - `gen_order` - Order finding
   - `gen_PH_log` - Pohlig-Hellman discrete log

4. **Miller's Algorithm** (for pairings)
   - Line evaluation
   - Vertical line evaluation

---

## 4. Implementation Order

### Phase 1: Foundation
1. **Finite field arithmetic** (Fp)
   - Addition, subtraction, multiplication, division
   - Square root (Tonelli-Shanks)
   - Quadratic residue testing

2. **Basic curve structure**
   - `ellinit` for Fp curves (short Weierstrass)
   - Coefficient computation (b2, b4, b6, b8, c4, c6, disc, j)

### Phase 2: Core Point Operations
3. **Point arithmetic** (affine coordinates)
   - `ellisoncurve`
   - `elladd`
   - `ellneg`
   - `ellsub`

4. **Jacobian coordinates** (for efficiency)
   - `FpJ_dbl`, `FpJ_add`
   - Conversion functions

5. **Scalar multiplication**
   - `ellmul` using double-and-add
   - NAF (Non-Adjacent Form) optimization

### Phase 3: Group Operations
6. **Random point generation**
   - `ellrandom`

7. **Point order**
   - `ellorder` using factorization

8. **Point counting**
   - Naive method for small p
   - Baby-step giant-step for medium p

9. **Group structure**
   - `ellgroup`
   - `ellgenerators`

### Phase 4: Advanced Features
10. **Discrete logarithm**
    - `elllog` (Pohlig-Hellman)

11. **Pairings**
    - Miller's algorithm
    - `elltatepairing`
    - `ellweilpairing`

12. **Extension fields**
    - FpXQ arithmetic
    - Curves over Fq

### Phase 5: Optimization
13. **Montgomery ladder** (constant-time)
14. **Precomputation tables** (for fixed-base multiplication)
15. **SEA algorithm** (for very large primes)

---

## 5. Test Vectors

### 5.1 Basic Operations over Fp

From PARI test files (`reference/pari/src/test/32/ell`):

```typescript
// Curve: y^2 = x^3 + a4*x + a6 over Fp
// Test: ellisoncurve with point (0, 0) on y^2 = x^3
const testCurve1 = { a4: 0n, a6: 0n, p: 1009n };
// Point (0, 0) is on curve: 0^2 = 0^3 + 0*0 + 0

// Scalar multiplication: [0] * P = O (point at infinity)
// ellmul(E, [0,0], 0) -> [0] (infinity)

// Point negation on curve y^2 = x^3 over Fp
// ellneg(E, [0,0]) -> [0, 0] (since a1=a3=0)
```

### 5.2 Group Structure

From test/32/ell (lines 167-200):
```typescript
// ellgroup tests with various primes
const testVectors = [
  { p: 1031n, group: [504n, 2n] },      // Non-cyclic
  { p: 2053n, group: [1008n, 2n] },     // Non-cyclic
  { p: 4099n, group: [4196n] },          // Cyclic
  { p: 8209n, group: [8291n] },          // Cyclic
  { p: 16411n, group: [8280n, 2n] },    // Non-cyclic
];
```

### 5.3 Well-Known Test Curves

**secp256k1** (Bitcoin):
```typescript
const secp256k1 = {
  p: 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2Fn,
  a4: 0n,
  a6: 7n,
  // Generator
  Gx: 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798n,
  Gy: 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8n,
  // Order
  n: 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141n,
};
```

**P-256 (NIST)**:
```typescript
const P256 = {
  p: 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFFn,
  a4: -3n, // mod p
  a6: 0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604Bn,
};
```

### 5.4 Simple Arithmetic Tests

```typescript
// Small prime for easy verification
const p = 23n;
// Curve: y^2 = x^3 + x + 1 over F_23
const E = { a4: 1n, a6: 1n, p: 23n };

// Points on curve (verify y^2 = x^3 + x + 1 mod 23):
// (0, 1): 1 = 0 + 0 + 1 = 1
// (0, 22): 484 mod 23 = 1
// (1, 7): 49 mod 23 = 3 = 1 + 1 + 1
// (1, 16): 256 mod 23 = 3

// Point addition test:
// P = (0, 1), Q = (1, 7)
// slope = (7 - 1) / (1 - 0) = 6 mod 23
// x_R = 6^2 - 0 - 1 = 35 mod 23 = 12
// y_R = 6*(0 - 12) - 1 = -73 mod 23 = -73 + 4*23 = 19
// P + Q = (12, 19)

// Point doubling test:
// [2]P where P = (0, 1)
// slope = (3*0^2 + 1) / (2*1) = 1/2 = 12 mod 23 (since 2*12 = 24 = 1 mod 23)
// x = 12^2 - 2*0 = 144 mod 23 = 6
// y = 12*(0 - 6) - 1 = -73 mod 23 = 19
// [2]P = (6, 19)
```

---

## 6. References

### 6.1 PARI/GP Source Files

| File | Description |
|------|-------------|
| `elliptic.c` | Main elliptic curve functions |
| `FpE.c` | Curves over prime fields Fp |
| `Fle.c` | Small prime (ulong) optimizations |
| `FlxqE.c` | Curves over Fq with ulong characteristic |
| `F2xqE.c` | Curves over F_{2^n} |
| `ellsea.c` | Schoof-Elkies-Atkin algorithm |
| `bb_group.c` | Generic group algorithms |
| `parigen.h` | Type definitions |
| `pariinl.h` | Inline accessor functions |
| `paridecl.h` | Function declarations |

### 6.2 Key Constants and Types

From `paridecl.h:3340`:
```c
enum { t_ELL_Rg = 0, t_ELL_Q, t_ELL_Qp, t_ELL_Fp, t_ELL_Fq, t_ELL_NF };
```

From `paripriv.h:187`:
```c
// Cached data indices for Fp/Fq curves
enum { FF_CARD = 1, FF_GROUP, FF_GROUPGEN, FF_O };
```

### 6.3 External References

1. PARI/GP User's Guide, Chapter on Elliptic Curves
2. J.H. Silverman, "The Arithmetic of Elliptic Curves" (GTM 106)
3. Hyperelliptic.org EFD (Explicit-Formulas Database): https://hyperelliptic.org/EFD/
4. IEEE P1363: Standard for Public-Key Cryptography

### 6.4 Efficient Formula References

**Jacobian Doubling (dbl-2007-bl):**
- Cost: 1M + 8S + 1*a + 10add + 1*8 + 2*2 + 1*3
- Source: http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-2007-bl

**Jacobian Addition (add-2007-bl):**
- Cost: 11M + 5S + 9add + 4*2
- Source: http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-2007-bl

---

## Appendix A: Formula Summary

### A.1 General Weierstrass Form

Curve equation: `y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6`

**Derived quantities:**
```
b2 = a1^2 + 4*a2
b4 = a1*a3 + 2*a4
b6 = a3^2 + 4*a6
b8 = a1^2*a6 + 4*a2*a6 + a2*a3^2 - a4*(a4 + a1*a3)

c4 = b2^2 - 24*b4
c6 = 36*b2*b4 - b2^3 - 216*b6

disc = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
     = (c4^3 - c6^2) / 1728

j = c4^3 / disc
```

### A.2 Point Addition (Affine)

For points P = (x1, y1) and Q = (x2, y2):

**Case P != Q:**
```
lambda = (y2 - y1) / (x2 - x1)
x3 = lambda^2 + a1*lambda - a2 - x1 - x2
y3 = -(y1 + a1*x3 + a3 + lambda*(x3 - x1))
```

**Case P = Q (doubling):**
```
lambda = (3*x1^2 + 2*a2*x1 + a4 - a1*y1) / (2*y1 + a1*x1 + a3)
x3 = lambda^2 + a1*lambda - a2 - 2*x1
y3 = -(y1 + a1*x3 + a3 + lambda*(x3 - x1))
```

### A.3 Short Weierstrass (char > 3)

Curve equation: `y^2 = x^3 + a4*x + a6`

**Point addition (P != Q):**
```
lambda = (y2 - y1) / (x2 - x1)
x3 = lambda^2 - x1 - x2
y3 = lambda*(x1 - x3) - y1
```

**Point doubling:**
```
lambda = (3*x1^2 + a4) / (2*y1)
x3 = lambda^2 - 2*x1
y3 = lambda*(x1 - x3) - y1
```

### A.4 Jacobian Coordinates

Representation: (X : Y : Z) corresponds to affine (X/Z^2, Y/Z^3)

Point at infinity: (1 : 1 : 0)

**Doubling (2007-bl formula):**
```
XX = X1^2
YY = Y1^2
YYYY = YY^2
ZZ = Z1^2
S = 2*((X1+YY)^2 - XX - YYYY)
M = 3*XX + a4*ZZ^2
T = M^2 - 2*S
X3 = T
Y3 = M*(S-T) - 8*YYYY
Z3 = (Y1+Z1)^2 - YY - ZZ
```

**Addition (2007-bl formula):**
```
Z1Z1 = Z1^2
Z2Z2 = Z2^2
U1 = X1*Z2Z2
U2 = X2*Z1Z1
S1 = Y1*Z2*Z2Z2
S2 = Y2*Z1*Z1Z1
H = U2 - U1
r = 2*(S2 - S1)
I = (2*H)^2
J = H*I
V = U1*I
X3 = r^2 - J - 2*V
Y3 = r*(V - X3) - 2*S1*J
Z3 = ((Z1+Z2)^2 - Z1Z1 - Z2Z2)*H
```
