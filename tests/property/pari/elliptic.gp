\\ Property tests for elliptic curve operations
\\ This script outputs results in a parseable format for comparison with parigp-ts
\\
\\ Usage: gp -q < elliptic.gp
\\ Or with seed: gp -q --seed 12345 < elliptic.gp
\\
\\ Output format:
\\ # <operation>
\\ RESULT: <value>
\\
\\ Note: We use small primes to keep tests fast and deterministic

\\ Set random seed for reproducibility
setrand(12345);

\\ Helper function to format points
format_point(P) = {
  if(P == [0],
    return("INFINITY"),
    return(Str("[", P[1], ", ", P[2], "]"))
  )
};

\\ Helper function to format curve info
format_curve(E) = {
  return(Str("a1=", E.a1, ", a2=", E.a2, ", a3=", E.a3, ", a4=", E.a4, ", a6=", E.a6, ", disc=", E.disc, ", j=", E.j))
};

\\ Helper function to format polynomial coefficients
format_poly(P) = {
  local(coeffs);
  if(P == 0, return("[]"));
  coeffs = Vec(Polrev(P));
  return(Str(coeffs))
};

\\ =============================================================================
\\ SECTION 1: ellinit tests
\\ =============================================================================

\\ Test 1: ellinit with short Weierstrass form [a4, a6] over small prime
print("# ellinit([1, 2], 101)");
E1 = ellinit([1, 2], 101);
print("RESULT: ", format_curve(E1));
print();

\\ Test 2: ellinit with general Weierstrass form [a1, a2, a3, a4, a6]
print("# ellinit([0, 0, 0, 1, 1], 23)");
E2 = ellinit([0, 0, 0, 1, 1], 23);
print("RESULT: ", format_curve(E2));
print();

\\ Test 3: ellinit with non-zero a1, a2, a3
print("# ellinit([1, 2, 1, 3, 2], 37)");
E3 = ellinit([1, 2, 1, 3, 2], 37);
print("RESULT: ", format_curve(E3));
print();

\\ =============================================================================
\\ SECTION 2: ellcard tests (various primes)
\\ =============================================================================

\\ Test 4: ellcard basic
print("# ellcard(E) for E: y^2 = x^3 + x + 1 over F_23");
E4 = ellinit([0, 0, 0, 1, 1], 23);
card = ellcard(E4);
print("RESULT: ", card);
print();

\\ Test 5: ellcard for another curve
print("# ellcard(E) for E: y^2 = x^3 + 3x + 5 over F_101");
E5 = ellinit([0, 0, 0, 3, 5], 101);
card5 = ellcard(E5);
print("RESULT: ", card5);
print();

\\ Test 6: ellcard for secp256k1-like curve over small field
print("# ellcard(E) for E: y^2 = x^3 + 7 over F_97");
E6 = ellinit([0, 0, 0, 0, 7], 97);
card6 = ellcard(E6);
print("RESULT: ", card6);
print();

\\ Test 7: ellcard for larger prime
print("# ellcard for E: y^2 = x^3 + 2x + 3 over F_997");
E7 = ellinit([0, 0, 0, 2, 3], 997);
card7 = ellcard(E7);
print("RESULT: ", card7);
print();

\\ Test 8: ellcard for prime p=509
print("# ellcard for E: y^2 = x^3 + 5x + 7 over F_509");
E8 = ellinit([0, 0, 0, 5, 7], 509);
card8 = ellcard(E8);
print("RESULT: ", card8);
print();

\\ Test 9: ellcard for prime p=1009
print("# ellcard for E: y^2 = x^3 + x + 1 over F_1009");
E9 = ellinit([0, 0, 0, 1, 1], 1009);
card9 = ellcard(E9);
print("RESULT: ", card9);
print();

\\ Test 10: ellcard for curve with j=0 (a4=0)
print("# ellcard for E: y^2 = x^3 + 1 over F_67");
E10 = ellinit([0, 0, 0, 0, 1], 67);
card10 = ellcard(E10);
print("RESULT: ", card10);
print();

\\ Test 11: ellcard for curve with j=1728 (a6=0)
print("# ellcard for E: y^2 = x^3 + x over F_73");
E11 = ellinit([0, 0, 0, 1, 0], 73);
card11 = ellcard(E11);
print("RESULT: ", card11);
print();

\\ =============================================================================
\\ SECTION 3: Point addition and scalar multiplication
\\ =============================================================================

\\ Test 12: Point addition - elladd
print("# elladd on E: y^2 = x^3 + x + 1 over F_23");
E12 = ellinit([0, 0, 0, 1, 1], 23);
P1 = [0, 1];
P2 = [1, 7];
P3 = elladd(E12, P1, P2);
print("# P1 = ", format_point(P1));
print("# P2 = ", format_point(P2));
print("RESULT: ", format_point(P3));
print();

\\ Test 13: Point doubling
print("# elladd(P, P) - point doubling on E: y^2 = x^3 + x + 1 over F_23");
E13 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
P_doubled = elladd(E13, P, P);
print("# P = ", format_point(P));
print("RESULT: ", format_point(P_doubled));
print();

\\ Test 14: Point negation
print("# ellneg on E: y^2 = x^3 + x + 1 over F_23");
E14 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
neg_P = ellneg(E14, P);
print("# P = ", format_point(P));
print("RESULT: ", format_point(neg_P));
print();

\\ Test 15: Scalar multiplication with n=5
print("# ellmul on E: y^2 = x^3 + x + 1 over F_23 with n=5");
E15 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
nP = ellmul(E15, P, 5);
print("# P = ", format_point(P));
print("# n = 5");
print("RESULT: ", format_point(nP));
print();

\\ Test 16: Scalar multiplication with n=12
print("# ellmul with n=12 on E: y^2 = x^3 + x + 1 over F_23");
E16 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
nP = ellmul(E16, P, 12);
print("# P = ", format_point(P));
print("# n = 12");
print("RESULT: ", format_point(nP));
print();

\\ Test 17: Scalar multiplication with negative scalar
print("# ellmul with n=-3 on E: y^2 = x^3 + x + 1 over F_23");
E17 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
nP = ellmul(E17, P, -3);
print("# P = ", format_point(P));
print("# n = -3");
print("RESULT: ", format_point(nP));
print();

\\ Test 18: Large scalar multiplication
print("# ellmul with large scalar n=1000 over F_101");
E18 = ellinit([0, 0, 0, 3, 5], 101);
P = [1, 3];
result = ellmul(E18, P, 1000);
print("# P = ", format_point(P));
print("# n = 1000");
print("RESULT: ", format_point(result));
print();

\\ =============================================================================
\\ SECTION 4: Point order tests
\\ =============================================================================

\\ Test 19: ellorder for P=[0, 1]
print("# ellorder on E: y^2 = x^3 + x + 1 over F_23");
E19 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
ord = ellorder(E19, P);
print("# P = ", format_point(P));
print("RESULT: ", ord);
print();

\\ Test 20: ellorder for P=[1, 7]
print("# ellorder for P=[1, 7] on E: y^2 = x^3 + x + 1 over F_23");
E20 = ellinit([0, 0, 0, 1, 1], 23);
P = [1, 7];
ord = ellorder(E20, P);
print("# P = ", format_point(P));
print("RESULT: ", ord);
print();

\\ Test 21: Verify n*P = O
print("# ellmul(P, order) should give infinity");
E21 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
ord = ellorder(E21, P);
nP = ellmul(E21, P, ord);
print("# P = ", format_point(P));
print("# order = ", ord);
print("RESULT: ", format_point(nP));
print();

\\ =============================================================================
\\ SECTION 5: Identity element tests
\\ =============================================================================

\\ Test 22: elladd(P, O) = P
print("# elladd(P, O) should give P");
E22 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
O = [0];
result = elladd(E22, P, O);
print("# P = ", format_point(P));
print("RESULT: ", format_point(result));
print();

\\ Test 23: elladd(O, P) = P
print("# elladd(O, P) should give P");
E23 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
O = [0];
result = elladd(E23, O, P);
print("# P = ", format_point(P));
print("RESULT: ", format_point(result));
print();

\\ Test 24: elladd(P, -P) = O
print("# elladd(P, -P) should give O");
E24 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
neg_P = ellneg(E24, P);
result = elladd(E24, P, neg_P);
print("# P = ", format_point(P));
print("# -P = ", format_point(neg_P));
print("RESULT: ", format_point(result));
print();

\\ Test 25: ellmul with n=0 gives infinity
print("# ellmul with n=0 should give infinity");
E25 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
result = ellmul(E25, P, 0);
print("# P = ", format_point(P));
print("RESULT: ", format_point(result));
print();

\\ Test 26: ellmul with n=1 gives P
print("# ellmul with n=1 should give P");
E26 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
result = ellmul(E26, P, 1);
print("# P = ", format_point(P));
print("RESULT: ", format_point(result));
print();

\\ =============================================================================
\\ SECTION 6: Discriminant and j-invariant tests
\\ =============================================================================

\\ Test 27: discriminant
print("# discriminant of E: y^2 = x^3 + 2x + 3 over F_101");
E27 = ellinit([0, 0, 0, 2, 3], 101);
print("RESULT: ", E27.disc);
print();

\\ Test 28: j-invariant
print("# j-invariant of E: y^2 = x^3 + 2x + 3 over F_101");
E28 = ellinit([0, 0, 0, 2, 3], 101);
print("RESULT: ", E28.j);
print();

\\ =============================================================================
\\ SECTION 7: ellisoncurve tests
\\ =============================================================================

\\ Test 29: ellisoncurve for valid point
print("# ellisoncurve for P=[0,1] on E: y^2 = x^3 + x + 1 over F_23");
E29 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
result = ellisoncurve(E29, P);
print("# P = ", format_point(P));
print("RESULT: ", result);
print();

\\ Test 30: ellisoncurve for invalid point
print("# ellisoncurve for invalid P=[0,5] on E: y^2 = x^3 + x + 1 over F_23");
E30 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 5];
result = ellisoncurve(E30, P);
print("# P = ", format_point(P));
print("RESULT: ", result);
print();

\\ =============================================================================
\\ SECTION 8: ellgroup tests
\\ =============================================================================

\\ Test 31: ellgroup
print("# ellgroup for E: y^2 = x^3 + x + 1 over F_23");
E31 = ellinit([0, 0, 0, 1, 1], 23);
grp = ellgroup(E31);
print("RESULT: ", grp);
print();

\\ Test 32: ellgroup for larger curve
print("# ellgroup for E: y^2 = x^3 + 2x + 3 over F_101");
E32 = ellinit([0, 0, 0, 2, 3], 101);
grp = ellgroup(E32);
print("RESULT: ", grp);
print();

\\ Test 33: ellgroup for curve over F_997
print("# ellgroup for E: y^2 = x^3 + 5x + 7 over F_997");
E33 = ellinit([0, 0, 0, 5, 7], 997);
grp = ellgroup(E33);
print("RESULT: ", grp);
print();

\\ =============================================================================
\\ SECTION 9: Associativity and commutativity tests
\\ =============================================================================

\\ Test 34: Associativity check
print("# (P + Q) + R = P + (Q + R) associativity check");
E34 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
Q = [1, 7];
R = elladd(E34, P, Q);
S = elladd(E34, R, P);
T = elladd(E34, Q, P);
U = elladd(E34, P, T);
print("# P = ", format_point(P));
print("# Q = ", format_point(Q));
print("# (P+Q)+P = ", format_point(S));
print("RESULT: ", format_point(U));
print();

\\ Test 35: Verify -3P = -(3P)
print("# Verify -3P = -(3P)");
E35 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
neg3P = ellmul(E35, P, -3);
pos3P = ellmul(E35, P, 3);
neg_pos3P = ellneg(E35, pos3P);
print("# -3P = ", format_point(neg3P));
print("RESULT: ", format_point(neg_pos3P));
print();

\\ =============================================================================
\\ SECTION 10: elllog tests (discrete logarithm)
\\ =============================================================================

\\ Test 36: elllog basic
print("# elllog on E: y^2 = x^3 + x + 1 over F_23");
E36 = ellinit([0, 0, 0, 1, 1], 23);
G = [0, 1];
\\ Compute P = 5*G and then find log
P = ellmul(E36, G, 5);
dlog = elllog(E36, P, G);
print("# G = ", format_point(G));
print("# P = 5*G = ", format_point(P));
print("RESULT: ", dlog);
print();

\\ Test 37: elllog with different multiplier
print("# elllog on E: y^2 = x^3 + x + 1 over F_23 for 12*G");
E37 = ellinit([0, 0, 0, 1, 1], 23);
G = [0, 1];
P = ellmul(E37, G, 12);
dlog = elllog(E37, P, G);
print("# G = ", format_point(G));
print("# P = 12*G = ", format_point(P));
print("RESULT: ", dlog);
print();

\\ Test 38: elllog with order provided
print("# elllog on E: y^2 = x^3 + 3x + 5 over F_101");
E38 = ellinit([0, 0, 0, 3, 5], 101);
G = [1, 3];
ordG = ellorder(E38, G);
P = ellmul(E38, G, 17);
dlog = elllog(E38, P, G, ordG);
print("# G = ", format_point(G));
print("# P = 17*G = ", format_point(P));
print("# order(G) = ", ordG);
print("RESULT: ", dlog);
print();

\\ Test 39: elllog for identity
print("# elllog for identity point");
E39 = ellinit([0, 0, 0, 1, 1], 23);
G = [0, 1];
P = [0];  \\ Point at infinity
dlog = elllog(E39, P, G);
print("# G = ", format_point(G));
print("# P = O");
print("RESULT: ", dlog);
print();

\\ Test 40: elllog on larger field
print("# elllog on E over F_509");
E40 = ellinit([0, 0, 0, 2, 5], 509);
G = ellpointorder(E40, ellcard(E40))[1];  \\ Get a generator
ordG = ellorder(E40, G);
P = ellmul(E40, G, 73);
dlog = elllog(E40, P, G, ordG);
print("# G = ", format_point(G));
print("# P = 73*G = ", format_point(P));
print("RESULT: ", dlog);
print();

\\ =============================================================================
\\ SECTION 11: elldivpol tests (division polynomials)
\\ =============================================================================

\\ Test 41: elldivpol for n=2
print("# elldivpol(E, 2) for E: y^2 = x^3 + x + 1 over F_23");
E41 = ellinit([0, 0, 0, 1, 1], 23);
psi2 = elldivpol(E41, 2);
print("RESULT: ", format_poly(psi2));
print();

\\ Test 42: elldivpol for n=3
print("# elldivpol(E, 3) for E: y^2 = x^3 + x + 1 over F_23");
E42 = ellinit([0, 0, 0, 1, 1], 23);
psi3 = elldivpol(E42, 3);
print("RESULT: ", format_poly(psi3));
print();

\\ Test 43: elldivpol for n=4
print("# elldivpol(E, 4) for E: y^2 = x^3 + x + 1 over F_23");
E43 = ellinit([0, 0, 0, 1, 1], 23);
psi4 = elldivpol(E43, 4);
print("RESULT: ", format_poly(psi4));
print();

\\ Test 44: elldivpol for n=5
print("# elldivpol(E, 5) for E: y^2 = x^3 + x + 1 over F_23");
E44 = ellinit([0, 0, 0, 1, 1], 23);
psi5 = elldivpol(E44, 5);
print("RESULT: ", format_poly(psi5));
print();

\\ Test 45: elldivpol for different curve
print("# elldivpol(E, 3) for E: y^2 = x^3 + 2x + 3 over F_101");
E45 = ellinit([0, 0, 0, 2, 3], 101);
psi3 = elldivpol(E45, 3);
print("RESULT: ", format_poly(psi3));
print();

\\ Test 46: elldivpol for n=6
print("# elldivpol(E, 6) for E: y^2 = x^3 + x + 1 over F_23");
E46 = ellinit([0, 0, 0, 1, 1], 23);
psi6 = elldivpol(E46, 6);
print("RESULT: ", format_poly(psi6));
print();

\\ =============================================================================
\\ SECTION 12: elltatepairing tests
\\ =============================================================================

\\ Test 47: elltatepairing basic
print("# elltatepairing on E: y^2 = x^3 + x + 1 over F_23");
E47 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
Q = [1, 7];
ordP = ellorder(E47, P);
tate = elltatepairing(E47, P, Q, ordP);
print("# P = ", format_point(P));
print("# Q = ", format_point(Q));
print("# m = order(P) = ", ordP);
print("RESULT: ", lift(tate));
print();

\\ Test 48: elltatepairing with same point
print("# elltatepairing(P, P) on E: y^2 = x^3 + x + 1 over F_23");
E48 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
ordP = ellorder(E48, P);
tate = elltatepairing(E48, P, P, ordP);
print("# P = ", format_point(P));
print("# m = ", ordP);
print("RESULT: ", lift(tate));
print();

\\ Test 49: elltatepairing with infinity
print("# elltatepairing(O, Q)");
E49 = ellinit([0, 0, 0, 1, 1], 23);
P = [0];  \\ Point at infinity
Q = [0, 1];
tate = elltatepairing(E49, P, Q, 28);
print("# P = O");
print("# Q = ", format_point(Q));
print("RESULT: ", lift(tate));
print();

\\ Test 50: elltatepairing on larger field
print("# elltatepairing on E over F_101");
E50 = ellinit([0, 0, 0, 2, 3], 101);
P = [12, 43];  \\ Check if this is on curve: 43^2 = 1849 = 30 mod 101; 12^3 + 2*12 + 3 = 1728 + 24 + 3 = 1755 = 38 mod 101. Let's find a valid point
\\ Try x=0: y^2 = 3 mod 101, sqrt(3) doesn't exist in F_101
\\ Try x=1: y^2 = 1 + 2 + 3 = 6 mod 101
\\ Try x=2: y^2 = 8 + 4 + 3 = 15 mod 101
\\ Try x=5: y^2 = 125 + 10 + 3 = 138 = 37 mod 101
\\ Try x=10: y^2 = 1000 + 20 + 3 = 1023 = 13 mod 101
\\ Try x=11: y^2 = 1331 + 22 + 3 = 1356 = 42 mod 101
\\ Use random point finding instead
P = [73, 35];  \\ 35^2 = 1225 = 13 mod 101; 73^3 + 2*73 + 3 = 389017 + 146 + 3 = 389166 = 13 mod 101 - yes!
Q = [28, 60];  \\ 60^2 = 3600 = 66 mod 101; 28^3 + 2*28 + 3 = 21952 + 56 + 3 = 22011 = 66 mod 101 - yes!
ordP = ellorder(E50, P);
tate = elltatepairing(E50, P, Q, ordP);
print("# P = ", format_point(P));
print("# Q = ", format_point(Q));
print("# m = ", ordP);
print("RESULT: ", lift(tate));
print();

\\ =============================================================================
\\ SECTION 13: ellweilpairing tests
\\ =============================================================================

\\ Test 51: ellweilpairing basic
print("# ellweilpairing on E: y^2 = x^3 + x + 1 over F_23");
E51 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
Q = [1, 7];
\\ For Weil pairing, both points must have the same order
ordP = ellorder(E51, P);
ordQ = ellorder(E51, Q);
m = if(ordP <= ordQ, ordP, ordQ);
\\ Scale points to have order m
if(ordP > m, P = ellmul(E51, P, ordP/m));
if(ordQ > m, Q = ellmul(E51, Q, ordQ/m));
weil = ellweilpairing(E51, P, Q, m);
print("# P = ", format_point(P));
print("# Q = ", format_point(Q));
print("# m = ", m);
print("RESULT: ", lift(weil));
print();

\\ Test 52: ellweilpairing symmetry check e(P,Q) * e(Q,P) = 1
print("# ellweilpairing symmetry: e(P,Q) * e(Q,P)");
E52 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
Q = [1, 7];
\\ Use a common order
ord = gcd(ellorder(E52, P), ellorder(E52, Q));
w1 = ellweilpairing(E52, P, Q, ord);
w2 = ellweilpairing(E52, Q, P, ord);
product = lift(w1 * w2);
print("# e(P,Q) = ", lift(w1));
print("# e(Q,P) = ", lift(w2));
print("RESULT: ", product);
print();

\\ Test 53: ellweilpairing(P, P) = 1
print("# ellweilpairing(P, P) should be 1");
E53 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
ordP = ellorder(E53, P);
weil = ellweilpairing(E53, P, P, ordP);
print("# P = ", format_point(P));
print("RESULT: ", lift(weil));
print();

\\ =============================================================================
\\ SECTION 14: Known Curves Tests (secp256k1-like, P-256-like with small primes)
\\ =============================================================================

\\ Test 54: secp256k1-like curve y^2 = x^3 + 7 over small prime
print("# secp256k1-like: y^2 = x^3 + 7 over F_2063");
E54 = ellinit([0, 0, 0, 0, 7], 2063);
card = ellcard(E54);
print("RESULT: ", card);
print();

\\ Test 55: ellcard for a curve similar structure to P-256 (a=-3)
print("# P-256-like: y^2 = x^3 - 3x + b over F_1997");
\\ P-256 has a = -3, let's use small prime
E55 = ellinit([0, 0, 0, -3, 5], 1997);
card = ellcard(E55);
print("RESULT: ", card);
print();

\\ Test 56: Curve with a=-3 (common in standards)
print("# Curve y^2 = x^3 - 3x + 1 over F_127");
E56 = ellinit([0, 0, 0, -3, 1], 127);
card = ellcard(E56);
grp = ellgroup(E56);
print("# card = ", card);
print("RESULT: ", grp);
print();

\\ =============================================================================
\\ SECTION 15: Random Curve Tests (deterministic with seed)
\\ =============================================================================

\\ Test 57: Random curve 1
print("# Random curve 1: a4=17, a6=31 over F_229");
E57 = ellinit([0, 0, 0, 17, 31], 229);
card = ellcard(E57);
print("RESULT: ", card);
print();

\\ Test 58: Random curve 2
print("# Random curve 2: a4=42, a6=13 over F_307");
E58 = ellinit([0, 0, 0, 42, 13], 307);
card = ellcard(E58);
print("RESULT: ", card);
print();

\\ Test 59: Random curve 3 with point operations
print("# Random curve 3 point operation: a4=7, a6=11 over F_167");
E59 = ellinit([0, 0, 0, 7, 11], 167);
\\ Find a point on the curve
P = [2, 0];  \\ Check: 0 = 8 + 14 + 11 = 33 mod 167. Not on curve
\\ Try to find a valid point
found = 0;
for(x = 0, 166,
  y2 = Mod(x^3 + 7*x + 11, 167);
  if(issquare(y2),
    y = sqrt(y2);
    P = [x, lift(y)];
    found = 1;
    break
  )
);
if(found,
  Q = ellmul(E59, P, 7);
  print("# P = ", format_point(P));
  print("# 7*P = ", format_point(Q));
  print("RESULT: ", format_point(Q))
,
  print("RESULT: NO_POINT_FOUND")
);
print();

\\ Test 60: Random curve with elllog
print("# Random curve elllog test: a4=5, a6=9 over F_251");
E60 = ellinit([0, 0, 0, 5, 9], 251);
\\ Find a generator
card = ellcard(E60);
G = ellpointorder(E60, card)[1];
ordG = ellorder(E60, G);
\\ Compute 23*G and find dlog
P = ellmul(E60, G, 23);
dlog = elllog(E60, P, G, ordG);
print("# G = ", format_point(G));
print("# P = 23*G");
print("RESULT: ", dlog);
print();

\\ =============================================================================
\\ SECTION 16: Edge Cases
\\ =============================================================================

\\ Test 61: Very small prime p=5
print("# Edge case: very small prime p=5");
E61 = ellinit([0, 0, 0, 1, 1], 5);
card = ellcard(E61);
print("RESULT: ", card);
print();

\\ Test 62: Prime p=7
print("# Edge case: prime p=7");
E62 = ellinit([0, 0, 0, 0, 1], 7);  \\ y^2 = x^3 + 1 (non-singular)
card = ellcard(E62);
print("RESULT: ", card);
print();

\\ Test 63: Prime p=11
print("# Edge case: prime p=11");
E63 = ellinit([0, 0, 0, 3, 4], 11);
card = ellcard(E63);
grp = ellgroup(E63);
print("# card = ", card);
print("RESULT: ", grp);
print();

\\ Test 64: Curve where group might be non-cyclic
print("# Curve with potential non-cyclic group over F_31");
E64 = ellinit([0, 0, 0, 1, 0], 31);
card = ellcard(E64);
grp = ellgroup(E64);
print("# card = ", card);
print("RESULT: ", grp);
print();

\\ =============================================================================
\\ SECTION 17: More elldivpol tests
\\ =============================================================================

\\ Test 65: elldivpol for n=7
print("# elldivpol(E, 7) for E: y^2 = x^3 + 2x + 3 over F_101");
E65 = ellinit([0, 0, 0, 2, 3], 101);
psi7 = elldivpol(E65, 7);
print("RESULT: ", format_poly(psi7));
print();

\\ Test 66: elldivpol for n=8
print("# elldivpol(E, 8) for E: y^2 = x^3 + 1x + 2 over F_47");
E66 = ellinit([0, 0, 0, 1, 2], 47);
psi8 = elldivpol(E66, 8);
print("RESULT: ", format_poly(psi8));
print();

\\ Test 67: elldivpol degree check for n=5
\\ For odd n, degree should be (n^2-1)/2 = 12
print("# elldivpol(E, 5) degree check");
E67 = ellinit([0, 0, 0, 1, 1], 23);
psi5 = elldivpol(E67, 5);
deg = poldegree(psi5);
print("# degree = ", deg);
print("RESULT: ", deg);
print();

\\ =============================================================================
\\ SECTION 18: Additional pairing tests
\\ =============================================================================

\\ Test 68: elltatepairing bilinearity check
print("# elltatepairing bilinearity: e(2P, Q) = e(P, Q)^2");
E68 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
Q = [1, 7];
ordP = ellorder(E68, P);
tP2Q = elltatepairing(E68, ellmul(E68, P, 2), Q, ordP);
tPQ = elltatepairing(E68, P, Q, ordP);
tPQ2 = tPQ^2;
print("# e(2P, Q) = ", lift(tP2Q));
print("# e(P, Q)^2 = ", lift(tPQ2));
print("RESULT: ", lift(tP2Q) == lift(tPQ2));
print();

\\ Test 69: ellweilpairing bilinearity check
print("# ellweilpairing bilinearity check");
E69 = ellinit([0, 0, 0, 1, 1], 23);
P = [0, 1];
Q = [1, 7];
ord = gcd(ellorder(E69, P), ellorder(E69, Q));
P2 = ellmul(E69, P, 2);
wP2Q = ellweilpairing(E69, P2, Q, ord);
wPQ = ellweilpairing(E69, P, Q, ord);
wPQ2 = wPQ^2;
print("# e(2P, Q) = ", lift(wP2Q));
print("# e(P, Q)^2 = ", lift(wPQ2));
print("RESULT: ", lift(wP2Q) == lift(wPQ2));
print();

\\ =============================================================================
\\ SECTION 19: Verify elllog consistency
\\ =============================================================================

\\ Test 70: elllog verify with ellmul
print("# elllog verification: ellmul(G, elllog(P,G)) == P");
E70 = ellinit([0, 0, 0, 3, 5], 101);
G = [1, 3];
ordG = ellorder(E70, G);
P = ellmul(E70, G, 37);
dlog = elllog(E70, P, G, ordG);
P_verify = ellmul(E70, G, dlog);
print("# P = ", format_point(P));
print("# dlog = ", dlog);
print("# G^dlog = ", format_point(P_verify));
print("RESULT: ", P == P_verify);
print();

\\ =============================================================================
\\ SECTION 20: Large field tests
\\ =============================================================================

\\ Test 71: ellcard for prime p=4999
print("# ellcard for E over F_4999");
E71 = ellinit([0, 0, 0, 7, 11], 4999);
card = ellcard(E71);
print("RESULT: ", card);
print();

\\ Test 72: ellorder on larger field
print("# ellorder on E over F_4999");
E72 = ellinit([0, 0, 0, 7, 11], 4999);
\\ Find a point
P = [0, 0];
found = 0;
for(x = 0, 100,
  y2 = Mod(x^3 + 7*x + 11, 4999);
  if(issquare(y2),
    y = sqrt(y2);
    P = [x, lift(y)];
    found = 1;
    break
  )
);
if(found,
  ord = ellorder(E72, P);
  print("# P = ", format_point(P));
  print("RESULT: ", ord)
,
  print("RESULT: NO_POINT_FOUND")
);
print();

\\ =============================================================================
\\ SECTION 21: More comprehensive elllog tests
\\ =============================================================================

\\ Test 73: elllog with k=1 (P = G)
print("# elllog when P = G (expect 1)");
E73 = ellinit([0, 0, 0, 1, 1], 23);
G = [0, 1];
dlog = elllog(E73, G, G);
print("# G = ", format_point(G));
print("# P = G");
print("RESULT: ", dlog);
print();

\\ Test 74: elllog with full group order
print("# elllog with various multipliers on F_67");
E74 = ellinit([0, 0, 0, 2, 5], 67);
G = [0, 0];
\\ Find a generator
card = ellcard(E74);
G = ellpointorder(E74, card)[1];
ordG = ellorder(E74, G);
\\ Test with k = ordG - 1
k = ordG - 1;
P = ellmul(E74, G, k);
dlog = elllog(E74, P, G, ordG);
print("# k = ", k);
print("# ordG = ", ordG);
print("RESULT: ", dlog);
print();

\\ =============================================================================
\\ SECTION 22: elldivpol evaluation tests
\\ =============================================================================

\\ Test 75: Evaluate psi_n at an n-torsion point (should be 0)
print("# psi_n evaluated at n-torsion point");
E75 = ellinit([0, 0, 0, 1, 1], 23);
\\ Find a 7-torsion point (if exists)
card = ellcard(E75);
\\ card = 28 = 4 * 7, so 7-torsion points exist
G = [0, 1];
ordG = ellorder(E75, G);
cofactor = ordG / gcd(ordG, 7);
P7 = ellmul(E75, G, cofactor);  \\ This should have order dividing 7
ordP7 = ellorder(E75, P7);
if(ordP7 == 7,
  psi7 = elldivpol(E75, 7);
  val = subst(psi7, x, P7[1]);
  print("# P7 has order ", ordP7);
  print("# psi_7(P7.x) = ", lift(val));
  print("RESULT: ", lift(val) == 0)
,
  print("RESULT: NO_7_TORSION")
);
print();

\\ =============================================================================
\\ SECTION 23: Additional ellcard tests for comprehensive coverage
\\ =============================================================================

\\ Test 76: ellcard for supersingular-like curve (trace = 0 over some extension)
print("# ellcard for E: y^2 = x^3 + 1 over F_13");
E76 = ellinit([0, 0, 0, 0, 1], 13);
card = ellcard(E76);
print("RESULT: ", card);
print();

\\ Test 77: ellcard for prime p = 1021
print("# ellcard for E: y^2 = x^3 + 2x + 5 over F_1021");
E77 = ellinit([0, 0, 0, 2, 5], 1021);
card = ellcard(E77);
print("RESULT: ", card);
print();

\\ Test 78: ellcard for curve with specific trace
print("# ellcard for E: y^2 = x^3 + 11x + 17 over F_257");
E78 = ellinit([0, 0, 0, 11, 17], 257);
card = ellcard(E78);
print("RESULT: ", card);
print();

\\ =============================================================================
\\ SECTION 24: Additional elllog edge cases
\\ =============================================================================

\\ Test 79: elllog where base point has composite order
print("# elllog with composite order base point");
E79 = ellinit([0, 0, 0, 3, 5], 101);
G = [1, 3];
ordG = ellorder(E79, G);
\\ Use a point that is a multiple of G
P = ellmul(E79, G, ordG - 5);
dlog = elllog(E79, P, G, ordG);
print("# ordG = ", ordG);
print("# expected = ", ordG - 5);
print("RESULT: ", dlog);
print();

\\ =============================================================================
\\ SECTION 25: Verify ellmul(P, ellorder(P)) = O for various points
\\ =============================================================================

\\ Test 80: Order verification for multiple points
print("# Verify ellmul(P, order(P)) = O");
E80 = ellinit([0, 0, 0, 5, 7], 127);
\\ Test with random-ish point
found = 0;
for(x = 1, 50,
  y2 = Mod(x^3 + 5*x + 7, 127);
  if(issquare(y2),
    y = sqrt(y2);
    P = [x, lift(y)];
    ord = ellorder(E80, P);
    Q = ellmul(E80, P, ord);
    if(Q == [0],
      print("# P = ", format_point(P));
      print("# order = ", ord);
      print("RESULT: ", format_point(Q));
      found = 1;
      break
    )
  )
);
if(!found, print("RESULT: NO_POINT_FOUND"));
print();

quit;
