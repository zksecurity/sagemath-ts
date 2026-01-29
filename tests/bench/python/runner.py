#!/usr/bin/env sage
"""
Benchmark runner for SageMath.

Usage:
    sage tests/bench/python/runner.py --input tests/bench/inputs/arith.bench.inputs.json
    cat tests/bench/inputs/arith.bench.inputs.json | sage tests/bench/python/runner.py
"""

import json
import sys
import time
from sage.all import *
from sage.env import SAGE_VERSION


def hex_to_int(value):
    return int(value.replace(" ", "").replace("\n", ""), 16)


# Curve parameters (prime field short Weierstrass form)
SECP256K1 = {
    "p": hex_to_int("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F"),
    "a": 0,
    "b": 7,
    "gx": hex_to_int("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798"),
    "gy": hex_to_int("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8"),
}

P256 = {
    "p": hex_to_int("FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF"),
    "a": hex_to_int("FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC"),
    "b": hex_to_int("5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B"),
    "gx": hex_to_int("6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296"),
    "gy": hex_to_int("4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5"),
}

BLS12_381_G1 = {
    "p": hex_to_int(
        "1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F624"
        "1EABFFFEB153FFFFB9FEFFFFFFFFAAAB"
    ),
    "a": 0,
    "b": 4,
    "gx": int(
        "3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507"
    ),
    "gy": int(
        "1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569"
    ),
}

SECP256K1_FIELD = GF(SECP256K1["p"])
SECP256K1_CURVE = EllipticCurve(SECP256K1_FIELD, [SECP256K1["a"], SECP256K1["b"]])
SECP256K1_G = SECP256K1_CURVE([SECP256K1["gx"], SECP256K1["gy"]])

P256_FIELD = GF(P256["p"])
P256_CURVE = EllipticCurve(P256_FIELD, [P256["a"], P256["b"]])
P256_G = P256_CURVE([P256["gx"], P256["gy"]])

BLS12_381_FIELD = GF(BLS12_381_G1["p"])
BLS12_381_CURVE = EllipticCurve(BLS12_381_FIELD, [BLS12_381_G1["a"], BLS12_381_G1["b"]])
BLS12_381_G1_POINT = BLS12_381_CURVE([BLS12_381_G1["gx"], BLS12_381_G1["gy"]])


def parse_arg(arg):
    if isinstance(arg, list):
        return [Integer(value) for value in arg]
    return Integer(arg)


def percentile(values, p):
    if not values:
        return None
    values = sorted(values)
    index = (len(values) - 1) * p
    lower = int(index)
    upper = int(index + 1)
    if upper >= len(values):
        return values[lower]
    weight = index - lower
    return values[lower] * (1 - weight) + values[upper] * weight


def execute_function(module, function_name, args):
    function_map = {
        'arith': {
            'gcd': gcd,
            'lcm': lcm,
            'xgcd': xgcd,
            'factor': factor,
            'is_prime': is_prime,
            'is_prime_power': is_prime_power,
            'next_prime': next_prime,
            'previous_prime': previous_prime,
            'euler_phi': euler_phi,
            'radical': radical,
            'moebius': moebius,
            'kronecker_symbol': kronecker_symbol,
            'legendre_symbol': legendre_symbol,
            'jacobi_symbol': jacobi_symbol,
            'power_mod': power_mod,
            'inverse_mod': inverse_mod,
            'crt': crt,
            'isqrt': isqrt,
            'is_square': is_square,
            'is_squarefree': is_squarefree,
            'divisors': divisors,
            'number_of_divisors': number_of_divisors,
            'sigma': sigma,
            'prime_range': prime_range,
            'trial_division': trial_division,
            'squarefree_part': squarefree_part,
            'prime_factors': prime_factors,
            'valuation': valuation,
        },
        'bench_crypto': {
            'rsa_2048_pow': lambda m, e, n: power_mod(m % n, e, n),
            'secp256k1_mul': lambda k: SECP256K1_G * k,
            'p256_mul': lambda k: P256_G * k,
            'bls12_381_g1_mul': lambda k: BLS12_381_G1_POINT * k,
        },
    }

    if module not in function_map:
        raise ValueError(f"Unknown module: {module}")

    if function_name not in function_map[module]:
        raise ValueError(f"Unknown function: {module}.{function_name}")

    func = function_map[module][function_name]
    return func(*args)


def run_case(module_name, bench_case):
    inputs = [[parse_arg(arg) for arg in args] for args in bench_case['inputs']]
    warmup = int(bench_case.get('warmup', 0))
    iterations = int(bench_case.get('iterations', 1))
    repeats = int(bench_case.get('repeats', 1))

    total_ns = []
    per_call_ns = []
    error = None

    try:
        for i in range(warmup):
            args = inputs[i % len(inputs)]
            execute_function(module_name, bench_case['function'], args)

        for _ in range(repeats):
            start = time.perf_counter_ns()
            for i in range(iterations):
                args = inputs[i % len(inputs)]
                execute_function(module_name, bench_case['function'], args)
            end = time.perf_counter_ns()
            elapsed = end - start
            total_ns.append(str(elapsed))
            per_call_ns.append(elapsed / iterations)
    except Exception as exc:
        error = str(exc)
        total_ns = []
        per_call_ns = []

    return {
        'id': bench_case.get('id', bench_case['function']),
        'function': bench_case['function'],
        'inputs': len(inputs),
        'iterations': iterations,
        'warmup': warmup,
        'repeats': repeats,
        'total_ns': total_ns,
        'per_call_ns': per_call_ns,
        'p50_ns': percentile(per_call_ns, 0.5),
        'p90_ns': percentile(per_call_ns, 0.9),
        'error': error,
    }


def read_input(input_path):
    if input_path:
        with open(input_path, 'r', encoding='utf-8') as handle:
            return handle.read()
    return sys.stdin.read()


def main():
    input_path = None
    args = sys.argv[1:]
    for i, arg in enumerate(args):
        if arg == '--input' and i + 1 < len(args):
            input_path = args[i + 1]

    raw = read_input(input_path)
    if not raw.strip():
        print('Error: No input provided. Pass benchmark inputs via --input or stdin.', file=sys.stderr)
        sys.exit(1)

    suite = json.loads(raw)
    results = [run_case(suite['module'], bench_case) for bench_case in suite['cases']]

    output = {
        'module': suite['module'],
        'inputFile': input_path,
        'runtime': f"sage {SAGE_VERSION}",
        'timestamp': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'results': results,
    }

    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
