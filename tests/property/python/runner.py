#!/usr/bin/env sage
"""
Property test runner for SageMath.
Executes operations with deterministic seeding and outputs results.

This script reads test case definitions from stdin (JSON format) and
executes them in SageMath, outputting results in a format that can be
compared with the TypeScript implementation.

Usage:
    sage runner.py < test_cases.json
    sage runner.py --seed 42 < test_cases.json

Output format:
    {"function": "gcd", "args": [12, 8], "result": "4", "seed": 42}

Areas
-----
This runner contains NO area-specific code.  Every area lives in its own module
under ``areas/`` and is discovered by globbing that directory, so several agents
can add areas in parallel without ever editing this file.  See
``tests/property/README.md``.
"""

import glob
import importlib
import json
import os
import sys
import random as py_random
from sage.all import *

# Importing areas/ would otherwise leave a __pycache__ directory inside the repo
# (which .gitignore does not cover); the runner is short-lived so caching buys
# nothing.
sys.dont_write_bytecode = True

# Make ``areas`` importable as a package no matter which directory `sage` was
# invoked from (compare.ts runs us with cwd = project root).
_RUNNER_DIR = os.path.dirname(os.path.abspath(__file__))
if _RUNNER_DIR not in sys.path:
    sys.path.insert(0, _RUNNER_DIR)

AREAS_DIR = os.path.join(_RUNNER_DIR, 'areas')

# name -> loaded {'functions': ..., 'formatters': ...}; filled in on demand.
_LOADED_AREAS = {}


def discover_area_names():
    """
    Return the sorted names of the area modules present in ``areas/``.

    An area is any ``areas/<area>.py`` whose basename does not start with ``_``
    (underscore-prefixed modules are shared helpers, not areas).  The name must
    match the ``module`` field of ``cases/<area>.cases.json``.
    """
    names = []
    for path in sorted(glob.glob(os.path.join(AREAS_DIR, '*.py'))):
        name = os.path.splitext(os.path.basename(path))[0]
        if name.startswith('_'):
            continue
        names.append(name)
    return names


def load_area(name):
    """
    Import ``areas/<name>.py`` and return its dispatch tables.

    Every area module must export a ``FUNCTIONS`` dict mapping the ``function``
    names used in ``cases/<name>.cases.json`` to callables.  It may optionally
    export a ``FORMATTERS`` dict mapping a function name to a ``result -> str``
    callable, for results the generic ``format_result`` below cannot render.

    Areas are imported lazily and one at a time, so a broken area module can
    never take down the areas owned by other agents.

    Returns:
        dict: {'functions': dict, 'formatters': dict}
    """
    if name in _LOADED_AREAS:
        return _LOADED_AREAS[name]

    module = importlib.import_module('areas.' + name)
    functions = getattr(module, 'FUNCTIONS', None)
    if not isinstance(functions, dict):
        raise ValueError(
            f"Area module areas/{name}.py must export a FUNCTIONS dict "
            f"mapping function names to callables"
        )
    formatters = getattr(module, 'FORMATTERS', {})
    _LOADED_AREAS[name] = {'functions': functions, 'formatters': formatters}
    return _LOADED_AREAS[name]


def set_seed(seed):
    """Set random seed for reproducibility across both SageMath and Python random."""
    set_random_seed(seed)
    py_random.seed(seed)


def generate_arg(generator_spec, seed):
    """
    Generate an argument based on the generator specification.

    Supported generators:
    - randomBigint(min, max): Random integer in range [min, max]
    - randomPrime(min, max): Random prime in range [min, max]
    - fixedValue(value): Fixed value
    - randomList(generator, length): List of random values

    Note: Uses Python's random module (Mersenne Twister) to match TypeScript.
    """
    if generator_spec.startswith('randomBigint('):
        # Parse: randomBigint(min, max)
        params = generator_spec[13:-1]  # Remove 'randomBigint(' and ')'
        parts = params.split(',')
        min_val = int(parts[0].strip())
        max_val = int(parts[1].strip())
        # Use Python's random.randint to match TypeScript's Mersenne Twister
        return Integer(py_random.randint(min_val, max_val))

    elif generator_spec.startswith('randomPrime('):
        # Parse: randomPrime(min, max)
        params = generator_spec[12:-1]
        parts = params.split(',')
        min_val = int(parts[0].strip())
        max_val = int(parts[1].strip())
        # Use Python's random with rejection sampling to match TypeScript
        max_attempts = 10000
        for _ in range(max_attempts):
            candidate = py_random.randint(min_val, max_val)
            if is_prime(candidate):
                return Integer(candidate)
        # Fallback: find next prime from min
        p = next_prime(min_val - 1)
        return p

    elif generator_spec.startswith('fixedValue('):
        # Parse: fixedValue(value) - value can be integer or array
        value_str = generator_spec[11:-1].strip()
        if value_str.startswith('[') and value_str.endswith(']'):
            # Parse array: [1, 2, 3]
            inner = value_str[1:-1].strip()
            if not inner:
                return []
            parts = inner.split(',')
            return [Integer(p.strip()) for p in parts]
        return Integer(value_str)

    elif generator_spec.startswith('randomList('):
        # Parse: randomList(generator, length)
        params = generator_spec[11:-1]
        # Find the last comma that separates length
        last_comma = params.rfind(',')
        inner_generator = params[:last_comma].strip()
        length = int(params[last_comma+1:].strip())
        return [generate_arg(inner_generator, seed + i) for i in range(length)]

    else:
        # Assume it's a literal value
        return Integer(generator_spec)


def format_result(result):
    """
    Format a SageMath result as a string for comparison.

    This ensures consistent string representation across Python and TypeScript.
    """
    if result is None:
        return 'null'
    if isinstance(result, bool):
        return 'True' if result else 'False'
    elif isinstance(result, Factorization):
        # Format: "2^2 * 3" or "1" for empty
        if len(result) == 0:
            return '1'
        parts = []
        for (p, e) in result:
            if e == 1:
                parts.append(str(p))
            else:
                parts.append(f'{p}^{e}')
        return ' * '.join(parts)
    elif isinstance(result, (list, tuple)):
        # Format lists/tuples to match Python's str() output
        if isinstance(result, tuple):
            return '(' + ', '.join(format_result(x) for x in result) + ')'
        else:
            return '[' + ', '.join(format_result(x) for x in result) + ']'
    elif hasattr(result, '__iter__') and not isinstance(result, str):
        # Other iterables
        return '[' + ', '.join(format_result(x) for x in result) + ']'
    else:
        return str(result)


def execute_function(module, function_name, args):
    """
    Execute a SageMath function with the given arguments.

    Args:
        module: Area name (e.g., "arith", "elliptic_curves"), matching areas/<module>.py
        function_name: Function name (e.g., "gcd", "is_prime")
        args: List of arguments

    Returns:
        The result of the function call
    """
    if module not in discover_area_names():
        raise ValueError(f"Unknown module: {module}")

    functions = load_area(module)['functions']
    if function_name not in functions:
        raise ValueError(f"Unknown function: {module}.{function_name}")

    func = functions[function_name]
    return func(*args)


def format_for(module, function_name, result):
    """Format a result, honouring an area-specific FORMATTERS override if present."""
    if module not in discover_area_names():
        return format_result(result)
    formatter = load_area(module)['formatters'].get(function_name)
    if formatter is not None:
        return formatter(result)
    return format_result(result)


def run_test_case(case, seed):
    """
    Run a single test case with the given seed.

    Args:
        case: Test case definition dict
        seed: Random seed to use

    Returns:
        Result dict with function, args, result, and seed
    """
    set_seed(seed)

    module = case.get('module', 'arith')
    function_name = case['function']
    arg_generators = case.get('argGenerators', [])

    # Generate arguments
    args = []
    for i, gen in enumerate(arg_generators):
        arg = generate_arg(gen, seed + i * 1000)
        args.append(arg)

    # Execute function
    try:
        result = execute_function(module, function_name, args)
        formatted_result = format_for(module, function_name, result)
        error = None
    except Exception as e:
        formatted_result = None
        error = str(e)

    return {
        'function': function_name,
        'args': [str(a) for a in args],
        'result': formatted_result,
        'error': error,
        'seed': seed,
    }


def run_test_suite(test_suite):
    """
    Run all test cases in a test suite.

    Args:
        test_suite: Dict with 'module' and 'cases' keys

    Returns:
        List of result dicts
    """
    module = test_suite.get('module', 'arith')
    cases = test_suite.get('cases', [])
    results = []

    for case in cases:
        case['module'] = module
        seeds = case.get('seeds', [42])

        for seed in seeds:
            result = run_test_case(case, seed)
            results.append(result)

    return results


def main():
    """Main entry point."""
    # Parse command line arguments
    seed_override = None
    for i, arg in enumerate(sys.argv[1:], 1):
        if arg == '--seed' and i < len(sys.argv) - 1:
            seed_override = int(sys.argv[i + 1])

    # Read test cases from stdin
    input_data = sys.stdin.read()

    if not input_data.strip():
        print("Error: No input provided. Pass test cases as JSON via stdin.", file=sys.stderr)
        sys.exit(1)

    try:
        test_suite = json.loads(input_data)
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON: {e}", file=sys.stderr)
        sys.exit(1)

    # Override seeds if specified
    if seed_override is not None:
        for case in test_suite.get('cases', []):
            case['seeds'] = [seed_override]

    # Run tests
    results = run_test_suite(test_suite)

    # Output results as JSON
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
