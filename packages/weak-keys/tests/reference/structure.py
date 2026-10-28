"""Independent small-integer oracle for fresh property cases; no Sage/third-party deps."""
import json
import math
import sys


def classify(n):
    prime = n > 1 and all(n % divisor for divisor in range(2, math.isqrt(n) + 1))
    power = False
    for base in range(2, math.isqrt(n) + 1):
        value = base * base
        while value < n:
            value *= base
        if value == n:
            power = True
            break
    return {"prime": prime, "power": power}


print(json.dumps([classify(n) for n in json.load(sys.stdin)]))
