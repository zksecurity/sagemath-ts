"""State/copy fixtures driven by the bundled generic.py implementation."""
import json
from sage.all import QuadraticField


def iterator_state(g, n, indexed, copy_kind, mutate_before, mutate_op, throw_at):
    indexed = bool(indexed)
    copies = calls = 0
    class Box:
        def __init__(self, value): self.value = value % 11
        def __eq__(self, other): return self.value == other.value
    if copy_kind:
        def cp(self):
            nonlocal copies
            copies += 1
            if copy_kind == 2 or (copy_kind == 6 and self.value == 5):
                raise ValueError('copy blocked')
            if copy_kind == 3 or (copy_kind == 5 and self.value == 5):
                return None
            return Box(self.value)
        Box.__copy__ = cp
    p, p0 = Box(2), Box(5)
    def op(x, y):
        nonlocal calls
        calls += 1
        if calls == throw_at: raise ValueError('step blocked')
        if mutate_op:
            x.value = (x.value + y.value) % 11
            return x
        return Box(x.value + y.value)
    def state(): return [calls, copies, p.value, p0.value]
    events = []
    try:
        it = g.multiples(p, n, p0, indexed, operation='other', op=op)
        events.append(['construct', *state()])
    except Exception as e:
        events.append(['construct-error', type(e).__name__ + ': ' + str(e), *state()])
        return json.dumps(events, separators=(',', ':'))
    if mutate_before: p.value, p0.value = 9, 8
    for _ in range(max(n, 0) + 2):
        try:
            value = next(it)
            output = [int(value[0]), value[1].value] if indexed else value.value
            events.append(['value', output, *state()])
        except StopIteration: events.append(['done', *state()])
        except Exception as e:
            events.append(['next-error', type(e).__name__ + ': ' + str(e), *state()])
    return json.dumps(events, separators=(',', ':'))


def iterator_standard(g, n, mode, supplied, parent_kind, mutate):
    # parent_kind selects a TS property/method representation of Sage's parent().
    trace = []
    class Parent:
        def zero(self): trace.append(['zero']); return Box(0)
        def one(self): trace.append(['one']); return Box(1)
    ring = Parent()
    class Box:
        def __init__(self, value): self.value = value % 11
        def parent(self): return ring
        def __add__(self, other):
            trace.append(['add', self.value, other.value])
            return Box(self.value + other.value)
        def __mul__(self, other):
            trace.append(['mul', self.value, other.value])
            return Box(self.value * other.value)
    p, p0 = Box(2), Box(5)
    it = g.multiples(p, n, p0 if supplied else None, operation='+' if mode == 0 else '*')
    events = [['construct', list(trace)]]
    if mutate: p.value, p0.value = 9, 8
    for _ in range(n + 1):
        try:
            value = next(it)
            events.append(['value', value.value, list(trace), p.value, p0.value])
        except StopIteration: events.append(['done', list(trace), p.value, p0.value])
    return json.dumps(events, separators=(',', ':'))


def parent_field(g, d, mode, k):
    K = QuadraticField(d, 'a')
    a, b = K.gen(), K.gen() + 1
    try:
        oper, identity, inverse, op = g._parse_group_def(K, '+' if mode == 0 else '*', None, None, None)
        values = [identity, inverse(a), op(a, b), g._power_func(oper, identity, inverse, op)(a, k)]
        result = [[str(c) for c in value.list()] for value in values] + [identity == identity, a == identity]
        error = None
    except Exception as e:
        result, error = None, type(e).__name__ + ': ' + str(e)
    return json.dumps([result, error], separators=(',', ':'))
