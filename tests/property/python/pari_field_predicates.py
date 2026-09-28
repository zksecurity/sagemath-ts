"""Bundled PARI finite-field resultant/norm/square-predicate comparisons."""
import atexit
import hashlib
import json
from pathlib import Path
import select
import subprocess
import tempfile

from bundled_pari import bundled_pari_build

_processes = {}


class PariError(Exception):
    pass


def _call(kind, values, timeout=30):
    process = _processes.get(kind)
    if process is None:
        build = bundled_pari_build()
        obj = next(p.parent for p in build.glob('O*/pari.cfg'))
        source = Path(__file__).resolve().parents[1] / 'native' / 'pari_field_predicates.c'
        flags = ['-O3', '-fno-strict-aliasing']
        compiler = subprocess.check_output(['cc', '--version'])
        key = hashlib.sha256(source.read_bytes() + str(build).encode() + compiler
                             + repr(flags).encode()).hexdigest()[:16]
        folder = Path(tempfile.gettempdir()) / ('sage-bundled-pari_field_predicates-' + kind + '-' + key)
        folder.mkdir(exist_ok=True)
        executable = folder / 'oracle'
        if not executable.exists():
            library = next(p for p in obj.glob('libpari*') if p.suffix in ('.dylib', '.so'))
            subprocess.run(['cc', *flags, '-I' + str(obj), '-I' + str(build / 'src/headers'),
                            '-I' + str(build / 'src/basemath'), str(source), str(library),
                            '-Wl,-rpath,' + str(obj), '-o', str(executable)],
                           check=True, capture_output=True)
        process = subprocess.Popen([str(executable)], stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE, text=True, bufsize=1)
        _processes[kind] = process
        atexit.register(process.terminate)
    compact = lambda value: '[' + ','.join(map(str, value)) + ']' if isinstance(value, list) else str(value)
    process.stdin.write(' '.join(compact(value) for value in values) + '\n')
    process.stdin.flush()
    if not select.select([process.stdout], [], [], timeout)[0]:
        process.terminate()
        process.wait()
        del _processes[kind]
        raise RuntimeError('bundled PARI ' + kind + ' timed out; no comparison result')
    line = process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):
        # Match cypari2 rather than pari_err2str's added terminal punctuation.
        raise PariError(json.loads(line[6:]).removesuffix('.'))
    if not line or line == 'ERROR':
        raise RuntimeError('bundled PARI ' + kind + ' oracle failed: ' + line)
    return line


def pari_field_predicates(*args):
    return _call('pari_field_predicates', list(args))


def ff_extension_is_square(p,T,a,seed):
    from sage.all import GF,PolynomialRing
    R=PolynomialRing(GF(p),'x');F=GF(p**(len(T)-1),'a',modulus=R(T),impl='pari_ffelt')
    x=F(R(a));coeffs=[int(c)for c in x.polynomial().list()]
    result=pari_field_predicates(7,p,coeffs,T,seed)
    assert x.is_square()==json.loads(result)['value']
    return result


def pari_field_record_scalar(op,p,a,seed):
    if op==0:
        from pari_ff_square_root import pari_ff_square_root
        return pari_ff_square_root(p,[0,1],[a],seed)
    return pari_field_predicates(6 if op==1 else 7,p,[a],[0,1],seed)


def ff_extension_trace(p, T, a, seed):
    """Execute the bundled Sage trace body with an observed native PARI call."""
    import textwrap
    from types import SimpleNamespace
    from sage.all import GF, PolynomialRing, ZZ
    R = PolynomialRing(GF(p), 'x')
    F = GF(p**(len(T)-1), 'a', modulus=R(T), impl='pari_ffelt')
    x = F(R(a))
    source = (Path(__file__).resolve().parents[3] / 'reference/sage/src/sage/rings/finite_rings/element_base.pyx').read_text()
    body = source[source.index('    def trace(self):'):source.index('    def multiplicative_order(self):')]
    namespace = {}
    exec(compile(textwrap.dedent(body), 'element_base.pyx', 'exec'), namespace)
    calls = []
    native = None
    def trace():
        nonlocal native
        calls.append('FF_trace')
        native = json.loads(pari_field_predicates(8, p, list(map(int, x.polynomial().list())), T, seed))
        return SimpleNamespace(lift=lambda: ZZ(native['value']))
    proxy = SimpleNamespace(parent=lambda: F, __pari__=lambda: SimpleNamespace(trace=trace))
    result = namespace['trace'](proxy)
    assert result == x.trace()
    return json.dumps({'value': str(result), 'parent': str(result.parent()),
                       'identity': result.parent() is F.prime_subfield(),
                       'calls': calls, 'state': native['state']}, separators=(',', ':'))
