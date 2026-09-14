"""Bundled PARI finite-field scalar square-root comparisons."""
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
        source = Path(__file__).resolve().parents[1] / 'native' / 'pari_ff_square_root.c'
        flags = ['-O3', '-fno-strict-aliasing']
        compiler = subprocess.check_output(['cc', '--version'])
        key = hashlib.sha256(source.read_bytes() + str(build).encode() + compiler
                             + repr(flags).encode()).hexdigest()[:16]
        folder = Path(tempfile.gettempdir()) / ('sage-bundled-pari_ff_square_root-' + kind + '-' + key)
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


def pari_ff_square_root(*args):
    return _call('pari_ff_square_root', list(args))


def ff_extension_sqrt(p, T, a, seed, mode):
    # element_pari_ffelt.pyx:1032-1103: option guard, FF_issquareall, [x,-x].
    from sage.all import GF, PolynomialRing
    F = GF(p**(len(T)-1), 'a', modulus=PolynomialRing(GF(p), 'x')(T), impl='pari_ffelt')
    value = F(PolynomialRing(GF(p), 'x')(a))
    reduced = [int(c) for c in value.polynomial().list()]
    # Native option failures do not call PARI. Obtain its freshly seeded state
    # through zero's FF_issquareall, whose source preserves that state.
    if mode in (2,5,6):
        state=json.loads(pari_ff_square_root(p,T,[],seed))['state']
        options={'unknown': True} if mode==6 else {'extend': True, 'all': mode==5}
        try:
            value.sqrt(**options)
        except Exception as error:
            return json.dumps({'error':type(error).__name__, 'message':str(error),
                               'state':state},separators=(',',':'))
        raise AssertionError('native unsupported sqrt option succeeded')
    raw=json.loads(pari_ff_square_root(p,T,reduced,seed));root=raw['value']
    # Check the Sage wrapper's success/list/parent contract with real elements;
    # root choice comes from bundled PARI, rather than Sage's installed version.
    native_roots=value.sqrt(all=True, extend=False)
    expected=[] if root is None else [F(PolynomialRing(GF(p), 'x')(list(map(int,root))))]
    if root and p!=2:expected.append(-expected[0])
    assert len(native_roots)==len(expected) and set(native_roots)==set(expected)
    assert all(r.parent() is F for r in native_roots)
    if root is None and mode not in (3,4):
        try:
            value.sqrt()
        except Exception as error:
            return json.dumps({'error':type(error).__name__, 'message':str(error),
                               'state':raw['state']},separators=(',',':'))
        raise AssertionError('native nonsquare sqrt succeeded')
    if mode in (3,4):
        values=[] if root is None else [root]
        if root and p!=2: values.append([str((-int(c))%p) for c in root])
    else:values=root
    return json.dumps({'value':values,'state':raw['state']},separators=(',',':'))
