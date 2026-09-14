"""Bundled PARI Gram LLL, HNF transforms and Smith-group representatives."""
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


def _call(kind, values):
    process = _processes.get(kind)
    if process is None:
        build = bundled_pari_build()
        obj = next(p.parent for p in build.glob('O*/pari.cfg'))
        source = Path(__file__).resolve().parents[1] / 'native' / 'pari_lll_dependents.c'
        flags = ['-O3', '-fno-strict-aliasing']
        compiler = subprocess.check_output(['cc', '--version'])
        key = hashlib.sha256(source.read_bytes() + str(build).encode() + compiler
                             + repr(flags).encode()).hexdigest()[:16]
        folder = Path(tempfile.gettempdir()) / ('sage-bundled-pari-lll-dependents-' + kind + '-' + key)
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
    if not select.select([process.stdout], [], [], 30)[0]:
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


def pari_lll_dependents(op,n,flag,flat):
    if op==5 and not -(1<<63)<=n<(1<<63):
        import struct,math
        if not math.isfinite(struct.unpack('>d',struct.pack('>Q',int(flat[0])))[0]):
            raise PariError('overflow in dbltor [NaN or Infinity]')
        raise OverflowError('Python int too large to convert to C long')
    return _call('lll',[op,n,flag,flat])
