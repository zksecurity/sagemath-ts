"""Bundled NTL factor selection bodies compiled against Sage's native NTL primitives."""
import atexit
import hashlib
import json
from pathlib import Path
import re
import select
import subprocess
import tempfile
from sage.env import SAGE_LOCAL

_process = None


class Error(Exception):
    """NTL LogicErrorObject mapped to the port's existing JavaScript Error type."""


def _build():
    root = Path(__file__).resolve().parents[3]
    original = (root / 'reference/ntl/src/ZZXFactoring.cpp').read_text()
    names = ['mul', 'inplace_rev', 'RecordPattern', 'NumFactors', 'CalcPossibleDegrees', 'ConstTermTest', 'BalCopy', 'InvMul', 'RemoveFactors', 'unpack', 'SubPattern']
    masked = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"', lambda m: " " * len(m.group()), original)
    bodies = []
    for name in names:
        for match in re.finditer(r'(?:void|long) ' + name + r'\(', original):
            start = match.start()
            if name == "mul" and original[start:].startswith("void mul(ZZX&"):
                continue
            end = masked.index('{', start) + 1
            depth = 1
            while depth:
                depth += (masked[end] == '{') - (masked[end] == '}')
                end += 1
            bodies.append((start, original[start:end]))
    bodies = [body for _, body in sorted(bodies)]
    template = (root / 'tests/property/native/ntl_factor_selection.cpp').read_text()
    source = template.replace('// BUNDLED_NTL_SELECTION_ROUTINES', '\n'.join(bodies))
    local = Path(SAGE_LOCAL)
    compiler = subprocess.check_output(['c++', '--version'])
    flags = ['-std=c++17', '-O2']
    version = (local / 'include/NTL/version.h').read_bytes()
    key = hashlib.sha256(source.encode() + compiler + version + str(local).encode()
                         + repr(flags).encode()).hexdigest()[:16]
    folder = Path(tempfile.gettempdir()) / ('sage-bundled-ntl-selection-' + key)
    folder.mkdir(exist_ok=True)
    executable = folder / 'oracle'
    if not executable.exists():
        path = folder / 'oracle.cpp'
        path.write_text(source)
        subprocess.run(['c++', *flags, '-I' + str(local / 'include'), str(path),
                        '-L' + str(local / 'lib'), '-lntl', '-lgmp',
                        '-Wl,-rpath,' + str(local / 'lib'), '-o', str(executable)],
                       check=True, capture_output=True)
    return executable


def ntl_factor_selection(*args):
    global _process
    if _process is None:
        _process = subprocess.Popen([str(_build())], stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE, text=True, bufsize=1)
        atexit.register(_process.terminate)
    compact = lambda x: '[' + ' '.join(map(str, x)) + ']' if isinstance(x, list) else str(x)
    _process.stdin.write(' '.join(map(compact, args)) + '\n')
    _process.stdin.flush()
    if not select.select([_process.stdout], [], [], 30)[0]:
        _process.terminate()
        _process.wait()
        _process = None
        raise RuntimeError('bundled NTL selection oracle timed out; no comparison result')
    line = _process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):
        raise Error(json.loads(line[6:]))
    if not line:
        _process.wait()
        _process = None
        raise RuntimeError('bundled NTL selection oracle failed without a result')
    return json.dumps(json.loads(line), separators=(',', ':'))
