"""Compare Hensel lifting and retained state with the complete bundled NTL."""
import atexit
import hashlib
import json
from pathlib import Path
import select
import subprocess
import tempfile
from sage.env import SAGE_LOCAL

_process = None


class Error(Exception):
    """NTL LogicErrorObject mapped to the port's existing JavaScript Error type."""


def _build():
    from ntl_native_library import native_library
    root = Path(__file__).resolve().parents[3]
    prefix = native_library()
    source = root / 'tests/property/native/ntl_hensel_state.cpp'
    import re
    original = (root / 'reference/ntl/src/ZZXFactoring.cpp').read_text()
    names = ['HenselLift', 'HenselLift1', 'BuildTree', 'RecTreeLift', 'TreeLift', 'MultiLift', 'AdditionalLifting']
    masked = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"', lambda m: " " * len(m.group()), original)
    bodies = []
    for name in names:
        start = re.search(r'void ' + name + r'\(', original).start()
        end = masked.index('{', start) + 1
        depth = 1
        while depth:
            depth += (masked[end] == '{') - (masked[end] == '}')
            end += 1
        bodies.append(original[start:end])
    generated = source.read_text().replace('// BUNDLED_NTL_HENSEL_ROUTINES', '\n'.join(bodies))
    digest = hashlib.sha256(generated.encode() + str(prefix).encode()).hexdigest()[:16]
    folder = Path(tempfile.gettempdir()) / ('sage-bundled-ntl-hensel-state-' + digest)
    folder.mkdir(exist_ok=True)
    executable = folder / 'oracle'
    if not executable.exists():
        source = folder / 'oracle.cpp'
        source.write_text(generated)
        local = Path(SAGE_LOCAL)
        subprocess.run(['c++', '-std=c++17', '-O2', '-I' + str(prefix / 'include'),
                        '-I' + str(root / 'reference/ntl/src'), str(source),
                        str(prefix / 'lib/libntl.a'), '-L' + str(local / 'lib'), '-lgmp', '-lgf2x',
                        '-Wl,-rpath,' + str(local / 'lib'), '-o', str(executable)], check=True, capture_output=True)
    return executable


def ntl_hensel_state(*args):
    global _process
    if _process is None:
        _process = subprocess.Popen([str(_build())], stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE, text=True, bufsize=1)
        atexit.register(_process.terminate)
    compact = lambda x: '[' + ' '.join(compact(y) for y in x) + ']' if isinstance(x, list) else str(x)
    _process.stdin.write(' '.join(map(compact, args)) + '\n')
    _process.stdin.flush()
    if not select.select([_process.stdout], [], [], 30)[0]:
        _process.terminate()
        _process.wait()
        _process = None
        raise RuntimeError('bundled NTL polynomial state oracle timed out; no comparison result')
    line = _process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):
        raise Error(json.loads(line[6:]))
    if not line:
        _process.wait()
        _process = None
        raise RuntimeError('bundled NTL polynomial state oracle failed without a result')
    return json.dumps(json.loads(line), separators=(',', ':'))
