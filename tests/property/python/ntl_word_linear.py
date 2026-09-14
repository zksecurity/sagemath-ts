"""Bundled NTL word linear algebra modules compiled against Sage's native NTL primitives."""
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
    source = (root / 'tests/property/native/ntl_word_linear.cpp').read_text()
    originals = [root / 'reference/ntl/src/mat_lzz_p.cpp', root / 'reference/ntl/src/lzz_p.cpp', root / 'reference/ntl/include/NTL/lzz_p.h', root / 'reference/ntl/include/NTL/mat_lzz_p.h', root / 'reference/ntl/include/NTL/matrix.h']
    context_source = originals[1].read_text()
    body = context_source[context_source.index('zz_pInfoT::zz_pInfoT(long NewP, long maxroot)'):context_source.index('zz_pInfoT::zz_pInfoT(INIT_FFT_TYPE')]
    body += context_source[context_source.index('void zz_p::init(long p, long maxroot)'):context_source.index('void zz_p::FFTInit')]
    body += context_source[context_source.index('zz_pContext::zz_pContext(long p, long maxroot)'):context_source.index('zz_pContext::zz_pContext(INIT_FFT_TYPE')]
    source = source.replace('/*__NATIVE_CONTEXT_BODY__*/', body)
    local = Path(SAGE_LOCAL)
    context_class = lambda h: h[h.index('class zz_pInfoT {'):h.index('extern ')]
    if context_class(originals[2].read_text()) != context_class((local / 'include/NTL/lzz_p.h').read_text()):
        raise RuntimeError('native word linear-algebra context layout differs from bundled source')
    compiler = subprocess.check_output(['c++', '--version'])
    flags = ['-std=c++17', '-O2']
    version = b''.join((local / 'include/NTL' / name).read_bytes() for name in ['version.h','config.h','mach_desc.h'])
    key = hashlib.sha256(source.encode() + b"".join(p.read_bytes() for p in originals) + compiler + version + str(local).encode()
                         + repr(flags).encode()).hexdigest()[:16]
    folder = Path(tempfile.gettempdir()) / ('sage-bundled-ntl-word-linear-' + key)
    folder.mkdir(exist_ok=True)
    executable = folder / 'oracle'
    if not executable.exists():
        path = folder / 'oracle.cpp'
        path.write_text(source)
        subprocess.run(['c++', *flags, '-I' + str(local / 'include'), '-I' + str(root / 'reference/ntl/include'), str(path), str(originals[0]),
                        '-L' + str(local / 'lib'), '-lntl', '-lgmp',
                        '-Wl,-rpath,' + str(local / 'lib'), '-o', str(executable)],
                       check=True, capture_output=True)
    return executable


def ntl_word_linear(*args):
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
        raise RuntimeError('bundled NTL word linear-algebra oracle timed out; no comparison result')
    line = _process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):
        raise Error(json.loads(line[6:]))
    if not line:
        _process.wait()
        _process = None
        raise RuntimeError('bundled NTL word linear-algebra oracle failed without a result')
    return json.dumps(json.loads(line), separators=(',', ':'))
