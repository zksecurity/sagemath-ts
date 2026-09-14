"""Bundled NTL word root prerequisite modules compiled against Sage's native NTL primitives."""
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
    source = (root / 'tests/property/native/ntl_root_prerequisites.cpp').read_text()
    originals = [root / ('reference/ntl/src/' + name) for name in ['lzz_pX.cpp', 'lzz_pX1.cpp', 'lzz_pXFactoring.cpp', 'mat_lzz_p.cpp']]
    original = (root / 'reference/ntl/src/ZZ.cpp').read_text()
    sha = original[original.index('#if (NTL_BITS_PER_INT32 == 32)'):original.index('// ******************** ChaCha20 stuff')]
    chacha = original[original.index('#define LE(p)'):original.index('old_RandomStream::old_RandomStream')]
    generic_start = original.index('struct RandomStream_impl {\n   _ntl_uint32 state[16];\n   unsigned char buf[64];')
    generic = original[generic_start:original.index('\n#endif\n#endif // defined(NTL_RANDOMSTREAM_AES256CTR)', generic_start)]
    pimpl = original[original.index('// Boilerplate PIMPL code'):original.index('NTL_TLS_GLOBAL_DECL(UniquePtr<RandomStream>,  CurrentRandomStream);')]
    header = (root / 'reference/ntl/include/NTL/ZZ.h').read_text()
    interface = header[header.index('class RandomStream {'):header.index('// this is the number of bits we can pass through the set_nonce')]
    source = source.replace('/*__NATIVE_RANDOM_BODY__*/', sha + chacha + generic + pimpl)
    sampling = original[original.index('static inline\nunsigned long WordFromBytes'):original.index('// More prime generation stuff...')]
    source = source.replace('/*__NATIVE_SAMPLING_BODY__*/', sampling)
    local = Path(SAGE_LOCAL)
    installed_header = (local / 'include/NTL/ZZ.h').read_text()
    installed_interface = installed_header[installed_header.index('class RandomStream {'):installed_header.index('// this is the number of bits we can pass through the set_nonce')]
    cached_class = header[header.index('struct RandomBndGenerator {'):header.index('void VectorRandomWord(long k, unsigned long* x);')]
    installed_cached_class = installed_header[installed_header.index('struct RandomBndGenerator {'):installed_header.index('void VectorRandomWord(long k, unsigned long* x);')]
    if cached_class != installed_cached_class:
        raise RuntimeError('native cached sampler public interface differs from bundled source')
    if interface != installed_interface:
        raise RuntimeError('native RandomStream public interface differs from bundled source')
    compiler = subprocess.check_output(['c++', '--version'])
    flags = ['-std=c++17', '-O2']
    version = b''.join((local / 'include/NTL' / name).read_bytes() for name in ['version.h','config.h','mach_desc.h'])
    key = hashlib.sha256(source.encode() + b"".join(p.read_bytes() for p in originals) + (root / 'reference/ntl/src/ZZ.cpp').read_bytes() + (root / 'reference/ntl/include/NTL/ZZ.h').read_bytes() + compiler + version + str(local).encode()
                         + repr(flags).encode()).hexdigest()[:16]
    folder = Path(tempfile.gettempdir()) / ('sage-bundled-ntl-root-prerequisites-' + key)
    folder.mkdir(exist_ok=True)
    executable = folder / 'oracle'
    if not executable.exists():
        path = folder / 'oracle.cpp'
        path.write_text(source)
        subprocess.run(['c++', *flags, '-I' + str(local / 'include'), '-I' + str(root / 'reference/ntl/include'), str(path), *map(str, originals),
                        '-L' + str(local / 'lib'), '-lntl', '-lgmp',
                        '-Wl,-rpath,' + str(local / 'lib'), '-o', str(executable)],
                       check=True, capture_output=True)
    return executable


def ntl_root_prerequisites(*args):
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
        raise RuntimeError('bundled NTL root prerequisite oracle timed out; no comparison result')
    line = _process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):
        raise Error(json.loads(line[6:]))
    if not line:
        _process.wait()
        _process = None
        raise RuntimeError('bundled NTL root prerequisite oracle failed without a result')
    return json.dumps(json.loads(line), separators=(',', ':'))
