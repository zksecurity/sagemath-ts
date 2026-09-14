"""Bundled NTL probable-prime routines and deterministic byte streams compiled against Sage's native NTL primitives."""
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
    source = (root / 'tests/property/native/ntl_integer_reconstruction.cpp').read_text()
    originals = [root / 'reference/ntl/src/ZZ.cpp', root / 'reference/ntl/include/NTL/ZZ.h']
    original = originals[0].read_text()
    sha = original[original.index('#if (NTL_BITS_PER_INT32 == 32)'):original.index('// ******************** ChaCha20 stuff')]
    chacha = original[original.index('#define LE(p)'):original.index('old_RandomStream::old_RandomStream')]
    generic_start = original.index('struct RandomStream_impl {\n   _ntl_uint32 state[16];\n   unsigned char buf[64];')
    generic = original[generic_start:original.index('\n#endif\n#endif // defined(NTL_RANDOMSTREAM_AES256CTR)', generic_start)]
    pimpl = original[original.index('// Boilerplate PIMPL code'):original.index('NTL_TLS_GLOBAL_DECL(UniquePtr<RandomStream>,  CurrentRandomStream);')]
    header = originals[1].read_text()
    interface = header[header.index('class RandomStream {'):header.index('// this is the number of bits we can pass through the set_nonce')]
    source = source.replace('/*__NATIVE_RANDOM_BODY__*/', sha + chacha + generic + pimpl)
    source = source.replace('/*__NATIVE_RANDOM_CLASS__*/', interface)
    sampling = original[original.index('static inline\nunsigned long WordFromBytes'):original.index('// More prime generation stuff...')]
    source = source.replace('/*__NATIVE_SAMPLING_BODY__*/', sampling)
    source = source.replace('/*__NATIVE_SEED_BODY__*/', original[original.index('void SetSeed(const unsigned char *data, long dlen)'):original.index('static\nvoid InitRandomStream()')])
    prime_sequence = original[original.index('static Lazy< Vec<char> > lowsieve_storage;'):original.index('\nlong Jacobi(const ZZ& aa, const ZZ& nn)')]
    source = source.replace('/*__NATIVE_PRIME_SEQUENCE__*/', prime_sequence)
    primality = original[original.index('long PowerMod(long a, long ee, long n)'):original.index('void NextPrime(ZZ& n, const ZZ& m, long NumTrials)')]
    primality += original[original.index('long RandomPrime_long(long l, long NumTrials)'):original.index('static Lazy< Vec<char> > lowsieve_storage;')]
    primality += original[original.index('static\ndouble Log2(double x)'):original.index('void MultiThreadedGenGermainPrime(ZZ& n, long k, long err)')]
    source = source.replace('/*__NATIVE_PRIME_BODY__*/', primality)
    fft_source = root / 'reference/ntl/src/FFT.cpp'
    fft_header = root / 'reference/ntl/include/NTL/FFT.h'
    body = fft_source.read_text()
    body = body[body.index('long IsFFTPrime(long n, long& w)'):body.index('#ifdef NTL_FFT_LAZYMUL')]
    source = source.replace('/*__NATIVE_FFT_PRIMES__*/', body)
    originals += [fft_source, fft_header]
    reconstruction_sources = [root / 'reference/ntl/src' / n for n in ['mat_ZZ.cpp','mat_lzz_p.cpp','mat_ZZ_p.cpp','lzz_p.cpp','ZZ_p.cpp']]
    mat,word,big,ctx,bigctx = [p.read_text() for p in reconstruction_sources]
    cut = lambda s,start,end:s[s.index(start):s.index(end,s.index(start))]
    body = cut(original,'long CRT(ZZ& gg, ZZ& a, long G, long p)','void sub(ZZ& x, long a, const ZZ& b)')
    body += cut(big,'void determinant(ZZ_p& d, const mat_ZZ_p& M_in)','long IsIdent(const mat_ZZ_p& A, long n)')
    body += cut(ctx,'zz_pInfoT::zz_pInfoT(long NewP, long maxroot)','zz_pInfoT::zz_pInfoT(INIT_FFT_TYPE')
    body += cut(ctx,'void zz_p::init(long p, long maxroot)','void zz_p::UserFFTInit')
    body += cut(ctx,'zz_pContext::zz_pContext(long p, long maxroot)','zz_pContext::zz_pContext(INIT_USER_FFT_TYPE')
    body += cut(bigctx,'ZZ_pInfoT::ZZ_pInfoT(const ZZ& NewP)','// we use a lazy strategy')
    body += cut(bigctx,'void ZZ_p::init(const ZZ& p)','void ZZ_pContext::save()')
    source = source.replace('/*__NATIVE_RECONSTRUCTION_BODY__*/',body)
    originals += reconstruction_sources
    native_headers = ['ZZ_p.h','lzz_p.h','mat_ZZ.h','mat_ZZ_p.h','mat_lzz_p.h','matrix.h','vec_ZZVec.h']
    originals += [root / 'reference/ntl/include/NTL' / n for n in native_headers]
    local = Path(SAGE_LOCAL)
    for name,start in [('ZZ_p.h','class ZZ_pInfoT {'),('lzz_p.h','class zz_pInfoT {')]:
        if cut((root/'reference/ntl/include/NTL'/name).read_text(),start,'extern ') != cut((local/'include/NTL'/name).read_text(),start,'extern '):
            raise RuntimeError('native integer reconstruction context layout differs from bundled source')
    if fft_header.read_text() != (local / 'include/NTL/FFT.h').read_text():
        raise RuntimeError('native FFT-prime cache interface differs from bundled source')
    installed_header = (local / 'include/NTL/ZZ.h').read_text()
    installed_interface = installed_header[installed_header.index('class RandomStream {'):installed_header.index('// this is the number of bits we can pass through the set_nonce')]
    cached = lambda h: h[h.index('struct RandomBndGenerator {'):h.index('void VectorRandomWord(long k, unsigned long* x);')]
    prime = lambda h: h[h.index('class PrimeSeq {'):h.index('\n};',h.index('class PrimeSeq {'))+3]
    if interface != installed_interface or cached(header) != cached(installed_header) or prime(header) != prime(installed_header):
        raise RuntimeError('native prime stream, sampler or prime-sequence interface differs from bundled source')
    compiler = subprocess.check_output(['c++', '--version'])
    flags = ['-std=c++17', '-O2']
    version = b''.join((local / 'include/NTL' / name).read_bytes() for name in ['version.h', 'config.h', 'mach_desc.h', 'gmp_aux.h', *native_headers])
    key = hashlib.sha256(source.encode() + b"".join(p.read_bytes() for p in originals) + compiler + version + str(local).encode()
                         + repr(flags).encode()).hexdigest()[:16]
    folder = Path(tempfile.gettempdir()) / ('sage-bundled-ntl-integer-reconstruction-' + key)
    folder.mkdir(exist_ok=True)
    executable = folder / 'oracle'
    if not executable.exists():
        path = folder / 'oracle.cpp'
        path.write_text(source)
        subprocess.run(['c++', *flags, '-I' + str(local / 'include'), '-I' + str(root / 'reference/ntl/include'), str(path), str(reconstruction_sources[0]), str(reconstruction_sources[1]),
                        '-L' + str(local / 'lib'), '-lntl', '-lgmp',
                        '-Wl,-rpath,' + str(local / 'lib'), '-o', str(executable)],
                       check=True, capture_output=True)
    return executable


def ntl_integer_reconstruction(*args):
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
        raise RuntimeError('bundled NTL FFT-prime cache oracle timed out; no comparison result')
    line = _process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):
        raise Error(json.loads(line[6:]))
    if not line:
        _process.wait()
        _process = None
        raise RuntimeError('bundled NTL FFT-prime cache oracle failed without a result')
    return json.dumps(json.loads(line), separators=(',', ':'))
