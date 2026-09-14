"""Build the complete bundled NTL once, without mixing its globals with another version.

The fixed non-threaded 64-bit profile matches the existing Sage-backed oracles.
The source tree is copied to a content-addressed temporary cache before configuring.
"""
from pathlib import Path
import fcntl
import hashlib
import json
import shutil
import subprocess
import tempfile
from sage.env import SAGE_LOCAL


def native_library():
    root = Path(__file__).resolve().parents[3]
    original = root / 'reference/ntl'
    local = Path(SAGE_LOCAL)
    compiler = subprocess.check_output(['c++', '--version'])
    flags = [
        'CXX=c++', 'CXXFLAGS=-O2 -std=c++17', 'SHARED=off', 'TUNE=generic',
        'NTL_THREADS=off', 'NTL_THREAD_BOOST=off', 'NTL_TLS_HACK=on',
        'NTL_GMP_LIP=on', 'NTL_GF2X_LIB=on', 'NTL_FFT_BIGTAB=on',
        'NTL_FFT_LAZYMUL=on', 'NTL_SPMM_ULL=on', 'NTL_AVOID_BRANCHING=on',
        'NTL_TBL_REM=on', 'NTL_CRT_ALTCODE=on', 'NTL_GF2X_ALTCODE1=on',
        'NTL_ENABLE_AVX_FFT=off', 'NTL_RANDOM_AES256CTR=off',
        'GMP_PREFIX=' + str(local), 'GF2X_PREFIX=' + str(local),
    ]
    digest = hashlib.sha256(b'complete-ntl-oracle-v1' + compiler + repr(flags).encode())
    for path in sorted(p for folder in ['src', 'include'] for p in (original / folder).rglob('*') if p.is_file()):
        digest.update(str(path.relative_to(original)).encode() + b'\0' + path.read_bytes())
    digest.update((local / 'include/gmp.h').read_bytes())
    digest.update((local / 'include/gf2x.h').read_bytes())
    folder = Path(tempfile.gettempdir()) / ('sage-bundled-ntl-library-' + digest.hexdigest()[:16])
    folder.mkdir(exist_ok=True)
    prefix = folder / 'install'
    with (folder / 'build.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not (folder / 'complete.json').exists():
            source = folder / 'source'
            if not source.exists():
                shutil.copytree(original, source)
            with (folder / 'build.log').open('w') as log:
                commands = [
                    ['./configure', *flags, 'PREFIX=' + str(prefix)],
                    ['make', '-j2'], ['make', 'check'], ['make', 'install'],
                ]
                for command in commands:
                    log.write(json.dumps(command) + '\n')
                    log.flush()
                    try:
                        subprocess.run(command, cwd=source / 'src', stdout=log, stderr=subprocess.STDOUT, check=True)
                    except subprocess.CalledProcessError as error:
                        raise RuntimeError('bundled NTL build/check failed; inspect ' + str(folder / 'build.log')) from error
            if not (prefix / 'lib/libntl.a').exists():
                raise RuntimeError('bundled NTL build did not produce its static library')
            (folder / 'complete.json').write_text(json.dumps({'flags': flags, 'compiler': compiler.decode(), 'source_digest': digest.hexdigest(), 'native_checks': 'make check passed'}) + '\n')
    return prefix
