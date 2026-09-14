"""Compare original PARI DPE state with the bundled 64-bit native library.

Resource-boundary observations are intentionally a separate comparison dispatcher.
"""
import atexit
import json
import hashlib
from pathlib import Path
import select
import subprocess
import tempfile

from bundled_pari import bundled_pari_build

_process = None


class PariError(Exception):
    pass


def native_pari_lll_dpe(mode, m, n, dn, dd, en, ed, keep, track, want, flat):
    global _process
    if _process is None:
        build = bundled_pari_build()
        obj = next(p.parent for p in build.glob('O*/pari.cfg'))
        source = Path(__file__).resolve().parents[1] / 'native/pari_lll_dpe.c'
        flags = ['-O3', '-fno-strict-aliasing']
        compiler = subprocess.check_output(['cc', '--version'])
        key = hashlib.sha256(source.read_bytes() + (build / 'src/basemath/lll.c').read_bytes()
                             + str(build).encode() + compiler + repr(flags).encode()).hexdigest()[:16]
        folder = Path(tempfile.gettempdir()) / ('sage-bundled-pari-lll-dpe-' + key)
        folder.mkdir(exist_ok=True)
        executable = folder / 'oracle'
        if not executable.exists():
            library = next(p for p in obj.glob('libpari*') if p.suffix in ('.dylib', '.so'))
            subprocess.run(['cc', *flags, '-I' + str(obj), '-I' + str(build / 'src/headers'),
                            '-I' + str(build / 'src/basemath'), str(source), str(library),
                            '-Wl,-rpath,' + str(obj), '-o', str(executable)],
                           check=True, capture_output=True)
        _process = subprocess.Popen([str(executable)], stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE, text=True, bufsize=1)
        atexit.register(_process.terminate)
    values = [mode, m, n, dn, dd, en, ed, keep, track, want, '[' + ','.join(map(str, flat)) + ']']
    _process.stdin.write(' '.join(map(str, values)) + '\n')
    _process.stdin.flush()
    if not select.select([_process.stdout], [], [], 30)[0]:
        _process.terminate()
        _process.wait()
        _process = None
        raise RuntimeError('bundled PARI DPE LLL stage timed out; no comparison result')
    line = _process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):
        # cypari2 omits the final punctuation added by pari_err2str.
        raise PariError(json.loads(line[6:]).removesuffix('.'))
    if not line:
        raise RuntimeError('bundled PARI DPE LLL oracle failed: ' + line)
    return line


def native_pari_lll_dpe_resource(*args):
    """Separate resource-boundary comparison; do not count this as exact error parity."""
    try:
        native_pari_lll_dpe(*args)
    except PariError as error:
        if str(error).startswith('the PARI stack overflows !'):
            return 'resource_failure'
        raise
    return 'unexpected_success'
