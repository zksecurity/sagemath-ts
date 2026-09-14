"""Compare the original PARI fast stage using the audited AArch64 Clang FP profile.

The native stage can loop on dependent keep-first inputs. Such probes live in the
separate timeout inventory, not in the ordinary return-value comparisons.
"""
import atexit
import hashlib
from pathlib import Path
import select
import subprocess
import tempfile

from bundled_pari import bundled_pari_build

_process = None


def native_pari_lll_fast(op, m, n, dn, dd, en, ed, keep, track, flat):
    global _process
    if _process is None:
        build = bundled_pari_build()
        obj = next(p.parent for p in build.glob('O*/pari.cfg'))
        source = Path(__file__).resolve().parents[1] / 'native/pari_lll_fast.c'
        flags = ['-O3', '-fno-strict-aliasing']
        compiler = subprocess.check_output(['cc', '--version'])
        key = hashlib.sha256(source.read_bytes() + (build / 'src/basemath/lll.c').read_bytes()
                             + str(build).encode() + compiler + repr(flags).encode()).hexdigest()[:16]
        folder = Path(tempfile.gettempdir()) / ('sage-bundled-pari-lll-fast-' + key)
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
    values = [op, m, n, dn, dd, en, ed, keep, track, '[' + ','.join(map(str, flat)) + ']']
    _process.stdin.write(' '.join(map(str, values)) + '\n')
    _process.stdin.flush()
    if not select.select([_process.stdout], [], [], 30)[0]:
        _process.terminate()
        _process.wait()
        _process = None
        raise RuntimeError('bundled PARI fast LLL stage timed out; no comparison result')
    line = _process.stdout.readline().rstrip('\n')
    if not line or line.startswith('ERROR'):
        raise RuntimeError('bundled PARI fast LLL oracle failed: ' + line)
    return line
