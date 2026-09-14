"""Bundled NTL word transposed projection modules compiled against Sage's native NTL primitives."""
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
    source = (root / 'tests/property/native/ntl_word_projection.cpp').read_text()
    originals = [root / ('reference/ntl/src/' + name) for name in ['lzz_pX.cpp', 'lzz_pX1.cpp', 'lzz_pXFactoring.cpp', 'mat_lzz_p.cpp']]
    local = Path(SAGE_LOCAL)
    compiler = subprocess.check_output(['c++', '--version'])
    flags = ['-std=c++17', '-O2']
    profile = b"".join((local / 'include/NTL' / name).read_bytes()
                       for name in ['version.h', 'config.h', 'mach_desc.h'])
    key = hashlib.sha256(source.encode() + b"".join(p.read_bytes() for p in originals) + compiler + profile + str(local).encode()
                         + repr(flags).encode()).hexdigest()[:16]
    folder = Path(tempfile.gettempdir()) / ('sage-bundled-ntl-word-projection-' + key)
    folder.mkdir(exist_ok=True)
    executable = folder / 'oracle'
    if not executable.exists():
        path = folder / 'oracle.cpp'
        path.write_text(source)
        subprocess.run(['c++', *flags, '-I' + str(local / 'include'), '-I' + str(root / 'reference/ntl/include'), str(path), *map(str, originals),
                        '-L' + str(local / 'lib'), '-lntl', '-lgmp',
                        '-Wl,-rpath,' + str(local / 'lib'), '-o', str(executable)],
                       check=True, capture_output=True)
    # Revalidate the undefined native boundary whenever sources/compiler/profile
    # change, rather than silently keeping a guard after an upstream repair.
    probe_file = folder / 'prepared-zero-probe-v1.json'
    if not probe_file.exists():
        guard = ' if(param[1]==0 && H.mat.NumRows()>0 && a.length()<=F.n)throw std::runtime_error("ProjectPowers: prepared argument requires a positive count");\n'
        if source.count(guard) != 1:
            raise RuntimeError('native prepared-zero guard marker changed')
        unguarded = source.replace(guard, '')
        path = folder / 'unguarded.cpp'
        path.write_text(unguarded)
        probe = folder / 'unguarded'
        subprocess.run(['c++', *flags, '-I' + str(local / 'include'), '-I' + str(root / 'reference/ntl/include'), str(path), *map(str, originals),
                        '-L' + str(local / 'lib'), '-lntl', '-lgmp',
                        '-Wl,-rpath,' + str(local / 'lib'), '-o', str(probe)],
                       check=True, capture_output=True)
        observations = []
        for count in [1, 0]:
            data = f'4 5 [2 {count} 1 1] [5 3 1 0 1 0 2 1 2 2 0 1 0]\n'
            result = subprocess.run([str(probe)], input=data, text=True,
                                    capture_output=True, timeout=3)
            observations.append({'input': data, 'returncode': result.returncode,
                                 'stdout': result.stdout, 'stderr': result.stderr})
        if (observations[0]['returncode'] != 0 or
                json.loads(observations[0]['stdout']) != ['1'] or
                observations[1]['returncode'] != -11):
            raise RuntimeError('native prepared-zero behavior changed; review the explicit adapter guard')
        probe_file.write_text(json.dumps({'unguarded_wrapper_sha256': hashlib.sha256(unguarded.encode()).hexdigest(),
                                         'observations': observations}, indent=2) + '\n')
    return executable


def ntl_word_projection(*args):
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
        raise RuntimeError('bundled NTL word projection oracle timed out; no comparison result')
    line = _process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):
        raise Error(json.loads(line[6:]))
    if not line:
        _process.wait()
        _process = None
        raise RuntimeError('bundled NTL word projection oracle failed without a result')
    return json.dumps(json.loads(line), separators=(',', ':'))
