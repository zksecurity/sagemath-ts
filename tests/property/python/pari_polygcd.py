"""Direct native polynomial GCD and Bezout comparisons."""
import atexit
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
from bundled_pari import bundled_pari_build
_process=None

def native_pari_polygcd(op,p,a,b):
    global _process
    class RangeError(Exception):pass
    class PariError(Exception):pass
    class PariInvError(Exception):pass
    if op==15 and not (b and b[0]):raise RangeError('gen_pow_i requires a nonzero exponent')
    if op in (3,4,5,10):
        if not 1<=p<2**64:raise RangeError('modulus must be a positive word integer')
        if any(not 0<=c<p for c in list(a)+list(b)):raise RangeError('word polynomial coefficients must be reduced')
    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_polygcd.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-polygcd-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*')if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),
                            '-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    compact=lambda v:'['+','.join(str(c)for c in v)+']'
    _process.stdin.write(f'{op} {p} {compact(a)} {compact(b)}\n');_process.stdin.flush()
    line=_process.stdout.readline().rstrip('\n')
    if op in (6,14) and (line.startswith('ERROR impossible inverse in FpXQ_inv:') or line.startswith('ERROR impossible inverse in Flxq_inv:')):
        # Existing quotient-kernel error adapter (DEVIATIONS.md, Extension Arithmetic and PARI Quotient Kernels).
        raise PariInvError('impossible inverse in FpXQ_inv')
    if line.startswith('ERROR '):raise PariError(json.loads(line[6:]))
    if not line.startswith('OK '):raise RuntimeError('bundled PARI polynomial GCD oracle failed: '+line)
    return line[3:]
