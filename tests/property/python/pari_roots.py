"""Direct native finite polynomial roots/counts comparisons."""
import atexit
import hashlib
from pathlib import Path
import subprocess
import tempfile
from bundled_pari import bundled_pari_build
_process=None

def native_pari_roots(op,p,a):
    global _process
    class RangeError(Exception):pass
    class PariError(Exception):pass
    if op in (1,4,5):
        if not 1<=p<2**64:raise RangeError('modulus must be a positive word integer')
        if any(not 0<=c<p for c in a):raise RangeError('word polynomial coefficients must be reduced')
    if op in (0,1) and p<2:raise RangeError('polynomial roots require a modulus greater than one')
    if op==2 and p<2:
        degree=len(a)-1
        while degree>=0 and a[degree]==0:degree-=1
        if degree>=2 and degree<=abs(p):raise RangeError('total splitting requires a modulus greater than one')

    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_roots.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-roots-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*')if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),
                            '-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    compact=lambda v:'['+','.join(str(c)for c in v)+']'
    _process.stdin.write(f'{op} {p} {compact(a)}\n');_process.stdin.flush()
    line=_process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):raise PariError(line[6:])
    if not line.startswith('OK '):raise RuntimeError('bundled PARI polynomial roots oracle failed: '+line)
    return line[3:]
