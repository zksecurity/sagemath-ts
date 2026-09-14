"""Direct native minimal-polynomial and random-state comparisons."""
import atexit
import hashlib
from pathlib import Path
import subprocess
import tempfile
from bundled_pari import bundled_pari_build
_process=None

def native_pari_minpoly(op,seed,p,a,b):
    global _process
    class RangeError(Exception):pass
    class PariError(Exception):pass
    if op==1:
        if not 1<=p<2**64:raise RangeError('modulus must be a positive word integer')
        if any(not 0<=c<p for c in list(a)+list(b)):raise RangeError('word polynomial coefficients must be reduced')
    word=op==1 or 0<p<2**64
    aa=[int(c)%int(p)for c in a]if word else list(a)
    bb=[int(c)%int(p)for c in b]if word else list(b)
    while aa and aa[-1]==0:aa.pop()
    while bb and bb[-1]==0:bb.pop()
    if p<2:raise RangeError('minimal polynomial requires a modulus greater than one')
    if len(bb)<2:raise RangeError('minimal polynomial modulus must have positive degree')
    if word and len(aa)>=len(bb):raise RangeError('power polynomial exceeds modulus degree')
    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_minpoly.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-minpoly-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*')if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),
                            '-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    compact=lambda v:'['+','.join(str(c)for c in v)+']'
    _process.stdin.write(f'{op} {seed} {p} {compact(a)} {compact(b)}\n');_process.stdin.flush()
    line=_process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):raise PariError(line[6:])
    if not line.startswith('OK '):raise RuntimeError('bundled PARI minimal-polynomial oracle failed: '+line)
    return line[3:]
