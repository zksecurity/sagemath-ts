"""Direct native modular square roots comparisons."""
import atexit
import hashlib
from pathlib import Path
import subprocess
import tempfile
from bundled_pari import bundled_pari_build
_process=None

def native_pari_sqrt(op,p,a,y):
    global _process
    class RangeError(Exception):pass
    class PariError(Exception):pass
    if op==3 and y==0:raise RangeError('gen_pow_fold requires a nonzero exponent')
    if op==2:
        if not 1<=p<2**64:raise RangeError('modulus must be a positive word integer')
        if not 0<=a<p:raise RangeError('word square-root argument must be reduced')
    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_sqrt.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-sqrt-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*')if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),
                            '-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    compact=lambda v:'['+','.join(str(c)for c in v)+']'
    _process.stdin.write(f'{op} {p} {a} {y}\n');_process.stdin.flush()
    line=_process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):raise PariError(line[6:])
    if not line.startswith('OK '):raise RuntimeError('bundled PARI modular square-root oracle failed: '+line)
    return line[3:]
