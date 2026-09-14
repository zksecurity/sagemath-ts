"""Direct native polynomial division and reciprocal comparisons."""
import atexit
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
from bundled_pari import bundled_pari_build
_process=None

def native_pari_polydiv(op,p,a,b):
    global _process
    class RangeError(Exception):pass
    class PariError(Exception):pass
    if op in (3,4,5):
        if not 1<=p<2**64:raise RangeError('modulus must be a positive word integer')
        if any(not 0<=c<p for c in list(a)+(list(b)if op in (3,4) else [])):raise RangeError('word polynomial coefficients must be reduced')
        length=max((i+1 for i,c in enumerate(a)if c),default=0)
        if op==4 and p==1 and not any(b) and length+3<(159 if p<=3037000493 else 89):raise RangeError('remainder divisor must be nonzero')
    if op==1 and p==1:
        al=max((i+1 for i,c in enumerate(a)if c),default=0)
        bl=max((i+1 for i,c in enumerate(b)if c),default=0)
        if (al>=bl and bl>1) or al-bl+3>=111:
            raise RangeError('remainder divisor must be nonzero')
    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_polydiv.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-polydiv-'+key);folder.mkdir(exist_ok=True)
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
    if line.startswith('ERROR '):raise PariError(json.loads(line[6:]))
    if not line.startswith('OK '):raise RuntimeError('bundled PARI polynomial division oracle failed: '+line)
    return line[3:]
