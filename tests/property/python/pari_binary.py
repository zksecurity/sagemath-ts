"""Compare packed polynomial/matrix adapters to the original PARI library."""
import atexit
import hashlib
from pathlib import Path
import subprocess
import tempfile
from bundled_pari import bundled_pari_build

_process=None


def native_pari_binary(op,a,b,c,d,e):
    global _process
    class RangeError(Exception):pass
    class PariError(Exception):pass
    if any(x<0 for x in (a,b,c)):raise RangeError('polynomial bits must be nonnegative')
    if op==2 and any(a>>i&1 for i in range(1,int(a).bit_length(),2)):raise RangeError('polynomial must be a square')
    if (op in (5,14) and not b) or (op==13 and not c) or (op==7 and not a) or (op==10 and d>=2 and not b):raise RangeError('polynomial divisor must be nonzero')
    if op==8 and a<2:raise RangeError('modulus must have positive degree')
    if op==10 and not 0<=d<=2**53-1:raise RangeError('power count must be nonnegative')
    if op==18 and not a:raise RangeError('polynomial divisor must be nonzero')
    if op==9:
        if not 0<=d<=2**53-1:raise RangeError('row count must be nonnegative')
        if a>>int(d*b):raise RangeError('column exceeds matrix row count')
    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_binary.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-binary-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*')if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),
                            '-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    _process.stdin.write(' '.join(map(str,(op,a,b,c,d,e)))+'\n');_process.stdin.flush()
    line=_process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):raise PariError(line[6:])
    if not line.startswith('OK '):raise RuntimeError('bundled PARI binary oracle failed: '+line)
    return line[3:]
