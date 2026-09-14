"""Persistent original-PARI process; isolated from Sage's older linked PARI."""
import atexit
import hashlib
from pathlib import Path
import subprocess
import tempfile
from bundled_pari import bundled_pari_build

_process=None


def native_pari_rng(op,seed,limit,count,length=7):
    global _process
    class RangeError(Exception):pass
    class PariError(Exception):pass
    if op==1 and count>0 and not 1<=limit<=64:raise RangeError('bits must be between 1 and 64')
    if (op==2 or op==6 and length>0) and count>0 and not 1<=limit<2**64:raise RangeError('limit must be a positive unsigned word')
    if (op==3 or op==7 and length>0) and count>0 and limit<=0:raise RangeError('limit must be positive')
    if op in (4,5) and count>0 and not 0<=limit<=2**53-1:raise RangeError('length must be nonnegative')
    if op in (6,7) and count>0 and not 0<=length<=2**53-1:raise RangeError('length must be nonnegative')
    if _process is None:
        build=bundled_pari_build()
        obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_random.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-rng-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*')if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),
                            '-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    _process.stdin.write(' '.join(map(str,(op,seed,limit,count,length)))+'\n');_process.stdin.flush()
    line=_process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):raise PariError(line[6:])
    if not line.startswith('OK '):raise RuntimeError('bundled PARI RNG oracle failed: '+line)
    return line[3:]
