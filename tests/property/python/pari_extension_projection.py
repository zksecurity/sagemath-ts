"""Direct bundled native PARI extension comparisons."""
import atexit,json,hashlib,subprocess,tempfile
from pathlib import Path
from bundled_pari import bundled_pari_build
_process=None
class RangeError(Exception):pass
class PariError(Exception):pass

def guard(mode,p,T,polys):
    if mode==2:
        if T<2:raise RangeError('extension modulus must have positive degree')
        if any(c<0 for a in polys for c in a):raise RangeError('polynomial bits must be nonnegative')
    else:
        if p<2:raise RangeError('extension arithmetic requires a modulus greater than one')
        t=list(T)
        while t and t[-1]==0:t.pop()
        if len(t)<2:raise RangeError('extension modulus must have positive degree')
        if mode==1:
            if p>=2**64:raise RangeError('modulus must be a positive word integer')
            if any(c<0 or c>=p for c in T+[v for a in polys for f in a for v in f]):raise RangeError('word polynomial coefficients must be reduced')

def call(values):
    global _process
    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_extension_projection.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-extension-projection-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*')if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),
                            '-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    compact=lambda v:'['+','.join(compact(c)for c in v)+']' if isinstance(v,list) else str(v)
    _process.stdin.write(' '.join(compact(v)for v in values)+'\n');_process.stdin.flush()
    line=_process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):raise PariError(json.loads(line[6:]))
    if not line.startswith('OK '):raise RuntimeError('bundled PARI extension oracle failed: '+line)
    return line[3:]

def native_pari_extension_projection(mode,code,p,seed,n,T,a,b):
    op=code%10
    if seed<=0:raise PariError('domain error in setrand: n <= 0')
    if op==0:
        if n<0 or n>2**53-1:raise RangeError('length must be nonnegative')
        if mode==1 and (p<0 or p>=2**64):raise RangeError('modulus must be an unsigned word integer')
        t=list(T)
        while t and t[-1]==0:t.pop()
        if mode==1 and any(c<0 or c>=2**64 for c in t):raise RangeError('modulus coefficients must be unsigned words')
        d=len(t)-1
        if n and d<0:raise RangeError('coefficient length must be nonnegative')
        if n and d and p<=0:raise RangeError('limit must be a positive unsigned word'if mode==1 else 'limit must be positive')
    else:
        if op>=2 and (n<0 or n>2**53-1):raise RangeError('truncation length must be nonnegative')
        guard(mode,p,T,[a]if op==3 else[a,b])
    return call([mode,code,p,seed,n,T,a,b])
