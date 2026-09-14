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
        source=Path(__file__).resolve().parents[1]/'native/pari_extension_minpoly.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-extension-minpoly-'+key);folder.mkdir(exist_ok=True)
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

def native_pari_extension_minpoly(mode,code,p,seed,T,S,x):
    if seed<=0:raise PariError('domain error in setrand: n <= 0')
    guard(mode,p,T,[x,S])
    def canonical(P):
        P=[list(c)if isinstance(c,list)else c for c in P]
        for c in P:
            if isinstance(c,list):
                while c and c[-1]==0:c.pop()
        while P and (P[-1]==0 or P[-1]==[]):P.pop()
        return P
    n=len(canonical(S))-1
    if n<1:raise RangeError('minimal polynomial modulus must have positive degree')
    if len(canonical(x))>n:raise RangeError('power polynomial exceeds modulus degree')
    return call([mode,code,p,seed,T,S,x])
