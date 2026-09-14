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
        source=Path(__file__).resolve().parents[1]/'native/pari_inner_cache.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-inner-cache-'+key);folder.mkdir(exist_ok=True)
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

def native_pari_inner_cache(op,p,n,Q,x,V,T):
    import math
    def trim(a):
        a=list(a)
        while a and a[-1]==0:a.pop()
        return a
    if op==1 and not 1<=p<2**64:raise RangeError('modulus must be a positive word integer')
    if op<2 and not 0<=n<=2**53-1:raise RangeError('power count must be nonnegative')
    if op>=2:
        q,t,v=trim(Q),trim(T),list(map(trim,V));xx=trim(x)
        word=op in (3,5)or op==4 and 0<p<2**64
        if op==4 and word:
            q=trim([c%p for c in q]);t=trim([c%p for c in t]);xx=trim([c%p for c in xx])
        if q:
            generated=op>=4
            count=math.isqrt(len(q)-1)if generated else 0
            length=count+1 if generated else len(v)
            columns=length if len(q)<=length else length-1
            message=None
            if not t:message='polynomial modulus must be nonzero'
            elif not length or length==1 and len(q)>1:message='power table is too short'
            elif word:
                used=([[1],xx]if columns>=2 else [[1]])if generated else v[:columns]
                if any(len(a)>len(t)-1 for a in used):message='power polynomial exceeds modulus degree'
            if message:
                if generated:call([1 if word else 0,p,count,[],xx,[],t])
                raise RangeError(message)
            if word and not 1<=p<2**64:raise RangeError('modulus must be a positive word integer')
            if not word and p<1:raise RangeError('modulus must be positive')
    return call([op,p,n,Q,x,V,T])
