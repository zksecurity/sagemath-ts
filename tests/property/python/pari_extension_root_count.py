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
        source=Path(__file__).resolve().parents[1]/'native/pari_extension_root_count.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-extension-root-count-'+key);folder.mkdir(exist_ok=True)
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

def canonical(f):
    f=[list(c)if isinstance(c,list)else c for c in f]
    for c in f:
        if isinstance(c,list):
            while c and c[-1]==0:c.pop()
    while f and (f[-1]==0 or f[-1]==[]):f.pop()
    return f

def word_guard(p,f,T=()):
    if p<1 or p>=2**64:raise RangeError('modulus must be a positive word integer')
    if any(c<0 or c>=p for c in list(T)+[v for a in f for v in a]):raise RangeError('word polynomial coefficients must be reduced')

def native_pari_extension_root_count(mode,op,p,T,f):
    if mode==3:
        if any(isinstance(c,list)for c in f):raise TypeError('prime-field polynomial coefficients must be integers')
        return call([mode,op,p,T,f])
    if op==3:
        if mode==1:word_guard(p,f)
        return call([mode,op,p,T,f])
    if op==2:
        if mode==1:word_guard(p,f)
        elif p==0:call([mode,3,p,T,f]) # derivative errors precede GCD guards
        guard(mode,p,T,[f])
        return call([mode,op,p,T,f])
    raw=canonical(f)
    m=0 if mode==4 else mode
    if m==0 and 0<abs(p)<2**64:
        q=abs(p);t=[c%q for c in T]
        reduced=canonical([[v%q for v in c]if isinstance(c,list)else[c%q]for c in raw])
        if len(reduced)>2:guard(1,q,t,[reduced])
    else:
        if m==1:word_guard(p,raw,T)
        elif m==2 and (T<0 or any(c<0 for c in raw)):raise RangeError('polynomial bits must be nonnegative')
        if len(raw)>2:guard(m,p,T,[raw])
    return call([mode,op,p,T,f])
