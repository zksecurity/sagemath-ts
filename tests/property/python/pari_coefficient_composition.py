"""Direct bundled native PARI extension comparisons."""
import atexit,json,hashlib,subprocess,tempfile
from pathlib import Path
from bundled_pari import bundled_pari_build
_process=None
class RangeError(Exception):pass
class PariError(Exception):pass

def call(values):
    global _process
    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_coefficient_composition.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-coefficient-composition-'+key);folder.mkdir(exist_ok=True)
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

def native_pari_coefficient_composition(mode,op,p,T,P,x,V):
    if mode==1 and not 1<=p<2**64:raise RangeError('modulus must be a positive word integer')
    if mode==3 and P[0]<0:raise RangeError('polynomial bits must be nonnegative')
    if mode<2:
        def trim(a):
            a=list(a)
            while a and a[-1]==0:a.pop()
            return a
        polys=[c if isinstance(c,int) else trim(c) for c in P]
        while polys and (polys[-1]==0 or polys[-1]==[]):polys.pop()
        table=list(map(trim,V))
        if op%10==0:
            # Native power preparation precedes even empty substitution. Return
            # the actual table so unsafe matrix-copy guards use its real shape.
            table=json.loads(call([mode,9+(10 if op>=10 else 0),p,T,P,x,V]))
        t=trim(T)
        for i,c in enumerate(polys):
            if isinstance(c,int) or not c:continue
            message=None;length=len(table)
            if not t:message='polynomial modulus must be nonzero'
            elif not length or length==1 and len(c)>1:message='power table is too short'
            elif mode==1:
                used=table[:length if len(c)<=length else length-1]
                if any(len(a)>len(t)-1 for a in used):message='power polynomial exceeds modulus degree'
            if message is None and p<1:message='modulus must be positive'
            if message:
                # Evaluate the safe prefix with the same table/cache first;
                # earlier native inverse errors must still win over this guard.
                call([mode,11 if op>=10 else 1,p,T,polys[:i],x,table])
                raise RangeError(message)
    return call([mode,op,p,T,P,x,V])
