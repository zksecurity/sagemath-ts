"""Native PARI composition, with explicit guards for unsafe C input contracts."""
import atexit
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
from bundled_pari import bundled_pari_build
_process=None

def native_pari_composition(op,p,n,Q,x,T):
    global _process
    class RangeError(Exception):pass
    class PariError(Exception):pass
    if op==4 and any(not -2**53+1<=v<=2**53-1 for v in (p,n,Q[0]if Q else 0)):raise RangeError('parameters must be safe integers')
    if op in (8,9) and Q:
        if not T:raise RangeError('polynomial modulus must be nonzero')
        raise RangeError('power table is too short')
    if op==5 and not 0<n<2**64:raise RangeError('gen_powu_i requires a nonzero unsigned word')
    if op in (1,3) and not 0<=n<=2**53-1:raise RangeError('power count must be nonnegative')
    if op==0:
        # The native entry converts word coefficients before its zero shortcut.
        word=0<abs(p)<2**64
        def canonical(v):
            v=[int(c)%abs(p) for c in v] if word else list(v)
            while v and v[-1]==0:v.pop()
            return v
        q,t,xx=canonical(Q),canonical(T),canonical(x)
        if q:
            if not t:raise RangeError('polynomial modulus must be nonzero')
            length=__import__('math').isqrt(len(q)-1)+1
            columns=length if len(q)<=length else length-1
            if word and (len(t)==1 or (columns>=2 and len(xx)>=len(t))):raise RangeError('power polynomial exceeds modulus degree')
    if op in (1,2,3) and Q:
        if not T:raise RangeError('polynomial modulus must be nonzero')
        count=int(n) if op in (1,3) else __import__('math').isqrt(len(Q)-1)
        length=count+1
        if length==1 and len(Q)>1:raise RangeError('power table is too short')
        columns=length if len(Q)<=length else length-1
        if op in (2,3)or(op==0 and 0<abs(p)<2**64):
            if len(T)==1 or(columns>=2 and len(x)>=len(T)):raise RangeError('power polynomial exceeds modulus degree')
    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_composition.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-composition-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*')if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),
                            '-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    compact=lambda v:json.dumps(list(map(int,v)),separators=(',',':'))
    _process.stdin.write(f'{op} {p} {n} {compact(Q)} {compact(x)} {compact(T)}\n');_process.stdin.flush()
    line=_process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):raise PariError(json.loads(line[6:]))
    if not line.startswith('OK '):raise RuntimeError('bundled PARI composition oracle failed: '+line)
    return line[3:]
