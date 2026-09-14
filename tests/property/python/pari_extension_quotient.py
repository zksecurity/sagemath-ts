"""Direct native extension-polynomial quotient comparisons."""
import atexit
import json
import hashlib
from pathlib import Path
import subprocess
import tempfile
from bundled_pari import bundled_pari_build
_process=None

def native_pari_extension_quotient(mode,op,p,n,T,a,b,S):
    global _process
    class RangeError(Exception):pass
    class PariError(Exception):pass
    if mode==2:
        if T<2:raise RangeError('extension modulus must have positive degree')
        if any(c<0 for c in a+b+S):raise RangeError('polynomial bits must be nonnegative')
    else:
        if p<2:raise RangeError('extension arithmetic requires a modulus greater than one')
        t=list(T)
        while t and t[-1]==0:t.pop()
        if len(t)<2:raise RangeError('extension modulus must have positive degree')
        if mode==1:
            if p>=2**64:raise RangeError('modulus must be a positive word integer')
            if any(c<0 or c>=p for c in T+[v for f in a+b+S for v in f]):raise RangeError('word polynomial coefficients must be reduced')

    operation=op%10
    if operation==6 and (n<0 or n>=0xffffffff):raise RangeError('power count must be a nonnegative array length')
    if operation==7 and (n<0 or n>=2**64):raise RangeError('exponent must be an unsigned word integer')
    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_extension_quotient.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-extension-quotient-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*')if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),
                            '-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    raw=list(S)
    def zero(c):return not any(c) if isinstance(c,list) else c==0
    while raw and zero(raw[-1]):raw.pop()
    cache=(op%100)//10
    cached=cache==1 or cache==2 and len(raw)+2>12
    if mode==0 and p<2**64 and operation==5 and n not in (-1,0,1) and cached and any(isinstance(c,list) for c in raw):
        native_pari_extension_quotient(mode,op-operation+8,p,0,T,a,b,S)
        raise RangeError('generic word powering cannot convert cached polynomial coefficients')
    compact=lambda v:'['+','.join(compact(c)for c in v)+']' if isinstance(v,list) else str(v)
    _process.stdin.write(f'{mode} {op} {p} {n} {compact(T)} {compact(a)} {compact(b)} {compact(S)}\n');_process.stdin.flush()
    line=_process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):raise PariError(json.loads(line[6:]))
    if not line.startswith('OK '):raise RuntimeError('bundled PARI extension polynomial quotient oracle failed: '+line)
    return line[3:]
