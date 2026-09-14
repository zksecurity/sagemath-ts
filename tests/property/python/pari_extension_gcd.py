"""Direct native extension-polynomial GCD comparisons."""
import atexit
import json
import hashlib
from pathlib import Path
import subprocess
import tempfile
from bundled_pari import bundled_pari_build
_process=None

def native_pari_extension_gcd(mode,op,p,T,a,b):
    global _process
    class RangeError(Exception):pass
    class PariError(Exception):pass
    if mode==2:
        if T<2:raise RangeError('extension modulus must have positive degree')
        if any(c<0 for c in a+b):raise RangeError('polynomial bits must be nonnegative')
    else:
        if p<2:raise RangeError('extension arithmetic requires a modulus greater than one')
        t=list(T)
        while t and t[-1]==0:t.pop()
        if len(t)<2:raise RangeError('extension modulus must have positive degree')
        if mode==1:
            if p>=2**64:raise RangeError('modulus must be a positive word integer')
            if any(c<0 or c>=p for c in T+[v for f in a+b for v in f]):raise RangeError('word polynomial coefficients must be reduced')

    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_extension_gcd.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-extension-gcd-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*')if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),
                            '-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    compact=lambda v:'['+','.join(compact(c)for c in v)+']' if isinstance(v,list) else str(v)
    _process.stdin.write(f'{mode} {op} {p} {compact(T)} {compact(a)} {compact(b)}\n');_process.stdin.flush()
    line=_process.stdout.readline().rstrip('\n')
    if line.startswith('ERROR '):raise PariError(json.loads(line[6:]))
    if not line.startswith('OK '):raise RuntimeError('bundled PARI extension polynomial GCD oracle failed: '+line)
    return line[3:]

def native_binary_extension_gcd_termination(degree,op):
    # Warm only the native executable before starting the operation deadline.
    if _process is None:native_pari_extension_gcd(2,0,2,7,[],[])
    a='['+','.join(map(str,[1]+[0]*(degree-1)+[1]))+']'
    try:
        r=subprocess.run(_process.args,input=f'2 {op} 2 7 {a} [0,0,1]\n',
                         text=True,capture_output=True,timeout=2,check=True)
    except subprocess.TimeoutExpired:return 'timeout'
    return r.stdout.strip()
