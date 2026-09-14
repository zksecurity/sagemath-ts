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
        source=Path(__file__).resolve().parents[1]/'native/pari_extension_aut.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-extension-aut-'+key);folder.mkdir(exist_ok=True)
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

def native_pari_extension_aut(mode,code,p,n,T,S,phi,B,a):
    from pari_extension_matrix import guard
    op=code%10;unsigned=mode==1 and op==1
    if n==0 or n<(1 if unsigned else -(2**63)) or n>=(2**64 if unsigned else 2**63):
        raise RangeError('automorphism exponent must be a nonzero '+('unsigned' if unsigned else 'signed')+' word')
    if op==1 and mode!=2:guard(mode,p,T,[B,a])
    else:
        guard(mode,p,T,[[phi]])
        guard(mode,p,T,[B,a if op else []])
    guard(mode,p,T,[S])
    return call([mode,code,p,n,T,S,phi,B,a])
