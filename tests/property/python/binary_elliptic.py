"""Live bundled PARI binary elliptic point kernels; no saved outputs."""
import atexit
import hashlib
import json
from pathlib import Path
import select
import subprocess
import tempfile
from bundled_pari import bundled_pari_build
_process=None

def pari_f2_elliptic(T,a,P,Q,n,ch,op):
    global _process
    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_binary_elliptic.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()+subprocess.check_output(['cc','--version'])).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-pari-binary-elliptic-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*') if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),'-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    _process.stdin.write(json.dumps([T,a,P,Q,n,ch,op],default=int)+'\n');_process.stdin.flush()
    if not select.select([_process.stdout],[],[],30)[0]:
        _process.terminate()
        try: _process.wait(timeout=2)
        except subprocess.TimeoutExpired: _process.kill();_process.wait()
        _process=None
        raise RuntimeError('live binary elliptic comparison timed out')
    line=_process.stdout.readline().strip()
    if line.startswith('ERROR '): return json.dumps(dict(error='PariError',message=json.loads(line[6:])),separators=(',',':'))
    if not line.startswith('OK '): raise RuntimeError('native binary elliptic oracle failed: '+line)
    def strings(x): return None if x is None else list(map(strings,x)) if isinstance(x,list) else str(x)
    return json.dumps(dict(value=strings(json.loads(line[3:]))),separators=(',',':'))

def pari_f2_curve(T,mode,x,y,a3,a4,n,m,op):
    return pari_f2_elliptic(T,[mode,a3,a4],[x,y],m,n,[],op+10)
