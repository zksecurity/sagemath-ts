"""Native PARI coordinates plus bundled Sage scalar caller/dependency checks."""
import atexit
import ast
import hashlib
import json
from pathlib import Path
import select
import subprocess
import tempfile
from bundled_pari import bundled_pari_build
from sage.all import EllipticCurve, GF, Integer, ZZ, pari, Sequence
from sage.structure.coerce_actions import IntegerMulAction
from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_finite_field
from cypari2.handle_error import PariError

_process=None

def pari_elliptic_model(p,cs,P,ch,n,op):
    global _process
    if _process is None:
        build=bundled_pari_build()
        obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_elliptic_model.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()+subprocess.check_output(['cc','--version'])).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-pari-elliptic-model-'+key)
        folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*') if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),'-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    _process.stdin.write(json.dumps([p,cs,P,ch,n,op],default=int)+'\n');_process.stdin.flush()
    if not select.select([_process.stdout],[],[],30)[0]:
        _process.terminate()
        try: _process.wait(timeout=2)
        except subprocess.TimeoutExpired: _process.kill(); _process.wait()
        _process=None
        raise RuntimeError('live PARI elliptic model comparison timed out')
    line=_process.stdout.readline().strip()
    if not line or line.startswith('ERROR'): raise RuntimeError(line)
    def strings(x): return list(map(strings,x)) if isinstance(x,list) else str(x)
    return json.dumps(strings(json.loads(line)),separators=(',',':'))

_source=Path(__file__).resolve().parents[3]/'reference/sage/src/sage/schemes/elliptic_curves/ell_point.py'
_tree=ast.parse(_source.read_text())
_method=next(n for c in _tree.body if isinstance(c,ast.ClassDef) and c.name=='EllipticCurvePoint_finite_field'
             for n in c.body if isinstance(n,ast.FunctionDef) and n.name=='_acted_upon_')

def ec_pari_scalar(p,cs,target,n,known,encoding=0,family=0):
    from sage.all import Zmod
    E=EllipticCurve(Zmod(p) if family==2 else GF(p),cs)
    P=E(target) if target else E(0)
    if known: P._order=P.order()
    else:
        try: del P._order
        except AttributeError: pass
    calls=[]
    class Pari:
        @staticmethod
        def ellmul(curve,point,k):
            calls.append(['ellmul',list(map(str,curve.a_invariants())),str(point),str(k)])
            return pari.ellmul(curve,point,k)
    namespace=dict(ZZ=ZZ,pari=Pari,Sequence=Sequence,PariError=PariError,IntegerMulAction=IntegerMulAction,
                   EllipticCurvePoint_finite_field=EllipticCurvePoint_finite_field)
    exec(compile(ast.Module(body=[_method],type_ignores=[]),str(_source),'exec'),namespace)
    Q=namespace['_acted_upon_'](P,float(n)+0.5 if encoding==2 else Integer(n),False)
    return json.dumps(dict(value=str(Q),order=str(Q._order) if hasattr(Q,'_order') else None,calls=calls,
                          model_cached=E.pari_curve() is E.pari_curve(),alias_cached=E.__pari__() is E.pari_curve()),separators=(',',':'))


def ec_pari_transformed(p,r,s,t,n,known,family):
    # Substitute old x=x+r, old y=y+s*x+t in y^2=x^3+x+1.
    cs=[2*s,3*r-s*s,2*t,3*r*r+1-2*s*t,1+r*r*r+r-t*t]
    P=[-r,1+s*r-t]
    return ec_pari_scalar(p,cs,P,n,known,0,family)

def pari_elliptic_transformed(p,r,s,t,n):
    cs=[2*s,3*r-s*s,2*t,3*r*r+1-2*s*t,1+r*r*r+r-t*t]
    return pari_elliptic_model(p,cs,[-r,1+s*r-t],[1,0,0,0],n,2)
