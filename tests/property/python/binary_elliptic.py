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


def pari_f2_model(T,cs,P,n,mode,encoding):
    return pari_f2_elliptic(T,cs,P,[],n,[encoding],20+mode)


def ec_f2_scalar(T,cs,point,n,known,encoding,family,mode):
    # Run the bundled caller, retaining native Sage parents and PARI conversion.
    import ast
    from sage.all import GF, Zmod, PolynomialRing, Integer, ZZ, pari, Sequence
    from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_finite_field
    from sage.structure.coerce_actions import IntegerMulAction
    from cypari2.handle_error import PariError
    from curve_pari_model import _method, _source
    from sage.all import EllipticCurve
    degree=int(T).bit_length()-1
    bits=lambda v:[(int(v)>>i)&1 for i in range(max(1,int(v).bit_length()))]
    K=GF(2**degree,'a',modulus=PolynomialRing(GF(2),'t')(bits(T)),impl='pari_ffelt') if degree>1 else Zmod(2) if family==2 else GF(2)
    decode=lambda v:K(bits(v)) if degree>1 else K(v)
    pack=lambda v:sum(int(c)<<i for i,c in enumerate(v.polynomial().list())) if degree>1 else int(v)
    coefficients=list(map(decode,cs)); target=list(map(decode,point))
    if mode==3: coefficients[0]=K.zero()
    if mode>=2:
        a1,a2,a3,a4,_=coefficients;x,y=target
        coefficients[4]=y*y+a1*x*y+a3*y-x**3-a2*x*x-a4*x
    try: E=EllipticCurve(K,coefficients)
    except ArithmeticError: return json.dumps(dict(singular=True),separators=(',',':'))
    P=E(target) if target else E(0)
    if known: P._order=P.order()
    else:
        try: del P._order
        except AttributeError: pass
    calls=[]
    def coordinates(Q): return [] if Q.is_zero() else list(map(lambda a:str(pack(a)),Q[:2]))
    class Pari:
        @staticmethod
        def ellmul(curve,point,k):
            calls.append(['ellmul',list(map(lambda a:str(pack(a)),curve.a_invariants())),coordinates(point),str(k)])
            return pari.ellmul(curve,point,k)
    namespace=dict(ZZ=ZZ,pari=Pari,Sequence=Sequence,PariError=PariError,IntegerMulAction=IntegerMulAction,
                   EllipticCurvePoint_finite_field=EllipticCurvePoint_finite_field)
    exec(compile(ast.Module(body=[_method],type_ignores=[]),str(_source),'exec'),namespace)
    try:
        Q=namespace['_acted_upon_'](P,float(n)+0.5 if encoding==2 else float(n) if encoding==3 else Integer(n),False)
        return json.dumps(dict(value=coordinates(Q),order=str(Q._order) if hasattr(Q,'_order') else None,calls=calls,
            model_cached=E.pari_curve() is E.pari_curve(),alias_cached=E.__pari__() is E.pari_curve(),
            same_curve=Q.curve() is E,same_field=all(c.parent() is K for c in Q)),separators=(',',':'))
    except Exception as e:
        return json.dumps(dict(error=type(e).__name__,message=str(e),calls=calls),separators=(',',':'))
