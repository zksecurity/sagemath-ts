"""Live bundled PARI odd-extension elliptic point kernels; no saved outputs."""
import atexit
import hashlib
import json
from pathlib import Path
import select
import subprocess
import tempfile
from bundled_pari import bundled_pari_build
_process=None

def pari_fq_elliptic(p,T,a,P,Q,ch,n,ordinary,backend,op):
    global _process
    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_odd_elliptic.c'
        key=hashlib.sha256(source.read_bytes()+str(build).encode()+subprocess.check_output(['cc','--version'])).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-pari-odd-elliptic-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            library=next(p for p in obj.glob('libpari*') if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),'-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    _process.stdin.write(json.dumps([p,T,a,P,Q,ch,n,ordinary,backend,op],default=int)+'\n');_process.stdin.flush()
    if not select.select([_process.stdout],[],[],30)[0]:
        _process.terminate()
        try: _process.wait(timeout=2)
        except subprocess.TimeoutExpired: _process.kill();_process.wait()
        _process=None
        raise RuntimeError('live odd elliptic comparison timed out')
    line=_process.stdout.readline().strip()
    if line.startswith('ERROR '): return json.dumps(dict(error='PariError',message=json.loads(line[6:])),separators=(',',':'))
    if not line.startswith('OK '): raise RuntimeError('native odd elliptic oracle failed: '+line)
    def strings(x): return None if x is None else list(map(strings,x)) if isinstance(x,list) else str(x)
    return json.dumps(dict(value=strings(json.loads(line[3:]))),separators=(',',':'))


def pari_fq_curve(p,T,a,P,m,ch,n,ordinary,backend,op):
    return pari_fq_elliptic(p,T,a,P,m,ch,n,ordinary,backend,op+10)


def pari_fq_model(p,T,cs,P,n,mode,encoding):
    return pari_fq_elliptic(p,T,cs,P,[],[encoding],n,mode,0,20)


def ec_fq_scalar(p,T,cs,point,n,known,encoding,family,mode):
    # Run the bundled caller, retaining native Sage parents and PARI conversion.
    import ast
    from sage.all import GF, Zmod, PolynomialRing, Integer, ZZ, pari, Sequence
    from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_finite_field
    from sage.structure.coerce_actions import IntegerMulAction
    from cypari2.handle_error import PariError
    from curve_pari_model import _method, _source
    from sage.all import EllipticCurve
    degree=len(T)-1
    def digits(v):
        out=[]
        while v: out.append(int(v%p));v//=p
        return out
    K=GF(p**degree,'a',modulus=PolynomialRing(GF(p),'t')(T),impl='pari_ffelt') if degree>1 else Zmod(p) if family==2 else GF(p)
    decode=lambda v:K(digits(v)) if degree>1 else K(v)
    pack=lambda v:sum(int(c)*int(p)**i for i,c in enumerate(v.polynomial().list())) if degree>1 else int(v)
    coefficients=list(map(decode,cs)); target=list(map(decode,point))
    if mode==3: coefficients[1]=-coefficients[0]**2/K(4)
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


def pari_generic_order(modulus,a,order,encoding,factors):
    return pari_fq_elliptic(modulus,[0,1],a,[factors[i:i+2] for i in range(0,len(factors),2)],0,[encoding],order,0,0,30)


def pari_prime_order_bound(p,cs,P,n):
    return pari_fq_elliptic(p,[0,1],cs,P,0,[0],n,0,0,31)


def pari_extension_trace(t,n,q):
    return pari_fq_elliptic(q,[0,1],t,[],0,[0],n,0,0,32)

def pari_extension_card(p,cs,n):
    return pari_fq_elliptic(p,[0,1],cs,[],0,[0],n,0,0,33)
