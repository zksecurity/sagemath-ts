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
        native_source=(Path(__file__).resolve().parents[3]/'reference/pari/src/basemath/FpE.c').read_text()
        start=native_source.index('\nGEN\nFp_ellcard(')+1
        end=native_source.index('\n}\n',start)+3
        dispatch=native_source[start:end].replace('Fp_ellcard(', 'audit_Fp_ellcard(',1)
        static_cards=[]
        for filename,name in [('FpE.c','FpXQ_ellcardj'),('FlxqE.c','Flxq_ellcardj')]:
            body=(Path(__file__).resolve().parents[3]/'reference/pari/src/basemath'/filename).read_text()
            card_start=body.index('\nstatic GEN\n'+name+'(')+1
            card_end=body.index('\n}\n',card_start)+3
            static_cards.append(body[card_start:card_end])
        cards='\n'.join(static_cards)
        body=(Path(__file__).resolve().parents[3]/'reference/pari/src/basemath/FlxqE.c').read_text()
        shanks_parts=[]
        for name in ['closest_lift','_FlxqE_order_multiple','_FlxqE_order','Flxq_kronecker','Flxq_ellpoint','Flxq_ellcard_Shanks']:
            part_start=body.rfind('static ',0,body.index('\n'+name+'('))
            brace=body.index('{',part_start); depth=1; part_end=brace+1
            while depth:
                depth+=(body[part_end]=='{')-(body[part_end]=='}');part_end+=1
            shanks_parts.append(body[part_start:part_end])
        shanks='\n'.join(shanks_parts)
        key=hashlib.sha256(source.read_bytes()+dispatch.encode()+cards.encode()+shanks.encode()+str(build).encode()+subprocess.check_output(['cc','--version'])).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-pari-odd-elliptic-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            (folder/'pari_cardinality_dispatch.h').write_text(dispatch)
            (folder/'pari_constant_j.h').write_text(cards)
            (folder/'pari_extension_shanks.h').write_text(shanks)
            library=next(p for p in obj.glob('libpari*') if p.suffix in ('.dylib','.so'))
            subprocess.run(['cc','-O2','-I'+str(folder),'-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),'-Wl,-rpath,'+str(obj),'-o',str(executable)],check=True,capture_output=True)
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


def pari_cardinality_dispatch(p,cs,cm):
    return pari_fq_elliptic(p,[0,1],cs,[],0,[0],cm,0,0,34)

def pari_word_elliptic(p,a,P,Q,n,op):
    return pari_fq_elliptic(p,[0,1],a,P,Q,[0],n,op,0,35)

def pari_prime_card_state(p,cs,seed,direct=0):
    return pari_fq_elliptic(p,[0,1],cs,[],[],[0],seed,0,0,36)

def pari_extension_sqrt(p,T,z,seed,mode=0,degree=0):
    return pari_fq_elliptic(p,T,z,[],[],[degree],seed,mode,0,37)

def pari_word_extension_card(p,T,a4,a6,seed,ordinary=0):
    return pari_fq_elliptic(p,T,[a4,a6],[],[],[0],seed,ordinary,0,38)

def pari_constant_j_card(p,T,z,kind,seed,backend):
    result=pari_fq_elliptic(p,T,z,[],[],[0],seed,kind,backend,39)
    # Independently check the native twist formulas on small generated fields.
    if int(p)**(len(T)-1)<=1000:
        from sage.all import GF,PolynomialRing
        p=int(p);n=len(T)-1
        K=GF(p**n,'a',modulus=PolynomialRing(GF(p),'x')(T)) if n>1 else GF(p)
        x=K.gen() if n>1 else K(0)
        u=sum((K(c)*x**i for i,c in enumerate(z)),K(0)) or K(1)
        if kind==0:a,b=K(0),u
        elif kind==1:a,b=u,K(0)
        else:
            j=int(kind)%p
            while j==0 or j==1728%p:j=(j+1)%p
            g=K(j)/K((1728-j)%p);a=3*g*u*u;b=2*g*u**3
        count=1
        for t in K:
            rhs=t**3+a*t+b
            count+=1 if rhs==0 else 2 if rhs.is_square() else 0
        assert int(json.loads(result)['value'][0])==count, 'native constant-j count disagrees with independent enumeration'
    return result


def pari_extension_shanks(p,T,a4,a6,seed):
    return pari_fq_elliptic(p,T,[a4,a6],[],[],[0],seed,0,1,40)
