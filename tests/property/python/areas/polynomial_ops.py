"""SageMath side of the ``polynomial_ops`` property-test area (ZZ[x] and GF(p)[x]).

Cases: tests/property/cases/polynomial_ops.cases.json
"""

from sage.all import *

from ._helpers import (
    poly_derivative,
    poly_eval,
    poly_factor,
    poly_gcd,
    poly_is_irreducible,
    poly_mul,
    poly_pow,
    poly_quo_rem,
    poly_roots,
)


def _make_poly_zz(coeffs):
    """Create a polynomial over ZZ from coefficient list [c0, c1, ..., cn]."""
    R = PolynomialRing(ZZ, 'x')
    return R(list(map(Integer, coeffs)))


def poly_gcd_zz(coeffs1, coeffs2):
    """Compute GCD of two polynomials over ZZ."""
    f = _make_poly_zz(coeffs1)
    g = _make_poly_zz(coeffs2)
    result = gcd(f, g)
    return [int(c) for c in result.list()]


def poly_eval_zz(coeffs, x):
    """Evaluate polynomial over ZZ at a point."""
    f = _make_poly_zz(coeffs)
    return int(f(Integer(x)))


FUNCTIONS = {
    'poly_gcd_zz': poly_gcd_zz,
    'poly_eval_zz': poly_eval_zz,
    'poly_gcd_ff': poly_gcd,
    'poly_eval_ff': poly_eval,
    'poly_factor_ff': poly_factor,
    'poly_roots_ff': poly_roots,
    'poly_derivative_ff': poly_derivative,
    'poly_is_irreducible_ff': poly_is_irreducible,
    'poly_mul_ff': poly_mul,
    'poly_quo_rem_ff': poly_quo_rem,
    'poly_pow_ff': poly_pow,
}

# Call FLINT itself, bypassing Sage's zero/one polynomial GCD shortcuts.
# Resolve symbols through Sage's loaded extension so no hardcoded library path
# or optional system FLINT installation is required.
import ctypes
import sage.rings.polynomial.polynomial_integer_dense_flint as _flint_extension
_flint = ctypes.CDLL(_flint_extension.__file__)
class _FmpzPoly(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.POINTER(ctypes.c_long)), ('alloc', ctypes.c_long), ('length', ctypes.c_long)]
_PolyPtr = ctypes.POINTER(_FmpzPoly)
for _name in ('fmpz_poly_init', 'fmpz_poly_clear'):
    getattr(_flint, _name).argtypes = [_PolyPtr]
    getattr(_flint, _name).restype = None
_flint.fmpz_poly_set_str.argtypes = [_PolyPtr, ctypes.c_char_p]
_flint.fmpz_poly_set_str.restype = ctypes.c_int
_flint.fmpz_poly_get_str.argtypes = [_PolyPtr]
_flint.fmpz_poly_get_str.restype = ctypes.c_void_p
_flint.flint_free.argtypes = [ctypes.c_void_p]
_flint.flint_free.restype = None

def flint_poly_gcd(a, b, method):
    polys = [_FmpzPoly() for _ in range(3)]
    for poly in polys: _flint.fmpz_poly_init(ctypes.byref(poly))
    try:
        for poly, coeffs in zip(polys, (a,b)):
            raw = (str(len(coeffs)) + '  ' + ' '.join(map(str,coeffs))).encode()
            if _flint.fmpz_poly_set_str(ctypes.byref(poly), raw): raise RuntimeError('invalid FLINT polynomial fixture')
        f = getattr(_flint, 'fmpz_poly_gcd' + method)
        f.argtypes = [_PolyPtr, _PolyPtr, _PolyPtr]
        f.restype = ctypes.c_int if method == '_heuristic' else None
        result = f(ctypes.byref(polys[2]), ctypes.byref(polys[0]), ctypes.byref(polys[1]))
        if method == '_heuristic' and not result: return None
        ptr = _flint.fmpz_poly_get_str(ctypes.byref(polys[2]))
        try: return [Integer(c) for c in ctypes.string_at(ptr).decode().split()[1:]]
        finally: _flint.flint_free(ptr)
    finally:
        for poly in polys: _flint.fmpz_poly_clear(ctypes.byref(poly))
for _method in ('', '_subresultant', '_heuristic', '_modular'):
    FUNCTIONS['flint_poly_gcd' + _method] = lambda a,b,method=_method: flint_poly_gcd(a,b,method)

def flint_nmod_poly_gcd(p, a, b):
    R = PolynomialRing(GF(p), 'x')
    g = R(a).gcd(R(b))
    # Sage returns the operand on zero shortcuts; FLINT's public wrapper is monic.
    return list(g.monic()) if g else []
FUNCTIONS['flint_nmod_poly_gcd'] = flint_nmod_poly_gcd

_flint.fmpz_poly_fit_length.argtypes = [_PolyPtr, ctypes.c_long]
_flint.fmpz_poly_fit_length.restype = None
_flint.fmpz_poly_divrem.argtypes = [_PolyPtr] * 4
_flint.fmpz_poly_divrem.restype = None
_flint._fmpz_poly_divrem.argtypes = [ctypes.POINTER(ctypes.c_long)] * 3 + [ctypes.c_long, ctypes.POINTER(ctypes.c_long), ctypes.c_long, ctypes.c_int]
_flint._fmpz_poly_divrem.restype = ctypes.c_int

def flint_poly_divrem(a, b, exact):
    polys = [_FmpzPoly() for _ in range(4)]
    for poly in polys: _flint.fmpz_poly_init(ctypes.byref(poly))
    try:
        for poly, coeffs in zip(polys, (a,b)):
            raw = (str(len(coeffs)) + '  ' + ' '.join(map(str,coeffs))).encode()
            if _flint.fmpz_poly_set_str(ctypes.byref(poly),raw): raise RuntimeError('invalid FLINT polynomial fixture')
        A,B,Q,R = polys
        if not B.length: raise ValueError('nonzero divisor required by FLINT fixture')
        if exact and A.length >= B.length:
            _flint.fmpz_poly_fit_length(ctypes.byref(Q), A.length-B.length+1)
            _flint.fmpz_poly_fit_length(ctypes.byref(R), A.length)
            if not _flint._fmpz_poly_divrem(Q.coeffs,R.coeffs,A.coeffs,A.length,B.coeffs,B.length,1): return None
            Q.length, R.length = A.length-B.length+1, A.length
        else:
            _flint.fmpz_poly_divrem(ctypes.byref(Q),ctypes.byref(R),ctypes.byref(A),ctypes.byref(B))
        result = []
        for poly in (Q,R):
            ptr = _flint.fmpz_poly_get_str(ctypes.byref(poly))
            try: values = [Integer(c) for c in ctypes.string_at(ptr).decode().split()[1:]]
            finally: _flint.flint_free(ptr)
            while values and not values[-1]: values.pop()
            result.append(values)
        return result
    finally:
        for poly in polys: _flint.fmpz_poly_clear(ctypes.byref(poly))
import json
def flint_poly_divrem_json(*args):
    result = flint_poly_divrem(*args)
    return json.dumps(None if result is None else [[str(c) for c in poly] for poly in result], separators=(',',':'))
FUNCTIONS['flint_poly_divrem'] = flint_poly_divrem_json

def flint_nmod_poly_divrem(p, a, b):
    R = PolynomialRing(Zmod(p), 'x')
    q, r = R(a).quo_rem(R(b))
    return json.dumps([[str(c) for c in poly.list()] for poly in (q,r)], separators=(',',':'))
FUNCTIONS['flint_nmod_poly_divrem'] = flint_nmod_poly_divrem

def flint_fmpq_poly_gcd(a, b):
    R = PolynomialRing(QQ, 'x')
    g = R(a).gcd(R(b))
    d = g.denominator()
    return json.dumps([[str(c*d) for c in g.list()], str(d)], separators=(',', ':'))
FUNCTIONS['flint_fmpq_poly_gcd'] = flint_fmpq_poly_gcd

_native_xgcd_module = None

def flint_nmod_poly_xgcd(n, a, b):
    global _native_xgcd_module
    import os, sys
    from pathlib import Path
    sys.stdout.flush()
    libc = ctypes.CDLL(None)
    libc.fflush(None)
    saved = os.dup(1)
    try:
        os.dup2(2, 1)
        if _native_xgcd_module is None:
            from sage.misc.cython import cython_import
            _native_xgcd_module = cython_import(str(Path(__file__).parent.parent / 'polynomial_xgcd_flint.pyx'))
        try:
            value = _native_xgcd_module.native_xgcd(a, b, n)
            result = {'value': [[str(c) for c in poly] for poly in value]}
        except RuntimeError:
            result = {'error': 'noninvertible'}
        return json.dumps(result, separators=(',', ':'))
    finally:
        sys.stdout.flush()
        libc.fflush(None)
        os.dup2(saved, 1)
        os.close(saved)
FUNCTIONS['flint_nmod_poly_xgcd'] = flint_nmod_poly_xgcd

_native_resultant_module = None

def native_resultant_comparison(a, b, n=0, same=0):
    global _native_resultant_module
    import os, sys
    from pathlib import Path
    sys.stdout.flush()
    libc = ctypes.CDLL(None)
    libc.fflush(None)
    saved = os.dup(1)
    try:
        os.dup2(2, 1)
        if _native_resultant_module is None:
            from sage.misc.cython import cython_import
            _native_resultant_module = cython_import(str(Path(__file__).parent.parent / 'polynomial_resultant_flint.pyx'))
        try:
            value = _native_resultant_module.native_modular_resultant(a,b,n,bool(same)) if n else _native_resultant_module.native_integer_resultant(a,b)
            result = {'value': str(value)}
        except RuntimeError:
            result = {'error': 'noninvertible'}
        return json.dumps(result, separators=(',', ':'))
    finally:
        sys.stdout.flush()
        libc.fflush(None)
        os.dup2(saved, 1)
        os.close(saved)
FUNCTIONS['flint_poly_resultant'] = lambda a,b: native_resultant_comparison(a,b)
FUNCTIONS['flint_nmod_poly_resultant'] = lambda n,a,b,same: native_resultant_comparison(a,b,n,same)

def flint_poly_xgcd(a, b):
    global _native_xgcd_module
    import os, sys
    from pathlib import Path
    sys.stdout.flush()
    saved = os.dup(1)
    try:
        os.dup2(2, 1)
        if _native_xgcd_module is None:
            from sage.misc.cython import cython_import
            _native_xgcd_module = cython_import(str(Path(__file__).parent.parent / 'polynomial_xgcd_flint.pyx'))
        return json.dumps(_native_xgcd_module.native_integer_xgcd(a,b), separators=(',', ':'))
    finally:
        sys.stdout.flush()
        os.dup2(saved, 1)
        os.close(saved)
FUNCTIONS['flint_poly_xgcd'] = flint_poly_xgcd

def flint_fmpq_poly_xgcd(a, den_a, b, den_b):
    R = PolynomialRing(QQ, 'x')
    result = R([QQ(c)/den_a for c in a]).xgcd(R([QQ(c)/den_b for c in b]))
    return json.dumps([[[str(c*f.denominator()) for c in f.list()], str(f.denominator())] for f in result], separators=(',', ':'))
FUNCTIONS['flint_fmpq_poly_xgcd'] = flint_fmpq_poly_xgcd

def flint_fmpq_poly_resultant(a, den_a, b, den_b):
    R = PolynomialRing(QQ, 'x')
    r = R([QQ(c)/den_a for c in a]).resultant(R([QQ(c)/den_b for c in b]))
    return json.dumps([str(r.numerator()),str(r.denominator())],separators=(',', ':'))
FUNCTIONS['flint_fmpq_poly_resultant'] = flint_fmpq_poly_resultant

def native_real_result(run):
    import struct
    try:
        value = float(run())
        return json.dumps({'value':str(int.from_bytes(struct.pack('>d',value),'big'))},separators=(',', ':'))
    except Exception as e:
        if 'overflow' in str(e).lower(): return json.dumps({'error':'overflow'},separators=(',', ':'))
        raise

def real_from_bits(bits):
    import struct
    return struct.unpack('>d',int(bits%(2**64)).to_bytes(8,'big'))[0]
FUNCTIONS['pari_double_roundtrip'] = lambda a: native_real_result(lambda: pari(real_from_bits(a)))
FUNCTIONS['pari_real_resultant'] = lambda a,b: native_real_result(lambda: pari([real_from_bits(x) for x in a]).Polrev('x').polresultant(pari([real_from_bits(x) for x in b]).Polrev('x'),'x'))

def pari_real_arithmetic(operation, a, b):
    x,y = pari(real_from_bits(a)),pari(real_from_bits(b))
    try:
        z = x+y if operation == 0 else x-y if operation == 1 else x*y if operation == 2 else x/y if operation == 3 else (x/y)+(x*y) if operation == 4 else (x+y)/(x-y)
        q = z.sage().exact_rational()
        return json.dumps([str(z.sign()),str(z.exponent()),str(q.numerator()),str(q.denominator())],separators=(',', ':'))
    except Exception as e:
        if 'inverse' in str(e): return json.dumps({'error':'inverse'},separators=(',', ':'))
        raise
FUNCTIONS['pari_real_arithmetic'] = pari_real_arithmetic

# QQ._derivative delegates to native fmpq_poly_derivative, which calls the
# original _fmpz_poly_derivative; denominator-one inputs exercise that kernel.
def flint_poly_derivative(a):
    g=PolynomialRing(QQ,'x')(a)._derivative()
    return json.dumps(list(map(str,g.list())),separators=(',', ':'))
def flint_fmpq_poly_derivative(a,den):
    g=PolynomialRing(QQ,'x')([QQ(c)/den for c in a])._derivative()
    return json.dumps([[str(c*g.denominator()) for c in g.list()],str(g.denominator())],separators=(',', ':'))
FUNCTIONS['flint_poly_derivative']=flint_poly_derivative
FUNCTIONS['flint_fmpq_poly_derivative']=flint_fmpq_poly_derivative

_native_product_module = None

def native_product_module():
    global _native_product_module
    if _native_product_module is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _native_product_module=cython_import(str(Path(__file__).parent.parent/'polynomial_mul_flint.pyx'))
        finally:
            sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return _native_product_module
FUNCTIONS['flint_nmod_poly_product']=lambda p,a,b,method,same:json.dumps(native_product_module().native_product(p,a,b,method,bool(same)),separators=(',',':'))
FUNCTIONS['flint_nmod_poly_power']=lambda p,a,e:json.dumps(native_product_module().native_power(p,a,e),separators=(',',':'))

_native_power_module=None
def native_power_module():
    global _native_power_module
    if _native_power_module is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _native_power_module=cython_import(str(Path(__file__).parent.parent/'polynomial_power_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return _native_power_module
def native_power_frame(coeffs):
    # Compare every coefficient without storing multi-megabyte decimal vectors.
    from hashlib import sha256
    return [len(coeffs),sha256(','.join(coeffs).encode()).hexdigest()]
FUNCTIONS['flint_fmpz_poly_power']=lambda a,e,method:json.dumps(native_power_frame(native_power_module().native_integer_power(a,e,method)),separators=(',',':'))
def flint_fmpq_poly_power(a,den,e):
    coeffs,d=native_power_module().native_rational_power(a,den,e)
    return json.dumps([native_power_frame(coeffs),d],separators=(',',':'))
FUNCTIONS['flint_fmpq_poly_power']=flint_fmpq_poly_power

def ntl_gf2x_power(a,e):
    h=ntl.GF2X(list(ZZ(a).digits(2)))**e
    return format(int(ZZ([int(str(c)) for c in h.list()],2)),'x')
def ntl_gf2x_product(a,b,same):
    x=ntl.GF2X(list(ZZ(a).digits(2)));y=x if same else ntl.GF2X(list(ZZ(b).digits(2)))
    return format(int(ZZ([int(str(c)) for c in (x*y).list()],2)),'x')
def ntl_zzpx_power(p,a,e):
    h=ntl.ZZ_pX(a,p)**e
    # ntl_ZZ_p.__int__ returns a C int; its decimal representation is exact.
    return json.dumps(native_power_frame([format(int(str(c)),'x') for c in h.list()]),separators=(',',':'))
def ntl_zzpex_power(p,f,a,e):
    context=ntl.ZZ_pEContext(ntl.ZZ_pX(f,p));h=ntl.ZZ_pEX([ntl.ZZ_pE(list(ZZ(c).digits(p)),context) for c in a],context)**e
    coefficients=[[str(c) for c in str(v).strip('[]').split()] for v in h.list()]
    return json.dumps(coefficients,separators=(',',':'))
FUNCTIONS.update(ntl_gf2x_power=ntl_gf2x_power,ntl_gf2x_product=ntl_gf2x_product,ntl_zzpx_power=ntl_zzpx_power,ntl_zzpex_power=ntl_zzpex_power)

def flint_polynomial_series(kind,p,a,da,b,db,e,n):
    coefficients,d=native_power_module().native_polynomial_series(kind,p,a,da,b,db,e,n)
    return json.dumps([native_power_frame(coefficients),d],separators=(',',':'))
FUNCTIONS['flint_polynomial_series']=flint_polynomial_series

def ntl_zzpx_dense_xgcd(p,a,b):
    result=ntl.ZZ_pX(a,p).xgcd(ntl.ZZ_pX(b,p),plain=False)
    return json.dumps([[str(c) for c in h.list()] for h in result],separators=(',',':'))
def ntl_zzpex_inverse_series(p,f,a,n):
    context=ntl.ZZ_pEContext(ntl.ZZ_pX(f,p))
    h=ntl.ZZ_pEX([ntl.ZZ_pE(list(ZZ(c).digits(p)),context) for c in a],context).invert_and_truncate(n)
    return json.dumps([[str(c) for c in str(v).strip('[]').split()] for v in h.list()],separators=(',',':'))
FUNCTIONS.update(ntl_zzpx_dense_xgcd=ntl_zzpx_dense_xgcd,ntl_zzpex_inverse_series=ntl_zzpex_inverse_series)

def flint_full_product(a,da,b,db,rational,same):
    coeffs,d=native_power_module().native_full_product(a,da,b,db,int(rational),int(same))
    return json.dumps([native_power_frame(coeffs),d],separators=(',',':'))
def ntl_zzpx_full_product(p,a,b,same):
    x=ntl.ZZ_pX(a,p);y=x if same else ntl.ZZ_pX(b,p)
    return json.dumps(native_power_frame([format(int(str(c)),'x') for c in (x*y).list()]),separators=(',',':'))
def ntl_zzpex_full_product(p,f,a,b,same):
    context=ntl.ZZ_pEContext(ntl.ZZ_pX(f,p))
    def poly(v):return ntl.ZZ_pEX([ntl.ZZ_pE(list(ZZ(c).digits(p)),context) for c in v],context)
    x=poly(a);y=x if same else poly(b)
    return json.dumps([[str(c) for c in str(v).strip('[]').split()] for v in (x*y).list()],separators=(',',':'))
FUNCTIONS.update(flint_full_product=flint_full_product,ntl_zzpx_full_product=ntl_zzpx_full_product,ntl_zzpex_full_product=ntl_zzpex_full_product)

_native_modular_module=None
def native_modular_module():
    global _native_modular_module
    if _native_modular_module is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _native_modular_module=cython_import(str(Path(__file__).parent.parent/'polynomial_modular_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return _native_modular_module

def flint_modular_power(p,a,m,e,method,precision):
    import os,sys
    native=native_power_module();sys.stdout.flush();saved=os.dup(1)
    try:
        os.dup2(2,1)
        try:return json.dumps(native_power_frame(native.native_modular_power(p,a,m,e,int(method),int(precision))),separators=(',',':'))
        except RuntimeError as err:
            # Native FLINT aborts for zero precision; the dense-array boundary
            # rejects it as a RangeError before allocation. Compare rejection.
            if method==3 and precision==0 and str(err)=='Aborted':return json.dumps({'invalid_precision':True},separators=(',',':'))
            raise
    finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
def ntl_extension_modular(p,f,a,b,e,method):
    context=ntl.ZZ_pEContext(ntl.ZZ_pX(f,p))
    def poly(v):return ntl.ZZ_pEX([ntl.ZZ_pE(list(ZZ(c).digits(p)),context) for c in v],context)
    A=poly(a);B=poly(b)
    def frame(h):return [[str(c) for c in str(v).strip('[]').split()] for v in h.list()]
    try:
        result=[frame(h) for h in A.xgcd(B)] if method==2 else frame(native_modular_module().native_extension_modular_power(A,e,B,int(method)))
        return json.dumps({'value':result},separators=(',',':'))
    except Exception as err:return json.dumps({'error':type(err).__name__,'message':str(err)},separators=(',',':'))
FUNCTIONS.update(flint_modular_power=flint_modular_power,ntl_extension_modular=ntl_extension_modular)

# Exercise the product shared by modular Newton division and half-GCD directly.
FUNCTIONS['flint_nmod_gcd_product']=lambda p,a,b,same:json.dumps(native_product_module().native_product(p,a,b,0,bool(same)),separators=(',',':'))

_native_evaluation_module=None
_native_evaluation_module_error=None
def native_evaluation_module():
    global _native_evaluation_module,_native_evaluation_module_error
    if _native_evaluation_module_error is not None:raise _native_evaluation_module_error
    if _native_evaluation_module is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            # Compile the bundled internal routines against native FLINT primitives.
            # Installed FLINT canonicalizes these results; bundled FLINT retains
            # raw denominators. The tested C bodies are read verbatim, not reimplemented.
            import tempfile
            root=Path(__file__).resolve().parents[4]/'reference'/'flint'/'src'
            routines=[
                ('fmpz_poly/evaluate_horner_fmpq.c','_fmpz_poly_evaluate_horner_fmpq'),
                ('fmpz_poly/evaluate_divconquer_fmpq.c','_fmpz_poly_evaluate_divconquer_fmpq'),
                ('fmpz_poly/evaluate_fmpq.c','_fmpz_poly_evaluate_fmpq'),
                ('fmpq_poly/evaluate.c','_fmpq_poly_evaluate_fmpq'),
                ('fmpq_poly/evaluate.c','_fmpq_poly_evaluate_fmpz'),
                ('fmpz/get.c','fmpz_get_d'),
            ]
            bodies=[]
            for filename,name in routines:
                source=(root/filename).read_text();start=source.index(('double' if name=='fmpz_get_d' else 'void')+'\n'+name+'(')
                opening=source.index('{',start);depth=1;end=opening+1
                while depth:
                    depth+=(source[end]=='{')-(source[end]=='}');end+=1
                bodies.append(source[start:end])
            original='\n'.join(bodies)
            for _,name in routines:original=original.replace(name,'audit'+name)
            headers='\n'.join('#include <flint/'+name+'.h>' for name in ['fmpz','fmpz_poly','fmpz_vec','fmpq','longlong','mpn_extras'])
            block='cdef extern from *:\n    """\n'+headers+'\ntypedef long slong;\ntypedef unsigned long ulong;\n'+'#define DOUBLE_MAX WORD(9007199254740992)\n#define DOUBLE_MIN WORD(-9007199254740992)\n'+original+'\n    """\n'
            block+='    void audit_fmpz_poly_evaluate_fmpq(fmpz_t,fmpz_t,const fmpz*,long,const fmpz_t,const fmpz_t)\n'
            block+='    void audit_fmpq_poly_evaluate_fmpq(fmpz_t,fmpz_t,const fmpz*,const fmpz_t,long,const fmpz_t,const fmpz_t)\n'
            block+='    void audit_fmpq_poly_evaluate_fmpz(fmpz_t,fmpz_t,const fmpz*,const fmpz_t,long,const fmpz_t)\n'
            block+='    double auditfmpz_get_d(const fmpz_t)\n'
            template=(Path(__file__).parent.parent/'polynomial_evaluation_native.pyx').read_text()
            with tempfile.TemporaryDirectory(prefix='sage-bundled-evaluation-') as directory:
                generated=Path(directory)/'bundled_evaluation.pyx'
                generated.write_text(template.replace('# BUNDLED_EVALUATION_ROUTINES',block))
                _native_evaluation_module=cython_import(str(generated))
        except Exception as err:
            _native_evaluation_module_error=err
            raise
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return _native_evaluation_module
FUNCTIONS['flint_polynomial_evaluation']=lambda kind,a,den,num,xden:json.dumps(native_evaluation_module().native_evaluation(int(kind),a,den,num,xden),separators=(',',':'))
def ntl_extension_evaluation(p,f,a,x):
    context=ntl.ZZ_pEContext(ntl.ZZ_pX(f,p))
    def coefficient(c):return ntl.ZZ_pE(list(ZZ(c).digits(p)),context)
    h=native_evaluation_module().native_extension_evaluation(ntl.ZZ_pEX(list(map(coefficient,a)),context),coefficient(x))
    return json.dumps(str(h).strip('[]').split(),separators=(',',':'))
FUNCTIONS['ntl_extension_evaluation']=ntl_extension_evaluation
FUNCTIONS['flint_nmod_evaluation']=lambda p,a,x:str(PolynomialRing(Zmod(p),'x')(a)(Zmod(p)(x)))
FUNCTIONS['flint_polynomial_composition']=lambda kind,a,da,b,db,p:json.dumps(native_evaluation_module().native_composition(int(kind),a,da,b,db,int(p)),separators=(',',':'))

FUNCTIONS['ntl_prime_evaluation']=lambda p,a,x:str(PolynomialRing(Zmod(p),'x',implementation='NTL')(a)(Zmod(p)(x)))

def flint_integer_double(value):
    import struct
    return str(int.from_bytes(struct.pack('>d',native_evaluation_module().native_integer_double(value)),'big'))
FUNCTIONS['flint_integer_double']=flint_integer_double
