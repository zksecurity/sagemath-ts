"""Original MPFR decimal conversion and Sage real-literal oracles."""
from sage.all import *
import json
_native=None

def native_module():
    global _native
    if _native is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _native=cython_import(str(Path(__file__).parent.parent/'mpfr_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return _native

TEXTS=['','x',' ',' +.5','1.2junk','1e','1e+','-nan','@NaN@','nan(foo)','@Inf@','Inf','-Infinity','+0','-0','1.',' .1 ','1_0','1.2e3 ','0e999999999999999999999','1e4096','1e-4096','@nan','nan@','infinityx','-@NaN@']
def conversion_text(kind,a,b):
    if kind==0:return str(a)
    if kind==1:return str(a)+'e'+str(b)
    if kind==2:return repr(float(a)/float(b))
    if kind==3:return TEXTS[int(a)]
    if kind==4:return ('-' if a<0 else '')+'0.'+'0'*int(b)+str(abs(a))
    return str(int(a)*2**int(b)) if b>=0 else str(int(a)*5**int(-b))+'e'+str(b)

def mpfr_conversion(kind,a,b,precision,digits):
    from sage.rings.real_mpfr import create_RealNumber
    text=conversion_text(kind,a,b)
    p=int(precision) if precision else int(create_RealNumber(text).precision())
    return json.dumps(native_module().native_conversion(text,p,int(digits)),separators=(',',':'))
FUNCTIONS={'mpfr_conversion':mpfr_conversion}

def real_literal_frame(x):
    import math,struct
    from sage.rings.real_mpfr import RealLiteral
    try:
        q=x.exact_rational();rational=[str(q.numerator()),str(q.denominator())]
    except Exception as e:rational=[type(e).__name__,str(e)]
    # Bundled RealLiteral.__float__ reparses its original text at 53 bits.
    # Installed Sage 10.3 lacks this override; execute the bundled expression.
    v=float(x.numerical_approx(53)) if isinstance(x,RealLiteral) else float(x)
    s,m,e=x.sign_mantissa_exponent()
    return [str(x),str(x.parent()),int(x.precision()),rational,[int(s),str(m),str(e)],'NaN' if math.isnan(v) else struct.pack('>d',v).hex(),type(x).__name__,x.literal if isinstance(x,RealLiteral) else None,int(x.base) if isinstance(x,RealLiteral) else None]

OBJECTS=[None,True,False,float('nan'),float('inf'),-float('inf'),-0.0,1.0,1e-5,1e-4,1e16,1e20,1.133759543500045e153]
def real_literal_conversion(kind,a,b,precision,pad,min_prec,action,arg):
    from sage.rings.real_mpfr import RealLiteral,create_RealNumber
    value=OBJECTS[int(a)] if kind==6 else float(a)/float(b) if kind==2 else conversion_text(kind,a,b)
    x=RealLiteral(RealField(precision),str(value),10) if precision else create_RealNumber(value,pad=int(pad),min_prec=int(min_prec))
    if action==1:x=-x
    elif action==2:x=x.numerical_approx(int(arg))
    elif action==3:x=x.numerical_approx(digits=int(arg))
    elif action==4:x=abs(x)
    elif action==5:x=-x.numerical_approx(int(arg))
    return json.dumps(real_literal_frame(x),separators=(',',':'))
FUNCTIONS['real_literal_conversion']=real_literal_conversion

SEQUENCES=[['-1','','nan','-nan','@nan','-inf','nan(foo)','1e+',''],['1.234567890123456789','1.2junk','-0','nan','@NaN@','-0.25',''],['-1','nan','-nan','+nan','-inf','nan(payload)','0','-NaN']]
def mpfr_sequence(kind,precision,digits):
    return json.dumps(native_module().native_sequence(SEQUENCES[int(kind)],int(precision),int(digits)),separators=(',',':'))
FUNCTIONS['mpfr_sequence']=mpfr_sequence

def real_literal_format(a,exponent,precision,digits,flags):
    from sage.rings.real_mpfr import RealLiteral
    x=RealLiteral(RealField(precision),str(a)+'e'+str(exponent))
    return x.str(digits=int(digits),truncate=bool(flags&1),skip_zeroes=bool(flags&2),no_sci=[None,False,True,2][int(flags>>2)],e='E')
FUNCTIONS['real_literal_format']=real_literal_format


def mpfr_double_conversion(raw,precision,digits):
    return json.dumps(native_module().native_double(int(raw),int(precision),int(digits)),separators=(',',':'))
FUNCTIONS['mpfr_double_conversion']=mpfr_double_conversion

MPFR_DOUBLE_SEQUENCES = [[13830554455654793216, 9221120237041090560], [9223372036854775808, 18444492273895866368, 9221120237041090560], [18442240474082181120, 9221120237041090560, 0], [4608308318706860032, 13832806255468478464, 1, 9223372036854775809], [9218868437227405311, 18442240474082181119, 4503599627370496], [9221120237041090560, 18442240474082181120, 9218868437227405312, 9223372036854775808], [4607182418800017409, 4607182418800017407, 0], [9218868437227405313, 18442240474082181121, 13830554455654793216]]

def mpfr_double_sequence(kind,precision,digits):
    return json.dumps(native_module().native_double_sequence(MPFR_DOUBLE_SEQUENCES[int(kind)],int(precision),int(digits)),separators=(',',':'))
FUNCTIONS['mpfr_double_sequence']=mpfr_double_sequence


def real_format53(raw):
    import struct
    value=RealField(53)(struct.unpack('>d',int(raw).to_bytes(8,'big'))[0])
    return json.dumps([repr(value),format(value,'.6f')],separators=(',',':'))
FUNCTIONS['real_format53']=real_format53


def real_numeric_constructor(kind,value,precision,action):
    import struct
    if kind==0:source=struct.unpack('>d',int(value).to_bytes(8,'big'))[0]
    elif kind==1:source=Integer(value)
    else:source=conversion_text(3,value,0)
    x=RealField(int(precision))(source)
    if action==1:x=-x
    elif action==2:x=abs(x)
    return json.dumps([real_literal_frame(x),x.str(),x.str(digits=6)],separators=(',',':'))
FUNCTIONS['real_numeric_constructor']=real_numeric_constructor


def mpfr_integer_conversion(value,precision,digits):
    return json.dumps(native_module().native_integer(value,int(precision),int(digits)),separators=(',',':'))
FUNCTIONS['mpfr_integer_conversion']=mpfr_integer_conversion


REAL_OBSERVER_TEXTS=['0','-0','NaN','+infinity','-infinity','1','-1','1.00000000000000000001','-1.00000000000000000001','9007199254740993','1e400','-1e400','1e-400','-1e-400','0.5','-0.5', '-NaN', '-@NaN@']
REAL_OBSERVER_OPS=['sign','is_NaN','is_positive_infinity','is_negative_infinity','is_infinity','is_integer','is_square','multiplicative_order','floor','ceil','round','trunc']
def real_observer(kind,a,b,precision,op):
    text=REAL_OBSERVER_TEXTS[int(a)] if kind==0 else str(a)+'e'+str(b)
    x=RealField(int(precision))(text)
    return str(getattr(x,REAL_OBSERVER_OPS[int(op)])())
FUNCTIONS['real_observer']=real_observer

def mpfr_observer(kind,a,b,precision,target,operation,mode,integer):
    text=REAL_OBSERVER_TEXTS[int(a)] if kind==0 else str(a)+'e'+str(b)
    return json.dumps(native_module().native_observer(text,int(precision),int(target),int(operation),int(mode),int(integer)),separators=(',',':'))
FUNCTIONS['mpfr_observer']=mpfr_observer


def mpfr_rint_extreme(sign,exponent,precision,target,mode):
    return json.dumps(native_module().native_rint_extreme(int(sign),int(exponent),int(precision),int(target),int(mode)),separators=(',',':'))
FUNCTIONS['mpfr_rint_extreme']=mpfr_rint_extreme


def real_fraction(kind,a,b,precision,literal):
    from sage.rings.real_mpfr import RealLiteral
    text=REAL_OBSERVER_TEXTS[int(a)] if kind==0 else str(a)+'e'+str(b)
    R=RealField(int(precision))
    x=RealLiteral(R,text) if literal else R(text)
    return json.dumps(real_literal_frame(x.frac()),separators=(',',':'))
FUNCTIONS['real_fraction']=real_fraction

def real_compare(kind,a,b,precision,other_kind,other_a,other_b,other_precision,flags):
    from sage.rings.real_mpfr import RealLiteral
    left=REAL_OBSERVER_TEXTS[int(a)] if kind==0 else str(a)+'e'+str(b)
    right=REAL_OBSERVER_TEXTS[int(other_a)] if other_kind==0 else str(other_a)+'e'+str(other_b)
    R=RealField(int(precision));S=RealField(int(other_precision))
    x=RealLiteral(R,left) if flags&1 else R(left)
    y=float(right) if flags&4 else RealLiteral(S,right) if flags&2 else S(right)
    if flags&4 and flags&1 and precision>53:
        # Sage 10.3 lacks bundled RealLiteral.__float__. RDF's ToRDF calls
        # that override, which reparses the literal at 53 bits before float().
        x=RDF(float(x.numerical_approx(53)))
    return json.dumps([int(x>y)-int(x<y),bool(x==y)],separators=(',',':'))
FUNCTIONS['real_compare']=real_compare

def mpfr_copy_fraction(kind,a,b,precision,target,alias):
    text=REAL_OBSERVER_TEXTS[int(a)] if kind==0 else str(a)+'e'+str(b)
    return json.dumps(native_module().native_copy_fraction(text,int(precision),int(target),bool(alias)),separators=(',',':'))
FUNCTIONS['mpfr_copy_fraction']=mpfr_copy_fraction


BINARY_TEXTS = ['0', '-0', 'nan', '-nan', 'inf', '-inf', '1', '-1', '1.5', '-1.5', '1.25', '-1.25', '1.75', '-1.75', '0.1', '-0.1', '1e1000', '-1e1000', '1e-1000', '-1e-1000', '1.00000000000000000000000000000000000001', '-1.00000000000000000000000000000000000001']
def mpfr_binary(left,right,lp,rp,target,operation,alias,le,re):
    return json.dumps(native_module().native_binary(BINARY_TEXTS[left],BINARY_TEXTS[right],int(lp),int(rp),int(target),int(operation),int(alias),int(le),int(re)),separators=(',',':'))
FUNCTIONS['mpfr_binary'] = mpfr_binary

def mpfr_binary_halfway(precision,extra,parity,sign,other_sign,gap,swapped,alias):
    # Exact integer mantissa with one halfway bit; no decimal oracle approximation.
    mantissa=(2**(int(precision)-1)+int(parity))*2**int(extra)+2**(int(extra)-1)
    left=str(int(sign)*mantissa);right=str(int(other_sign))
    lp=int(precision+extra);rp=int(precision+3);le=int(precision+extra);re=le-int(gap)
    if swapped:left,right,lp,rp,le,re=right,left,rp,lp,re,le
    return json.dumps(native_module().native_binary(left,right,lp,rp,int(precision),0,int(alias),le,re),separators=(',',':'))
FUNCTIONS['mpfr_binary_halfway']=mpfr_binary_halfway
