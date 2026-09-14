"""Execute the bundled FLINT factor kernels through a minimal C ABI adapter."""
import ctypes
import hashlib
from pathlib import Path
import subprocess
import tempfile
from bundled_flint import bundled_flint_build

_native = None


def native_factor(op, p, d, a, b, c, seed1=-1, seed2=0):
    global _native
    op, p, d = int(op), int(p), int(d)
    if not 0 <= op <= 20: raise ValueError('unknown native factor operation')
    if not 2 <= p < 2**64: raise ValueError('modulus outside native word contract')
    class RangeError(Exception): pass
    def normalized(v):
        v=[int(x)%p for x in v]
        while v and not v[-1]:v.pop()
        return v
    A,B,C=map(normalized,(a,b,c))
    # Only undefined/aborting native calls are guarded. Valid calls below run C.
    if op == 1 and not 1<=d<=2**53-1:raise RangeError('deflation must be positive')
    if op == 2 and not 0<=d<=2**53-1:raise RangeError('inflation must be nonnegative')
    if op == 3 and len(B)<=1:raise RangeError('factor must have positive degree')
    if op in (4,5,6) and not 0<=d<=2**53-1:raise RangeError('length must be nonnegative')
    if op == 7 and (len(A)<3 or not 1<=d<=2**53-1):raise RangeError('equal-degree splitting requires degree at least two and positive d')
    if op == 8 and (not 1<=d<=2**53-1 or len(A)<=1 or (len(A)-1)%d):raise RangeError('invalid equal-degree factorization input')
    if op == 9 and len(A)<=1:raise RangeError('polynomial division by zero' if not A else 'polynomial must have positive degree')
    if op in (10,11) and not A:raise RangeError('polynomial division by zero')
    if (op == 15 and not B) or (op == 16 and not C):raise RangeError('polynomial division by zero')
    if op == 16 and len(A)>=len(C):raise RangeError('outer degree must be smaller than modulus degree')
    if op == 17 and (not 0<=d<=2**53-1):raise RangeError('invalid polynomial count')
    if op == 17 and len(a) != d * (len(C)-1):raise ValueError('invalid native composition vector')
    if op == 18 and (len(a) != seed1*d or len(b) != d*seed2):raise ValueError('invalid native matrix shape')
    if op == 19 and d<1:raise RangeError('matrix must have at least one row')
    if op == 19 and (len(C)<2 or len(a)!=d*(len(C)-1) or len(B)>=len(C)):raise ValueError('invalid native row evaluation')
    if op == 20 and len(C)<=1:raise RangeError('modulus must have positive degree')
    if op == 20 and not 2<=len(C)<len(B):raise ValueError('invalid native matrix reduction')
    if _native is None:
        build = bundled_flint_build()
        source = Path(__file__).resolve().parents[1]/'native/flint_factor.c'
        key = hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16]
        folder = Path(tempfile.gettempdir())/('sage-bundled-flint-factor-'+key)
        folder.mkdir(exist_ok=True)
        library = folder/'oracle.so'
        if not library.exists():
            subprocess.run(['cc','-shared','-fPIC','-O2','-I'+str(build/'src'),str(source),
                            '-L'+str(build),'-lflint','-Wl,-rpath,'+str(build),'-o',str(library)],
                           check=True, capture_output=True)
        _native = ctypes.CDLL(str(library)).audit_flint_factor
        word, length = ctypes.c_ulong, ctypes.c_long
        ptr = ctypes.POINTER(word)
        _native.argtypes = [ctypes.c_int,word,length,ptr,length,ptr,length,ptr,length,
                           ctypes.c_int,word,word,ptr,length]
        _native.restype = length
    arrays = [(ctypes.c_ulong*len(v))(*(int(x)%p for x in v)) for v in (a,b,c)]
    cap = max(64, 8*(len(a)+len(b)+len(c)+abs(d)+1))
    while True:
        out = (ctypes.c_ulong*cap)()
        n = _native(op,p,d,arrays[0],len(a),arrays[1],len(b),arrays[2],len(c),
                    int(seed1>=0),int(seed1),int(seed2),out,cap)
        if n >= 0: return ','.join(map(str,out[:n]))
        cap *= 2
