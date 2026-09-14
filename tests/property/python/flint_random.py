"""Compile the bundled FLINT random kernels, independent of the installed FLINT version."""
import ctypes
import hashlib
from pathlib import Path
import subprocess
import tempfile

_native = None

def native_rng(op, seed1, seed2, limit, count):
    global _native
    # Out-of-range native shifts have no portable result. Match the documented guard.
    if op in (3,4) and count>0 and not 0<=limit<=64:
        class RangeError(Exception):pass
        raise RangeError('bits must be between 0 and 64')
    if _native is None:
        root = Path(__file__).resolve().parents[3]
        header = (root/'reference/flint/src/flint.h.in').read_text()
        extras = (root/'reference/flint/src/ulong_extras.h').read_text()
        source = (root/'reference/flint/src/ulong_extras/randomisation.c').read_text()
        init = header[header.index('FLINT_INLINE\nvoid flint_rand_init('):header.index('FLINT_DEPRECATED void _flint_rand_init_gmp_state')]
        limb = extras[extras.index('ULONG_EXTRAS_INLINE ulong _n_randlimb('):extras.index('#else', extras.index('ULONG_EXTRAS_INLINE ulong _n_randlimb('))]
        randint = extras[extras.index('ULONG_EXTRAS_INLINE ulong _n_randint('):extras.index('/* Basic arithmetic')]
        kernels = source[source.index('ulong n_randlimb('):source.index('ulong n_randprime(')]
        # Only platform definitions are supplied here; algorithm bodies above are
        # extracted verbatim from the original source, not translated into Python.
        shim = r'''
#include <stdint.h>
#include <stddef.h>
#include <limits.h>
typedef unsigned long ulong;
typedef struct { void *__gmp_state; ulong __randval; ulong __randval2; } flint_rand_struct;
typedef flint_rand_struct flint_rand_t[1];
#define FLINT_UNUSED(x) x
#define FLINT_INLINE static inline
#define ULONG_EXTRAS_INLINE static inline
#define FLINT64 1
#define FLINT_BITS 64
#define UWORD(x) ((ulong)(x))
#define WORD(x) ((long)(x))
#define UWORD_MAX ULONG_MAX
#define WORD_MAX LONG_MAX
#define COEFF_MAX (LONG_MAX >> 1)
#define umul_ppmm(hi,lo,a,b) do { __uint128_t t=(__uint128_t)(a)*(b); (hi)=(ulong)(t>>64); (lo)=(ulong)t; } while(0)
ulong n_randint(flint_rand_t, ulong);
'''
        wrapper = r'''
void native_sequence(ulong *out, int op, int seeded, ulong s1, ulong s2, ulong limit, long count) {
    flint_rand_t state;
    flint_rand_init(state);
    if (seeded) flint_rand_set_seed(state,s1,s2);
    for(long i=0;i<count;i++) {
        switch(op) {
          case 0: out[i]=n_randlimb(state);break;
          case 1: out[i]=n_randint(state,limit);break;
          case 2: out[i]=n_urandint(state,limit);break;
          case 3: out[i]=n_randbits(state,limit);break;
          case 4: out[i]=n_randtest_bits(state,limit);break;
          case 5: out[i]=n_randtest(state);break;
          default:out[i]=n_randtest_not_zero(state);
        }
    }
    flint_rand_clear(state);
    flint_rand_get_seed(out+count,out+count+1,state);
}
'''
        code = shim+init+limb+randint+kernels+wrapper
        folder = Path(tempfile.gettempdir())/('sage-bundled-flint-rng-'+hashlib.sha256(code.encode()).hexdigest()[:16])
        folder.mkdir(exist_ok=True)
        cfile=folder/'oracle.c'; library=folder/'oracle.so'
        if not library.exists():
            cfile.write_text(code)
            subprocess.run(['cc','-shared','-fPIC','-O2',str(cfile),'-o',str(library)],check=True,capture_output=True)
        _native=ctypes.CDLL(str(library)).native_sequence
        _native.argtypes=[ctypes.POINTER(ctypes.c_ulong),ctypes.c_int,ctypes.c_int,ctypes.c_ulong,ctypes.c_ulong,ctypes.c_ulong,ctypes.c_long]
        _native.restype=None
    output=(ctypes.c_ulong*(int(count)+2))()
    _native(output,int(op),int(seed1>=0),int(seed1),int(seed2),int(limit),int(count))
    return ','.join(map(str,output))
