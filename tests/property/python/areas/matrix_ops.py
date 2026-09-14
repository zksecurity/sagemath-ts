"""SageMath side of the ``matrix_ops`` property-test area.

Cases: tests/property/cases/matrix_ops.cases.json
"""

from sage.all import *


def matrix_determinant_2x2(a, b, c, d):
    """Compute determinant of 2x2 matrix [[a,b],[c,d]]."""
    M = matrix(ZZ, [[a, b], [c, d]])
    return M.determinant()


def matrix_determinant_3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33):
    """Compute determinant of 3x3 matrix."""
    M = matrix(ZZ, [[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])
    return M.determinant()


def matrix_determinant_4x4(*args):
    """Compute determinant of 4x4 matrix."""
    M = matrix(ZZ, 4, 4, list(args))
    return M.determinant()


def matrix_rank_2x3(a11, a12, a13, a21, a22, a23):
    """Compute rank of 2x3 matrix."""
    M = matrix(ZZ, [[a11, a12, a13], [a21, a22, a23]])
    return M.rank()


def matrix_rank_3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33):
    """Compute rank of 3x3 matrix."""
    M = matrix(ZZ, [[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])
    return M.rank()


def matrix_hnf_2x2(a, b, c, d):
    """Compute HNF of 2x2 matrix, return as list of lists."""
    M = matrix(ZZ, [[a, b], [c, d]])
    H = M.hermite_form()
    return [[int(H[i,j]) for j in range(2)] for i in range(2)]


def matrix_hnf_3x3(*args):
    """Compute HNF of 3x3 matrix, return as list of lists."""
    M = matrix(ZZ, 3, 3, list(args))
    H = M.hermite_form()
    return [[int(H[i,j]) for j in range(3)] for i in range(3)]


def matrix_snf_2x2(a, b, c, d):
    """Compute SNF of 2x2 matrix, return diagonal elements."""
    M = matrix(ZZ, [[a, b], [c, d]])
    D, _, _ = M.smith_form()
    return [int(D[i,i]) for i in range(2)]


def matrix_snf_3x3(*args):
    """Compute SNF of 3x3 matrix, return diagonal elements."""
    M = matrix(ZZ, 3, 3, list(args))
    D, _, _ = M.smith_form()
    return [int(D[i,i]) for i in range(3)]


def matrix_lll_2x2(a, b, c, d):
    """Compute LLL-reduced basis of 2x2 matrix, return as list of lists."""
    M = matrix(ZZ, [[a, b], [c, d]])
    L = M.LLL()
    return [[int(L[i,j]) for j in range(2)] for i in range(2)]


def matrix_lll_3x3(*args):
    """Compute LLL-reduced basis of 3x3 matrix, return as list of lists."""
    M = matrix(ZZ, 3, 3, list(args))
    L = M.LLL()
    return [[int(L[i,j]) for j in range(3)] for i in range(3)]


def matrix_elementary_divisors_2x2(a, b, c, d):
    """Compute elementary divisors of 2x2 matrix."""
    M = matrix(ZZ, [[a, b], [c, d]])
    divs = M.elementary_divisors()
    return [int(d) for d in divs]


def matrix_elementary_divisors_3x3(*args):
    """Compute elementary divisors of 3x3 matrix."""
    M = matrix(ZZ, 3, 3, list(args))
    divs = M.elementary_divisors()
    return [int(d) for d in divs]


FUNCTIONS = {
    # Determinant
    'determinant_2x2': matrix_determinant_2x2,
    'determinant_3x3': matrix_determinant_3x3,
    'determinant_4x4': matrix_determinant_4x4,
    # Rank
    'rank_2x3': matrix_rank_2x3,
    'rank_3x3': matrix_rank_3x3,
    # HNF
    'hnf_2x2': matrix_hnf_2x2,
    'hnf_3x3': matrix_hnf_3x3,
    # SNF
    'snf_2x2': matrix_snf_2x2,
    'snf_3x3': matrix_snf_3x3,
    # LLL
    'lll_2x2': matrix_lll_2x2,
    'lll_3x3': matrix_lll_3x3,
    # Elementary divisors
    'elementary_divisors_2x2': matrix_elementary_divisors_2x2,
    'elementary_divisors_3x3': matrix_elementary_divisors_3x3,
}

_m4ri_native=None
_m4ri_native_error=None
def m4ri_native_module():
    global _m4ri_native,_m4ri_native_error
    if _m4ri_native_error is not None:raise _m4ri_native_error
    if _m4ri_native is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _m4ri_native=cython_import(str(Path(__file__).parent.parent/'m4ri_native.pyx'))
        except Exception as e:_m4ri_native_error=e;raise
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return _m4ri_native

def m4ri_packed_fixture(m,n,seed):
    state=int(seed);rows=[]
    for i in range(m):
        row=0
        for j in range(0,n,64):
            state=(state*6364136223846793005+1442695040888963407)&((1<<64)-1)
            row|=state<<j
        rows.append(row&((1<<n)-1))
    return rows

def m4ri_product(method,m,k,n,seed,parameter,square):
    import json
    a=m4ri_packed_fixture(m,k,seed)
    b=a if square else m4ri_packed_fixture(m if method==3 else k,k if method==3 else n,seed+1)
    return json.dumps(m4ri_native_module().native_m4ri(int(method),int(m),int(k),int(n),a,b,int(parameter),bool(square)),separators=(',',':'))
FUNCTIONS['m4ri_product']=m4ri_product

def m4ri_strassen_boundary(m,k,n,seed,parameter,square):
    import json
    from m4ri_pinned import strassen
    return json.dumps(strassen(m,k,n,seed,parameter,square),separators=(',',':'))
FUNCTIONS['m4ri_strassen_boundary']=m4ri_strassen_boundary

def binary_matrix_product(method,m,k,right_rows,n,seed,parameter,square):
    import json
    aa=m4ri_packed_fixture(m,k,seed);bb=aa if square else m4ri_packed_fixture(right_rows,n,seed+1)
    A=matrix(GF(2),m,k,[(row>>j)&1 for row in aa for j in range(k)])
    B=A if square else matrix(GF(2),right_rows,n,[(row>>j)&1 for row in bb for j in range(n)])
    if method==0:C=A*B
    elif method==1:
        if m and k and n and k!=right_rows:return json.dumps('invalid-native-classical-dimensions')
        C=A._multiply_classical(B)
    elif method==2:C=A._multiply_m4rm(B,parameter)
    elif k==right_rows and m and k and n and 0<parameter<=128:
        # Isolate native empty-window faults at small cutoffs. The pinned
        # library computes the result; reconstruct the original return parent.
        from m4ri_pinned import strassen
        result=strassen(m,k,n,seed,parameter,int(square))
        C=matrix(GF(2),m,n,[(int(row)>>j)&1 for row in result[2] for j in range(n)])
    else:C=A._multiply_strassen(B,parameter)
    values=[str(sum(ZZ(C[i,j])<<j for j in range(C.ncols()))) for i in range(C.nrows())]
    return json.dumps([str(C.nrows()),str(C.ncols()),values,C is A,C is B],separators=(',',':'))
FUNCTIONS['binary_matrix_product']=binary_matrix_product


def binary_matrix_access(action,m,n,seed,args):
    import json
    A=matrix(GF(2),m,n,[(i*n+j+seed)%2 for i in range(m) for j in range(n)])
    if action==0:result=[int(x) for x in A.row(*args)]
    elif action==1:
        B=A.submatrix(*args)
        result=[str(B.nrows()),str(B.ncols()),[str(x) for x in B.list()]]
    elif action==2:result=str(A[args[0],args[1]])
    else:
        A[args[0],args[1]]=args[2]
        result=[str(x) for x in A.list()]
    return json.dumps(result,separators=(',',':'))
FUNCTIONS['binary_matrix_access']=binary_matrix_access


def binary_matrix_index_conversion(action,m,n,kind,position,value):
    import json
    if kind in (0,1):x=int(value)
    elif kind==2:x=ZZ(value)
    elif kind==3:x=bool(value)
    elif kind==4:x=float(value)+0.5
    elif kind==5:x=float('nan')
    elif kind==6:x=float('inf')
    elif kind==7:x=str(value)
    elif kind==8:x=None
    else:x=QQ(value)/2
    A=matrix(GF(2),m,n,[(i*n+j)%2 for i in range(m) for j in range(n)])
    if action==0:result=[int(c) for c in A.row(x)]
    elif action==1:
        args=[0,0];args[position]=x;result=str(A[args[0],args[1]])
    elif action==2:
        args=[0,0];args[position]=x;A[args[0],args[1]]=1;result=[str(c) for c in A.list()]
    else:
        args=[0,0,0,0];args[position]=x;B=A.submatrix(*args)
        result=[str(B.nrows()),str(B.ncols()),[str(c) for c in B.list()]]
    return json.dumps(result,separators=(',',':'))
FUNCTIONS['binary_matrix_index_conversion']=binary_matrix_index_conversion

def binary_matrix_slice_words(m,n,seed,r0,c0,r1,c1):
    import json
    rows=m4ri_packed_fixture(m,n,seed)
    A=matrix(GF(2),m,n,[(row>>j)&1 for row in rows for j in range(n)])
    B=A.submatrix(r0,c0,r1-r0,c1-c0)
    frame=lambda X:[str(X.nrows()),str(X.ncols()),[str(sum(ZZ(X[i,j])<<j for j in range(X.ncols()))) for i in range(X.nrows())]]
    return json.dumps([frame(B),m4ri_native_module().native_submatrix(m,n,rows,r0,c0,r1,c1),B is A],separators=(',',':'))
FUNCTIONS['binary_matrix_slice_words']=binary_matrix_slice_words

def binary_matrix_density(method,m,n,pattern,seed,resolution):
    import json,struct,math
    from m4ri_pinned import density
    if method==2:
        value=density(m,n,pattern,seed,resolution)
        return json.dumps('NaN' if math.isnan(value) else struct.pack('>d',value).hex())
    rows=m4ri_packed_fixture(m,n,seed)
    entries=[((rows[i]>>j)&1) if pattern==5 else 0 if pattern==0 else 1 if pattern==1 else int(j>=n//2) if pattern==2 else int(j<n//2) if pattern==3 else (i*n+j)%2 if pattern==4 else int(i*n+j<seed) for i in range(m) for j in range(n)]
    A=matrix(GF(2),m,n,entries)
    if method==0:
        value=A.density();return json.dumps([type(value).__name__,str(value)],separators=(',',':'))
    if m and n==0:density(m,n,pattern,seed,1) # isolate the unguarded original native fault
    from areas.real_literals import real_literal_frame
    return json.dumps(real_literal_frame(A.density(approx=True)),separators=(',',':'))
FUNCTIONS['binary_matrix_density']=binary_matrix_density

BINARY_ALGORITHMS=['heuristic','m4ri','pluq','classical','ple','linbox','default','bogus']
def binary_matrix_elimination(method,m,n,seed,algorithm,reduced,sequence):
    import json
    rows=m4ri_packed_fixture(m,n,seed)
    A=matrix(GF(2),m,n,[(row>>j)&1 for row in rows for j in range(n)])
    algo=BINARY_ALGORITHMS[int(algorithm)]
    frame=lambda B:[str(x) for x in B.list()]
    def observe(call):
        try:return ['ok',call()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():
        if method==0:return int(A.rank(algo))
        if method==1:
            A.echelonize(algo,0,bool(reduced));return frame(A)
        return frame(A.echelon_form(algo,0,reduced=bool(reduced)))
    if sequence==1:A.rank()
    elif sequence==2:A.echelonize('m4ri',0,False)
    elif sequence==3:A.echelon_form('m4ri',0,reduced=False)
    def invalid_echelon():
        x=A.echelonize('bogus')
        # The empty-matrix self return is observed in the dedicated return test.
        return None
    return json.dumps([observe(call),frame(A),observe(lambda:int(A.rank('bogus'))),observe(invalid_echelon),frame(A)],separators=(',',':'))
FUNCTIONS['binary_matrix_elimination']=binary_matrix_elimination

def m4ri_echelon(m,n,seed,algorithm,full,k):
    import json
    from m4ri_pinned import echelon
    return json.dumps(echelon(m,n,seed,algorithm,full,k),separators=(',',':'))
FUNCTIONS['m4ri_echelon']=m4ri_echelon

def m4ri_ple(m,n,seed,kind,k):
    import json
    from m4ri_pinned import ple
    return json.dumps(ple(m,n,seed,kind,k),separators=(',',':'))
FUNCTIONS['m4ri_ple']=m4ri_ple

def m4ri_trsm(m,n,seed,upper,cutoff):
    import json
    from m4ri_pinned import trsm
    return json.dumps(trsm(m,n,seed,upper,cutoff),separators=(',',':'))
FUNCTIONS['m4ri_trsm']=m4ri_trsm

def binary_matrix_cache(m,n,seed,prepare,operation):
    import json
    rows=m4ri_packed_fixture(m,n,seed)
    A=matrix(GF(2),m,n,[(row>>j)&1 for row in rows for j in range(n)])
    if prepare==1:A.rank()
    elif prepare==2:A.echelonize('m4ri',0,False)
    elif prepare==3:A.echelon_form('m4ri',0,reduced=False)
    elif prepare==4:A=A.echelon_form('m4ri',0,reduced=False)
    elif prepare==5:A.set_immutable()
    elif prepare==6:A.echelonize('classical')
    before=[A.is_immutable(),A.is_mutable()]
    def observe(call):
        try:return ['ok',call()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():
        nonlocal A
        # Installed M4RI faults before the bundled release's empty-row guard.
        # Run the pinned kernel and the original cache invalidation instead.
        if A.is_mutable() and ((m==0 and n>1 and operation in (4,7)) or (n==0 and m>1 and operation in (2,6))):
            from m4ri_pinned import empty_swap
            A._clear_cache();axis=int(m==0);size=n if axis else m
            empty_swap(m,n,axis,0,size-1)
            return None
        if operation==0:A[0,0]=1-int(A[0,0])
        elif operation==1:A[m,n]=1
        elif operation==2:A.swap_rows(0,max(0,m-1))
        elif operation==3:A.swap_rows(-1,m)
        elif operation==4:A.swap_columns(0,max(0,n-1))
        elif operation==5:A.swap_columns(-1,n)
        elif operation==6:A.permute_rows(SymmetricGroup(m)([1+(i+1)%m for i in range(m)]))
        elif operation==7:A.permute_columns(SymmetricGroup(n)([1+(i+1)%n for i in range(n)]))
        elif operation==8:A.randomize(0)
        elif operation==9:A.randomize(-1)
        elif operation==10:A.echelonize('bogus')
        elif operation==11:A.echelonize('m4ri')
        elif operation==12:return [str(x) for x in A.echelon_form('bogus').list()]
        elif operation==13:A._clear_cache()
        elif operation==14:A.set_immutable()
        elif operation==15:A=A.__copy__()
        elif operation==16:binary_ordering_module().doubly_lexical_ordering(A,inplace=True)
        elif operation==17:binary_ordering_module().doubly_lexical_ordering(A,inplace=False)
        return None
    result=observe(call)
    def invalid_echelon():A.echelonize('bogus')
    return json.dumps([before,result,[str(x) for x in A.list()],[A.is_immutable(),A.is_mutable()],observe(lambda:int(A.rank('bogus'))),observe(invalid_echelon)],separators=(',',':'))
FUNCTIONS['binary_matrix_cache']=binary_matrix_cache

def binary_matrix_echelon_options(m,n,seed,algorithm,reduced,k,prepare):
    import json
    rows=m4ri_packed_fixture(m,n,seed)
    A=matrix(GF(2),m,n,[(row>>j)&1 for row in rows for j in range(n)])
    if prepare==1:A.echelonize()
    elif prepare==2:A.set_immutable()
    elif prepare==3:A.rank()
    def observe(call):
        try:return ['ok',call()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():
        value=A.echelonize(BINARY_ALGORITHMS[int(algorithm)],0,bool(reduced),k=int(k))
        return 'self' if value is A else 'None'
    result=observe(call)
    return json.dumps([result,[str(x) for x in A.list()],observe(lambda:int(A.rank('bogus')))],separators=(',',':'))
FUNCTIONS['binary_matrix_echelon_options']=binary_matrix_echelon_options

def m4ri_echelon_pattern(m,n,seed,pattern,algorithm,full,k):
    import json
    from m4ri_pinned import echelon_pattern
    return json.dumps(echelon_pattern(m,n,seed,pattern,algorithm,full,k),separators=(',',':'))
FUNCTIONS['m4ri_echelon_pattern']=m4ri_echelon_pattern

_binary_ordering=None
def binary_ordering_module():
    global _binary_ordering
    if _binary_ordering is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _binary_ordering=cython_import(str(Path(__file__).parent.parent/'binary_ordering_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return _binary_ordering

def binary_factorization(method,m,n,seed,algorithm,param,kind):
    import json
    from sage.matrix.matrix_mod2_dense import ple,pluq
    rows=m4ri_packed_fixture(m,n,seed)
    A=matrix(GF(2),m,n,[(row>>j)&1 for row in rows for j in range(n)])
    name=['standard','mmpf' if method else 'russian','naive','bogus'][int(algorithm)]
    value=int(param) if kind in (0,1) else ZZ(param)
    if algorithm<3 and (not m or not n or (algorithm==1 and 9<=param<=16)):
        # Project the bundled wrapper's copied result from its pinned dependency.
        # Empty native faults and oversized Russian panels differ across releases.
        from m4ri_pinned import factorization
        rank,P,Q,packed=factorization(m,n,seed,3*method+algorithm,int(param))
        B=matrix(GF(2),m,n,[(int(row)>>j)&1 for row in packed for j in range(n)])
    else:B,P,Q=(pluq if method else ple)(A,name,value)
    return json.dumps([[str(x) for x in B.list()],P,Q,[str(x) for x in A.list()],B is A,B.is_mutable()],separators=(',',':'))
FUNCTIONS['binary_factorization']=binary_factorization

def m4ri_factorization(m,n,seed,kind,k):
    import json
    from m4ri_pinned import factorization
    return json.dumps(factorization(m,n,seed,kind,k),separators=(',',':'))
FUNCTIONS['m4ri_factorization']=m4ri_factorization

def binary_factorization_conversion(method,algorithm,kind,value):
    import json
    from sage.matrix.matrix_mod2_dense import ple,pluq
    param=float(value)/2 if kind==0 else float('nan') if kind==1 else float('inf') if kind==2 else -float('inf') if kind==3 else str(value) if kind==4 else None if kind==5 else bool(value) if kind==6 else QQ(value)/2 if kind==7 else ZZ(value) if kind==8 else int(value)
    A=matrix(GF(2),[[0,1,0,1,0],[0,1,1,0,1],[0,0,1,1,1]])
    B,P,Q=(pluq if method else ple)(A,['standard','mmpf' if method else 'russian','naive','bogus'][int(algorithm)],param)
    return json.dumps([[str(x) for x in B.list()],P,Q],separators=(',',':'))
FUNCTIONS['binary_factorization_conversion']=binary_factorization_conversion

def binary_swap_conversion(axis,m,n,prepare,kind,value,position):
    import json
    rows=m4ri_packed_fixture(m,n,42)
    A=matrix(GF(2),m,n,[(r>>j)&1 for r in rows for j in range(n)])
    if prepare==1:A.rank()
    elif prepare==2:A.echelonize('m4ri')
    elif prepare==3:A.set_immutable()
    elif prepare==4:A=A.echelon_form('m4ri')
    v=int(value) if kind==0 else ZZ(value) if kind==1 else QQ(value)/2 if kind==2 else float(value)/2 if kind==3 else float('nan') if kind==4 else float('inf') if kind==5 else str(value) if kind==6 else None if kind==7 else bool(value)
    if kind==3 and v.is_integer():v=int(v) # documented JS numeric-index adapter
    args=(v,0) if position==0 else (0,v) if position==1 else (v,v)
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():
        # Older installed M4RI crashes valid swaps with an empty other axis.
        import operator
        i,j=(operator.index(x) for x in args)
        if all(-(1<<63)<=x<(1<<63) for x in (i,j)) and A.is_mutable() and ((axis and m==0 and 0<=i<n and 0<=j<n) or (not axis and n==0 and 0<=i<m and 0<=j<m)):
            from m4ri_pinned import empty_swap
            A._clear_cache();empty_swap(m,n,int(axis),i,j)
            return None
        return (A.swap_columns if axis else A.swap_rows)(*args)
    result=observe(call)
    return json.dumps([result,[str(x) for x in A.list()],observe(lambda:int(A.rank('bogus'))),A.is_immutable()],separators=(',',':'))
FUNCTIONS['binary_swap_conversion']=binary_swap_conversion

def binary_inverse(m,n,seed,pattern,prepare):
    import json
    rows=m4ri_packed_fixture(m,n,seed)
    if pattern==1:rows=[(r&(-1<<i))|(1<<i) for i,r in enumerate(rows)]
    elif pattern==2:rows=[0]*m
    A=matrix(GF(2),m,n,[(r>>j)&1 for r in rows for j in range(n)])
    if prepare==1:A.rank()
    elif prepare==2:A.set_immutable()
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():
        B=~A
        return [B.nrows(),B.ncols(),[str(x) for x in B.list()],B is A,B.is_mutable()]
    result=observe(call)
    return json.dumps([result,[str(x) for x in A.list()],observe(lambda:int(A.rank('bogus')))],separators=(',',':'))
FUNCTIONS['binary_inverse']=binary_inverse

def m4ri_inverse(n,seed,pattern,k):
    import json
    from m4ri_pinned import inverse
    return json.dumps(inverse(n,seed,pattern,k),separators=(',',':'))
FUNCTIONS['m4ri_inverse']=m4ri_inverse

def binary_permutation(axis,m,n,prepare,perm):
    import json
    rows=m4ri_packed_fixture(m,n,42)
    A=matrix(GF(2),m,n,[(r>>j)&1 for r in rows for j in range(n)])
    if prepare==1:A.rank()
    elif prepare==2:A.echelonize('m4ri')
    elif prepare==3:A.set_immutable()
    elif prepare==4:A=A.echelon_form('m4ri')
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():
        p=SymmetricGroup(len(perm))([int(x) for x in perm])
        if A.is_mutable() and ((axis and m==0) or (not axis and n==0)):
            # Bundled native swaps have an empty-row guard absent in installed M4RI.
            from m4ri_pinned import empty_swap
            A._clear_cache()
            for cycle in p.cycle_tuples():
                cycle=[x-1 for x in reversed(cycle)]
                for elt in cycle:
                    if not(0<=cycle[0]<(n if axis else m) and 0<=elt<(n if axis else m)):
                        raise IndexError('matrix '+('column' if axis else 'row')+' index out of range')
                    if cycle[0]!=elt:empty_swap(m,n,int(axis),cycle[0],elt)
            return None
        return (A.permute_columns if axis else A.permute_rows)(p)
    result=observe(call)
    return json.dumps([result,A.nrows(),A.ncols(),[str(x) for x in A.list()],observe(lambda:int(A.rank('bogus')))],separators=(',',':'))
FUNCTIONS['binary_permutation']=binary_permutation

def binary_solve(m,n,brows,bcols,seed,pattern,check,prepare):
    import json
    rows=m4ri_packed_fixture(m,n,seed)
    if pattern==1:rows=[0]*m
    elif pattern==2:rows=[(1<<i) if i<n else 0 for i in range(m)]
    A=matrix(GF(2),m,n,[(r>>j)&1 for r in rows for j in range(n)])
    bits=m4ri_packed_fixture(brows,bcols,seed+1)
    B=matrix(GF(2),brows,bcols,[(r>>j)&1 for r in bits for j in range(bcols)])
    if pattern==3 and brows==m:
        x=m4ri_packed_fixture(n,bcols,seed+2)
        B=A*matrix(GF(2),n,bcols,[(r>>j)&1 for r in x for j in range(bcols)])
    if prepare==1:A.rank()
    elif prepare==2:A.set_immutable();B.set_immutable()
    elif prepare==3:A.echelon_form()
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():
        # Installed 10.3 has no specialized binary _solve_right_general.
        # Dispatch the bundled method with the bundled public shape/rank checks.
        if A.nrows()!=B.nrows():raise ValueError('number of rows of self must equal number of rows of right-hand side')
        if A.is_square() and A.rank()==A.nrows():X=A._solve_right_nonsingular_square(B,check_rank=False)
        else:X=binary_solve_module()._solve_right_general(A,B,check=True if check==2 else bool(check))
        return [X.nrows(),X.ncols(),[str(v) for v in X.list()],X.is_mutable()]
    result=observe(call)
    return json.dumps([result,[str(x) for x in A.list()],[str(x) for x in B.list()],observe(lambda:int(A.rank('bogus'))),observe(lambda:[str(x) for x in A.echelon_form('bogus').list()])],separators=(',',':'))
FUNCTIONS['binary_solve']=binary_solve

def binary_kernel(m,n,seed,pattern,basis,algorithm,prepare):
    import json
    rows=m4ri_packed_fixture(m,n,seed)
    if pattern==1:rows=[0]*m
    elif pattern==2:rows=[(1<<(n-1)) if n else 0 for _ in range(m)]
    A=matrix(GF(2),m,n,[(r>>j)&1 for r in rows for j in range(n)])
    if prepare==1:A.rank()
    elif prepare==2:A.set_immutable()
    elif prepare==3:A.echelon_form()
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():
        if algorithm==5:
            # Bundled matrix2.pyx adds linbox to accepted names, then rejects GF(2).
            raise ValueError("'linbox' matrix kernel algorithm only available over the rationals, not over Finite Field of size 2")
        B=A.right_kernel_matrix(basis=['echelon','pivot','computed','default','bogus','LLL'][int(basis)],algorithm=['default','generic','pluq','bogus','flint','linbox','pari','padic'][int(algorithm)])
        return [B.nrows(),B.ncols(),[str(v) for v in B.list()],B.is_mutable()]
    result=observe(call)
    return json.dumps([result,[str(x) for x in A.list()],observe(lambda:int(A.rank('bogus'))),observe(lambda:[str(x) for x in A.echelon_form('bogus').list()])],separators=(',',':'))
FUNCTIONS['binary_kernel']=binary_kernel

def m4ri_solve(m,n,bc,seed,pattern,cutoff,check):
    import json
    from m4ri_pinned import solve
    return json.dumps(solve(m,n,bc,seed,pattern,cutoff,check),separators=(',',':'))
FUNCTIONS['m4ri_solve']=m4ri_solve

def m4ri_kernel(m,n,seed,pattern,cutoff):
    import json
    from m4ri_pinned import kernel
    return json.dumps(kernel(m,n,seed,pattern,cutoff),separators=(',',':'))
FUNCTIONS['m4ri_kernel']=m4ri_kernel

_binary_solve=None
def binary_solve_module():
    global _binary_solve
    if _binary_solve is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _binary_solve=cython_import(str(Path(__file__).parent.parent/'binary_solve_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return _binary_solve

def binary_kernel_proof(m,n,basis,algorithm,proof):
    import json
    A=matrix(GF(2),m,n)
    value=[None,True,False,0,1,2,'bad'][int(proof)]
    if algorithm==2:raise ValueError("'linbox' matrix kernel algorithm only available over the rationals, not over Finite Field of size 2")
    B=A.right_kernel_matrix(basis=['default','computed','bogus','LLL'][int(basis)],algorithm=['default','bogus','linbox'][int(algorithm)],proof=value)
    return json.dumps([B.nrows(),B.ncols(),[str(x) for x in B.list()],B.is_mutable()],separators=(',',':'))
FUNCTIONS['binary_kernel_proof']=binary_kernel_proof

def binary_basic(m,n,r,c,seed,operation,prepare):
    import json
    def make(m,n,seed):
        rows=m4ri_packed_fixture(m,n,seed)
        return matrix(GF(2),m,n,[(v>>j)&1 for v in rows for j in range(n)])
    A=make(m,n,seed);B=make(r,c,seed+1)
    if prepare==1:A.rank()
    elif prepare==2:A.set_immutable()
    elif prepare==3:A.echelon_form()
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():
        if operation==0:X=A+B
        elif operation==1:X=A-B
        elif operation==2:X=A.transpose()
        elif operation==3:X=A.augment(B)
        elif operation==4:return int(A.determinant())
        elif operation==5:X=-A
        elif operation==6:X=A.__copy__()
        elif operation==7:return [[str(v) for v in x] for x in A.columns()]
        elif operation==8:return [str(v) for v in A.list()]
        elif operation==9:return A==B
        elif operation==10:return str(A)
        return [X.nrows(),X.ncols(),[str(v) for v in X.list()],X.is_mutable()]
    result=observe(call)
    return json.dumps([result,[str(v) for v in A.list()],observe(lambda:int(A.rank('bogus'))),observe(lambda:[str(v) for v in A.echelon_form('bogus').list()])],separators=(',',':'))
FUNCTIONS['binary_basic']=binary_basic

def binary_columns_state(m,n,copy,operation,prepare):
    import json
    rows=m4ri_packed_fixture(m,n,42)
    A=matrix(GF(2),m,n,[(r>>j)&1 for r in rows for j in range(n)])
    if prepare==1:A.rank()
    elif prepare==2:A.echelon_form()
    C=A.columns() if copy==2 else A.columns(copy=bool(copy))
    D=A.columns(copy=False)
    before=[C is D, bool(C and D and C[0] is D[0])]
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():
        if operation==1 and C and m:C[0][0]=1-int(C[0][0])
        elif operation==2 and C:C[0]=[0]*m
        elif operation==3:C.append([1]*m)
        elif operation==4 and C:C.pop()
        elif operation==5 and m and n:A[0,0]=1-int(A[0,0])
        elif operation==6:A[m,n]=1
        elif operation==7:A._clear_cache()
        return None
    result=observe(call);E=A.columns(copy=False);F=A.columns()
    return json.dumps([before,result,[E is D,F is E,bool(E and F and E[0] is F[0])],[[str(v) for v in row] for row in C],[[str(v) for v in row] for row in E],[str(v) for v in A.list()],observe(lambda:int(A.rank('bogus')))],separators=(',',':'))
FUNCTIONS['binary_columns_state']=binary_columns_state

def m4ri_basic(m,n,right,seed,op):
    import json
    from m4ri_pinned import basic
    return json.dumps(basic(m,n,right,seed,op),separators=(',',':'))
FUNCTIONS['m4ri_basic']=m4ri_basic

_binary_format=None
def binary_format_module():
    global _binary_format
    if _binary_format is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _binary_format=cython_import(str(Path(__file__).parent.parent/'binary_format_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return _binary_format

def binary_format(m,n,seed,mapping,zero,one,minus):
    import json
    rows=m4ri_packed_fixture(m,n,seed)
    A=matrix(GF(2),m,n,[(r>>j)&1 for r in rows for j in range(n)])
    texts=[None,'','.', 'zero','1','one','🙂','a\nb']
    counter=0
    def callback(x):
        nonlocal counter
        counter+=1
        return str((counter+int(x))%3)
    mp=None if mapping==0 else {} if mapping==1 else {GF(2)(0):'xx'} if mapping==2 else {GF(2)(1):'y'} if mapping==3 else {GF(2)(0):'xx',GF(2)(1):'y'} if mapping==4 else callback
    try:result=['ok',A.str(mp,texts[int(zero)],texts[int(one)],texts[int(minus)])]
    except Exception as e:result=['error',type(e).__name__,str(e)]
    state=sorted([[str(k),v] for k,v in mp.items()]) if isinstance(mp,dict) else counter
    return json.dumps([result,state],separators=(',',':'),ensure_ascii=False)
FUNCTIONS['binary_format']=binary_format

def binary_subdivision(m,n,rows,cols,operation,prepare):
    import json
    bits=m4ri_packed_fixture(m,n,42)
    A=matrix(GF(2),m,n,[(r>>j)&1 for r in bits for j in range(n)])
    if prepare==1:A.rank()
    elif prepare==2:A.echelon_form()
    elif prepare==3:A.set_immutable()
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    result=observe(lambda:A.subdivide(list(rows),list(cols)))
    if operation==1:A=A.__copy__()
    elif operation==2:
        # Bundled transpose preserves subdivisions; installed 10.3 omits this step.
        r,c=A.subdivisions();A=A.transpose()
        if r or c:A.subdivide(c,r)
    elif operation==3:A=-A
    elif operation==4:A=A+A
    elif operation==5:A=A.augment(A,subdivide=True)
    elif operation==6:A=A.augment(A,subdivide=False)
    elif operation==7:A=A.echelon_form()
    elif operation==8:
        r,c=A.get_subdivisions();r.append(0);c.append(0)
    elif operation==9:result=observe(lambda:A.subdivide())
    elif operation==10:result=observe(lambda:A.subdivide([1],None))
    return json.dumps([result,[[str(x) for x in row] for row in A.subdivisions()],A.is_mutable(),observe(lambda:binary_format_module().str(A)),observe(lambda:A.str(zero='.')),observe(lambda:int(A.rank('bogus')))],separators=(',',':'))
FUNCTIONS['binary_subdivision']=binary_subdivision

def binary_subdivision_conversion(kind,value,side,prepare):
    import json
    A=matrix(GF(2),[[1,0,1],[0,1,0]]);A.subdivide(1,1)
    if prepare==1:A.rank()
    elif prepare==2:A.set_immutable()
    elif prepare==3:A.echelon_form()
    def convert(x):
        return float(x)/2 if kind==0 else QQ(x)/2 if kind==1 else ZZ(x) if kind==2 else int(x) if kind==3 else str(x) if kind==4 else None if kind==5 else bool(x) if kind==6 else float('nan') if kind==7 else float('inf')
    v=convert(value)
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():
        if side==0:A.subdivide(v,None)
        elif side==1:A.subdivide(None,v)
        elif side==2:A.subdivide([v,convert(0)],[])
        elif side==3:A.subdivide([],[v,convert(2)])
        elif side==4:A.subdivide([v,convert(2)],['bad',0])
        else:A.subdivide([v,0],[])
        return None
    return json.dumps([observe(call),[[str(x) for x in r] for r in A.subdivisions()],observe(lambda:int(A.rank('bogus')))],separators=(',',':'))
FUNCTIONS['binary_subdivision_conversion']=binary_subdivision_conversion

def binary_format_invalid(m,n,kind):
    import json
    A=matrix(GF(2),m,n,[(i+j)%2 for i in range(m) for j in range(n)])
    mp=[{0:None},{1:1},lambda x:None,lambda x:1,lambda x:['x'],None,None,1,'0',[],[0],[1],{0:{'x':1}}][int(kind)]
    zero=1 if kind==5 else None;one=False if kind==6 else None
    return A.str(mp,zero,one)
FUNCTIONS['binary_format_invalid']=binary_format_invalid

def binary_augmentation_subdivisions(m,n,c,left_rows,left_cols,right_rows,right_cols,flag,operation):
    import json
    def make(m,n,seed):
        rows=m4ri_packed_fixture(m,n,seed)
        return matrix(GF(2),m,n,[(r>>j)&1 for r in rows for j in range(n)])
    A=make(m,n,42);B=make(m,c,43)
    A.subdivide(list(left_rows),list(left_cols));B.subdivide(list(right_rows),list(right_cols))
    A.set_immutable();B.set_immutable()
    C=A.augment(B,subdivide=bool(flag))
    if operation==1:C=C.__copy__()
    elif operation==2:
        r,s=C.subdivisions();C=C.transpose()
        if r or s:C.subdivide(s,r) # bundled transpose metadata fix absent in installed release
    elif operation==3:C=-C
    elif operation==4:C=C.submatrix(0,0,C.nrows(),C.ncols())
    return json.dumps([C.nrows(),C.ncols(),[[str(x) for x in r] for r in C.subdivisions()],C.str(),C.str(zero='.'),C.is_mutable(),[[str(x) for x in r] for r in A.subdivisions()],[[str(x) for x in r] for r in B.subdivisions()]],separators=(',',':'))
FUNCTIONS['binary_augmentation_subdivisions']=binary_augmentation_subdivisions

# Cached versus freshly extracted binary rows, including mutation and invalidation.
def binary_row_state(m,n,index,from_list,operation,prepare):
    import json
    bits=m4ri_packed_fixture(m,n,42)
    A=matrix(GF(2),m,n,[(r>>j)&1 for r in bits for j in range(n)])
    if prepare==1:A.rank()
    elif prepare==2:A.echelon_form()
    elif prepare==3:A.set_immutable()
    C=A.row(index,from_list=bool(from_list));D=A.row(index,from_list=True)
    before=C is D
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():
        if operation==1 and n:C[0]=1-int(C[0])
        elif operation==2 and m and n:A[0,0]=1-int(A[0,0])
        elif operation==3:A[m,n]=1
        elif operation==4:A._clear_cache()
        elif operation==5:A.swap_rows(0,m-1)
        elif operation==6:A.subdivide([0],[0])
        elif operation==7:A.echelonize()
        return None
    result=observe(call);E=A.row(index,from_list=bool(from_list));F=A.row(index,from_list=True)
    return json.dumps([before,result,[C is E,D is F,E is F],[str(v) for v in C],[str(v) for v in E],[str(v) for v in A.list()],observe(lambda:int(A.rank('bogus')))],separators=(',',':'))
FUNCTIONS['binary_row_state']=binary_row_state

def binary_entry_value(kind,value):
    return (int(value) if kind==0 else ZZ(value) if kind==1 else QQ(value)/3 if kind==2 else QQ(value)/2 if kind==3 else (int(float(value)/2) if (float(value)/2).is_integer() else float(value)/2) if kind==4 else bool(value) if kind==5 else str(value) if kind==6 else None if kind==7 else GF(2)(value) if kind==8 else GF(3)(value) if kind==9 else Mod(value,4) if kind==10 else float('nan') if kind==11 else float('inf'))

def binary_constructor(m,n,layout,kind,value):
    import json
    v=binary_entry_value(kind,value)
    flat=[v]*(m*n);nested=[[v]*n for _ in range(m)]
    entries=None if layout<2 else v if layout==2 else flat if layout==3 else nested if layout==4 else flat[:-1] if layout==5 else flat+[v] if layout==6 else nested[:-1] if layout==7 else nested+[[v]*n] if layout==8 else [row[:-1] for row in nested] if layout==9 else [row+[v] for row in nested] if layout==10 else [nested[0] if nested else [],v]
    # Bundled args.pyx recognizes strings as scalars; installed 10.3 predates
    # issue 34821. These fixtures use valid numeric strings, so conversion may
    # precede flat-list validation without changing any competing error.
    if kind==6 and layout==2:
        if m!=n:raise TypeError('nonzero scalar matrix must be square')
        scalar=GF(2)(v)
        entries=[scalar if i==j else GF(2)(0) for i in range(m) for j in range(n)]
    elif kind==6 and n==1 and isinstance(entries,list) and entries and isinstance(entries[0],str):
        entries=[GF(2)(x) for x in entries]
    A=MatrixSpace(GF(2),int(m),int(n))(entries)
    return json.dumps([A.nrows(),A.ncols(),[str(x) for x in A.list()],A.is_mutable()],separators=(',',':'))
FUNCTIONS['binary_constructor']=binary_constructor

def binary_set_entry(m,n,index,kind,value,prepare):
    import json
    A=matrix(GF(2),m,n,[(i+1)%2 for i in range(m*n)])
    if prepare==1:A.rank()
    elif prepare==2:A.echelon_form()
    elif prepare==3:A.set_immutable()
    v=binary_entry_value(kind,value)
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    def call():A[index,index]=v
    return json.dumps([observe(call),[str(x) for x in A.list()],observe(lambda:int(A.rank('bogus')))],separators=(',',':'))
FUNCTIONS['binary_set_entry']=binary_set_entry

def binary_dimensions(kind,m,n):
    import json
    def cv(v):
        if kind==0:return int(v)
        if kind==1:return int(v)
        if kind==2:return ZZ(v)
        if kind==3:return QQ(v)/2
        if kind==4:return float(v)/2
        if kind==5:return float('nan')
        return float('inf')
    A=MatrixSpace(GF(2),cv(m),cv(n))()
    return json.dumps([A.nrows(),A.ncols(),[str(x) for x in A.list()]],separators=(',',':'))
FUNCTIONS['binary_dimensions']=binary_dimensions

def binary_random_state(m,n,seed,density,nonzero,prepare,repeats):
    import json
    from sage.misc.randstate import random as sage_random
    A=matrix(GF(2),m,n,[(i+1)%2 for i in range(m*n)])
    if prepare==1:A.rank()
    elif prepare==2:A.echelon_form()
    elif prepare==3:A.set_immutable()
    densities=[1.0,-1.0,0.0,1/100,1/3,0.5,0.999,1.0,2.0,float('nan'),float('inf'),-float('inf')]
    set_random_seed(seed)
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    states=[]
    for _ in range(repeats):
        result=observe(lambda:A.randomize(nonzero=bool(nonzero)) if density==0 else A.randomize(densities[int(density)],bool(nonzero)))
        states.append([result,[str(v) for v in A.list()]])
    return json.dumps([states,[str(sage_random()) for _ in range(3)],observe(lambda:int(A.rank('bogus'))),A.is_mutable()],separators=(',',':'))
FUNCTIONS['binary_random_state']=binary_random_state

def binary_factory(op,m,n,kind,seed,density):
    import json
    from sage.misc.randstate import random as sage_random
    def cv(v):return int(v) if kind<2 else ZZ(v) if kind==2 else QQ(v)/2 if kind==3 else float(v)/2
    rows=cv(m);cols=None if n==-999 else cv(n)
    densities=[None,-1.0,0.0,1/100,1/3,0.5,0.999,1.0,2.0,float('nan'),float('inf'),-float('inf')]
    set_random_seed(seed)
    if op==0:A=zero_matrix(GF(2),rows,cols)
    elif op==1:A=identity_matrix(GF(2),rows)
    else:A=random_matrix(GF(2),rows,cols,density=densities[int(density)])
    return json.dumps([A.nrows(),A.ncols(),[str(v) for v in A.list()],A.is_mutable(),[str(sage_random()) for _ in range(3)]],separators=(',',':'))
FUNCTIONS['binary_factory']=binary_factory

def binary_from_entries(m,n,layout,kind,value):
    import json
    v=binary_entry_value(kind,value);entries=[[v]*n for _ in range(m)]
    if layout==1:entries=entries[:-1]
    elif layout==2:entries.append([v]*n)
    elif layout==3 and entries:entries[-1]=entries[-1][:-1]
    elif layout==4 and entries:entries[-1].append(v)
    elif layout==5 and entries:entries[0]=[]
    A=matrix(GF(2),entries)
    return json.dumps([A.nrows(),A.ncols(),[str(v) for v in A.list()],A.is_mutable()],separators=(',',':'))
FUNCTIONS['binary_from_entries']=binary_from_entries


BINARY_DENSITY_STRINGS=['bad', '', '0.5', ' 1_0e-2 ', 'inf', '-Infinity', 'nan', 'NaN', '0x1', '1e9999', '-1e9999', '-0.0', '1__0', '1_.0', '١.٥', '１.０', '\xa00.5\u2003', '0.5\x1c', '\x1c0.5', '0.5\ufeff', '\ufeff0.5', '1\x00', "a'b", 'a"b', 'a\'"b', '🙂', '\ud800']
def binary_density_value(kind):
    values=[None,False,True,2**2048,ZZ(2)**2048,-(2**2048),-(ZZ(2)**2048),QQ(1)/3,QQ(2)**2048,-QQ(1)/(ZZ(2)**2048),RealField(80)('0.333333333333333333333333333'),RealField(53)('0.5'),GF(2)(1),GF(3)(2),Mod(0,4),[],{},b'0.5',b'bad',GF(2)(1)]
    return values[int(kind)] if kind<20 else BINARY_DENSITY_STRINGS[int(kind)-20]

def binary_density_input(m,n,kind,prepare,nonzero):
    import json
    from sage.misc.randstate import random as sage_random
    density=binary_density_value(kind)
    A=matrix(GF(2),m,n,[(i+1)%2 for i in range(m*n)])
    if prepare==1:A.rank()
    elif prepare==2:A.echelon_form()
    elif prepare==3:A.set_immutable()
    set_random_seed(42)
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    result=observe(lambda:A.randomize(density,bool(nonzero)))
    return json.dumps([result,[str(x) for x in A.list()],[str(sage_random()) for _ in range(3)],observe(lambda:int(A.rank('bogus')))],separators=(',',':'),ensure_ascii=False)
FUNCTIONS['binary_density_input']=binary_density_input

def python_float_text(codes,as_bytes):
    import json
    from python_float_native import text_float
    return json.dumps(text_float(tuple(int(c) for c in codes), bool(as_bytes)),separators=(',',':'),ensure_ascii=False)
FUNCTIONS['python_float_text']=python_float_text


def binary_predicates(m,n,mask,certificate,prepare):
    import json
    A=matrix(GF(2),m,n,[(mask>>i)&1 for i in range(m*n)])
    if prepare==1:A.rank()
    elif prepare==2:A.set_immutable()
    def observe(f):
        try:return ['ok',f()]
        except Exception as e:return ['error',type(e).__name__,str(e)]
    result=binary_ordering_module().is_Gamma_free(A,bool(certificate))
    return json.dumps([result,A==A.__copy__(),A==A.transpose(),A==matrix(GF(2),m,n,0),[str(v) for v in A.list()],A.is_mutable(),observe(lambda:int(A.rank('bogus')))],separators=(',',':'))
FUNCTIONS['binary_predicates']=binary_predicates

def binary_png_data(m,n,mask,mode):
    import json,os,tempfile
    from PIL import Image
    from sage.matrix.matrix_mod2_dense import to_png,from_png
    A=matrix(GF(2),m,n,[(mask>>i)&1 for i in range(m*n)])
    fd,path=tempfile.mkstemp(suffix='.png');os.close(fd)
    try:
        if mode==0:
            # Bundled to_png lowercases this message; installed 10.3 capitalizes it.
            if not m or not n:raise TypeError(f"cannot write image with dimensions {n} x {m}")
            to_png(A,path)
        else:
            image=Image.new('P',(n,m));image.putpalette([0,0,0,255,255,255]+[0]*762)
            image.putdata([1-int(x) for x in A.list()]);image.save(path,bits=1)
        with Image.open(path) as image:
            pixels=list(image.convert('L').getdata());width,height=image.size
        B=from_png(path)
        return json.dumps([width,height,pixels,[str(x) for x in B.list()],B.is_mutable()],separators=(',',':'))
    finally:os.unlink(path)
FUNCTIONS['binary_png_data']=binary_png_data


def binary_density_factory(m,n,kind):
    import json
    from sage.misc.randstate import random as sage_random
    density=binary_density_value(kind)
    set_random_seed(42)
    try:
        A=random_matrix(GF(2),m,n,density=density)
        result=['ok',[str(x) for x in A.list()],A.is_mutable()]
    except Exception as e:result=['error',type(e).__name__,str(e)]
    return json.dumps([result,[str(sage_random()) for _ in range(3)]],separators=(',',':'),ensure_ascii=False)
FUNCTIONS['binary_density_factory']=binary_density_factory
