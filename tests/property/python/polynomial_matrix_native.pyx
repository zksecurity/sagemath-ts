from sage.matrix.matrix_gf2e_dense cimport Matrix_gf2e_dense
from sage.libs.m4rie cimport mzed_read_elem

def native_binary_matrix_rows(Matrix_gf2e_dense a):
    # matrix_gf2e_dense.get_unsafe reads these exact words before calling the
    # unavailable Givaro _cache on a PARI-backed coefficient field.
    return [[mzed_read_elem(a._entries,i,j) for j in range(a._ncols)] for i in range(a._nrows)]
