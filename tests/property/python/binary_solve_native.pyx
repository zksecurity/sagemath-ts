# Bundled Sage solve method, absent in installed Sage 10.3.
# Only dedentation, typed self and the missing old .pxd declaration are adapted.
from sage.matrix.matrix_mod2_dense cimport Matrix_mod2_dense
from sage.libs.m4ri cimport mzd_t, rci_t, mzd_init, mzd_copy, mzd_free
from cysignals.signals cimport sig_on, sig_off
cdef extern from "m4ri/solve.h":
    int mzd_solve_left(mzd_t *, mzd_t *, int, int)

def _solve_right_general(Matrix_mod2_dense self, B, check=True):
    """
    Solve the matrix equation AX = B for X using the M4RI library.

    INPUT:

    - ``B`` -- a matrix
    - ``check`` -- boolean (default: ``True``); whether to check if the
      matrix equation has a solution

    EXAMPLES::

        sage: A = matrix(GF(2), [[1, 0], [0, 1], [1, 1]])
        sage: A.solve_right(vector([1, 1, 0]))
        (1, 1)
        sage: A.solve_right(vector([1, 1, 1]))
        Traceback (most recent call last):
        ...
        ValueError: matrix equation has no solutions

    TESTS::

        sage: n = 128
        sage: m = 128
        sage: A = random_matrix(GF(2), n, m)
        sage: B = A * random_vector(GF(2), m)
        sage: A * A.solve_right(B) == B
        True
        sage: m = 64
        sage: A = random_matrix(GF(2), n, m)
        sage: B = A * random_vector(GF(2), m)
        sage: A * A.solve_right(B) == B
        True
        sage: m = 256
        sage: A = random_matrix(GF(2), n, m)
        sage: B = A * random_vector(GF(2), m)
        sage: A * A.solve_right(B) == B
        True
        sage: matrix(GF(2), 2, 0).solve_right(matrix(GF(2), 2, 2)) == matrix(GF(2), 0, 2)
        True
        sage: matrix(GF(2), 2, 0).solve_right(matrix(GF(2), 2, 2) + 1)
        Traceback (most recent call last):
        ...
        ValueError: matrix equation has no solutions
        sage: matrix(GF(2), 2, 0).solve_right(matrix(GF(2), 2, 2) + 1, check=False) == matrix(GF(2), 0, 2)
        True
        sage: matrix(GF(2), 2, 2).solve_right(matrix(GF(2), 2, 0)) == matrix(GF(2), 2, 0)
        True

    Check that it can be interrupted::

        sage: set_random_seed(12345)
        sage: n, m = 20000, 19968
        sage: A = random_matrix(GF(2), n, m)
        sage: x = random_vector(GF(2), m)
        sage: B = A*x
        sage: from sage.doctest.util import ensure_interruptible_after
        sage: with ensure_interruptible_after(0.5): sol = A.solve_right(B)
    """
    cdef mzd_t *B_entries = (<Matrix_mod2_dense>B)._entries

    cdef Matrix_mod2_dense X  # the solution
    X = self.new_matrix(nrows=self._entries.ncols, ncols=B_entries.ncols)
    if self._entries.ncols == 0 or B_entries.ncols == 0:
        # special case: empty matrix
        if check and B != 0:
            raise ValueError("matrix equation has no solutions")
        return X
    cdef rci_t rows = self._entries.nrows
    if self._entries.nrows < self._entries.ncols:
        rows = self._entries.ncols  # mzd_solve_left requires ncols <= nrows

    cdef mzd_t *lhs = mzd_init(rows, self._entries.ncols)
    mzd_copy(lhs, self._entries)
    cdef mzd_t *rhs = mzd_init(rows, B_entries.ncols)
    mzd_copy(rhs, B_entries)

    cdef int ret
    try:
        sig_on()
        # although it is called mzd_solve_left, it does the same thing as solve_right
        ret = mzd_solve_left(lhs, rhs, 0, check)
        sig_off()

        if ret == 0:
            # solution is placed in rhs
            rhs.nrows = self._entries.ncols
            mzd_copy(X._entries, rhs)
            return X
        else:
            raise ValueError("matrix equation has no solutions")
    finally:
        mzd_free(lhs)
        mzd_free(rhs)

