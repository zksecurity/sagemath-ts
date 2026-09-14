# Original bundled Sage method, unavailable in installed Sage 10.3.
# Only dedentation and a typed self argument adapt the method into a Cython function.
from sage.matrix.matrix_mod2_dense cimport Matrix_mod2_dense
from sage.libs.m4ri cimport mzd_read_bit

def doubly_lexical_ordering(Matrix_mod2_dense self, inplace=False):
    r"""
    Return a doubly lexical ordering of the matrix.

    A doubly lexical ordering of a matrix is an ordering of the rows
    and of the columns of the matrix so that both the rows and the
    columns, as vectors, are lexically increasing. See [Lub1987]_.
    A lexical ordering of vectors is the standard dictionary ordering,
    except that vectors will be read from highest to lowest coordinate.
    Thus row vectors will be compared from right to left, and column
    vectors from bottom to top.

    INPUT:

    - ``inplace`` -- boolean (default: ``False``); using ``inplace=True``
      will permute the rows and columns of the current matrix
      according to a doubly lexical ordering; this will modify the matrix

    OUTPUT:

    A pair ``(row_ordering, col_ordering)`` of
    :class:`~sage.groups.perm_gps.constructor.PermutationGroupElement`
    that represents a doubly lexical ordering of the rows or columns.

    .. SEEALSO::

        :meth:`~sage.matrix.matrix2.Matrix.permutation_normal_form`;
        a similar matrix normal form

    ALGORITHM:

    The algorithm is adapted from section 3 of [HAM1985]_. The time
    complexity of this algorithm is `O(n \cdot m^2)` for a `n \times m`
    matrix.

    EXAMPLES::

        sage: A = Matrix(GF(2), [[0, 1],
        ....:                    [1, 0]])
        sage: r, c = A.doubly_lexical_ordering()
        sage: r
        (1,2)
        sage: c
        ()
        sage: A.permute_rows_and_columns(r, c); A
        [1 0]
        [0 1]

    ::

        sage: A = Matrix(GF(2), [[0, 1],
        ....:                    [1, 0]])
        sage: r, c = A.doubly_lexical_ordering(inplace=True); A
        [1 0]
        [0 1]

    TESTS:

    This algorithm works correctly for the matrix in
    Example 3.7 in [HAM1985]_::

        sage: A = Matrix(GF(2), [[1, 1, 0, 0, 0, 0, 0],
        ....:                    [1, 1, 0, 0, 0, 0, 0],
        ....:                    [1, 1, 0, 1, 0, 0, 0],
        ....:                    [0, 0, 1, 1, 0, 0, 0],
        ....:                    [0, 1, 1, 1, 1, 0, 0],
        ....:                    [0, 0, 0, 0, 0, 1, 1],
        ....:                    [0, 0, 0, 0, 0, 1, 1],
        ....:                    [0, 0, 0, 0, 1, 1, 1],
        ....:                    [0, 0, 0, 1, 1, 1, 0]])
        sage: r, c = A.doubly_lexical_ordering()
        sage: B = A.with_permuted_rows_and_columns(r, c)
        sage: for i in range(B.ncols()):
        ....:     for j in range(i):
        ....:         for k in reversed(range(B.nrows())):
        ....:             assert B[k][j] <= B[k][i]
        ....:             if B[k][j] < B[k][i]:
        ....:                 break
        sage: for i in range(B.nrows()):
        ....:     for j in range(i):
        ....:         for k in reversed(range(B.ncols())):
        ....:             assert B[j][k] <= B[i][k]
        ....:             if B[j][k] < B[i][k]:
        ....:                 break
        sage: r, c = A.doubly_lexical_ordering(inplace=True)
        sage: A == B
        True

    An immutable matrix calling with ``inplace=True`` will raise an error::

        sage: A = Matrix(GF(2), [[0, 1], [1, 0]], immutable=True)
        sage: r, c = A.doubly_lexical_ordering(inplace=True)
        Traceback (most recent call last):
        ...
        TypeError: this matrix is immutable;
         use inplace=False or apply to a mutable copy.

    The algorithm works collectly for a matrix with nrows=0 or ncols=0::

        sage: A = Matrix(GF(2), 0, 2, [])
        sage: A.doubly_lexical_ordering()
        ((), ())
        sage: B = Matrix(GF(2), 2, 0, [])
        sage: B.doubly_lexical_ordering()
        ((), ())
    """
    if inplace and self.is_immutable():
        raise TypeError("this matrix is immutable;"
                        " use inplace=False or apply to a mutable copy.")

    from sage.groups.perm_gps.permgroup_named import SymmetricGroup
    from sage.groups.perm_gps.permgroup_element import make_permgroup_element_v2
    symmetric_group_nrows = SymmetricGroup(self._nrows)
    symmetric_group_ncols = SymmetricGroup(self._ncols)

    if self._nrows == 0 or self._ncols == 0:
        return symmetric_group_nrows.one(), symmetric_group_ncols.one()

    cdef list partition_rows = [False for _ in range(self._nrows - 1)]
    cdef int partition_num = 1
    cdef list row_swapped = list(range(1, self._nrows + 1))
    cdef list col_swapped = list(range(1, self._ncols + 1))

    cdef Matrix_mod2_dense A = self if inplace else self.__copy__()

    cdef int i, col, row, partition_i, partition_start, partition_end
    cdef int largest_col, row_start, row_end
    for i in reversed(range(1, A._ncols + 1)):
        # count 1 for each partition and column
        count1 = [[0]*partition_num for _ in range(i)]
        for col in range(i):
            parition_i = 0
            for row in reversed(range(1, A._nrows)):
                count1[col][parition_i] += mzd_read_bit(A._entries, row, col)
                if partition_rows[row - 1]:
                    parition_i += 1
            count1[col][parition_i] += mzd_read_bit(A._entries, 0, col)  # special case of row == 0

        # calculate largest_col = col s.t. count1[col] is lexicographically largest (0 <= col < i)
        _, largest_col = max((c, i) for i, c in enumerate(count1))

        # We refine each partition of rows according to the value of A[:][largest_col].
        # and also move down rows that satisfy A[row][largest_col] = 1 in each partition.
        partition_start = 0
        for _ in range(partition_num):
            partition_end = partition_start
            while partition_end < A._nrows - 1 and not partition_rows[partition_end]:
                partition_end += 1
            row_start = partition_start
            row_end = partition_end
            while row_start < row_end:
                while row_start < row_end and not mzd_read_bit(A._entries, row_start, largest_col):
                    row_start += 1
                while row_start < row_end and mzd_read_bit(A._entries, row_end, largest_col):
                    row_end -= 1
                if row_start < row_end:  # swap row
                    A.swap_rows_c(row_start, row_end)
                    row_swapped[row_start], row_swapped[row_end] = row_swapped[row_end], row_swapped[row_start]
            partition_start = partition_end + 1  # for next partition

        for row in range(A._nrows - 1):
            if mzd_read_bit(A._entries, row, largest_col) != mzd_read_bit(A._entries, row + 1, largest_col):
                if not partition_rows[row]:
                    partition_rows[row] = True
                    partition_num += 1

        # swap column
        A.swap_columns_c(largest_col, i - 1)
        col_swapped[largest_col], col_swapped[i - 1] = col_swapped[i - 1], col_swapped[largest_col]

    row_ordering = make_permgroup_element_v2(symmetric_group_nrows, row_swapped, symmetric_group_nrows.domain())
    col_ordering = make_permgroup_element_v2(symmetric_group_ncols, col_swapped, symmetric_group_ncols.domain())

    return row_ordering, col_ordering

def is_Gamma_free(Matrix_mod2_dense self, certificate=False):
    r"""
    Return True if the matrix is `\Gamma`-free.

    A matrix is `\Gamma`-free if it does not contain a 2x2 submatrix
    of the form:

    .. MATH::

        \begin{pmatrix}
            1 & 1 \\
            1 & 0
        \end{pmatrix}

    INPUT:

    - ``certificate`` -- boolean (default: ``False``); whether to return a
      certificate for no-answers (see OUTPUT section)

    OUTPUT:

    When ``certificate`` is set to ``False`` (default) this method only
    returns ``True`` or ``False`` answers. When ``certificate`` is set to
    ``True``, the method either returns ``(True, None)`` or ``(False,
    (r1, c1, r2, c2))`` where ``r1``, ``r2``-th rows and ``c1``,
    ``c2``-th columns of the matrix constitute the `\Gamma`-submatrix.

    ALGORITHM:

    For each 1 entry, the algorithm finds the next 1 in the same row and
    the next 1 in the same column, and check the 2x2 submatrix that contains
    these entries forms `\Gamma` submatrix. The time complexity of
    this algorithm is `O(n \cdot m)` for a `n \times m` matrix.

    EXAMPLES::

        sage: A = Matrix(GF(2), [[1, 1],
        ....:                    [0, 0]])
        sage: A.is_Gamma_free()
        True
        sage: B = Matrix(GF(2), [[1, 1],
        ....:                    [1, 0]])
        sage: B.is_Gamma_free(certificate=True)
        (False, (0, 0, 1, 1))

    TESTS:

    The algorithm works collectly for larger matrices::

        sage: A = Matrix(GF(2), [[1, 0, 1],
        ....:                    [0, 0, 0],
        ....:                    [1, 0, 0]])
        sage: A.is_Gamma_free(certificate=True)
        (False, (0, 0, 2, 2))
        sage: B = Matrix(GF(2), [[1, 0, 1],
        ....:                    [0, 0, 0],
        ....:                    [1, 0, 1]])
        sage: B.is_Gamma_free(certificate=True)
        (True, None)
    """
    cdef int i, j, i_bottom, j_right

    for i in range(self._nrows):
        j = 0
        while j < self._ncols:
            if mzd_read_bit(self._entries, i, j):  # if A[i][j] == 1
                # find the next 1 in the row
                j_right = j + 1
                while j_right < self._ncols and not mzd_read_bit(self._entries, i, j_right):
                    j_right += 1
                if j_right < self._ncols:
                    # find the next 1 in the column
                    i_bottom = i + 1
                    while i_bottom < self._nrows and not mzd_read_bit(self._entries, i_bottom, j):
                        i_bottom += 1
                    if i_bottom < self._nrows and not mzd_read_bit(self._entries, i_bottom, j_right):
                        # A[i_bottom][j_right] == 0
                        if certificate:
                            return False, (i, j, i_bottom, j_right)
                        else:
                            return False
                j = j_right
            else:
                j += 1

    if certificate:
        return True, None
    else:
        return True
