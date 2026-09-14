# Bundled str method; only dedentation and a typed self adapt this function.
from sage.matrix.matrix_mod2_dense cimport Matrix_mod2_dense
from sage.libs.m4ri cimport mzd_read_bit
import sage.matrix.matrix_dense as matrix_dense

def str(Matrix_mod2_dense self, rep_mapping=None, zero=None, plus_one=None, minus_one=None,
        *, unicode=False, shape=None, character_art=False,
        left_border=None, right_border=None,
        top_border=None, bottom_border=None):
    r"""
    Return a nice string representation of the matrix.

    INPUT:

    - ``rep_mapping`` -- dictionary or callable used to override
      the usual representation of elements.  For a dictionary,
      keys should be elements of the base ring and values the
      desired string representation.

    - ``zero`` -- string (default: ``None``); if not ``None`` use
      the value of ``zero`` as the representation of the zero
      element.

    - ``plus_one`` -- string (default: ``None``); if not ``None``
      use the value of ``plus_one`` as the representation of the
      one element.

    - ``minus_one`` -- ignored.  Only for compatibility with
      generic matrices.

    - ``unicode`` -- boolean (default: ``False``);
      whether to use Unicode symbols instead of ASCII symbols
      for brackets and subdivision lines

    - ``shape`` -- one of ``'square'`` or ``'round'`` (default: ``None``).
      Switches between round and square brackets.
      The default depends on the setting of the ``unicode`` keyword
      argument. For Unicode symbols, the default is round brackets
      in accordance with the TeX rendering,
      while the ASCII rendering defaults to square brackets.

    - ``character_art`` -- boolean (default: ``False``); if ``True``, the
      result will be of type :class:`~sage.typeset.ascii_art.AsciiArt` or
      :class:`~sage.typeset.unicode_art.UnicodeArt` which support line
      breaking of wide matrices that exceed the window width

    - ``left_border``, ``right_border`` -- sequence (default: ``None``);
      if not ``None``, call :func:`str` on the elements and use the
      results as labels for the rows of the matrix. The labels appear
      outside of the parentheses.

    - ``top_border``, ``bottom_border`` -- sequence (default: ``None``);
      if not ``None``, call :func:`str` on the elements and use the
      results as labels for the columns of the matrix. The labels appear
      outside of the parentheses.

    EXAMPLES::

        sage: B = matrix(GF(2), 3, 3, [0, 1, 0, 0, 1, 1, 0, 0, 0])
        sage: B  # indirect doctest
        [0 1 0]
        [0 1 1]
        [0 0 0]
        sage: block_matrix([[B, 1], [0, B]])
        [0 1 0|1 0 0]
        [0 1 1|0 1 0]
        [0 0 0|0 0 1]
        [-----+-----]
        [0 0 0|0 1 0]
        [0 0 0|0 1 1]
        [0 0 0|0 0 0]
        sage: B.str(zero='.')
        '[. 1 .]\n[. 1 1]\n[. . .]'

        sage: M = matrix.identity(GF(2), 3)
        sage: M.subdivide(None, 2)
        sage: print(M.str(unicode=True, shape='square'))
        ⎡1 0│0⎤
        ⎢0 1│0⎥
        ⎣0 0│1⎦
        sage: print(unicode_art(M))  # indirect doctest
        ⎛1 0│0⎞
        ⎜0 1│0⎟
        ⎝0 0│1⎠
    """
    # Set the mapping based on keyword arguments
    # We ignore minus_one (it's only there for compatibility with Matrix)
    if (rep_mapping is not None or zero is not None or plus_one is not None
            or unicode or shape is not None or character_art
            or left_border is not None or right_border is not None
            or top_border is not None or bottom_border is not None):
        # Shunt mappings off to the generic code since they might not be
        # single characters
        return matrix_dense.Matrix_dense.str(self, rep_mapping=rep_mapping,
                                             zero=zero, plus_one=plus_one,
                                             unicode=unicode, shape=shape,
                                             character_art=character_art,
                                             left_border=left_border,
                                             right_border=right_border,
                                             top_border=top_border,
                                             bottom_border=bottom_border)

    if self._nrows == 0 or self._ncols == 0:
        return "[]"

    cdef Py_ssize_t i,j, last_i
    cdef list s = []
    empty_row = b' '*(self._ncols*2-1)
    cdef char *row_s
    cdef char *div_s

    cdef list row_div, col_div
    if self._subdivisions is not None:
        row_s = empty_row
        div_s = row_divider = b'[' + (b'-' * (self._ncols*2-1)) + b']'
        row_div, col_div = self.subdivisions()
        last_i = 0
        for i in col_div:
            if i == last_i or i == self._ncols:
                # Adjacent column divisions messy, use generic code
                return matrix_dense.Matrix_dense.str(self, rep_mapping)
            row_s[2*i-1] = c'|'
            div_s[2*i] = c'+'
            last_i = i

    for i in range(self._nrows):
        row_s = row = b'[' + empty_row + b']'
        for j in range(self._ncols):
            row_s[1+2*j] = c'0' + mzd_read_bit(self._entries, i, j)
        s.append(row)

    if self._subdivisions is not None:
        for i in reversed(row_div):
            s.insert(i, row_divider)

    return (b"\n".join(s)).decode()

