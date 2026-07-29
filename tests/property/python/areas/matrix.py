"""SageMath side of the ``matrix`` property-test area.

Matrix results are transported as nested integer lists so the transcript tests
the entries without depending on SageMath's pretty-printer.

Cases: tests/property/cases/matrix.cases.json
"""

from sage.all import ZZ, matrix


def _matrix(nrows, ncols, entries):
    return matrix(ZZ, nrows, ncols, entries)


def _rows(M):
    return [[ZZ(x) for x in row] for row in M.rows()]


def determinant_2x2(a, b, c, d):
    return _matrix(2, 2, [a, b, c, d]).determinant()


def determinant_3x3(a, b, c, d, e, f, g, h, i):
    return _matrix(3, 3, [a, b, c, d, e, f, g, h, i]).determinant()


def rank_2x3(a, b, c, d, e, f):
    return _matrix(2, 3, [a, b, c, d, e, f]).rank()


def rank_3x3(a, b, c, d, e, f, g, h, i):
    return _matrix(3, 3, [a, b, c, d, e, f, g, h, i]).rank()


def hnf_2x2(a, b, c, d):
    return _rows(_matrix(2, 2, [a, b, c, d]).hermite_form())


def hnf_3x3(a, b, c, d, e, f, g, h, i):
    return _rows(_matrix(3, 3, [a, b, c, d, e, f, g, h, i]).hermite_form())


def snf_2x2(a, b, c, d):
    return _rows(_matrix(2, 2, [a, b, c, d]).smith_form(transformation=False))


def snf_3x3(a, b, c, d, e, f, g, h, i):
    return _rows(
        _matrix(3, 3, [a, b, c, d, e, f, g, h, i]).smith_form(
            transformation=False
        )
    )


def lll_2x2(a, b, c, d):
    return _rows(_matrix(2, 2, [a, b, c, d]).LLL())


def lll_3x3(a, b, c, d, e, f, g, h, i):
    return _rows(_matrix(3, 3, [a, b, c, d, e, f, g, h, i]).LLL())


FUNCTIONS = {
    'determinant_2x2': determinant_2x2,
    'determinant_3x3': determinant_3x3,
    'rank_2x3': rank_2x3,
    'rank_3x3': rank_3x3,
    'hnf_2x2': hnf_2x2,
    'hnf_3x3': hnf_3x3,
    'snf_2x2': snf_2x2,
    'snf_3x3': snf_3x3,
    'lll_2x2': lll_2x2,
    'lll_3x3': lll_3x3,
}
