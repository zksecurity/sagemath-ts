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
