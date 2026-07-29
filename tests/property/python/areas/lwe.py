"""SageMath side of the ``lwe`` property-test area.

The cases compare deterministic structural facts.  The sample-vector case
constructs and calls the real oracle, but compares only the vector dimension
because SageMath and sagemath-ts intentionally use different random streams.

Cases: tests/property/cases/lwe.cases.json
"""

from sage.crypto.lwe import LindnerPeikert, Regev, UniformNoiseLWE


def regev_q(n):
    return Regev(n).K.order()


def regev_dimension(n):
    return Regev(n).n


def lindner_peikert_m(n):
    return LindnerPeikert(n).m


def sample_vector_length(n):
    a, _ = Regev(n)()
    return len(a)


def uniform_noise_min_n(n):
    UniformNoiseLWE(n)
    return 'valid'


FUNCTIONS = {
    'regev_q': regev_q,
    'regev_dimension': regev_dimension,
    'lindner_peikert_m': lindner_peikert_m,
    'sample_vector_length': sample_vector_length,
    'uniform_noise_min_n': uniform_noise_min_n,
}
