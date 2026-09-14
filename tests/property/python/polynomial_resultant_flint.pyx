# Expected results come from native FLINT, including dispatch and alias behavior.
from cysignals.signals cimport sig_on, sig_off
from sage.rings.integer cimport Integer
from sage.libs.flint.types cimport nmod_poly_t, fmpz_poly_t, fmpz_t
from sage.libs.flint.nmod_poly cimport (
    nmod_poly_init, nmod_poly_clear, nmod_poly_set_coeff_ui, nmod_poly_resultant,
)
from sage.libs.flint.fmpz_poly cimport (
    fmpz_poly_init, fmpz_poly_clear, fmpz_poly_set_str, fmpz_poly_resultant,
)
from sage.libs.flint.fmpz cimport fmpz_init, fmpz_clear, fmpz_get_mpz

def native_integer_resultant(a, b):
    cdef fmpz_poly_t A, B
    cdef fmpz_t r
    cdef Integer answer = Integer.__new__(Integer)
    fmpz_poly_init(A)
    fmpz_poly_init(B)
    fmpz_init(r)
    try:
        raw = (str(len(a)) + '  ' + ' '.join(map(str,a))).encode()
        if fmpz_poly_set_str(A, raw): raise ValueError('invalid polynomial fixture')
        raw = (str(len(b)) + '  ' + ' '.join(map(str,b))).encode()
        if fmpz_poly_set_str(B, raw): raise ValueError('invalid polynomial fixture')
        sig_on()
        try:
            fmpz_poly_resultant(r, A, B)
        finally:
            sig_off()
        fmpz_get_mpz(answer.value, r)
        return answer
    finally:
        fmpz_clear(r)
        fmpz_poly_clear(A)
        fmpz_poly_clear(B)

def native_modular_resultant(a, b, unsigned long n, same=False):
    cdef nmod_poly_t A, B
    cdef unsigned long r
    nmod_poly_init(A, n)
    nmod_poly_init(B, n)
    try:
        for i, c in enumerate(a): nmod_poly_set_coeff_ui(A, i, c % n)
        for i, c in enumerate(b): nmod_poly_set_coeff_ui(B, i, c % n)
        sig_on()
        try:
            if same:
                r = nmod_poly_resultant(A, A)
            else:
                r = nmod_poly_resultant(A, B)
        finally:
            sig_off()
        return r
    finally:
        nmod_poly_clear(A)
        nmod_poly_clear(B)
