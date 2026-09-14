# Direct native FLINT oracle: coefficients and errors come from the original library.
from cysignals.signals cimport sig_on, sig_off
from sage.libs.flint.types cimport nmod_poly_t
from sage.libs.flint.nmod_poly cimport (
    nmod_poly_init, nmod_poly_clear, nmod_poly_set_coeff_ui,
    nmod_poly_get_coeff_ui, nmod_poly_length, nmod_poly_xgcd,
)

def native_xgcd(a, b, unsigned long n):
    cdef nmod_poly_t A, B, G, S, T
    nmod_poly_init(A, n)
    nmod_poly_init(B, n)
    nmod_poly_init(G, n)
    nmod_poly_init(S, n)
    nmod_poly_init(T, n)
    try:
        for i, c in enumerate(a): nmod_poly_set_coeff_ui(A, i, c % n)
        for i, c in enumerate(b): nmod_poly_set_coeff_ui(B, i, c % n)
        sig_on()
        try:
            nmod_poly_xgcd(G, S, T, A, B)
        finally:
            sig_off()
        return [
            [nmod_poly_get_coeff_ui(G, i) for i in range(nmod_poly_length(G))],
            [nmod_poly_get_coeff_ui(S, i) for i in range(nmod_poly_length(S))],
            [nmod_poly_get_coeff_ui(T, i) for i in range(nmod_poly_length(T))],
        ]
    finally:
        nmod_poly_clear(A)
        nmod_poly_clear(B)
        nmod_poly_clear(G)
        nmod_poly_clear(S)
        nmod_poly_clear(T)

from sage.rings.integer cimport Integer
from sage.libs.flint.types cimport fmpz_poly_t, fmpz_t
from sage.libs.flint.flint cimport flint_free
from sage.libs.flint.fmpz_poly cimport (
    fmpz_poly_init, fmpz_poly_clear, fmpz_poly_set_str, fmpz_poly_get_str, fmpz_poly_xgcd,
)
from sage.libs.flint.fmpz cimport fmpz_init, fmpz_clear, fmpz_get_mpz

cdef poly_strings(fmpz_poly_t poly):
    cdef char * raw = fmpz_poly_get_str(poly)
    try:
        return (<bytes>raw).decode().split()[1:]
    finally:
        flint_free(raw)

def native_integer_xgcd(a, b):
    cdef fmpz_poly_t A, B, S, T
    cdef fmpz_t r
    cdef Integer answer = Integer.__new__(Integer)
    fmpz_poly_init(A)
    fmpz_poly_init(B)
    fmpz_poly_init(S)
    fmpz_poly_init(T)
    fmpz_init(r)
    try:
        raw = (str(len(a)) + '  ' + ' '.join(map(str,a))).encode()
        if fmpz_poly_set_str(A, raw): raise ValueError('invalid polynomial fixture')
        raw = (str(len(b)) + '  ' + ' '.join(map(str,b))).encode()
        if fmpz_poly_set_str(B, raw): raise ValueError('invalid polynomial fixture')
        sig_on()
        try:
            fmpz_poly_xgcd(r, S, T, A, B)
        finally:
            sig_off()
        fmpz_get_mpz(answer.value, r)
        return [str(answer), poly_strings(S), poly_strings(T)]
    finally:
        fmpz_clear(r)
        fmpz_poly_clear(A)
        fmpz_poly_clear(B)
        fmpz_poly_clear(S)
        fmpz_poly_clear(T)
