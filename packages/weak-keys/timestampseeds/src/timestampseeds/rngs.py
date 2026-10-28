"""Named algorithms, not a claim to emulate any library's key generator."""
import random
from collections.abc import Iterator

WIDTHS = {"lrand48": 4, "drand48-state": 6, "glibc-rand": 4, "ansi-lcg31": 4,
          "python-mt19937": 4, "java-random": 4, "msvc-rand": 4, "cpp-mt19937": 4}
RNGS = tuple(WIDTHS)


def words(name: str, seed: int) -> Iterator[tuple[int, int]]:
    """Yield (integer, serialized byte width); no native byte order or float rounding."""
    if type(seed) is not int or seed < 0:
        raise ValueError("Seed must be a nonnegative integer")
    if name in ("lrand48", "drand48-state"):
        # POSIX srand48: low 32 seed bits become high 32 state bits.
        state = ((seed & 0xFFFFFFFF) << 16) | 0x330E
        while True:
            state = (0x5DEECE66D * state + 0xB) & ((1 << 48) - 1)
            yield (state >> 17, 4) if name == "lrand48" else (state, 6)
    elif name == "glibc-rand":
        # glibc TYPE_3 random()/rand(): degree 31, separation 3.
        seed &= 0xFFFFFFFF
        if seed == 0:
            seed = 1
        state = [seed]
        word = seed if seed < (1 << 31) else seed - (1 << 32)
        for _ in range(1, 31):
            # Park-Miller signed division must truncate toward zero, unlike //.
            hi = abs(word) // 127773 * (-1 if word < 0 else 1)
            lo = word - hi * 127773
            word = 16807 * lo - 2836 * hi
            if word < 0:
                word += 2147483647
            state.append(word & 0xFFFFFFFF)
        front, rear = 3, 0
        for i in range(310):
            state[front] = (state[front] + state[rear]) & 0xFFFFFFFF
            front, rear = (front + 1) % 31, (rear + 1) % 31
        while True:
            state[front] = (state[front] + state[rear]) & 0xFFFFFFFF
            result = state[front] >> 1
            front, rear = (front + 1) % 31, (rear + 1) % 31
            yield result, 4
    elif name == "ansi-lcg31":
        # Explicit common example constants, NOT a portable C rand() definition.
        state = seed & 0x7FFFFFFF
        while True:
            state = (1103515245 * state + 12345) & 0x7FFFFFFF
            yield state, 4
    elif name == "python-mt19937":
        # CPython integer seed expansion; NOT std::mt19937(seed).
        rng = random.Random(seed)
        while True:
            yield rng.getrandbits(32), 4
    elif name == "java-random":
        # java.util.Random(long).nextInt(), treated as an unsigned 32-bit word.
        state = (seed ^ 0x5DEECE66D) & ((1 << 48) - 1)
        while True:
            state = (0x5DEECE66D * state + 0xB) & ((1 << 48) - 1)
            yield state >> 16, 4
    elif name == "msvc-rand":
        # Microsoft CRT srand(unsigned int)/rand(), 15-bit output in four bytes.
        state = seed & 0xFFFFFFFF
        while True:
            state = (214013 * state + 2531011) & 0xFFFFFFFF
            yield (state >> 16) & 0x7FFF, 4
    elif name == "cpp-mt19937":
        # std::mt19937(integer), not seed_seq or CPython integer expansion.
        state = [seed & 0xFFFFFFFF]
        for i in range(1, 624):
            previous = state[-1]
            state.append((1812433253 * (previous ^ (previous >> 30)) + i) & 0xFFFFFFFF)
        while True:
            # In-place twist deliberately uses already updated words after wraparound.
            for i in range(624):
                joined = (state[i] & 0x80000000) | (state[(i + 1) % 624] & 0x7FFFFFFF)
                state[i] = state[(i + 397) % 624] ^ (joined >> 1) ^ (0x9908B0DF if joined & 1 else 0)
            for value in state:
                value ^= value >> 11
                value ^= (value << 7) & 0x9D2C5680
                value ^= (value << 15) & 0xEFC60000
                value ^= value >> 18
                yield value, 4
    else:
        raise ValueError(f"Unknown RNG: {name}")
