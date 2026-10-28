// Original platform APIs, deliberately independent of the TypeScript generators.
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
int main(int argc, char** argv) {
    if (argc != 4) return 2;
    const std::string name = argv[1];
    const auto seed = std::strtoull(argv[2], nullptr, 10);
    const auto count = std::strtoul(argv[3], nullptr, 10);
#ifdef _WIN32
    if (name != "msvc-rand") return 2;
    std::srand(static_cast<unsigned int>(seed));
    for (unsigned long i = 0; i < count; ++i) std::cout << std::rand() << '\n';
#else
    if (name == "lrand48" || name == "drand48-state") {
        srand48(static_cast<long>(seed));
        for (unsigned long i = 0; i < count; ++i)
            std::cout << (name == "lrand48" ? static_cast<std::uint64_t>(lrand48()) : static_cast<std::uint64_t>(drand48() * 281474976710656.0)) << '\n';
    } else if (name == "glibc-rand") {
#ifndef __GLIBC__
        return 3; // Not a glibc reference on this platform.
#else
        std::srand(static_cast<unsigned int>(seed));
        for (unsigned long i = 0; i < count; ++i) std::cout << std::rand() << '\n';
#endif
    } else return 2;
#endif
}
